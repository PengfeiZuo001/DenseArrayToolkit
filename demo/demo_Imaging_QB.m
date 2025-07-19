%% DenseArrayToolkit Common Conversion Point (CCP) Stacking Main Function Example
% This script demonstrates the complete workflow for seismic imaging using receiver
% functions and Common Conversion Point stacking method. The workflow processes
% seismic data from dense arrays to image subsurface structures.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata and geometry
%   4. Create velocity model - Set up 3D velocity structure for migration
%   5. CCP imaging - Perform Common Conversion Point stacking
%   6. Visualization - Display results using various plotting methods
%   7. Save results - Store processed data and images
%
% This script is designed for analyzing receiver functions from dense seismic arrays
% to image subsurface discontinuities and structures.

clear; clc; close all;
load roma;
%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions to MATLAB path
% and loading configuration parameters for various processing steps. The config file
% contains crucial parameters that control data processing, imaging, and analysis.
% Call custom function setupPaths() to add dependency packages and functions to search path
setupPaths();

% Load various parameters from config file (such as data path, processing parameters, etc.)
config = loadConfig();

% Extract processing parameters from configuration structure
% - dataFolder: Directory containing seismic data files
% - PreprocessingParam: Parameters for filtering, windowing, and quality control
% - MigParam: Parameters controlling migration process
% - RadonParam: Settings for Radon transform array processing
% - DeconvParam: Parameters for receiver function deconvolution
% - CCPParam: Settings for Common Conversion Point stacking
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;

% Define data directories for two phases of the experiment (QBI and QBII)

dataFolder1 = './data/event_waveforms_QBI';
dataFolder2 = './data/event_waveforms_QBII';
%% 1. Read data
% Load seismic waveform data from SAC format files. The data is organized in two
% directories (QBI and QBII) containing event-based waveforms. The read_SAC function
% creates a structured array containing waveform data and metadata for each recording.
% Read SAC format seismic data from dataFolder and encapsulate into DataStruct
DataStruct1 = read_SAC(dataFolder1);
DataStruct2 = read_SAC(dataFolder2);
DataStruct = [DataStruct1 DataStruct2];
%% 2. Preprocessing
% Apply standard seismic data preprocessing steps defined in PreprocessingParam:
% - Filtering: Remove unwanted frequency components
% - Demeaning: Remove DC offset from signals
% - Time window selection: Extract relevant portion of seismograms
% - Quality control: Remove noisy or incomplete recordings
% These steps prepare the data for receiver function analysis and imaging.
% Process DataStruct according to configured preprocessing parameters (filtering, demeaning, slicing, etc.)
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Get array and event information
% Extract and organize metadata about the seismic array geometry and earthquake
% source locations. This information is crucial for:
% - Spatial analysis of the seismic array
% - Ray path calculations
% - Migration of receiver functions
% - Quality control based on source-receiver geometry
% Extract station list and event list from DataStruct
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract latitude and longitude information of stations and events for subsequent consistency filtering
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% Event ID list
eventid = {eventList.evid};

% Generate event-station correspondence table for quick indexing later
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Generate a 3D velocity model for migration imaging:
% 1. Create an imaging grid based on array geometry using principal component
%    analysis to optimize grid orientation
% 2. Set grid spacing (dx, dy) for horizontal dimensions
% 3. Define vertical sampling (npts) for depth profiles
% 4. Interpolate velocity structure onto the imaging grid
% This velocity model will be used to migrate receiver functions in the CCP stacking.
% Perform principal component analysis based on station locations to automatically create imaging grid
dx = 10;
dy = 10;
dz = 1;
zmax = 100;
xpad = 50;
ypad = 50;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);
% Create or obtain velocity model for use in subsequent imaging
npts = 10;
gridStruct = getVelocityModel('3D',gridStruct,npts);
%% 5. Migration imaging
% Perform Common Conversion Point (CCP) stacking for each event that meets quality
% criteria. This process includes:
% 1. Extract event gathers (recordings of the same event at different stations)
% 2. Quality control: Skip events with insufficient station coverage (<60 stations)
% 3. Optional Radon Transform for array processing (currently commented out)
% 4. Receiver function calculation through deconvolution
% 5. CCP stacking to create 3D image volume
%
% The results are accumulated across all events to enhance signal-to-noise ratio
% and improve imaging quality. Two key matrices are maintained:
% - dimg: Stores the stacked amplitudes
% - count: Tracks the number of traces contributing to each image point
% Initialize arrays for accumulating migration results across all events
dimg = [];    % 4D array to store CCP stacking results (X, Y, Z, Event)
count = [];   % 4D array tracking number of traces contributing to each point
nMigratedEvents = 1;    % Counter for successfully processed events
minTrace = 60; % Minimum number of traces required for CCP imaging
minSNR = 5;    % Minimum SNR of the RF required for CCP imaging
% Iterate through all events that meet the filtering criteria
for iEvent = 1:length(eventid)
    evid = eventid{iEvent};
    % Extract the subset of seismic records corresponding to the current event (Common Event Gather)
    gather = getCommonEventGather(DataStruct, evid);

    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));
    % Skip events with fewer than 60 valid stations to ensure imaging quality
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end

    % Configure and perform deconvolution to compute receiver functions
    % - radonfilter: Optional array processing step (disabled)
    % - verbose: Control detailed output during processing
    DeconvParam.radonfilter = false;  % Disable additional Radon filtering
    DeconvParam.verbose = false;       % Suppress detailed processing output
    % Transform seismic recordings into receiver functions through deconvolution
    % This highlights converted phases from subsurface discontinuities
    gather = deconv(gather, DeconvParam);

    % Call CCPCommonEventGather for common conversion point stacking imaging
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    dimg(:,:,:,nMigratedEvents) = ccpResult.img;
    count(:,:,:,nMigratedEvents) = ccpResult.count;

    % Close all figure windows to avoid generating too many windows during batch processing
    close all;

    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Visualization
% Create visualizations of the 3D CCP stacking results:
% 1. Normalize the stacked image by the hit count
% 2. Generate a 3D volume plot with orthogonal slices
% 3. Apply custom colormap (roma) for optimal visualization
% 4. Scale color limits based on RMS amplitude
% 5. Option to select profile lines for cross-sectional views
%
% Note: The profile selection section (ginput) is dependent on user interaction
% and will create cross-sections along specified profiles using plotCCPXsectionCartesian
% Prepare final image volume for visualization
X = ccpResult.X;  % X coordinates of the imaging grid
Y = ccpResult.Y;  % Y coordinates of the imaging grid
Z = ccpResult.Z;  % Depth coordinates
% Compute final image by stacking across events and normalizing by hit count
% This reduces artifacts and enhances true structural features
V = sum(dimg,4)./max(sum(count,4),1);
% Create 3D visualization with orthogonal slices
figure;
% Set up figure window with appropriate size and white background
set(gcf,'Position',[50 50 800 800],'Color','w')
% Display orthogonal slices through the volume at X=90km, Y=90km, Z=40km
h = slice(X,Y,Z,V,90,90,40); hold on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
colormap(flipud(roma));
cmax = rms(V(:));
caxis([-cmax cmax]);
view(0,90)
axis xy;
axis equal;
% [xpoints,ypoints] = ginput;
% define profile location
profilePoints = [
91.461	37.47
91.948	38.168
93.599	38.782
93.994	38.681
95.53531 39.21347];
% coordinates conversion from lat, lon to the projected coordinate
[xpoints,ypoints] = latlonToProjectedCoords(profilePoints(:,1), profilePoints(:,2), gridStruct);

plot(xpoints, ypoints, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(xpoints, ypoints, 'r-', 'LineWidth', 2);

profileAll={};
for n = 1:length(xpoints)-1
    profile=[xpoints(n),ypoints(n);xpoints(n+1),ypoints(n+1)];
    profileAll{n} = profile;
end
plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profileAll)
%% 7. Save results
% Save the final CCP imaging results to a MAT file for future reference and analysis.
% The saved results include:
% - 3D image volume
% - Grid coordinates
% - Processing parameters
% - Other relevant metadata
% Write ccpResult to the specified file
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);
