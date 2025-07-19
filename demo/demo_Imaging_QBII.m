%% DenseArrayToolkit Common Conversion Point (CCP) Stacking Main Function Example - QBII Array
% This script demonstrates the complete workflow for seismic imaging using receiver
% functions and Common Conversion Point stacking method. The workflow processes
% seismic data from the QBII array deployment to image subsurface structures.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata and geometry
%   4. Create velocity model - Set up 3D velocity structure for migration
%   5. CCP imaging - Perform Common Conversion Point stacking
%   6. Visualization - Display results with multiple cross-sections
%   7. Save results - Store processed data and images
%
% This script focuses on processing data from the QBII array deployment,
% providing detailed subsurface imaging through CCP stacking analysis.

clear; clc; close all;

%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions to MATLAB path
% and loading configuration parameters for various processing steps. This ensures
% access to all required processing routines in the DenseArrayToolkit.
setupPaths();

% Load configuration file containing essential parameters for data processing,
% including paths, preprocessing settings, and imaging parameters
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

dataFolder = './data/event_waveforms_QBII';
%% 1. Read data
% Load seismic waveform data from QBII array in SAC format
% The data is encapsulated into a structured array (DataStruct) containing
% waveforms and metadata for each recording
DataStruct = read_SAC(dataFolder);

%% 2. Preprocessing
% Apply standard seismic data preprocessing steps defined in PreprocessingParam:
% - Filtering: Remove unwanted frequency components
% - Demeaning: Remove DC offset from signals
% - Time window selection: Extract relevant portion of seismograms
% - Quality control: Remove noisy or incomplete recordings
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Get array and event information
% Extract and organize metadata about the seismic array geometry and earthquake
% source locations. This information is crucial for:
% - Spatial analysis of the seismic array
% - Ray path calculations
% - Migration of receiver functions
% - Quality control based on source-receiver geometry
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis and filtering
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% List of unique event identifiers
eventid = {eventList.evid};

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Generate a 3D velocity model for migration imaging:
% 1. Create an imaging grid based on array geometry using principal component
%    analysis to optimize grid orientation
% 2. Set grid spacing (dx, dy) for horizontal dimensions
% 3. Define vertical sampling for depth profiles
% 4. Generate 3D velocity structure for accurate imaging
dx = 5;  % Horizontal grid spacing (km)
dy = 5;  % Vertical grid spacing (km)
gridStruct = createGrid(DataStruct, dx, dy);
% Create 3D velocity model with 5-point sampling for detailed structure
gridStruct = getVelocityModel('3D',gridStruct,5);
%% 5. Migration imaging
% Perform Common Conversion Point (CCP) stacking for each event that meets quality
% criteria. This process includes:
% 1. Extract event gathers (recordings of the same event at different stations)
% 2. Quality control: Skip events with insufficient station coverage (<60 stations)
% 3. Optional Radon Transform for array processing
% 4. Receiver function calculation through deconvolution
% 5. CCP stacking to create 3D image volume

% Initialize arrays for accumulating migration results across all events:
dimg = [];    % 4D array to store CCP stacking results (X, Y, Z, Event)
count = [];   % 4D array tracking number of traces contributing to each point
nMigratedEvents = 1;    % Counter for successfully processed events

% Process each event that meets the quality criteria
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    % Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    
    % Quality control: Skip events with insufficient station coverage
    % Minimum 60 stations required to ensure reliable imaging results
    if length(gather) < 60
        continue
    end
  
    % Optional Radon Transform filtering (currently disabled)
    % Can be enabled to enhance signal coherency across the array
%     gather = radonTransform(gather, RadonParam);
    
    % Configure deconvolution parameters for receiver function computation
    DeconvParam.radonfilter = false;  % Disable Radon filtering
    DeconvParam.verbose = false;       % Suppress detailed processing output
    % Compute receiver functions through deconvolution
    % This isolates converted phases from the P-wave coda
    gather = deconv(gather, DeconvParam);
    
    % Perform Common Conversion Point stacking
    % This maps receiver function amplitudes to subsurface points
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    dimg(:,:,:,nMigratedEvents) = ccpResult.img;        % Store image volume
    count(:,:,:,nMigratedEvents) = ccpResult.count;     % Store hit count

    % Clear figure windows to avoid memory issues during batch processing
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Visualization
% Create visualizations of the 3D CCP stacking results:
% 1. Extract grid coordinates and compute normalized stacked volume
% 2. Generate 3D volume plot with orthogonal slices
% 3. Create cross-sections along predefined profiles
% 4. Apply custom colormap for optimal visualization

% Prepare final image volume
X = ccpResult.X;  % X coordinates of the imaging grid
Y = ccpResult.Y;  % Y coordinates of the imaging grid
Z = ccpResult.Z;  % Depth coordinates
% Compute final image by stacking across events and normalizing by hit count
V = sum(dimg,4)./max(sum(count,4),1);

% Create 3D visualization with orthogonal slices
figure;
set(gcf,'Position',[50 50 800 800],'Color','w')
h = slice(X,Y,Z,V,90,90,40);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
colormap(flipud(roma));
cmax = rms(V(:));
caxis([-cmax cmax]);

% Define and plot multiple cross-sections through the volume
% These profiles are chosen to highlight key structural features
x1 = 0;
y1 = 0;
x2 = 30;
y2 = 90;
x3 = 90;
y3 = 100;
x4 = 130;
y4 = 160;
x5 = 310;
y5 = 30;
profile1=[x1,y1;x2,y2];
profile2=[x2,y2;x3,y3];
profile3=[x3,y3;x4,y4];
profile4=[x4,y4;x5,y5];
profile = {profile1,profile2,profile3,profile4};
plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile)
%% 7. Save results
% Save the final CCP imaging results to a MAT file for future analysis
% The saved results include:
% - 3D image volume
% - Grid coordinates
% - Processing parameters
% - Hit count distribution
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);
