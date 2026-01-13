%% DenseArrayToolkit Common Conversion Point (CCP) Stacking Main Function Example
% This script demonstrates the complete workflow for seismic imaging using receiver
% functions and Common Conversion Point stacking method. The workflow processes
% seismic data from the Sichuan-Longmenshan array to image subsurface structures.
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
% This script focuses on processing data from the QBI array deployment for
% subsurface imaging and structural analysis.

clear; clc; close all;
cd ..
%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions to MATLAB path
% and loading configuration parameters for various processing steps.
setupPaths();

% Load configuration file containing essential parameters for data processing
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

dataFolder = '../DenseArrayToolkit-data/event_waveforms_SL';
%% 1. Read data
% Load seismic waveform data in SAC format from the QBI array deployment
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
% 3. Define velocity structure using a 3D model
dx = 4;
dy = 4; 
dz = 1;
zmax = 100;
xpad = 30;
ypad = 30;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);
% Create or obtain 3D velocity model with specified sampling (10 points)
% This model will be used for ray tracing in migration
npts = 5;
gridStruct = getVelocityModel('3D',gridStruct,npts);

% Set up 3D migration parameters based on grid structure
MigParam.paramMig = setMigParam3D(gridStruct);
%% 5. Migration imaging
% Perform Common Conversion Point (CCP) stacking for each event that meets quality
% criteria. This process includes:
% 1. Extract event gathers (recordings of the same event at different stations)
% 2. Quality control: Skip events with insufficient station coverage (<60 stations)
% 3. Optional Radon Transform for array processing
% 4. Receiver function calculation through deconvolution
% 5. CCP stacking to create 3D image volume

% Initialize arrays for accumulating migration results across all events
nMigratedEvents = 1;   % Counter for successfully processed events
minTrace = 60; % Minimum number of traces required for CCP imaging
minSNR = 5;    % Minimum SNR of the RF required for CCP imaging
DeconvParam.radonfilter = 0; % Apply radon transform to seismograms
DeconvParam.gauss = 2.5;
CCPParam.smoothLength = 3;
% Process each event that meets the filtering criteria
ccpResults = [];
ccpResultsRadon = [];
migResults = [];

for iEvent = 1:length(eventid)
    evid = eventid{iEvent};
    % Extract seismic records for current event (Common Event Gather)
    gather = getCommonEventGather(DataStruct, evid);
    DeconvParam.radonfilter = false;
    DeconvParam.verbose = 0;
    
    gather = deconv(gather, DeconvParam);
%     plotCommonEventGather(gather, [], 'trace', 'wigb')

    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));
    % Skip events with fewer than 60 valid stations to ensure imaging quality
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end
    % Check if Radon filtering is enabled for enhanced signal quality
    DeconvParam.radonfilter = true;
    if DeconvParam.radonfilter
        RadonParam.highs = 1.2;
        RadonParam.pmax = 0.02;
        RadonParam.pmin = -0.02;
        % Apply Radon Transform for array processing
        gatherRadon = radonTransform2D(gather, gridStruct, RadonParam);
        export_fig(['./figures/SL_radon_post_rfs_',num2str(evid),'.png'],'-r300')

        CCPParam.smoothLength = 3;
        ccpResultRadon = CCPCommonEventGather(gatherRadon, gridStruct, CCPParam);
        ccpResultsRadon = [ccpResultsRadon; ccpResultRadon];
    end    
    
    % Apply Common Conversion Point stacking for the current event gather
    CCPParam.smoothLength = 3;
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    % Apply migration imaging
    MigParam.itermax = 10;
    MigParam.tmax = 80;
    migResult = leastSquaresMig3D(gather, gridStruct, MigParam);
    migResults = [migResults; migResult];
%      pause;
    % Clear figure windows to avoid memory issues during batch processing
    close all;
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Stacking
% Combine imaging results from all processed events using stackImagingResults
% This function averages migration results and normalizes CCP results by count
% to create consolidated 3D image volumes for visualization
smoothLength = 3;
ccpImage = stackImagingResults(ccpResults,smoothLength);
migImage = stackImagingResults(migResults);
ccpImageRadon = stackImagingResults(ccpResultsRadon,smoothLength);
%% 7. Visualization
% Display 3D imaging results using interactive visualization tools
% This includes volume slicing and cross-section profiling for both
% migration and CCP stacking results

% Configure visualization options for predefined profiles
options = struct();
options.profileType = 'interactive';  % Use predefined profile paths
options.dem = load('./visualization/SichuanLongmenshan_DEM.mat');

% Visualize stacked migration and CCP results
visualizeImage(migImage, gridStruct, options);
% visualizeImage(ccpImageRadon, gridStruct, options);