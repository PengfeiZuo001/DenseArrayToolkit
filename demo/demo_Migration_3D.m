%% DenseArrayToolkit 3D Migration Main Function Example
% This script demonstrates a complete workflow for 3D seismic imaging using
% least-squares migration methods. It processes seismic array data to create
% volumetric subsurface images in 3D space.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata for 3D processing
%   4. Create velocity model - Set up 3D velocity structure for migration
%   5. Migration imaging - Perform 3D least-squares migration
%   6. Stacking - Combine results from multiple events
%   7. Visualization - Display 3D imaging results
%
% This script implements 3D least-squares migration with rank reduction
% preprocessing for improved imaging quality in dense array applications.

clear; clc; close all;

%% 0. Setup paths and parameters
% Load configuration file containing essential parameters for data processing,
% including paths, processing parameters, and imaging settings
config = loadConfig();

% Extract processing parameters from configuration structure
% - dataFolder: Directory containing seismic data files
% - PreprocessingParam: Parameters for filtering, windowing, and quality control
% - MigParam: Parameters controlling migration algorithms
% - RadonParam: Settings for Radon transform array processing
% - DeconvParam: Parameters for receiver function deconvolution
% - CCPParam: Settings for Common Conversion Point stacking
% - RankReductionParam: Parameters for rank reduction preprocessing
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;
RankReductionParam = config.RankReductionParam;

%% 1. Read data
% Load seismic waveform data in SAC format from the specified directory
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
% Extract metadata about the seismic array geometry and earthquake sources
% This information is crucial for:
% - Spatial analysis and quality control
% - 3D migration volume definition
% - Processing optimization for dense arrays
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% Create filtered list of event IDs for processing
eventid = {eventList.evid};

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Set up the 3D imaging grid and velocity model for migration:
% 1. Define grid spacing in x, y, and z directions
% 2. Create 3D imaging volume based on array geometry
% 3. Generate 3D velocity model for migration
% Note: Using 3D velocity model for accurate ray tracing in volumetric imaging
dx = 10;    % Horizontal x-direction grid spacing (km)
dy = 10;    % Horizontal y-direction grid spacing (km)
dz = 1;     % Vertical grid spacing (km)
zmax = 100; % Maximum imaging depth (km)
xpad = 20;  % Horizontal padding in x-direction (km)
ypad = 20;  % Horizontal padding in y-direction (km)

% Create 3D imaging grid with specified parameters
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 3D velocity model for migration imaging
nptsx = gridStruct.nx;
nptsy = gridStruct.ny;
% nptsx = 5;
% nptsy = 4;

gridStruct = getVelocityModel('3D', gridStruct, nptsx,nptsy);

%% Update rank reduction parameters for 3D processing
% Configure rank reduction parameters based on the 3D grid structure
% This ensures the rank reduction processing aligns with the imaging volume
RankReductionParam.nx = gridStruct.nx;    % Number of grid points in x-direction
RankReductionParam.ny = gridStruct.ny;    % Number of grid points in y-direction
RankReductionParam.ox = min(gridStruct.x); % Minimum x-coordinate of grid
RankReductionParam.oy = min(gridStruct.y); % Minimum y-coordinate of grid
RankReductionParam.mx = max(gridStruct.x); % Maximum x-coordinate of grid
RankReductionParam.my = max(gridStruct.y); % Maximum y-coordinate of grid
RankReductionParam.rank = 10;             % Initial rank for reduction
RankReductionParam.rank = 5;
RankReductionParam.fhigh = 2.4;
%% 5. Migration imaging
% Perform 3D least-squares migration with rank reduction preprocessing:
% - Standard 3D migration
% - 3D least-squares migration for improved resolution
% - Rank reduction preprocessing for noise suppression
%
% Initialize structures to store migration results from all events
migResults = [];
ccpResults = [];
nMigratedEvents = 0;    % Counter for successfully processed events

% Set up 3D migration parameters based on grid structure
MigParam.paramMig = setMigParam3D(gridStruct);

minTrace = 100; % Minimum number of traces required for migration imaging
minSNR = 3;    % Minimum SNR of the RF required for migration imaging

% Process each event in the dataset
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    
    % Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    
    % Extract signal-to-noise ratios for quality control
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    
    % Quality control: Skip events with insufficient station coverage or low SNR
    % Minimum minTrace stations required to ensure reliable 3D imaging results
    if length(gather) < minTrace || mean(snrAll) < minSNR
        continue
    end
    
    % Compute receiver functions through deconvolution
    % This isolates converted phases from the P-wave coda
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);

    % Apply rank reduction preprocessing to improve signal quality
    % This helps suppress noise and enhance coherent signals for better imaging
    RankReductionParam.plotRankReduction = 1;
    MigParam.tmax = 80;

    [gatherDRR, dout_regular] = rankReduction_new(gather, gridStruct, RankReductionParam);
    MigParam.paramMig.dout_regular = dout_regular;
    
    % Perform 3D least-squares migration
    % This method provides improved resolution compared to standard migration
    migResult = leastSquaresMig3D(gatherDRR, gridStruct, MigParam);

    % Store migration results for current event
    % mig - Standard 3D migration results
    % migls - 3D least-squares migration results
    migResults = [migResults; migResult];
    
    % Perform CCP stacking for comparison with migration results
    ccpResult = CCPCommonEventGather(gatherDRR, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    % Optional pause for visualization (can be removed for batch processing)
%     pause;
    
    % Close all figure windows to avoid memory issues during batch processing
    close all;
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Stacking
% Combine imaging results from all processed events using stackImagingResults
% This function averages migration results and normalizes CCP results by count
% to create consolidated 3D image volumes for visualization
ccpImage = stackImagingResults(ccpResults);
migImage = stackImagingResults(migResults);

%% 7. Visualization
% Display 3D imaging results using interactive visualization tools
% This includes volume slicing and cross-section profiling for both
% migration and CCP stacking results

% Configure visualization options for predefined profiles
options = struct();
options.profileType = 'predefined';  % Use predefined profile paths

% Define North-South profile crossing Baiyan Obo mining area
options.profilePoints(:,1) = [76.6927   76.6927 nan 0       200];
options.profilePoints(:,2) = [0         150     nan 66.9016 66.9016];
options.dem = load('../visualization/Baiyanebo_DEM.mat');

% Visualize stacked migration and CCP results
visualizeImage(migImage, gridStruct, options);
visualizeImage(ccpImage, gridStruct, options);