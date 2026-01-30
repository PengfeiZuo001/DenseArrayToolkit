%% DenseArrayToolkit 3D CCP Imaging Example - BaiyanEbo Region
% This script demonstrates 3D Common Conversion Point (CCP) stacking for
% seismic imaging in the BaiyanEbo region. It processes receiver functions
% to create volumetric subsurface images showing crustal and upper mantle
% structures beneath the dense seismic array.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata for 3D processing
%   4. Create velocity model - Set up 3D velocity structure for CCP stacking
%   5. Compute receiver functions - Extract P-to-S converted phases
%   6. CCP stacking - Perform 3D Common Conversion Point imaging
%   7. Results output - Visualize and save 3D imaging results
%
% This script implements 3D CCP stacking with rank reduction preprocessing
% for improved signal-to-noise ratio in dense array applications.

clear; clc; close all;
%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions and
% dependencies to the MATLAB path. This ensures access to all required
% processing routines in the DenseArrayToolkit.
cd ../
setupPaths();

% Load configuration file containing essential parameters for data processing,
% including paths, processing parameters, and imaging settings
config = loadConfig();

% Extract processing parameters from configuration structure
% - dataFolder: Directory containing seismic data files
% - PreprocessingParam: Parameters for filtering, windowing, and quality control
% - MigParam: Parameters controlling migration algorithms (not used in CCP)
% - RadonParam: Settings for Radon transform array processing (optional)
% - DeconvParam: Parameters for receiver function deconvolution
% - CCPParam: Settings for Common Conversion Point stacking
% - RankReductionParam: Parameters for rank reduction preprocessing
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;

%% 1. Read data
% Load seismic waveform data in SAC format from the specified directory
% The data is encapsulated into a structured array (DataStruct) containing
% waveforms and metadata for each recording
dataFolder1= './data/event_waveforms_BY';
DataStruct1 = read_SAC(dataFolder1);

dataFolder2= '../DenseArrayToolkit-data/event_waveforms_BY_local';
DataStruct2 = read_SAC(dataFolder2);

DataStruct = [DataStruct1 DataStruct2];
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
% - 3D CCP volume definition
% - Processing optimization for dense arrays
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% Optional azimuthal filtering (commented out for this region)
% idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% Create filtered list of event IDs for processing
eventid = {eventList.evid};
% eventid = eventid(idxConsistentEQ);

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Set up the 3D imaging grid and velocity model for CCP stacking:
% 1. Define grid spacing in x, y, and z directions
% 2. Create 3D imaging volume based on array geometry
% 3. Generate 3D velocity model for accurate depth conversion
% Note: Using 3D velocity model for precise ray tracing and depth conversion
dx = 5;    % Horizontal x-direction grid spacing (km)
dy = 5;    % Horizontal y-direction grid spacing (km)
dz = 0.5;   % Vertical grid spacing (km) - finer resolution for CCP
zmax = 100; % Maximum imaging depth (km)
xpad = 40;
ypad = 40;
% Create 3D imaging grid with specified parameters
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 3D velocity model for CCP imaging
% npts = 10 specifies the number of interpolation points for velocity model
npts = 5;
gridStruct = getVelocityModel('3D', gridStruct, npts);

%% 5. Compute receiver functions
% Apply deconvolution to extract receiver functions from seismic waveforms
% This isolates P-to-S converted phases that contain information about
% subsurface discontinuities beneath the array
DeconvParam.verbose = 0;
DataStruct = deconv(DataStruct, DeconvParam);

%% 6. CCP stacking
% Perform 3D Common Conversion Point stacking to create volumetric images:
% - Stack receiver functions at their theoretical conversion points
% - Account for 3D ray geometry and velocity structure
% - Apply rank reduction preprocessing for noise suppression
%
% Initialize arrays to store CCP results from all events:
ccpResults = [];    % 3D CCP image results
nMigratedEvents = 0;    % Counter for successfully processed events

minTrace = 100; % Minimum number of traces required for CCP imaging
minSNR = 1;    % Minimum SNR of the RF required for CCP imaging

% Process each event in the dataset
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    
    % Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    
    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    % Quality control: Skip events with insufficient station coverage
    % Minimum minTrace stations required to ensure reliable 3D imaging results
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end

    % Apply deconvolution to extract receiver functions
    % This isolates P-to-S converted phases from the P-wave coda
    gather = deconv(gather, DeconvParam);

    % Perform 3D CCP stacking using Common Conversion Point method
    % This maps receiver functions to their theoretical conversion points
    CCPParam.gauss=DeconvParam.gauss;
    ccpResult = FresnelCCPCommonEventGather(gather, gridStruct, CCPParam);
    
    % Store CCP results for current event
    ccpResults = [ccpResults; ccpResult];
    
    % Close all figure windows to avoid memory issues during batch processing
%     close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end

%% 7. Results output
% Create final 3D CCP image by stacking all events and normalizing by hit count
% This produces the final volumetric image showing subsurface structure
ccpImage = stackImagingResults(ccpResults);

% Configure visualization options
options = struct();
options.profileType = 'predefined';  % Enable interactive profile selection
% options.smoothingParams = struct(...
%     'radius', 3, ...        % Smoothing radius
%     'eps', 0.01, ...       % Regularization parameter
%     'order', 2);           % Smoothing order
% N-S profile crossing Baiyan Obo minning area
options.profilePoints(:,1) = [76.6927 76.6927 nan 0 200];
options.profilePoints(:,2) = [0 150 nan 66.9016 66.9016];
options.dem = load('./visualization/Baiyanebo_DEM_small.mat');
% E-W profile crossing Baiyan Obo minning area
% options.profilePoints(:,1) = [0; 200];
% options.profilePoints(:,2) = [66.9016; 66.9016];
visualizeImage(ccpImage, gridStruct, options);

% Optional: Save results to file for future analysis
% save './results/BaiyanEbo_ccp.mat' 'X' 'Y' 'Z' 'V' 'gridStruct'
