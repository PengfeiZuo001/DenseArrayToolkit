%% DenseArrayToolkit Common Conversion Point (CCP) Stacking Main Function Example
% This script demonstrates the complete workflow for seismic imaging using receiver
% functions and Common Conversion Point stacking method. The workflow processes
% seismic data from the QBI array to image subsurface structures.
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

dataFolder = './data/event_waveforms_QBI';
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
dx = 5;
dy = 5;
gridStruct = createGrid(DataStruct, dx, dy);
% Create or obtain 3D velocity model with specified sampling (10 points)
% This model will be used for ray tracing in migration
npts = 10;
gridStruct = getVelocityModel('3D',gridStruct,npts);
%% 5. Migration imaging
% Perform Common Conversion Point (CCP) stacking for each event that meets quality
% criteria. This process includes:
% 1. Extract event gathers (recordings of the same event at different stations)
% 2. Quality control: Skip events with insufficient station coverage (<60 stations)
% 3. Optional Radon Transform for array processing
% 4. Receiver function calculation through deconvolution
% 5. CCP stacking to create 3D image volume

% Initialize arrays for accumulating migration results across all events
dimg = [];    % 4D array to store CCP stacking results (X, Y, Z, Event)
count = [];   % 4D array tracking number of traces contributing to each point
nMigratedEvents = 1;   % Counter for successfully processed events
minTrace = 60; % Minimum number of traces required for CCP imaging
minSNR = 3;    % Minimum SNR of the RF required for CCP imaging
DeconvParam.radonfilter = 0; % Apply radon transform to seismograms
% Process each event that meets the filtering criteria
for iEvent = 1:length(eventid)
    evid = eventid{iEvent};
    % Extract seismic records for current event (Common Event Gather)
    gather = getCommonEventGather(DataStruct, evid);
    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));
    % Skip events with fewer than 60 valid stations to ensure imaging quality
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end

    % Check if Radon filtering is enabled for enhanced signal quality
    if DeconvParam.radonfilter
        RadonParam.highs = 1.2;
        RadonParam.pmax = 0.05;
        RadonParam.pmin = -0.05;
        % Apply Radon Transform for array processing
        gather = radonTransform(gather, RadonParam);
    end

    % Perform deconvolution to obtain receiver functions
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);

    % Apply Common Conversion Point stacking for the current event gather
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    dimg(:,:,:,nMigratedEvents) = ccpResult.img;
    count(:,:,:,nMigratedEvents) = ccpResult.count;

    % Clear figure windows to avoid memory issues during batch processing
    close all;
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Visualization
% Create visualizations of the CCP stacking results using the visualizeCCPResults
% function, which provides:
% 1. 3D volume visualization with topography overlay
% 2. Interactive profile selection
% 3. Structure-oriented filtering of cross-sections
% 4. Comparative display of original and filtered profiles

% Prepare the normalized image volume
ccpResult.V = sum(dimg,4)./max(sum(count,4),1);

% Configure visualization options
options = struct();
options.profileType = 'predefined';  % Enable interactive profile selection
options.smoothingParams = struct(...
    'radius', 3, ...        % Smoothing radius
    'eps', 0.01, ...       % Regularization parameter
    'order', 2);           % Smoothing order
% define profile location
profilePoints = [
91.461	37.47
91.948	38.168
93.599	38.782
93.994	38.681
95.53531 39.21347];
[px,py] = latlonToProjectedCoords(profilePoints(:,1), profilePoints(:,2), gridStruct);
options.profilePoints(:,1) = px;
options.profilePoints(:,2) = py;
% Generate visualizations and get profile data structure
profileStruct = visualizeCCPResults(ccpResult, gridStruct, options);
%% 7. Save results
% Save the final CCP imaging results to a MAT file for future reference and analysis.
% The saved results include:
% - 3D image volume
% - Grid coordinates
% - Processing parameters
% - Other relevant metadata
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);
