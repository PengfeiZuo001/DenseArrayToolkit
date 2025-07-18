d%% DenseArrayToolkit 2D Migration Main Function Example
% This script demonstrates a complete workflow for 2D seismic imaging using both
% Common Conversion Point (CCP) stacking and migration methods. It processes
% seismic array data to create subsurface images along 2D profiles.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata and filter by azimuth
%   4. Create velocity model - Set up 2D velocity structure for migration
%   5. Migration imaging - Perform CCP stacking and least-squares migration
%   6. Visualization - Display and compare different imaging results
%   7. Save results - Store processed data and images
%
% This script implements both conventional and least-squares migration approaches
% alongside CCP stacking for comparative analysis of imaging methods.

clear; clc; close all;

%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions and
% dependencies to the MATLAB path. This ensures access to all required
% processing routines in the DenseArrayToolkit.
setupPaths();

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
% - Azimuthal consistency filtering
% - Migration path calculations
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis and filtering
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% Filter events based on azimuthal consistency
% Only keep events within specified maximum azimuth difference (config.max_angle_diff)
% This ensures 2D approximation validity for the imaging profile
idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% Create filtered list of event IDs that meet azimuthal criteria
eventid = {eventList.evid};
eventid = eventid(idxConsistentEQ);

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Set up the imaging grid and velocity model for migration:
% 1. Define grid spacing in horizontal (dx) and vertical (dy) directions
% 2. Create imaging grid based on array geometry
% 3. Generate 1D velocity model for migration
% Note: Using 1D velocity model for simplified ray tracing in 2D migration
dx = 5;  % Horizontal grid spacing (km)
dy = 5;  % Vertical grid spacing (km)
gridStruct = createGrid(DataStruct, dx, dy);
% Create 1D velocity model for migration imaging
gridStruct = getVelocityModel('1D',gridStruct);
%% 5. Migration imaging
% Perform both conventional migration and CCP stacking for comparison:
% - Standard Kirchhoff migration
% - Least-squares migration for improved resolution
% - CCP stacking as a reference method
%
% Initialize arrays to store results from all events:
dmig   = [];  % Standard migration results
dmigls = [];  % Least-squares migration results
dimg   = [];  % CCP stacking results
nMigratedEvents = 1;    % Counter for successfully processed events

% Process each event that meets the azimuthal filtering criteria
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    % Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    
    % Quality control: Skip events with insufficient station coverage
    % Minimum 50 stations required to ensure reliable imaging results
    if length(gather) < 50
        continue
    end
  
    % Optional Radon Transform filtering for enhanced signal quality
    % This helps suppress noise and improve coherency in the data
    if DeconvParam.radonfilter
        % Apply Radon Transform for array-based signal enhancement
        gather = radonTransform(gather, RadonParam);
    end

    % Compute receiver functions through deconvolution
    % This isolates converted phases from the P-wave coda
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);
    
    % Perform 2D CCP stacking
    CCPParam.imagingType = '2D';  % Set imaging mode to 2D
    % Apply Common Conversion Point stacking
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    % Normalize CCP image by hit count to account for uneven sampling
    dimg(:,:,nMigratedEvents) = ccpResult.img./max(ccpResult.count,1);

    % Perform least-squares migration
    % This method provides improved resolution compared to standard migration
    migResult = leastSquaresMig(gather, gridStruct, MigParam);
    
    % Store migration results for current event
    dmig(:, :, nMigratedEvents)   = migResult.mig;
    dmigls(:, :, nMigratedEvents) = migResult.migls;
    
    % Clear figure windows to avoid memory issues during batch processing
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end

% Extract spatial coordinates from final migration result for visualization
x = migResult.x;    % Horizontal distance coordinate
z = migResult.z;    % Depth coordinate

%% 6. Visualization
% Create comparative display of different imaging methods:
% - CCP stacking results (dimg)
% - Standard migration results (dmig)
% - Least-squares migration results (dmigls)
% This allows direct comparison of the different imaging approaches
plotCCPMigrationResults(dimg,dmig,dmigls,x,z)

%% 7. Save results
% Store final migration and CCP results to files for future analysis
% Results include:
% - Migration results (conventional and least-squares)
% - CCP stacking results
% - Grid coordinates and parameters
% - Processing configuration
write_MigResult([config.outputFolder,'/migResult.mat'], migResult);
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);
