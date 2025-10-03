%% DenseArrayToolkit 2D Migration Main Function Example
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

dataFolder = './data/event_waveforms_QBI';
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
dx = 4;  % Horizontal grid spacing (km)
dy = 4;  % Vertical grid spacing (km)
dz = 1;
zmax = 100;
xpad = 40;
ypad = 40;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);
% Create 1D velocity model for migration imaging
gridStruct = getVelocityModel('2D',gridStruct);
%% 5. Migration imaging
% Perform both conventional migration and CCP stacking for comparison:
% - Standard Kirchhoff migration
% - Least-squares migration for improved resolution
% - CCP stacking as a reference method
%
% Initialize arrays to store results from all events: 
migResults = []; % Migration results
ccpResults = []; % CCP stacking results

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
  
    % Compute receiver functions through deconvolution
    DeconvParam.gauss=2.5;
    DeconvParam.verbose = false;
    DeconvParam.radonfilter = false;
    gather = deconv(gather, DeconvParam);   

    % Optional Radon Transform filtering for enhanced signal quality
    % This helps suppress noise and improve coherency in the data
    % Check if Radon filtering is enabled for enhanced signal quality
    DeconvParam.radonfilter = true;
    if DeconvParam.radonfilter
        RadonParam.highs = 1.2;
        RadonParam.pmax = 0.02;
        RadonParam.pmin = -0.02;
        % Apply Radon Transform for array processing
        gatherRadon = radonTransform(gather, gridStruct, RadonParam);
% 
%         figure;
%         set(gcf,'Position',[100 100 800 400],'Color','w')
%         ax1 = subplot(121);
%         plotCommonEventGather(gather, [], 'trace', 'imagesc',ax1)
%         set(ax1,'YLim',[0 20])
%         ax2 = subplot(122);
%         plotCommonEventGather(gatherRadon, [], 'trace', 'imagesc',ax2)
%         set(ax2,'YLim',[0 20])
    end    
    
    % Perform 2D CCP stacking
    CCPParam.imagingType = '2D';  % Set imaging mode to 2D
    CCPParam.smoothLength = 0;
    % Apply Common Conversion Point stacking
    ccpResult = CCPCommonEventGather(gatherRadon, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    % Perform least-squares migration
    % This method provides improved resolution compared to standard migration
    MigParam.itermax = 30;
    MigParam.gauss = DeconvParam.gauss;
    migResult = leastSquaresMig(gatherRadon, gridStruct, MigParam);
    
    % Store migration results for current event
    migResults = [migResults; migResult];
    
    % Clear figure windows to avoid memory issues during batch processing
%     close all; 
    nMigratedEvents = nMigratedEvents + 1;
end

%% 6. Visualization
% Create comparative display of different imaging methods:
% - CCP stacking results (ccpResults)
% - Standard migration results (migResults)
% - Least-squares migration results (migResults)
% This allows direct comparison of the different imaging approaches
smoothLength = 3;
plotCCPMigrationResults(ccpResults,migResults,gridStruct,smoothLength)

%% 7. Save results
% Store final migration and CCP results to files for future analysis
% Results include:
% - Migration results (conventional and least-squares)
% - CCP stacking results
% - Grid coordinates and parameters
% - Processing configuration
write_MigResult([config.outputFolder,'/migResults.mat'], migResults);
write_MigResult([config.outputFolder,'/ccpResults.mat'], ccpResults);
