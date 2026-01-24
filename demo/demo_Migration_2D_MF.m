%% DenseArrayToolkit 2D Migration Example Script
% This script demonstrates a complete workflow for 2D seismic migration imaging
% using dense array data. It processes teleseismic receiver functions to create
% high-resolution subsurface images along 2D profiles.
%
% The workflow includes two complementary imaging methods:
%   - Common Conversion Point (CCP) stacking: Traditional method for imaging
%     seismic discontinuities by stacking receiver functions along conversion points
%   - Least-squares migration: Advanced method that improves resolution by
%     solving an inverse problem to account for limited aperture and uneven coverage
%
% Main processing steps:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for receiver function analysis
%   3. Get array and event information - Extract metadata and filter by azimuth
%   4. Create velocity model - Set up 2D velocity structure for migration
%   5. Migration imaging - Perform CCP stacking and least-squares migration
%   6. Visualization - Display and compare different imaging results
%   7. Save results - Store processed data and images for further analysis
%
% This implementation allows direct comparison between conventional CCP stacking
% and advanced migration techniques, providing insights into their relative
% strengths for subsurface imaging.

clear; clc; close all;

%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary toolbox functions
% and dependencies to the MATLAB path. This ensures access to all required
% processing routines in the DenseArrayToolkit.
cd ..
setupPaths();

% Load configuration file containing essential parameters for data processing
% The configuration includes paths, processing parameters, and imaging settings
% that control the entire migration workflow
config = loadConfig();

% Extract processing parameters from configuration structure for clarity
% These parameters control different stages of the seismic imaging workflow:
% - dataFolder: Directory containing seismic data files in SAC format
% - PreprocessingParam: Parameters for data filtering, windowing, and quality control
% - MigParam: Parameters controlling migration algorithms and convergence criteria
% - RadonParam: Settings for Radon transform array processing and noise suppression
% - DeconvParam: Parameters for receiver function deconvolution and Gaussian filtering
% - CCPParam: Settings for Common Conversion Point stacking and binning parameters
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;

% Override data folder path for this specific demonstration
% This points to the Qaidam Basin dataset for 2D migration example
dataFolder = '../DenseArrayToolkit-data/event_waveforms_MF';
%% 1. Read data
% Load seismic waveform data in SAC (Seismic Analysis Code) format from the 
% specified directory. The read_SAC function reads both waveform data and 
% metadata (headers) for each recording station.
%
% DataStruct output contains:
% - Waveform traces (seismograms) for each station-event pair
% - SAC header information including timing, coordinates, and event parameters
% - Station metadata for array geometry analysis
DataStruct = read_SAC(dataFolder);

%% 2. Preprocessing
% Apply standard seismic data preprocessing steps to prepare waveforms for 
% receiver function analysis and migration. The preprocessing function:
%
% - Applies bandpass filtering to isolate frequency bands of interest
% - Removes DC offset and linear trends from signals
% - Selects appropriate time windows around seismic phases
% - Performs quality control to remove noisy or corrupted traces
% - Normalizes amplitudes for consistent processing
%
% This step ensures clean, consistent data suitable for deconvolution and imaging
PreprocessingParam.highs = 5.0;
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Get array and event information
% Extract metadata about the seismic array geometry and earthquake sources
% This information is essential for spatial analysis, quality control, and
% ensuring the validity of the 2D migration approximation
%
% stationList: Contains station coordinates and metadata
% eventList: Contains earthquake source parameters and locations
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis and filtering
% These coordinates are used to calculate azimuths and distances for
% quality control and 2D profile selection
stlo = [stationList.stlo]';  % Station longitude (degrees)
stla = [stationList.stla]';  % Station latitude (degrees)
evla = [eventList.evla]';    % Event epicenter latitude (degrees)
evlo = [eventList.evlo]';    % Event epicenter longitude (degrees)

% Filter events based on azimuthal consistency to ensure 2D approximation validity
% The 2D migration approach requires events to be approximately aligned along
% a single azimuth. This filtering removes events with scattered back-azimuths
% that would violate the 2D imaging assumption
idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% Create filtered list of event IDs that meet azimuthal criteria
% Only events with consistent back-azimuths (within config.max_angle_diff) are retained
eventid = {eventList.evid};
eventid = eventid(idxConsistentEQ);

% Generate event-station correspondence table for efficient data access
% This table maps each event to the stations that recorded it, facilitating
% rapid data retrieval during the migration loop
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Set up the imaging grid and velocity model for 2D migration processing
% This step defines the spatial domain where subsurface structures will be imaged
% and establishes the velocity structure used for ray tracing and migration
%
% Grid parameters:
dx = 0.5;  % Horizontal grid spacing in x-direction (km)
dy = 0.5;  % Horizontal grid spacing in y-direction (km) 
dz = 0.2;  % Vertical grid spacing (km)
zmax = 30; % Maximum depth for imaging (km)
xpad = 5;  % Padding distance beyond array extent in x-direction (km)
ypad = 5;  % Padding distance beyond array extent in y-direction (km)
%
% The createGrid function generates a 2D/3D grid structure based on:
% - Array geometry from station coordinates
% - Specified grid spacing and depth range
% - Padding to ensure complete coverage of the imaging domain
GridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 1D velocity model for migration imaging
% For 2D migration, a 1D velocity model is typically sufficient and provides
% computational efficiency while maintaining reasonable accuracy for teleseismic
% imaging. The velocity model includes P-wave and S-wave velocities as a function
% of depth, which are essential for calculating travel times and migration operators
GridStruct = getVelocityModel('3D',GridStruct,5);
%% 5. Migration imaging
% Core processing loop: Perform CCP stacking for each event
% This section implements the main imaging algorithms to create subsurface
% images from receiver function data.
%
% - CCP Stacking: Traditional method that stacks receiver functions along
%   theoretical conversion points using 1D ray tracing

% Results from all events are accumulated for final stacking and comparison
ccpResults = []; % Array to store CCP stacking results from all events

nMigratedEvents = 0; % Counter for successfully processed events

minTrace = 60;
minSNR = 1;    % Minimum SNR of the RF required for CCP imaging
goodEvent = {};
% Process each event that passed azimuthal filtering
% The loop iterates through all events with consistent back-azimuths
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; % Current event ID
    
    % Extract Common Event Gather: all recordings of this event across stations
    % This groups seismograms from the same earthquake recorded at different
    % stations, enabling event-based processing
    gather = getCommonEventGather(DataStruct, evid);

    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    % Quality control: Skip events with insufficient station coverage
    % Minimum minTrace stations required to ensure reliable 3D imaging results

    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end
    goodEvent = [goodEvent;evid];

    % Compute receiver functions through iterative deconvolution
    % Receiver functions isolate the P-to-S converted phases that reveal
    % subsurface discontinuities. Gaussian filtering (gauss=2.5) controls
    % the frequency content and resolution of the resulting images
    DeconvParam.gauss = 10;    % Gaussian width parameter for frequency filtering
    DeconvParam.verbose = false; % Suppress verbose output during processing
    gather = deconv(gather, DeconvParam);

    % Apply Radon Transform for enhanced signal-to-noise ratio
    % The Radon transform helps suppress coherent noise and improve
    % signal coherency across the array by focusing energy along moveout curves
    RadonParam.highs = DeconvParam.gauss/2;  % High-slowness cutoff (s/km)
    RadonParam.pmax = 0.06;   % Maximum slowness (s/km)
    RadonParam.pmin = -0.06;  % Minimum slowness (s/km)
    RadonParam.N1 = 10;
    RadonParam.plotRadon = 0; 
%     gatherRadon = radonTransform2D(gather, GridStruct, RadonParam);
    gatherRadon = gather;
    
    % Perform 2D Common Conversion Point (CCP) stacking
    % CCP stacking bins receiver functions based on their theoretical conversion
    % points in the subsurface and stacks them to create a migrated image
    CCPParam.imagingType = '3D';    % Set imaging mode to 2D
    CCPParam.smoothLength = 0;      % No additional smoothing applied
    ccpResult = CCPCommonEventGather(gatherRadon, GridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult]; % Accumulate CCP results

    % Update event counter for progress tracking
    nMigratedEvents = nMigratedEvents + 1;
end

%% 6. Visualization
% Create final 3D CCP image by stacking all events and normalizing by hit count
% This produces the final volumetric image showing subsurface structure
ccpImage = stackImagingResults(ccpResults,3);

% Configure visualization options
options = struct();
options.profileType = 'interactive';  % Enable interactive profile selection
options.smoothingParams = struct(...
    'radius', 3, ...        % Smoothing radius
    'eps', 0.01, ...       % Regularization parameter
    'order', 2);           % Smoothing order

% options.profilePoints(:,1) = [76.6927 76.6927 nan 0 200];
% options.profilePoints(:,2) = [0 150 nan 66.9016 66.9016];
options.dem = load('./visualization/Qaidam_DEM.mat');

% E-W profile crossing Baiyan Obo minning area
% options.profilePoints(:,1) = [0; 200];
% options.profilePoints(:,2) = [66.9016; 66.9016];
visualizeImage(ccpImage, GridStruct, options);