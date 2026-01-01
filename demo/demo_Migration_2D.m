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
dataFolder = '../data/event_waveforms_QBI';
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
dx = 4;  % Horizontal grid spacing in x-direction (km)
dy = 4;  % Horizontal grid spacing in y-direction (km) 
dz = 1;  % Vertical grid spacing (km)
zmax = 100; % Maximum depth for imaging (km)
xpad = 40;  % Padding distance beyond array extent in x-direction (km)
ypad = 40;  % Padding distance beyond array extent in y-direction (km)
%
% The createGrid function generates a 2D/3D grid structure based on:
% - Array geometry from station coordinates
% - Specified grid spacing and depth range
% - Padding to ensure complete coverage of the imaging domain
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 1D velocity model for migration imaging
% For 2D migration, a 1D velocity model is typically sufficient and provides
% computational efficiency while maintaining reasonable accuracy for teleseismic
% imaging. The velocity model includes P-wave and S-wave velocities as a function
% of depth, which are essential for calculating travel times and migration operators
gridStruct = getVelocityModel('2D',gridStruct);
%% 5. Migration imaging
% Core processing loop: Perform 2D migration and CCP stacking for each event
% This section implements the main imaging algorithms to create subsurface
% images from receiver function data. Two complementary methods are applied:
%
% - CCP Stacking: Traditional method that stacks receiver functions along
%   theoretical conversion points using 1D ray tracing
% - Least-squares Migration: Advanced inverse method that improves resolution
%   by accounting for limited aperture and uneven station coverage
%
% Results from all events are accumulated for final stacking and comparison
migResults = []; % Array to store migration results from all events
ccpResults = []; % Array to store CCP stacking results from all events

nMigratedEvents = 0; % Counter for successfully processed events

minTrace = 50;
% Process each event that passed azimuthal filtering
% The loop iterates through all events with consistent back-azimuths
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; % Current event ID
    
    % Extract Common Event Gather: all recordings of this event across stations
    % This groups seismograms from the same earthquake recorded at different
    % stations, enabling event-based processing
    gather = getCommonEventGather(DataStruct, evid);

    % Quality control: Skip events with insufficient station coverage
    % A minimum of 50 stations is required to ensure reliable imaging results
    % and adequate spatial sampling for migration
    if length(gather) < minTrace
        continue % Skip to next event if insufficient stations
    end

    % Compute receiver functions through iterative deconvolution
    % Receiver functions isolate the P-to-S converted phases that reveal
    % subsurface discontinuities. Gaussian filtering (gauss=2.5) controls
    % the frequency content and resolution of the resulting images
    DeconvParam.gauss = 2.5;    % Gaussian width parameter for frequency filtering
    DeconvParam.verbose = false; % Suppress verbose output during processing
    gather = deconv(gather, DeconvParam);

    % Apply Radon Transform for enhanced signal-to-noise ratio
    % The Radon transform helps suppress coherent noise and improve
    % signal coherency across the array by focusing energy along moveout curves
    RadonParam.highs = 1.2;   % High-slowness cutoff (s/km)
    RadonParam.pmax = 0.02;   % Maximum slowness (s/km)
    RadonParam.pmin = -0.02;  % Minimum slowness (s/km)
    gatherRadon = radonTransform2D(gather, gridStruct, RadonParam);

    % Perform 2D Common Conversion Point (CCP) stacking
    % CCP stacking bins receiver functions based on their theoretical conversion
    % points in the subsurface and stacks them to create a migrated image
    CCPParam.imagingType = '2D';    % Set imaging mode to 2D
    CCPParam.smoothLength = 0;      % No additional smoothing applied
    ccpResult = CCPCommonEventGather(gatherRadon, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult]; % Accumulate CCP results

    % Perform least-squares migration for improved resolution
    % This method solves an inverse problem to account for limited aperture
    % and uneven station coverage, providing sharper images than standard migration
    MigParam.itermax = 30;          % Maximum iterations for convergence
    MigParam.gauss = DeconvParam.gauss; % Use same Gaussian parameter as deconvolution
    migResult = leastSquaresMig2D(gatherRadon, gridStruct, MigParam);

    % Store migration results for current event
    migResults = [migResults; migResult]; % Accumulate migration results

    % Update event counter for progress tracking
    nMigratedEvents = nMigratedEvents + 1;
end

%% 6. Visualization
% Create comprehensive comparative display of different imaging methods
% This visualization allows direct assessment of the relative performance
% and characteristics of CCP stacking versus least-squares migration
%
% The plotCCPMigrationResults function generates side-by-side comparisons of:
% - CCP stacking results: Traditional binning and stacking method
% - Least-squares migration: Advanced inverse imaging with improved resolution
%
% Key visualization features:
% - Common color scale for direct amplitude comparison
% - Grid coordinates for spatial reference
% - Smoothing applied to enhance interpretability while preserving features
smoothLength = 3; % Smoothing parameter for display (grid points)
plotCCPMigrationResults(ccpResults, migResults, gridStruct, smoothLength);

%% 7. Save results
% Store final migration and CCP results to files for future analysis and sharing
% This preserves the processed data, enabling:
% - Reproducible research and validation
% - Further analysis without reprocessing
% - Comparison with other datasets or methods
%
% Saved data includes:
% - Migration results: Complete least-squares migration outputs for all events
% - CCP stacking results: Traditional CCP images for all events
% - Grid coordinates: Spatial reference system for interpretation
% - Processing parameters: Configuration used for reproducibility
% - Quality metrics: Processing statistics and event counts
write_MigResult([config.outputFolder, '/migResults.mat'], migResults);
write_MigResult([config.outputFolder, '/ccpResults.mat'], ccpResults);
