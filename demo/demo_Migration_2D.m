%% DenseArrayToolkit 2D Migration Example Script
% This script demonstrates a complete workflow for 2D seismic migration imaging 
% using dense array data. 
%
% The workflow includes two complementary imaging methods:
%   - Common Conversion Point (CCP) stacking: 
%   - Least-squares migration: 

% Copyright (C) Zhejiang University, seismology group
% References:
%   
clear; clc; close all;

% input your elevation data if needed
dem = load('Qaidam_DEM.mat');
%% 0. Setup paths and parameters
% if this is first time to run the code, please
% cd to main directory, and run:
%--------------------------------------------------------------------------
% setupPaths();
%--------------------------------------------------------------------------

%% step1. Load configuration file containing essential parameters for data processing
% Almost all parameters are included in this file
config = loadConfig();

%% step2. Extract processing parameters from configuration structure for clarity
% These parameters control different stages of the seismic imaging workflow:
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;

%% step3. Read data
% Load seismic waveform data in SAC.
% DataStruct output contains:
% - Waveform traces (seismograms) for each station-event pair
% - SAC header information
% - Station metadata for array geometry analysis
DataStruct = read_SAC(dataFolder);

%% step4. Preprocessing
% The preprocessing function:
% - Applies bandpass filtering to isolate frequency bands of interest
% - Removes DC offset and linear trends from signals
% - Selects appropriate time windows around seismic phases
% - Performs quality control to remove noisy or corrupted traces
% - Normalizes amplitudes for consistent processing
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% step5. Get array and event information
% Extract metadata about the seismic array geometry and earthquake sources
% - stationList: 
% - eventList: 
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis and filtering
% These coordinates are used to calculate azimuths and distances for
% quality control and 2D profile selection
stlo = [stationList.stlo]';                                                 % Station longitude (degrees)
stla = [stationList.stla]';                                                 % Station latitude (degrees)
evla = [eventList.evla]';                                                   % Event epicenter latitude (degrees)
evlo = [eventList.evlo]';                                                   % Event epicenter longitude (degrees)

%% step6. Filter events based on azimuthal consistency to ensure 2D approximation validity
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

%% step7. Create meshgrid
% Set up the imaging grid and velocity model for 2D migration processing
% This step defines the spatial domain where subsurface structures will be imaged
% and establishes the velocity structure used for ray tracing and migration
% Grid parameters:
dx = 4;                                                                     % Horizontal grid spacing in x-direction (km)
dy = 4;                                                                     % Horizontal grid spacing in y-direction (km) 
dz = 1;                                                                     % Vertical grid spacing (km)
zmax = 100;                                                                 % Maximum depth for imaging (km)
xpad = 40;                                                                  % Padding distance beyond array extent in x-direction (km)
ypad = 40;                                                                  % Padding distance beyond array extent in y-direction (km)

% The createGrid function generates a 2D/3D grid structure based on:
% - Array geometry
% - Padding 
GridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

%% step8. Create velocity model for migration imaging
% For 2D migration, velocity models are extracted from 3D models. So both 1D and 2D
% are available
nxpts = 5;                                                                  % points along x
nypts = 2;                                                                  % points along y
GridStruct = getVelocityModel('2D',GridStruct,nxpts,nypts);

%% step9. RF imaging for each event
migResults = [];
ccpResults = []; 
nMigratedEvents = 0; 
smoothLength = 2;                                                           % Smoothing parameter for display (grid points)
for iEvent = 1:length(eventid)
    disp(['==========> Running No. ', nMigratedEvents])
    evid = eventid{iEvent}; % Current event ID
    
    %% Extract Common Event Gather: all recordings of this event across stations
    gather = getCommonEventGather(DataStruct, evid);

    % the number of stations are 
    if length(gather) < MigParam.minRatio*length(stla)
        continue 
    end

    %% Compute receiver functions through iterative deconvolution
    DeconvParam.verbose = false; 
    gather = deconv(gather, DeconvParam);

    %% Apply Radon Transform for enhanced signal-to-noise ratio
    % The Radon transform helps suppress coherent noise and improve
    % signal coherency across the array by focusing energy along moveout curves
    if MigParam.is_radon
        RadonParam.highs = 1.2;                                             % High-slowness cutoff (s/km)
        RadonParam.pmax = 0.06;                                             % Maximum slowness (s/km)
        RadonParam.pmin = -0.06;                                            % Minimum slowness (s/km) (also see loadConfig())
        RadonParam.N1 = 10;
        RadonParam.plotRadon = 0;
        gather = radonTransform2D(gather, GridStruct, RadonParam);
    end

    %% Perform 2D Common Conversion Point (CCP) stacking
    CCPParam.imagingType = '2D';                                            % Set imaging mode to 2D
    CCPParam.smoothLength = 0;                                              % No additional smoothing applied
    ccpResult = CCPCommonEventGather(gather, GridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult]; 

    %% Perform least-squares migration for improved resolution
    migResult = leastSquaresMig2D(gather, GridStruct, MigParam);
    migResults = [migResults; migResult]; 

    %% plot single event
    % plotCCPMigrationResults(ccpResult, migResult, GridStruct, smoothLength, dem);

    % Update event counter for progress tracking
    nMigratedEvents = nMigratedEvents + 1;
end

%% step10. Stacking
ccpImage = stackImagingResults(ccpResults);
migImage = stackImagingResults(migResults);

%% step10. Visualization
plotCCPMigrationResults(ccpImage, migImage, GridStruct, smoothLength, dem);

%% step12. Save results
% write_MigResult([config.outputFolder, '/migResults.mat'], migResults);
% write_MigResult([config.outputFolder, '/ccpResults.mat'], ccpResults);
