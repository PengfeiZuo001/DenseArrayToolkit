%% DenseArrayToolkit 2D Migration Example Script
% This script demonstrates a complete workflow for 3D seismic migration imaging 
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
%% step1. Setup paths and parameters
config = loadConfig();

%% step2. Extract processing parameters from configuration structure
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

%% step3. Read data
% Load seismic waveform data in SAC format from the specified directory
DataStruct = read_SAC(dataFolder);

%% step4. Preprocessing
% Apply standard seismic data preprocessing steps defined in PreprocessingParam:
% - Filtering: 
% - Demeaning:
% - Time window selection: 
% - Quality control:
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% step5. Get array and event information
% Extract metadata about the seismic array geometry and earthquake sources
% This information is crucial for:
% - Spatial analysis and quality control
% - 3D migration volume definition
% - Processing optimization for dense arrays
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis
stlo = [stationList.stlo]'; 
stla = [stationList.stla]'; 
evla = [eventList.evla]'; 
evlo = [eventList.evlo]'; 

% Create filtered list of event IDs for processing
eventid = {eventList.evid};

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% step6. Create 3D imaging grid with specified parameters
% Set up the 3D imaging grid and velocity model for migration:
dx = 4;                                                                     % Horizontal x-direction grid spacing (km)
dy = 4;                                                                     % Horizontal y-direction grid spacing (km)
dz = 1;                                                                     % Vertical grid spacing (km)
zmax = 100;                                                                 % Maximum imaging depth (km)
xpad = 40;                                                                  % Horizontal padding in x-direction (km)
ypad = 40;                                                                  % Horizontal padding in y-direction (km)

gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

%% step7. Create 3D velocity model
% to be updated because this is too slow
% nptsx = gridStruct.nx; 
% nptsy = gridStruct.ny;
nxpts = 5;
nypts = 2;
gridStruct = getVelocityModel('3D', gridStruct, nxpts,nypts);

%% step8. Update rank reduction parameters for 3D processing (not use)
RankReductionParam.nx = gridStruct.nx;   
RankReductionParam.ny = gridStruct.ny;    
RankReductionParam.ox = min(gridStruct.x);
RankReductionParam.oy = min(gridStruct.y); 
RankReductionParam.mx = max(gridStruct.x); 
RankReductionParam.my = max(gridStruct.y);
RankReductionParam.rank = 10;            
RankReductionParam.rank = 5;
RankReductionParam.fhigh = 2.4;

%% step9. Perform 3D least-squares migration
migResults = [];
ccpResults = [];
nMigratedEvents = 0;  
smoothLength = 2;                                                           % Smoothing parameter for display (grid points)
for iEvent = 1:length(eventid)
    disp(['==========> Running No. ', nMigratedEvents])
    evid = eventid{iEvent}; 
    
    %% Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    
    % Extract signal-to-noise ratios for quality control
    % snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    
    % Quality control: Skip events with insufficient station coverage or low SNR
    % Minimum minTrace stations required to ensure reliable 3D imaging results
    if length(gather) < MigParam.minRatio*length(stla)
        continue
    end
    
    %% Compute receiver functions through deconvolution
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);

    %% Apply rank reduction preprocessing to improve signal quality or Radon 
    % This helps suppress noise and enhance coherent signals
    if MigParam.ssa
        [gather, reconed_rf] = rankReduction_new(gather, gridStruct, RankReductionParam);
        MigParam.reconed_rf = reconed_rf;
    end

    if MigParam.is_radon
        RadonParam.highs = 1.2;                                             % High-slowness cutoff (s/km)
        RadonParam.pmax = 0.06;                                             % Maximum slowness (s/km)
        RadonParam.pmin = -0.06;                                            % Minimum slowness (s/km)
        RadonParam.N1 = 10;
        RadonParam.plotRadon = 0;
        gather = radonTransform2D(gather, gridStruct, RadonParam);
    end
    
    %% Perform 3D least-squares migration
    migResult = leastSquaresMig3D(gather, gridStruct, MigParam);
    migResults = [migResults; migResult];
    
    %% Perform CCP stacking for comparison with migration results
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    %% single
    % plotCCPMigrationResults3D(ccpResult,migResult, gridStruct,smoothLength,dem)

    nMigratedEvents = nMigratedEvents + 1;
end

%% step10. Stacking
ccpImage = stackImagingResults(ccpResults);
migImage = stackImagingResults(migResults);

%% step11. Visualization
plotCCPMigrationResults3D(ccpImage,migImage, gridStruct,smoothLength,dem)

%% step12. Save results
% write_MigResult([config.outputFolder, '/migResults3D.mat'], migResults);
% write_MigResult([config.outputFolder, '/ccpResults3D.mat'], ccpResults);
