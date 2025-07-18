%% demo_rankReduction.m
%  This script demonstrates a complete workflow for receiver function computation 
%  and "Rank Reduction". The process includes reading data, computing receiver 
%  functions, multi-event reconstruction and stacking, and final visualization.
%  Main steps include:
%    0. Setup paths and parameters
%    1. Read data
%    2. Preprocessing
%    3. Compute receiver functions
%    4. Perform rank reduction (DRR-OTG) for each event
%    5. Stack receiver functions

clear; clc; close all;

%% 0. Setup paths and parameters
% --------------------------------------------------
try
    setupPaths();  % Add project toolboxes or codes to MATLAB path
catch ME
    warning('setupPaths() not found or failed: %s\nUsing current path instead.', ME.message);
end

% Load configuration
if exist('loadConfig','file') == 2
    config = loadConfig();
else
    error('No loadConfig() found. Please implement or provide config structure.');
end

% Extract commonly used fields from config
% dataFolder = config.dataFolder;
dataFolder = './data/event_waveforms_BY';
PreprocessingParam = config.PreprocessingParam;
DeconvParam        = config.DeconvParam;
RankReductionParam = config.RankReductionParam;

%% 1. Read data
% --------------------------------------------------
fprintf('\n[Step 1] Reading SAC data from folder: %s\n', dataFolder);
DataStruct = read_SAC(dataFolder);
if isempty(DataStruct)
    error('No data read from %s. Check if files exist.', dataFolder);
end

%% 2. Preprocessing
% --------------------------------------------------
fprintf('\n[Step 2] Preprocessing data...\n');
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Compute receiver functions
% --------------------------------------------------
fprintf('\n[Step 3] Computing receiver functions (deconvolution)...\n');
DeconvParam.verbose = false;
DataStruct = deconv(DataStruct, DeconvParam);

%% 4. Perform rank reduction (DRR-OTG) for each event
% --------------------------------------------------
fprintf('\n[Step 4] Doing rank reduction (DRR-OTG) for each event...\n');

gridStruct = createGrid(DataStruct);

eventList = getEvents(DataStruct);  % Returns struct array of event info
eventIDs  = {eventList.evid};
DataStructDRR = [];  % Store reconstruction results

minStationCount = 50; % Threshold: skip if event has fewer stations than this

for iEvent = 1:length(eventIDs)
    evid = eventIDs{iEvent};
    [gather, matchIndex] = getCommonEventGather(DataStruct, evid);

    if length(gather) < minStationCount
        fprintf('Event %s only has %d stations (<%d). Skipped.\n', ...
            evid, length(gather), minStationCount);
        continue;
    end

    % Call rankReduction for 3D data reconstruction of this event
    % gather => [ gather(i).RF.itr ...], gather(i).StationInfo ...
    [gatherReconstructed, d1_otg] = rankReduction(gather, gridStruct, RankReductionParam);

    % Merge results into DataStruct_drr
    DataStructDRR = [DataStructDRR; gatherReconstructed(:)];

    % For debugging/visualization of d1_otg, implement here or in rankReduction
    close all;  % To prevent too many figure windows
end

% Transpose results back to DataStruct-like shape (if needed)
DataStructDRR = DataStructDRR';

%% 5. Stack receiver functions
% --------------------------------------------------
fprintf('\n[Step 5] Stacking receiver functions by station...\n');

% stackCommonStationGather() assumes stacking multiple events for the same station
% Output: seisout => stacked data, depth0 => depth axis, mohoStruct => optional Moho depth
[seisout, depth0, mohoStruct] = stackCommonStationGather(DataStructDRR);

fprintf('\nDone! All steps completed.\n');

%% 6. Structural oriented filter (post-processing)
% Get station info
stationList = getStations(DataStructDRR);
% stationList should have fields .stlo, .stla
% Flatten stlo, stla to vector if needed
stlo = [stationList.stlo]';  
stla = [stationList.stla]';

% transform lat lon to x y (relative to param.lonmin, param.latmin)
 [rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct);
% shift to ensure min coordinate=0
rx = rx - min(rx);
ry = ry - min(ry);
param = RankReductionParam;
% store for reference
param.x = rx;  % unregular location 
param.y = ry;

% Define grid
param.ox = 0;  % origin x
param.oy = 0;  % origin y
param.mx = ceil(max(rx));  % max value of x
param.my = ceil(max(ry));  % max value of y


dt = 0.1;

d0 = seisout;
[d1_otg, d1] = drr3drecon_otg(...
    d0, rx, ry, ...
    param.nx, param.ny, ...
    param.ox, param.oy, param.mx, param.my, ...
    param.flow, param.fhigh, dt, ...
    param.rank, param.K, param.niter, param.eps, ...
    param.verb, param.mode);

% 3D slope calculation (inline and xline)
% default parameter values are suitable for most cases
[dipi,dipx] = str_dip3d(d1_otg);

% Structural smoothing
r1=2;
r2=2;
eps=0.01;
order=2;

d1_otg_str=str_pwsmooth_lop3d(d1_otg,dipi,dipx,r1,r2,eps,order);

figure;
subplot(211)
imagesc(reshape(d1_otg,201,11*14))
subplot(212)
imagesc(reshape(d1_otg_str,201,11*14))