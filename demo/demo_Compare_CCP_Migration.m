%% DenseArrayToolkit CCP vs Migration Comparison Example
% This script provides a direct comparison between 3D CCP stacking and 
% 3D least-squares migration imaging methods using identical preprocessing
% and data processing workflows. The goal is to enable fair comparison
% between the two imaging approaches for dense array applications.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply identical filters and data preparation
%   3. Get array and event information - Extract metadata for 3D processing
%   4. Create velocity model - Set up 3D velocity structure for both methods
%   5. Compute receiver functions - Extract P-to-S converted phases
%   6. Parallel imaging - Apply both CCP stacking and migration to same data
%   7. Results comparison - Visualize and compare imaging results
%
% This script ensures identical preprocessing for both methods to enable
% fair comparison of imaging quality and resolution characteristics.

clear; clc; close all;
cd ../

%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions and
% dependencies to the MATLAB path. This ensures access to all required
% processing routines in the DenseArrayToolkit.
setupPaths();

% Load configuration file containing essential parameters for data processing
config = loadConfig();

% Extract processing parameters from configuration structure
% Both CCP and migration will use identical preprocessing parameters
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;
RankReductionParam = config.RankReductionParam;

%% 1. Read data
% Load seismic waveform data in SAC format from the specified directory
% The same dataset will be used for both CCP and migration imaging
DataStruct = read_SAC(dataFolder);

%% 2. Preprocessing
% Apply identical preprocessing steps to ensure fair comparison:
% - Filtering: Remove unwanted frequency components
% - Demeaning: Remove DC offset from signals
% - Time window selection: Extract relevant portion of seismograms
% - Quality control: Remove noisy or incomplete recordings
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Get array and event information
% Extract metadata about the seismic array geometry and earthquake sources
% This information is used consistently for both imaging methods
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
% Set up the 3D imaging grid and velocity model for both methods:
% Using identical grid parameters ensures spatial consistency between results
dx = 10;    % Horizontal x-direction grid spacing (km)
dy = 10;    % Horizontal y-direction grid spacing (km)
dz = 0.5;   % Vertical grid spacing (km) - fine resolution for both methods
zmax = 100; % Maximum imaging depth (km)

% Create 3D imaging grid with identical parameters for both methods
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax);

% Create 3D velocity model for both CCP and migration imaging
npts = 10;
gridStruct = getVelocityModel('3D', gridStruct, npts);

%% 5. Compute receiver functions
% Apply deconvolution to extract receiver functions from seismic waveforms
% Identical receiver functions will be used for both imaging methods
DataStruct = deconv(DataStruct, DeconvParam);

%% 6. Parallel imaging - CCP and Migration
% Apply both imaging methods to identical preprocessed data:
% - CCP stacking: Maps receiver functions to conversion points
% - Migration: Least-squares migration for improved resolution
% Both methods use identical preprocessing and rank reduction
MigParam.paramMig = setMigParam3D(gridStruct);e
% Initialize storage arrays for both methods
ccp_img = [];      % 3D CCP image results
ccp_count = [];    % 3D CCP hit count for normalization
mig_img = [];      % 3D migration image results
mig_imgls = [];    % 3D least-squares migration results
nProcessedEvents = 1;  % Counter for successfully processed events

minTrace = 100; % Minimum number of traces required for CCP imaging
minSNR = 5;    % Minimum SNR of the RF required for CCP imaging

% Process each event with both imaging methods
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    
    % Extract seismic records for current event
    gather = getCommonEventGather(DataStruct, evid);
    
    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    % Quality control: Skip events with insufficient station coverage
    % Minimum minTrace stations required to ensure reliable 3D imaging results
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end

    % Apply deconvolution to extract receiver functions
    gather = deconv(gather, DeconvParam);
    
    % Apply identical rank reduction preprocessing to both methods
    RankReductionParam.rank = 10;
    [gatherReconstructed, d1_otg] = rankReduction(gather, gridStruct, RankReductionParam);

    % === CCP STACKING ===
    % Perform 3D Common Conversion Point stacking
    ccpResult = CCPCommonEventGather(gatherReconstructed, gridStruct, CCPParam);
    ccp_img(:,:,:,nProcessedEvents) = ccpResult.img;
    ccp_count(:,:,:,nProcessedEvents) = ccpResult.count;
    
    % === MIGRATION IMAGING ===
    % Perform 3D least-squares migration
    migResult = leastSquaresMig3D(gatherReconstructed, gridStruct, MigParam);
    mig_img(:,:,:,nProcessedEvents) = migResult.mig;
    mig_imgls(:,:,:,nProcessedEvents) = migResult.migls;
    
    % Progress indicator
    fprintf('Processed event %d/%d: %s\n', iEvent, length(eventid), evid);
    
    % Close figure windows to avoid memory issues
    close all;
    
    nProcessedEvents = nProcessedEvents + 1;
end

%% 7. Results comparison
% Create normalized final images for both methods and perform comparison

% Extract grid coordinates (identical for both methods)
X = ccpResult.X;  % X-coordinate grid (km)
Y = ccpResult.Y;  % Y-coordinate grid (km)
Z = ccpResult.Z;  % Z-coordinate grid (km)

% === CCP Results ===
% Create normalized CCP image by stacking all events
ccp_final = sum(ccp_img,4)./max(sum(ccp_count,4),1);

% === Migration Results ===
% Create normalized migration images by stacking all events
mig_final = sum(mig_img,4)./nProcessedEvents;
migls_final = sum(mig_imgls,4)./nProcessedEvents;

%% 8. Visualization and comparison
% Create comprehensive comparison plots

% === 3D CCP Image ===
figure('Name', 'CCP Stacking Results', 'Position', [100, 100, 800, 600]);
subplot(1,3,1);
V_ccp = mean(ccp_img,4);
h1 = slice(X,Y,Z,V_ccp,90,90,40);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
set(h1(:),'EdgeColor','none'); set(gca,'ZDir','reverse');
colormap(flipud(roma)); cmax = rms(V_ccp(:)); caxis([-cmax cmax]);
title('CCP Stacking - Average Image');
colorbar;

% === 3D Migration Image ===
subplot(1,3,2);
V_mig = mean(mig_img,4);
V_mig = permute(V_mig,[3,2,1]);
h2 = slice(X,Y,Z,V_mig,90,90,40);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
set(h2(:),'EdgeColor','none'); set(gca,'ZDir','reverse');
colormap(flipud(roma)); cmax = rms(V_mig(:)); caxis([-cmax cmax]);
title('Migration - Average Image');
colorbar;

% === 3D Least-Squares Migration ===
subplot(1,3,3);
V_migls = mean(mig_imgls,4);
V_migls = permute(V_migls,[3,2,1]);
h3 = slice(X,Y,Z,V_migls,90,90,40);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
set(h3(:),'EdgeColor','none'); set(gca,'ZDir','reverse');
colormap(flipud(roma)); cmax = rms(V_migls(:)); caxis([-cmax cmax]);
title('Least-Squares Migration - Average Image');
colorbar;

%% 9. Save comparison results
% Store final comparison results for further analysis
results.CCP = ccp_final;
results.Migration = mig_final;
results.LSMigration = migls_final;
results.X = X;
results.Y = Y;
results.Z = Z;
results.Fold = fold_map;
results.Parameters = struct('dx', dx, 'dy', dy, 'dz', dz, 'zmax', zmax, ...
                           'nEvents', nProcessedEvents-1);

% Save comparison results
save([config.outputFolder, '/CCP_Migration_Comparison.mat'], 'results');
fprintf('Comparison results saved to: %s/CCP_Migration_Comparison.mat\n', config.outputFolder);

%% 10. Summary statistics
% Display processing summary
fprintf('\n=== Processing Summary ===\n');
fprintf('Total events processed: %d\n', nProcessedEvents-1);
fprintf('Grid dimensions: %d x %d x %d\n', length(X), length(Y), length(Z));
fprintf('Grid spacing: dx=%.1f km, dy=%.1f km, dz=%.1f km\n', dx, dy, dz);
fprintf('Maximum depth: %.1f km\n', zmax);
fprintf('Results saved to: %s/CCP_Migration_Comparison.mat\n', config.outputFolder);
