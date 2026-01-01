%% demo_rankReduction.m
% DEMO_RANKREDUCTION - Demonstration of Rank Reduction Method for Receiver Function Processing
%
% This script demonstrates a complete workflow for receiver function computation 
% and rank reduction using the Damped Rank Reduction with Off-the-Grid (DRR-OTG)
% method. The process includes reading seismic data, preprocessing, computing 
% receiver functions, applying rank reduction for signal enhancement, and 
% stacking the results for imaging.
%
% Key Features:
% - Reads SAC format seismic waveforms
% - Preprocesses data (filtering, resampling, alignment)
% - Computes receiver functions using iterative deconvolution
% - Applies 3D rank reduction for signal enhancement and noise suppression
% - Stacks receiver functions by station for final imaging
%
% Methodology:
% Rank reduction is a signal processing technique that exploits the low-rank
% structure of seismic data in the frequency-space domain to separate coherent
% signals from random noise. The DRR-OTG method handles irregular station
% geometries by reconstructing data on a regular grid while preserving the
% original signal characteristics.

clear; clc; close all;

%% 0. Setup Paths and Parameters
% --------------------------------------------------
% Load configuration parameters for all processing steps
config = loadConfig();

% Define data folder containing SAC waveform files
dataFolder = '../data/event_waveforms_BY';

% Extract processing parameters from configuration
PreprocessingParam = config.PreprocessingParam;  % Data preprocessing settings
DeconvParam        = config.DeconvParam;         % Receiver function computation
RankReductionParam = config.RankReductionParam;  % Rank reduction parameters

fprintf('=== Rank Reduction Demo for Receiver Function Processing ===\n');

%% 1. Read Seismic Data
% --------------------------------------------------
% Read SAC format waveform data from specified directory
% DataStruct contains organized seismic traces with station and event metadata
fprintf('\n[Step 1] Reading SAC waveform data from: %s\n', dataFolder);
DataStruct = read_SAC(dataFolder);
if isempty(DataStruct)
    error('No seismic data found in %s. Please check the data folder path.', dataFolder);
end
fprintf('   Successfully loaded %d seismic traces\n', length(DataStruct));

%% 2. Data Preprocessing
% --------------------------------------------------
% Apply standard preprocessing steps to prepare data for receiver function analysis:
% - Bandpass filtering (0.1-2.0 Hz typical for receiver functions)
% - Time window selection around P-wave arrival
% - Resampling to consistent sampling rate
% - Rotation to radial and transverse components
fprintf('\n[Step 2] Preprocessing seismic data...\n');
DataStruct = preprocessing(DataStruct, PreprocessingParam);
fprintf('   Preprocessing completed successfully\n');

%% 3. Receiver Function Computation
% --------------------------------------------------
% Compute receiver functions using iterative deconvolution:
% - Deconvolves vertical component from radial component
% - Estimates Earth's impulse response beneath each station
% - Uses water-level stabilization and Gaussian filtering
fprintf('\n[Step 3] Computing receiver functions using iterative deconvolution...\n');
DeconvParam.verbose = false;  % Suppress detailed output for cleaner demo
DataStruct = deconv(DataStruct, DeconvParam);
fprintf('   Receiver function computation completed\n');

%% 4. Rank Reduction Processing (DRR-OTG Method)
% --------------------------------------------------
% Apply rank reduction to enhance signal-to-noise ratio and reconstruct
% coherent signals from irregular station geometries
fprintf('\n[Step 4] Applying rank reduction (DRR-OTG) for signal enhancement...\n');

% Define 3D grid parameters for spatial reconstruction
dx = 10;    % Grid spacing in x-direction (km)
dy = 10;    % Grid spacing in y-direction (km) 
dz = 1;     % Grid spacing in z-direction (depth, km)
zmax = 100; % Maximum depth for imaging (km)
xpad = 50;  % Padding in x-direction for boundary handling (km)
ypad = 50;  % Padding in y-direction for boundary handling (km)

% Create regular grid structure for spatial reconstruction
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Get list of unique seismic events in the dataset
eventList = getEvents(DataStruct);
eventIDs  = {eventList.evid};
fprintf('   Processing %d seismic events...\n', length(eventIDs));

% Initialize structure for storing rank-reduced results
DataStructDRR = [];
minStationCount = 50;  % Minimum stations per event for reliable rank reduction

% Process each event individually
for iEvent = 1:length(eventIDs)
    evid = eventIDs{iEvent};
    
    % Extract all receiver functions for current event
    [gather, matchIndex] = getCommonEventGather(DataStruct, evid);
    
    % Skip events with insufficient station coverage
    if length(gather) < minStationCount
        fprintf('   Event %s: Only %d stations (minimum %d required) - Skipping\n', ...
            evid, length(gather), minStationCount);
        continue;
    end
    
    fprintf('   Processing event %s with %d stations...\n', evid, length(gather));
    
    % Apply 3D rank reduction using DRR-OTG method
    % - reconstructs coherent signals while suppressing random noise
    % - handles irregular station geometries through off-the-grid reconstruction
    [gatherReconstructed, d1_otg] = rankReduction3D(gather, gridStruct, RankReductionParam);
    
    % Store reconstructed receiver functions
    DataStructDRR = [DataStructDRR; gatherReconstructed(:)];
end

% Transpose results to match original DataStruct format
DataStructDRR = DataStructDRR';
fprintf('   Rank reduction completed for %d events\n', length(DataStructDRR));

%% 5. Stack Receiver Functions
% --------------------------------------------------
% Stack rank-reduced receiver functions by station to enhance coherent signals
% and create final images of subsurface structure
fprintf('\n[Step 5] Stacking receiver functions by station...\n');

% Stack common station gathers to improve signal-to-noise ratio
% - Aligns receiver functions from different events for each station
% - Averages aligned traces to suppress random noise
% - Produces high-quality station stacks for imaging
[seisout, depth0, mohoStruct] = stackCommonStationGather(DataStructDRR);

fprintf('   Stacking completed successfully\n');
fprintf('   Output dimensions: %d stations x %d depth points\n', ...
    size(seisout, 2), size(seisout, 1));

%% 6. Results Summary and Output
% --------------------------------------------------
fprintf('\n=== Processing Summary ===\n');
fprintf('Input data:  %d seismic traces from %d events\n', length(DataStruct), length(eventIDs));
fprintf('Output data: %d rank-reduced receiver functions\n', length(DataStructDRR));
fprintf('Final stacks: %d station stacks\n', size(seisout, 2));
fprintf('Depth range: 0 to %d km\n', max(depth0));

% Optional: Save results for further analysis
if config.saveResults
    outputFile = fullfile(config.outputFolder, 'rankReduction_results.mat');
    save(outputFile, 'DataStructDRR', 'seisout', 'depth0', 'mohoStruct', 'config');
    fprintf('Results saved to: %s\n', outputFile);
end

fprintf('\n=== Rank Reduction Demo Completed Successfully ===\n');
