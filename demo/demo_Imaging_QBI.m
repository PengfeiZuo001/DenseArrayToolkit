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
% eventid = {eventList.evid};

% load good events
fid = fopen('good_evt_new_id.txt','r');
C = textscan(fid,'%s');
eventid = C{1};

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
dz = 1;
zmax = 100;
xpad = 30;
ypad = 30;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);
% Create or obtain 3D velocity model with specified sampling (10 points)
% This model will be used for ray tracing in migration
npts = 5;
gridStruct = getVelocityModel('3D',gridStruct,npts);

% Set up 3D migration parameters based on grid structure
MigParam.paramMig = setMigParam3D(gridStruct);
%% 5. Migration imaging
% Perform Common Conversion Point (CCP) stacking for each event that meets quality
% criteria. This process includes:
% 1. Extract event gathers (recordings of the same event at different stations)
% 2. Quality control: Skip events with insufficient station coverage (<60 stations)
% 3. Optional Radon Transform for array processing
% 4. Receiver function calculation through deconvolution
% 5. CCP stacking to create 3D image volume

% Initialize arrays for accumulating migration results across all events
nMigratedEvents = 1;   % Counter for successfully processed events
minTrace = 60; % Minimum number of traces required for CCP imaging
minSNR = 5;    % Minimum SNR of the RF required for CCP imaging
DeconvParam.radonfilter = 1; % Apply radon transform to seismograms
DeconvParam.gauss = 2.5;
CCPParam.smoothLength = 3;
% Process each event that meets the filtering criteria
ccpResults = [];
ccpResultsRadon = [];
migResults = [];

for iEvent = 1:length(eventid)
    evid = eventid{iEvent};
%     evid='20231221145556';
    % Extract seismic records for current event (Common Event Gather)
    gather = getCommonEventGather(DataStruct, evid);
    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));
    % Skip events with fewer than 60 valid stations to ensure imaging quality
%     if length(gather) < minTrace || mean(snrAll)< minSNR || length(gather) > 100
%         continue
%     end
    DeconvParam.radonfilter = false;
    gather = deconv(gather, DeconvParam);
    % Check if Radon filtering is enabled for enhanced signal quality
    DeconvParam.radonfilter = true;
    if DeconvParam.radonfilter
        RadonParam.highs = 1.2;
        RadonParam.pmax = 0.02;
        RadonParam.pmin = -0.02;
        % Apply Radon Transform for array processing
        gatherRadon = radonTransform(gather, gridStruct, RadonParam);
%         gatherRadon = deconv(gather, DeconvParam);
%         export_fig(['./figures/radon_post_rfs_',num2str(evid),'.png'],'-r300')
    end    
    
    % compare two RFs
%     figure;
%     set(gcf,'Position',[100 100 400 800],'Color','w')
%     ax1 = subplot(211);
%     plotCommonEventGather(gather, [], 'trace', 'imagesc',ax1)
%     set(ax1,'YLim',[0 20])
%     ax2 = subplot(212);
%     plotCommonEventGather(gatherRadon, [], 'trace', 'imagesc',ax2)
%     set(ax2,'YLim',[0 20])
%     export_fig(['./figures/radon_postdecon_rfs_',num2str(evid),'.png'],'-r300')

    % Apply Common Conversion Point stacking for the current event gather
    CCPParam.smoothLength = 3;
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    CCPParam.smoothLength = 3;
    ccpResultRadon = CCPCommonEventGather(gatherRadon, gridStruct, CCPParam);
    ccpResultsRadon = [ccpResultsRadon; ccpResultRadon];

    % Apply migration imaging
    MigParam.itermax = 10;
    MigParam.tmax = 80;
    migResult = leastSquaresMig3D(gatherRadon, gridStruct, MigParam);
    migResults = [migResults; migResult];

%      pause;
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
count = 0;
V = zeros(size(ccpResults(1).img));
for n=1:length(ccpResults)
    V = V+ccpResults(n).img;
    count = count+ccpResults(n).count;
end
ccpResult.V = V./max(count,1);

count = 0;
V = zeros(size(ccpResults(1).img));
for n=1:length(ccpResultsRadon)
    V = V+ccpResultsRadon(n).img;
    count = count+ccpResultsRadon(n).count;
end
ccpResultRadon.V = V./max(count,1);

V = zeros(size(migResults(1).migls));
for n=1:length(migResults)
    V = V+migResults(n).migls;
end
V = V/length(migResults);
V = permute(V,[3,2,1]);
migResult.V = V-mean(V(:));

% ccpResult.V = sum(dimgRadon,4)./max(sum(countRadon,4;),1);
% Configure visualization options
options = struct();
options.profileType = 'predefined';  % Enable interactive profile selection
options.smoothingParams = struct(...
    'radius', 3, ...        % Smoothing radius
    'eps', 0.01, ...       % Regularization parameter
    'order', 2);           % Smoothing order
% define profile location
% profilePoints = [
% 91.461	37.47
% 91.948	38.168
% 93.599	38.782
% 93.994	38.681
% 95.53531 39.21347];
% [px,py] = latlonToProjectedCoords(profilePoints(:,1), profilePoints(:,2), gridStruct);
profilePoints = [
    -20 -10;
    107 48;
    212 47;
    345 -5];
px = profilePoints(:,1);
py = profilePoints(:,2);
options.profilePoints(:,1) = px;
options.profilePoints(:,2) = py;
% Generate visualizations and get profile data structure
% profileStruct = visualizeCCPResults(migResult, gridStruct, options);
visualizeCCPResults(ccpResult, gridStruct, options);
visualizeCCPResults(ccpResultRadon, gridStruct, options);
visualizeCCPResults(migResult, gridStruct, options);
%% 7. Save results
% Save the final CCP imaging results to a MAT file for future reference and analysis.
% The saved results include:
% - 3D image volume
% - Grid coordinates
% - Processing parameters
% - Other relevant metadata
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);
