%% DenseArrayToolkit 3D Migration Main Function Example
% This script demonstrates a complete workflow for 3D seismic imaging using
% least-squares migration methods. It processes seismic array data to create
% volumetric subsurface images in 3D space.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata for 3D processing
%   4. Create velocity model - Set up 3D velocity structure for migration
%   5. Migration imaging - Perform 3D least-squares migration
%   6. Save results - Store processed 3D imaging results
%
% This script implements 3D least-squares migration with rank reduction
% preprocessing for improved imaging quality in dense array applications.

clear; clc; close all;

%% 0. Setup paths and parameters
% Initialize the processing environment by adding necessary functions and
% dependencies to the MATLAB path. This ensures access to all required
% processing routines in the DenseArrayToolkit.
cd ../
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
% - RankReductionParam: Parameters for rank reduction preprocessing
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;
RankReductionParam = config.RankReductionParam;

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
% - 3D migration volume definition
% - Processing optimization for dense arrays
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
% Set up the 3D imaging grid and velocity model for migration:
% 1. Define grid spacing in x, y, and z directions
% 2. Create 3D imaging volume based on array geometry
% 3. Generate 3D velocity model for migration
% Note: Using 3D velocity model for accurate ray tracing in volumetric imaging
dx = 10;    % Horizontal x-direction grid spacing (km)
dy = 10;    % Horizontal y-direction grid spacing (km)
dz = 1;     % Vertical grid spacing (km)
zmax = 100; % Maximum imaging depth (km)
xpad = 20;  % Horizontal padding in x-direction (km)
ypad = 20;  % Horizontal padding in y-direction (km)

% Create 3D imaging grid with specified parameters
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 3D velocity model for migration imaging
npts = 5;
gridStruct = getVelocityModel('3D', gridStruct, npts);

%% update param of rank reduction
% note: regular grid is greater than stations (rx,ry)
% This part needs to be updated to make it more concise!!!
RankReductionParam.nx = gridStruct.nx;
RankReductionParam.ny = gridStruct.ny;
RankReductionParam.ox = min(gridStruct.x); 
RankReductionParam.oy = min(gridStruct.y); 
RankReductionParam.mx = max(gridStruct.x);
RankReductionParam.my = max(gridStruct.y);
RankReductionParam.rank = 10;

%% 5. Migration imaging
% Perform 3D least-squares migration with rank reduction preprocessing:
% - Standard 3D migration
% - 3D least-squares migration for improved resolution
% - Rank reduction preprocessing for noise suppression
%
% Initialize structure to store migration results farom all events
% migResults = [];
ccpResults = [];
nMigratedEvents = 1;    % Counter for successfully processed events

% Set up 3D migration parameters based on grid structure
MigParam.paramMig = setMigParam3D(gridStruct);

minTrace = 100; % Minimum number of traces required for CCP imaging
minSNR = 3;    % Minimum SNR of the RF required for CCP imaging

migResults = [];
ccpResults = [];
% Process each event in the dataset
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    
    % Extract seismic records for current event (Common Event Gather)
    % This groups all recordings of the same event across different stations
    gather = getCommonEventGather(DataStruct, evid);
    % Extract the SNR
    snrAll = cell2mat(cellfun(@(rf) rf.snr, {gather.RF}, 'UniformOutput', false));    
    % Quality control: Skip events with insufficient station coverage
    % Minimum minTrace stations required to ensure reliable 3D imaging results
    if length(gather) < minTrace || mean(snrAll)< minSNR
        continue
    end
    
    % Compute receiver functions through deconvolution
    % This isolates converted phases from the P-wave coda
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);

    % Apply rank reduction preprocessing to improve signal quality
    % This helps suppress noise and enhance coherent signals
    MigParam.paramMig.isReconRFs=1;
    MigParam.gauss = 5;
    MigParam.tmax = 80;
    if MigParam.paramMig.isReconRFs
        RankReductionParam.rank = 5;
        RankReductionParam.fhigh = 2.4;
        [gatherDRR, d1_otg] = rankReduction_new(gather, gridStruct, RankReductionParam);
        MigParam.paramMig.dotg = d1_otg;
        itrCell = {gather.RF};
        d0 = cell2mat(cellfun(@(rf) rf.itr, itrCell,'UniformOutput', false));
        itrCell = {gatherDRR.RF};
        d1 = cell2mat(cellfun(@(rf) rf.itr, itrCell,'UniformOutput', false));
        t = gather(1).RF.ittime;
        figure;
        set(gcf,'Position',[100 100 1200 600],'Color','w')
        ax1=subplot(121);
        imagesc(1:length(gather),t(30:300),d0(30:300,:))
        caxis([-0.1 0.1])
        colormap(seismic(1))
        xlabel('Trace Number')
        ylabel('Time (s)')
        set(ax1,'fontsize',14)
        ax2=subplot(122);
        imagesc(1:length(gather),t(30:300),d1(30:300,:))
        caxis([-0.1 0.1])
        xlabel('Trace Number')
        ylabel('Time (s)')
        set(ax2,'fontsize',14)
        export_fig(['./figures/BY_rank_',num2str(evid),'.png'], '-r300');

        figure;
        d3d = permute(d1_otg,[2,3,1]);
        [nx,ny,nt] = size(d3d);
        h=slice(d3d,round(ny/2),round(nx/2),80);
        set(h,'EdgeColor','none');
        colormap(seismic(1))
        set(gca,'ZDir','reverse')
        caxis([-0.1 0.1])
        zlim([0 300])
        export_fig(['./figures/BY_rank_3D_',num2str(evid),'.png'], '-r300');
    end
    % Perform 3D least-squares migration
    % This method provides improved resolution compared to standard migration
    migResult = leastSquaresMig3D(gatherDRR, gridStruct, MigParam);

    % Store migration results for current event
    % mig - Standard 3D migration results
    % migls - 3D least-squares migration results
    migResults = [migResults; migResult];
 

    ccpResult = CCPCommonEventGather(gatherDRR, gridStruct, CCPParam);
    ccpResults = [ccpResults; ccpResult];

    % Pause for visualization (can be removed for batch processing)
%     pause;
    
    % Close all figure windows to avoid memory issues during batch processing
    close all;
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. Stacking
% ccpImage = stackImagingResults(ccpResults);
% migImage = stackImagingResults(migResults);

%% 7. Visualization
V = zeros(size(migResults(1).mig));
for n=1:length(migResults)
    V = V+migResults(n).migls;
end
V = V/length(migResults);
V = permute(V,[3,2,1]);
migResult.V = V;

count = 0;
for n=1:length(ccpResults)
    V = V+ccpResults(n).img;
    count = count+ccpResults(n).count;
end
ccpResult.V = V./max(count,1);

% Configure visualization options
options = struct();
options.profileType = 'predefined';  % Enable interactive profile selection
% options.smoothingParams = struct(...
%     'radius', 3, ...        % Smoothing radius
%     'eps', 0.01, ...       % Regularization parameter
%     'order', 2);           % Smoothing order
% N-S profile crossing Baiyan Obo minning area
options.profilePoints(:,1) = [76.6927 76.6927 nan 0 200];
options.profilePoints(:,2) = [0 150 nan 66.9016 66.9016];
% E-W profile crossing Baiyan Obo minning area
% options.profilePoints(:,1) = [0; 200];
% options.profilePoints(:,2) = [66.9016; 66.9016];

visualizeImage(migResult, gridStruct, options);
visualizeImage(ccpResult, gridStruct, options);
% pointA = [76.6927, 66.9016, 0]; % 点A (x1,y1,z1)
% pointB = [76.6927, 66.9016, 100]; % 点B (x2,y2,z2)
% lx = [pointA(1), pointB(1)];
% ly = [pointA(2), pointB(2)];
% lz = [pointA(3), pointB(3)];
% hold on;
% plot3(lx, ly, lz, '--','LineWidth', 2,'Color', [0.5 0.5 0.5]);