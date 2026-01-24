%% DenseArrayToolkit 3D CCP Imaging Example - BaiyanEbo Region
% This script demonstrates 3D Common Conversion Point (CCP) stacking for
% seismic imaging in the BaiyanEbo region. It processes receiver functions
% to create volumetric subsurface images showing crustal and upper mantle
% structures beneath the dense seismic array.
%
% Main processing steps include:
%   0. Setup paths and parameters - Initialize environment and load configurations
%   1. Read data - Import seismic waveform data in SAC format
%   2. Preprocessing - Apply filters and prepare data for analysis
%   3. Get array and event information - Extract metadata for 3D processing
%   4. Create velocity model - Set up 3D velocity structure for CCP stacking
%   5. Compute receiver functions - Extract P-to-S converted phases
%   6. CCP stacking - Perform 3D Common Conversion Point imaging
%   7. Results output - Visualize and save 3D imaging results
%
% This script implements 3D CCP stacking with rank reduction preprocessing
% for improved signal-to-noise ratio in dense array applications.

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
% - MigParam: Parameters controlling migration algorithms (not used in CCP)
% - RadonParam: Settings for Radon transform array processing (optional)
% - DeconvParam: Parameters for receiver function deconvolution
% - CCPParam: Settings for Common Conversion Point stacking
% - RankReductionParam: Parameters for rank reduction preprocessing
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;

%% 1. Read data
% Load seismic waveform data in SAC format from the specified directory
% The data is encapsulated into a structured array (DataStruct) containing
% waveforms and metadata for each recording
dataFolder1= './data/event_waveforms_BY';
DataStruct1 = read_SAC(dataFolder1);

dataFolder2= '/Users/yunfeng/30_40/research/BY_local/event_waveforms';
DataStruct2 = read_SAC(dataFolder2);

DataStruct = [DataStruct1 DataStruct2];
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
% - 3D CCP volume definition
% - Processing optimization for dense arrays
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% Extract station and event coordinates for spatial analysis
stlo = [stationList.stlo]';  % Station longitude
stla = [stationList.stla]';  % Station latitude
evla = [eventList.evla]';    % Event epicenter latitude
evlo = [eventList.evlo]';    % Event epicenter longitude

% Optional azimuthal filtering (commented out for this region)
% idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% Create filtered list of event IDs for processing
eventid = {eventList.evid};
% eventid = eventid(idxConsistentEQ);

% Generate event-station correspondence table for efficient data access
EventStationTable = getEventStationTable(DataStruct);

%% 4. Create velocity model
% Set up the 3D imaging grid and velocity model for CCP stacking:
% 1. Define grid spacing in x, y, and z directions
% 2. Create 3D imaging volume based on array geometry
% 3. Generate 3D velocity model for accurate depth conversion
% Note: Using 3D velocity model for precise ray tracing and depth conversion
dx = 10;    % Horizontal x-direction grid spacing (km)
dy = 10;    % Horizontal y-direction grid spacing (km)
dz = 0.5;   % Vertical grid spacing (km) - finer resolution for CCP
zmax = 100; % Maximum imaging depth (km)
xpad = 40;
ypad = 40;
% Create 3D imaging grid with specified parameters
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad);

% Create 3D velocity model for CCP imaging
% npts = 10 specifies the number of interpolation points for velocity model
npts = 5;
gridStruct = getVelocityModel('3D', gridStruct, npts);

%% 5. Compute receiver functions
% Apply deconvolution to extract receiver functions from seismic waveforms
% This isolates P-to-S converted phases that contain information about
% subsurface discontinuities beneath the array
DeconvParam.verbose = 0;
DataStruct = deconv(DataStruct, DeconvParam);

%% 6. CCP stacking
% Perform 3D Common Conversion Point stacking to create volumetric images:
% - Stack receiver functions at their theoretical conversion points
% - Account for 3D ray geometry and velocity structure
% - Apply rank reduction preprocessing for noise suppression
%
% Initialize arrays to store CCP results from all events:
ccpResults = [];    % 3D CCP image results
nMigratedEvents = 0;    % Counter for successfully processed events

minTrace = 100; % Minimum number of traces required for CCP imaging
minSNR = 5;    % Minimum SNR of the RF required for CCP imaging

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

    % Apply deconvolution to extract receiver functions
    % This isolates P-to-S converted phases from the P-wave coda
    gather = deconv(gather, DeconvParam);
%     RadonParam.N1 = 5;
%     gatherRadon = radonTransform3D(gather, gridStruct, RadonParam);
    % Set up Radon transform parameters
    itrCell = {gather.RF};
    d1 = cell2mat(cellfun(@(rf) rf.itr, itrCell,'UniformOutput', false));
    % px/py - Slowness parameters (1/velocity) in x/y directions
    px = linspace(-0.01,0.01,20); % 20 slowness values from -0.01 to 0.01 s/m
    py = linspace(-0.01,0.01,20); % Same range for y direction
    dt = 0.1; % Time sampling interval (s)
    t = (0:size(d1,1)-1)*dt; % Time vector
    
    stationList = getStations(gather);
    stlo = [stationList.stlo]';  % station longitude
    stla = [stationList.stla]';  % station latitude
    
    % Convert to projected coordinates (2D)
    [hx, hy] = latlonToProjectedCoords(stlo, stla, gridStruct);    
    % Create parameter structure for Radon transform
    Param.hx  = hx;       % Receiver x-coordinates (required by radon_op)
    Param.hy  = hy;       % Receiver y-coordinates
    Param.px  = px;       % Slowness parameters in x-direction
    Param.py  = py;       % Slowness parameters in y-direction
    Param.nt = length(t); % Number of time samples
    Param.dt = dt;        % Time sampling interval
    Param.type = 1;       % Transform type (1 = linear Radon transform)
    Param.isGrid = 0;

    % Get number of slowness parameters
    npx = length(px);
    npy = length(py);
    
    % Preallocate transform model matrix
    ma = zeros(Param.nt, npx, npy); % Initial model (all zeros)
    
    % Set PCG (Preconditioned Conjugate Gradient) parameters
    N1 = 10; % Maximum number of iterations
    N2 = 1;  % Number of restart iterations
    
    % Apply 3D Radon Transform using PCG algorithm

    % Perform inverse Radon transform using PCG
    % yc_pcg - Preconditioned Conjugate Gradient solver
    % @radon3d_op - Handle to Radon transform operator
    % d1_otg - Input data (time gather after rank reduction)
    % ma - Initial model (zeros)
    % N1, N2 - PCG parameters
    % 1 - Flag indicating forward/inverse transform
    mi = yc_pcg(@radon3d_op, Param, d1, ma, N1, N2, 1);
    
    % Perform forward Radon transform to reconstruct data from model
%     d1_radon = radon3d_op(mi_z, Param, 1);
    Param.hx = gridStruct.x;      % x-offset coordinates of receivers
    Param.hy = gridStruct.y;      % y-offset coordinates of receivers
    Param.isGrid = 1;
    d1_reg = radon3d_op(mi, Param, 1);  % forward modeling from the found model
    % Visualize original and Radon-transformed data side by side
    figure;
    % Reshape and concatenate original and transformed data for display
    imagesc([d1 d1_radon]);
    caxis([-0.1 0.1]); % Set color axis limits
    colormap(seismic(1))
    colorbar;
    title('Original (left) vs Radon-transformed (right) data');
    xlabel('Trace number');
    ylabel('Time sample');
    % Apply rank reduction preprocessing to improve signal quality
    % This helps suppress noise and enhance coherent signals
%     RankReductionParam = config.RankReductionParam;
%     RankReductionParam.rank = 3;  % Set rank reduction parameter
%     [gatherReconstructed, d1_otg,~] = rankReduction3D(gather, gridStruct, RankReductionParam);
%     RankReductionParam.rank = 10;  % Set rank reduction parameter
%     [gatherReconstructed1, d1_otg,~] = rankReduction3D(gather, gridStruct, RankReductionParam);

%     figure;
%     set(gcf,'Position',[0 0 1500 400],'Color','w')
%     ax1 = subplot(131);
%     plotCommonEventGather(gather,evid,'trace','imagesc',ax1)
%     ax1.Title.String='Raw RFs';
%     ax2 = subplot(132);
%     plotCommonEventGather(gatherReconstructed1,evid,'trace','imagesc',ax2)
%     ax2.Title.String='Rank = 10';
%     ax3 = subplot(133);
%     plotCommonEventGather(gatherReconstructed,evid,'trace','imagesc',ax3)
%     ax3.Title.String='Rank = 3'; 

%     plotZRandRF(gather)
%     [gatherReconstructed, d1_otg] = radonTransform3D(gather, gridStruct, RankReductionParam);


    % Perform 3D CCP stacking using Common Conversion Point method
    % This maps receiver functions to their theoretical conversion points
    ccpResult = CCPCommonEventGather(gatherReconstructed, gridStruct, CCPParam);
    
    % Store CCP results for current event
    ccpResults = [ccpResults; ccpResult];
    
    % Close all figure windows to avoid memory issues during batch processing
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end

%% 7. Results output
% Create final 3D CCP image by stacking all events and normalizing by hit count
% This produces the final volumetric image showing subsurface structure
ccpImage = stackImagingResults(ccpResults,7);

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
options.dem = load('./visualization/Baiyanebo_DEM_small.mat');
% E-W profile crossing Baiyan Obo minning area
% options.profilePoints(:,1) = [0; 200];
% options.profilePoints(:,2) = [66.9016; 66.9016];
visualizeImage(ccpImage, gridStruct, options);

% Optional: Save results to file for future analysis
% save './results/BaiyanEbo_ccp.mat' 'X' 'Y' 'Z' 'V' 'gridStruct'
