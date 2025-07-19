% test_radon3d.m - Test script for 3D Radon transform operations
%
% This script tests the 3D Radon transform implementation by:
% 1. Loading test data (seismic gathers and parameters)
% 2. Applying rank reduction to the input data
% 3. Setting up Radon transform parameters
% 4. Performing forward and inverse Radon transforms
% 5. Visualizing the results

% Load test data files:
%   gather.mat - Contains seismic gather data
%   gridStruct.mat - Contains grid structure information
%   RankReductionParam.mat - Contains parameters for rank reduction
load gather.mat 
load gridStruct.mat 
load RankReductionParam.mat 

% Apply rank reduction to input data
% Inputs:
%   gather - Input seismic data
%   gridStruct - Grid structure containing spatial information
%   RankReductionParam - Parameters for rank reduction
% Outputs:
%   gatherReconstructed - Reconstructed gather after rank reduction
%   d1_otg - Output time gather after rank reduction
%   reconGrid - Grid structure after reconstruction
[gatherReconstructed, d1_otg, reconGrid] = rankReduction(gather, gridStruct, RankReductionParam);

% Set up Radon transform parameters
% px/py - Slowness parameters (1/velocity) in x/y directions
px = linspace(-0.01,0.01,20); % 20 slowness values from -0.01 to 0.01 s/m
py = linspace(-0.01,0.01,20); % Same range for y direction
dt = 0.1; % Time sampling interval (s)
t = (0:size(d1_otg,1)-1)*dt; % Time vector

% Extract spatial coordinates from reconstructed grid
hx = reconGrid.x; % x-coordinates of receivers
hy = reconGrid.y; % y-coordinates of receivers

% Create parameter structure for Radon transform
Param.hx  = hx;       % Receiver x-coordinates (required by radon_op)
Param.hy  = hy;       % Receiver y-coordinates
Param.px  = px;       % Slowness parameters in x-direction
Param.py  = py;       % Slowness parameters in y-direction
Param.nt = length(t); % Number of time samples
Param.dt = dt;        % Time sampling interval
Param.type = 1;       % Transform type (1 = linear Radon transform)

% Get number of slowness parameters
npx = length(px);
npy = length(py);

% Preallocate transform model matrix
ma = zeros(Param.nt, npx, npy); % Initial model (all zeros)

% Set PCG (Preconditioned Conjugate Gradient) parameters
N1 = 10; % Maximum number of iterations
N2 = 1;  % Number of restart iterations

%% Apply 3D Radon Transform using PCG algorithm
try
    % Perform inverse Radon transform using PCG
    % yc_pcg - Preconditioned Conjugate Gradient solver
    % @radon3d_op - Handle to Radon transform operator
    % d1_otg - Input data (time gather after rank reduction)
    % ma - Initial model (zeros)
    % N1, N2 - PCG parameters
    % 1 - Flag indicating forward/inverse transform
    mi_z = yc_pcg(@radon3d_op, Param, d1_otg, ma, N1, N2, 1);
    
    % Perform forward Radon transform to reconstruct data from model
    d1_otg_radon = radon3d_op(mi_z, Param, 1);  
    
    % Visualize original and Radon-transformed data side by side
    figure;
    % Reshape and concatenate original and transformed data for display
    imagesc([reshape(d1_otg,550,11*14) reshape(d1_otg_radon,550,11*14)]);
    caxis([-0.02 0.02]); % Set color axis limits
    colorbar;
    title('Original (left) vs Radon-transformed (right) data');
    xlabel('Trace number');
    ylabel('Time sample');
    
catch ME
    % Error handling for Radon transform failure
    warnMsg = sprintf('[RadonTransform] Event %s, Z-component failed: %s', ...
        eventID, ME.message);
    warning(warnMsg);
    % Record error in processing history
    commonEventGather = appendHistory(commonEventGather, warnMsg);
    DataStruct(matchIndex) = commonEventGather;
end
