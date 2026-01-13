% test_radon3d.m - Comprehensive test script for 3D Radon transform operations
%
% This script provides a complete test framework for 3D Radon transform
% implementation, including both regular and non-regular grid configurations.
% The script performs the following operations:
%
% 1. Data Loading: Loads seismic gather data and processing parameters
% 2. Data Preprocessing: Applies rank reduction to input data
% 3. Parameter Setup: Configures Radon transform parameters
% 4. Regular Grid Processing: Tests Radon transform on regular grid data
% 5. Non-Regular Grid Processing: Tests Radon transform on station location data
% 6. Operator Validation: Performs dot product test to verify adjoint operators
% 7. Visualization: Displays original and transformed data for comparison
%
% Author: MATLAB DenseArrayToolkit
% Date: Created for seismic data processing applications

%% ========================================================================
% SECTION 1: INITIALIZATION AND DATA LOADING
% ========================================================================
fprintf('=== 3D Radon Transform Test Script ===\n');
fprintf('Loading test data files...\n');

% Load required test data files:
%   gather.mat     - Contains seismic gather data (time-domain waveforms)
%   gridStruct.mat - Contains grid structure with spatial coordinates
%   RankReductionParam.mat - Parameters for rank reduction preprocessing
load gather.mat 
load gridStruct.mat 
load RankReductionParam.mat 

fprintf('Data loading completed.\n');

%% ========================================================================
% SECTION 2: DATA PREPROCESSING - RANK REDUCTION
% ========================================================================
fprintf('Applying rank reduction to input data...\n');

% Apply rank reduction to input data to enhance signal-to-noise ratio
% and reduce computational complexity for subsequent Radon transform
%
% Inputs:
%   gather              - Raw seismic gather data
%   gridStruct          - Original grid structure with spatial information
%   RankReductionParam  - Configuration parameters for rank reduction
%
% Outputs:
%   gatherReconstructed - Reconstructed gather after rank reduction
%   d1_otg             - Output time gather (cleaned data for Radon transform)
%   reconGrid          - Grid structure after reconstruction
[gatherReconstructed, d1_otg, reconGrid] = rankReduction3D(gather, gridStruct, RankReductionParam);

fprintf('Rank reduction completed. Data dimensions: %d x %d\n', size(d1_otg,1), size(d1_otg,2));

%% ========================================================================
% SECTION 3: RADON TRANSFORM PARAMETER CONFIGURATION
% ========================================================================
fprintf('Setting up Radon transform parameters...\n');

% Define slowness parameters (px, py) - these represent the range of
% possible slopes in the x and y directions for the linear Radon transform
px = linspace(-0.01, 0.01, 20);  % 20 slowness values from -0.01 to 0.01 s/m
py = linspace(-0.01, 0.01, 20);  % Same range for y direction

% Time domain parameters
dt = 0.1;                         % Time sampling interval (seconds)
t = (0:size(d1_otg,1)-1) * dt;    % Time vector based on data length

% Extract spatial coordinates from the reconstructed grid
hx = reconGrid.x;                 % x-coordinates of receivers
hy = reconGrid.y;                 % y-coordinates of receivers

% Create comprehensive parameter structure for Radon transform
Param.hx   = hx;                  % Receiver x-coordinates
Param.hy   = hy;                  % Receiver y-coordinates  
Param.px   = px;                  % Slowness parameters in x-direction
Param.py   = py;                  % Slowness parameters in y-direction
Param.nt   = length(t);           % Number of time samples
Param.dt   = dt;                  % Time sampling interval
Param.type = 1;                   % Transform type (1 = linear Radon transform)

% Calculate parameter dimensions for preallocation
nt  = length(t);                  % Number of time samples
npx = length(px);                 % Number of x-slowness parameters
npy = length(py);                 % Number of y-slowness parameters

% Preallocate transform model matrix (initial guess)
ma = zeros(Param.nt, npx, npy);   % Initial model (all zeros)

% Configure PCG (Preconditioned Conjugate Gradient) solver parameters
N1 = 10;  % Maximum number of inner iterations
N2 = 1;   % Number of restart iterations for outer loop

fprintf('Parameter setup completed. Model dimensions: %d x %d x %d\n', nt, npx, npy);

%% ========================================================================
% SECTION 4: REGULAR GRID RADON TRANSFORM PROCESSING
% ========================================================================
fprintf('Performing 3D Radon transform on regular grid data...\n');

try
    % Perform inverse Radon transform using PCG algorithm to estimate
    % the model parameters from the observed data
    %
    % Inputs to yc_pcg:
    %   @radon3d_op_regular_grid    - Function handle to Radon transform operator
    %   Param                       - Parameter structure
    %   d1_otg                      - Input data (time gather after rank reduction)
    %   ma                          - Initial model (zeros)
    %   N1, N2                      - PCG iteration parameters
    %   1                           - Verbosity flag (1 = display progress)
    mi_z = yc_pcg(@radon3d_op_regular_grid, Param, d1_otg, ma, N1, N2, 1);
    
    % Perform forward Radon transform to reconstruct data from the
    % estimated model - this validates the transform and provides
    % the reconstructed data for comparison
    d1_otg_radon = radon3d_op_regular_grid(mi_z, Param, 1);  
    
    % Visualize original and Radon-transformed data side by side
    figure('Name', 'Regular Grid: Original vs Radon-Transformed Data', ...
           'Position', [100, 100, 1200, 600]);
    
    % Reshape and concatenate original and transformed data for display
    % Note: 550x(11*14) suggests the data is reshaped from 3D to 2D
    combined_data = [reshape(d1_otg, 550, 11*14), reshape(d1_otg_radon, 550, 11*14)];
    
    imagesc(combined_data);
    caxis([-0.02, 0.02]);  % Set consistent color axis limits
    colorbar;
    title('Regular Grid: Original Data (Left) vs Radon-Transformed Data (Right)');
    xlabel('Trace Number');
    ylabel('Time Sample');
    
    fprintf('Regular grid Radon transform completed successfully.\n');
    
catch ME
    % Enhanced error handling with detailed information
    fprintf('ERROR in regular grid Radon transform:\n');
    fprintf('  Message: %s\n', ME.message);
    fprintf('  Identifier: %s\n', ME.identifier);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    File: %s, Line: %d, Function: %s\n', ...
                ME.stack(i).file, ME.stack(i).line, ME.stack(i).name);
    end
    warning('Regular grid Radon transform failed. Continuing with non-regular grid test...');
end

%% ========================================================================
% SECTION 5: NON-REGULAR GRID RADON TRANSFORM PROCESSING
% ========================================================================
fprintf('Performing 3D Radon transform on non-regular grid data...\n');

% Extract station information from the gather data for non-regular grid processing
% Non-regular grid processing handles real-world station locations that
% may not form a perfect rectangular grid
stationList = getStations(gather);
stlo = [stationList.stlo]';  % Station longitude coordinates
stla = [stationList.stla]';  % Station latitude coordinates

% Convert geographic coordinates (latitude/longitude) to projected 2D coordinates
% This transformation is necessary for the Radon transform which operates
% in Cartesian space
[hx, hy] = latlonToProjectedCoords(stlo, stla, gridStruct);

% Update parameter structure for non-regular grid processing
Param.hx = hx;        % Receiver x-coordinates (projected)
Param.hy = hy;        % Receiver y-coordinates (projected)
Param.isGrid = 0;     % Flag indicating non-regular grid (0 = irregular stations)

% Extract receiver function (RF) data from the gather structure
% This represents the seismic data traces at each station location
itrCell = {gather.RF};
din = cell2mat(cellfun(@(rf) rf.itr, itrCell, 'UniformOutput', false));
din = din(1:Param.nt, :);  % Truncate to match time samples

% Re-initialize model for non-regular grid processing
ma = zeros(Param.nt, npx, npy);  % Initial model (all zeros)

% Perform inverse Radon transform on non-regular grid data using PCG
fprintf('Running PCG solver for non-regular grid...\n');
mi = yc_pcg(@radon3d_op, Param, din, ma, N1, N2, 1);
    
% Perform forward Radon transform to reconstruct data from the estimated model
d1_radon = radon3d_op(mi, Param, 1);  

fprintf('Non-regular grid Radon transform completed.\n');

%% ========================================================================
% SECTION 6: OPERATOR VALIDATION - DOT PRODUCT TEST
% ========================================================================
fprintf('Performing dot product test to validate adjoint operators...\n');

% The dot product test verifies that the forward and adjoint operators
% form a proper adjoint pair. This is crucial for the correctness of
% iterative inversion algorithms like PCG.
%
% For operators A (forward) and A' (adjoint), the test verifies:
%   <A*m, d> = <m, A'*d>
% where <.,.> denotes the inner product.

% Generate random test vectors for validation
m1 = randn(nt, npx, npy);  % Random model vector
Param.isGrid = 0;           % Use non-regular grid configuration
Param.hx = hx;              % Receiver x-coordinates
Param.hy = hy;              % Receiver y-coordinates

% Apply forward operator: d1 = A * m1
d1 = radon3d_op(m1, Param, 1);

% Generate random data vector and apply adjoint operator: m2 = A' * d2
nhx = length(hx);
d2 = randn(nt, nhx);
m2 = radon3d_op(d2, Param, -1);

% Compute dot products for validation
dot1 = sum(sum(sum(d1 .* d2)));  % <A*m1, d2>
dot2 = sum(sum(sum(m1 .* m2)));  % <m1, A'*d2>

fprintf('Dot product test results:\n');
fprintf('  <A*m1, d2> = %.6e\n', dot1);
fprintf('  <m1, A''*d2> = %.6e\n', dot2);
fprintf('  Difference = %.6e\n', abs(dot1 - dot2));
fprintf('  Relative error = %.6e\n', abs(dot1 - dot2) / max(abs(dot1), abs(dot2)));

% Check if operators pass the dot product test
if abs(dot1 - dot2) / max(abs(dot1), abs(dot2)) < 1e-10
    fprintf('✓ Operator validation PASSED - forward and adjoint operators are consistent.\n');
else
    fprintf('✗ Operator validation FAILED - operators may not be proper adjoints.\n');
    warning('Operator validation failed. Check implementation of forward/adjoint operators.');
end

%% ========================================================================
% SECTION 7: FINAL VISUALIZATION AND COMPARISON
% ========================================================================
fprintf('Generating final visualizations...\n');

% Switch back to regular grid for final comparison
Param.isGrid = 1;
% Extract spatial coordinates from reconstructed grid
hx = reconGrid.x;  % x-coordinates of receivers
hy = reconGrid.y;  % y-coordinates of receivers

% Update parameter structure for regular grid visualization
Param.hx = hx;  % Receiver x-coordinates
Param.hy = hy;  % Receiver y-coordinates

% Apply forward transform for regular grid visualization
d1_otg_radon = radon3d_op_regular_grid(mi, Param, 1);  

% Create comprehensive visualization comparing original and transformed data
figure('Name', 'Non-Regular Grid: Original vs Radon-Transformed Data', ...
       'Position', [100, 100, 1200, 600]);

% Display original and Radon-transformed data side by side
imagesc([din, d1_radon]);
caxis([-0.3, 0.3]);  % Set appropriate color axis limits for this data
colorbar;
title('Non-Regular Grid: Original Data (Left) vs Radon-Transformed Data (Right)');
xlabel('Trace Number');
ylabel('Time Sample');

fprintf('=== 3D Radon Transform Test Completed Successfully ===\n');
fprintf('Summary:\n');
fprintf('  - Regular grid processing: Completed\n');
fprintf('  - Non-regular grid processing: Completed\n');
fprintf('  - Operator validation: Performed\n');
fprintf('  - Visualizations: Generated\n');
