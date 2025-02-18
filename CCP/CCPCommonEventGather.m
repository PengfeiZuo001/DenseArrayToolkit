function ccpResult = CCPCommonEventGather(gather, velocityModel, profileStruct, param)
% CCPCommonEventGather  Perform Common Conversion Point stacking (CCP) on seismic data.
%
% Usage:
%   ccpResult = CCPCommonEventGather(gather, velocityModel, profileStruct, param)
%
% Inputs:
%   gather         : Struct array of seismic traces for one event or gather.
%                    Each element includes fields:
%                       .RF.itr            - Iterative decon RF [Nt x 1]
%                       .TravelInfo.rayParam, .TravelInfo.baz - Ray parameters
%                       .TimeAxis.t_resample, .TimeAxis.dt_resample - time axis
%   velocityModel  : Struct describing velocity model with fields:
%                       .x, .z, .nx, .nz, .dx, .dz, .vp, .vs, ...
%   profileStruct  : Struct with profile-related info:
%                       .line_points   [N x 2] (lon, lat) or (x, y)
%                       .center        [lon0, lat0] or [x0, y0]
%                       .direction     'NE-SW' or other orientation
%   param          : Parameters for migration and reconstruction:
%       .dz         (double) - Depth sampling for migration (km)
%       .zmax       (double) - Maximum depth for imaging (km)
%       .plotBinned (bool)   - If true, show binned data wiggle plot
%       .xpad       (double) - Horizontal padding
%       etc.
%
% Outputs:
%   ccpResult : Struct containing CCP results with fields:
%       .x, .z    - Horizontal and depth axes
%       .img      - CCP image
%       .param    - A copy of the param struct for reference
%
% Author: Yunfeng Chen
% Date: Feb. 18, 2025

%% Input Validation and Parameter Defaults
if nargin < 4
    param = struct();
end

% Validate gather
if isempty(gather) || ~isstruct(gather)
    error('CCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end

% Validate velocityModel fields
requiredVMfields = {'x', 'z', 'vp', 'vs', 'dx', 'dz', 'nx', 'nz'};
for f = requiredVMfields
    if ~isfield(velocityModel, f{1})
        error('CCPCommonEventGather:MissingVelocityModelField', 'Field "%s" is missing in velocityModel.', f{1});
    end
end

% Fill default parameters
defaultParams = struct('dz', 1, 'zmax', 200, 'xpad', 0);
param = fillDefaultParams(param, defaultParams);

% Display event info if available
if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    disp(['Processing event: ', gather(1).EventInfo.evid]);
end

%% Unpack Profile and Velocity Model Information
% Extract line endpoints from profileStruct
lon1 = profileStruct.line_points(1, 1);
lat1 = profileStruct.line_points(1, 2);
lon2 = profileStruct.line_points(end, 1);
lat2 = profileStruct.line_points(end, 2);

% Unpack velocity model fields
vp   = velocityModel.vp(:, 1);  % Assuming vp is 1D (velocity at top layer)
vs   = velocityModel.vs(:, 1);  % Assuming vs is 1D (velocity at top layer)
xpad = param.xpad;
zmax = param.zmax;

%% Calculate Conversion Points
nrf = length(gather);

% Extract RFs and times from gather
rfsAll = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
timeAll = cellfun(@(rf) rf.ittime, {gather.RF}, 'UniformOutput', false);

% Extract ray parameters and back-azimuths from gather
raypAll = cellfun(@(ti) ti.rayParam / 6371, {gather.TravelInfo}, 'UniformOutput', false);
raypAll = cell2mat(raypAll);
bazAll = cellfun(@(ti) ti.baz, {gather.TravelInfo}, 'UniformOutput', false);
bazAll = cell2mat(bazAll);

% Extract station coordinates from gather
latAll = cellfun(@(si) si.stla, {gather.StationInfo}, 'UniformOutput', false);
latAll = cell2mat(latAll);
lonAll = cellfun(@(si) si.stlo, {gather.StationInfo}, 'UniformOutput', false);
lonAll = cell2mat(lonAll);

% Perform ray tracing for CCP points
model_type = 'flat';
disp('Ray tracing started');
tic;
[cp, ~, ~] = rf_ccp(raypAll, bazAll, velocityModel.dz, zmax, velocityModel.z, vp, vs, latAll, lonAll, model_type);
toc;
disp('Ray tracing completed');

% Perform time-to-depth conversion for RFs
disp('Time-to-depth conversion started');
tic;
[~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, velocityModel.dz, param.zmax, velocityModel.z, vp, vs);
toc;
disp('Time-to-depth conversion completed');

% Attach RF amplitudes to CCP points
for k = 1:nrf
    cp(k).amp = rfsAll_depth{k};
end

% Project CCP points to profile
ndep = length(cp(1).zpos);
dist_projected = zeros(ndep, nrf);

% Calculate projected distances for each RF
for k = 1:nrf
    slon_tmp = [cp(k).lonb];
    slat_tmp = [cp(k).latb];
    data = [slon_tmp(:), slat_tmp(:)];
    centered_data = data - profileStruct.center;
    projected_points = profileStruct.center + (centered_data * profileStruct.direction) * profileStruct.direction';
    slat_projected = projected_points(:, 2);
    slon_projected = projected_points(:, 1);
    [deg, ~] = distance(lat1, lon1, slat_projected, slon_projected);
    dist = deg * 2 * pi * 6371 / 360;
    dist_projected(:, k) = dist;
end

%% CCP Stacking Process
% Set grid sizes (distance direction: dx, depth direction: dz)
[X, Z] = meshgrid(velocityModel.x, velocityModel.z);
nx = length(velocityModel.x);
nz = length(velocityModel.z);

% Initialize stacking result matrices
V = zeros(nz, nx); % Store accumulated amplitude values
count = zeros(nz, nx); % Store sample count for each grid

depth = [cp.zpos];
amp = [cp.amp];
amp = amp(2:end,:);
[nt,ntr] = size(amp);

% Perform parallel processing for speed
parfor i = 1:nz
    for j = 1:nx
        for n = 1:ntr
            xi = dist_projected(:, n);
            zi = depth(:, n);
            vi = amp(:, n);
            keep = xi >= velocityModel.x(j) - 2 * velocityModel.dx & xi <= velocityModel.x(j) + 2 * velocityModel.dx & ...
                   zi >= velocityModel.z(i) - 2 * velocityModel.dz & zi <= velocityModel.z(i) + 2 * velocityModel.dz;
            V(i, j) = V(i, j) + sum(vi(keep));
            count(i, j) = count(i, j) + sum(keep);
        end
    end
end

%% Generate CCP Image
% Calculate average amplitude for each grid if there are valid data in the grid
VI = V ./ max(count, 1);  % 避免除以0

%% Plot CCP Image
try
    load roma;
    cmap = flipud(roma);
catch
    cmap = parula;
end

figure;
set(gcf, 'Position', [100 100 800 400], 'color', 'w');
imagesc(velocityModel.x, velocityModel.z, VI);
caxis([-0.1 0.1]);
colormap(cmap);
colorbar;
xlabel('Distance (km)');
ylabel('Depth (km)');
title('CCP image');
set(gca, 'fontsize', 14);

%% Save CCP results
ccpResult = struct('depth', depth, 'dist_projected', dist_projected, 'amp', amp, ...
                   'x', velocityModel.x, 'z', velocityModel.z, 'img', VI);

end

function param = fillDefaultParams(param, defaultParams)
% FILLDEFAULTPARAMS  Fill default parameters in param struct if they are not already set.
%
% Inputs:
%   param         : Parameters struct to be filled.
%   defaultParams : Struct with default parameter values.
%
% Outputs:
%   param         : Parameters struct with default values filled.

% Fill missing fields in param with defaultParams
fields = fieldnames(defaultParams);
for i = 1:numel(fields)
    if ~isfield(param, fields{i})
        param.(fields{i}) = defaultParams.(fields{i});
    end
end

end
