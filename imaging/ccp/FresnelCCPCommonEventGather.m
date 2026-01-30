function ccpResult = FresnelCCPCommonEventGather(gather, gridStruct, param)
% FresenlCCPCommonEventGather_FresnelZone  Perform CCP stacking with depth-dependent weighting.
%
% This version accounts for the expansion of the Fresnel zone with depth and
% addresses sparse ray coverage by using a Gaussian weighting kernel.
%
% Usage:
%   ccpResult = FresnelCCPCommonEventGather(gather, gridStruct, param)
%
% Inputs:
%   gather      : Struct array of seismic traces (RFs and travel info)
%   gridStruct  : Struct containing imaging grid and velocity model
%   param       : Parameter struct:
%       .imagingType  - '2D' or '3D'
%       .gauss        - Gaussian filter parameter (used to estimate frequency)
%       .plotCCP      - Boolean to toggle plotting
%
% Outputs:
%   ccpResult   : Result struct (X, Y, Z coordinates, img, and weights)

%% 1. Input Validation and Parameter Defaults
if isempty(gather) || ~isstruct(gather)
    error('CCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end

if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    fprintf('Processing event: %s\n', gather(1).EventInfo.evid);
end

%% 2. Unpack Grid and Velocity Information
% Obtain a reference 1D velocity profile for ray tracing and Fresnel zone calculation
if strcmp(gridStruct.ModelType ,'1D')
    vp = gridStruct.vp(:, 1);
    vs = gridStruct.vs(:, 1);
    z  = gridStruct.z;
elseif strcmp(gridStruct.ModelType ,'2D')
    % Use mean profile for initial ray tracing
    vp = mean(gridStruct.vp, 'all'); 
    vs = mean(gridStruct.vs, 'all');
    z  = gridStruct.z;
elseif strcmp(gridStruct.ModelType ,'3D')
    vp   = mean(mean(gridStruct.VP,3),2);  % Average of 3D model
    vs   = mean(mean(gridStruct.VS,3),2);  % Average of 3D model
    z = gridStruct.z;
else
    [z, ~, vp, vs, ~, ~] = ak135('cont');
end

dz   = gridStruct.dz;
zmax = max(gridStruct.z);
zout = 0:dz:zmax;

%% 3. Ray Tracing and Time-to-Depth Migration
nrf     = length(gather);
rfsAll  = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
timeAll = cellfun(@(rf) rf.ittime, {gather.RF}, 'UniformOutput', false);
raypAll = cell2mat(cellfun(@(ti) ti.rayParam / 6371, {gather.TravelInfo}, 'UniformOutput', false));
bazAll  = cell2mat(cellfun(@(ti) ti.baz, {gather.TravelInfo}, 'UniformOutput', false));
latAll  = cell2mat(cellfun(@(si) si.stla, {gather.StationInfo}, 'UniformOutput', false));
lonAll  = cell2mat(cellfun(@(si) si.stlo, {gather.StationInfo}, 'UniformOutput', false));

disp('Starting ray tracing...');
tic;
[cp, ~, MidPoints] = rf_ccp(raypAll, bazAll, dz, zmax, z, vp, vs, latAll, lonAll, 'flat');
toc;

% Apply 3D/2D velocity corrections if model is available
if ismember(gridStruct.ModelType, {'2D', '3D'})
    RayDepths = (1*dz:dz:zmax)';
    [TimeCorrections, ~, ~] = correct_RFs(MidPoints, RayDepths, gridStruct.Fvp, gridStruct.Fvs, z, vp, vs);
else
    TimeCorrections = zeros(length(zout), nrf);
end

disp('Migrating RFs from time to depth...');
[~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, dz, zmax, z, vp, vs, TimeCorrections);

% Assign migrated amplitudes and projected coordinates to the CP structure
for k = 1:nrf
    cp(k).amp = rfsAll_depth{k};
    [rx, ry]  = latlonToProjectedCoords([cp(k).lonb], [cp(k).latb], gridStruct);
    cp(k).rx  = rx;
    cp(k).ry  = ry;
end

%% 4. Pre-calculate Fresnel Zone Parameters
% Estimate dominant frequency from Gaussian parameter (a). Approx f ~ a/2.
f_center = param.gauss / 2; 

% Interpolate Vs to the imaging depth grid for Fresnel radius calculation
vs_interp = interp1(z, vs, gridStruct.z, 'linear', 'extrap');
% Physical Fresnel Zone Radius: R = sqrt(Vs * z / (2 * f))
R_fresnel = sqrt((vs_interp(:) .* gridStruct.z(:)) ./ (2 * f_center));

% Set a minimum radius to prevent artifacts in the shallow crust and fill gaps
R_min = gridStruct.dx * 1.5;
R_fresnel(R_fresnel < R_min) = R_min;

% Gaussian decay factor (sigma = R / sigma_factor)
sigma_factor = 2; 

%% 5. CCP Stacking Process
switch param.imagingType
    case '2D'
        [X, Z] = meshgrid(gridStruct.x, gridStruct.z);
        nx = length(gridStruct.x);
        nz = length(gridStruct.z);
        V = zeros(nz, nx);
        count = zeros(nz, nx);
        dx = gridStruct.dx;

        for n = 1:length(cp)
            xx = cp(n).rx;   % Horizontal distance along profile
            zz = cp(n).zpos; % Depth
            amp_vec = cp(n).amp;
            
            for k = 1:length(zz)
                if isnan(amp_vec(k)), continue; end
                
                % Map depth to grid index
                z_idx = round((zz(k) - gridStruct.z(1)) / dz) + 1;
                if z_idx < 1 || z_idx > nz, continue; end
                
                R = R_fresnel(z_idx);
                % Determine horizontal index range within 2 * Radius
                x_idx_range = max(1, floor((xx(k) - 2*R - gridStruct.x(1)) / dx) + 1) : ...
                              min(nx, ceil((xx(k) + 2*R - gridStruct.x(1)) / dx) + 1);
                
                for ix = x_idx_range
                    dist_x = abs(gridStruct.x(ix) - xx(k));
                    if dist_x <= 2*R
                        % Apply Gaussian weight based on distance
                        weight = exp(-(dist_x^2) / (2 * (R/sigma_factor)^2));
                        V(z_idx, ix) = V(z_idx, ix) + amp_vec(k) * weight;
                        count(z_idx, ix) = count(z_idx, ix) + weight;
                    end
                end
            end
        end
        % Normalize by total weight to maintain amplitude fidelity
%         img_final = V ./ (count + eps);
        ccpResult = struct('X', X, 'Z', Z, 'img', V, 'count', count);

    case '3D'
        nx = gridStruct.nx; ny = gridStruct.ny; nz = gridStruct.nz;
        [X, Y, Z] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);
        V = zeros(ny, nx, nz);
        count = zeros(ny, nx, nz);
        dx = gridStruct.dx; dy = gridStruct.dy;

        for n = 1:length(cp)
            if mod(n, 50) == 0, fprintf('Stacking trace %d/%d...\n', n, length(cp)); end
            xx = cp(n).rx; yy = cp(n).ry; zz = cp(n).zpos; amp_vec = cp(n).amp;

            for k = 1:length(zz)
                if isnan(amp_vec(k)), continue; end
                
                z_idx = round((zz(k) - gridStruct.z(1)) / dz) + 1;
                if z_idx < 1 || z_idx > nz, continue; end
                
                R = R_fresnel(z_idx);
                % Determine 2D bounding box for the Gaussian kernel
                ix_range = max(1, floor((xx(k)-2*R-gridStruct.x(1))/dx)+1) : min(nx, ceil((xx(k)+2*R-gridStruct.x(1))/dx)+1);
                iy_range = max(1, floor((yy(k)-2*R-gridStruct.y(1))/dy)+1) : min(ny, ceil((yy(k)+2*R-gridStruct.y(1))/dy)+1);
                
                for ix = ix_range
                    for iy = iy_range
                        dist_sq = (gridStruct.x(ix)-xx(k))^2 + (gridStruct.y(iy)-yy(k))^2;
                        if dist_sq <= (2*R)^2
                            weight = exp(-dist_sq / (2 * (R/sigma_factor)^2));
                            V(iy, ix, z_idx) = V(iy, ix, z_idx) + amp_vec(k) * weight;
                            count(iy, ix, z_idx) = count(iy, ix, z_idx) + weight;
                        end
                    end
                end
            end
        end
%         img_final = V ./ (count + eps);
        ccpResult = struct('X', X, 'Y', Y, 'Z', Z, 'img', V, 'count', count);
end
end