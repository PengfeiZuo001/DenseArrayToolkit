function ccpResult = CCPCommonEventGather(gather, gridStruct, param)
% CCPCommonEventGather  Perform Common Conversion Point stacking (CCP) on seismic data.
%
% Usage:
%   ccpResult = CCPCommonEventGather(gather, gridStruct)
%
% Inputs:
%   gather         : Struct array of seismic traces for one event or gather.
%                    Each element includes fields:
%                       .RF.itr            - Iterative decon RF [Nt x 1]
%                       .TravelInfo.rayParam, .TravelInfo.baz - Ray parameters
%                       .TimeAxis.t_resample, .TimeAxis.dt_resample - time axis
%   gridStruct  : Struct with grid-related info:
%
%   param          : CCP parameters:
%       .imagingType (string) - conduct 2D or 3D imaging
%       .plotCCP     (bool) - if true, plot CCP results
%       .smoothLength (int) - if greater than 0, apply smoothing to CCP
%                             image

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
% Validate input data structure and parameters
% --------------------------------------------------
% Validate seismic gather structure
% - Check if gather is non-empty struct array
% - Each gather element should contain RF data and travel time information
if isempty(gather) || ~isstruct(gather)
    error('CCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end

% Display event info if available
if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    disp(['Processing event: ', gather(1).EventInfo.evid]);
end

%% Unpack Grid Information
% Extract grid configuration
% --------------------------------------------------
% gridStruct contains either 1D or 3D velocity information and
% defines the imaging grid parameters
if strcmp(gridStruct.ModelType ,'1D')
    vp   = gridStruct.vp(:, 1);  % Assuming vp is 1D
    vs   = gridStruct.vs(:, 1);  % Assuming vs is 1D
    z = gridStruct.z;
else
    [z, r, vp, vs, ~, ~] = ak135('cont');
end
dz = gridStruct.dz;
zmax = max(gridStruct.z);
zout = 0:dz:zmax;
%% Calculate Conversion Points and Ray Tracing
% --------------------------------------------------
% Key Steps:
% 1. Extract receiver function (RF) data and time axes
% 2. Calculate ray parameters from travel time information
% 3. Perform ray tracing to find conversion points
% 4. Apply 3D velocity corrections if needed
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
[cp, ~ ,MidPoints] = rf_ccp(raypAll, bazAll, gridStruct.dz, zmax, z, vp, vs, latAll, lonAll, model_type);
toc;
disp('Ray tracing completed');

% make time correction if a regional 3D velocity model is available
if strcmp(gridStruct.ModelType ,'3D')
    Fvp = gridStruct.Fvp;
    Fvs = gridStruct.Fvs;
    % correct for heterogeneity
    RayDepths = (1*dz:dz:zmax)';
    [TimeCorrections, ~, ~] = correct_RFs(MidPoints, RayDepths, Fvp, Fvs, z, vp, vs);
else
    TimeCorrections = zeros(length(zout),length(rfsAll));
end

% Time-to-Depth Conversion
% --------------------------------------------------
% Convert RFs from time domain to depth domain using:
% - Ray parameters (raypAll)
% - Velocity model (vp, vs)
% - Time corrections from 3D velocity model (if applicable)
% This converts the time-based RF measurements to depth coordinates
% matching the imaging grid
disp('Time-to-depth conversion started');
tic;
[~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, gridStruct.dz, zmax, z, vp, vs, TimeCorrections);
toc;
disp('Time-to-depth conversion completed');


% Attach RF amplitudes to CCP points
for k = 1:nrf
    cp(k).amp = rfsAll_depth{k};
end

% Project CCP points to profile
for k = 1:nrf
    slon_tmp = [cp(k).lonb];
    slat_tmp = [cp(k).latb];
    [rx, ry] = latlonToProjectedCoords(slon_tmp, slat_tmp, gridStruct);
    cp(k).rx = rx;
    cp(k).ry = ry;
end

%% CCP Stacking Process
% --------------------------------------------------
% Core algorithm steps:
% 1. Initialize imaging grid based on gridStruct parameters
% 2. Bin RF amplitudes into spatial grid cells
% 3. Stack (sum) amplitudes in each cell
% 4. Normalize by sample count in each cell
%
% Handles both 1D and 3D velocity models differently:
% - 1D: Simple 2D (distance-depth) stacking
% - 3D: Full 3D spatial binning with progress tracking
switch param.imagingType
    % project to profile for 2D imaging
    case '2D'
        % Set grid sizes (distance direction: dx, depth direction: dz)
        [X, Z] = meshgrid(gridStruct.x, gridStruct.z);
        nx = length(gridStruct.x);
        nz = length(gridStruct.z);

        % Initialize stacking result matrices
        V = zeros(nz, nx); % Store accumulated amplitude values
        count = zeros(nz, nx); % Store sample count for each grid
        dx = gridStruct.dx;      

        for n = 1:length(cp)
            % 2D stacking
            xx=cp(n).rx;
            zz=cp(n).zpos;
            for j=1:length(zz)
                col_idx=floor((xx(j)-gridStruct.x(1))/dx)+1;
                % Calculate z-layer index
                row_idx = floor((zz(j) - gridStruct.z(1)) / gridStruct.dz) + 1;
                % Validate array indices
                if col_idx>=1 && col_idx<=nx && row_idx>=1 && row_idx<=nz
                    amp = cp(n).amp(j);
                    if ~isnan(amp)
                        count(row_idx,col_idx) = count(row_idx,col_idx) + 1;
                        V(row_idx,col_idx) = V(row_idx,col_idx) + amp;
                    end
                end
            end
        end

        if param.smoothLength > 0 
            K = (1/param.smoothLength^2)*ones(param.smoothLength,param.smoothLength);
            V = conv2(V,K,'same');
        end

        % Generate CCP Image
%         V = V ./ max(count, 1);  % 避免除以0

        %% Plot CCP Image
        if param.plotCCP
            try
                load roma;
                cmap = flipud(roma);
            catch
                cmap = parula;
            end
    
            figure;
            set(gcf, 'Position', [100 100 800 400], 'color', 'w');
            imagesc(gridStruct.x, gridStruct.z, V);
            caxis([-0.1 0.1]);
            colormap(cmap);
            colorbar;
            xlabel('Distance (km)');
            ylabel('Depth (km)');
            title('CCP image');
            set(gca, 'fontsize', 14);
        end

        %% Save CCP results
        ccpResult = struct('X', X, 'Z', Z, 'img', V, 'count', count);
    case '3D'
        nx = gridStruct.nx;
        ny = gridStruct.ny;
        nz = gridStruct.nz;
        % Set grid sizes (distance direction: dx, depth direction: dz)
        [X, Y, Z] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);

        % Initialize stacking result matrices
        V = zeros(ny, nx, nz); % Store accumulated amplitude values
        count = zeros(ny, nx, nz); % Store sample count for each grid
        dx = gridStruct.dx;
        dy = gridStruct.dy;

        for n=1:length(cp)
            if mod(n,10) == 0
                disp(['Binning ',num2str(n),'/',num2str(length(cp)),' traces']);
            end
            xx=cp(n).rx;
            yy=cp(n).ry;
            zz=cp(n).zpos;

            % 3D stacking
            for k=1:length(zz)
                yy_idx=floor((yy(k)-gridStruct.y(1))/dy)+1;
                xx_idx=floor((xx(k)-gridStruct.x(1))/dx)+1;
                % Calculate z-layer index
                zz_idx = floor((zz(k) - gridStruct.z(1)) / gridStruct.dz) + 1;
                % Validate array indices
                if yy_idx>=1 && yy_idx<=ny && xx_idx>=1 && xx_idx<=nx && zz_idx>=1 && zz_idx<=nz
                    amp = cp(n).amp(k);
                    if ~isnan(amp)
                        count(yy_idx,xx_idx,zz_idx) = count(yy_idx,xx_idx,zz_idx) + 1;
                        V(yy_idx,xx_idx,zz_idx) = V(yy_idx,xx_idx,zz_idx) + amp;
                    end
                end
            end

        end
%         % apply smoothing to the CCP image
        if param.smoothLength > 0
            V = smooth3(V,'box',param.smoothLength);
            count = smooth3(count,'box',param.smoothLength);
        end

        ccpResult = struct('X', X, 'Y', Y, 'Z', Z, 'img', V, 'count', count);
end

end
