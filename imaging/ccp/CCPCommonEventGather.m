function ccpResult = CombinedCCPCommonEventGather(gather, gridStruct, param)
% CombinedCCPCommonEventGather  CCP stacking with selectable uniform or Fresnel-weighted mode.
%
% Usage:
%   ccpResult = CombinedCCPCommonEventGather(gather, gridStruct, param)
%
% Inputs:
%   gather      : Struct array of seismic traces (RFs and travel info)
%   gridStruct  : Struct containing imaging grid and velocity model
%   param       : Parameter struct:
%       .imagingType  - '2D' or '3D'
%       .stackMode    - 'uniform' (均匀网格) or 'fresnel' (Fresnel带加权)
%       .gauss        - (仅fresnel模式) Gaussian滤波参数，估算主频
%       .plotCCP      - 是否绘图
%       .smoothLength - (仅uniform模式) 平滑长度
%
% Outputs:
%   ccpResult   : Result struct (X, Y, Z coordinates, img, and weights/count)
%
% Author: 合并自CCPCommonEventGather与FresnelCCPCommonEventGather
% Date: Feb. 2026

%% 1. 输入检查与参数
if isempty(gather) || ~isstruct(gather)
    error('CombinedCCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end
if ~isfield(param, 'stackMode')
    param.stackMode = 'uniform';
end
if ~isfield(param, 'imagingType')
    param.imagingType = '2D';
end

%% 2. 网格与速度模型
if strcmp(gridStruct.ModelType ,'1D')
    vp = gridStruct.vp(:, 1);
    vs = gridStruct.vs(:, 1);
    z  = gridStruct.z;
elseif strcmp(gridStruct.ModelType ,'2D')
    vp = mean(gridStruct.vp, 'all'); 
    vs = mean(gridStruct.vs, 'all');
    z  = gridStruct.z;
elseif strcmp(gridStruct.ModelType ,'3D')
    vp = mean(mean(gridStruct.VP,3),2);
    vs = mean(mean(gridStruct.VS,3),2);
    z  = gridStruct.z;
else
    [z, ~, vp, vs, ~, ~] = ak135('cont');
end
dz   = gridStruct.dz;
zmax = max(gridStruct.z);
zout = 0:dz:zmax;

%% 3. 射线追踪与时深变换
nrf     = length(gather);
rfsAll  = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
timeAll = cellfun(@(rf) rf.ittime, {gather.RF}, 'UniformOutput', false);
raypAll = cell2mat(cellfun(@(ti) ti.rayParam / 6371, {gather.TravelInfo}, 'UniformOutput', false));
bazAll  = cell2mat(cellfun(@(ti) ti.baz, {gather.TravelInfo}, 'UniformOutput', false));
latAll  = cell2mat(cellfun(@(si) si.stla, {gather.StationInfo}, 'UniformOutput', false));
lonAll  = cell2mat(cellfun(@(si) si.stlo, {gather.StationInfo}, 'UniformOutput', false));

[cp, ~, MidPoints] = rf_ccp(raypAll, bazAll, dz, zmax, z, vp, vs, latAll, lonAll, 'flat');

if ismember(gridStruct.ModelType, {'2D', '3D'})
    RayDepths = (1*dz:dz:zmax)';
    [TimeCorrections, ~, ~] = correct_RFs(MidPoints, RayDepths, gridStruct.Fvp, gridStruct.Fvs, z, vp, vs);
else
    TimeCorrections = zeros(length(zout), nrf);
end

[~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, dz, zmax, z, vp, vs, TimeCorrections);

for k = 1:nrf
    cp(k).amp = rfsAll_depth{k};
    [rx, ry]  = latlonToProjectedCoords([cp(k).lonb], [cp(k).latb], gridStruct);
    cp(k).rx  = rx;
    cp(k).ry  = ry;
end

%% 4. 预计算Fresnel带参数（仅fresnel模式）
if strcmpi(param.stackMode, 'fresnel')
    if ~isfield(param, 'gauss')
        param.gauss = 2; % 默认主频估算
    end
    f_center = param.gauss / 2;
    vs_interp = interp1(z, vs, gridStruct.z, 'linear', 'extrap');
    R_fresnel = sqrt((vs_interp(:) .* gridStruct.z(:)) ./ (2 * f_center));
    R_min = gridStruct.dx * 1.5;
    R_fresnel(R_fresnel < R_min) = R_min;
    sigma_factor = 2;
end

%% 5. CCP叠加主流程
switch param.imagingType
    case '2D'
        [X, Z] = meshgrid(gridStruct.x, gridStruct.z);
        nx = length(gridStruct.x);
        nz = length(gridStruct.z);
        V = zeros(nz, nx);
        count = zeros(nz, nx);
        dx = gridStruct.dx;

        for n = 1:length(cp)
            xx = cp(n).rx;
            zz = cp(n).zpos;
            amp_vec = cp(n).amp;
            for k = 1:length(zz)
                if isnan(amp_vec(k)), continue; end
                z_idx = round((zz(k) - gridStruct.z(1)) / dz) + 1;
                if z_idx < 1 || z_idx > nz, continue; end

                if strcmpi(param.stackMode, 'fresnel')
                    R = R_fresnel(z_idx);
                    x_idx_range = max(1, floor((xx(k) - 2*R - gridStruct.x(1)) / dx) + 1) : ...
                                  min(nx, ceil((xx(k) + 2*R - gridStruct.x(1)) / dx) + 1);
                    for ix = x_idx_range
                        dist_x = abs(gridStruct.x(ix) - xx(k));
                        if dist_x <= 2*R
                            weight = exp(-(dist_x^2) / (2 * (R/sigma_factor)^2));
                            V(z_idx, ix) = V(z_idx, ix) + amp_vec(k) * weight;
                            count(z_idx, ix) = count(z_idx, ix) + weight;
                        end
                    end
                else % 均匀网格
                    col_idx = floor((xx(k)-gridStruct.x(1))/dx)+1;
                    if col_idx>=1 && col_idx<=nx
                        count(z_idx,col_idx) = count(z_idx,col_idx) + 1;
                        V(z_idx,col_idx) = V(z_idx,col_idx) + amp_vec(k);
                    end
                end
            end
        end

        % 平滑（仅uniform模式）
        if strcmpi(param.stackMode, 'uniform') && isfield(param, 'smoothLength') && param.smoothLength > 0
            K = (1/param.smoothLength^2)*ones(param.smoothLength,param.smoothLength);
            V = conv2(V,K,'same');
            count = conv2(count,K,'same');
        end

        ccpResult = struct('X', X, 'Z', Z, 'img', V, 'count', count);

    case '3D'
        nx = gridStruct.nx; ny = gridStruct.ny; nz = gridStruct.nz;
        [X, Y, Z] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);
        V = zeros(ny, nx, nz);
        count = zeros(ny, nx, nz);
        dx = gridStruct.dx; dy = gridStruct.dy;

        for n = 1:length(cp)
            xx = cp(n).rx; yy = cp(n).ry; zz = cp(n).zpos; amp_vec = cp(n).amp;
            for k = 1:length(zz)
                if isnan(amp_vec(k)), continue; end
                z_idx = round((zz(k) - gridStruct.z(1)) / dz) + 1;
                if z_idx < 1 || z_idx > nz, continue; end

                if strcmpi(param.stackMode, 'fresnel')
                    R = R_fresnel(z_idx);
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
                else % 均匀网格
                    xx_idx=floor((xx(k)-gridStruct.x(1))/dx)+1;
                    yy_idx=floor((yy(k)-gridStruct.y(1))/dy)+1;
                    if yy_idx>=1 && yy_idx<=ny && xx_idx>=1 && xx_idx<=nx
                        count(yy_idx,xx_idx,z_idx) = count(yy_idx,xx_idx,z_idx) + 1;
                        V(yy_idx,xx_idx,z_idx) = V(yy_idx,xx_idx,z_idx) + amp_vec(k);
                    end
                end
            end
        end

        % 平滑（仅uniform模式）
        if strcmpi(param.stackMode, 'uniform') && isfield(param, 'smoothLength') && param.smoothLength > 0
            V = smooth3(V,'box',param.smoothLength);
            count = smooth3(count,'box',param.smoothLength);
        end

        ccpResult = struct('X', X, 'Y', Y, 'Z', Z, 'img', V, 'count', count);
end

end