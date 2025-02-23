function ccpResult = CCPCommonEventGather(gather, velocityModel, gridStruct)
% CCPCommonEventGather  Perform Common Conversion Point stacking (CCP) on seismic data.
%
% Usage:
%   ccpResult = CCPCommonEventGather(gather, velocityModel, gridStruct)
%
% Inputs:
%   gather         : Struct array of seismic traces for one event or gather.
%                    Each element includes fields:
%                       .RF.itr            - Iterative decon RF [Nt x 1]
%                       .TravelInfo.rayParam, .TravelInfo.baz - Ray parameters
%                       .TimeAxis.t_resample, .TimeAxis.dt_resample - time axis
%   velocityModel  : Struct describing velocity model with fields:
%                       .x, .z, .nx, .nz, .dx, .dz, .vp, .vs, ...
%   gridStruct  : Struct with profile-related info:

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
% Validate gather
if isempty(gather) || ~isstruct(gather)
    error('CCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end

% Validate velocityModel fields
% requiredVMfields = {'x', 'z', 'vp', 'vs', 'dx', 'dz', 'nx', 'nz'};
% for f = requiredVMfields
%     if ~isfield(velocityModel, f{1})
%         error('CCPCommonEventGather:MissingVelocityModelField', 'Field "%s" is missing in velocityModel.', f{1});
%     end
% end


gridType = '3D';

% Display event info if available
if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    disp(['Processing event: ', gather(1).EventInfo.evid]);
end

%% Unpack Profile and Velocity Model Information
% Unpack velocity model fields
% vp   = velocityModel.vp(:, 1);  % Assuming vp is 1D
% vs   = velocityModel.vs(:, 1);  % Assuming vs is 1D
[z, r, vp, vs, ~, ~] = ak135('cont');

dz = gridStruct.dz;
zmax = max(gridStruct.z);

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
[cp, ~ ,MidPoints] = rf_ccp(raypAll, bazAll, gridStruct.dz, zmax, z, vp, vs, latAll, lonAll, model_type);
toc;
disp('Ray tracing completed');

Fvp = velocityModel.vp;
Fvs = velocityModel.vs;
% correct for heterogeneity
RayDepths = (1*dz:dz:zmax)';
[TimeCorrections, ~, ~] = correct_RFs(MidPoints, RayDepths, Fvp, Fvs, z, vp, vs);

% Perform time-to-depth conversion for RFs
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
switch gridType
    % project to profile for 2D imaging
    case '2D'
        % Set grid sizes (distance direction: dx, depth direction: dz)
        [X, Z] = meshgrid(gridStruct.x, gridStruct.z);
        nx = length(gridStruct.x);
        nz = length(gridStruct.z);

        % Initialize stacking result matrices
        V = zeros(nz, nx); % Store accumulated amplitude values
        count = zeros(nz, nx); % Store sample count for each grid
        % Perform parallel processing for speed
        parfor i = 1:nz
            for j = 1:nx
                for n = 1:length(cp)
                    xi = cp(n).rx;
                    zi = cp(n).zpos;
                    vi = cp(n).amp;
                    keep = xi >= gridStruct.x(j) - 2 * gridStruct.dx & xi <= gridStruct.x(j) + 2 * gridStruct.dx & ...
                        zi >= gridStruct.z(i) - 2 * gridStruct.dz & zi <= gridStruct.z(i) + 2 * gridStruct.dz;
                    V(i, j) = V(i, j) + sum(vi(keep));
                    count(i, j) = count(i, j) + sum(keep);
                end
            end
        end

        % Generate CCP Image
        V = V ./ max(count, 1);  % 避免除以0

        %% Plot CCP Image
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
        xmin = min(gridStruct.x);
        ymin = min(gridStruct.y);
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
                i=floor((yy(k)-ymin)/dy)+1;
                j=floor((xx(k)-xmin)/dx)+1;
                if i>0 && j>0 && i<=ny && j<=nx
                    amp=cp(n).amp(k);
                    if ~isnan(amp)
                        count(i,j,k)=count(i,j,k)+1;
                        V(i,j,k)=V(i,j,k)+amp;
                    end
                end
            end

        end
        V = V./max(count,1);

        % V1 = zeros(nx, ny, nz); % Store accumulated amplitude values
        % count1 = zeros(nx, ny, nz); % Store sample count for each grid       
        % % Perform parallel processing for speed        
        %         parfor i = 1:nz
        %             for j = 1:nx
        %                 for k = 1:ny
        %                     for n = 1:length(cp)
        %                         xi = cp(n).rx;
        %                         yi = cp(n).ry;
        %                         zi = cp(n).zpos;
        %                         vi = cp(n).amp;
        %                         keep = xi >= gridStruct.x(j) - 2 * gridStruct.dx & xi <= gridStruct.x(j) + 2 * gridStruct.dx & ...
        %                             yi >= gridStruct.y(k) - 2 * gridStruct.dy & yi <= gridStruct.y(k) + 2 * gridStruct.dy & ...
        %                             zi >= gridStruct.z(i) - 2 * gridStruct.dz & zi <= gridStruct.z(i) + 2 * gridStruct.dz;
        %                         V1(j, k, i) = V1(j, k, i) + sum(vi(keep));
        %                         count1(j, k, i) = count1(j, k, i) + sum(keep);
        %                     end
        %                 end
        %             end
        %         end
        %         V1 = V1./max(count1,1);
        %         figure; 
        %         h = slice(permute(V,[2,3,1]),10,10,10);
        %         set(h(:),'EdgeColor','none')
        %         set(gca,'ZDir','reverse')
        %         colormap(flipud(roma));

        ccpResult = struct('X', X, 'Y', Y, 'Z', Z, 'img', V, 'count', count);
end
end
