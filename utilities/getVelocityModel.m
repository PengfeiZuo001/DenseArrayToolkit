function gridStruct = getVelocityModel(ModelType, gridStruct, nptsx,nptsy)
%GETVELOCITYMODEL Generate velocity model based on specified type and parameters
%
% Inputs:
%   ModelType  - String, type of velocity model: '1D', '2D', or '3D'
%   gridStruct - Structure containing grid information:
%                For 1D: nx, z
%                For 3D: x, y, z, LatMin, LatMax, LonMin, LonMax
%   npts       - Number of points for interpolation (used in 3D model)
%
% Outputs:
%   gridStruct - Updated structure with velocity information:
%                For 1D: adds vp, vs fields
%                For 3D: adds VP, VS, Fvp, Fvs fields
%
% Note: 2D implementation is currently placeholder only

% Input validation
if ~ischar(ModelType) || ~ismember(ModelType, {'1D', '2D', '3D'})
    error('ModelType must be one of: 1D, 2D, or 3D');
end


gridStruct.ModelType = ModelType;
switch ModelType
    case '1D'
        % Get 1D velocity model from CRUST1.0
        model = obtain_crust1_QB();
        nx = gridStruct.nx;
        
        % Interpolate P and S wave velocities to grid points
        vp = interp1(model(:,1), model(:,2), gridStruct.z, 'nearest', 'extrap');
        vs = interp1(model(:,1), model(:,3), gridStruct.z, 'nearest', 'extrap');
        vel = repmat(vp(:), 1, nx);
        vel_s = repmat(vs(:), 1, nx);
        
        % Apply smoothing to velocity models using moving average
        N = 5;  % Smoothing window size
        [vel,~] = moving_avg(vel, N, 'constant', 2);
        [vel,~] = moving_avg(vel, N, 'constant');
        [vel_s,~] = moving_avg(vel_s, N, 'constant', 2);
        [vel_s,~] = moving_avg(vel_s, N, 'constant');

        gridStruct.vp = vel;
        gridStruct.vs = vel_s;
        
    case '2D'
        % TODO: Implement 2D velocity model
%         warning('2D velocity model not implemented yet');
        % extract 2D velocity profile
%         X = gridStruct.XInOriginalCoord;
%         Y = gridStruct.YInOriginalCoord;
%         [LON, LAT] = xy2latlon(X, Y, gridStruct.originLon, gridStruct.originLat);

        % load regional velocity model
        filename = './velocity_model/Zhao2013_QB_PS2.2.txt';
        matrix = readmatrix(filename);
        lon = matrix(:,1); 
        lat = matrix(:,2); 
        depth = matrix(:,3);
        vp = matrix(:,4); 
        vs = matrix(:,5);
        
        % Create interpolants for the 3D velocity model
        Fvp = scatteredInterpolant(lon, lat, depth, vp, 'linear', 'nearest');
        Fvs = scatteredInterpolant(lon, lat, depth, vs, 'linear', 'nearest');
        gridStruct.Fvp = Fvp;
        gridStruct.Fvs = Fvs;

        % velocity model grid
        [xv,yv] = latlonToProjectedCoords(lon, lat, gridStruct);
        % migrition grid
        [XX, YY, ZZ] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);
        
        % 使用散点插值到规则网格
        F_vp = scatteredInterpolant(xv, yv, depth, vp, 'linear', 'none');
        F_vs = scatteredInterpolant(xv, yv, depth, vs, 'linear', 'none');
        
        % 在规则网格点上求值
        Vp = F_vp(XX, YY, ZZ);
        Vs = F_vs(XX, YY, ZZ);
        
        % 转置以匹配期望的维度顺序 [nz,nx,ny]
        Vp = permute(Vp, [3 2 1]);
        Vs = permute(Vs, [3 2 1]);

%         gridStruct.VP = Vp;
%         gridStruct.VS = Vs;
        vel = mean(Vp,3);
        vel_s = mean(Vs,3);
        % Apply smoothing to velocity models using moving average
        N = 5;  % Smoothing window size
        [vel,~] = moving_avg(vel, N, 'constant', 2);
        [vel,~] = moving_avg(vel, N, 'constant');
        [vel_s,~] = moving_avg(vel_s, N, 'constant', 2);
        [vel_s,~] = moving_avg(vel_s, N, 'constant');
        
        gridStruct.vp = vel;
        gridStruct.vs = vel_s;

        figure;
        imagesc(gridStruct.x, gridStruct.z, vel);
    case '3D'
        % Extract geographical boundaries
        LatMin = gridStruct.LatMin;
        LatMax = gridStruct.LatMax;
        LonMin = gridStruct.LonMin;
        LonMax = gridStruct.LonMax;

        % Get AK135 reference model
        [z, rho, vp, vs, qk, qm] = ak135('cont');
        zmax = 100;  % ccp 800km
        dz = 1;

        % Initialize arrays for storing 3D model data
        X = []; Y = []; Z = []; VP = []; VS = [];
        knode = 0;

        % Create lat/lon grid for interpolation
        lonall = linspace(LonMin, LonMax, nptsx);
        latall = linspace(LatMin, LatMax, nptsy);

        vpgrid = zeros(zmax/dz+1,nptsx,nptsy);
        vsgrid = zeros(zmax/dz+1,nptsx,nptsy);

        % Loop through each grid point to build 3D model
        for i = 1:length(latall)
            for j = 1:length(lonall)
                knode = knode + 1;
                disp(['Processing node: ',num2str(knode), ' [',num2str(lonall(j)),',',num2str(latall(i)), '] of ', num2str(nptsx*nptsy)]);

                % Get velocity model at current location
                lat = latall(i);
                lon = lonall(j);
% note that                 m0 starts at 0 km depth. The elevation infomration is
                % missing, to include topography set if_topo to 1
                % if_topo = 1;
                % [m0,~] = obtain_crust1_v2(lat,lon,[],if_topo);
                m0 = obtain_ustc_litho(lat,lon);
                % define the interface in the crust 1.0 model
                m_interface = [];
                for l = 1:size(m0,1)
                    if l == 1
                        m_interface(l,:) = m0(l,:);
                    else
                        m_interface(2*(l-1),:) = [m0(l,1) m0(l-1,2:4)];
                        m_interface(2*(l-1)+1,:) = m0(l,:);
                    end
                end

                EPS = 1e-6;
                ztemp = m_interface(:,1);
                idisc = find( ztemp(1:end-1) == ztemp(2:end) );
                ztemp(idisc) = ztemp(idisc) - EPS;

                zpos = 0:dz:zmax;
                vptemp = interp1( ztemp, m_interface(:,2), zpos(1:end-1)+0.5*dz, 'linear','extrap');
                vstemp = interp1( ztemp, m_interface(:,3), zpos(1:end-1)+0.5*dz, 'linear','extrap');
                Z = [Z; zpos(1:end-1)'];
                X = [X; ones(size(vptemp'))*lon];
                Y = [Y; ones(size(vptemp'))*lat];
                VP = [VP; vptemp'];
                VS = [VS; vstemp'];

                vpgrid(:,j,i) = [vptemp,vptemp(end)];
                vsgrid(:,j,i) = [vstemp,vstemp(end)];

            end
        end
%         % Create interpolants for the 3D velocity model (used in ccp)
        Fvp = scatteredInterpolant(X, Y, Z, VP);
        Fvs = scatteredInterpolant(X, Y, Z, VS);
        %-----------------3D model---------------------------------
        gridStruct.Fvp = Fvp;
        gridStruct.Fvs = Fvs;
        %----------------------------------------------------------

        %----------smooth---------------------
        N = 5;  % Smoothing window size
        vp0 = zeros(size(vpgrid));
        vs0 = zeros(size(vpgrid));
        for i = 1:size(vpgrid,1)
            v1 = squeeze(vpgrid(i,:,:));
            [v1,~] = moving_avg(v1, N, 'constant', 2);
            [vp0(i,:,:),~] = moving_avg(v1, N, 'constant');
            
            v2 = squeeze(vsgrid(i,:,:));
            [v2,~] = moving_avg(v2, N, 'constant', 2);
            [vs0(i,:,:),~] = moving_avg(v2, N, 'constant');
        end
        %------------
        gridStruct.vp = vp0;
        gridStruct.vs = vs0;
        %% plot velocity model at 40 km depth
        figure;
        set(gcf,'Position',[0 0 1000 1000],'Color','w')
        hold on;

        idx = gridStruct.z == 40;
        V = squeeze(gridStruct.vs(idx,:,:));
        hm = pcolor(gridStruct.XInOriginalCoord,gridStruct.YInOriginalCoord,V');
        set(hm,'EdgeColor','none')
        cm = colormap('jet');
        colormap(flipud(cm));
        % 绘制台站的位置
        scatter(gridStruct.stationX, gridStruct.stationY, 50,'r^', 'filled', 'DisplayName', 'Stations','MarkerEdgeColor','k');

        % 绘制网格点的位置
        scatter(gridStruct.XInOriginalCoord(:), gridStruct.YInOriginalCoord(:), 10, 'k', 'filled', ...
            'DisplayName', 'Grid Points');

        % 设置图形
        xlabel('X (km)');
        ylabel('Y (km)');
%         legend('show','Location','best');
        axis equal;
        grid on;
        title('Velocity model');
        hold off;
        set(gca,'fontsize',14)


end
end