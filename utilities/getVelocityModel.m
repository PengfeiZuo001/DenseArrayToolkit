function velocityModel = getVelocityModel(ModelType, varargin)

switch  ModelType
    case '2D'
        xpad = 50;
        xmax = 350;
        model = obtain_crust1_QB();
        dx = 4;
        x = 0-xpad:dx:xmax+xpad;
        nx = length(x);
        dz = 1;
        z = 0:dz:100;
        nz = length(z);
        vp = interp1(model(:,1),model(:,2),z,'nearest','extrap');
        vs = interp1(model(:,1),model(:,3),z,'nearest','extrap');
        vel = repmat(vp(:),1,nx);
        vel_s = repmat(vs(:),1,nx);
        
        % smooth the velocity model
        N=5;
        [vel,~]=moving_avg(vel,N,'constant',2);
        [vel,~]=moving_avg(vel,N,'constant');
        [vel_s,~]=moving_avg(vel_s,N,'constant',2);
        [vel_s,~]=moving_avg(vel_s,N,'constant');
        
        velocityModel.x = x;
        velocityModel.z = z;
        velocityModel.xmax = xmax;
        velocityModel.nx = nx;
        velocityModel.nz = nz;
        velocityModel.dx = dx;
        velocityModel.dz = dz;
        velocityModel.xpad = xpad;
        velocityModel.vp = vel;
        velocityModel.vs = vel_s;
    
    case '3D'

        LatMin = varargin{1, 1}.LatMin;
        LatMax = varargin{1, 1}.LatMax;
        LonMin = varargin{1, 1}.LonMin;
        LonMax = varargin{1, 1}.LonMax;

        npts = varargin{1, 1}.npts;

        [z, rho, vp, vs, qk, qm] = ak135( 'cont' );
        zmax = 800;
        dz = 0.5;

        X = [];
        Y = [];
        Z = [];
        VP = [];
        VS = [];
        knode = 0;
        latall = linspace(LatMin,LatMax,npts);
        lonall = linspace(LonMin,LonMax,npts);
        % plot it on the CCP map
        % [lontemp,lattemp] = meshgrid(lonall,latall);
        % scatter(lontemp(:),lattemp(:),'ko','filled')
        for i = 1:length(latall)
            for j = 1:length(lonall)
                knode = knode + 1;
                disp(['Now extracting velocity model at node #',num2str(knode)]);
                lat = latall(i);
                lon = lonall(j);
                % note that m0 starts at 0 km depth. The elevation infomration is
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
                dmoho(knode) = m0(end,1);
                dmax = 1000;
                keepz=z>dmoho(knode) & z<dmax;
                % deal with the upper mantle, use the PREM model
                m_interface = [m_interface; z(keepz) vp(keepz) vs(keepz) rho(keepz)];
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
                % plot(vptemp,zpos(1:end-1)); hold on;
                % ylim([-5 100])
                % axis ij;
            end
        end
        % plot the model
        % imagesc(1:knode,zpos(1:end-1),reshape(VS,length(zpos)-1,knode));
        % interpolate 
        Fvp = scatteredInterpolant(X,Y,Z,VP);
        Fvs = scatteredInterpolant(X,Y,Z,VS);

        velocityModel.vp = Fvp;
        velocityModel.vs = Fvs;

end
% figure;
% imagesc(x,z,vel); % hold on; scatter(rx,zeros(size(rx)),100,'r^','filled')
end
