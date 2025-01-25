function velocityModel = getVelocityModel(ModelParam)

    minlat = ModelParam.LatMin;
    maxlat = ModelParam.LatMax;
    minlon = ModelParam.LonMin;
    maxlon = ModelParam.LonMax;
    npt = ModelParam.spacing;

    disp('Create velocity model starts')
    X = [];
    Y = [];
    Z = [];
    VP = [];
    VS = [];
    knode = 0;
    latall = linspace(minlat,maxlat,npt);
    lonall = linspace(minlon,maxlon,npt);
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
            if_topo = 1;
            % [m0,nsedi] = obtain_crust1_QB(lat,lon,[],if_topo);
            [m0,nsedi] = obtain_crust1_v2(lat,lon,[],if_topo);
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
            % find the depth beneath the moho from prem
    %         prem_model = prem;
    %         z_prem = prem_model.depth;
    %         vp_prem = prem_model.vp;
    %         vs_prem = prem_model.vs;
    %         rho_prem = prem_model.rho;
    %         keepz = z_prem > dmoho(knode) & z_prem < 1000;
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
    disp('Create velocity model completes')

    velocityModel.vp = Fvp;
    velocityModel.vs = Fvs;

end