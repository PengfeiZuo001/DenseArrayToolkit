function ccpResult = CCPStacking(DataStruct, velocityModel, CCPParam)

    ccp_data_directory = './matfiles/CCPData/';

    [z, rho, vp, vs, qk, qm] = ak135( 'cont' );
    zmax = 800;
    dz = 0.5;
    Fvp = velocityModel.vp;
    Fvs = velocityModel.vs;

    LatMin = CCPParam.LatMin;
    LatMax = CCPParam.LatMax;
    LonMin = CCPParam.LonMin;
    LonMax = CCPParam.LonMax;
    BinSpacing = CCPParam.BinSpacing;
    BinSize = CCPParam.BinSize;

    stalist = CCPParam.StationCode;

    % call function to set up grid nodes
    CCPGrid = ccp_setup_grid(LatMin,LatMax,LonMin,LonMax,BinSpacing,BinSize);

    %% 1D raytracing and time to depth conversion
    disp('Migration starts')
    start_index = 1;
    cp_lat = [];
    cp_lon = [];

    for n = 1:length(stalist)
        disp(['Now processing station: ', stalist{n}]);
        seis = {};
        time = {};
        [CommonStationGather, ~] = getCommonStationGather(DataStruct,stalist{n});
        nevt = length(CommonStationGather);
        backaz = zeros(1,nevt);
        p = zeros(1,nevt);
        parfor ii = 1:nevt
            seis{ii} = CommonStationGather(ii).RF.itr;
            time{ii} = CommonStationGather(ii).RF.ittime;
            backaz(1,ii) = CommonStationGather(ii).TravelInfo.baz;
            p(1,ii) = CommonStationGather(ii).TravelInfo.rayParam;   % s/rad
        end
        nseis=length(seis);
        lat = CommonStationGather(1).StationInfo.stla;
        lon = CommonStationGather(1).StationInfo.stlo;

        % convert to s/km
        if mean(p)>100
            p = p./6371;
        end
        % 1D ray tracing
        [cp, RayMatrix, MidPoints] = findConversionPoints(p, backaz, dz, zmax, z, vp, vs, lat, lon, 'spherical');
        % save the conversion point location for GMT plot
        cp_lat = [cp_lat; cp.latb(60,:)'];
        cp_lon = [cp_lon; cp.lonb(60,:)'];
        RayDepths = 1*dz:dz:zmax;
        RayDepths = RayDepths(:); 
        % corret for heterogenity
        [TimeCorrections, Tpds3D, Tpds1D] = correct_RFs(MidPoints, RayDepths, Fvp, Fvs, z, vp, vs);
        % reindex the RFs
        end_index = start_index + nseis - 1;
        index = start_index:1:end_index;
        index_matrix = repmat(index,size(RayMatrix,1),1);
        RayMatrix(:,:,7) = index_matrix;
        % update the index
        start_index = end_index + 1;
        % map RFs to depth
        [timeout, seisout, depth0] = migrate_RFs( time, seis, p, dz, zmax, z, vp, vs, TimeCorrections);
        % normalization
        seisout_norm = normalize_RFs(seisout);
        % read in migrated RFs and depth0
        RayMatrix(:,:,1) = cell2mat(seisout_norm);
        RayMatrix(:,:,2) = repmat(depth0',1,nseis);
        %% save the Raymatrix
        if ~exist(ccp_data_directory,'dir')
            mkdir(ccp_data_directory)
        end
        CCPMatFile = [ccp_data_directory,'/CCPData',stalist{n},'.mat'];
        save(CCPMatFile,'RayMatrix','MidPoints');
    end
    disp('Migration completes')

    %% find the RFs fall within each bin
    disp('Assign RF starts')
    RayMatrix_all = [];
    for n = 1:length(stalist)
        %
        filename=[ccp_data_directory,'/CCPData',stalist{n},'.mat'];
        if exist(filename,'file')
            load(filename,'RayMatrix');
        else
            continue
        end
        % save matrix for individual station into a big matrix
        RayMatrix_all = [RayMatrix_all, RayMatrix];
    %     nptsTotal = size(CCPGrid,1) * size(RayMatrix,1);
    %     nptsUpdate = floor(nptsTotal / 100);
        kk = 0;
        for m = 1:size(CCPGrid,1) % Each Bin
            for k = 1:size(RayMatrix,1) % Each Depth
                lons = RayMatrix(k,:,4) * pi/180;
                lats = RayMatrix(k,:,3) * pi/180;
                tlon = CCPGrid{m,2} * pi/180;
                tlat = CCPGrid{m,1} * pi/180;
                dlon = lons - tlon;
                dlat = lats - tlat;
                a = (sin(dlat/2)).^2 + cos(lats) .* cos(tlat) .* (sin(dlon/2)).^2;
                angles = 2 .* atan2(sqrt(a),sqrt(1-a));
                dist = 6371 * angles;
                % It would be really easy to adjust the size of your bin as you
                % go to deeper depths.
                %--------------------------------------------------------------
                Indx = find(dist <= CCPGrid{m,3}); % This ONLY takes RFs falling into the bin
                Temp{k,1} = RayMatrix(k,Indx,7); % Record IDs
                TempDist{k,1} = dist(Indx);
                if ~isempty(Indx)
                    disp('')
                end
            end
            CCPGrid{m,4} = [CCPGrid{m,4} Temp];
            CCPGrid{m,5} = [CCPGrid{m,5} TempDist];
            disp(['Processed bin: ',num2str(m),'/',num2str(size(CCPGrid,1))])
        end
        
    end
    disp('Assign RF completes')

    %% ccp stacking
    disp('CCP stacking starts')
    for i = 1:size(CCPGrid,1) % For each bin
        disp(['Now stacking CCP bin #',num2str(i)]);
        RRFAmps = [];
        Weights = [];
        MaxHits(i) = 0;
        for k = 1:length(CCPGrid{i,4})
            temp(k) = numel([CCPGrid{i,4}{k,:}]);
        end
        MaxHits(i) = max(temp);
        clear temp
        if MaxHits(i) > 0
            % Build 2 matrices with MaxHits # columns and length(DepthAxis) rows
            Rtemp = NaN * ones(length(CCPGrid{i,5}),MaxHits(i));
            Wtemp = NaN * ones(length(CCPGrid{i,5}),MaxHits(i));
            % Add the value in the RayMatrix using the RF indices to the
            % matrices we just built.
            for k = 1:length(CCPGrid{i,4})
                Ids = cell2mat(CCPGrid{i,4}(k,:));
                %                 Dist = CCPGrid{n,4}{k};
                if ~isempty(Ids)
                    for l = 1:length(Ids)
                        temp = find(RayMatrix_all(k,:,7) == Ids(l));
                        if ~isempty(temp)   % 06/12/2011: ADDED IF STATEMENT TO CHECK IF THE RECORD IDS EXIST IN THIS MATRIX
                            Rtemp(k,l) = RayMatrix_all(k,temp,1);
                            
                            % Ttemp(k,l) = RayMatrix(k,temp,2);
                            %                         Wtemp(k,l) = exp(-(Dist(l)).^2./(2.*CCPGrid{n,3}.^2));
                            Wtemp(k,l) = 1; % weight of the matrix
                        end
                    end
                end
            end
            RRFAmps = [RRFAmps Rtemp];
            Weights = [Weights Wtemp];
            DepthAxis = depth0;
            ResampleNumber = 50;
            RandomSelectionCount = 0.8; %resample 80% of the original traces
            % apply bootstrap to each bin
            [RRFBootMean,RRFBootSigma,Peaks410,Peaks660,Amps410,Amps660] = ...
                ccp_bootstrap(RRFAmps,Weights,DepthAxis,ResampleNumber,round(size(RRFAmps,2)*RandomSelectionCount));
        else
            RRFBootMean = NaN * ones(length(CCPGrid{i,5}),1);
            RRFBootSigma = NaN * ones(length(CCPGrid{i,5}),1);
        end
        % save the bin
        BootAxis = depth0;
        LatCCP = CCPGrid{i,1};
        LonCCP = CCPGrid{i,2};
        CCPBinLead = ['./matfiles/CCPData', filesep 'Bin_'];
        CCPPath = [CCPBinLead, num2str(sprintf('%0.0f',i)), '.mat'];
        save([CCPPath],'RRFBootMean','RRFBootSigma','BootAxis','LatCCP','LonCCP');
    end
    disp('CCP stacking completes')

    %% rough output
    rf = {};
    nbin  = size(CCPGrid,1);
    lat = [];
    lon = [];
    dep = [];
    k = 0;
    dep0 = 0:dz:zmax;
    
    for n = 1:nbin
        % check if bin exists
        binfile = ['./matfiles/CCPData','/Bin_',num2str(n),'.mat'];
        if exist(binfile,'file')
            k = k + 1;
            load(binfile);
            rf{k} = RRFBootMean;
            lat(k) = LatCCP;
            lon(k) = LonCCP;
        end
        % scatter3(ones(size(rf{k}))*lon(k), ones(size(rf{k}))*lat(k), -dep0,30, rf{k}, 'filled' ); hold on;
    end
    rf1 = cell2mat(rf);
    %% create the interpolation volume
    k = 0;
    V = [];
    for i = 1:length(lat)
        for j = 1:length(dep0)
            k = k + 1;
            V(k,:) = [lon(i) lat(i) dep0(j) rf1(j,i)];
        end
    end
    F = scatteredInterpolant(V(:,1),V(:,2),V(:,3),V(:,4));

    ccpResult.rayMatrix = RayMatrix_all;
    ccpResult.cp = cp;
    ccpResult.rf = F;

end
