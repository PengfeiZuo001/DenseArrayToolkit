function ccpResult = CCPStacking(DataStruct, velocityModel, CCPParam)

    % 设置CCP数据目录
    ccp_data_directory = './matfiles/CCPData/';
    % 读取参考模型
    [z, rho, vp, vs, qk, qm] = ak135('cont');
    % 设置最大深度和深度间隔
    zmax = 800;
    dz = 0.5;
    % 获取速度模型的vp和vs
    Fvp = velocityModel.vp;
    Fvs = velocityModel.vs;

    % 从CCP参数中获取经纬度范围和网格参数
    LatMin = CCPParam.LatMin;
    LatMax = CCPParam.LatMax;
    LonMin = CCPParam.LonMin;
    LonMax = CCPParam.LonMax;
    BinSpacing = CCPParam.BinSpacing;
    BinSize = CCPParam.BinSize;

    % 获取站点列表
    stalist = CCPParam.StationCode;

    % 调用函数设置网格节点
    CCPGrid = ccp_setup_grid(LatMin, LatMax, LonMin, LonMax, BinSpacing, BinSize);
    
    % 提取网格中心和半径，并绘制圆圈
    centers = [cell2mat(CCPGrid(:,2)) cell2mat(CCPGrid(:,1))];
    radii = km2deg(cell2mat(CCPGrid(:,3)));
    viscircles(centers,radii); hold on;
    % 绘制站点位置
    stationList = getStations(DataStruct);
    stlo=[stationList.stlo];
    stla=[stationList.stla];
    scatter(stlo,stla,100,'^','filled')
    disp('Create CCP bin completes')

    %% 1D射线追踪和时间到深度转换
    disp('Migration starts')
    start_index = 1;

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
        RayDepths = [1*dz:dz:zmax]';
        % correct for heterogeneity
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
        filename=fullfile(ccp_data_directory,['CCPData',stalist{n},'.mat']);
        if exist(filename,'file')
            load(filename,'RayMatrix');
        else
            continue
        end
        % save matrix for individual station into a big matrix
        RayMatrix_all = [RayMatrix_all, RayMatrix];
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
                if ~isempty(Ids)
                    for l = 1:length(Ids)
                        temp = find(RayMatrix_all(k,:,7) == Ids(l));
                        if ~isempty(temp)   % 06/12/2011: ADDED IF STATEMENT TO CHECK IF THE RECORD IDS EXIST IN THIS MATRIX
                            Rtemp(k,l) = RayMatrix_all(k,temp,1);
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

    %% create the interpolation volume
    amp = {};
    nbin  = size(CCPGrid,1);
    lat = [];
    lon = [];
    k = 0;
    for n = 1:nbin
        % check if bin exists
        binfile = ['./matfiles/CCPData','/Bin_',num2str(n),'.mat'];
        if exist(binfile,'file')
            k = k + 1;
            load(binfile);
            amp{k} = RRFBootMean;
            lat(k) = LatCCP;
            lon(k) = LonCCP;
        end
    end
    amp = cell2mat(amp);
    dep = BootAxis;
    k = 0;
    V = [];
    for i = 1:length(lat)
        for j = 1:length(dep)
            k = k + 1;
            V(k,:) = [lon(i) lat(i) dep(j) amp(j,i)];
        end
    end
    F = scatteredInterpolant(V(:,1),V(:,2),V(:,3),V(:,4));

    ccpResult.CCPGrid = CCPGrid;
    ccpResult.rf = F;

end
