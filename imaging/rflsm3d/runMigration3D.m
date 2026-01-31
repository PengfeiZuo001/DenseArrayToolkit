function [mig,dp] = runMigration3D(rfshift,take_off,back_az,src_func,save_wavefield,gridStruct,param)

    % grid
    dx = gridStruct.dx;
    dy = gridStruct.dy;
    dz = gridStruct.dz;
    vp = gridStruct.vp;
    vs = gridStruct.vs;
    
    x = gridStruct.x;
    y = gridStruct.y;
    xo = gridStruct.XInOriginalCoord(:);
    yo = gridStruct.YInOriginalCoord(:);
    

    dt = param.dt;
    nt = param.nt;
    
    src_type = param.src_type;
    fpeak = param.fpeak;
    ispred = param.ispred;

    % for migration

    paramMig.x = gridStruct.x;
    paramMig.y = gridStruct.y;
    paramMig.z = gridStruct.z;

    paramMig.vp = vp;
    paramMig.vs = vs;
    paramMig.dx = dx;
    paramMig.dy = dy;
    paramMig.dz = dz;
    paramMig.flow = param.flow;
    paramMig.fhigh = param.fhigh;
    paramMig.bc = param.bc;
    paramMig.nt = param.nt;
    paramMig.dt = param.dt;

    paramMig.rx = gridStruct.rx(:,1);
    paramMig.ry = gridStruct.ry(:,2);
    paramMig.plot = 0;

    %% 
    dsrc = genPlaneWave(src_func,take_off,back_az,xo,yo,x,y,vp,nt,dt,src_type,fpeak);
    
    disp('==========>')
    disp('do migration')
    mig = ssfm_adj_3D(rfshift,dsrc,save_wavefield,paramMig);

    % predict data
    if ispred
        save_wavefield = 0;  % 1: only to calculate time diff
        [dp,~,~] = ssfm_fd_3D(mig,dsrc,save_wavefield,paramMig);
        S = repmat(any(rfshift),size(rfshift,1),1);

        dp = S.*dp;
    else
        dp = [];
    end


end