function [lsmig,pre_rf] = runLSM3D(rfshift,take_off,back_az,src_func,save_wavefield,gridStruct,param)

    % grid
    dx = gridStruct.dx;
    dy = gridStruct.dy;
    dz = gridStruct.dz;
    nx = gridStruct.nx;
    ny = gridStruct.ny;
    nz = gridStruct.nz;

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

    % for cg
    tol = param.tol;
    mu = param.mu;
    itermax = param.itermax;

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
    
    %%  dot product test
    % S=repmat(squeeze(any(dsrc)),1,1,size(dsrc,1));
    % S=permute(S,[3,1,2]);
    % S=S(:);
    % itermax = 10;
    % if_cg = 1;
    % L = @(m) S.*ssfm_fd_3D(m,dsrc,save_wavefield,param,if_cg);
    % Lt= @(d) ssfm_adj_3D(S.*d,dsrc,save_wavefield,param,if_cg);  
    % m1=randn(nz,nx,ny);
    % d2=randn(nt,nx,ny);
    % m1=m1(:);
    % d2=d2(:);
    % d1=L(m1);
    % m2= Lt(d2);
    % dot1=sum(sum(sum(d1.*d2)))
    % dot2=sum(sum(sum(m1.*m2)))

%% ---------------------------------------------------------------
    disp('Preconditioned least-squares migration begins')
    din = rfshift(:);
    S = repmat(any(rfshift),size(rfshift,1),1);
    S = S(:);
    if_cg = 1;
    P = @(u) precon_x_3d(u,nx,ny,nz,dx,dy,dz);
    Pt= @(d) preconT_x_3d(d,nx,ny,nz,dx,dy,dz);

    L = @(m) S.*ssfm_fd_3D(m,dsrc,save_wavefield,paramMig,if_cg);
    
    Lt = @(d) ssfm_adj_3D(S.*d,dsrc,save_wavefield,paramMig,if_cg);    
    LP = @(u) L(P(u));
    PtLt = @(d) Pt(Lt(d));

    %-----------------------------
    %% dot product test
    % m1=randn(nz,nx,ny);
    % d2=randn(nt,nx,ny);
    % m1=m1(:);
    % d2=d2(:);
    % 
    % d1=LP(m1);
    % m2= PtLt(d2);
    % 
    % dot1=sum(sum(sum(d1.*d2)));
    % dot2=sum(sum(sum(m1.*m2)));
    % 
    % fprintf('dot1: %3f \n',dot1);
    % fprintf('dot2: %3f \n',dot2);
    % if abs(dot1 - dot2) < 10^-5
    %     disp('pass the dot product  test! ')
    % else
    %     disp('try again...')
    % end
    %--------------------------------

    A = @(u) PtLt(LP(u)) + mu * u;
    b = PtLt(din);

    tic;
    disp('PCG begins...')
    [utmp1,flag,relres,iter,resvec] = pcg(A,b,tol,itermax);
    disp(['=====> Iteration: ,',num2str(itermax),',flag: ',num2str(flag),', min value:',num2str(min(resvec))])
    toc;

    mtmp1 = P(utmp1(:));
    lsmig = reshape(mtmp1,nz,nx,ny);

    % predict data
    if ispred
        Hd = L(lsmig);
        pre_rf = reshape(Hd/max(Hd(:)),nt,nx,ny);
    else
        pre_rf = [];
    end

end