function [param] = genVelocityModel(param)

    filename = './velocity_model/Zhao2013_QB_PS2.2.txt';
    x = param.x;
    yv = param.y;
    dz = param.dz;
    nx = param.nx;
    ny = param.ny;
    nz = param.nz;

    slon_ref = param.reflon;
    slat_ref = param.reflat;

    %%-------------------
    matrix = readmatrix(filename);
    lon = matrix(:,1); lat = matrix(:,2); depth = matrix(:,3);
    vp = matrix(:,4); vs = matrix(:,5);

    % velocity model grid
    [xv,yv] = latlon2xy(lon,lat,slon_ref,slat_ref);
    % migrition grid
    [XX, YY, ZZ] = meshgrid(param.x, param.y, param.z);
    
    % 使用散点插值到规则网格
    F_vp = scatteredInterpolant(xv, yv, depth, vp, 'linear', 'none');
    F_vs = scatteredInterpolant(xv, yv, depth, vs, 'linear', 'none');
    
    % 在规则网格点上求值
    Vp = F_vp(XX, YY, ZZ);
    Vs = F_vs(XX, YY, ZZ);
    
    % 转置以匹配期望的维度顺序 [nz,nx,ny]
    Vp = permute(Vp, [3 2 1]);
    Vs = permute(Vs, [3 2 1]);

    % plot
    % figure
    % set(gcf,'Position',[100,100,1200,400],'Color','w')
    % subplot(211)
    % slice(XX,YY,ZZ,permute(Vp,[3 2 1]),[200],[100],[]);
    % set(gca,'ZDir','reverse','FontSize',22);
    % set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on','TickDir','in');
    % shading interp
    % title('Vp model')
    % xlim([0 300]);ylim([0 200]);zlim([0 100])
    % xlabel('X (km)');ylabel('Y (km)');zlabel('Depth (km)');
    % view(-36,25)
    % 
    % subplot(212)
    % slice(XX,YY,ZZ,permute(Vs,[3 2 1]),[200],[100],[]);
    % set(gca,'ZDir','reverse','FontSize',22);
    % set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on','TickDir','in');
    % shading interp
    % title('Vs model')
    % colormap(("jet"))
    % % axis equal
    % xlim([0 300]);ylim([0 200]);zlim([0 100])
    % xlabel('X (km)');ylabel('Y (km)');zlabel('Depth (km)');
    % view(-36,25)
    % export_fig('./velocity-model-QB.png','-r300')


    % load ./matfiles/velocity_model_1D.mat
    % Vp = repmat(v1,1,1,ny);
    % Vs = repmat(v2,1,1,ny);
    %%-------------------
    param.vp = Vp;      % km/s  [nz,nx,ny]
    param.vs = Vs;      % km/s  [nz,nx,ny]

    % figure
    % subplot(131)
    % imagesc(v2);
    % hold on;
    % clim([1.2 4])
    % subplot(132)
    % v3s = Vs(:,:,11);
    % imagesc(v3s);
    % clim([1.2 4])
    % colormap('jet')
    % 
    % subplot(133)
    % plot(v2(:,41),1:nz,'Color','r','LineStyle','-','LineWidth',3)
    % hold on
    % plot(v3s(:,41),1:nz,'Color','b','LineStyle','-','LineWidth',3)
    % set(gca,'YDir','reverse')
    % legend('1D','3D')
    % grid on;
    % xlim([0 6]);
    % ylim([0 100])

end
