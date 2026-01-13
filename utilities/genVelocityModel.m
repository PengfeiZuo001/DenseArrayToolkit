function [gridStruct] = genVelocityModel(gridStruct)

    filename = './velocity_model/Zhao2013_QB_PS2.2.txt';
    x = gridStruct.x;
    yv = gridStruct.y;
    dz = gridStruct.dz;
    nx = gridStruct.nx;
    ny = gridStruct.ny;
    nz = gridStruct.nz;

    slon_ref = gridStruct.originLon;
    slat_ref = gridStruct.originLat;

    %%-------------------
    matrix = readmatrix(filename);
    lon = matrix(:,1); lat = matrix(:,2); depth = matrix(:,3);
    vp = matrix(:,4); vs = matrix(:,5);

    % velocity model grid
%     [xv,yv] = latlon2xy(lon,lat,slon_ref,slat_ref);
    [xv,yv] = latlonToProjectedCoords(lon, lat, gridStruct);
    % migrition grid
    [XX, YY, ZZ] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);
    
    % 使用散点插值到规则网格
    Fvp = scatteredInterpolant(xv, yv, depth, vp, 'linear', 'none');
    Fvs = scatteredInterpolant(xv, yv, depth, vs, 'linear', 'none');
    
    % 在规则网格点上求值
    Vp = Fvp(XX, YY, ZZ);
    Vs = Fvs(XX, YY, ZZ);
    
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
    gridStruct.Fvp = Fvp;
    gridStruct.Fvs = Fvs;
    gridStruct.VP = Vp;      % km/s  [nz,nx,ny]
    gridStruct.VS = Vs;      % km/s  [nz,nx,ny]
    gridStruct.ModelType = '3D';
    %% plot velocity model at 40 km depth
    figure;
    set(gcf,'Position',[0 0 1000 1000],'Color','w')
    hold on;
    
    idx = gridStruct.z == 40;
    V = squeeze(gridStruct.VS(idx,:,:));
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
