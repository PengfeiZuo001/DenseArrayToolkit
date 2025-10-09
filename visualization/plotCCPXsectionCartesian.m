function [distAll,depthAll,VAll] = plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile,dem)
%% create 2D cross-section
% function [distAll,depthAll,VAll] = plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile)

F = scatteredInterpolant(X(:),Y(:),Z(:),V(:));
depth0 = 0:0.5:max(Z(:));
% load colormap
cmap = load('./visualization/colormap/roma.mat');
% load DEM
if ~isempty(dem)
    demX = X(:,:,1);
    demY = Y(:,:,1);
    [nx,ny] = size(demX);
    [LON, LAT] = ProjectedCoordsTolatlon(demX(:),demY(:), gridStruct);
    demZ = interp2(dem.demLon,dem.demLat,dem.demZ,LON,LAT)/1000.;
    demZ = reshape(demZ,nx,ny);
end
figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
ax1 = subplot(1,1,1);
hold on;
h = slice(X,Y,Z,V,[],[],100);
colormap(ax1, flipud(cmap.roma));  % 为波形数据设置不同的色标
cmax = 2*rms(V(:));  % 计算色标的最大值
caxis([-cmax, cmax]);  % 设置波形数据的色标范围
% colorbar;  % 为波形数据显示colorbar


xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(ax1,'ZDir','reverse','FontSize',16)
% % 绘制网格点的位置
% scatter(gridStruct.X(:), gridStruct.Y(:), 10, 'k', 'filled', ...
%     'DisplayName', 'Grid Points');

VAll=[];
distAll = [];
depthAll = [];
totalDist = 0;
elevAll = [];
for n = 1:length(profile)
    x1 = profile{n}(1,1);
    y1 = profile{n}(1,2);
    x2 = profile{n}(2,1);
    y2 = profile{n}(2,2);
    if x1 == x2
        x2 = x1+1e-5;
    end
    l = (y2-y1)/(x2-x1);
    npoints = 100;
    xp = linspace(x1,x2,npoints)';
    yp = y1+(xp-x1)*l;
    dist0 = sqrt((xp-x1).^2 + (yp-y1).^2);
    % plot the geometry of RFs
    Y_profile = repmat(yp',length(depth0),1);
    X_profile = repmat(xp',length(depth0),1);
    Z_profile = repmat(depth0',1,length(xp));
    DIST_profile = repmat(dist0',length(depth0),1);
    Vprofile = F(X_profile,Y_profile,Z_profile);
    Vprofile(isnan(Vprofile)) = 0;

    scatter(xp,yp,50,'ko','filled');
    surface(X_profile,Y_profile,Z_profile,Vprofile,'EdgeColor','none')

    colormap(flipud(cmap.roma))
    cmax = 2*rms(Vprofile(:));
    caxis([-cmax,cmax]);
    
    grid on; box on;
    view(140,20)
    if n > 1
        totalDist = max(distAll(:));
    end
    distAll = [distAll, DIST_profile+totalDist];
    depthAll = [depthAll,Z_profile];
    VAll = [VAll, Vprofile];
    % extract surface elevation from DEM file
%     [lonp, latp] = ProjectedCoordsTolatlon(xp, yp, gridStruct);
%     elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
%     elevAll = [elevAll; elevp];
end
zlim([-50 100])
distElev = distAll(1,:)';

ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...       % 背景透明
           'XTick', [], 'YTick', [], 'ZTick', [], ...
           'View', ax1.View); % 隐藏坐标刻度
if ~isempty(dem)
    demZ = -30*demZ/max(demZ(:));
    hdem = surface(ax2,demX,demY,demZ); hold on;
    % extract station elevation
    rx = gridStruct.rx(:,1);
    ry = gridStruct.ry(:,2);
    rz = interp2(demX,demY,demZ,rx,ry);
    scatter3(ax2,rx,ry,rz,100,'^','MarkerFaceColor','r','MarkerEdgeColor','k');
end
colormap(ax2,'parula')
set(hdem,'EdgeColor','none','FaceAlpha',0.8)
zlim([ax1.ZLim])
xlim([ax1.XLim])
ylim([ax1.YLim])
set(ax2,'ZDir','reverse','Visible','off')
linkaxes([ax1,ax2])