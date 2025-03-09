function plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile)
%% create 2D cross-section
F = scatteredInterpolant(X(:),Y(:),Z(:),V(:));
depth0 = 0:0.5:max(Z(:));
% load colormap
cmap = load('./visualization/colormap/roma.mat');
% load DEM
dem = load('./visualization/Qaidam_DEM.mat');
demX = X(:,:,1);
demY = Y(:,:,1);
[nx,ny] = size(demX);
[LON, LAT] = ProjectedCoordsTolatlon(demX(:),demY(:), gridStruct);
demZ = interp2(dem.demLon,dem.demLat,dem.demZ,LON,LAT)/1000.;
demZ = reshape(demZ,nx,ny);

figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
ax1 = subplot(1,1,1);
hold on;
h = slice(X,Y,Z,V,[],[],100);
colormap(ax1, flipud(cmap.roma));  % 为波形数据设置不同的色标
cmax = rms(V(:));  % 计算色标的最大值
caxis([-cmax, cmax]);  % 设置波形数据的色标范围
% colorbar;  % 为波形数据显示colorbar

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(ax1,'ZDir','reverse','YDir','reverse','FontSize',16)
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
    Vprofile(isnan(Vprofile)) = NaN;

    scatter(xp,yp,50,'ko','filled');
    surface(X_profile,Y_profile,Z_profile,Vprofile,'EdgeColor','none')

    colormap(flipud(cmap.roma))
    cmax = rms(Vprofile(:));
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
    [lonp, latp] = ProjectedCoordsTolatlon(xp, yp, gridStruct);
    elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
    elevAll = [elevAll; elevp];
end
zlim([-30 100])
distElev = distAll(1,:)';

ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...       % 背景透明
           'XTick', [], 'YTick', [], 'ZTick', [], ...
           'View', ax1.View); % 隐藏坐标刻度
linkaxes([ax1,ax2])
demZ = -30*demZ/max(demZ(:));
hdem = surface(ax2,demX,demY,demZ); hold on;
% extract station elevation
rx = gridStruct.rx(:,1);
ry = gridStruct.ry(:,2);
rz = interp2(demX,demY,demZ,rx,ry);
scatter3(ax2,rx,ry,rz,100,'^','MarkerFaceColor','r','MarkerEdgeColor','k');
% scatter3(ax2,demX(:),demY(:),-demZ(:),30,demZ(:),'filled')
colormap(ax2,'parula')
set(hdem,'EdgeColor','none','FaceAlpha',0.8)
zlim([ax1.ZLim])
xlim([ax1.XLim])
ylim([ax1.YLim])
set(ax2,'ZDir','reverse','YDir','reverse','Visible','off')

fig = figure();
set(gcf,'Position',[50 50 1200 600],'Color','w')
subplot(511,'Parent',fig)
plot(distElev,elevAll,'k','linewidth',2);
patch([distElev(1);distElev;distElev(end)],[0;elevAll;0],[0.8 0.8 0.8]);
ylabel('Elev. (km)')
xlim([0 max(distAll(1,:))]);
set(gca,'fontsize',14)
ylim([1 6]);
subplot(5,1,2:5,'Parent',fig);
h = pcolor(distAll,depthAll,VAll);
ylim([0 100]);
colormap(flipud(cmap.roma))
cmax = rms(VAll(:));
caxis([-cmax,cmax]);
set(h,'EdgeColor','none')
axis ij
xlabel('Distance (km)')
ylabel('Depth (km)')
set(gca,'fontsize',14)