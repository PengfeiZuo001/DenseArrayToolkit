function plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile)
%% create 2D cross-section
F = scatteredInterpolant(X(:),Y(:),Z(:),V(:));
depth0 = 0:0.5:200;
% load colormap
cmap = load('./visualization/colormap/roma.mat');
% load DEM
dem = load('./visualization/Qaidam_DEM.mat');

figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
hold on;
h = slice(X,Y,Z,V,[],[],40);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
% % 绘制台站的位置
% scatter(gridStruct.stationX, gridStruct.stationY, 'r^', 'filled', 'DisplayName', 'Stations');
% scatter(gridStruct.rxInOriginalCoord(:,1),gridStruct.rxInOriginalCoord(:,2),'b^','DisplayName','Projected X location')
% scatter(gridStruct.ryInOriginalCoord(:,1),gridStruct.ryInOriginalCoord(:,2),'g^','DisplayName','Projected Y location')
% 
% % 绘制主轴和次轴方向
% % quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (x_max-x_min)/2*principal_axis(1), (x_max-x_min)/2*principal_axis(2), ...
% %     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Principal Axis');
% % quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (y_max-y_min)/2*secondary_axis(1), (y_max-y_min)/2*secondary_axis(2), ...
% %     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Secondary Axis');
% plot(gridStruct.principleAxisInOriginalCoord(:,1),gridStruct.principleAxisInOriginalCoord(:,2),'b','linewidth',2,'DisplayName', 'Principal Axis')
% plot(gridStruct.secondaryAxisInOriginalCoord(:,1),gridStruct.secondaryAxisInOriginalCoord(:,2),'g','linewidth',2,'DisplayName', 'Secondary Axis')
% 
% % 绘制网格点的位置
% scatter(gridStruct.XInOriginalCoord(:), gridStruct.YInOriginalCoord(:), 10, 'k', 'filled', ...
%     'DisplayName', 'Grid Points');
% 
% % 设置图形
% xlabel('X (km)');
% ylabel('Y (km)');
% legend('show','Location','best');
% axis equal;
% grid on;
% title('Station Positions, PCA Axes, and Grid Points in Original Coordinate System');
% hold off;
% set(gca,'fontsize',14)
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
    xlabel(['Longitude' ,char(176)])
    ylabel(['Latitude' ,char(176)])
    zlabel('Depth (km)')
    set(gca,'FontSize',16)

    surface(X_profile,Y_profile,Z_profile,Vprofile,'EdgeColor','none')
    zlim([0 100])
    colormap(flipud(cmap.roma))
    caxis([-0.01 0.01])
    colorbar
    grid on; box on;
    view(140,20)
    if n > 1
        totalDist = max(distAll(:));
    end
    distAll = [distAll, DIST_profile+totalDist];
    depthAll = [depthAll,Z_profile];
    VAll = [VAll, Vprofile];
    % saveas(gcf, 'example_ccp.png','png')
    % extract surface elevation from DEM file
    [lonp, latp] = ProjectedCoordsTolatlon(xp, yp, gridStruct);
    elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
    elevAll = [elevAll; elevp];
end
distElev = distAll(1,:)';
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
caxis([-0.005,0.005]);
set(h,'EdgeColor','none')
axis ij
xlabel('Distance (km)')
ylabel('Depth (km)')
set(gca,'fontsize',14)