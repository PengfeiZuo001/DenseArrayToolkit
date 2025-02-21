function  gridStruct = createGrid(DataStruct, dx, dy)
% 从数据结构中提取台站信息
stationList = getStations(DataStruct);
stlo = [stationList.stlo]';  % 台站经度
stla = [stationList.stla]';  % 台站纬度
originLon=min(stlo);
originLat=min(stla);
% 将经纬度转换为笛卡尔坐标（相对最小经纬度）
[stationX, stationY] = latlon2xy(stlo, stla, originLon, originLat);
% stationX = stationX-mean(stationX);
% stationY = stationY-mean(stationY);

% % 计算所有相邻台站之间的欧几里得距离
% n = length(stationX);
% distances = zeros(n-1, 1);  % 存储相邻台站之间的距离
% for i = 1:n-1
%     distances(i) = sqrt((stationX(i) - stationX(i+1))^2 + (stationY(i) - stationY(i+1))^2);
% end
% 
% % 计算相邻台站的平均距离
% avgDistance = mean(distances);  % 计算相邻台站的平均距离

% 执行主成分分析（PCA）
coords = [stationX(:), stationY(:)];
[coeff, ~, ~] = pca(coords);  % coeff 包含主成分方向

% 主轴（最大方差方向）和次轴（最小方差方向）
principal_axis = coeff(:,1);  % 主轴
secondaryAxis = coeff(:,2);  % 次轴

% 计算主成分和次成分的范围 
projection_on_principal_axis = coords * principal_axis;
projection_on_secondary_axis = coords * secondaryAxis;

% 台站的投影位置
rx = [projection_on_principal_axis,zeros(size(projection_on_principal_axis))];
ry = [zeros(size(projection_on_secondary_axis)),projection_on_secondary_axis];
% 台站投影位置在原始坐标系中的位置
rxInOriginalCoord = rx * coeff'; 
ryInOriginalCoord = ry * coeff'; 

% for n = 1:length(DataStruct)
%     idx = strcmp({stationList.sta},DataStruct.StationInfo.sta);
%     DataStruct(n).StationInfo.rx = rx(1,)
% end
% 获取台阵的范围
xpad = 5*dx;
ypad = 5*dy;

x_min = min(projection_on_principal_axis)-xpad;
x_max = round(max(projection_on_principal_axis)/dx)*dx+xpad;
y_min = min(projection_on_secondary_axis)-ypad;
y_max = round(max(projection_on_secondary_axis)/dy)*dy+ypad;

% 创建网格
nx = floor((x_max - x_min) / dx) + 1;
ny = floor((y_max - y_min) / dy) + 1;

% 生成网格点的坐标
x = linspace(x_min, x_max, nx);
y = linspace(y_min, y_max, ny);
[X, Y] = meshgrid(x, y);

% 将网格点从主成分坐标系变换回原始坐标系
gridPoints = [X(:), Y(:)];
gridPointsInOriginalCoord = gridPoints * coeff';  % 逆变换
XInOriginalCoord = reshape(gridPointsInOriginalCoord(:,1),size(X));
YInOriginalCoord = reshape(gridPointsInOriginalCoord(:,2),size(Y));

% generate points on principle and secondary axes
principleAxisInOriginalCoord =  x' * coeff(:, 1)';
secondaryAxisInOriginalCoord =  y' * coeff(:, 2)';
% shift principle axis to the center of the grid
shiftVector=[0,mean(y)] * coeff';
xshift=shiftVector(1);
yshift=shiftVector(2);
principleAxisInOriginalCoord = principleAxisInOriginalCoord+[xshift,yshift];
rxInOriginalCoord(:,1)=rxInOriginalCoord(:,1)+xshift;
rxInOriginalCoord(:,2)=rxInOriginalCoord(:,2)+yshift;
% shift rx accordingly
rx(:,2)=rx(:,2)+mean(y);

[tmplon,tmplat]=xy2latlon(principleAxisInOriginalCoord(:,1),...
    principleAxisInOriginalCoord(:,2),originLon,originLat);
principleAxisLatLon = [tmplon,tmplat];
[tmplon,tmplat]=xy2latlon(secondaryAxisInOriginalCoord(:,1),...
    secondaryAxisInOriginalCoord(:,2),originLon,originLat);
secondaryAxisLatLon = [tmplon,tmplat];

% calculate amount of shift in original coordinates
% shiftVector=[mean(x),mean(y)] * coeff';
% xshift=shiftVector(1);
% yshift=shiftVector(2);
% stationX = stationX-xshift;
% stationY = stationY-yshift;
% rxInOriginalCoord(:,1)=rxInOriginalCoord(:,1)-xshift;
% rxInOriginalCoord(:,2)=rxInOriginalCoord(:,2)-yshift;
% ryInOriginalCoord(:,1)=ryInOriginalCoord(:,1)-xshift;
% ryInOriginalCoord(:,2)=ryInOriginalCoord(:,2)-yshift;
% points_on_principle_axis(:,1)=points_on_principle_axis(:,1)-xshift;
% points_on_principle_axis(:,2)=points_on_principle_axis(:,2)-yshift;
% points_on_secondary_axis(:,1)=points_on_secondary_axis(:,1)-xshift;
% points_on_secondary_axis(:,2)=points_on_secondary_axis(:,2)-yshift;
% XInOriginalCoord=XInOriginalCoord-xshift;
% YInOriginalCoord=YInOriginalCoord-yshift;

% 绘制图形
figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
hold on;

% 绘制台站的位置
scatter(stationX, stationY, 'r^', 'filled', 'DisplayName', 'Stations');
scatter(rxInOriginalCoord(:,1),rxInOriginalCoord(:,2),'b^','DisplayName','Projected X location')
scatter(ryInOriginalCoord(:,1),ryInOriginalCoord(:,2),'g^','DisplayName','Projected Y location')

% 绘制主轴和次轴方向
% quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (x_max-x_min)/2*principal_axis(1), (x_max-x_min)/2*principal_axis(2), ...
%     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Principal Axis');
% quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (y_max-y_min)/2*secondary_axis(1), (y_max-y_min)/2*secondary_axis(2), ...
%     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Secondary Axis');
plot(principleAxisInOriginalCoord(:,1),principleAxisInOriginalCoord(:,2),'b','linewidth',2,'DisplayName', 'Principal Axis')
plot(secondaryAxisInOriginalCoord(:,1),secondaryAxisInOriginalCoord(:,2),'g','linewidth',2,'DisplayName', 'Secondary Axis')

% 绘制网格点的位置
scatter(XInOriginalCoord(:), YInOriginalCoord(:), 10, 'k', 'filled', ...
    'DisplayName', 'Grid Points');

% 设置图形
xlabel('X (km)');
ylabel('Y (km)');
legend('show','Location','best');
axis equal;
grid on;
title('Station Positions, PCA Axes, and Grid Points in Original Coordinate System');
hold off;
set(gca,'fontsize',14)

% 返回网格的参数
gridStruct = struct( ...
    'originLon', originLon, ...
    'originLat', originLat, ...
    'originX', 0, ...
    'originY', 0, ...
    'X', X, ...
    'Y', Y, ...
    'XInOriginalCoord',XInOriginalCoord, ...
    'YInOriginalCoord', YInOriginalCoord, ...
    'dx', dx, ...
    'dy', dy, ...
    'nx', nx, ...
    'ny', ny, ...
    'rx', rx, ...
    'ry', ry, ...
    'rxInOriginalCoord',rxInOriginalCoord, ...
    'ryInOriginalCoord',ryInOriginalCoord, ...
    'coeff', coeff, ...
    'principalAxisLatLon', principleAxisLatLon, ...
    'secondaryAxisLatLon', secondaryAxisLatLon, ...
    'principleAxisInOriginalCoord', principleAxisInOriginalCoord, ...
    'secondaryAxisInOriginalCoord', secondaryAxisInOriginalCoord, ...
    'type', 'PCA' ...
);
end
