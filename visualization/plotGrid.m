function plotGrid(gridStruct)
% 绘制图形
figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
hold on;

% 绘制台站的位置
scatter(gridStruct.stationX, gridStruct.stationY, 'r^', 'filled', 'DisplayName', 'Stations');
scatter(gridStruct.rxInOriginalCoord(:,1),gridStruct.rxInOriginalCoord(:,2),'b^','DisplayName','Projected X location')
scatter(gridStruct.ryInOriginalCoord(:,1),gridStruct.ryInOriginalCoord(:,2),'g^','DisplayName','Projected Y location')

% 绘制主轴和次轴方向
plot(gridStruct.principleAxisInOriginalCoord(:,1),gridStruct.principleAxisInOriginalCoord(:,2),'b','linewidth',2,'DisplayName', 'Principal Axis')
plot(gridStruct.secondaryAxisInOriginalCoord(:,1),gridStruct.secondaryAxisInOriginalCoord(:,2),'g','linewidth',2,'DisplayName', 'Secondary Axis')

% 绘制网格点的位置
scatter(gridStruct.XInOriginalCoord(:), gridStruct.YInOriginalCoord(:), 10, 'k', 'filled', ...
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