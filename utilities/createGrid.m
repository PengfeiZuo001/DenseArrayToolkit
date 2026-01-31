function  gridStruct = createGrid(DataStruct, dx, dy, dz, zmax, xpad, ypad)
if nargin < 4
    dz = 1;
    zmax = 100;
    xpad = dx;
    ypad = dy;
end

stationList = getStations(DataStruct);
stlo = [stationList.stlo]'; 
stla = [stationList.stla]';  
originLon=min(stlo);
originLat=min(stla);

[stationX, stationY] = latlon2xy(stlo, stla, originLon, originLat);

if nargin < 2

    n = length(stationX);
    distances = zeros(n-1, 1); 
    for i = 1:n-1
        distances(i) = sqrt((stationX(i) - stationX(i+1))^2 + (stationY(i) - stationY(i+1))^2);
    end
    %
    avgDistance = mean(distances);  
    dx = avgDistance;
    dy = avgDistance;
end

coords = [stationX(:), stationY(:)];
[coeff, ~, ~] = pca(coords);  % coeff 

threshold = cosd(5); 
v1 = coeff(:,1);     
v2 = coeff(:,2);   

dot_east = v1(1);
dot_north_v1 = v1(2);

if abs(dot_east) > abs(dot_north_v1) 
    if abs(dot_east) >= threshold
        new_v1 = sign(dot_east)*[1; 0];
        
        new_v2 = sign(v2(2))*[0; 1];
        
        coeff_adjusted = [new_v1, new_v2];
        if det(coeff_adjusted) < 0
            coeff_adjusted(:,2) = -coeff_adjusted(:,2);
        end
        coeff = coeff_adjusted;
    end
else 
    if abs(dot_north_v1) >= threshold
        new_v1 = sign(dot_north_v1)*[0; 1];
        new_v2 = sign(v2(1))*[1; 0];
        
        coeff_adjusted = [new_v1, new_v2];
        
        if det(coeff_adjusted) < 0
            coeff_adjusted(:,2) = -coeff_adjusted(:,2);
        end
        coeff = coeff_adjusted;
    end
end

principal_axis = coeff(:,1);  % main axis
secondaryAxis = coeff(:,2);  % secondary axis

projection_on_principal_axis = coords * principal_axis;
projection_on_secondary_axis = coords * secondaryAxis;

rx = [projection_on_principal_axis,zeros(size(projection_on_principal_axis))];
ry = [zeros(size(projection_on_secondary_axis)),projection_on_secondary_axis];
rxInOriginalCoord = rx * coeff'; 
ryInOriginalCoord = ry * coeff'; 

% for n = 1:length(DataStruct)
%     idx = strcmp({stationList.sta},DataStruct.StationInfo.sta);
%     DataStruct(n).StationInfo.rx = rx(1,)
% end
% xpad = 1*dx;
% ypad = 1*dy;

x_min = min(projection_on_principal_axis)-xpad;
x_max = round(max(projection_on_principal_axis)/dx)*dx+xpad;
y_min = min(projection_on_secondary_axis)-ypad;
y_max = round(max(projection_on_secondary_axis)/dy)*dy+ypad;

nx = floor((x_max - x_min) / dx) + 1;
ny = floor((y_max - y_min) / dy) + 1;

x = linspace(x_min, x_max, nx);
y = linspace(y_min, y_max, ny);
[X, Y] = meshgrid(x, y);

gridPoints = [X(:), Y(:)];
gridPointsInOriginalCoord = gridPoints * coeff';  
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

figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
hold on;

scatter(stationX, stationY, 'r^', 'filled', 'DisplayName', 'Stations');
scatter(rxInOriginalCoord(:,1),rxInOriginalCoord(:,2),'b^','DisplayName','Projected X location')
scatter(ryInOriginalCoord(:,1),ryInOriginalCoord(:,2),'g^','DisplayName','Projected Y location')


% quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (x_max-x_min)/2*principal_axis(1), (x_max-x_min)/2*principal_axis(2), ...
%     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Principal Axis');
% quiver(mean(gridPointsInOriginalCoord(:,1)), mean(gridPointsInOriginalCoord(:,2)), (y_max-y_min)/2*secondary_axis(1), (y_max-y_min)/2*secondary_axis(2), ...
%     'MaxHeadSize', 2, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Secondary Axis');
plot(principleAxisInOriginalCoord(:,1),principleAxisInOriginalCoord(:,2),'b','linewidth',2,'DisplayName', 'Principal Axis')
plot(secondaryAxisInOriginalCoord(:,1),secondaryAxisInOriginalCoord(:,2),'g','linewidth',2,'DisplayName', 'Secondary Axis')

scatter(XInOriginalCoord(:), YInOriginalCoord(:), 10, 'k', 'filled', ...
    'DisplayName', 'Grid Points');

xlabel('X (km)');
ylabel('Y (km)');
legend('show','Location','best');
axis equal;
grid on;
title('Station Positions, PCA Axes, and Grid Points in Original Coordinate System');
hold off;
set(gca,'fontsize',14)
% construct depth axis
z = 0:dz:zmax;
nz = length(z);
% calculate the area limit in lat,lon
[LonMin,LatMin] = xy2latlon(min(XInOriginalCoord(:)),min(YInOriginalCoord(:)),originLon,originLat);
[LonMax,LatMax] = xy2latlon(max(XInOriginalCoord(:)),max(YInOriginalCoord(:)),originLon,originLat);

% return
gridStruct = struct( ...
    'stationLon', stlo, ...
    'stationLat', stla, ...
    'stationX', stationX, ...
    'stationY', stationY, ...
    'originLon', originLon, ...
    'originLat', originLat, ...
    'originX', 0, ...
    'originY', 0, ...
    'LatMin', LatMin, ...
    'LonMin', LonMin, ...
    'LatMax', LatMax, ...
    'LonMax', LonMax, ...
    'X', X, ...
    'Y', Y, ...
    'x', x, ...
    'y', y, ...
    'z', z, ...
    'XInOriginalCoord', XInOriginalCoord, ...
    'YInOriginalCoord', YInOriginalCoord, ...
    'dx', dx, ...
    'dy', dy, ...
    'dz', dz, ...
    'nx', nx, ...
    'ny', ny, ...
    'nz', nz, ...
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
