function plotStations(DataStruct, demFile)
% PLOTSTATIONS - Using MATLAB Mapping Toolbox to plot station locations
%                on a background DEM (digital elevation model).
%
% Usage:
%   plotStations(DataStruct, demFile)
%
% Inputs:
%   DataStruct - struct array with .StationInfo.stla, .StationInfo.stlo, .StationInfo.sta
%   demFile    - path to a .mat file containing DEM data:
%                demLat, demLon, demZ
%
% Prerequisite: MATLAB Mapping Toolbox
%
% Example:
%   plotStations(DataStruct, 'myDEM.mat')

%% 1. Load DEM data
if nargin < 2 || isempty(demFile)
    demFile = 'Qaidam_DEM.mat';  % 你可以把默认的 DEM 文件名写在这里
end
S = load(demFile);
demLat = S.demLat;  % [M x 1]
demLon = S.demLon;  % [N x 1]
demZ   = S.demZ;    % [M x N], Z(i,j) 与 (demLat(i), demLon(j)) 对应

%% 2. 从 DataStruct 中提取台站信息
lats = [];
lons = [];
staname = {};
for n = 1:length(DataStruct)
    lats(n)    = DataStruct(n).StationInfo.stla;
    lons(n)    = DataStruct(n).StationInfo.stlo;
    staname{n} = DataStruct(n).StationInfo.sta;
end

% 去重（若有同名站）
[staname, idx] = unique(staname);
lats = lats(idx);
lons = lons(idx);

%% 3. 确定地理范围
%   可以根据台站位置设置一个合适的范围，也可以手动指定。
latMin = min(lats) - 0.25;  % 留一点边界
latMax = max(lats) + 0.25;
lonMin = min(lons) - 0.25;
lonMax = max(lons) + 0.25;

% 若数据很小，也要保证 latMin < latMax, lonMin < lonMax
if latMin == latMax, latMax = latMax + 0.5; end
if lonMin == lonMax, lonMax = lonMax + 0.5; end

latlim = [latMin, latMax];
lonlim = [lonMin, lonMax];

%% 4. 创建映射轴并设置范围 (worldmap or axesm)
figure('Position',[10 10 1000 800],'Name','StationMap_MappingToolbox','Color','w');
% worldmap 会根据 latlim/lonlim 自动选取投影并设置坐标
worldmap(latlim, lonlim);
set(gca,'FontSize',10,'FontWeight','normal');

% 或者使用 axesm('MapProjection','mercator'); 然后 setm(ax,'MapLatLimit',latlim,'MapLonLimit',lonlim);

%% 5. 绘制DEM背景
% geoshow 可以将 (demLat, demLon, demZ) 以 'DisplayType','texturemap' 方式显示
% 要保证 demLat 为列向量(升序)或行向量(升序)，demLon 同理
% 如果 lat 降序，需要翻转 demZ；请根据实际情况做检查

geoshow(demLat, demLon, demZ, ...
    'DisplayType','texturemap');

% 给地图设置地形色表 (demcmap) 或自定义
demcmap(demZ);   % 需要 Mapping Toolbox; 若无可用 colormap(parula), colormap('jet'), etc.
colorbar;        % 显示色带
hold on;

%% 6. 绘制台站位置
% 这里使用 geoshow(...) 的 point 模式，也可以使用 plotm, scatterm, textm
geoshow(lats, lons, ...
    'DisplayType','point', ...
    'Marker','^', 'MarkerFaceColor','r','MarkerEdgeColor','k', ...
    'MarkerSize',15);

% 如果想在点旁边标注台站名，可以用 textm:
for i = 1:length(lats)
    textm(lats(i), lons(i), staname{i}, ...
        'VerticalAlignment','bottom','HorizontalAlignment','right', ...
        'FontSize',12,'Color','k','FontWeight','normal');
end

%% 7. 其他修饰
title('Station Map with DEM', 'FontSize',12,'FontWeight','bold');

% gridm, framem, etc. 你可以根据需要再添加
gridm('on');    % 绘制经纬网
framem('on');   % 绘制地图框
mlabel on;      % 标注纬线
plabel on;      % 标注经线

end
