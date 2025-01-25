%% 接收函数计算与可视化主函数
% 本脚本示例了从读入数据到计算接收函数，再到可视化叠加结果的典型流程。
% 包括以下步骤：
%   0. 设置路径和参数
%   1. 读入数据
%   2. 预处理
%   3. 计算接收函数
%   4. 叠加接收函数
%   5. 可视化

clear; clc; close all;

%% 0. 设置路径和参数
% 0.1 设置路径
%   这里会调用自定义函数 setupPaths()，它通常负责将项目中所需的工具箱或代码库
%   添加到 MATLAB 的搜索路径中。
setupPaths()

% 0.2 加载配置
%   使用自定义函数 loadConfig() 来加载所需的配置文件，例如数据路径、处理参数等。
config = loadConfig();

% 0.3 使用配置参数
%   这里将配置文件中的重要参数提取为局部变量，以便后续函数使用。
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;        % 若本脚本不需要，可忽略或删除
DeconvParam        = config.DeconvParam;

%% 1. 读入数据
%   使用 read_SAC() 函数从 dataFolder 目录读取 SAC 格式的地震数据，并将
%   数据封装到 DataStruct 结构体中。
DataStruct = read_SAC(dataFolder);

%% 2. 预处理
%   调用 preprocessing()，对数据进行预处理操作：
%   例如去均值、去趋势、带通滤波、剪切时窗等（具体操作由 PreprocessingParam 决定）。
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. 计算接收函数
%   调用 deconv()，通过逆滤波等方法计算地震记录的接收函数。
%   可以在 DeconvParam 中配置不同的去卷积方法、稳态因子、平滑参数等。
DataStruct = deconv(DataStruct, DeconvParam);

%% 4. 叠加接收函数
%   针对同一个台站的所有事件接收函数进行叠加，得到更加稳定和清晰的结果。
%   stackCommonStationGather() 返回叠加后的数据 seisout、深度坐标 depth0，以及
%   可选的莫霍面信息 mohoStruct。
[seisout, depth0, mohoStruct] = stackCommonStationGather(DataStruct);

%% 5. 可视化
%   在完成接收函数计算及叠加后，对结果进行可视化和检查。

% 5.1 可视化单条地震记录示例
trace_index = 100;  % 指定要绘制的记录索引
%   plotWaveforms() 会绘制特定地震记录在时域或频域的波形
plotWaveforms(DataStruct, trace_index);

% 5.2 可视化某一台站的接收函数
station = getStations(DataStruct);
stationList = {station.sta};
station = stationList{1};  % 选取第一个台站来演示
plotCommonStationGather(DataStruct, station);

% 5.3 可视化所有事件的接收函数
%   首先获取事件列表，然后针对每个事件循环绘图并保存为 PNG 文件。
event     = getEvents(DataStruct);
eventList = {event.evid};
for iEvent = 1:length(eventList)
    EventID = eventList{iEvent};
    
    % 提取该事件对应的地震记录
    gather = getCommonEventGather(DataStruct, EventID);
    
    % 使用 plotCommonEventGather() 进行绘图
    %   如果指定 'trace'，则绘制每道波形；也可使用其他方式，如 'wiggle'、'stack' 等
    plotCommonEventGather(DataStruct, EventID, 'trace');
    
    % 将当前图形保存为 PNG 文件，分辨率设为 150 dpi
    export_fig(['./figures/', EventID, '.png'], '-r150');
    
    % 关闭当前所有图窗，防止在批量处理时生成过多窗口
    close all
end

% 5.4 可视化台站和事件在地图上的分布
%   plotStations() 根据 DataStruct 中的台站信息绘制台站位置分布，并将外部地图文件
%   传递给函数绘制背景。
plotStations(DataStruct, 'Baiyanebo_DEM.mat');

%   plotEvents() 用于绘制事件位置分布
plotEvents(DataStruct);