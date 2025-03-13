%% demo_rankReduction.m
%  本脚本示例了一个完整的地震接收函数计算与 "Rank Reduction" 流程，从读入数据到
%  计算接收函数、再到多事件重建和叠加，以及最终可视化。主要步骤包括：
%    0. 设置路径和参数
%    1. 读入数据
%    2. 预处理
%    3. 计算接收函数
%    4. 对每个事件做 rank reduction (DRR-OTG) 重建
%    5. 叠加接收函数
%    6. 可视化结果

clear; clc; close all;

%% 0. 设置路径和参数
% --------------------------------------------------
try
    setupPaths();  % 添加项目中的工具箱或代码到 MATLAB 路径
catch ME
    warning('setupPaths() not found or failed: %s\nUsing current path instead.', ME.message);
end

% 加载配置
if exist('loadConfig','file') == 2
    config = loadConfig();
else
    error('No loadConfig() found. Please implement or provide config structure.');
end

% 从 config 中提取常用字段
% dataFolder = config.dataFolder;
dataFolder = './data/event_waveforms_BY';
PreprocessingParam = config.PreprocessingParam;
DeconvParam        = config.DeconvParam;
RankReductionParam = config.RankReductionParam;

%% 1. 读入数据
% --------------------------------------------------
fprintf('\n[Step 1] Reading SAC data from folder: %s\n', dataFolder);
DataStruct = read_SAC(dataFolder);
if isempty(DataStruct)
    error('No data read from %s. Check if files exist.', dataFolder);
end

%% 2. 预处理
% --------------------------------------------------
fprintf('\n[Step 2] Preprocessing data...\n');
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. 计算接收函数
% --------------------------------------------------
fprintf('\n[Step 3] Computing receiver functions (deconvolution)...\n');
DeconvParam.verbose = false;
DataStruct = deconv(DataStruct, DeconvParam);

%% 4. 对每个事件做 rank reduction (DRR-OTG) 重建
% --------------------------------------------------
fprintf('\n[Step 4] Doing rank reduction (DRR-OTG) for each event...\n');

dx = 10;
dy = 10;
gridStruct = createGrid(DataStruct, dx, dy);

eventList = getEvents(DataStruct);  % 例如返回 struct array of event info
eventIDs  = {eventList.evid};
DataStruct_drr = [];  % 存储重建后的结果

minStationCount = 50; % 阈值: 若事件台站数小于此数, 则跳过

for iEvent = 1:length(eventIDs)
    evid = eventIDs{iEvent};
    [gather, matchIndex] = getCommonEventGather(DataStruct, evid);

    if length(gather) < minStationCount
        fprintf('Event %s only has %d stations (<%d). Skipped.\n', ...
            evid, length(gather), minStationCount);
        continue;
    end

    % 调用rankReduction对该事件进行3D数据重构
    % gather => [ gather(i).RF.itr ...], gather(i).StationInfo ...
    [gatherReconstructed, d1_otg] = rankReduction(gather, gridStruct, RankReductionParam);

    % 将结果并入 DataStruct_drr
    DataStruct_drr = [DataStruct_drr; gatherReconstructed(:)];
    pause;
    % 若需要调试/可视化 d1_otg，可在 rankReduction 内部或此处实现
    % close all;  % 如果怕产生过多图窗
end

% 将结果转置回与 DataStruct 类似的形状(若需要)
DataStruct_drr = DataStruct_drr';

%% 5. 叠加接收函数
% --------------------------------------------------
fprintf('\n[Step 5] Stacking receiver functions by station...\n');

% stackCommonStationGather() 假设对同一台站的多个事件进行叠加
% 输出：seisout => 叠加后数据, depth0 => 深度轴, mohoStruct => 可选莫霍深度
[seisout, depth0, mohoStruct] = stackCommonStationGather(DataStruct_drr);

%% 6. 可视化
% --------------------------------------------------
fprintf('\n[Step 6] Visualizing some results...\n');

% 6.1 可视化单条记录示例
trace_index = 100;  
if trace_index <= length(DataStruct_drr)
    plotWaveforms(DataStruct_drr, trace_index);
else
    warning('trace_index=%d > data length. Skip plotWaveforms.', trace_index);
end

% 6.2 可视化某一台站的接收函数
stationStruct = getStations(DataStruct_drr);  % eg. returns struct array .sta, .stlo, .stla
stationNames  = {stationStruct.sta};
if ~isempty(stationNames)
    stationToPlot = stationNames{1}; 
    plotCommonStationGather(DataStruct_drr, stationToPlot);
end

% 6.3 批量绘制每个事件的接收函数
eventForPlot = getEvents(DataStruct_drr);
eventIDsPlot = {eventForPlot.evid};

for iEvt = 1:length(eventIDsPlot)
    eID = eventIDsPlot{iEvt};
    gatherEvt = getCommonEventGather(DataStruct_drr, eID);
    if isempty(gatherEvt)
        continue;
    end
    plotCommonEventGather(DataStruct_drr, eID, 'trace','wigb');
    outFig = fullfile('./figures',[eID, '.png']);
    export_fig(outFig, '-r150');
    close all;
end

% 6.4 可视化台站分布
plotStations(DataStruct_drr, 'Baiyanebo_DEM.mat');

% 6.5 可视化事件分布
plotEvents(DataStruct_drr);

fprintf('\nDone! All steps completed.\n');
