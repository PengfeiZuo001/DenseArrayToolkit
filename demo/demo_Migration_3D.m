
clear; clc; close all;

%% 0. 设置路径和参数
cd ../
setupPaths();
config = loadConfig();

% 使用配置文件中的参数
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;
CCPParam           = config.CCPParam;
%% 1. 读入数据
% 读取 dataFolder 中的 SAC 格式地震数据，并封装到 DataStruct 中
DataStruct = read_SAC(dataFolder);

%% 2. 预处理
% 根据配置的预处理参数（如滤波、去均值、切片等）对 DataStruct 进行处理
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. 获取台阵和事件信息
% 从 DataStruct 中提取台站列表和事件列表
stationList = getStations(DataStruct);
eventList   = getEvents(DataStruct);

% 提取台站和事件的经纬度信息，便于后续一致性筛选
stlo = [stationList.stlo]';  % 台站经度
stla = [stationList.stla]';  % 台站纬度
evla = [eventList.evla]';    % 事件震中纬度
evlo = [eventList.evlo]';    % 事件震中经度

% 筛选后的事件 ID 列表
eventid = {eventList.evid};

% 生成事件-台站对应表，便于后续快速索引
EventStationTable = getEventStationTable(DataStruct);

%% 4. 创建速度模型
% 根据数据和配置文件中的测线长度信息，获取成像剖面的中心和方向等参数
% profileStruct = getImagingProfile(DataStruct, config.profile_length);
dx = 10;
dy = 10;
dz = 1;
zmax = 100;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax);
% 创建或获取速度模型，在后续偏移成像中使用
gridStruct = getVelocityModel('3D',gridStruct);

%% 5. 偏移成像
% 准备保存偏移结果的矩阵，这里将所有事件的成像结果进行累积存储
mig = struct();
nMigratedEvents  = 1;   % 用于计数成功完成成像的事件数

%%  
MigParam.paramMig = setMigParam3D(gridStruct);

% 遍历所有符合筛选条件的事件
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    % 提取当前事件对应的地震记录子集（Common Event Gather）
    gather = getCommonEventGather(DataStruct, evid);
    
    % 如果该事件的有效台站数小于 50，则跳过，以保证成像质量
    if length(gather) < 50
        continue
    end
  
    % 对数据进行反褶积处理，以获取接收函数
    DeconvParam.verbose = false;
    gather = deconv(gather, DeconvParam);
    
    % 调用 leastSquaresMig 进行偏移成像
    migResult = leastSquaresMig3D(gather, gridStruct, MigParam);
    
    % 将本次事件的成像结果累积到 dmig 和 dmigls 中
    mig(nMigratedEvents).mig   = migResult.mig;
    mig(nMigratedEvents).migls = migResult.migls;
    
    % 关闭所有图窗，避免在批量处理时生成过多窗口
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end


