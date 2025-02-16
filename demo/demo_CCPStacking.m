% 该脚本主要用于计算CCP

clear; clc; close all;
cd ../
%% 0. 设置路径和参数
% 调用自定义函数 setupPaths() 来添加依赖包和函数的搜索路径
setupPaths();
% 加载配置文件中的各种参数（如数据路径、处理参数等）
config = loadConfig();

% 使用配置文件中的参数
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;
RadonParam         = config.RadonParam;
DeconvParam        = config.DeconvParam;

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

% 根据 config.max_angle_diff 对事件进行方位角一致性的筛选
% idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% 筛选后的事件 ID 列表
eventid = {eventList.evid};
% eventid = eventid(idxConsistentEQ);

% 生成事件-台站对应表，便于后续快速索引
EventStationTable = getEventStationTable(DataStruct);

%% 4. 创建速度模型
% 根据数据和配置文件中的测线长度信息，获取成像剖面的中心和方向等参数
% profileStruct = getImagingProfile(DataStruct, config.profile_length);
ModelParam.LatMin = min(stla)-0.2;
ModelParam.LatMax = max(stla)+0.2;
ModelParam.LonMin = min(stlo)-0.2;
ModelParam.LonMax = max(stlo)+0.2;
ModelParam.npts = 5;

ModelType = '3D';
velocityModel = getVelocityModel(ModelType,ModelParam);

%% 5. 计算接收函数
DataStruct = deconv(DataStruct, DeconvParam);

%% 6. CCP叠加
CCPParam.LatMin = min(stla)-0.2;
CCPParam.LatMax = max(stla)+0.2;
CCPParam.LonMin = min(stlo)-0.2;
CCPParam.LonMax = max(stlo)+0.2;
CCPParam.BinSpacing = 25;
CCPParam.BinSize = 25;
CCPParam.StationCode = {stationList.sta}';

ccpResult = CCPStacking(DataStruct, velocityModel, CCPParam);

%% 7. 结果输出
% plotCCPResults(ccpResult,stlo,stla)

lon1 = 109.9250;
lat1 = 42.3421;
lon2 = 109.9250;
lat2 = 41.1693;
profile1=[lon1,lat1;lon2,lat2];

lon1 = 109.0;
lat1 = 41.7666;
lon2 = 111.5;
lat2 = 41.7666;
profile2=[lon1,lat1;lon2,lat2];
profile = {profile1,profile2};
plotCCPXsection(ccpResult,stlo,stla,profile)

