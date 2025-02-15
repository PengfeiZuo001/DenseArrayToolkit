%% DenseArrayToolkit Migration主函数示例
% 该脚本演示了从数据读取到成像结果可视化的基本流程。
% 包含以下主要步骤：
%   0. 设置路径和参数
%   1. 读入数据
%   2. 预处理
%   3. 获取台阵和事件信息
%   4. 创建速度模型
%   5. 偏移成像
%   6. 可视化
%   7. 保存结果

clear; clc; close all;

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
idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% 筛选后的事件 ID 列表
eventid = {eventList.evid};
eventid = eventid(idxConsistentEQ);

% 生成事件-台站对应表，便于后续快速索引
EventStationTable = getEventStationTable(DataStruct);

%% 4. 创建速度模型
% 根据数据和配置文件中的测线长度信息，获取成像剖面的中心和方向等参数
profileStruct = getImagingProfile(DataStruct, config.profile_length);

% 创建或获取速度模型，在后续偏移成像中使用
velocityModel = getVelocityModel();

%% 5. 偏移成像
%   以下注释掉的是其他可能的成像方案示例，可根据实际需要进行调用
%   ccpResult = ccp(DataStruct, velocityModel, CCPParam);
%   hkResult  = hk(DataStruct, HKParam);

% 准备保存偏移结果的矩阵，这里将所有事件的成像结果进行累积存储
dmig   = [];  % 存储常规偏移结果
dmigls = [];  % 存储最小二乘偏移结果
nMigratedEvents  = 1;   % 用于计数成功完成成像的事件数

% 遍历所有符合筛选条件的事件
for iEvent = 5
    evid = eventid{iEvent}; 
    % 提取当前事件对应的地震记录子集（Common Event Gather）
    gather = getCommonEventGather(DataStruct, evid);
    
    % 如果该事件的有效台站数小于 50，则跳过，以保证成像质量
    if length(gather) < 50
        continue
    end
    
    % 根据成像剖面，将台站投影到 profile 上，便于后续空间处理
    gather = project_stations(gather, profileStruct.center, profileStruct.direction, 1);
    
    % 进行射线变换（RadonTransform）以实现数组处理，如围绕某一时距-速度域的滤波等
    gather = radonTransform(gather, RadonParam);
    
    % 设置反褶积参数，此处启用 radonfilter，进一步提升信号质量
    DeconvParam.radonfilter = true;
    
    % 对数据进行反褶积处理，以获取更清晰的波形特征
    gather = deconv(gather, DeconvParam);
    
    % 调用 leastSquaresMig 进行偏移成像，得到常规和最小二乘偏移结果
    migResult = leastSquaresMig(gather, velocityModel, profileStruct, MigParam);
    
    % 将本次事件的成像结果累积到 dmig 和 dmigls 中
    dmig(:, :, nMigratedEvents)   = migResult.mig;
    dmigls(:, :, nMigratedEvents) = migResult.migls;
    
    % 关闭所有图窗，避免在批量处理时生成过多窗口
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end

% 在循环结束后，获取最后一次偏移结果的横坐标 x 和深度坐标 z 用于可视化
x = migResult.x;
z = migResult.z;

%% 6. 可视化
% 将所有事件的成像结果 dmig / dmigls 汇总并进行绘图
plotMigrationResults(dmig, dmigls, x, z);

%% 7. 保存结果
% 将最后一次得到的 migResult（也可以改为需要保存的汇总结果）写入到指定文件中
write_MigResult([config.outputFolder,'/migResult.mat'], migResult);
