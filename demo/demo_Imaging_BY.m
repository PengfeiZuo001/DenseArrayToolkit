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

% 根据 config.max_angle_diff 对事件进行方位角一致性的筛选
% idxConsistentEQ = filter_earthquakes_by_azimuth(stlo, stla, evlo, evla, config.max_angle_diff);

% 筛选后的事件 ID 列表
eventid = {eventList.evid};
% eventid = eventid(idxConsistentEQ);

% 生成事件-台站对应表，便于后续快速索引
EventStationTable = getEventStationTable(DataStruct);

%% 4. 创建速度模型
% 根据数据和配置文件中的测线长度信息，获取成像剖面的中心和方向等参数
dx = 10;
dy = 10;
dz = 0.5;
zmax = 100;
gridStruct = createGrid(DataStruct, dx, dy, dz, zmax);

npts = 10;
gridStruct = getVelocityModel('3D',gridStruct,npts);

%% 5. 计算接收函数
DataStruct = deconv(DataStruct, DeconvParam);

%% 6. CCP叠加
% 准备保存偏移结果的矩阵，这里将所有事件的成像结果进行累积存储

dimg = [];    % 存储最CCP结果
count = [];
nMigratedEvents  = 1;   % 用于计数成功完成成像的事件数

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
    gather = deconv(gather, DeconvParam);
    
    % 调用rankReduction对该事件进行3D数据重构
    RankReductionParam = config.RankReductionParam;
    RankReductionParam.rank = 10;
    [gatherReconstructed, d1_otg] = rankReduction(gather, RankReductionParam);

    % 调用 CCPCommonEventGather 进行共转换点叠加成像
    ccpResult = CCPCommonEventGather(gatherReconstructed, gridStruct, CCPParam);
    dimg(:,:,:,nMigratedEvents) = ccpResult.img;
    count(:,:,:,nMigratedEvents) = ccpResult.count;
    % 关闭所有图窗，避免在批量处理时生成过多窗口
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end


%% 7. 结果输出
% smooth the image
X = ccpResult.X;
Y = ccpResult.Y;
Z = ccpResult.Z;
V = sum(dimg,4)./max(sum(count,4),1);

figure;
V = mean(dimg,4);
h = slice(X,Y,Z,V,90,90,40);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
colormap(flipud(roma));
cmax = rms(V(:));
caxis([-cmax cmax]);

figure;
fold_map = sum(count,4);
h = slice(X,Y,Z,fold_map,90,90,50);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
cm=colormap('hot');
colormap(flipud(cm));
caxis([0 20])
% save './results/BaiyanEbo_ccp.mat' 'X' 'Y' 'Z' 'V_smooth' 'gridStruct'



