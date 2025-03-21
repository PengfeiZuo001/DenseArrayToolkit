%% DenseArrayToolkit Common Conversion Point Stacking主函数示例
% 该脚本演示了从数据读取到成像结果可视化的基本流程。
% 包含以下主要步骤：
%   0. 设置路径和参数
%   1. 读入数据
%   2. 预处理
%   3. 获取台阵和事件信息
%   4. 创建速度模型
%   5. CCP成像
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
CCPParam           = config.CCPParam;

dataFolder = './data/event_waveforms_QBII';
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

% 事件 ID 列表
eventid = {eventList.evid};

% 生成事件-台站对应表，便于后续快速索引
EventStationTable = getEventStationTable(DataStruct);

%% 4. 创建速度模型
% 根据台站位置，进行主成分分析，自动创建成像网格
dx = 5;
dy = 5;
gridStruct = createGrid(DataStruct, dx, dy);
% 创建或获取速度模型，在后续成像中使用
gridStruct = getVelocityModel('3D',gridStruct,5);
%% 5. 偏移成像
% 准备保存偏移结果的矩阵，这里将所有事件的成像结果进行累积存储
dimg = [];    % 存储最CCP结果
count = [];
nMigratedEvents  = 1;   % 用于计数成功完成成像的事件数

% 遍历所有符合筛选条件的事件
for iEvent = 1:length(eventid)
    evid = eventid{iEvent}; 
    % 提取当前事件对应的地震记录子集（Common Event Gather）
    gather = getCommonEventGather(DataStruct, evid);
    
    % 如果该事件的有效台站数小于 60，则跳过，以保证成像质量
    if length(gather) < 60
        continue
    end
  
    % 进行拉东变换（RadonTransform）以实现台阵处理
%     gather = radonTransform(gather, RadonParam);
    
    % 设置反褶积参数，此处启用 radonfilter，进一步提升信号质量
    DeconvParam.radonfilter = false;
    DeconvParam.verbose = false;
    % 对数据进行反褶积处理，以获取接收函数
    gather = deconv(gather, DeconvParam);
    
    % 调用 CCPCommonEventGather 进行共转换点叠加成像
    ccpResult = CCPCommonEventGather(gather, gridStruct, CCPParam);
    dimg(:,:,:,nMigratedEvents) = ccpResult.img;
    count(:,:,:,nMigratedEvents) = ccpResult.count;

    % 关闭所有图窗，避免在批量处理时生成过多窗口
    close all;
    
    nMigratedEvents = nMigratedEvents + 1;
end
%% 6. 可视化
% smooth the image
X = ccpResult.X;
Y = ccpResult.Y;
Z = ccpResult.Z;
V = sum(dimg,4)./max(sum(count,4),1);

% 绘制切片
figure;
set(gcf,'Position',[50 50 800 800],'Color','w')
h = slice(X,Y,Z,V,90,90,40);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
colormap(flipud(roma));
cmax = rms(V(:));
caxis([-cmax cmax]);

% plot user defined profiles
x1 = 0;
y1 = 0;
x2 = 30;
y2 = 90;
x3 = 90;
y3 = 100;
x4 = 130;
y4 = 160;
x5 = 310;
y5 = 30;
profile1=[x1,y1;x2,y2];
profile2=[x2,y2;x3,y3];
profile3=[x3,y3;x4,y4];
profile4=[x4,y4;x5,y5];
profile = {profile1,profile2,profile3,profile4};
plotCCPXsectionCartesian(X,Y,Z,V_smooth,gridStruct,profile)
%% 7. 保存结果
% 将ccpResult写入到指定文件中
write_MigResult([config.outputFolder,'/ccpResult.mat'], ccpResult);