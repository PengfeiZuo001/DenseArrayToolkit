clear; clc; close all;

% addpath 
addpath ./processRFmatlab-master/
processRFmatlab_startup;
addpath ~/MATLAB/MatSAC;
addpath ./data_io/
addpath ./preprocessing/
addpath ./deconvolution/
addpath(genpath('./utilities/'));
javaaddpath('./utilities/FMI/lib/FMI.jar');
addpath('./utilities/FMI/matTaup');

% 1. 读入数据
dataFolder = './data/event_waveforms_BY';
DataStruct = read_SAC(dataFolder);

% 2. 预处理：带通滤波等
DataStruct = preprocessing(DataStruct);
 
% 3. 拉东变换与结构导向滤波
% DataStruct = radonTransform(DataStruct);
% DataStruct = structureOrientedMedianFilter(DataStruct, ...);

% 4. 计算接收函数 
DeconvParam.gauss = 2.5;
DeconvParam.waterlevel = 0.01;
DeconvParam.itmax = 100;
DeconvParam.minderr = 1e-5;
DeconvParam.phaseshift = 5;
DeconvParam.verbose = false;
DeconvParam.radonfilter = false;
DataStruct = deconv(DataStruct,DeconvParam);

% 5. 获取台阵和事件信息
stationList = getStations(DataStruct);
eventList = getEvents(DataStruct);
EventStationTable = getEventStationTable(DataStruct);

% 6. 偏移成像
% ccpResult = ccp(DataStruct, velocityModel, CCPParam);
% hkResult = hk(DataStruct, HKParam)
% migResult = leastSquaresMig(DataStruct, velocityModel, MigParam);
 
% 7. 可视化
trace_index = 100;
plotWaveforms(DataStruct,trace_index);

station = stationList{1};
plotCommonStationGather(DataStruct,station)

% plot all common event gathers 
for n = 1:length(eventList)
    EventID = eventList{n};
    gather = getCommonEventGather(DataStruct,EventID);
    plotCommonEventGather(DataStruct,EventID,'trace');
    export_fig(['./figures/',EventID,'.png'],'-r150');
    close all
end

% plot event and station distributions
plotStations(DataStruct,'Baiyanebo_DEM.mat');
plotEvents(DataStruct)

% plotMigrationResults(migResult);
 
% 8. 保存结果
% write_MigResult('migResult.mat', migResult);
