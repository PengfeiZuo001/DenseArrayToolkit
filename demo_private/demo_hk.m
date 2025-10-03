%% 接收函数HK叠加计算主函数
% 本脚本示例了从读入数据到计算接收函数，再到HK叠加结果的流程。
% 包括以下步骤：
%   0. 设置路径和参数
%   1. 读入数据
%   2. 预处理
%   3. 计算接收函数
%   4. HK叠加
clear; clc; close all;
%% 0. 设置路径和参数
% 0.1 设置路径
%   这里会调用自定义函数 setupPaths()，它通常负责将项目中所需的工具箱或代码库
%   添加到 MATLAB 的搜索路径中。
% addpath ('..')
cd ../
setupPaths()
% 0.2 加载配置
%   使用自定义函数 loadConfig() 来加载所需的配置文件，例如数据路径、处理参数等。

config = loadConfig();

% 0.3 使用配置参数
%   这里将配置文件中的重要参数提取为局部变量，以便后续函数使用。
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
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
%% 4. H-k叠加
%   利用接收函数H-k叠加，分别计算台站下方的地壳厚度和地壳平均Vp/Vs。
%   getCommonStationGather()统计了同一台站的所有接收函数数据。
%   HKstacking()中以前三个台站为例进行计算。
%   HK的计算结果保存在"./HK_bootstrap_tmp.txt",HK叠加结果图和接收函数图
%   保存在'./figures/'中。
HKParam.H = 30:0.1:70;
HKParam.K = 1.6:0.01:2.0;
HKParam.W = [1/3 1/3 1/3];
HKresults=HKstacking(DataStruct,HKParam);