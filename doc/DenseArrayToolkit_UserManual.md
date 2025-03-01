# DenseArrayToolkit User Manual

<!-- TOC -->
- [DenseArrayToolkit User Manual](#densearraytoolkit-user-manual)
  - [1. Introduction](#1-introduction)
  - [2. Installation \& Requirements](#2-installation--requirements)
  - [3. Quick Start](#3-quick-start)
  - [4. Workflow \& Modules](#4-workflow--modules)
    - [4.1 Data Reading](#41-data-reading)
    - [4.2 Preprocessing](#42-preprocessing)
    - [4.3 Receiver Function Calculation](#43-receiver-function-calculation)
    - [4.4 Rank Reduction / Other Seismic Processing](#44-rank-reduction--other-seismic-processing)
    - [4.5 Stacking \& Visualization](#45-stacking--visualization)
  - [5. GUI Usage](#5-gui-usage)
    - [5.1 Main Interface Overview](#51-main-interface-overview)
    - [5.2 File Menu](#52-file-menu)
    - [5.3 Tools Menu (Rank Reduction)](#53-tools-menu-rank-reduction)
    - [5.4 Common Operations](#54-common-operations)
  - [6. Advanced Features \& Customization](#6-advanced-features--customization)
  - [7. Troubleshooting](#7-troubleshooting)
  - [8. Appendices](#8-appendices)
    - [8.1 Data Structure Reference](#81-data-structure-reference)
    - [8.2 Parameter Lists](#82-parameter-lists)
    - [8.3 Example Scripts](#83-example-scripts)
  - [9. Revision History](#9-revision-history)
  - [10. Authors](#10-authors)
  - [11. License](#11-license)
  - [12. Contributing](#12-contributing)
  - [13. FAQ](#13-faq)
  - [14. Contact Information](#14-contact-information)
<!-- /TOC -->

---

## 1. Introduction

**DenseArrayToolkit** 是一套面向地震台阵数据处理与成像的 MATLAB 工具包，提供了多种常见的数据读取、预处理、接收函数计算以及台阵处理和成像方法（例如 **Rank Reduction**、**Radon Transform**、**CCP**、**Migration Imaging** 等）功能。

- **主要功能：**
  - 快速加载 SAC 格式数据
  - 多种预处理手段（滤波、去均值、去趋势等）
  - 接收函数计算与叠加
  - 台阵方法：**Rank Reduction (DRR)**、**Radon Transform**、**Migration** 等
  - 友好的 **GUI** 支持常用流程可视化与操作

- **使用人群：**
  - 地震学研究人员、地球物理专业学生
  - 需要针对大规模台阵数据进行半自动/自动处理的工程师

---

## 2. Installation & Requirements

1. **MATLAB 版本要求**  
   - 推荐使用 R2020b 或更高版本  
   - 如果要使用 **GeographicAxes** / **UIAxes** 的地图功能，需要 R2020b+

2. **依赖的工具箱：**  
   - *Mapping Toolbox*（若使用与地理投影相关的绘图函数）
   - *Signal Processing Toolbox*

3. **安装方法：**  
   - 将 `DenseArrayToolkit` 文件夹放置在 MATLAB 搜索路径，或在 MATLAB 命令行执行:
     ```matlab
     addpath('path_to_DenseArrayToolkit');
     savepath;
     ```
   - 对于 GUI，仅需保证 `.mlapp` 文件也在搜索路径下或同一文件夹中。

4. **外部依赖**
   - 若某些脚本调用第三方库，需要说明如何获取与安装。

---

## 3. Quick Start

以最简要的方式示范从**读入数据**到**可视化**结果的流程，见demo_Stacking.m：

1. **Setup Paths**
   ```matlab
   setupPaths();  % 添加 DenseArrayToolkit 到 MATLAB path

2. **Load Config**
   ```matlab
   config = loadConfig();  % 包含 dataFolder, PreprocessingParam
3. **Read Data**
    ```matlab
    DataStruct = read_SAC(config.dataFolder); %读入SAC数据
   
4. **Preprocess**
   ```matlab
   DataStruct = preprocessing(DataStruct, config.PreprocessingParam); % 数据预处理

5. **Deconvolution**
   ```matlab
   DataStruct = deconv(DataStruct, config.DeconvParam); % 计算接收函数
6. **Stack & Visualization**
   ```matlab
   [seisout, depth0] = stackCommonStationGather(DataStruct); % 叠加接收函数
   plotCommonStationGather(DataStruct, 'STATION_ID'); % 绘制台站道集
若想使用 GUI，在 MATLAB 命令行执行    
   ```matlab
   DenseArrayToolkit_GUI
   ```
即可打开主界面并按照相应菜单完成同样的处理流程。


## 4. Workflow & Modules
### 4.1 Data Reading
`read_SAC(dataFolder)`:

从指定目录读取 SAC 文件；返回 DataStruct。
DataStruct 结构示例：
DataStruct(i).TraceData: 波形数组
DataStruct(i).StationInfo: 台站信息
DataStruct(i).EventInfo: 事件信息
常见错误：SAC 文件头缺失、文件格式不兼容等。
### 4.2 Preprocessing
`preprocessing(DataStruct, PreprocessingParam)`:

可配置去均值、去趋势、带通滤波、重采样等
常用参数：
参数名	默认值	说明
tstart, tend	-60,600	读波形相对事件起始时间窗
lows, highs	0.1,2.0	带通滤波截止频率 (Hz)
resample_period	0.1	重采样周期 (s)
### 4.3 Receiver Function Calculation
`deconv(DataStruct, DeconvParam)`:

通过逆滤波(去卷积)计算接收函数
常用参数包括 gauss、waterlevel、itmax、minderr 等
输出仍是 DataStruct 形式，新添加 .RF 字段包含了接收函数结果
### 4.4 Rank Reduction / Other Seismic Processing
`rankReduction(gather, RankReductionParam)`:

实现矩阵降秩(去噪+插值)；需要 latlon2xy、drr3drecon_otg 等内部函数
示例参数： nx, ny, rank, K, niter 等
如果只想处理指定事件，可先用 getCommonEventGather(DataStruct, evid)
radonTransform, migrationImaging 等同理。
### 4.5 Stacking & Visualization
- `stackCommonStationGather(DataStruct)`:
将同一台站的多事件接收函数叠加返回seisout（叠加后波形/数据）与深度轴depth0
常用可视化函数：
- `plotCommonStationGather(DataStruct, 'STATION_ID')`: 绘制同一台站不同事件的波形
- `plotCommonEventGather(DataStruct, 'EVID')`: 绘制同一事件不同台站的波形
- `plotStations(DataStruct, demFile)`, `plotEvents(DataStruct)`: 显示台站/事件在地图上的分
## 5. GUI Usage
### 5.1 Main Interface Overview
当您在 MATLAB 命令行执行:

matlab
Copy
Edit
DenseArrayToolkit_GUI
将出现主界面，包括顶栏菜单和一系列按钮：

File: “Load Config”, “Load Data”, “Save Results”
Compute: “Preprocess”, “Deconvolution”, “Stack” 等
Tools: “Rank Reduction”, “Radon Transform”, “Migration” (可选)
View: 用于可视化操作（波形查看、事件/台站分布）
### 5.2 File Menu
Load Config: 读取 config，其中包含 dataFolder, PreprocessingParam, ...
Load Data: 调用 read_SAC，将结果保存到 app.DataStruct
Save Results: 弹出对话框，保存处理后的 DataStruct 或 DataStructDRR 等
### 5.3 Tools Menu
Rank Reduction：
打开 RankReductionApp，可编辑相关参数（rank, nx, ny, etc.）
执行对当前 DataStruct 中某些或全部事件做降秩去噪
处理结果保存在 DataStructDRR
Radon Transform:
打开XXX
### 5.4 Common Operations
- 读取数据：File -> Load Data
- 进行预处理：Compute -> Preprocess
- 接收函数：Compute -> Deconvolution
- 台阵降噪：Tools -> Rank Reduction -> 选定参数并 Run
- 查看处理结果：View -> PlotStation / PlotEvent / PlotWaveforms
## 6. Advanced Features & Customization
自定义去卷积
修改 deconv.m 中的核函数或在 DeconvParam 里切换不同算法
添加新菜单或功能
在 GUI 中可在 createComponents(app) 里添加 uimenu，并在回调中注入新功能
修改默认配置
在 loadConfig.m 中更改各模块的默认参数，如 PreprocessingParam、RankReductionParam 等
## 7. Troubleshooting
无数据读入

检查 config.dataFolder 路径是否正确
SAC 文件是否有 .SAC 后缀并具备合法头段
Mapping Toolbox 错误

如果 plotStations 或地理投影函数报错，确保安装了 Mapping Toolbox；或使用 GeographicAxes 进行替代
内存不足

当处理超大量台阵数据时，可能出现内存压力；可尝试 -v7.3 保存分段数据或减少一次性加载
GUI 回调中报错

大多数都可在命令行查看 ME.message；也可在 app.statusLabel 查看
## 8. Appendices
### 8.1 Data Structure Reference
DataStruct(i) 常见字段：

.TraceData: 波形数组 [Nt x 1]
.StationInfo:
.stla, .stlo, .sta
.EventInfo:
.evla, .evlo, .evdp, .mag
.RF (若已计算接收函数)
.itr, .ittime, ...
### 8.2 Parameter Lists
| 组件                   | 参数名             | 默认值  | 说明                           |
|------------------------|--------------------|---------|---------------------------------|
| **PreprocessingParam** | `tstart`, `tend`  | -60,600 | 读波形相对时间窗 (sec)          |
|                        | `lows`, `highs`   | 0.1,2.0 | 带通滤波边界 (Hz)               |
|                        | `resample_period` | 0.1     | 重采样周期 (sec)               |
| **DeconvParam**        | `gauss`           | 5       | Gauss 参数                      |
|                        | `waterlevel`      | 0.01    | 去卷积水位                      |
|                        | `itmax`           | 100     | 迭代最大次数                    |
|                        | `minderr`         | 1e-5    | 最小误差                        |
| **RankReductionParam** | `nx`, `ny`        | 11,14   | 网格点数 (X, Y方向)            |
|                        | `rank`            | 5       | DRR 矩阵秩                     |
|                        | `K`               | 4       | OTG 算法用额外参数             |
|                        | `mode`            | 1       | DRR 算法模式 (1=默认)          |
|                        | `niter`           | 20      | 最大迭代次数                   |
|                        | `eps`             | 0.001   | 收敛阈值/正则化系数           |

### 8.3 Example Scripts
demo_main.m: 展示最基础的数据读取与处理流程
demo_rankReduction.m: 展示如何对单事件或多事件做降秩去噪
demo_radon.m: 演示使用 RadonParam 进行数据透射变换
## 9. Revision History
v1.0 (2025-01-30):
初始发布，包含基础预处理与 Rank Reduction GUI

v1.1 (2025-02-15):
新增 Radon Transform 方法
优化 GUI 交互（新增保存功能）

v1.2 (2025-03-01):
修复 SAC 读取若干Bug，完善配置文件结构

## 10. Authors
- Author 1: Name, Affiliation, Email
- Author 2: Name, Affiliation, Email
- Author 3: Name, Affiliation, Email

© 2025 DenseArrayToolkit Authors. 如果您在使用中遇到问题或有功能改进需求，欢迎与我们联系或在 GitHub Issues 中提交反馈。

## 11. License
本软件遵循 MIT 许可证。详细信息请参见 LICENSE 文件。

## 12. Contributing
欢迎贡献代码和报告问题！请参阅 CONTRIBUTING.md 了解详细的贡献指南。

## 13. FAQ
- **Q**: 如何安装 DenseArrayToolkit？
  **A**: 请参阅第 2 章 [安装与要求](#2-installation--requirements)。

- **Q**: 如何报告一个 Bug？
  **A**: 请在 GitHub 上提交一个 Issue，并提供详细的错误信息和复现步骤。

## 14. Contact Information
如果您有任何问题或建议，请通过以下方式联系我们：
- Email: support@densearraytoolkit.com
- GitHub: [DenseArrayToolkit](https://github.com/DenseArrayToolkit)