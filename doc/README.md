<div style="display: flex; justify-content: center; align-items: center; height: 100vh;">
  <div align="center">
    <h1>DenseArrayToolkit User Manual</h1>
    <p><strong>Zhejiang University, Hangzhou, China</strong></p>
    <p><strong>Version 0.1</strong></p>
    <p><strong>March 1, 2025</strong></p>
  </div>
</div>

<div style="page-break-after: always;"></div>

## Contents

<!-- TOC -->
  - [1. Introduction](#1-introduction)
  - [2. Installation & Requirements](#2-installation--requirements)
  - [3. Quick Start](#3-quick-start)
  - [4. Modules/Functions](#4-modulesfunctions)
    - [4.1 `read_SAC()`](#41-read_sac)
    - [4.2 `preprocessing()`](#42-preprocessing)
    - [4.3 `deconv()`](#43-deconv)
    - [4.4 `radonTransform()`](#44-radontransform)
    - [4.5 `rankReduction()`](#45-rankreduction)
    - [4.6 `CCPStacking()`](#46-ccpstacking)
    - [4.7 `leastSquaresMig()`](#47-leastsquaresmig)
    - [4.8 `leastSquaresMig3D()`](#48-leastsquaresmig3d)
    - [4.9 `HKstacking()`](#49-hkstacking)
    - [4.10 `plotCommonEventGather()`](#410-plotcommoneventgather)
    - [4.11 `plotCommonStationGather()`](#411-plotcommonstationgather)
  - [5. Examples](#5-examples)
    - [5.1 Rank Reduction Method](#51-rank-reduction-method)
    - [5.2 CCP Stacking](#52-ccp-stacking)
    - [5.3 Least Squares Migration](#53-least-squares-migration)
  - [6. Troubleshooting](#6-troubleshooting)
  - [7. Appendices](#7-appendices)
    - [7.1 Data Structure Reference](#71-data-structure-reference)
    - [7.2 Parameter Lists](#72-parameter-lists)
  - [8. Revision History](#8-revision-history)
  - [9. Authors](#9-authors)
  - [10. License](#10-license)
  - [11. Contributing](#11-contributing)
  - [12. FAQ](#12-faq)
  - [13. Contact Information](#13-contact-information)
  - [14. Acknowledgments](#14-acknowledgments)
  - [15. References](#15-references)
<!-- /TOC -->

<div style="page-break-after: always;"></div>

## 1. Introduction

**DenseArrayToolkit** 是一套面向地震密集台阵数据处理与成像的 MATLAB 工具包，提供了多种基本的数据读取、预处理、接收函数计算以及台阵处理和成像方法（例如 **Rank Reduction**、**Radon Transform**、**CCP**、**Migration Imaging** 等）功能。

- **主要功能：**
  - 快速加载 SAC 格式数据
  - 多种预处理手段（滤波、去均值、去趋势等）
  - 接收函数计算与成像
  - 台阵数据处理方法：**Rank Reduction (DRR)**、**Radon Transform**、**Migration** 等

- **使用人群：**
  - 地震学研究人员、地球物理专业学生
  - 需要针对大规模台阵数据进行半自动/自动处理的工程师

---

## 2. Installation & Requirements

1. **下载**
   - 在 GitHub 上克隆或下载该仓库，使用命令
   ```bash
   git clone https://github.com/DenseArrayToolkit/DenseArrayToolkit.git
   ```
   即可下载。
   - 或者直接从[这里](https://github.com/DenseArrayToolkit/DenseArrayToolkit/releases/download/v0.1/DenseArrayToolkit.zip)下载 `DenseArrayToolkit.zip` 。

3. **安装方法：**
   - 将 `DenseArrayToolkit` 文件夹放置在 MATLAB 搜索路径，或在 MATLAB 命令行执行:
     ```matlab
     addpath('path_to_DenseArrayToolkit');
     savepath;
     ```

4. **使用示例**
   - 在**主目录**下运行`demo/`文件夹下的脚本，即可运行示例。

5. **MATLAB 版本要求**
   - 推荐使用 R2020b 或更高版本
   - 如果要使用 **GeographicAxes** / **UIAxes** 的地图功能，需要 R2020b+

6. **依赖的工具箱：**
   - *Mapping Toolbox*（若使用与地理投影相关的绘图函数）
   - *Signal Processing Toolbox*

---

## 3. Quick Start

以下是地震数据处理的基本流程示例，展示了从**数据读取**到**接收函数可视化**的简要步骤，具体可参考`demo_Stacking.m`。

1. **Setup Paths**
   ```matlab
   setupPaths();  % 添加 DenseArrayToolkit 所需要的函数文件到 MATLAB 路径中
   ```

2. **Load Config**
   ```matlab
   config = loadConfig();  % 包含 dataFolder, PreprocessingParam
   ```

3. **Read Data**
   ```matlab
   DataStruct = read_SAC(config.dataFolder); % 读入SAC数据
   ```

4. **Preprocess**
   ```matlab
   DataStruct = preprocessing(DataStruct, config.PreprocessingParam); % 数据预处理，包括去均值、去趋势、带通滤波、重采样等。
   ```

5. **Deconvolution**
   ```matlab
   DataStruct = deconv(DataStruct, config.DeconvParam); % 计算接收函数
   ```

6. **Stack & Visualization**
   ```matlab
   [seisout, depth0] = stackCommonStationGather(DataStruct); % 叠加接收函数
   plotCommonStationGather(DataStruct, 'STATION_ID'); % 绘制台站道集
   ```

---

## 4. Modules/Functions

本小节主要介绍工具包中常用的各个模块的函数及其使用方法。

### 4.1 `read_SAC()`

该函数用于读取`dataFolder`目录下的 SAC 文件，并将其组装成包含头文件、三分量波形及相应台站、事件和时间信息的结构数组。

- **用法**：
  ```matlab
  DataStruct = read_SAC(dataFolder)
  DataStruct = read_SAC(dataFolder, maxFiles)
  ```

- **输入**：
  - `dataFolder`：包含 SAC 文件的文件夹路径。
  - `maxFiles`（可选）：要读取的最大文件数，默认为读取所有可用文件。

- **输出**：
  - `DataStruct`：结构数组，每个元素包含一个地震事件的波形和头文件数据，包括台站和事件信息。

### 4.2 `preprocessing()`

该函数对 `DataStruct` 进行一系列预处理步骤，包括信噪比（SNR）计算、去均值和趋势、ENZ 到 TRZ 的旋转、带通滤波、以及重采样等。

- **用法**：
  ```matlab
  DataStruct = preprocessing(DataStruct, PreprocessingParam)
  ```

- **输入**：
  - `DataStruct`：结构数组，每个元素通常对应一个三分量记录（E, N, Z）。
  - `PreprocessingParam`：预处理参数的结构体，包含默认字段如 `tstart`、`tend`、`lows`、`highs`、`resample_period` 等。

- **输出**：
  - `DataStruct`：预处理后的结构数组，包含处理后的波形数据和相关信息。

### 4.3 `deconv()`

该函数对单台三分量波形进行反褶积计算，并将结果存储在 `DataStruct(n).RF` 字段中。

- **用法**：
  ```matlab
  DataStruct = deconv(DataStruct, DeconvParam)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含预处理后的波形数据及相关信息。
  - `DeconvParam`：反褶积参数的结构体，包含默认字段如 `gauss`、`waterlevel`、`itmax`、`minderr` 等。

- **输出**：
  - `DataStruct`：更新后的结构数组，包含 `DataStruct(n).RF` 字段，即反褶积结果。

### 4.4 `radonTransform()`

该函数对接收函数进行 Radon 变换，利用 Radon 参数进行噪声去除或分离，最终生成 Radon 变换后的波形数据。

- **用法**：
  ```matlab
  DataStruct = radonTransform(DataStruct, RadonParam)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含处理后的波形数据及相关信息。
  - `RadonParam`：Radon 变换参数的结构体，包含默认字段如 `lows`、`highs`、`pmax`、`pmin` 等。

- **输出**：
  - `DataStruct`：更新后的结构数组，包含 Radon 变换后的波形数据。

### 4.5 `rankReduction()`

该函数对接收函数进行非规则网格重建（DRR-OTG），通过减少数据的秩来实现去噪和插值，最终生成降秩后的波形数据。

- **用法**：
  ```matlab
  [gather, d1_otg] = rankReduction(gather, gridStruct, RankReductionParam)
  ```

- **输入**：
  - `gather`：结构数组，包含接收函数数据及相关信息。
  - `gridStruct`：网格结构体，便于后续的降秩重建。
  - `RankReductionParam`：降秩重建参数的结构体，包含默认字段如 `nx`、`ny`、`rank`、`K` 等。

- **输出**：
  - `gather`：更新后的结构数组，包含降秩后的波形数据。
  - `d1_otg`：重建后的三维数据体。

### 4.6 `CCPStacking()`

该函数对接收函数进行CCP叠加，利用速度模型和CCP参数进行射线追踪和时间到深度的转换，最终生成CCP叠加结果。

- **用法**：
  ```matlab
  ccpResult = CCPStacking(DataStruct, velocityModel, CCPParam)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含接收函数数据及相关信息。
  - `velocityModel`：速度模型，用于射线追踪和时间到深度转换。
  - `CCPParam`：CCP参数的结构体，包含经纬度范围、网格参数等。

- **输出**：
  - `ccpResult`：包含CCP叠加结果的结构体。

### 4.7 `leastSquaresMig()`

该函数对接收函数进行2D最小二乘偏移，利用速度模型和最小二乘偏移参数进行偏移成像，最终生成偏移成像结果。

- **用法**：
  ```matlab
  MigResult = leastSquaresMig(gather, gridStruct, MigrationParam)
  ```

- **输入**：
  - `gather`：结构数组，包含地震记录及相关信息。
  - `gridStruct`：网格结构体，包含网格相关信息。
  - `MigrationParam`：偏移参数的结构体，包含默认字段如 `dz`、`zmax`、`itermax`、`mu` 等。

- **输出**：
  - `MigResult`：包含偏移成像结果的结构体，包括水平和深度轴、偏移图像和最小二乘偏移结果。

### 4.8 `leastSquaresMig3D()`

该函数对接收函数进行3D最小二乘偏移，利用速度模型和最小二乘偏移参数进行偏移成像，最终生成偏移成像结果。

- **用法**：
  ```matlab
  MigResult = leastSquaresMig3D(gather, gridStruct, MigrationParam)
  ```

- **输入**：
  - `gather`：结构数组，包含地震记录及相关信息。
  - `gridStruct`：网格结构体，包含网格相关信息。
  - `MigrationParam`：偏移参数的结构体，包含默认字段如 `dz`、`zmax`、`itermax`、`mu` 等。

- **输出**：
  - `MigResult`：包含偏移成像结果的结构体，包括水平和深度轴、偏移图像和最小二乘偏移结果。

### 4.9 `HKstacking()`

该函数利用接收函数进行HK叠加，计算台站下方地壳厚度和Vp/Vs比值。通过射线追踪和时间到深度的转换，最终生成HK叠加结果。

- **用法**：
  ```matlab
  HKresults = HKstacking(DataStruct, HKParam)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含台站和事件数据。
  - `HKParam`：HK叠加参数的结构体，包含地壳厚度和Vp/Vs比值的搜索范围等。

- **输出**：
  - `HKresults`：包含HK叠加结果的结构体，包括台站下方地壳厚度、Vp/Vs比值及其误差估计。

### 4.10 `plotCommonEventGather()`

该函数用于绘制单个事件的接收函数，展示不同台站的波形。

- **用法**：
  ```matlab
  plotCommonEventGather(DataStruct, EventID)
  ```

- **输入**：
  - `DataStruct`：结构数组，每个元素包含事件信息和接收函数数据。
  - `EventID`：事件ID，指定要绘制的事件。

- **输出**：
  - 无直接输出，函数在新图形中绘制指定事件的接收函数。

### 4.11 `plotCommonStationGather()`

该函数用于绘制同一台站的接收函数道集，展示不同事件的波形。

- **用法**：
  ```matlab
  plotCommonStationGather(DataStruct, 'STATION_ID')
  ```

- **输入**：
  - `DataStruct`：结构数组，每个元素包含台站信息和接收函数数据。
  - `station`：字符串，指定要绘制的台站名。

- **输出**：
  - 无直接输出，函数在新图形中绘制指定台站的接收函数道集。

---

## 5. Examples

### 5.1 Rank Reduction Method

#### 概述

`demo_rankReduction.m` 脚本提供了一个完整的工作流程示例，用于介绍阻尼秩约简方法。该脚本旨在引导用户从初始数据获取到高级数据处理和可视化的全过程。

#### 详细步骤

1. **设置路径和参数**：
   - 脚本首先通过 `setupPaths()` 和 `loadConfig()` 设置路径并加载配置参数。

2. **读取数据**：
   - 使用 `read_SAC()` 从指定文件夹读取地震波形数据。脚本检查数据是否成功加载。

3. **预处理**：
   - 使用 `preprocessing()` 对加载的地震波形数据进行预处理。

4. **计算接收函数**：
   - 使用 `deconv()` 函数进行去卷积以计算地震接收函数。

5. **阻尼降秩**：
   - 对每个地震事件，脚本使用 `rankReduction()` 执行重建。它为每个事件重建三维地震数据并存储结果。

6. **叠加接收函数**：
   - 脚本使用 `stackCommonStationGather()` 按地震台站叠加接收函数，生成叠加后的地震数据和深度轴。

7. **构造导向滤波**：
   - 后处理涉及结构平滑和倾角计算以提高地震数据质量。

#### 使用方法

要运行脚本，请确保正确设置所有必要的路径和配置。该脚本设计在 MATLAB 环境中执行，要求所需的函数和地震数据可访问。

### 5.2 CCP Stacking

#### 概述

`demo_CCPStacking.m` 脚本提供了一个完整的工作流程示例，用于介绍CCP叠加方法。该脚本旨在引导用户从初始数据获取到高级数据处理和可视化的全过程。

#### 详细步骤

1. **设置路径和参数**：
   - 使用 `setupPaths()` 和 `loadConfig()` 设置路径并加载配置参数。

2. **读取数据**：
   - 使用 `read_SAC()` 从指定文件夹读取SAC格式的地震数据。

3. **预处理**：
   - 使用 `preprocessing()` 根据配置的预处理参数对数据进行处理。

4. **获取台阵和事件信息**：
   - 使用 `getStationInfo()` 和 `getEventInfo()` 提取台站和事件的经纬度信息，并进行方位角一致性筛选。

5. **创建速度模型**：
   - 使用 `createVelocityModel()` 根据配置文件创建导入或生成速度模型。

6. **计算接收函数**：
   - 使用 `deconv()` 函数计算接收函数。

7. **CCP叠加**：
   - 使用 `CCPStacking()` 函数进行CCP叠加，生成叠加结果。

8. **结果输出**：
   - 使用 `plotCCPResults()` 和 `plotCCPXsection()` 可视化CCP叠加结果。

#### 使用方法

要运行该脚本，请确保正确设置所有必要的路径和配置。该脚本设计在 MATLAB 环境中执行，要求所需的函数和数据可访问。

### 5.3 Least Squares Migration

#### 概述

`demo_Migration_2D.m` 脚本提供了一个完整的工作流程示例，用于介绍2D地震成像方法。该脚本旨在引导用户从初始数据获取到高级数据处理和可视化的全过程。

#### 详细步骤

1. **设置路径和参数**：
   - 使用 `setupPaths()` 和 `loadConfig()` 设置路径并加载配置参数。

2. **读取数据**：
   - 使用 `read_SAC()` 从指定文件夹读取SAC格式的地震波形数据。

3. **预处理**：
   - 使用 `preprocessing()` 根据配置的预处理参数对数据进行处理。

4. **获取台阵和事件信息**：
   - 使用 `getStationInfo()` 和 `getEventInfo()` 提取台站和事件的经纬度信息，并进行方位角一致性筛选。

5. **创建规则网格**：
   - 使用 `createRegularGrid()` 根据配置文件创建规则网格。

6. **创建速度模型**：
   - 使用 `createVelocityModel()` 根据配置文件创建导入或生成速度模型。

7. **成像处理**：
   - 使用 `CCPCommonEventGather()` 和 `leastSquaresMig()` 函数进行CCP叠加和最小二乘偏移成像处理，生成成像结果。

8. **可视化**：
   - 使用 `plotCCPMigrationResults()` 可视化成像结果。

9. **保存结果**：
   - 使用 `write_MigResult()` 存储成像结果。

#### 使用方法

要运行该脚本，请确保正确设置所有必要的路径和配置。该脚本设计在 MATLAB 环境中执行，要求所需的函数和数据可访问。

---

## 6. Troubleshooting

- **无数据读入**
  - 检查 `config.dataFolder` 路径是否正确。
  - 确保 SAC 文件有 `.SAC` 后缀并具备合法头段。

- **路径配置不成功**
  - 确保 `setupPaths()` 函数已正确实现并在 MATLAB 路径中。
  - 检查配置文件路径是否正确，确保所有依赖包和函数文件都在指定路径中。
  - 如果路径中包含特殊字符或空格，可能导致路径配置失败，建议使用绝对路径。
  - 使用 `addpath` 手动添加路径时，确保路径拼写正确并使用 MATLAB 的完整路径格式。

- **Mapping Toolbox 错误**
  - 如果 `plotStations` 或地理投影函数报错，确保安装了 Mapping Toolbox；或使用 GeographicAxes 进行替代。

- **内存不足**
  - 当处理超大量台阵数据时，可能出现内存压力；可尝试 -v7.3 保存分段数据或减少一次性加载。

---

## 7. Appendices

### 7.1 Data Structure Reference

DataStruct(i) 常见字段：

.TraceData: 波形数组 [Nt x 1]
.StationInfo:
.stla, .stlo, .sta
.EventInfo:
.evla, .evlo, .evdp, .mag
.RF (若已计算接收函数)
.itr, .ittime, ...

### 7.2 Parameter Lists

<table style="width:100%; border-collapse: collapse;">
  <tr style="text-align: center; vertical-align: middle;">
    <th>组件</th>
    <th>参数名</th>
    <th>默认值</th>
    <th>说明</th>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="5"><strong>PreprocessingParam</strong></td>
    <td>`tstart`, `tend`</td>
    <td>-60,600</td>
    <td>读波形相对时间窗 (sec)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`sig_leader`</td>
    <td>30</td>
    <td>P 波到时前的时间 (sec)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`record_len`</td>
    <td>120</td>
    <td>P 波到时后的记录长度 (sec)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`lows`, `highs`</td>
    <td>0.1,2.0</td>
    <td>带通滤波边界 (Hz)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`resample_period`</td>
    <td>0.1</td>
    <td>重采样周期 (sec)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="6"><strong>DeconvParam</strong></td>
    <td>`gauss`</td>
    <td>5.0</td>
    <td>Gauss 参数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`waterlevel`</td>
    <td>0.01</td>
    <td>去卷积水位</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`itmax`</td>
    <td>100</td>
    <td>迭代最大次数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`minderr`</td>
    <td>1e-5</td>
    <td>最小误差</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`phaseshift`</td>
    <td>5</td>
    <td>相位偏移</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`verbose`</td>
    <td>true</td>
    <td>是否输出详细信息</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="5"><strong>RadonParam</strong></td>
    <td>`lows`, `highs`</td>
    <td>0.1,1.2</td>
    <td>带通滤波边界 (Hz)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`pmax`, `pmin`</td>
    <td>0.05,-0.05</td>
    <td>最大/最小慢度 (s/km)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`minTraces`</td>
    <td>60</td>
    <td>每个事件的最小道数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`N1`, `N2`</td>
    <td>30,1</td>
    <td>CG 迭代次数，外循环次数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`plotRadon`</td>
    <td>true</td>
    <td>是否绘制 Radon 结果</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="8"><strong>RankReductionParam</strong></td>
    <td>`nx`, `ny`</td>
    <td>11,14</td>
    <td>网格点数 (X, Y方向)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`rank`</td>
    <td>5</td>
    <td>DRR 矩阵秩</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`niter`</td>
    <td>5</td>
    <td>最大迭代次数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`mode`</td>
    <td>1</td>
    <td>DRR 算法模式</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`verb`</td>
    <td>true</td>
    <td>是否输出详细信息</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`eps`</td>
    <td>0.00001</td>
    <td>收敛阈值/正则化系数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`K`</td>
    <td>4</td>
    <td>DRR 参数 K</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`flow`, `fhigh`</td>
    <td>0.1,1.2</td>
    <td>频率范围 (Hz)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="7"><strong>MigParam</strong></td>
    <td>`is_ssa`</td>
    <td>0</td>
    <td>是否使用 SSA 方法</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`flow`, `fhigh`</td>
    <td>0.1,1.2</td>
    <td>迁移频率范围 (Hz)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`rank_p`</td>
    <td>8</td>
    <td>SSA 秩参数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`alpha`</td>
    <td>0.9</td>
    <td>SSA 权重参数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`n_iter`</td>
    <td>20</td>
    <td>迁移迭代次数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`gauss`</td>
    <td>2.5</td>
    <td>迁移 Gauss 参数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`phaseshift`</td>
    <td>5</td>
    <td>迁移相位参数</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="3"><strong>CCPParam</strong></td>
    <td>`imagingType`</td>
    <td>'3D'</td>
    <td>使用 2D 或 3D 成像</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`plotCCP`</td>
    <td>false</td>
    <td>是否绘制 CCP 结果</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`smoothLength`</td>
    <td>3</td>
    <td>CCP 图像平滑长度</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="5"><strong>Global Parameters</strong></td>
    <td>`max_angle_diff`</td>
    <td>15</td>
    <td>最大方位角差异 (deg)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`profile_length`</td>
    <td>4</td>
    <td>剖面长度 (degree)</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`dataFolder`</td>
    <td>'./data/event_waveforms_BY'</td>
    <td>数据文件夹路径</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`outputFolder`</td>
    <td>'./results'</td>
    <td>输出文件夹路径</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`saveResults`</td>
    <td>true</td>
    <td>是否保存结果</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td rowspan="1"><strong>visualizeResults</strong></td>
    <td>true</td>
    <td>是否可视化结果</td>
  </tr>
</table>

---

## 8. Revision History

v1.0 (2025-01-30):
- 初始发布，包含基础预处理与 Rank Reduction GUI

v1.1 (2025-02-15):
- 完善预处理模块

v1.2 (2025-03-01):
- 完善接收函数计算模块

---

## 9. Authors

- Author 1: Name, Affiliation, Email
- Author 2: Name, Affiliation, Email
- Author 3: Name, Affiliation, Email

© 2025 DenseArrayToolkit Authors. 如果您在使用中遇到问题或有功能改进需求，欢迎与我们联系或在 GitHub Issues 中提交反馈。

---

## 10. License

本软件遵循 MIT 许可证。详细信息请参见 LICENSE 文件。

---

## 11. Contributing

欢迎贡献代码和报告问题！如果你有更好的想法，欢迎加入我们，一起完善这个项目，共同为地震密集台阵数据处理提供更好的工具。

---

## 12. FAQ

- **Q**: 如何安装 DenseArrayToolkit？
  **A**: 请参阅第 2 章 [安装与要求](#2-installation--requirements)。

- **Q**: 如何报告一个 Bug？
  **A**: 请在 GitHub 上提交一个 Issue，并提供详细的错误信息和复现步骤。

---

## 13. Contact Information

如果您有任何问题或建议，请通过以下方式联系我们：
- Email: yunfeng_chen@zju.edu.cn, pengfei.zuo@zju.edu.cn
- GitHub: [DenseArrayToolkit](https://github.com/DenseArrayToolkit)

---

## 14. Acknowledgments

- 感谢 MATLAB 社区提供的强大工具和资源。
- 感谢所有贡献者对 DenseArrayToolkit 项目的支持和帮助。

---

## 15. References

- [DenseArrayToolkit](https://github.com/DenseArrayToolkit)
- [DenseArrayToolkit](https://github.com/DenseArrayToolkit)