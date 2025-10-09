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
- [Contents](#contents)
- [1. Introduction](#1-introduction)
- [2. Installation \& Requirements](#2-installation--requirements)
- [3. Quick Start](#3-quick-start)
- [4. Modules/Functions](#4-modulesfunctions)
  - [4.1 `read_SAC()`](#41-read_sac)
  - [4.2 `preprocessing()`](#42-preprocessing)
  - [4.3 `deconv()`](#43-deconv)
  - [4.4 Array Processing Module](#44-array-processing-module)
    - [4.4.1 `radonTransform2D()`](#441-radontransform2d)
    - [4.4.2 `radonTransform3D()`](#442-radontransform3d)
    - [4.4.3 `rankReduction2D()`](#443-rankreduction2d)
    - [4.4.4 `rankReduction3D()`](#444-rankreduction3d)
    - [4.4.5 `fkFilter()`](#445-fkfilter)
  - [4.5 Imaging Module](#45-imaging-module)
    - [4.5.1 `CCPCommonEventGather()`](#451-ccpcommoneventgather)
    - [4.5.2 `leastSquaresMig2D()`](#452-leastsquaresmig2d)
    - [4.5.3 `leastSquaresMig3D()`](#453-leastsquaresmig3d)
    - [4.5.4 `HKstacking()`](#454-hkstacking)
  - [4.6 Visualization Module](#46-visualization-module)
    - [4.6.1 `plotCommonEventGather()`](#461-plotcommoneventgather)
    - [4.6.2 `plotCommonStationGather()`](#462-plotcommonstationgather)
    - [4.6.3 `plotCCPResults()`](#463-plotccpresults)
    - [4.6.4 `plotMigrationResults()`](#464-plotmigrationresults)
    - [4.6.5 `plotStations()`](#465-plotstations)
    - [4.6.6 `plotEvents()`](#466-plotevents)
- [5. Examples](#5-examples)
  - [5.1 demo\_rankReduction.m —— 接收函数阻尼降秩方法重建](#51-demo_rankreductionm--阵列接收函数阻尼秩约简处理)
  - [5.2 demo\_Stacking.m —— 接收函数共转换点叠加](#52-demo_stackingm--接收函数叠加与可视化)
  - [5.3 demo\_Migration\_2D.m —— 二维地震处理与成像](#53-demo_migration_2dm--二维地震处理与成像)
  - [5.4 demo\_Migration\_3D.m —— 三维地震处理与成像](#54-demo_migration_3dm--三维地震处理与成像)
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

**DenseArrayToolkit** 是一套面向地震密集台阵数据处理与成像的 MATLAB 工具包。它提供了一系列基础与高级功能，包括数据 I/O、预处理、接收函数计算、台阵数据处理和成像方法（例如阻尼降秩方法、拉东变换方法、共转换点叠加方法、接收函数偏移方法等）功能。

- **主要功能：**
  - 快速读取 SAC 格式地震波形数据
  - 常规数据预处理手段（滤波、去均值、去趋势、重采样、计算信噪比等）
  - 接收函数计算与成像（hk叠加、共转换点叠加、最小二乘偏移）
  - 台阵数据处理方法（拉东变换、阻尼降秩方法等）

- **使用人群：**
  - 地震学研究人员、地球物理专业本科生与研究生或对地震数据处理感兴趣的初学者
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
     addpath(genpath('path_to_DenseArrayToolkit'));
     savepath;
     ```

4. **使用示例**
   - 需要在**主目录**下运行`demo/`文件夹下的脚本，即可运行示例。

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
   config = loadConfig();  % 加载默认参数，包含 dataFolder, PreprocessingParam等
   ```

3. **Read Data**
   ```matlab
   DataStruct = read_SAC(config.dataFolder); % 输入SAC文件所在的文件路径，读入SAC数据
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

该函数对 `DataStruct` 结构体下的原始波形进行一系列预处理步骤，包括信噪比（SNR）计算、去均值和趋势、ENZ 到 TRZ 的旋转、带通滤波、以及重采样等。

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
  - `DataStruct`：更新后的结构数组，包含 `DataStruct(n).RF` 变量，即反褶积结果。

### 4.4 Array Processing Module

array_processing 模块是 DenseArrayToolkit 的核心处理模块，专门针对密集台阵数据的特点（低信噪比，分布不规则）开发了一系列先进的信号处理算法。该模块主要用于提高地震数据的信噪比、分离相干噪声与有效信号、以及改善数据的空间连续性。这些方法特别适用于处理不规则分布的台阵数据和压制面波等相干噪声。

**模块特点：**
- **多维处理能力**：支持 2D 和 3D 数据处理，适应不同维度的台阵布局
- **噪声压制**：专门针对地震数据中的相干噪声和随机噪声设计
- **信号增强**：通过数学变换增强微弱的地下结构信号
- **适应性**：适用于规则和不规则台阵分布的数据处理

**应用场景：**
- 压制随机噪声和散射噪声
- 分离不同传播路径的地震波
- 改善接收函数数据的空间连续性
- 提高后续成像处理的质量

#### 4.4.1 `radonTransform2D()`

该函数对接收函数进行二维 Radon 变换，利用慢度域（slowness domain）分离噪声与有效信号，适用于规则或近似规则线性台阵数据的噪声抑制与信号增强。

- **用法**：
  ```matlab
  DataStruct = radonTransform2D(DataStruct, gridStruct, param)
  ```

- **输入参数说明**：
  - `DataStruct`：结构数组，包含事件、波形、行程等信息，需含 `.EventInfo.evid`, `.Waveforms.dataProcessed`, `.TravelInfo.distDeg` 等字段。
  - `gridStruct`：离散网格结构体。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | lows           | 0.1         | 带通滤波低频 (Hz)            |
    | highs          | 1.2         | 带通滤波高频 (Hz)            |
    | pmax/pmin      | 0.05/-0.05  | 最大/最小慢度 (s/km)         |
    | minTraces      | 60          | 最小道数                     |
    | N1/N2          | 30/1        | 共轭梯度迭代次数/稀疏权重    |
    | plotRadon      | false       | 是否绘制 Radon 结果          |
    | order          | postdecon   | 滤波顺序（反褶积前/后）      |

- **输出**：
  - `DataStruct`：每个事件的 `.Waveforms.dataRadonFiltered` 字段填充为 Radon 处理后的三分量波形。

- **算法简要说明**：
  1. 对每个事件，按投影距离对台站排序，提取并带通滤波三分量波形。
  2. 通过共轭梯度法（CG）在慢度域进行 Radon 变换与反变换，分离噪声与有效信号。
  3. 支持反褶积前/后处理（order），反褶积前对原始数据垂直和径向分量地震记录进行滤波，反褶积后则直接对接收函数进行滤波，可选绘图。
  4. 适用于台阵分布较规则、噪声具有一定空间相干性的场景。

- **适用条件**：
  - 适合规则或近似规则线性台阵
  - 事件道数较多（>minTraces）
  - 主要用于随机噪声、散射噪声抑制及信号增强

#### 4.4.2 `radonTransform3D()`

该函数对接收函数进行三维 Radon 变换，利用二维慢度域（px, py）分离噪声与有效信号，适用于面状密集台阵的噪声抑制与信号增强。

- **用法**：
  ```matlab
  DataStruct = radonTransform3D(DataStruct, gridStruct, param)
  ```

- **输入参数说明**：
  - `DataStruct`：结构数组，包含事件ID、波形、时间等信息，需含 `.EventInfo.evid`, `.Waveforms.dataProcessed`, `.TravelInfo.distDeg` 等字段。
  - `gridStruct`：生成的离散网格结构体。
  - `param`：参数变量结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | lows           | 0.1         | 带通滤波低频 (Hz)            |
    | highs          | 1.2         | 带通滤波高频 (Hz)            |
    | pxmax/pxmin    | 0.05/-0.05  | x方向最大/最小慢度 (s/km)    |
    | pymax/pymin    | 0.05/-0.05  | y方向最大/最小慢度 (s/km)    |
    | minTraces      | 60          | 最小道数                     |
    | N1/N2          | 30/1        | 共轭梯度迭代次数/稀疏权重    |
    | plotRadon      | false       | 是否绘制 Radon 结果          |
    | order          | postdecon   | 滤波顺序（反褶积前/后）      |
    | type           | 1           | Radon类型(1=线性,2=抛物线,3=双曲线) |

- **输出**：
  - `DataStruct`：每个事件的 `.Waveforms.dataRadonFiltered` 变量为 Radon 处理后的三分量波形。

- **算法简要说明**：
  1. 对每个事件，按投影x/y坐标对台站排序，提取并带通滤波三分量波形。
  2. 在二维慢度域(px, py)上通过共轭梯度法（CG）进行Radon变换与反演，分离噪声与有效信号。
  3. 支持线性（Linear）、抛物线（Parabolic）、双曲线（Hyperbolic）等多种Radon类型。
  4. 适用于面状密集台阵、空间分布较均匀、噪声具有空间相干性的场景。

- **适用条件**：
  - 适合面状密集台阵（二维分布）
  - 事件道数较多（>minTraces）
  - 主要用于随机噪声、散射噪声抑制及信号增强

#### 4.4.3 `rankReduction2D()`

该函数可对接收函数进行二维阻尼降秩（Damped Rank Reduction, DRR-OTG）处理，通过矩阵降秩方法实现噪声抑制和信号空间重建，适用于规则或近似规则线性台阵。

- **用法**：
  ```matlab
  [gather, reconGrid] = rankReduction2D(gather, gridStruct, param)
  ```

- **输入参数说明**：
  - `gather`：结构数组，包含台站、接收函数等信息，需含 `.StationInfo.stlo`, `.StationInfo.stla`, `.RF.itr`, `.RF.ittime` 等字段。
  - `gridStruct`：离散网格结构体，包含规则网格坐标和不规则网格坐标。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | nx/ny          | 必填        | 网格点数（x/y方向）          |
    | rank           | 10          | DRR矩阵的秩                    |
    | K              | 5           | DRR阻尼参数                  |
    | niter          | 20          | 最大迭代次数                 |
    | eps            | 1e-3        | 正则化/收敛阈值              |
    | verb           | true        | 是否输出详细信息             |
    | mode           | 1           | 算法模式                     |
    | flow/fhigh     | 0.1/1.2     | 频率范围 (Hz)                |
    | tmax           | 30          | 最大时间窗 (s)               |
    | plotRankReduction | false    | 是否绘制降秩结果             |

- **输出**：
  - `gather`：更新后的结构数组，`.RF.itr`字段为降秩重建后的接收函数。
  - `reconGrid`：重建后的二维空间网格（x）。

- **算法简要说明**：
  1. 将所有台站的接收函数投影到x轴，构建时间-x的二维数据体。
  2. 采用DRR-OTG方法对数据体进行矩阵降秩重建，抑制噪声、填补缺失。
  3. 支持不规则台阵的空间插值与重建，参数可调节降秩强度与收敛速度。
  4. 适用于线性或近线性台阵、空间连续性要求高的场景。

- **适用条件**：
  - 适合规则或近似规则线性台阵
  - 台站空间分布较均匀
  - 主要用于噪声抑制、缺失数据插值、信号增强

#### 4.4.4 `rankReduction3D()`

该函数对接收函数进行三维阻尼降秩（Damped Rank Reduction, DRR-OTG）处理，通过三维张量降秩实现噪声抑制和空间重建，适用于面状密集台阵。

- **用法**：
  ```matlab
  [gather, reconGrid] = rankReduction3D(gather, gridStruct, param)
  ```

- **输入参数说明**：
  - `gather`：结构数组，包含台站、接收函数等信息，需含 `.StationInfo.stlo`, `.StationInfo.stla`, `.RF.itr`, `.RF.ittime` 等字段。
  - `gridStruct`：离散网格结构体，包含规则网格坐标和不规则网格坐标。
  - `param`：降秩所需变量结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | nx/ny          | 必填        | 网格点数（x/y方向）          |
    | rank           | 10          | DRR张量秩                    |
    | K              | 5           | DRR阻尼参数                  |
    | niter          | 20          | 最大迭代次数                 |
    | eps            | 1e-3        | 正则化/收敛阈值              |
    | verb           | true        | 是否输出详细信息             |
    | mode           | 1           | 算法模式                     |
    | flow/fhigh     | 0.1/1.2     | 频率范围 (Hz)                |
    | tmax           | 30          | 最大时间窗 (s)               |
    | plotRankReduction | false    | 是否绘制降秩结果             |

- **输出**：
  - `gather`：更新后的结构数组，`.RF.itr`字段为降秩重建后的接收函数。
  - `reconGrid`：重建后的三维空间网格（x, y）。

- **算法简要说明**：
  1. 将所有台站的接收函数投影到x/y平面，构建时间-x-y的三维数据体。
  2. 采用DRR-OTG方法对三维数据体进行张量降秩重建，抑制噪声、填补空间缺失。
  3. 支持不规则台阵的空间插值与重建，参数可调节降秩强度与收敛速度。
  4. 适用于面状密集台阵、空间连续性要求高的场景。

- **适用条件**：
  - 适合面状密集台阵（二维分布）
  - 台站空间分布较均匀
  - 主要用于噪声抑制、缺失数据插值、信号增强

#### 4.4.5 `fkFilter()`

该函数对接收函数进行频率-波数（f-k）锥形滤波，抑制空间一致性强的面波、散射噪声等，适用于规则或近似规则台阵。

- **用法**：
  ```matlab
  DataStruct = fkFilter(DataStruct, gridStruct, param)
  ```

- **输入参数说明**：
  - `DataStruct`：结构数组，包含事件、波形、行程等信息，需含 `.EventInfo.evid`, `.Waveforms.dataProcessed`, `.TravelInfo.distDeg` 等字段。
  - `gridStruct`：台阵投影坐标下网格结构体，包含规则网格坐标和不规则网格坐标。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | lows           | 0.1         | 带通滤波低频 (Hz)            |
    | highs          | 1.2         | 带通滤波高频 (Hz)            |
    | minTraces      | 60          | 最小道数                     |
    | w              | 0.1         | 锥形滤波半宽度（百分比）     |
    | plotFK         | false       | 是否绘制f-k滤波结果          |
    | plotFKspectrum | false       | 是否绘制f-k谱                |
    | order          | postdecon   | 滤波顺序（反褶积前/后）      |

- **输出**：
  - `DataStruct`：每个事件的 `.Waveforms.dataFKFiltered` 字段填充为f-k滤波后的三分量波形。

- **算法简要说明**：
  1. 对每个事件，按投影距离对台站排序，提取并带通滤波三分量波形。
  2. 对不规则台阵插值到规则网格，进行二维f-k变换，应用特定波数范围设计滤波器来抑制噪声。
  3. 滤波后插值回原始台站位置，恢复物理空间波形。
  4. 支持反褶积前/后处理（order），可选绘图和谱分析。
  5. 适用于随机噪声、散射噪声空间一致性强的场景。

- **适用条件**：
  - 适合规则或近似规则线性台阵
  - 事件道数较多（>minTraces）
  - 主要用于随机噪声、空间一致性噪声抑制、信号增强

### 4.5 Imaging Module

imaging 模块是 DenseArrayToolkit 的核心成像模块，实现了多种先进的地震成像算法，用于将接收函数数据转换为地下结构的图像。该模块结合了传统的 CCP 叠加方法和最新的最小二乘偏移技术，能够提供高分辨率的地下结构图像。

**应用场景：**
- 地壳和上地幔结构的精细成像
- 莫霍面深度和形态研究
- 地壳内部不连续面探测
- 地壳厚度和 Vp/Vs 比值估算
- 构造边界和断层系统的成像

#### 4.5.1 `CCPCommonEventGather()`

该函数实现共转换点叠加（Common Conversion Point stacking, CCP），通过射线追踪和时深转换，将接收函数能量投影到地下空间网格，实现2D或3D地壳结构成像。

- **用法**：
  ```matlab
  ccpResult = CCPCommonEventGather(gather, gridStruct, param)
  ```

- **输入参数说明**：
  - `gather`：结构数组，包含单事件或道集的接收函数数据，需含 `.RF.itr`, `.RF.ittime`, `.TravelInfo.rayParam`, `.TravelInfo.baz`, `.StationInfo.stla`, `.StationInfo.stlo` 等字段。
  - `gridStruct`：生成的网格参数结构体。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | imagingType    | '2D'        | 成像类型（'2D'或'3D'）       |
    | plotCCP        | false       | 是否绘制CCP成像结果          |
    | smoothLength   | 0           | CCP图像平滑核长度            |

- **输出**：
  - `ccpResult`：CCP结果结构体，包含成像网格（X/Z或X/Y/Z）、CCP叠加图像`img`、每格样本数`count`等。

- **算法简要说明**：
  1. 提取每道接收函数的射线参数、方位角、台站坐标等信息。
  2. 根据速度模型和射线参数，进行射线追踪，计算每道接收函数在地下的转换点位置。
  3. 若为2D成像，将所有转换点投影到剖面，按空间网格进行Bining并累加RF振幅，归一化获得CCP图像。
  4. 若为3D成像，将转换点投影到三维网格，进行空间Bining与累加。
  5. 支持3D速度模型下的时深校正与空间插值，提升成像精度。
  6. 可选平滑处理与可视化，输出空间分布的CCP能量图像。

- **适用条件**：
  - 适合二维剖面或三维台阵区域
  - 速度模型已知（支持1D/2D/3D）
  - 主要用于地壳界面、莫霍面等转换点成像

#### 4.5.2 `leastSquaresMig2D()`

该函数对接收函数进行二维最小二乘偏移（Least-squares Migration），通过最小化观测数据与正演模拟数据的残差，获得高分辨率地下结构成像，适用于二维台阵剖面。

- **用法**：
  ```matlab
  MigResult = leastSquaresMig2D(gather, gridStruct, param)
  ```

- **输入参数说明**：
  - `gather`：结构数组，包含台站、接收函数等信息，需含 `.StationInfo.stlo`, `.StationInfo.stla`, `.RF.itr`, `.RF.ittime` 等字段。
  - `gridStruct`：二维成像网格结构体，包含水平和深度轴信息。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | dz             | 1           | 成像深度步长 (km)            |
    | zmax           | 100         | 最大成像深度 (km)            |
    | itermax        | 20          | 最大迭代次数                 |
    | mu             | 0.1         | 正则化参数                   |
    | flow/fhigh     | 0.1/1.2     | 频率范围 (Hz)                |
    | gauss          | 2.5         | 震源子波主频参数               |
    | plotMig        | false       | 是否绘制成像结果             |

- **输出**：
  - `MigResult`：偏移结构体，包含水平和深度轴、偏移图像、最小二乘偏移结果等。

- **算法简要说明**：
  1. 构建二维成像网格，利用速度模型进行射线追踪，建立观测与模型的映射关系。
  2. 采用最小二乘反演，迭代优化地下结构模型，使正演模拟数据与观测接收函数残差最小。
  3. 支持正则化与高斯平滑，提升成像分辨率与稳定性。
  4. 适用于二维台阵剖面、地壳厚度和界面成像等场景。

- **适用条件**：
  - 适合二维台阵剖面
  - 速度模型已知或可估算
  - 主要用于高分辨率地壳/上地幔结构成像

#### 4.5.3 `leastSquaresMig3D()`

该函数对接收函数进行三维最小二乘偏移（Least-squares Migration），通过最小化观测数据与三维正演模拟数据的残差，获得高分辨率三维地下结构成像，适用于面状密集台阵。

- **用法**：
  ```matlab
  MigResult = leastSquaresMig3D(gather, gridStruct, param)
  ```

- **输入参数说明**：
  - `gather`：结构数组，包含台站、接收函数等信息，需含 `.StationInfo.stlo`, `.StationInfo.stla`, `.RF.itr`, `.RF.ittime` 等字段。
  - `gridStruct`：三维成像网格结构体，包含水平（x, y）和深度轴信息。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | dz             | 1           | 成像深度步长 (km)            |
    | zmax           | 100         | 最大成像深度 (km)            |
    | itermax        | 20          | 最大迭代次数                 |
    | mu             | 0.1         | 正则化参数                   |
    | flow/fhigh     | 0.1/1.2     | 频率范围 (Hz)                |
    | gauss          | 2.5         | 震源子波主频参数               |
    | plotMig        | false       | 是否绘制成像结果             |

- **输出**：
  - `MigResult`：结构体，包含水平（x, y）和深度轴、三维偏移图像、最小二乘偏移结果等。

- **算法简要说明**：
  1. 构建三维成像网格，利用速度模型进行三维射线追踪，建立观测与模型的映射关系。
  2. 采用三维最小二乘反演，迭代优化地下结构模型，使正演模拟数据与观测接收函数残差最小。
  3. 支持正则化与高斯平滑，提升三维成像分辨率与稳定性。
  4. 适用于面状密集台阵、地壳厚度和界面三维成像等场景。

- **适用条件**：
  - 适合面状密集台阵（三维分布）
  - 速度模型已知或可估算
  - 主要用于高分辨率地壳/上地幔三维结构成像

#### 4.5.4 `HKstacking()`

该函数利用接收函数进行HK叠加（H-K stacking），通过分析Ps及其多次波的到时，反演台站下方地壳厚度（H）和Vp/Vs比值（K），适用于单台站或台阵的地壳参数估算。

- **用法**：
  ```matlab
  HKresults = HKstacking(DataStruct, param)
  ```

- **输入参数说明**：
  - `DataStruct`：结构数组，包含台站、事件、接收函数等信息，需含 `.StationInfo.stla`, `.StationInfo.stlo`, `.RF.itr`, `.RF.ittime` 等字段。
  - `param`：结构体，主要参数如下：

    | 参数名         | 默认值      | 说明                         |
    | -------------- | ----------- | ---------------------------- |
    | Hrange         | [20 60]     | 地壳厚度搜索范围 (km)        |
    | Krange         | [1.6 2.0]   | Vp/Vs比值搜索范围            |
    | dH             | 0.1         | H步长 (km)                   |
    | dK             | 0.01        | K步长                        |
    | Vp             | 6.3         | 地壳P波速度 (km/s)           |
    | weight         | [0.7 0.2 0.1]| Ps, PpPs, PsPs+PpSs权重      |
    | plotHK         | false       | 是否绘制HK叠加曲面           |

- **输出**：
  - `HKresults`：结构体，包含每个台站的最优地壳厚度H、Vp/Vs比值K、误差估计、HK叠加曲面等。

- **算法简要说明**：
  1. 对每个台站，遍历H-K参数网格，计算理论Ps、PpPs、PsPs+PpSs到时。
  2. 在接收函数上提取对应振幅，按权重累加，获得HK叠加能量曲面。
  3. 取能量最大点为最优H、K，支持误差估计与可视化。
  4. 适用于单台站或台阵地壳厚度、Vp/Vs比值反演。

- **适用条件**：
  - 适合单台站或台阵
  - 事件方位和距离分布较均匀
  - 主要用于地壳厚度、Vp/Vs比值估算

### 4.6 Visualization Module

visualization 模块提供了一系列数据可视化函数，用于展示地震数据、接收函数和成像结果，帮助用户直观地分析和解释数据。

#### 4.6.1 `plotCommonEventGather()`

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

#### 4.6.2 `plotCommonStationGather()`

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

#### 4.6.3 `plotCCPResults()`

该函数用于可视化CCP叠加结果，显示地下结构的图像。

- **用法**：
  ```matlab
  plotCCPResults(ccpResult)
  ```

- **输入**：
  - `ccpResult`：包含CCP叠加结果的结构体。

- **输出**：
  - 无直接输出，函数在新图形中绘制CCP叠加结果。

#### 4.6.4 `plotMigrationResults()`

该函数用于可视化偏移成像结果，显示地下结构的反射特征。

- **用法**：
  ```matlab
  plotMigrationResults(MigResult)
  ```

- **输入**：
  - `MigResult`：包含偏移成像结果的结构体。

- **输出**：
  - 无直接输出，函数在新图形中绘制偏移成像结果。

#### 4.6.5 `plotStations()`

该函数用于绘制台站分布图，显示研究区域的台站位置。

- **用法**：
  ```matlab
  plotStations(DataStruct)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含台站位置信息。

- **输出**：
  - 无直接输出，函数在新图形中绘制台站分布图。

#### 4.6.6 `plotEvents()`

该函数用于绘制事件分布图，显示研究区域的地震事件位置。

- **用法**：
  ```matlab
  plotEvents(DataStruct)
  ```

- **输入**：
  - `DataStruct`：结构数组，包含事件位置信息。

- **输出**：
  - 无直接输出，函数在新图形中绘制事件分布图。

---

## 5. Examples

本节详细介绍 demo 文件夹下的四个核心示例脚本，涵盖 DenseArrayToolkit 的主要功能模块。每个示例均包含完整的处理流程、输入输出说明和典型运行方法，适合新用户快速上手和理解工具包的实际应用。

**示例程序特点：**
- **完整性**：每个示例均包含从数据读取、预处理、信号处理到结果可视化的全流程
- **实用性**：基于真实地震数据处理需求设计
- **可重复性**：参数配置标准化，便于复现
- **教育性**：脚本内含详细注释，便于学习

**学习建议：**
- 建议按顺序运行示例程序，从基础到高级逐步学习
- 可根据实际数据和需求调整参数
- 示例脚本可作为实际项目的模板

---

### 5.1 demo_rankReduction.m —— 阵列接收函数阻尼秩约简处理

**功能简介：**
本示例演示如何对地震接收函数进行阻尼降秩（Damped Rank Reduction, DRR-OTG）处理，实现信号增强与噪声抑制。适用于不规则台阵数据的信号重建。

**主要流程：**
1. 路径与参数设置：setupPaths()、loadConfig()，指定数据目录和参数
2. 数据读取：read_SAC() 读取 SAC 格式地震数据
3. 预处理：preprocessing() 进行滤波、去均值、重采样等
4. 接收函数计算：deconv() 反褶积获得接收函数
5. 台阵数据降秩处理：rankReduction3D() 对每个事件进行三维重建
6. 叠加：stackCommonStationGather() 按台站叠加接收函数
7. 结果输出：可保存处理结果至 mat 文件

**输入输出说明：**
- 输入：SAC 格式地震数据目录
- 输出：降秩后的接收函数、叠加结果 seisout/depth0、可选 Moho 结构 mohoStruct

**运行方法：**
在 MATLAB 命令行运行：
```matlab
demo_rankReduction
```
确保配置文件和数据路径正确。

---

### 5.2 demo_Stacking.m —— 接收函数叠加与可视化

**功能简介：**
本示例演示从数据读取、接收函数计算到按台站叠加和多种可视化的完整流程，适合初学者快速体验接收函数分析。

**主要流程：**
1. 路径与参数设置：setupPaths()、loadConfig()
2. 数据读取：read_SAC()
3. 预处理：preprocessing()
4. 接收函数计算：deconv()
5. 叠加：stackCommonStationGather()
6. 可视化：
   - plotWaveforms() 单道波形
   - plotCommonStationGather() 台站叠加结果
   - plotCommonEventGather() 事件道集
   - plotStations()/plotEvents() 台站与事件分布

**输入输出说明：**
- 输入：SAC 格式地震数据目录
- 输出：叠加接收函数、可视化图像（可自动保存 PNG）

**运行方法：**
在 MATLAB 命令行运行：
```matlab
demo_Stacking
```
可根据需要修改脚本中的参数和可视化选项。

---

### 5.3 demo_Migration_2D.m —— 二维地震处理与成像

**功能简介：**
本示例演示二维台阵地震成像的完整流程，包括 CCP 叠加与最小二乘偏移（Least-squares Migration）两种方法的对比。适用于二维剖面区域的高分辨率成像。

**主要流程：**
1. 路径与参数设置：setupPaths()、loadConfig()
2. 数据读取：read_SAC()
3. 预处理：preprocessing()
4. 台阵与事件信息提取：getStations()/getEvents()
5. filter_earthquakes_by_azimuth() 筛选台阵走向方向的地震
6. 网格与速度模型构建：createGrid()、getVelocityModel()
7. 成像处理：
   - radonTransform2D() 增强信噪比
   - CCPCommonEventGather() 传统 CCP 叠加
   - leastSquaresMig2D() 最小二乘偏移
8. 可视化：plotCCPMigrationResults() 对比显示
9. 结果保存：write_MigResult() 输出成像结果

**输入输出说明：**
- 输入：SAC 格式地震数据目录
- 输出：CCP 叠加与偏移成像结果、可视化图像、mat 文件

**运行方法：**
在 MATLAB 命令行运行：
```matlab
demo_Migration_2D
```
可根据实际剖面和数据调整参数。

---

### 5.4 demo_Migration_3D.m —— 三维地震处理与成像

**功能简介：**
本示例演示三维台阵地震成像流程，结合三维最小二乘偏移和三维阻尼降秩重建，适用于大规模面状密集台阵成像。

**主要流程：**
1. 路径与参数设置：setupPaths()、loadConfig()
2. 数据读取：read_SAC()
3. 预处理：preprocessing()
4. 台阵与事件信息提取：getStations()/getEvents()
5. 三维网格与速度模型构建：createGrid()、getVelocityModel()
6. 三维秩约简：rankReduction3D()
7. 三维成像处理：leastSquaresMig3D()
8. 叠加：stackImagingResults()
9. 可视化：visualizeImage() 

**输入输出说明：**
- 输入：SAC 格式地震数据目录
- 输出：三维成像体、可视化剖面、mat 文件

**运行方法：**
在 MATLAB 命令行运行：
```matlab
demo_Migration_3D
```
可根据实际台阵和数据规模调整参数。

---

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
    <td>震源函数 Gauss 参数（与接收函数计算一致）</td>
  </tr>
  <tr style="text-align: center; vertical-align: middle;">
    <td>`phaseshift`</td>
    <td>5</td>
    <td>相位时移参数（与接收函数计算一致）</td>
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

- Yunfeng Chen, Zhejiang University, yunfeng_chen@zju.edu.cn
- Pengfei Zuo, Zhejiang University, pengfei.zuo@zju.edu.cn
- Yangkang Chen, The University of Texas at Austin, yangkang.chen@beg.utexas.edu

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

- [DenseArrayToolkit](https://github.com/PengfeiZuo001/DenseArrayToolkit)
