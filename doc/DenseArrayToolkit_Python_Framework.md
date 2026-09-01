# DenseArrayToolkit Python 版本 —— 项目架构与开发框架

> **版本**: v0.1 (Draft)  
> **日期**: 2026-06-27  
> **基于**: DenseArrayToolkit MATLAB 版本 (commit a8f2d27)  
> **目标**: 构建高性能、易扩展、生态兼容的 Python 密集台阵地震学工具包

---

## 目录

1. [项目背景与动机](#1-项目背景与动机)
2. [Python 生态优势](#2-python-生态优势)
3. [总体架构设计](#3-总体架构设计)
4. [模块划分与实现计划](#4-模块划分与实现计划)
5. [核心数据模型设计](#5-核心数据模型设计)
6. [Pipeline API 设计](#6-pipeline-api-设计)
7. [计算引擎层设计](#7-计算引擎层设计)
8. [开发路线图](#8-开发路线图)
9. [工作量评估（修订版）](#9-工作量评估修订版)
10. [测试与验证策略](#10-测试与验证策略)
11. [项目文件结构](#11-项目文件结构)
12. [依赖清单](#12-依赖清单)
13. [编码规范与贡献指南](#13-编码规范与贡献指南)

---

## 1. 项目背景与动机

### 1.1 MATLAB 版本现状

DenseArrayToolkit 是一套面向密集台阵数据的接收函数处理与成像 MATLAB 工具包，具备以下核心能力：

| 类别 | 功能 | 成熟度 |
|------|------|--------|
| 预处理 | 波形预处理、滤波、SAC 读写 | ✅ 成熟 |
| 反褶积 | 迭代反褶积提取接收函数 | ✅ 成熟 |
| 台阵处理 | DRR-OTG 阻尼秩约简 (2D/3D) | ✅ 成熟 |
| 台阵处理 | Radon 变换 (2D/3D) | ✅ 成熟 |
| 台阵处理 | FK 滤波 | ✅ 成熟 |
| 台阵处理 | 结构导向滤波 (MATseistr) | ✅ 成熟 |
| 台阵处理 | 自适应匹配滤波 (AMF) | ✅ 成熟 |
| 成像 | CCP 叠加 (2D/3D, 含 Fresnel 加权) | ✅ 成熟 |
| 成像 | 最小二乘偏移 LSM (2D/3D) | ✅ 成熟 |
| 成像 | HK 叠加 | ✅ 成熟 |
| 可视化 | CCP 剖面、台站/事件分布、波形展示 | ✅ 成熟 |

### 1.2 迁移到 Python 的核心驱动力

1. **完全免费开源**: MATLAB 许可证费用 (~$2,150/年 + 工具箱) 限制了许多高校和发展中国家研究者的使用
2. **与地震学生态无缝对接**: ObsPy 是地震学数据处理的事实标准
3. **现代化开发体验**: 类型检查、包管理 (pip/conda)、CI/CD、容器化
4. **大规模计算**: Python 原生支持 HPC (Dask, Ray, mpi4py)，远超 MATLAB Parallel Server 的灵活性
5. **社区协作**: 开源社区更容易为 Python 项目贡献代码
6. **可复现性**: Docker 镜像 + Jupyter Notebook 实现完整可复现研究流程

### 1.3 已有 Python 生态基础

以下 Python 包已实现相关功能，可极大降低开发工作量：

| 功能 | 已有 Python 包 | 状态 |
|------|---------------|------|
| 阻尼秩约简 (DRR/DRR3D) | **Pydrr** | ✅ 可直接整合 |
| Radon 变换 2D | **PyRadon** / 相关实现 | ✅ 可直接整合 |
| SAC 读写 | **ObsPy** | ✅ 工业标准 |
| 坐标投影 | **Cartopy / pyproj** | ✅ 成熟 |
| 地震学通用处理 | **ObsPy** | ✅ 工业标准 |

**需要全新开发的核心模块**：
- 最小二乘偏移 LSM (2D/3D) —— 无现有 Python 实现
- CCP/HK 叠加的专用实现 —— 可借鉴 `rf` 和 `seispy` 包
- 台阵处理工作流集成 —— 胶水层

---

## 2. Python 生态优势

### 2.1 与 MATLAB 的技术对比

| 维度 | MATLAB | Python |
|------|--------|--------|
| **许可证** | 收费 (个人/学术/商业) | 完全免费 |
| **地震学核心库** | 无标准生态 | ObsPy (覆盖 20+ 格式、走时计算、事件处理) |
| **矩阵计算** | 内置 | NumPy + SciPy (功能超集) |
| **并行计算** | Parallel Toolbox (付费) | Dask / Ray / mpi4py (免费) |
| **GPU 加速** | 需要额外配置 | CuPy / PyTorch / JAX |
| **深度学习整合** | 受限 | PyTorch / TensorFlow 生态 |
| **包管理** | 无官方方案 | pip / conda / poetry |
| **容器化** | 因许可证难以分发 | Docker 镜像一键部署 |
| **交互式文档** | Live Script | Jupyter / Quarto (功能更强) |
| **CI/CD** | 受限 | GitHub Actions (免费) |
| **社区** | 封闭商业 | 全球最大开源社区 |

### 2.2 地震学 Python 生态全景

```
                    ┌──────────────┐
                    │  DenseArray  │  ← 本工具包
                    │   Toolkit    │
                    │   (Python)   │
                    └──────┬───────┘
                           │
          ┌────────────┬───┴───┬────────────┐
          │            │       │            │
    ┌─────▼─────┐ ┌───▼───┐ ┌─▼──────┐ ┌───▼────┐
    │   ObsPy   │ │ Pydrr │ │ PyGMT  │ │Cartopy │
    │ (SAC/波形)│ │(DRR)  │ │(绘图)  │ │(地图)  │
    └───────────┘ └───────┘ └────────┘ └────────┘
          │
    ┌─────▼─────────────────────────┐
    │  NumPy / SciPy / Numba / JAX  │
    │        (计算基础设施)           │
    └───────────────────────────────┘
```

---

## 3. 总体架构设计

### 3.1 设计哲学

> **OOP 做调度，NumPy 做计算。**
>
> 核心数值算法保持纯函数形式（NumPy 数组入 → NumPy 数组出），  
> 面向对象封装仅用于数据管理、工作流编排、参数校验和可视化。

### 3.2 三层架构

```
┌──────────────────────────────────────────────────────────┐
│                    User Layer                             │
│   Demo Scripts  │  Jupyter Notebooks  │  CLI / Web API   │
├──────────────────────────────────────────────────────────┤
│                 Pipeline Layer                            │
│   SeismicPipeline  │  RFProcessor  │  ImagingProcessor   │
│   (Builder Pattern, 链式调用, 参数管理, 进度日志)          │
├──────────────────────────────────────────────────────────┤
│                  Domain Layer                             │
│   Trace │ Station │ Event │ RF │ Gather │ VelocityModel   │
│   ImageGrid │ CCPResult │ HKResult │ MigrationResult     │
│   (@dataclass + pydantic 校验, .plot() 方法)              │
├──────────────────────────────────────────────────────────┤
│               Computation Engine                          │
│   rank_reduction_3d() │ radon_transform_3d()             │
│   ccp_stack() │ hk_stack() │ least_squares_migration()   │
│   fk_filter() │ iter_deconv()                             │
│   (纯函数, 输入/输出 ndarray, Numba 可选加速)              │
└──────────────────────────────────────────────────────────┘
```

### 3.3 设计原则

| 原则 | 说明 |
|------|------|
| **关注点分离** | 数据 (Domain) ≠ 算法 (Engine) ≠ 流程 (Pipeline) |
| **纯函数优先** | 核心计算不依赖外部状态，输入输出明确 |
| **组合优于继承** | 管道通过组合多个处理步骤构建 |
| **类型安全** | 全链路类型注解 + mypy 静态检查 |
| **渐进增强** | 纯 NumPy 实现 → Numba JIT 加速 → GPU (CuPy/JAX) |
| **可测试性** | 每层独立可测，核心算法对比 MATLAB 输出验证 |

---

## 4. 模块划分与实现计划

### 4.1 模块总览

```
densearray_toolkit/
│
├── densearray/
│   ├── __init__.py
│   │
│   ├── core/                       # Domain Layer — 核心数据模型
│   │   ├── __init__.py
│   │   ├── trace.py                # Trace, Station, Event
│   │   ├── rf.py                   # ReceiverFunction
│   │   ├── gather.py               # CommonEventGather, CommonStationGather
│   │   ├── velocity.py             # VelocityModel1D/2D/3D
│   │   ├── grid.py                 # ImageGrid, ProfileGrid
│   │   └── config.py               # Config, ProcessingParams
│   │
│   ├── io/                         # I/O 层
│   │   ├── __init__.py
│   │   ├── sac.py                  # SAC 读写 (基于 ObsPy)
│   │   ├── events.py               # 事件目录读写
│   │   ├── stations.py             # 台站文件读写
│   │   └── velocity_model.py       # 速度模型文件读写
│   │
│   ├── preprocessing/              # 预处理
│   │   ├── __init__.py
│   │   ├── filters.py              # 带通/高通/低通滤波
│   │   ├── resample.py             # 重采样
│   │   ├── snr.py                  # 信噪比计算与筛选
│   │   └── rotation.py             # RTZ 旋转
│   │
│   ├── deconvolution/              # 反褶积
│   │   ├── __init__.py
│   │   ├── iter_decon.py           # 迭代反褶积
│   │   └── water_level.py          # 水位反褶积
│   │
│   ├── array_processing/           # 台阵处理 (Phase 2 重点)
│   │   ├── __init__.py
│   │   ├── fk_filter.py            # FK 滤波
│   │   ├── radon2d.py              # 2D Radon 变换 [整合已有包]
│   │   ├── radon3d.py              # 3D Radon 变换
│   │   ├── rank_reduction_2d.py    # 2D 秩约简 [整合 Pydrr]
│   │   ├── rank_reduction_3d.py    # 3D 秩约简(DRR-OTG) [整合 Pydrr]
│   │   ├── soi_filter.py           # 结构导向滤波
│   │   └── amf.py                  # 自适应匹配滤波
│   │
│   ├── imaging/                    # 成像 (Phase 3 重点)
│   │   ├── __init__.py
│   │   ├── ccp/                    # CCP 叠加
│   │   │   ├── __init__.py
│   │   │   ├── raytracing.py       # 射线追踪
│   │   │   ├── conversion_points.py # 转换点计算
│   │   │   ├── time_depth.py       # 时深变换
│   │   │   ├── fresnel.py          # Fresnel 带加权
│   │   │   └── stacking.py         # 叠加计算
│   │   ├── lsm/                    # 最小二乘偏移 ⭐新开发
│   │   │   ├── __init__.py
│   │   │   ├── forward.py          # 正演算子
│   │   │   ├── adjoint.py          # 伴随算子
│   │   │   ├── solver.py           # 迭代求解器
│   │   │   └── regularization.py   # 正则化
│   │   ├── hk/                     # HK 叠加
│   │   │   ├── __init__.py
│   │   │   └── stacking.py         # HK 网格搜索
│   │   └── corrections.py          # 时深校正、地形校正
│   │
│   ├── visualization/              # 可视化
│   │   ├── __init__.py
│   │   ├── ccp.py                  # CCP 剖面/3D 图
│   │   ├── stations.py             # 台站/事件分布
│   │   ├── waveforms.py            # 波形图
│   │   ├── rf.py                   # 接收函数显示
│   │   ├── hk.py                   # HK 图
│   │   └── map_utils.py            # 底图工具
│   │
│   └── utilities/                  # 工具函数
│       ├── __init__.py
│       ├── coordinates.py          # 坐标转换 (latlon ↔ xy)
│       ├── geometry.py             # 几何计算 (距离、方位角)
│       ├── time_utils.py           # 时间处理
│       ├── stats.py                # 统计工具
│       └── signal.py               # 信号处理辅助
│
├── tests/                          # 测试
│   ├── test_core/
│   ├── test_io/
│   ├── test_preprocessing/
│   ├── test_deconvolution/
│   ├── test_array_processing/
│   ├── test_imaging/
│   └── validation/                 # MATLAB 交叉验证脚本
│
├── examples/                       # 示例与 Demo
│   ├── notebooks/                  # Jupyter Notebooks
│   │   ├── 01_read_and_preprocess.ipynb
│   │   ├── 02_deconvolution.ipynb
│   │   ├── 03_array_processing.ipynb
│   │   ├── 04_ccp_imaging.ipynb
│   │   ├── 05_lsm_imaging.ipynb
│   │   └── 06_hk_stacking.ipynb
│   └── scripts/                    # 命令行脚本
│
├── docs/                           # 文档
│   ├── user_guide.md
│   ├── api_reference.md
│   └── DenseArrayToolkit_Python_Framework.md  # 本文档
│
├── pyproject.toml                  # 项目配置 (PEP 621)
├── README.md
├── LICENSE
└── .github/
    └── workflows/
        └── tests.yml               # CI/CD
```

### 4.2 模块实现状态与策略

| 模块 | 实现策略 | 优先级 | 预估工期 |
|------|---------|--------|---------|
| `core/` (数据模型) | 全新开发 | 🔴 P0 | 2–3 周 |
| `io/` (数据读写) | 基于 ObsPy 封装 | 🔴 P0 | 1–2 周 |
| `preprocessing/` | 基于 ObsPy/Scipy 封装 | 🔴 P0 | 1–2 周 |
| `deconvolution/` | 参考 MATLAB 重写 | 🟡 P1 | 2 周 |
| `array_processing/fk_filter.py` | 参考 MATLAB 重写 | 🟡 P1 | 1 周 |
| `array_processing/rank_reduction_*.py` | **整合 Pydrr** | 🟡 P1 | 2 周 |
| `array_processing/radon2d.py` | **整合已有 PyRadon** | 🟡 P1 | 1 周 |
| `array_processing/radon3d.py` | 参考 MATLAB 重写 | 🟡 P1 | 3 周 |
| `array_processing/soi_filter.py` | 参考 MATLAB 重写 | 🟢 P2 | 2–3 周 |
| `array_processing/amf.py` | 参考 MATLAB 重写 | 🟢 P2 | 2–3 周 |
| `imaging/ccp/` | 参考 rf/seispy + MATLAB | 🟡 P1 | 4–5 周 |
| `imaging/lsm/` | ⭐ **完全新开发** | 🔴 P0 | 6–8 周 |
| `imaging/hk/` | 参考 rf/seispy + MATLAB | 🟡 P1 | 1–2 周 |
| `visualization/` | Matplotlib + Cartopy 新写 | 🟢 P2 | 3–4 周 |
| `utilities/` | 参考 MATLAB 重写 | 🟡 P1 | 2–3 周 |

---

## 5. 核心数据模型设计

### 5.1 设计原则

- 使用 `@dataclass` 定义数据结构，配合 `pydantic` 进行运行时校验
- 类型注解覆盖所有字段，支持 IDE 智能提示和 mypy 静态检查
- 数据类包含常见的运算方法和可视化方法
- 支持与 ObsPy 数据结构的双向转换

### 5.2 基础数据结构

```python
from dataclasses import dataclass, field
from typing import Optional, List, Tuple
import numpy as np
from obspy import UTCDateTime

@dataclass
class Station:
    """地震台站"""
    name: str
    network: str = ""
    latitude: float = 0.0
    longitude: float = 0.0
    elevation_m: float = 0.0
    # 投影坐标 (自动计算)
    projected_x: Optional[float] = None
    projected_y: Optional[float] = None

    def __repr__(self) -> str:
        return f"Station({self.name}, lon={self.longitude:.3f}, lat={self.latitude:.3f})"


@dataclass
class Event:
    """地震事件"""
    id: str
    origin_time: UTCDateTime
    latitude: float
    longitude: float
    depth_km: float
    magnitude: float = 0.0
    # 台站相对信息 (事件-台站配对后填充)
    distance_deg: Optional[float] = None
    back_azimuth: Optional[float] = None
    ray_param_s_deg: Optional[float] = None

    def __repr__(self) -> str:
        return f"Event({self.id}, M{self.magnitude:.1f}, {self.depth_km:.1f}km)"


@dataclass
class Trace:
    """单道三分量波形数据"""
    station: Station
    event: Event
    # 原始数据
    data: np.ndarray                     # [Nt x 3]: E, N, Z
    time_axis: np.ndarray                # [Nt]: 相对时间轴
    sampling_rate: float
    start_time: UTCDateTime
    # 预处理后
    data_rotated: Optional[np.ndarray] = None   # [Nt x 3]: R, T, Z
    data_filtered: Optional[np.ndarray] = None
    # 接收函数
    rf_radial: Optional[np.ndarray] = None      # 径向 RF
    rf_transverse: Optional[np.ndarray] = None  # 切向 RF
    rf_time: Optional[np.ndarray] = None         # RF 时间轴
    # 质量控制
    snr: Optional[float] = None
    quality_flag: bool = True
    # 处理历史
    history: List[str] = field(default_factory=list)

    @property
    def npts(self) -> int:
        return len(self.time_axis)

    def __repr__(self) -> str:
        return f"Trace({self.station.name}, {self.event.id})"


@dataclass
class CommonEventGather:
    """共事件点道集"""
    event: Event
    traces: List[Trace]
    # 属性
    n_stations: int = 0

    def __post_init__(self):
        self.n_stations = len(self.traces)

    @property
    def rf_matrix(self) -> np.ndarray:
        """返回 [Nt x Nstations] 的 RF 矩阵"""
        rfs = [t.rf_radial for t in self.traces if t.rf_radial is not None]
        return np.column_stack(rfs) if rfs else np.array([])

    def __repr__(self) -> str:
        return f"CommonEventGather({self.event.id}, {self.n_stations} stations)"


@dataclass
class CommonStationGather:
    """共台站道集"""
    station: Station
    traces: List[Trace]
    n_events: int = 0

    def __post_init__(self):
        self.n_events = len(self.traces)

    def __repr__(self) -> str:
        return f"CommonStationGather({self.station.name}, {self.n_events} events)"
```

### 5.3 速度模型 (继承与多态)

```python
from abc import ABC, abstractmethod

class VelocityModel(ABC):
    """速度模型抽象基类"""

    @abstractmethod
    def get_vp(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        ...

    @abstractmethod
    def get_vs(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        ...

    @abstractmethod
    def get_vp_vs_ratio(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        ...

    @property
    @abstractmethod
    def model_type(self) -> str:
        ...

    @abstractmethod
    def to_dict(self) -> dict:
        ...


class VelocityModel1D(VelocityModel):
    """1D 速度模型 (AK135, IASP91 等)"""

    def __init__(self, depths: np.ndarray, vp: np.ndarray, vs: np.ndarray,
                 name: str = "custom", rho: Optional[np.ndarray] = None):
        self.depths = np.asarray(depths)
        self.vp = np.asarray(vp)
        self.vs = np.asarray(vs)
        self.name = name
        self.rho = np.asarray(rho) if rho is not None else None

    def get_vp(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        return float(np.interp(depth, self.depths, self.vp))

    def get_vs(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        return float(np.interp(depth, self.depths, self.vs))

    def get_vp_vs_ratio(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        vp = self.get_vp(depth)
        vs = self.get_vs(depth)
        return vp / vs if vs > 0 else 1.73  # 默认泊松比

    @property
    def model_type(self) -> str:
        return "1D"

    @classmethod
    def from_ak135(cls) -> "VelocityModel1D":
        """从 AK135 模型创建"""
        ...

    @classmethod
    def from_file(cls, filepath: str) -> "VelocityModel1D":
        """从文件加载"""
        ...

    def to_dict(self) -> dict:
        return {"type": "1D", "name": self.name,
                "depths": self.depths.tolist(),
                "vp": self.vp.tolist(), "vs": self.vs.tolist()}

    def __repr__(self) -> str:
        return f"VelocityModel1D({self.name}, {len(self.depths)} layers)"


class VelocityModel3D(VelocityModel):
    """3D 速度模型 (网格化)"""

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray,
                 vp_3d: np.ndarray, vs_3d: np.ndarray, name: str = "custom_3d"):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.z = np.asarray(z)
        self.vp_3d = np.asarray(vp_3d)   # [nz, ny, nx]
        self.vs_3d = np.asarray(vs_3d)
        self.name = name

    def get_vp(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        # 3D 插值
        ...

    def get_vs(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        ...

    def get_vp_vs_ratio(self, depth: float, lat: float = 0, lon: float = 0) -> float:
        ...

    @property
    def model_type(self) -> str:
        return "3D"

    @classmethod
    def from_crust1(cls, region: Tuple[float, float, float, float]) -> "VelocityModel3D":
        """从 Crust1.0 全球模型提取区域"""
        ...

    def to_dict(self) -> dict:
        return {"type": "3D", "name": self.name,
                "x_range": [self.x[0], self.x[-1]],
                "y_range": [self.y[0], self.y[-1]],
                "z_range": [self.z[0], self.z[-1]]}

    def __repr__(self) -> str:
        return f"VelocityModel3D({self.name}, shape=({len(self.z)},{len(self.y)},{len(self.x)}))"
```

### 5.4 配置系统

```python
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class PreprocessingConfig:
    """预处理参数"""
    sampling_rate: float = 20.0
    bandpass_low: float = 0.05
    bandpass_high: float = 2.0
    time_before_p: float = 20.0
    time_after_p: float = 100.0
    quality_filter: bool = True
    snr_threshold: float = 3.0

@dataclass
class DeconvolutionConfig:
    """反褶积参数"""
    method: str = "iterative"        # "iterative" | "water_level"
    iterations: int = 400
    gauss_width: float = 2.5
    water_level: float = 0.01

@dataclass
class ArrayProcessingConfig:
    """台阵处理参数"""
    drr_rank: int = 10              # DRR 保留奇异值个数
    drr_k: float = 5.0              # 阻尼因子
    drr_niter: int = 20             # 最大迭代次数
    drr_tolerance: float = 1e-3
    drr_freq_min: float = 0.1
    drr_freq_max: float = 1.2
    radon_type: str = "parabolic"   # "linear" | "parabolic" | "hyperbolic"

@dataclass
class CCPConfig:
    """CCP 成像参数"""
    imaging_type: str = "2D"        # "2D" | "3D"
    stack_mode: str = "fresnel"     # "uniform" | "fresnel"
    gauss_factor: float = 1.0
    bin_size_km: float = 0.5
    max_depth_km: float = 80.0

@dataclass
class LSMConfig:
    """最小二乘偏移参数"""
    imaging_type: str = "2D"        # "2D" | "3D"
    max_iterations: int = 30
    tolerance: float = 1e-4
    regularization: str = "laplacian"  # "laplacian" | "tv" | None
    reg_weight: float = 0.01

@dataclass
class ProcessingConfig:
    """总配置"""
    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    deconv: DeconvolutionConfig = field(default_factory=DeconvolutionConfig)
    array_processing: ArrayProcessingConfig = field(default_factory=ArrayProcessingConfig)
    ccp: CCPConfig = field(default_factory=CCPConfig)
    lsm: LSMConfig = field(default_factory=LSMConfig)
    data_dir: str = "./data"
    output_dir: str = "./results"
    verbose: bool = True

    @classmethod
    def from_file(cls, filepath: str) -> "ProcessingConfig":
        """从 YAML/JSON 文件加载配置"""
        ...

    def save(self, filepath: str) -> None:
        """保存配置到文件"""
        ...

    def to_dict(self) -> dict:
        """转为字典 (用于序列化)"""
        ...
```

---

## 6. Pipeline API 设计

### 6.1 设计理念

采用 **Builder Pattern** (构建器模式) 实现链式调用，让代码具有高度可读性：

```python
# 用户可以一行链式调用完成整个流程
result = (SeismicPipeline(config)
    .read_sac("./data/")
    .preprocess()
    .deconvolve()
    .filter_by_quality()
    .to_common_event_gathers()
    .rank_reduction_3d()
    .ccp_imaging()
    .plot()
)
```

### 6.2 SeismicPipeline 核心类

```python
class SeismicPipeline:
    """主处理管道 —— 密集台阵接收函数处理全流程"""

    def __init__(self, config: ProcessingConfig):
        self.config = config
        self.traces: List[Trace] = []
        self.gathers: List[CommonEventGather] = []
        self._history: List[str] = []
        self._logger = logging.getLogger(__name__)

    # ==================== I/O ====================

    def read_sac(self, data_dir: str, station_file: Optional[str] = None) -> "SeismicPipeline":
        """读取 SAC 波形数据"""
        self.traces = sac_reader.read_directory(
            data_dir, station_file=station_file
        )
        self._log(f"Loaded {len(self.traces)} three-component traces")
        return self

    # ==================== 预处理 ====================

    def preprocess(self) -> "SeismicPipeline":
        """预处理: 滤波 → 重采样 → 截取 → 旋转 → 信噪比"""
        self.traces = preprocessing_pipeline(self.traces, self.config.preprocessing)
        self._log("Preprocessing complete")
        return self

    # ==================== 反褶积 ====================

    def deconvolve(self) -> "SeismicPipeline":
        """提取接收函数"""
        for trace in self.traces:
            rf_r, rf_t, rf_time = iter_deconvolve(
                trace.data_filtered,
                trace.sampling_rate,
                self.config.deconv,
            )
            trace.rf_radial = rf_r
            trace.rf_transverse = rf_t
            trace.rf_time = rf_time
        self._log("Deconvolution complete")
        return self

    # ==================== 质控 ====================

    def filter_by_quality(self, snr_threshold: Optional[float] = None) -> "SeismicPipeline":
        """按信噪比筛选"""
        threshold = snr_threshold or self.config.preprocessing.snr_threshold
        self.traces = [t for t in self.traces if (t.snr or 0) >= threshold]
        self._log(f"Kept {len(self.traces)} traces after QC (SNR >= {threshold})")
        return self

    # ==================== 道集 ====================

    def to_common_event_gathers(self) -> "SeismicPipeline":
        """构建共事件点道集"""
        self.gathers = build_common_event_gathers(self.traces)
        self._log(f"Built {len(self.gathers)} common-event gathers")
        return self

    # ==================== 台阵处理 ====================

    def rank_reduction_3d(self) -> "SeismicPipeline":
        """3D 阻尼秩约简 — 整合 Pydrr"""
        cfg = self.config.array_processing
        for gather in self.gathers:
            rf_mat = gather.rf_matrix            # [Nt x Nstations]
            station_x = [t.station.projected_x for t in gather.traces]
            station_y = [t.station.projected_y for t in gather.traces]

            # 调用 Pydrr 集成接口
            rf_denoised = pydrr_wrapper.drr_otg_3d(
                data=rf_mat,
                x=station_x, y=station_y,
                rank=cfg.drr_rank,
                K=cfg.drr_k,
                niter=cfg.drr_niter,
                eps=cfg.drr_tolerance,
                flow=cfg.drr_freq_min,
                fhigh=cfg.drr_freq_max,
                dt=1.0 / gather.traces[0].sampling_rate,
                verb=self.config.verbose,
                mode=1,  # denoising + reconstruction
            )
            gather._rf_processed = rf_denoised
        self._log("3D DRR-OTG rank reduction complete")
        return self

    def radon_transform_3d(self) -> "SeismicPipeline":
        """3D Radon 变换去噪"""
        ...

    # ==================== 成像 ====================

    def ccp_imaging(self, velocity_model: VelocityModel,
                    image_grid: ImageGrid) -> "SeismicPipeline":
        """CCP 叠加成像"""
        self._ccp_result = ccp_stacking(
            self.gathers, velocity_model, image_grid, self.config.ccp
        )
        self._log("CCP stacking complete")
        return self

    def lsm_imaging(self, velocity_model: VelocityModel,
                    image_grid: ImageGrid) -> "SeismicPipeline":
        """最小二乘偏移成像 ⭐新开发"""
        self._lsm_result = least_squares_migration(
            self.gathers, velocity_model, image_grid, self.config.lsm
        )
        self._log("LSM imaging complete")
        return self

    def hk_stacking(self, velocity_model_1d: VelocityModel1D) -> "SeismicPipeline":
        """HK 叠加"""
        ...

    # ==================== 可视化 ====================

    def plot_ccp(self, profile_idx: int = 0,
                 save_path: Optional[str] = None) -> None:
        """绘制 CCP 剖面"""
        ...

    def plot_results(self, save_dir: Optional[str] = None) -> "SeismicPipeline":
        """批量绘图"""
        ...

    # ==================== 结果输出 ====================

    def save_results(self, output_dir: Optional[str] = None) -> "SeismicPipeline":
        """保存处理结果"""
        ...

    @property
    def history(self) -> str:
        """处理历史摘要"""
        return "\n".join(f"  [{i+1}] {h}" for i, h in enumerate(self._history))

    def _log(self, message: str) -> None:
        """内部日志"""
        self._history.append(message)
        if self.config.verbose:
            self._logger.info(message)
```

---

## 7. 计算引擎层设计

### 7.1 核心原则：纯函数 + NumPy 数组

计算引擎层的所有函数遵循以下契约：

```python
def algorithm_name(
    data: np.ndarray,        # 输入数据
    *params: float,          # 标量参数
    **kwargs                  # 可选配置
) -> Union[np.ndarray, Tuple[np.ndarray, ...]]:
    """
    简短描述算法功能。

    Parameters
    ----------
    data : np.ndarray
        输入数据, shape (nt, nx) 或 (nt, nx, ny)
    ...

    Returns
    -------
    result : np.ndarray
        输出结果
    ...
    """
    # 纯计算逻辑，无副作用，无全局状态
    ...
```

### 7.2 关键算法接口

```python
# ---------- 秩约简 (整合 Pydrr) ----------

def rank_reduction_2d(
    data: np.ndarray,           # [nt, ntraces]
    rank: int = 10,
    K: float = 5.0,
    niter: int = 20,
    eps: float = 1e-3,
    dt: float = 0.1,
    flow: float = 0.1,
    fhigh: float = 1.2,
) -> np.ndarray:                # [nt, ntraces]
    """2D 阻尼秩约简 (DRR)"""
    from pydrr import drr2d
    return drr2d(data, rank=rank, K=K, niter=niter, eps=eps, dt=dt,
                 flow=flow, fhigh=fhigh)


def rank_reduction_3d(
    data: np.ndarray,           # [nt, ntraces]
    x: np.ndarray,              # [ntraces]
    y: np.ndarray,              # [ntraces]
    nx: int, ny: int,
    ox: float, oy: float,
    mx: float, my: float,
    rank: int = 10,
    K: float = 5.0,
    niter: int = 20,
    eps: float = 1e-3,
    dt: float = 0.1,
    flow: float = 0.1,
    fhigh: float = 1.2,
    mode: int = 1,              # 1=去噪+重建, 0=仅重建
) -> Tuple[np.ndarray, np.ndarray]:
    """3D 阻尼秩约简 (DRR-OTG) — 整合 Pydrr"""
    from pydrr import drr3d_otg
    return drr3d_otg(data, x, y, nx, ny, ox, oy, mx, my,
                     flow, fhigh, dt, rank, K, niter, eps, mode=mode)


# ---------- Radon 变换 ----------

def radon_transform_2d(
    data: np.ndarray,           # [nt, nx]
    p_min: float,
    p_max: float,
    n_p: int,
    radon_type: str = "parabolic",  # "linear" | "parabolic" | "hyperbolic"
    mu: float = 1.0,
    niter: int = 50,
    tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """2D Radon 变换 (稀疏约束 + 共轭梯度求解)"""
    ...


def radon_transform_3d(
    data: np.ndarray,           # [nt, nx, ny]
    dt: float,
    dx: float, dy: float,
    p_min: float, p_max: float, n_p: int,
    q_min: float, q_max: float, n_q: int,
    mu: float = 1.0,
    niter: int = 50,
    tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray]:
    """3D Radon 变换 (规则网格)"""
    ...


# ---------- FK 滤波 ----------

def fk_filter(
    data: np.ndarray,           # [nt, ntraces]
    dt: float,
    dx: float,
    v_min: float,
    v_max: float,
    taper_width: float = 0.1,
) -> np.ndarray:
    """FK 域锥形滤波"""
    ...


# ---------- 迭代反褶积 ----------

def iter_deconvolve(
    data_3c: np.ndarray,        # [nt, 3]   R, T, Z
    sampling_rate: float,
    config: DeconvolutionConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """迭代反褶积提取接收函数"""
    ...


# ---------- CCP 叠加 ----------

def ccp_stacking(
    gathers: List[CommonEventGather],
    velocity_model: VelocityModel,
    image_grid: "ImageGrid",
    config: CCPConfig,
) -> "CCPResult":
    """CCP 叠加成像 (支持 uniform / Fresnel 加权模式)"""
    ...


# ---------- 最小二乘偏移 ⭐核心新开发 ----------

def least_squares_migration(
    gathers: List[CommonEventGather],
    velocity_model: VelocityModel,
    image_grid: "ImageGrid",
    config: LSMConfig,
) -> "MigrationResult":
    """
    最小二乘偏移 (LSM) 成像。
    
    算法框架:
        m_{k+1} = m_k - α * L^T (L m_k - d) + λ R(m_k)
    
    其中:
        L   = 正演 (Kirchhoff) 算子
        L^T = 伴随 (偏移) 算子
        R   = 正则化项 (Laplacian / TV)
        α   = 步长 (线搜索或固定)
        λ   = 正则化权重
    """
    ...


# ---------- HK 叠加 ----------

def hk_stacking(
    rf_radial_list: List[np.ndarray],
    rf_time: np.ndarray,
    ray_params: np.ndarray,
    vp: float = 6.3,
    h_range: Tuple[float, float] = (20.0, 60.0),
    k_range: Tuple[float, float] = (1.5, 2.0),
    h_step: float = 0.5,
    k_step: float = 0.01,
    weight: Tuple[float, float, float] = (0.34, 0.33, 0.33),
) -> "HKResult":
    """H-κ 叠加"""
    ...
```

---

## 8. 开发路线图

### 8.1 阶段划分

```
               Q3 2026          Q4 2026          Q1 2027          Q2 2027
Phase 1  ████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  基础框架
Phase 2  ░░░░░░░░████████████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  台阵处理
Phase 3  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████████████████████░░░░░░  成像模块
Phase 4  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████  集成+文档
```

### 8.2 详细里程碑

#### Phase 1: 基础框架 (8–10 周)
- [ ] 项目脚手架：`pyproject.toml`、目录结构、CI/CD
- [ ] `core/` 数据模型：Trace, Station, Event, Gather, VelocityModel
- [ ] `io/` SAC 读写 (基于 ObsPy)
- [ ] `preprocessing/` 预处理管道
- [ ] `deconvolution/` 迭代反褶积
- [ ] `utilities/` 坐标转换、几何计算
- [ ] `config.py` 配置系统 (YAML 支持)
- [ ] 单元测试 + 与 MATLAB 输出交叉验证
- [ ] Demo Notebook 01 + 02

#### Phase 2: 台阵处理 (6–8 周) ⬇️ 工作量降低
- [ ] `fk_filter.py` FK 滤波
- [ ] `rank_reduction_2d.py` + `rank_reduction_3d.py` (整合 Pydrr)
- [ ] `radon2d.py` (整合已有 PyRadon)
- [ ] `radon3d.py` (参考 MATLAB 重写)
- [ ] `soi_filter.py` 结构导向滤波
- [ ] `amf.py` 自适应匹配滤波
- [ ] 交叉验证
- [ ] Demo Notebook 03

#### Phase 3: 成像模块 (10–14 周) ⭐ LSM 核心开发
- [ ] `imaging/ccp/` CCP 叠加 (2D/3D, uniform + Fresnel)
- [ ] `imaging/lsm/` 最小二乘偏移 ⭐
  - [ ] 正演算子 (Kirchhoff)
  - [ ] 伴随算子
  - [ ] 共轭梯度求解器
  - [ ] 正则化 (Laplacian, TV)
  - [ ] 2D + 3D 支持
- [ ] `imaging/hk/` HK 叠加
- [ ] `imaging/corrections.py` 时深校正
- [ ] 交叉验证 (CCP, HK 可与 rf/seispy 对比；LSM 与 MATLAB 对比)
- [ ] Demo Notebook 04 + 05 + 06

#### Phase 4: 集成、文档与发布 (6–8 周)
- [ ] `visualization/` 全部绘图功能
- [ ] 性能优化 (Numba JIT / 并行)
- [ ] 完整 API 文档 (Sphinx + ReadTheDocs)
- [ ] 用户手册更新
- [ ] 示例数据集 + 完整 Demo
- [ ] PyPI 发布
- [ ] Docker 镜像
- [ ] 教程视频 / 培训材料

---

## 9. 工作量评估（修订版）

> **相比初版评估大幅降低**：由于 Pydrr 和 PyRadon 等已有 Python 包可直接整合，Phase 2 工作量降低约 40%。

| 阶段 | 内容 | 预估人月 | 难度 | 与初版相比 |
|------|------|---------|------|-----------|
| Phase 1 | 基础框架 | **1.5–2** | ⭐⭐ | 持平 |
| Phase 2 | 台阵处理 | **2–2.5** | ⭐⭐⭐ | **↓ 25–37%** |
| Phase 3 | 成像模块 | **3–4** | ⭐⭐⭐⭐ | 持平 (LSM 需全新开发) |
| Phase 4 | 集成与文档 | **1.5–2** | ⭐⭐ | 持平 |
| **合计** | | **8–10.5 人月** | | **↓ 约 25%** |

### 团队配置建议

| 角色 | 人数 | 职责 |
|------|------|------|
| 核心算法开发 | 1–2 人 | Phase 2 (台阵处理) + Phase 3 (LSM 偏移) |
| 数据处理与可视化 | 1 人 | Phase 1 (基础框架) + Phase 4 (可视化) |
| 测试与文档 | 0.5 人 | 交叉验证 + API 文档 |

**推荐时间线**: 2–3 人团队，**5–7 个月**可完成 beta 版本。

---

## 10. 测试与验证策略

### 10.1 测试金字塔

```
         ┌──────────────┐
         │   E2E Tests  │  ← 完整管道测试 (少量)
         │ (demo scripts)│
        ┌┴──────────────┴┐
        │  Integration   │  ← 管道集成测试
        │    Tests       │
       ┌┴────────────────┴┐
       │   Unit Tests     │  ← 每个函数独立测试
       │  (pytest)        │
      ┌┴──────────────────┴┐
      │  MATLAB Validation │  ← 与 MATLAB 输出对比
      │  (交叉验证套件)     │
      └────────────────────┘
```

### 10.2 交叉验证策略

对于每个核心算法，使用相同输入数据对比 MATLAB 和 Python 输出：

```python
# tests/validation/test_rank_reduction_3d.py
import numpy as np
import h5py
import pytest
from densearray.array_processing import rank_reduction_3d

def test_rank_reduction_3d_against_matlab():
    """验证 3D DRR-OTG 与 MATLAB 输出一致"""
    # 加载 MATLAB 验证数据集
    with h5py.File("tests/validation/data/drr3d_validation.h5", "r") as f:
        input_data = f["input"][:]
        matlab_output = f["output_matlab"][:]

    # Python 计算
    python_output = rank_reduction_3d(
        input_data, x=f["x"][:], y=f["y"][:],
        nx=20, ny=20, ox=0, oy=0, mx=100, my=100,
        rank=10, K=5.0, niter=20, eps=1e-3,
        dt=0.1, flow=0.1, fhigh=1.2
    )[0]

    # 验证: 相对误差 < 1%
    rel_error = np.linalg.norm(python_output - matlab_output) / np.linalg.norm(matlab_output)
    assert rel_error < 0.01, f"Relative error {rel_error:.4f} exceeds 1%"
```

### 10.3 验证清单

| 算法 | 验证方法 | 容差 |
|------|---------|------|
| 迭代反褶积 | 与 MATLAB 输出逐点对比 | < 0.1% |
| DRR 2D/3D | 与 MATLAB + Pydrr 原生对比 | < 1% |
| Radon 2D/3D | 与 MATLAB 输出对比 | < 1% |
| FK 滤波 | 与 MATLAB 输出对比 | < 0.5% |
| CCP 叠加 | 与 MATLAB 输出对比 | < 2% |
| LSM 偏移 | 与 MATLAB 输出对比 | < 2% |
| HK 叠加 | 与 MATLAB + rf 包对比 | < 0.1% |

---

## 11. 项目文件结构

```
DenseArrayToolkit-Python/
│
├── pyproject.toml                    # PEP 621 项目配置
├── README.md                         # 项目介绍
├── LICENSE                           # GPL v3
├── CHANGELOG.md
├── .gitignore
├── .pre-commit-config.yaml           # pre-commit hooks
│
├── densearray/                       # 主包
│   ├── __init__.py
│   ├── _version.py                   # 版本号
│   ├── core/                         # 数据模型
│   ├── io/                           # 数据读写
│   ├── preprocessing/               # 预处理
│   ├── deconvolution/               # 反褶积
│   ├── array_processing/            # 台阵处理
│   ├── imaging/                      # 成像
│   │   ├── ccp/
│   │   ├── lsm/                      # ⭐ 新开发
│   │   └── hk/
│   ├── visualization/               # 可视化
│   └── utilities/                   # 工具函数
│
├── tests/                            # 测试
│   ├── conftest.py                   # pytest 配置
│   ├── test_core/
│   ├── test_io/
│   ├── test_preprocessing/
│   ├── test_deconvolution/
│   ├── test_array_processing/
│   ├── test_imaging/
│   │   ├── test_ccp.py
│   │   ├── test_lsm.py
│   │   └── test_hk.py
│   └── validation/                   # MATLAB 交叉验证
│       ├── data/                     # 验证数据集 (.h5)
│       ├── generate_validation_data.m # MATLAB 生成验证数据脚本
│       └── test_against_matlab.py
│
├── examples/                         # 示例
│   ├── notebooks/
│   │   ├── 01_read_and_preprocess.ipynb
│   │   ├── 02_deconvolution.ipynb
│   │   ├── 03_array_processing.ipynb
│   │   ├── 04_ccp_imaging.ipynb
│   │   ├── 05_lsm_imaging.ipynb
│   │   └── 06_hk_stacking.ipynb
│   ├── scripts/
│   └── data/                         # 示例数据
│
├── docs/                             # 文档
│   ├── index.md
│   ├── installation.md
│   ├── user_guide/
│   │   ├── quickstart.md
│   │   ├── preprocessing.md
│   │   ├── array_processing.md
│   │   ├── ccp_imaging.md
│   │   ├── lsm_imaging.md
│   │   └── hk_stacking.md
│   ├── api_reference/
│   ├── DenseArrayToolkit_Python_Framework.md  # 本文档
│   └── migration_guide.md            # MATLAB → Python 迁移指南
│
└── .github/
    ├── workflows/
    │   ├── tests.yml                 # 自动测试
    │   ├── docs.yml                  # 文档构建
    │   └── publish.yml               # PyPI 发布
    └── ISSUE_TEMPLATE/
```

---

## 12. 依赖清单

### 12.1 核心依赖

```toml
# pyproject.toml
[project]
name = "densearray-toolkit"
version = "0.1.0"
description = "Python toolkit for dense array receiver function processing and imaging"
requires-python = ">=3.10"
license = {text = "GPL-3.0-or-later"}

dependencies = [
    # 计算核心
    "numpy>=1.24",
    "scipy>=1.10",

    # 地震学 (ObsPy 生态)
    "obspy>=1.4",

    # 坐标与地图
    "pyproj>=3.5",
    "cartopy>=0.21",

    # 可视化
    "matplotlib>=3.7",

    # 数据处理
    "h5py>=3.8",              # 验证数据集
    "pyyaml>=6.0",            # 配置文件
    "tqdm>=4.65",             # 进度条

    # 数据校验
    "pydantic>=2.0",
]

[project.optional-dependencies]
# 台阵处理集成
array = [
    "pydrr>=1.0",             # 阻尼秩约简 (DRR)
    "pyradon>=0.1",           # Radon 变换 (如已发布)
]

# 性能加速
performance = [
    "numba>=0.57",
    "dask>=2023.6",
]

# 开发工具
dev = [
    "pytest>=7.0",
    "pytest-cov>=4.0",
    "mypy>=1.0",
    "ruff>=0.0.280",
    "pre-commit>=3.3",
    "sphinx>=7.0",
    "jupyter>=1.0",
]
```

### 12.2 与 MATLAB 工具箱的等价映射

| MATLAB 功能 | Python 等价 | 状态 |
|------------|------------|------|
| SAC 读写 (MatSAC) | ObsPy `obspy.io.sac` | ✅ |
| Signal Processing Toolbox | SciPy `scipy.signal` | ✅ |
| Mapping Toolbox | Cartopy + pyproj | ✅ |
| `fft()` / `ifft()` | NumPy `numpy.fft` | ✅ |
| `svds()` | SciPy `scipy.sparse.linalg.svds` | ✅ |
| `interp2()` / `interp1()` | SciPy `scipy.interpolate` | ✅ |
| `hankel()` | 自实现 (或 `scipy.linalg.hankel`) | ✅ |
| `struct` | `@dataclass` / `pydantic.BaseModel` | ✅ |
| `nargin` 默认参数 | Python 函数默认参数 | ✅ |
| DRR-OTG (MATdrr) | **Pydrr** | ✅ 已有 |
| Radon 2D (radon_op.m) | **PyRadon** (或自实现) | ✅ 已有 |
| LSM 偏移 | **全新开发** | ⬜ 待开发 |

---

## 13. 编码规范与贡献指南

### 13.1 代码风格

- **格式化**: Ruff (替代 Black + isort + flake8)
- **类型检查**: mypy (strict mode)
- **文档字符串**: NumPy docstring 格式
- **行宽**: 100 字符
- **命名**: PEP 8 (snake_case for functions/variables, PascalCase for classes)

### 13.2 示例代码风格

```python
"""Module docstring: brief description."""

from typing import Optional, Tuple

import numpy as np


def process_rf_data(
    rf_matrix: np.ndarray,
    sampling_rate: float,
    rank: int = 10,
    *,
    verbose: bool = False,
) -> Tuple[np.ndarray, dict]:
    """
    Process receiver function data using damped rank reduction.

    Parameters
    ----------
    rf_matrix : np.ndarray
        Input RF data matrix of shape (n_times, n_traces).
    sampling_rate : float
        Sampling rate in Hz.
    rank : int, optional
        Number of singular values to preserve. Default is 10.
    verbose : bool, optional
        If True, print progress information. Default is False.

    Returns
    -------
    rf_denoised : np.ndarray
        Denoised RF data, same shape as input.
    info : dict
        Dictionary containing processing metadata.

    Raises
    ------
    ValueError
        If rf_matrix is not 2D or rank exceeds matrix dimensions.

    Notes
    -----
    This function wraps the Pydrr implementation of the damped rank
    reduction (DRR) method [Chen et al., 2016].

    References
    ----------
    .. [1] Chen, Y., et al. (2016). Computers & Geosciences, 95, 59-66.
    """
    if rf_matrix.ndim != 2:
        raise ValueError(f"rf_matrix must be 2D, got shape {rf_matrix.shape}")

    # ... implementation ...
```

### 13.3 贡献流程

1. Fork 仓库 → 创建 feature 分支
2. 编写代码 + 测试 + 文档
3. 运行 `pre-commit run --all-files`
4. 提交 PR，描述变更内容
5. CI 自动运行测试 + 类型检查
6. 至少一位 reviewer 审核通过后合并

---

## 附录 A: 与现有 Python 包的关系

| 包名 | 功能 | 与本工具包关系 |
|------|------|--------------|
| **ObsPy** | 地震波形处理、格式转换、走时计算 | 基础依赖，用于 I/O 和通用处理 |
| **rf** (rf.readthedocs.io) | 接收函数计算 + CCP + HK | 参考实现，本工具包专注密集台阵增强 |
| **seispy** | 接收函数 + CCP + HK | 参考实现，本工具包差异化：DRR/Radon/LSM |
| **Pydrr** | 阻尼秩约简 DRR | **直接整合，作为核心依赖** |
| **PyRadon** | Radon 变换 | **直接整合，作为核心依赖** |
| **PyGMT** | GMT 绘图封装 | 可选依赖，用于高精度地图 |
| **pyfk** | FK 计算 | 参考实现 |

## 附录 B: 与 MATLAB 版本的 API 映射

| MATLAB 函数 | Python 等价 | 说明 |
|-------------|-----------|------|
| `read_SAC.m` | `densearray.io.sac.read_sac()` | 基于 ObsPy |
| `preprocessing.m` | `pipeline.preprocess()` | 管道方法 |
| `deconv.m` | `pipeline.deconvolve()` | 管道方法 |
| `rankReduction3D.m` | `densearray.array_processing.rank_reduction_3d()` | 整合 Pydrr |
| `radonTransform3D.m` | `densearray.array_processing.radon_transform_3d()` | 整合 PyRadon |
| `CCPCommonEventGather.m` | `densearray.imaging.ccp.stacking()` | 新实现 |
| `CCPStacking.m` | `densearray.imaging.ccp.stacking()` | 新实现 |
| `rf_ccp.m` | `densearray.imaging.ccp.raytracing()` | 新实现 |
| `hk_stacking.m` (demo) | `densearray.imaging.hk.stacking()` | 新实现 |
| `plotCCPXsection.m` | `ccp_result.plot_section()` | 数据类方法 |
| `getCommonEventGather.m` | `pipeline.to_common_event_gathers()` | 管道方法 |
| `skm2srad.m` | `densearray.utilities.geometry.km_to_radian()` | 工具函数 |

---

> **文档维护**: 本文档将随项目开发进展持续更新。  
> **反馈与贡献**: 欢迎通过 GitHub Issues 提出建议。  
> **许可证**: 本文档与项目代码共同遵循 GPL v3 协议。