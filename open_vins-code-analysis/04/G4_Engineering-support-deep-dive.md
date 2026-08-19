# 工程支撑工具集 精读报告 (G4)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/{utils,sim,plot}/`、`ov_eval/src/`、`ov_data/`
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 G4（Phase 3 函数级路线，阶段 G 工程支撑，P2-2 / P3-3/4/5）

> 本篇为**工程支撑层汇总**（了解接口/用途即可，不参与在线算法）。覆盖工具类、仿真器、绘图、评估工具集与数据集。

---

## Section 1: `ov_core/src/utils/` —— 通用工具

| 文件 | 主要接口 | 说明 |
|------|----------|------|
| `opencv_yaml_parse.h` | `YamlParser` 类：`parse_config` / `parse_external`（标量、矩阵、SE(3)；ROS 参数优先） | 参数读取（`readParam` 主实现） |
| `quat_ops.h` | `rot_2_quat`、`quat_2_Rot`、`quat_multiply`、`exp_se3`/`log_se3`、`Inv_se3`、`Jl_so3`/`Jr_so3`、`rot2rpy` | JPL 四元数/旋转工具（基于 Trawny 2005，全仓核心数学基础，见 S2） |
| `sensor_data.h` | `ImuData(t,w,a)`、`CameraData`（多相机图像+masks）及排序运算符 | 传感器数据结构 |
| `dataset_reader.h` | `DatasetReader`：`load_gt_file`（ASL/EuRoC 真值）、`get_gt_state`（17 维状态）、`load_simulated_trajectory` | 真值加载 |
| `print.h/.cpp` + `colors.h` | 带颜色日志 | 日志输出 |
| `opencv_lambda_body.h` | OpenCV 函数对象宏 | — |

> 注：仓库内**无独立 `math_helper` 文件**，数学工具集中在 `quat_ops.h`。

## Section 2: `ov_core/src/sim/` —— 仿真器

- `BsplineSE3.h/.cpp`：`BsplineSE3` 类，在 SE(3) 流形上做累积三次 B 样条插值，由轨迹点生成 C² 连续轨迹（可求速度/加速度），用于合成 IMU/相机数据。

## Section 3: `ov_core/src/plot/` —— 绘图工具

- `matplotlibcpp.h`：C++ 封装 matplotlib 的轻量头文件，供评估模块在 C++ 中直接画图（`plot_3errors` 等）。

## Section 4: `ov_eval/src/` —— 评估工具集

| 类 | 作用 |
|----|------|
| `AlignTrajectory` | SE3/SIM3/位置+偏航对齐 |
| `ResultTrajectory` | ATE/RPE/NEES 误差 |
| `ResultSimulation` | 全状态误差 + 3σ 界 |
| `Loader` / `Statistics` | 数据加载 / 统计 |

**关键可执行程序**：
- `error_singlerun` —— 单次运行误差
- `error_dataset` —— 数据集汇总
- `error_simulation` —— 仿真误差
- `live_align_trajectory` —— 实时对齐
- `plot_trajectories` —— 画轨迹
- `format_converter` / `pose_to_file` —— 格式转换
- `timing_*`（percentages/histogram/flamegraph）—— 耗时分析

## Section 5: `ov_data/` —— 数据集（非代码）

多数据集目录（euroc_mav、kaist、tum_vi、uzhfpv_*、rpng_*、sim 等），每个含 `.txt`/`.csv` 真值位姿与（部分）标定/配置数据，供测试与评估使用。

## Section 6: 在 pipeline 中的位置

- `utils`：被参数读取、数学运算（S2）、数据结构中转全面调用。
- `sim`：离线数据合成，不参与在线 pipeline。
- `plot`/`ov_eval`：离线评估，不参与在线。
- `ov_data`：测试输入，非代码。

## Section 7: 待详细补充项

- **G4a — `quat_ops.h` 数学函数**：与 S2（JPLQuat）对照的底层实现。
- **G4b — `BsplineSE3` 轨迹生成**：仿真器插值细节。
- **G4c — `ov_eval` 误差指标**：ATE/RPE/NEES 计算公式。

> 工程支撑层，按需查。遵循 SKILL.md "04/ 文档拆分规则"，子文档待需要时补。
