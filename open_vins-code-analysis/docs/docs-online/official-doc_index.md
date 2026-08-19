# OpenVINS 官方文档首页（缓存）

> 来源: https://docs.openvins.com/index.html
> 缓存日期: 2026-08-19

## 项目简介
OpenVINS 是特拉华大学 RPNG 团队开发的**滤波式视觉惯性估计研究平台**。核心滤波器为 **EKF**，融合惯性信息与稀疏视觉特征跟踪；采用 **MSCKF 滑动窗口** 公式（3D 特征不进状态向量，仅通过多帧约束更新位姿）。借鉴图优化思想，设计了模块化的 **type-based 状态系统** 以方便协方差管理。

- GitHub: https://github.com/rpng/open_vins
- 文档: https://docs.openvins.com/
- 主论文: https://pgeneva.com/downloads/papers/Geneva2020ICRA.pdf

## Project Features（特性清单）
- 滑动窗口视觉惯性 MSCKF
- 模块化协方差 type 系统
- 完备的推导文档
- 可扩展视觉惯性仿真器（SE(3) B 样条、任意相机数、任意传感器率、自动特征生成）
- **6 种特征参数化**：Global XYZ / Global inverse depth / Anchored XYZ / Anchored inverse depth / Anchored MSCKF inverse depth / Anchored single inverse depth
- 传感器内外参在线标定：相机-IMU 外参、相机-IMU 时间偏移、相机内参、IMU 内参（含 g-sensitivity）
- 环境 SLAM 特征：OpenCV ARUCO tag SLAM 特征 + 稀疏特征 SLAM 特征
- 视觉跟踪：单目 / 双目(同步) / 双相机(同步) / KLT 或描述子 / masked tracking
- **静态 + 动态状态初始化**
- **零速检测与更新 (ZUPT)**
- 开箱即用的 EuRoC / TUM-VI / UZH-FPV / KAIST 评估
- 完整评估套件（ATE, RPE, NEES, RMSE）

## Codebase Extensions（扩展）
- **ov_plane**：单目平面辅助 VIO
- **vicon2gt**：动捕真值轨迹生成
- **ov_maplab**：接入 maplab
- **ov_secondary**：基于 VINS-Fusion 的**松耦合**回环（不反哺估计器）

## License
GPL-3.0。引用见 `00_资料清单.md` §A 的 BibTeX（Geneva2020ICRA）。
