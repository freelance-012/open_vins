# OpenVINS 官方文档 — Zero Velocity Update (ZUPT 零速更新)

> 来源: https://docs.openvins.com/update-zerovelocity.html
> 缓存日期: 2026-08-19

## 概述

零速更新（ZUPT）的核心思想是：**利用系统处于静止这一运动知识来降低不确定性**。在以下场景尤其重要：
- 单目系统且无时序 SLAM 特征时，静止无法三角化特征，系统无法被更新；
- 典型自动驾驶场景：传感器在红灯处静止，而路口动态物体（其他车辆）会迅速污染系统。ZUPT + 跳过特征跟踪可应对。

## 恒定速度合成测量 (Constant Velocity Synthetic Measurement)

更新时构造一个合成"测量"：当前**真实**加速度与角速度均为零。相比直接说速度为零，这种方法能用 IMU 读数本身建模测量不确定性。

> 注意：这并非严格强制零速度，而是"恒定速度"约束。恒速（零加速度）时可能误检，需用速度幅值检查（velocity magnitude check）规避。

残差定义（合成测量减去测量函数），并对状态求雅可比。

## 零速检测 (Zero Velocity Detection)

检测本身是难题，已有大量工作（见下方引用）：
- 多数方法归结为简单阈值法，目标是找到最优阈值以最好地分类 ZUPT 段；
- 也有更复杂方法，处理阈值依赖运动类型（跑 vs 走）和平台特性（需忽略车辆引擎振动等非本质振动）。

### IMU 判据 (Inertial-based Detection)
基于上述测量模型调优 **卡方（chi-squared）阈值** 做检测；同时做速度幅值检查以防恒速非零误检。实际中常需将 IMU 噪声放大数倍才能获得正确检测（暗示 IMU 噪声被高估，或存在电机振动等额外频率注入噪声）。

### 视觉判据 (Disparity-based Detection)
利用连续两帧特征跟踪的视差：若特征跟踪视差变化极小，则判定静止。计算平均视差并对该值阈值化。在动态环境中可能失效，因此动态场景下建议 IMU 与视差方法联合使用。

## 对应代码模块
- `UpdaterZeroVelocity`（`ov_msckf/src/update/UpdaterZeroVelocity.*`）
- `VioManager` 中 `updaterZUPT->try_update(...)` 分支（见数据流 P2 报告 §7 #5）
- 配置项：`use_zupt`、`zupt_only_at_beginning`、`zupt_max_velocity`、`zupt_noise_multiplier`、视差检测相关参数

## 关联引用（OpenVINS citelist）
- [41] Wagstaff et al., 2017 IPIN — *Zero Velocity Detection: An Approach Comparing RNN and Traditional Signal Processing*
- [36] Ramanandan et al., 2011 TITS — 车载惯性静止检测
- [9] Davidson et al., 2009 — 惯性导航中静止段的检测与估计

> 注：以上三篇为 ZUPT 检测方法的参考文献，多为会议/期刊论文，未在本次收集中原地获取 PDF（部分付费），列为待补充（见 00_资料清单.md §E）。
