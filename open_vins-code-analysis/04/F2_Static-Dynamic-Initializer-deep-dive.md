# Static/Dynamic Initializer 精读报告 (F2)

> **仓库**: open_vins
> **模块路径**: `ov_init/src/static/StaticInitializer.cpp`、 `ov_init/src/dynamic/DynamicInitializer.cpp`
> **对应论文**: 静态初始化（静止启动，加速度均值对齐重力，gram_schmidt 求姿态 [11]）；动态初始化（Dong-Si & Mourikis 2012，CPI 连续预积分 + Ceres MLE 对齐 [11]）
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 F2（Phase 3 函数级路线，阶段 F 初始化，算法实现层）

> 上游：F1（`InertialInitializer` 调度本层）。本篇覆盖 **`StaticInitializer::initialize`（:37）** 与 **`DynamicInitializer::initialize`（:44）** 两种算法的内部实现。动态初始化是重点——它在内部构造 `ov_core::CpiV1` 预积分器（:243-307）并经 Ceres 因子（Factor_ImuCPIv1 等）做 MLE 对齐，详见 F3。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
InertialInitializer (F1)
├─ 静态 → StaticInitializer::initialize(timestamp, cov, order, t_imu, wait_for_jerk)  @ :37
│          假设: 载体静止启动, 仅 IMU
│          输出: 单个 IMU 状态 (q/p/v/bg/ba), 速度=0, yaw/位置不可观测固定
│
└─ 动态 → DynamicInitializer::initialize(timestamp, cov, order, _imu, _clones_IMU, _features_SLAM)  @ :44
           假设: 任意运动, IMU + 视觉特征库
           输出: IMU 状态 + IMU clones + SLAM 特征 + 协方差
           内部: 构造 CpiV1 (:243-307) → 线性系统求 v/g (:309-520) → Ceres MLE (:580-908) → F3
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 6 环节（初始化算法实现） |
| 输入 | 静态：仅 `imu_data`；动态：`imu_data` + `_db`（特征库） |
| 输出 | 初始状态 + 协方差；动态额外出 clones/SLAM 特征 |
| 调用方 | `InertialInitializer`（F1） |
| 被调用方 | `CpiV1`（F3）、`Factor_ImuCPIv1`/`Factor_ImageReprojCalib`/`Factor_GenericPrior`（ov_init/src/ceres/）、`InitializerHelper`（gram_schmidt/select_imu_readings） |

---

## Section 2: 静态初始化 (StaticInitializer::initialize :37)

### 2.1 步骤（带行号）

```cpp
:41-43   数据量检查 (窗口足够)
:46-53   取最新/最旧时间戳, 要求窗口 ≥ init_window_time
:57-64   切出两半窗 window_1to0 (最近) / window_2to1 (前一半)
:73-97   计算两窗加速度均值 a_avg 与样本标准差 a_var
:101-119 静止判定: 若 wait_for_jerk, 需新窗激励低且旧窗也低; 否则整体静止
:122-125 重力方向/姿态: z_axis = a_avg_2to1 归一化, InitializerHelper::gram_schmidt 构造 Ro (z 对齐 -g), q_GtoI = rot_2_quat(Ro)
         —— 只估计 roll/pitch, yaw 不可观测
:128-131 bias: gravity_inG=(0,0,gravity_mag); bg = w_avg_2to1; ba = a_avg_2to1 - R_GtoI*gravity_inG
:135-141 速度=0: 16维 imu_state 仅填 q/bg/ba, 位置/速度块为 0 (单目无尺度未知数, IMU 提供绝对尺度)
:144-146 初始协方差 (q:0.02², p:0.05², v:0.01² 静态)
```

### 2.2 算法要点

- **姿态**：静止时比力 `a ≈ R·g`，故加速度均值方向即重力反方向；用 gram_schmidt 把该方向对齐到 `-g`，得到 `R`（含 roll/pitch，yaw 自由）。
- **Bias**：陀螺零偏 = 角速度均值；加计零偏 = 观测加速度均值减去重力在机体系的分量。
- **速度/位置**：静止假设 → 速度恒 0；位置与 yaw 是 VIO 4-DoF 不可观测量，初始化时固定（位置=0，yaw 任意）。

---

## Section 3: 动态初始化 (DynamicInitializer::initialize :44)

### 3.1 前处理 (`:50-193`)

```cpp
:50-76   取 newest_cam_time, 清窗
:113-154 按 pose_dt_avg 选相机帧, 统计有效特征 (min_valid_features=8, min_num_meas_to_optimize)
:169-193 用 launch 给定 bias 初值, 检查角运动 ≥ init_dyn_min_deg (默认 45°)
```

### 3.2 构造 CpiV1 预积分器 (`:243-307`)  ★ 核心

遍历每个相机时间 `map_camera_times`：

```cpp
:248-253  首帧不预积分 (nullptr)
// I0→Ii (供线性系统)
:258  auto cpiI0toIi1 = std::make_shared<ov_core::CpiV1>(params.sigma_w, sigma_wb, sigma_a, sigma_ab, true);
:259  cpiI0toIi1->setLinearizationPoints(gyroscope_bias, accelerometer_bias);
:260-261  auto imu_sel = InitializerHelper::select_imu_readings(*imu_data, t0+calib_camimu_dt, t1+calib_camimu_dt);
:273-277  for 逐段: cpiI0toIi1->feed_IMU(t0, t1, wm0, am0, wm1, am1);
// Ii→Ii+1 (供 MLE)
:280-301 同样模式构造 cpiIitoIi1
:304-305 存入 map_camera_cpi_I0toIi / map_camera_cpi_IitoIi1
```
- **两组 CPI**：`I0→Ii`（相对首帧，用于线性系统求特征/速度/重力）与 `Ii→Ii+1`（逐帧，用于 Ceres MLE 残差）。
- 构造模式统一：`make_shared<CpiV1>(...)` + `setLinearizationPoints(bg, ba)` + 循环 `feed_IMU`（详见 F3）。
- 注意 `imu_avg=true`（:258/:282），即端点平均。

### 3.3 线性系统求 v/g (`:309-520`)

```cpp
:309-384  用 H_proj=[1 0 -u; 0 1 -v] 与 CPI 的 DT/R_k2tau/alpha_tau 构造 A x = b
          状态序 [features, v_I0, gravity]  (:365-376)
:395-472  重力模长约束 |g|=9.81: compute_dongsi_coeff 配分矩阵/特征值法 (最小实根) 解 state_grav / state_feat_vel
:479     校验 |g| 收敛
:491-520 由 CPI 反推各帧 p/v
:556-575 gram_schmidt(gravity_inI0) 对齐到全局重力系
```
- 先用**线性系统**快速得到特征位置、初始速度、重力方向（解析、快）。
- 重力方向由 `|g|=9.81` 约束从观测中解出（Dong-Si & Mourikis 2012 方法）。

### 3.4 Ceres MLE 对齐 (`:580-908`)

**因子类**（文件 `ov_init/src/ceres/`）：

| 因子 | 文件 | 作用 |
|------|------|------|
| `Factor_ImuCPIv1` | `Factor_ImuCPIv1.h:32` | CPI 预积分残差（误差=预积分测量 vs 当前状态，含 bias 线性化点） |
| `Factor_ImageReprojCalib` | `Factor_ImageReprojCalib.h:41` | 带标定的重投影残差（pose+feature+外参+内参） |
| `Factor_GenericPrior` | `Factor_GenericPrior.h:29` | 首帧固定先验（yaw+位置+零偏），消除 4-DoF 不可观测 |
| `State_JPLQuatLocal` | — | 四元数局部参数化 |

**调用**：
```cpp
:757-759  new Factor_ImuCPIv1(...) + AddResidualBlock (连接 k,k+1 帧 ori/pos/vel/bg/ba)
:893      Factor_ImageReprojCalib (重投影约束)
:726/:805/:837 各 prior
:908      ceres::Solve (DENSE_SCHUR + DOGLEG)
:992-1093 ceres::Covariance 恢复 IMU 协方差
```
- **估计量**：各相机帧相对位姿、速度、`_imu` 速度/重力、`bg/ba`；可选外参/内参（`init_dyn_mle_opt_calib`）。

---

## Section 4: 静态 vs 动态 差异

| | 静态 (Static) | 动态 (Dynamic) |
|---|---|---|
| 输入 | 仅 IMU | IMU + 视觉特征库 |
| 运动假设 | 静止启动 | 任意运动 |
| 重力/姿态 | 加速度均值直接对齐（gram_schmidt） | CPI+线性约束(\|g\|=9.81)求 g，再对齐 |
| 速度 | 恒 0 | MLE 估计非零初始速度 |
| 尺度 | IMU 物理尺度，无尺度未知数 | 同（IMU 提供绝对尺度） |
| Bias | 均值估计 | launch 初值 + MLE 精化 |
| 时间偏移 | 无 | 用 `calib_camimu_dt` 处理 cam→IMU（代码层面，未显式优化） |
| 输出 | 单个 `IMU` 状态 | `IMU` + IMU clones + SLAM 特征 + 协方差 |
| 视觉使用 | 无 | 重投影因子约束位姿/特征 |

---

## Section 5: 数学与公式对照

### 5.1 静态：重力对齐

静止时 $\mathbf{a}_{meas} \approx \mathbf{R}\,\mathbf{g}$，故 $\mathbf{R}$ 的 z 轴对齐 $-\mathbf{a}_{avg}$ 方向：
$$
\mathbf{z}_{axis} = \frac{\mathbf{a}_{avg}}{\|\mathbf{a}_{avg}\|},\quad \mathbf{R}_o = \text{gram\_schmidt}(\mathbf{z}_{axis}),\quad \mathbf{q}_{GtoI} = \text{rot\_2\_quat}(\mathbf{R}_o)
$$
只确定 roll/pitch，yaw 自由（VIO 4-DoF 不可观测）。

### 5.2 动态：CPI 预积分测量

逐相机帧对构造 `CpiV1`，输出相对位姿/速度/位置增量 $\{\mathbf{q}_{k\to i}, \boldsymbol{\beta}_{k\to i}, \boldsymbol{\alpha}_{k\to i}\}$ 及偏置雅可比（详见 F3）。

### 5.3 动态：重力约束

线性系统状态 $[\mathbf{f}, \mathbf{v}_{I0}, \mathbf{g}]$，加约束 $\|\mathbf{g}\| = 9.81$（Dong-Si & Mourikis 2012），用配分矩阵最小实根求解。

### 5.4 动态：Ceres 残差

- IMU-CPI 残差（F3 详）：姿态/速度/位置增量 vs 状态预测，含一阶 bias 修正。
- 重投影残差：观测像素 vs 投影（pose+feature+calib）。
- 先验：首帧 yaw+位置+零偏固定。

---

## Section 6: 关键实现细节与易错点

1. **两组 CPI 用途不同**（:243-307）：`I0→Ii` 给线性系统，`Ii→Ii+1` 给 MLE，不可混用。
2. **`imu_avg=true`**（:258/:282）：端点平均，IJRR 论文未做此平均（`CpiBase.h:114-115` 注释）。
3. **重力约束必做**（:395-472）：动态初始化靠 `|g|=9.81` 把尺度/重力解耦出来，跳过会退化。
4. **prior 消不可观测**（:726/:805/:837）：`Factor_GenericPrior` 固定首帧 yaw/位置/零偏，否则 MLE 尺度/偏航漂移。
5. **最小角运动检查**（:169-193）：`init_dyn_min_deg`（默认 45°），运动太小无法初始化动态。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **上游**：F1（`InertialInitializer` 调度）、S17（`try_to_initialize`）。
- **下游**：F3（`CpiV1` 预积分器实现）、Ceres 因子（`ov_init/src/ceres/`）。

### 待详细补充项

- **F2a — `compute_dongsi_coeff` 重力约束求解**：:395-472 的特征值/配分矩阵法细节。
- **F2b — `Factor_ImuCPIv1` 残差构造**：见 F3（Factor 消费 CPI 的部分）。
- **F2c — `Factor_ImageReprojCalib` 重投影因子**：带标定的视觉残差。
- **F2d — `InitializerHelper::gram_schmidt` / `select_imu_readings`**：静态姿态对齐与 IMU 选段辅助。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"）。
