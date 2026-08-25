# 动态初始化 DynamicInitializer 算法深度解析 (F1b)

> **生成日期**：2026-08-25
> **前置文档**：`F1_InertialInitializer-deep-dive.md`（调度层）、`F2_Static-Dynamic-Initializer-deep-dive.md`（算法层概览）、`F1a_StaticInitializer-algorithm-deep-dive.md`（静态初始化）
> **核心代码**：`ov_init/src/dynamic/DynamicInitializer.cpp`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Dong-Si & Mourikis "Estimator initialization in vision-aided inertial navigation with unknown camera-IMU calibration" IROS 2012; Eckenhoff et al. "Continuous Preintegration Theory for Graph-based Visual-Inertial Navigation" IJRR 2019

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
InertialInitializer::initialize()  (F1)
        │  策略判定: init_dyn_use && !is_still
        ▼
DynamicInitializer::initialize(timestamp, covariance, order, _imu, _clones_IMU, _features_SLAM)  ★ 本篇
   ├─ 1. 数据准备（窗口裁剪 + 特征复制）
   ├─ 2. 特征有效性验证（测量数 + 角运动检查）
   ├─ 3. CPI 预积分构造（I0→Ii 线性系统 + Ii→Ii+1 MLE）
   ├─ 4. 线性系统求解（特征位置 + 初速度 + 重力方向）
   ├─ 5. Gram-Schmidt 重力对齐
   ├─ 6. Ceres MLE 优化（位姿 + 速度 + bias + 特征）
   ├─ 7. 协方差恢复（ceres::Covariance）
   └─ 8. 状态输出（IMU + clones + SLAM 特征）
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | 动态初始化算法实现 |
| **输入** | `imu_data`（IMU 缓冲）、`_db`（特征库） |
| **输出** | `timestamp`、`covariance`、`order`、`_imu`（IMU 状态）、`_clones_IMU`（位姿）、`_features_SLAM`（特征） |
| **调用方** | `InertialInitializer`（F1） |
| **被调用方** | `CpiV1`（F3）、`InitializerHelper::select_imu_readings`（F2d）、Ceres 因子 |

### 1.2 一句话概括

`DynamicInitializer` 假设载体**任意运动**，利用 CPI 连续预积分构造线性系统求解特征位置/初速度/重力方向，再用 Ceres 图优化精化所有参数（位姿+速度+bias+特征），最后恢复协方差并输出完整初始化状态。

---

## 2. 数学模型

### 2.1 问题定义

给定初始化窗口 $[t_0, t_N]$ 内的：
- 相机观测：$\{(u_{j,k}, v_{j,k})\}$ —— 特征 $j$ 在时刻 $t_k$ 的像素坐标
- IMU 测量：$\{(\boldsymbol{\omega}_m^{(k)}, \mathbf{a}_m^{(k)})\}$

求解：
- IMU 状态：$\mathbf{x}_I = [{}^{I}_G\bar{q}^\top,\; {}^G\mathbf{p}_I^\top,\; {}^G\mathbf{v}_I^\top,\; \mathbf{b}_\omega^\top,\; \mathbf{b}_a^\top]^\top$
- 特征位置：$\{\mathbf{p}_{F_j}^G\}$
- 重力方向：$\mathbf{g}^G$（约束 $\|\mathbf{g}^G\| = 9.81$）

### 2.2 CPI 预积分测量

对每对连续相机时间戳 $(t_{k}, t_{k+1})$，构造两组 CPI：

**组 1：$I_0 \to I_k$（相对首帧）**

用于线性系统，输出：
- $\Delta\mathbf{R}_{0\to k}$：相对旋转
- $\boldsymbol{\alpha}_{0\to k}$：位置增量
- $\boldsymbol{\beta}_{0\to k}$：速度增量

**组 2：$I_k \to I_{k+1}$（逐帧）**

用于 Ceres MLE，输出相邻帧间的预积分测量。

### 2.3 线性系统构造

利用 CPI 测量，在 $I_0$ 系中写出 IMU 运动方程：

$$\mathbf{p}_{I_k}^{I_0} = \mathbf{v}_{I_0}^{I_0}\Delta t_{0k} - \frac{1}{2}\mathbf{g}^{I_0}\Delta t_{0k}^2 + \boldsymbol{\alpha}_{0\to k} \tag{F1b-1}$$

$$\mathbf{v}_{I_k}^{I_0} = \mathbf{v}_{I_0}^{I_0} - \mathbf{g}^{I_0}\Delta t_{0k} + \boldsymbol{\beta}_{0\to k} \tag{F1b-2}$$

特征重投影约束：

$$\begin{bmatrix} 1 & 0 & -u \\ 0 & 1 & -v \end{bmatrix} \mathbf{p}_{F_j}^{C_k} = \mathbf{0} \tag{F1b-3}$$

其中：

$$\mathbf{p}_{F_j}^{C_k} = \mathbf{R}_{I\to C}\left(\mathbf{R}_{I_0\to I_k}^\top(\mathbf{p}_{F_j}^{I_0} - \mathbf{p}_{I_k}^{I_0})\right) + \mathbf{p}_I^C \tag{F1b-4}$$

将式 (F1b-1) 代入 (F1b-4)，再代入 (F1b-3)，得到关于 $\mathbf{p}_{F_j}^{I_0}$、$\mathbf{v}_{I_0}^{I_0}$、$\mathbf{g}^{I_0}$ 的线性方程：

$$\mathbf{Y}\mathbf{p}_{F_j}^{I_0} - \Delta t_{0k}\mathbf{Y}\mathbf{v}_{I_0}^{I_0} + \frac{1}{2}\Delta t_{0k}^2\mathbf{Y}\mathbf{g}^{I_0} = \mathbf{Y}\boldsymbol{\alpha}_{0\to k} - \mathbf{H}_{\text{proj}}\mathbf{p}_I^C \tag{F1b-5}$$

其中 $\mathbf{Y} = \mathbf{H}_{\text{proj}}\mathbf{R}_{I\to C}\mathbf{R}_{I_0\to I_k}^\top$。

堆叠所有特征和所有时刻的观测，得到线性系统：

$$\mathbf{A}\mathbf{x} = \mathbf{b}, \quad \mathbf{x} = [\mathbf{p}_{F_1}^{I_0},\, \ldots,\, \mathbf{p}_{F_M}^{I_0},\, \mathbf{v}_{I_0}^{I_0},\, \mathbf{g}^{I_0}]^\top \tag{F1b-6}$$

### 2.4 重力约束

式 (F1b-6) 是超定线性系统，但 $\mathbf{g}^{I_0}$ 需满足物理约束：

$$\|\mathbf{g}^{I_0}\| = g = 9.81\,\text{m/s}^2 \tag{F1b-7}$$

**求解方法**（Dong-Si & Mourikis 2012）：

1. 对 $\mathbf{A}$ 做 QR 分解，分离出关于 $\mathbf{g}$ 的约束
2. 构造拉格朗日函数，导出关于 $\lambda$ 的 6 次多项式
3. 用伴随矩阵特征值分解求根
4. 选择使 $\|\mathbf{g}\|$ 最接近 9.81 的实根

详见 `F2a_compute_dongsi_coeff-deep-dive.md`。

### 2.5 Ceres MLE 优化

线性系统给出初值后，用 Ceres 图优化精化：

**优化变量**：
- 位姿：$\{(\mathbf{R}_{G\to I_k}, \mathbf{p}_{I_k}^G)\}_{k=0}^N$
- 速度：$\{\mathbf{v}_{I_k}^G\}_{k=0}^N$
- 零偏：$\mathbf{b}_\omega, \mathbf{b}_a$
- 特征：$\{\mathbf{p}_{F_j}^G\}$

**残差因子**：
- IMU CPI 因子（连接相邻帧）
- 重投影因子（特征观测约束）
- 先验因子（固定首帧 yaw+位置，消除 4-DoF 不可观性）

**求解**：DENSE_SCHUR + DOGLEG

---

## 3. 代码实现：逐段精读

### 3.1 数据准备

```cpp
// DynamicInitializer.cpp:48-94
// 1. 取最新/最老时间戳
double newest_cam_time = -1;
for (auto const &feat : _db->get_internal_data()) {
    for (auto const &camtimepair : feat.second->timestamps) {
        for (auto const &time : camtimepair.second) {
            newest_cam_time = std::max(newest_cam_time, time);
        }
    }
}
double oldest_time = newest_cam_time - params.init_window_time;

// 2. 裁剪窗口
_db->cleanup_measurements(oldest_time);

// 3. 深拷贝特征（线程安全）
std::unordered_map<size_t, std::shared_ptr<Feature>> features;
for (const auto &feat : _db->get_internal_data()) {
    auto feat_new = std::make_shared<Feature>();
    feat_new->featid = feat.second->featid;
    feat_new->uvs = feat.second->uvs;
    feat_new->uvs_norm = feat.second->uvs_norm;
    feat_new->timestamps = feat.second->timestamps;
    features.insert({feat.first, feat_new});
}
```

### 3.2 特征有效性验证

```cpp
// DynamicInitializer.cpp:100-164
const int min_num_meas_to_optimize = (int)params.init_window_time;
const int min_valid_features = 8;

// 验证每个特征的测量数
for (auto const &feat : features) {
    std::vector<double> times;
    for (auto const &camtime : feat.second->timestamps) {
        for (double time : camtime.second) {
            double time_dt = INFINITY;
            for (auto const &tmp : map_camera_times) {
                time_dt = std::min(time_dt, std::abs(time - tmp.first));
            }
            for (auto const &tmp : times) {
                time_dt = std::min(time_dt, std::abs(time - tmp));
            }
            // 采样：时间间隔 ≥ pose_dt_avg 或已是已知时刻
            if (time_dt >= pose_dt_avg || time_dt == 0.0) {
                times.push_back(time);
            }
        }
    }
    map_features_num_meas[feat.first] = (int)times.size();
    if (map_features_num_meas[feat.first] < min_num_meas_to_optimize)
        continue;  // 测量数不足，跳过
    count_valid_features++;
}

// 检查有效特征数
if (count_valid_features < min_valid_features) {
    return false;  // 特征不足
}
```

### 3.3 角运动检查

```cpp
// DynamicInitializer.cpp:172-194
double theta_inI_norm = 0.0;
for (size_t k = 0; k < readings.size() - 1; k++) {
    auto imu0 = readings.at(k);
    auto imu1 = readings.at(k + 1);
    double dt = imu1.timestamp - imu0.timestamp;
    Eigen::Vector3d wm = 0.5 * (imu0.wm + imu1.wm) - gyroscope_bias;
    theta_inI_norm += (-wm * dt).norm();  // 角增量
}

if (180.0 / M_PI * theta_inI_norm < params.init_dyn_min_deg) {
    return false;  // 角运动不足（默认 45°）
}
```

**理论意义**：足够的角运动是可观性的必要条件。没有旋转，无法区分重力和加速度。

### 3.4 CPI 预积分构造

```cpp
// DynamicInitializer.cpp:240-307
std::map<double, std::shared_ptr<ov_core::CpiV1>> map_camera_cpi_I0toIi, map_camera_cpi_IitoIi1;

for (auto const &timepair : map_camera_times) {
    double current_time = timepair.first;
    if (current_time == oldest_camera_time) {
        map_camera_cpi_I0toIi.insert({current_time, nullptr});
        map_camera_cpi_IitoIi1.insert({current_time, nullptr});
        continue;
    }
    
    // 组 1: I0 → Ii（线性系统）
    auto cpiI0toIi1 = std::make_shared<ov_core::CpiV1>(...);
    cpiI0toIi1->setLinearizationPoints(gyroscope_bias, accelerometer_bias);
    for (size_t k = 0; k < readings.size() - 1; k++) {
        cpiI0toIi1->feed_IMU(...);
    }
    map_camera_cpi_I0toIi.insert({current_time, cpiI0toIi1});
    
    // 组 2: Ii-1 → Ii（MLE）
    auto cpiIitoIi1 = std::make_shared<ov_core::CpiV1>(...);
    // ... 类似构造
    map_camera_cpi_IitoIi1.insert({current_time, cpiIitoIi1});
}
```

**两组 CPI 的区别**：
- `I0toIi`：相对首帧，用于线性系统（式 F1b-1/F1b-2）
- `IitoIi1`：相邻帧，用于 Ceres MLE（IMU 因子连接相邻位姿）

### 3.5 线性系统求解

```cpp
// DynamicInitializer.cpp:309-520
// 构造 A x = b
Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_measurements, system_size);
Eigen::VectorXd b = Eigen::VectorXd::Zero(num_measurements);

for (auto const &feat : features) {
    for (auto const &camtime : feat.second->timestamps) {
        for (size_t i = 0; i < camtime.second.size(); i++) {
            // 投影矩阵
            Eigen::MatrixXd H_proj = Eigen::MatrixXd::Zero(2, 3);
            H_proj << 1, 0, -uv_norm(0),
                      0, 1, -uv_norm(1);
            
            // Y = H_proj * R_ItoC * R_I0toI_k
            Eigen::MatrixXd Y = H_proj * R_ItoC * R_I0toIk;
            
            // 右侧常数项
            Eigen::MatrixXd b_i = Y * alpha_I0toIk - H_proj * p_IinC;
            
            // 填入 A 和 b
            H_i.block(0, 3*feat_idx, 2, 3) = Y;              // 特征位置
            H_i.block(0, 3*N_f + 0, 2, 3) = -DT * Y;         // 速度
            H_i.block(0, 3*N_f + 3, 2, 3) = 0.5 * DT*DT * Y; // 重力
        }
    }
}

// 重力约束求解
Eigen::Matrix<double, 7, 1> coeff = InitializerHelper::compute_dongsi_coeff(D, d, g);
// 伴随矩阵特征值分解 → 选最优实根
// 回代：g = (D - λ*I)^{-1} * d
// 回代：x1 = -A1A1_inv * A1^T * A2 * g + A1A1_inv * A1^T * b
```

### 3.6 Gram-Schmidt 重力对齐

```cpp
// DynamicInitializer.cpp:551-567
// 重力方向 → 全局坐标系
InitializerHelper::gram_schmidt(gravity_inI0, R_GtoI0);

// 转换到全局系
for (auto const &timepair : map_camera_times) {
    ori_GtoIi[timepair.first] = quat_multiply(ori_I0toIi.at(timepair.first), q_GtoI0);
    pos_IiinG[timepair.first] = R_GtoI0.transpose() * pos_IiinI0.at(timepair.first);
    vel_IiinG[timepair.first] = R_GtoI0.transpose() * vel_IiinI0.at(timepair.first);
}
for (auto const &feat : features_inI0) {
    features_inG[feat.first] = R_GtoI0.transpose() * feat.second;
}
```

### 3.7 Ceres MLE 优化

```cpp
// DynamicInitializer.cpp:572-908
ceres::Problem problem;

// 添加变量
for (auto const &timepair : map_camera_times) {
    // 位姿 + 速度 + bias
    problem.AddParameterBlock(var_ori, 4, ceres_jplquat);
    problem.AddParameterBlock(var_pos, 3);
    problem.AddParameterBlock(var_vel, 3);
    problem.AddParameterBlock(var_bias_g, 3);
    problem.AddParameterBlock(var_bias_a, 3);
}
for (auto const &feat : features) {
    problem.AddParameterBlock(var_feat, 3);
}

// 添加因子
// IMU CPI 因子
problem.AddResidualBlock(new Factor_ImuCPIv1(...), nullptr, factor_params);
// 重投影因子
problem.AddResidualBlock(new Factor_ImageReprojCalib(...), loss_function, factor_params);
// 先验因子（固定 yaw+位置）
problem.AddResidualBlock(new Factor_GenericPrior(...), nullptr, factor_params);

// 求解
ceres::Solve(options, &problem, &summary);
```

### 3.8 协方差恢复

```cpp
// DynamicInitializer.cpp:978-1073
ceres::Covariance problem_cov(options_cov);
problem_cov.Compute(covariance_blocks, &problem);

// 提取 IMU 协方差
CHECK(problem_cov.GetCovarianceBlockInTangentSpace(var_ori, var_ori, covtmp.data()));
covariance.block(0, 0, 3, 3) = covtmp.eval();
// ... 其他块

// 膨胀
covariance.block(0, 0, 3, 3) *= params.init_dyn_inflation_orientation;
covariance.block(6, 6, 3, 3) *= params.init_dyn_inflation_velocity;
covariance.block(9, 9, 3, 3) *= params.init_dyn_inflation_bias_gyro;
covariance.block(12, 12, 3, 3) *= params.init_dyn_inflation_bias_accel;

// 对称化
covariance = 0.5 * (covariance + covariance.transpose());
```

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{A}$ | 线性系统矩阵 | `A` | `:311` |
| $\mathbf{b}$ | 线性系统右侧 | `b` | `:312` |
| $\mathbf{Y}$ | 投影+旋转组合 | `Y` | `:363` |
| $\Delta\mathbf{R}_{0\to k}$ | CPI 相对旋转 | `R_I0toIk` | `:351` |
| $\boldsymbol{\alpha}_{0\to k}$ | CPI 位置增量 | `alpha_I0toIk` | `:352` |
| $\mathbf{g}^{I_0}$ | 重力方向 | `state_grav` | `:462` |
| $\mathbf{v}_{I_0}^{I_0}$ | 初速度 | `v_I0inI0` | `:469` |
| $\mathbf{P}_0$ | 初始协方差 | `covariance` | `:1022` |

---

## 5. 关键设计决策

### 5.1 为什么需要两组 CPI

| CPI 组 | 用途 | 数学意义 |
|--------|------|---------|
| `I0toIi` | 线性系统 | 式 (F1b-1)(F1b-2)，相对首帧的运动方程 |
| `IitoIi1` | Ceres MLE | IMU 因子连接相邻帧，约束位姿变化 |

### 5.2 为什么需要角运动检查

没有旋转时，重力方向与加速度方向无法区分（可观性缺失）。45° 是经验阈值，确保足够的旋转激励。

### 5.3 为什么协方差需要膨胀

线性化误差 + 有限观测窗口 → Hessian 逆过于乐观。膨胀系数：
- 姿态 10×：线性化误差显著
- 速度 10×：有限观测约束
- 零偏 100×：初值不准，观测少

---

## 6. 完整流程图

```
输入: imu_data, _db (特征库)
        │
        ├─ 1. 数据准备
        │     ├─ newest/oldest 时间戳
        │     ├─ 窗口裁剪
        │     └─ 特征深拷贝（线程安全）
        │
        ├─ 2. 特征验证
        │     ├─ 测量数 ≥ min_num_meas
        │     ├─ 有效特征 ≥ min_valid_features
        │     └─ 角运动 ≥ init_dyn_min_deg (45°)
        │
        ├─ 3. CPI 构造
        │     ├─ I0toIi (线性系统)
        │     ─ IitoIi1 (MLE)
        │
        ├─ 4. 线性系统
        │     ├─ 构造 A x = b
        │     ├─ compute_dongsi_coeff (重力约束)
        │     ├─ 伴随矩阵特征值分解
        │     └─ 回代：g, v, 特征位置
        │
        ├─ 5. 重力对齐
        │     ├─ gram_schmidt(gravity_inI0)
        │     └─ 转换到全局系
        │
        ├─ 6. Ceres MLE
        │     ├─ 添加变量（位姿+速度+bias+特征）
        │     ├─ 添加因子（IMU+重投影+先验）
        │     └─ 求解（DENSE_SCHUR + DOGLEG）
        │
        ├─ 7. 协方差恢复
        │     ├─ ceres::Covariance
        │     ├─ 膨胀（10×/100×）
        │     └─ 对称化
        │
        └─ 8. 状态输出
              ├─ _imu (IMU 状态)
              ├─ _clones_IMU (位姿)
              └─ _features_SLAM (特征)
```

---

## 7. 与其他文档的关联

- **上游**：`F1_InertialInitializer-deep-dive.md`（调度层）
- **下游**：`F3_CpiV1-continuous-preintegration-deep-dive.md`（CPI 实现）、`F2a_compute_dongsi_coeff-deep-dive.md`（重力约束求解）
- **相关**：`S12_StateHelper-initialize-deep-dive.md`（SLAM 特征初始化）

### 待详细补充项

- **F2a**：`compute_dongsi_coeff` 重力约束求解（伴随矩阵 + 特征值分解）
- **F2b**：`Factor_ImuCPIv1` 残差构造（CPI 测量 vs 状态预测）
- **F2c**：`Factor_ImageReprojCalib` 重投影因子（带标定的视觉残差）
- **F2d**：`InitializerHelper::gram_schmidt` / `select_imu_readings`（辅助函数）

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"）。
