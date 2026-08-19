# Propagator 单段预测 精读报告 (S7)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/Propagator.cpp`
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (离散噪声 Eq.129/130)、Eckenhoff 2018 TR CPI [12] / Eckenhoff 2019 IJRR [13] / Yang 2020 ACI2 [44] (连续时间 IMU 预积分思想)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S7（Phase 3 函数级路线，阶段 B 单段 F/Qd 计算）

> 上游：S5（编排入口 `:87` 调用）、S6（选段提供 `prop_data`）。本篇聚焦 **`Propagator::predict_and_compute` (:395)** —— 对**单段 IMU** 去 bias、转内参坐标系、按积分法（S8）算均值与状态转移 $\mathbf{F}$、构造离散噪声 $\mathbf{Q}_d=\mathbf{G}\mathbf{Q}_c\mathbf{G}^\top$，并把 IMU 均值前推 + **重设 FEJ**。这是预测步的数学核心，S8 是其下沉的三种积分路径。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
propagate_and_clone (S5 :87) 逐段调用
        │  (prop_data[i], prop_data[i+1])
        ▼
predict_and_compute(state, data_minus, data_plus, F, Qd)  @ :395  ★ 本篇
        ├─ 去 bias + 内参坐标变换 (Dw/Da/Tg, R_ACCtoIMU/R_GYROtoIMU)
        ├─ compute_Xi_sum (仅 RK4/ANALYTICAL)  → S8 :588
        ├─ predict_mean_* (按 integration_method)  → S8 :482/:507/:667
        ├─ compute_F_and_G_*  → S8 :830/:683
        ├─ Qd = G * Qc * Gᵀ  (Qc 用 1/dt, Trawny Eq.129/130)
        └─ state->_imu->set_value / set_fej  (均值前推 + FEJ 同步)
        ▼ 返回 F, Qd 给 S5 累加
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 预测第 2 步（单段积分） |
| 输入 | `data_minus`/`data_plus`（相邻两段 IMU）、`state` |
| 输出 | `F`（状态转移）、`Qd`（离散噪声），并改写 `_imu` 均值/FEJ |
| 调用方 | `propagate_and_clone` (S5 :87) |
| 被调用方 | `predict_mean_*`、`compute_F_and_G_*`、`compute_Xi_sum`（均属 S8）、`StateHelper` 无 |

### 1.2 一句话概括

`predict_and_compute` 对单段 IMU 先做**偏差校正**（减 bias、去重力灵敏度 `Tg·a`、转 IMU 内参坐标系），再据 `integration_method` 分发到三种均值积分（S8）与对应 F/G 雅可比，最后用连续噪声 $\mathbf{Q}_c=\sigma^2/\Delta t$（Trawny Eq.129/130）构造离散噪声 $\mathbf{Q}_d=\mathbf{G}\mathbf{Q}_c\mathbf{G}^\top$，并把前推后的 IMU 均值**同步写回 FEJ**（活跃态 FEJ 随传播刷新，见 S9/S10 讨论）。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:395-480
void Propagator::predict_and_compute(state,
    const ImuData &data_minus, const ImuData &data_plus,
    Eigen::MatrixXd &F,   // 输出: 状态转移 (N×N, N=imu_intrinsic_size+15)
    Eigen::MatrixXd &Qd); // 输出: 离散噪声 (N×N)
// 内部使用的噪声/内参:
//   _noises.sigma_w/a/wb/ab  → Qc 对角
//   _calib_imu_dw/da/tg, _calib_imu_ACCtoIMU/GYROtoIMU → 内参坐标变换
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 偏差校正与内参坐标变换

#### 推导说明

原始 IMU 测量 $\tilde{\boldsymbol{\omega}}_m, \tilde{\mathbf{a}}_m$ 含 bias 与内参误差。校正：
$$\hat{\mathbf{a}} = \mathbf{R}_{ACC\to IMU}\mathbf{D}_a(\tilde{\mathbf{a}}_m - \mathbf{b}_a),\quad \hat{\boldsymbol{\omega}} = \mathbf{R}_{GYRO\to IMU}\mathbf{D}_w(\tilde{\boldsymbol{\omega}}_m - \mathbf{b}_g - \mathbf{T}_g\hat{\mathbf{a}})$$

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:398-429
double dt = data_plus.timestamp - data_minus.timestamp;                     // :399
// IMU 内参 (静态标定)
Eigen::Matrix3d Dw = State::Dm(state->_options.imu_model, _calib_imu_dw->value()); // :403
Eigen::Matrix3d Da = State::Dm(state->_options.imu_model, _calib_imu_da->value()); // :404
Eigen::Matrix3d Tg = State::Tg(state->_calib_imu_tg->value());              // :405
// 去 bias (加速度取段首尾平均)
Eigen::Vector3d a_hat1 = data_minus.am - state->_imu->bias_a();             // :408
Eigen::Vector3d a_hat2 = data_plus.am  - state->_imu->bias_a();             // :409
Eigen::Vector3d a_hat_avg = .5*(a_hat1+a_hat2);                             // :410
// 转 ACCtoIMU 坐标系
Eigen::Matrix3d R_ACCtoIMU = _calib_imu_ACCtoIMU->Rot();                    // :414
a_hat1 = R_ACCtoIMU*Da*a_hat1;  a_hat2 = R_ACCtoIMU*Da*a_hat2;              // :415-416
a_hat_avg = R_ACCtoIMU*Da*a_hat_avg;                                        // :417
// 去 bias + 重力灵敏度 Tg·a (陀螺)
Eigen::Vector3d w_hat1 = data_minus.wm - _imu->bias_g() - Tg*a_hat1;        // :420
Eigen::Vector3d w_hat2 = data_plus.wm  - _imu->bias_g() - Tg*a_hat2;        // :421
Eigen::Vector3d w_hat_avg = .5*(w_hat1+w_hat2);                             // :422
// 转 GYROtoIMU 坐标系
Eigen::Matrix3d R_GYROtoIMU = _calib_imu_GYROtoIMU->Rot();                  // :426
w_hat1 = R_GYROtoIMU*Dw*w_hat1;  w_hat2 = R_GYROtoIMU*Dw*w_hat2;            // :427-428
w_hat_avg = R_GYROtoIMU*Dw*w_hat_avg;                                       // :429
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\hat{\mathbf{a}}=\mathbf{R}_{ACC\to IMU}\mathbf{D}_a(\tilde{\mathbf{a}}-\mathbf{b}_a)$ | `a_hat = R_ACCtoIMU*Da*(am - bias_a)` | `:415-417` |
| $\hat{\boldsymbol{\omega}}=\mathbf{R}_{GYRO\to IMU}\mathbf{D}_w(\tilde{\boldsymbol{\omega}}-\mathbf{b}_g-\mathbf{T}_g\hat{\mathbf{a}})$ | `w_hat = R_GYROtoIMU*Dw*(wm - bias_g - Tg*a)` | `:420-429` |

### 3.2 均值积分与 F/G 分发（下沉 S8）

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:431-457
// 解析积分预计算 Xi_sum (仅 RK4/ANALYTICAL 需要)
Eigen::Matrix<double,3,18> Xi_sum = Zero(3,18);                            // :432
if (RK4 || ANALYTICAL) compute_Xi_sum(state, dt, w_hat_avg, a_hat_avg, Xi_sum); // :433-436 → S8 :588

// 均值积分: 按 integration_method 选路
Eigen::Vector4d new_q; Eigen::Vector3d new_v, new_p;
if (ANALYTICAL)      predict_mean_analytic(state, dt, w_hat_avg, a_hat_avg, new_q, new_v, new_p, Xi_sum); // :441-442 → S8 :667
else if (RK4)        predict_mean_rk4(state, dt, w_hat1,a_hat1, w_hat2,a_hat2, new_q,new_v,new_p);        // :443-444 → S8 :507
else                 predict_mean_discrete(state, dt, w_hat_avg,a_hat_avg, new_q,new_v,new_p);            // :445-446 → S8 :482

// F/G 分配 (RK4/ANALYTICAL 用解析版, 否则离散版)
F = Zero(N, N);  Eigen::MatrixXd G = Zero(N, 12);                          // :450-451
if (RK4 || ANALYTICAL) compute_F_and_G_analytic(state, dt, w_hat_avg,a_hat_avg, w_uncorrected,a_uncorrected, new_q,new_v,new_p, Xi_sum, F, G); // :452-454 → S8 :683
else                  compute_F_and_G_discrete(state, dt, w_hat_avg,a_hat_avg, w_uncorrected,a_uncorrected, new_q,new_v,new_p, F, G);          // :455-456 → S8 :830
```

#### 对照注释

| 公式项 | 对应代码 | 行号 | 下沉 |
|--------|---------|------|------|
| 均值积分 | `predict_mean_*` | `:441-446` | **S8** |
| $\mathbf{F},\mathbf{G}$ 雅可比 | `compute_F_and_G_*` | `:452-456` | **S8** |

### 3.3 离散噪声构造（Trawny Eq.129/130）

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:459-479
// 连续噪声 Qc: 用 1/dt 缩放 (Trawny Eq.129/130)
Eigen::Matrix<double,12,12> Qc = Zero();
Qc.block(0,0,3,3) = pow(_noises.sigma_w,2)/dt * I;   // :463 角速度噪声
Qc.block(3,3,3,3) = pow(_noises.sigma_a,2)/dt * I;   // :464 加速度噪声
Qc.block(6,6,3,3) = pow(_noises.sigma_wb,2)/dt * I;  // :465 bias_g 随机游走
Qc.block(9,9,3,3) = pow(_noises.sigma_ab,2)/dt * I;  // :466 bias_a 随机游走
// 离散噪声注入: Qd = G * Qc * Gᵀ
Qd = Zero(N,N);
Qd = G * Qc * G.transpose();                          // :470
Qd = 0.5*(Qd+Qd.transpose());                         // :471 对称化

// 前推 IMU 均值 + 固化为 FEJ
Eigen::Matrix<double,16,1> imu_x = state->_imu->value();
imu_x.block(0,0,4,1) = new_q;  // :475 四元数
imu_x.block(4,0,3,1) = new_p;  // :476 位置
imu_x.block(7,0,3,1) = new_v;  // :477 速度
state->_imu->set_value(imu_x);                                // :478 估计移动
state->_imu->set_fej(imu_x);                                 // :479 FEJ 同步(每步重固)
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{Q}_c=\tfrac{\sigma^2}{\Delta t}\mathbf{I}$ | `Qc.block(..)=σ²/dt*I` | `:463-466` |
| $\mathbf{Q}_d=\mathbf{G}\mathbf{Q}_c\mathbf{G}^\top$ | `Qd = G*Qc*G.transpose()` | `:470` |
| 均值前推 + FEJ | `set_value` / `set_fej` | `:478-479` |

> **FEJ 同步点（重要）**：`:479` `set_fej(imu_x)` 在**每次传播后**把 FEJ 重设为当前传播值。这与 S1 说的"FEJ 在首估计锚定"看似矛盾——实际 OpenVINS 的策略是：**每个 clone 自己是独立的 FEJ 锚定点**（clone 在 S10 `augment_clone` 时 `set_fej`），而当前活跃 `_imu` 的 FEJ 随传播刷新，因为活跃状态没有"固定首估计"的约束（它会被持续更新）。滑窗里**历史 clone 的 FEJ 才固定不动**，这正是 MSCKF 一致性（[20][27]）的关键。

---

## Section 4: 函数调用链

```
predict_and_compute @ :395  ★ 本篇
  ├─ compute_Xi_sum(state, dt, w,a, Xi_sum)           @ :435 → S8 :588
  ├─ predict_mean_analytic / _rk4 / _discrete         @ :441-446 → S8 :667/:507/:482
  ├─ compute_F_and_G_analytic / _discrete             @ :452-456 → S8 :683/:830
  └─ state->_imu->set_value / set_fej                 @ :478-479 → S1 类型系统
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 | 下沉 |
|--------|-----------|------|------|
| `predict_and_compute` | `Propagator.cpp:395` | 单段 F/Qd + 均值 | **S7 (本篇)** |
| `compute_Xi_sum` | `Propagator.cpp:588` | 解析积分预计算 | **S8** |
| `predict_mean_discrete` | `Propagator.cpp:482` | 离散均值 | **S8** |
| `predict_mean_rk4` | `Propagator.cpp:507` | RK4 均值 | **S8** |
| `predict_mean_analytic` | `Propagator.cpp:667` | ACI2 解析均值 | **S8** |
| `compute_F_and_G_discrete` | `Propagator.cpp:830` | 离散 F/G | **S8** |
| `compute_F_and_G_analytic` | `Propagator.cpp:683` | 解析 F/G | **S8** |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `_noises.sigma_w/a/wb/ab` | yaml | `Propagator.h` | 精度 | IMU 噪声 | ★★★ 高 |
| `integration_method` | DISCRETE | `StateOptions` | 精度 | 积分法选型 | ★★ 中 |
| `imu_model` | 配置 | `StateOptions` | 标定 | KALIBR/标量内参 | ★★ 中 |

> `sigma_*` 直接进 `Qc`（`:463-466`）——噪声过大→协方差膨胀过快；过小→对异常不敏感。积分法选型见 S8。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `set_fej` 每步刷新活跃 IMU | ★ 低 | 理解一致性时 | 易与"FEJ 固定"混淆 | 明确：仅历史 clone FEJ 固定（S10） |
| 2 | 默认 DISCRETE 精度有限 | ★ 低 | 高频振动 | 均值漂移 | 切 RK4/ANALYTICAL (S8) |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `predict_and_compute` | VINS-Mono | OKVIS |
|---------|-------------------------------|-----------|-------|
| 偏差校正 | `Tg·a` 重力灵敏度 + 内参矩阵 | 简化 | 简化 |
| 离散噪声 | $\mathbf{Q}_c=\sigma^2/\Delta t$ (Trawny) | 解析预积分 | 欧拉 |
| 积分法 | 离散/RK4/解析(ACI2) | 解析 | 欧拉 |
| FEJ | 活跃态随传播刷新 | 无 | 部分 |

### 设计取舍分析

- 选择 **显式内参坐标变换（Dw/Da/Tg + R_ACCtoIMU/R_GYROtoIMU）** 是因为：支持 KALIBR 全模型标定（S5 `:124`），在线可标定 IMU 内参。
- 选择 **三种积分法分发** 是因为：DISCRETE 最快、RK4/ANALYTICAL（ACI2 [44]）精度高，按算力/动态性权衡（见 S8 对比）。
