# Propagator 积分路径 精读报告 (S8)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/Propagator.cpp`
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (离散积分 Eq.101/103)、Yang 2020 ACI2 [44] (连续时间解析积分 $\Xi$ 项)、Eckenhoff 2019 IJRR [13] (RK4 均值)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S8（Phase 3 函数级路线，阶段 B 三条积分路径 + F/G 雅可比）

> 上游：S7（`predict_and_compute` `:441-456` 按 `integration_method` 分发到本篇各函数）。本篇是**预测步数学最密处**——三种均值积分（离散/RK4/解析）与两套 F/G 雅可比（离散/解析）。对应 citelist [40]（离散）、[44]（ACI2 解析 $\Xi$ 积分）、[13]（RK4）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
predict_and_compute (S7 :441-456) 分发
   ├─ ANALYTICAL ─┬─ compute_Xi_sum    @ :588 (预计算 Ξ 项)  ┐
   │              ├─ predict_mean_analytic @ :667             ├─ 本篇 S8
   │              └─ compute_F_and_G_analytic @ :683          ┘
   ├─ RK4 ───────┬─ compute_Xi_sum    @ :588 (RK4 也用)      ┐
   │              └─ predict_mean_rk4 @ :507                  ├─ 本篇 S8
   └─ DISCRETE ──┬─ predict_mean_discrete @ :482             ┐
                  └─ compute_F_and_G_discrete @ :830          ┘
```

| 函数 | 文件:行号 | 被谁调用 |
|------|-----------|---------|
| `predict_mean_discrete` | `:482` | S7 `:446` |
| `predict_mean_rk4` | `:507` | S7 `:444` |
| `compute_Xi_sum` | `:588` | S7 `:435` (RK4/ANALYTICAL) |
| `predict_mean_analytic` | `:667` | S7 `:442` |
| `compute_F_and_G_analytic` | `:683` | S7 `:454` |
| `compute_F_and_G_discrete` | `:830` | S7 `:456` |

### 1.2 一句话概括

S8 是 IMU 积分的**底层数学实现**：`predict_mean_discrete` 用零阶+常加速度离散积分（[40] Eq.101/103）；`predict_mean_rk4` 用四阶龙格库塔（[13]）；`predict_mean_analytic` + `compute_Xi_sum` 用 ACI2 解析积分（[44]）预计算 $\Xi_1..\Xi_4$ 闭合项；两套 `compute_F_and_G_*` 给出误差态转移 $\mathbf{F}$ 与噪声耦合 $\mathbf{G}$。

---

## Section 2: 核心数据结构

```cpp
// 积分法由 state->_options.integration_method 决定:
//   DISCRETE  → predict_mean_discrete + compute_F_and_G_discrete
//   RK4       → predict_mean_rk4       + compute_Xi_sum + compute_F_and_G_analytic
//   ANALYTICAL→ predict_mean_analytic  + compute_Xi_sum + compute_F_and_G_analytic
// Xi_sum: 3×18 矩阵, 块存 [R_ktok1 | Ξ_1 | Ξ_2 | Jr_ktok1 | Ξ_3 | Ξ_4] (各 3×3)
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 离散均值积分（默认路径，[40] Eq.101/103）

#### 论文原文（Trawny 2005 [40]）
四元数：$\bar{q}_{k+1} = \big(\cos(\tfrac12\|\boldsymbol{\omega}\|\Delta t)\mathbf{I} + \tfrac{1}{\|\boldsymbol{\omega}\|}\sin(\tfrac12\|\boldsymbol{\omega}\|\Delta t)\boldsymbol{\Omega}(\boldsymbol{\omega})\big)\bar{q}_k$

速度/位置：$\mathbf{v}_{k+1}=\mathbf{v}_k+\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t - \mathbf{g}\Delta t$，$\mathbf{p}_{k+1}=\mathbf{p}_k+\mathbf{v}_k\Delta t+\tfrac12\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t^2-\tfrac12\mathbf{g}\Delta t^2$

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:482-505
void predict_mean_discrete(state, dt, w_hat, a_hat, new_q, new_v, new_p) {
  double w_norm = w_hat.norm();
  Eigen::Matrix4d I_4x4 = I;  Eigen::Matrix3d R_Gtoi = state->_imu->Rot();
  // 四元数: Ω 矩阵指数近似
  Eigen::Matrix4d bigO;
  if (w_norm > 1e-12)
    bigO = cos(0.5*w_norm*dt)*I_4x4 + 1/w_norm*sin(0.5*w_norm*dt)*Omega(w_hat);  // :493
  else bigO = I_4x4 + 0.5*dt*Omega(w_hat);                                       // :495
  new_q = quatnorm(bigO * state->_imu->quat());                                   // :497
  // 速度
  new_v = state->_imu->vel() + R_Gtoi.transpose()*a_hat*dt - _gravity*dt;        // :501
  // 位置
  new_p = state->_imu->pos() + state->_imu->vel()*dt
        + 0.5*R_Gtoi.transpose()*a_hat*dt*dt - 0.5*_gravity*dt*dt;               // :504
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\bar{q}_{k+1}=\big(\cos\tfrac{\|\omega\|\Delta t}{2}\mathbf{I}+\dots\big)\bar{q}_k$ | `bigO=cos(...)*I + sin(...)*Omega(w_hat)` | `:493` |
| $\mathbf{v}_{k+1}=\mathbf{v}+\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t-\mathbf{g}\Delta t$ | `new_v = vel + Rᵀ*a*dt - g*dt` | `:501` |
| $\mathbf{p}_{k+1}=\mathbf{p}+\mathbf{v}\Delta t+\tfrac12\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t^2-\tfrac12\mathbf{g}\Delta t^2$ | `:504` | `:504` |

### 3.2 RK4 均值积分（[13]）

#### 对应代码（结构）

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:507-586
void predict_mean_rk4(state, dt, w_hat1,a_hat1, w_hat2,a_hat2, new_q,new_v,new_p) {
  Eigen::Vector3d w_hat=w_hat1, a_hat=a_hat1;
  Eigen::Vector3d w_alpha=(w_hat2-w_hat1)/dt, a_jerk=(a_hat2-a_hat1)/dt;  // :514-515 线性变化率
  // y0 / k1 / k2 / k3 / k4 各求 q_dot=0.5*Ω(w)*dq, v_dot=Rᵀ*a - g, p_dot=v
  // 经典 RK4 加权:
  dq = quatnorm(dq_0 + 1/6*k1_q + 1/3*k2_q + 1/3*k3_q + 1/6*k4_q);       // :582
  new_q = quat_multiply(dq, q_0);                                          // :583
  new_p = p_0 + 1/6*k1_p + 1/3*k2_p + 1/3*k3_p + 1/6*k4_p;                // :584
  new_v = v_0 + 1/6*k1_v + 1/3*k2_v + 1/3*k3_v + 1/6*k4_v;                // :585
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| RK4 积分加权 $\tfrac16\mathbf{k}_1+\tfrac13\mathbf{k}_2+\tfrac13\mathbf{k}_3+\tfrac16\mathbf{k}_4$ | `:582-585` | `:582-585` |

> **要点**：RK4 用段首尾 `w_hat1/a_hat1` 与 `w_hat2/a_hat2` 线性插值中间点（`:514-515`），比 DISCRETE 的零阶近似更精确，适合高动态。

### 3.3 ACI2 解析积分：`compute_Xi_sum` + `predict_mean_analytic`（[44]）

#### 推导说明

ACI2（Analytical Continuous-time Integration, Yang 2020 [44]）在常角速度+常加速度假设下，把旋转积分闭式求出 $\Xi_1..\Xi_4$ 四项，避免数值积分误差。`compute_Xi_sum` 预计算这四项（及 $\mathbf{R}_{k\to k+1}$、$\mathbf{J}_r$），`predict_mean_analytic` 直接用它们闭合推进均值。

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:588-681
void compute_Xi_sum(state, dt, w_hat, a_hat, Xi_sum) {
  double w_norm=w_hat.norm(); double d_th=w_norm*dt;                       // :592-593
  Eigen::Vector3d k_hat = (w_norm>1e-12)? w_hat/w_norm : Zero();           // :594-597 角速度方向
  // 小角 / 大角两套闭合公式 (small_w = |w| < 0.5°/s)
  bool small_w = (w_norm < 1.0/180*M_PI/2);                                // :620
  if (!small_w) {
    Xi_1 = I*dt + (1-cos_dth)/w_norm*sK + (dt-sin_dth/w_norm)*sK2;         // :624 一阶旋转
    Xi_2 = 1/2*d_t2*I + (d_th-sin_dth)/w_norm2*sK + ...*sK2;              // :627 二阶旋转
    Xi_3 = ... (含 sA, sK, k_hat·a_hat 项) ...;                           // :630-633 一阶旋转×加速度
    Xi_4 = ... (含 d_t3, w_norm3 项) ...;                                 // :636-639 二阶旋转×加速度
  } else {
    Xi_1 = dt*(I + sin_dth*sK + (1-cos_dth)*sK2);                         // :644 小角近似
    Xi_2 = 1/2*dt*Xi_1;  Xi_3 = 1/2*d_t2*(...);  Xi_4 = 1/3*dt*Xi_3;      // :647-654
  }
  // 打包进 Xi_sum (3×18)
  Xi_sum.block(0,0,3,3)  = R_ktok1   = exp_so3(-w_hat*dt);                 // :659 旋转传递
  Xi_sum.block(0,3,3,3)  = Xi_1;                                          // :660
  Xi_sum.block(0,6,3,3)  = Xi_2;                                          // :661
  Xi_sum.block(0,9,3,3)  = Jr_ktok1  = Jr_so3(-w_hat*dt);                 // :662 右雅可比
  Xi_sum.block(0,12,3,3) = Xi_3;                                          // :663
  Xi_sum.block(0,15,3,3) = Xi_4;                                          // :664
}

void predict_mean_analytic(state, dt, w_hat, a_hat, new_q,new_v,new_p, Xi_sum) {
  Eigen::Matrix3d R_Gtok = state->_imu->Rot();
  Eigen::Matrix3d Xi_1 = Xi_sum.block(0,3,3,3);                           // :674
  Eigen::Matrix3d Xi_2 = Xi_sum.block(0,6,3,3);                           // :675
  new_q = quat_multiply(rot_2_quat(Xi_sum.block(0,0,3,3)), state->_imu->quat()); // :678
  new_v = state->_imu->vel() + R_Gtok.transpose()*Xi_1*a_hat - _gravity*dt;        // :679
  new_p = state->_imu->pos() + state->_imu->vel()*dt + R_Gtok.transpose()*Xi_2*a_hat - 0.5*_gravity*dt*dt; // :680
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\Xi_1,\Xi_2$ 旋转积分闭合项 | `Xi_1`/`Xi_2` `:624/:627` (大角) 或 `:644/:647` (小角) | `:624-647` |
| $\Xi_3,\Xi_4$ 旋转×加速度积分 | `Xi_3`/`Xi_4` `:630/:636` (大角) 或 `:650/:654` (小角) | `:630-654` |
| $\mathbf{R}_{k\to k+1}=\exp_{\mathfrak{so}(3)}(-\hat{\boldsymbol{\omega}}\Delta t)$ | `R_ktok1=exp_so3(-w_hat*dt)` | `:615/:659` |
| 右雅可比 $\mathbf{J}_r$ | `Jr_ktok1=Jr_so3(-w_hat*dt)` | `:616/:662` |

> **要点**：ACI2 [44] 的解析项直接对应连续时间 IMU 预积分的闭合解，与 CPI [12][13] 思想同源（见 F 阶段 CpiV1）。它比 RK4 更快且无离散化误差，是 OpenVINS 在精度/算力间的最优解。

### 3.4 F/G 雅可比（下沉：离散 vs 解析）

两套 `compute_F_and_G_*` 分别给出误差态转移矩阵 $\mathbf{F}=\partial\mathbf{x}_{k+1}/\partial\mathbf{x}_k$ 与噪声耦合 $\mathbf{G}=\partial\mathbf{x}_{k+1}/\partial\mathbf{w}$。离散版由 `predict_mean_discrete` 的解析微分得到（`:830` 起），解析版由 `Xi_sum` 的闭合微分得到（`:683` 起）。两者输出同维 $\mathbf{F}\in\mathbb{R}^{N\times N}$、$\mathbf{G}\in\mathbb{R}^{N\times 12}$（$N=$`imu_intrinsic_size+15`）。

> 详细 F/G 元素推导属 S7 调用的底层，本篇锚定其**位置与选型逻辑**：RK4/ANALYTICAL 共用 `compute_F_and_G_analytic`（`:683`），DISCRETE 用 `compute_F_and_G_discrete`（`:830`）。

---

## Section 4: 函数调用链

```
predict_and_compute (S7 :441-456)
  ├─ (ANALYTICAL/RK4) compute_Xi_sum      @ :588  ★ S8
  ├─ predict_mean_discrete                 @ :482  ★ S8 (DISCRETE)
  ├─ predict_mean_rk4                      @ :507  ★ S8 (RK4)
  ├─ predict_mean_analytic                 @ :667  ★ S8 (ANALYTICAL)
  ├─ compute_F_and_G_discrete              @ :830  ★ S8 (DISCRETE)
  └─ compute_F_and_G_analytic              @ :683  ★ S8 (RK4/ANALYTICAL)
```

### 关键函数速查

| 函数名 | 文件:行号 | 积分法 | 功能 |
|--------|-----------|--------|------|
| `predict_mean_discrete` | `:482` | DISCRETE | 零阶+常加速均值 [40] |
| `predict_mean_rk4` | `:507` | RK4 | 四阶龙格库塔 [13] |
| `compute_Xi_sum` | `:588` | RK4/ANALYTICAL | ACI2 Ξ 项预计算 [44] |
| `predict_mean_analytic` | `:667` | ANALYTICAL | ACI2 闭合均值 [44] |
| `compute_F_and_G_analytic` | `:683` | RK4/ANALYTICAL | 解析 F/G |
| `compute_F_and_G_discrete` | `:830` | DISCRETE | 离散 F/G |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `integration_method` | DISCRETE | `StateOptions` | 精度 | 积分法选型 | ★★ 中 |
| `imu_model` | 配置 | `StateOptions` | 标定 | 影响 F/G 内参块 | ★★ 中 |

> 选型权衡：DISCRETE 最快、高频振动误差大；RK4/ACI2 精度高、算力增。ACI2（[44]）为推荐高精度路径。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 默认 DISCRETE 精度有限 | ★ 低 | 高频振动 IMU | 均值漂移 | 切 RK4/ANALYTICAL |
| 2 | `small_w` 阈值硬编码 | ★ 低 | 临界角速度 | 大小角公式切换 | 可参数化 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS (S8) | VINS-Mono | OKVIS |
|---------|---------------|-----------|-------|
| 积分法 | 离散/RK4/**ACI2解析** | 解析预积分 | 欧拉 |
| 解析积分 | $\Xi_1..\Xi_4$ 闭合 (ACI2 [44]) | 流形预积分 | 无 |
| 噪声 Qd | $\mathbf{G}\mathbf{Q}_c\mathbf{G}^\top$ (Trawny) | 解析 | 欧拉 |

### 设计取舍分析

- 选择 **ACI2 解析积分（[44]）作为高精度路径** 是因为：常角速度+常加速度假设下的闭合解 $\Xi$ 项无离散误差，且比 RK4 少函数评估，是精度/算力最优。
- 保留 **DISCRETE 默认** 是因为：多数场景 IMU 频率高（$\Delta t$ 小），离散误差可接受，且实现最简、最稳。
- 提供 **RK4** 是因为：中等动态、不想引入 ACI2 闭合公式复杂度时的折中。
