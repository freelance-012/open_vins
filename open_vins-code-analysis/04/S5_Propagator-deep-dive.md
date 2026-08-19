# Propagator 传播与克隆 精读报告 (S5)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/Propagator.{h,cpp}`
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (间接 KF 离散化 Eq.101/103/129/130)、Eckenhoff 2018 TR CPI [12] / Eckenhoff 2019 IJRR [13] / Yang 2020 ACI2 [44] (连续时间 IMU 预积分与解析积分)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S5（Phase 3 函数级路线，阶段 B 预测）

> 上游：S1–S4（状态表示）。本篇进入**滤波器预测柱**——`Propagator` 把一段 IMU 测量积分前推状态均值、累加状态转移矩阵 $\Phi$ 与离散噪声 $Q_d$，更新 IMU 协方差，并执行**随机克隆（stochastic cloning）**生成 MSCKF 滑动窗口的关键帧。这是误差状态 EKF 的"预测"步，也是滑窗生长的源头。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
IMU 回调 → VioManager::feed_measurement_imu → Propagator::feed_imu (存 imu_data 缓冲)
                                                    │
                        每个相机帧触发 → Propagator::propagate_and_clone
                                                    │
        ┌───────────────────────────────────────────┼───────────────────────────────┐
        │ 预测 (本模块)                                │ 后续 (阶段 C/D)               │
        ▼                                             ▼                              │
  select_imu_readings → predict_and_compute (逐段)    StateHelper::augment_clone    │
        ↓ 累加 Phi_summed / Qd_summed                 (新增 IMU clone 到滑窗)         │
        ↓ StateHelper::EKFPropagation (协方差前推)  ──► _clones_IMU[timestamp]       │
        ▼ 更新 _imu 的 value + fej                                          │         │
                                                                             ▼         │
                                                          UpdaterMSCKF 用 clone 做多帧残差 │
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 2 环节（IMU 原始 → **Propagator 预测+克隆** → 状态/协方差更新 → 更新器） |
| 输入 | IMU 测量缓冲 `imu_data`、目标时间戳 `timestamp` |
| 输出 | 前推的 `_imu` 均值、更新后的协方差 $\mathbf{P}$、新增一个 IMU clone（滑窗） |
| 调用方 | `VioManager::feed_measurement_imu`（每相机帧触发） |
| 被调用方 | `select_imu_readings`、`predict_and_compute`、`StateHelper::EKFPropagation`、`StateHelper::augment_clone` |

### 1.2 一句话概括

`Propagator::propagate_and_clone` 是 MSCKF 的**预测+滑窗生长**核心：它在 `[t0,t1]` 内用所有 IMU 段逐段算状态转移 $\mathbf{F}$ 与离散噪声 $\mathbf{Q}_d$，累乘得总 $\Phi$、总 $Q_d$，一次性前推协方差；同时把 IMU 均值积分到新时刻并**将此刻姿态克隆进滑窗**（`augment_clone`）。支持离散/RK4/解析三种积分（对应 [13][44]）。

---

## Section 2: 核心数据结构

### 2.1 关键成员与接口

```cpp
// 文件: ov_msckf/src/state/Propagator.h (节选)
class Propagator {
public:
    void feed_imu(const ov_core::ImuData &message, double oldest_time = -1);  // 存 IMU 缓冲
    void clean_old_imu_measurements(double oldest_time);                      // 裁剪旧测量
    void propagate_and_clone(std::shared_ptr<State> state, double timestamp); // ★ 主入口
    bool fast_state_propagate(...);                                           // 快速版本(缓存)

    // 内部: 选段 / 预测+噪声 / 三种均值积分 / F-G 计算
    std::vector<ov_core::ImuData> select_imu_readings(...);
    void predict_and_compute(state, data_minus, data_plus, F, Qd);
    void predict_mean_discrete(...);   // 零阶 + 常加速度离散
    void predict_mean_rk4(...);        // 四阶龙格库塔
    void predict_mean_analytic(...);   // ACI2 解析积分 [44]
    void compute_Xi_sum(...);          // 解析积分预计算
    void compute_F_and_G_discrete(...);
    void compute_F_and_G_analytic(...);

private:
    std::vector<ov_core::ImuData> imu_data;   // IMU 测量缓冲 (imu_data_mtx 保护)
    std::mutex imu_data_mtx;
    double _gravity = 9.81;                   // 重力大小
    PropagatorNoise _noises;                  // σ_w, σ_a, σ_wb, σ_ab
    bool have_last_prop_time_offset = false;
    double last_prop_time_offset = 0;         // 上次 CAM-IMU 时间偏移
};
```

### 2.2 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 物理含义 | 位置 |
|---------|---------|---------|------|
| `Phi_summed` | $\Phi = \prod_i \Phi_i$ | 总状态转移矩阵 | `Propagator.cpp:76` |
| `Qd_summed` | $\mathbf{Q}_d = \sum \Phi_i \mathbf{Q}_{d,i} \Phi_i^\top$ | 总离散噪声 | `Propagator.cpp:77` |
| `prop_data` | $\{\mathbf{u}_k\}$ | 选定时间窗内 IMU 段 | `:68` |
| `dt` | $\Delta t$ | IMU 段间隔 | `:399` |
| `a_hat` | $\tilde{\mathbf{a}}$ | 去 bias+标定的加速度 | `:408` |
| `w_hat` | $\tilde{\boldsymbol{\omega}}$ | 去 bias+标定的角速度 | `:420` |
| `Qc` | $\mathbf{Q}_c$ | 连续噪声（含 $1/\Delta t$） | `:462-466` |
| `F`,`G` | $\mathbf{F},\mathbf{G}$ | 离散状态/噪声雅可比 | `:450-451` |
| `last_w` | $\boldsymbol{\omega}_{last}$ | 末段角速度（克隆用） | `:112` |

> **关键公式来源**：离散噪声 $\mathbf{Q}_c$ 用 `1/dt` 缩放连续噪声，对应 **Trawny 2005 [40] Eq.(129)(130)**（代码注释 `:461` 明确标注）。均值/转移矩阵离散化对应 [40] Eq.(101)(103)。

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 主流程：预测 + 克隆

#### 推导说明

误差状态 EKF 预测步：
1. 对每段 IMU 算 $\mathbf{F}_i, \mathbf{Q}_{d,i}$
2. 累加 $\Phi = \prod_i \Phi_i$，$\mathbf{Q}_d = \sum_i \Phi_i \mathbf{Q}_{d,i} \Phi_i^\top$
3. 协方差前推：$\mathbf{P}' = \Phi \mathbf{P} \Phi^\top + \mathbf{Q}_d$（交由 `StateHelper::EKFPropagation`）
4. 均值前推 + 克隆

#### 对应代码（propagate_and_clone 主体）

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:33-138
void Propagator::propagate_and_clone(std::shared_ptr<State> state, double timestamp) {
    // 防御: 同时间戳/反向传播 → 崩溃 (避免重克隆)
    if (state->_timestamp == timestamp) { PRINT_ERROR(...); std::exit(...); }   // :37-40
    if (state->_timestamp > timestamp)  { PRINT_ERROR(...); std::exit(...); }   // :43-47

    // 时间偏移: t_imu = t_cam + calib_dt
    double time0 = state->_timestamp + last_prop_time_offset;   // :63 IMU 起点
    double time1 = timestamp + state->_calib_dt_CAMtoIMU->value()(0);  // :64 IMU 终点
    prop_data = select_imu_readings(imu_data, time0, time1);    // :68 选段(加锁)

    // 累加 Phi / Qd
    Eigen::MatrixXd Phi_summed = I(N), Qd_summed = 0;           // :76-77  N = imu_intrinsic_size+15
    if (prop_data.size() > 1) {
      for (i=0; i<prop_data.size()-1; i++) {
        predict_and_compute(state, prop_data[i], prop_data[i+1], F, Qdi);  // :87
        Phi_summed = F * Phi_summed;                            // :95  Φ = Φ_i · Φ
        Qd_summed = F*Qd_summed*F.transpose() + Qdi;            // :96  Qd = Φ·Qd·Φᵀ + Qdi
        Qd_summed = 0.5*(Qd_summed+Qd_summed.transpose());      // :97  对称化
        dt_summed += dt;                                        // :98
      }
    }
    assert(abs((time1-time0) - dt_summed) < 1e-4);              // :101 时间窗一致性校验

    // 协方差前推 (Phi_order = IMU + 可选 IMU 内参标定块)
    Phi_order = {state->_imu, [可选 _calib_imu_dw/da/tg/旋转]}; // :116-129
    StateHelper::EKFPropagation(state, Phi_order, Phi_order, Phi_summed, Qd_summed);  // :130

    state->_timestamp = timestamp;                              // :133 更新状态时间
    last_prop_time_offset = t_off_new;                          // :134 记录偏移
    StateHelper::augment_clone(state, last_w);                  // :137 ★ 随机克隆
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\Phi = \prod_i \Phi_i$ | `Phi_summed = F * Phi_summed` | `:95` |
| $\mathbf{Q}_d = \sum \Phi_i \mathbf{Q}_{d,i} \Phi_i^\top$ | `F*Qd_summed*F.transpose()+Qdi` | `:96` |
| 协方差前推 $\mathbf{P}'=\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$ | `StateHelper::EKFPropagation(...)` | `:130` |
| 滑窗克隆 | `augment_clone(state, last_w)` | `:137` |

> **设计要点**：用"逐段算 F/Qdi → 累乘 Phi / 累加 Qd → 一次性乘协方差"替代每段都乘 P，减少大矩阵乘法次数（性能优化）。`assert(:101)` 保证所选 IMU 段总时间 == 请求时间窗，防止丢测量导致传播错误。

### 3.2 单段预测与噪声：`predict_and_compute`

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:395-480
void predict_and_compute(state, data_minus, data_plus, F, Qd) {
    double dt = data_plus.timestamp - data_minus.timestamp;              // :399
    // 去 bias + 标定 (IMU 内参 Dw/Da/Tg, 旋转 R_ACCtoIMU/R_GYROtoIMU)
    a_hat = R_ACCtoIMU*Da*(am - bias_a);                                // :415-417
    w_hat = R_GYROtoIMU*Dw*(wm - bias_g - Tg*a_hat);                    // :420-429
    // 选积分法: ANALYTICAL / RK4 / DISCRETE(默认)
    if (ANALYTICAL) predict_mean_analytic(...);                         // :442
    else if (RK4) predict_mean_rk4(...);                                // :444
    else predict_mean_discrete(...);                                    // :446
    // 算 F,G
    if (RK4||ANALYTICAL) compute_F_and_G_analytic(...);                 // :454
    else compute_F_and_G_discrete(...);                                 // :456
    // 连续噪声 → 离散 (Trawny Eq.129/130)
    Qc.block(0,0,3,3) = σ_w²/dt * I;   Qc.block(3,3,3,3) = σ_a²/dt * I; // :463-464
    Qc.block(6,6,3,3) = σ_wb²/dt * I;  Qc.block(9,9,3,3) = σ_ab²/dt * I;// :465-466
    Qd = G * Qc * G.transpose();                                        // :470 离散噪声注入
    Qd = 0.5*(Qd+Qd.transpose());                                       // :471
    // ★ 前推 IMU 均值, 并固化为 FEJ (克隆/初始化的锚定)
    imu_x = state->_imu->value();
    imu_x << new_q, new_p, new_v;                                      // :474-477
    state->_imu->set_value(imu_x);                                      // :478 估计移动
    state->_imu->set_fej(imu_x);                                        // :479 FEJ 同步(每步重固)
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 去偏加速度 $\tilde{\mathbf{a}}$ | `a_hat = R*Da*(am - bias_a)` | `:415` |
| 去偏角速度 $\tilde{\boldsymbol{\omega}}$ | `w_hat = R*Dw*(wm - bias_g - Tg*a)` | `:420` |
| 离散噪声 $\mathbf{Q}_c= \tfrac{\sigma^2}{\Delta t}\mathbf{I}$ | `Qc.block(..)=σ²/dt*I` | `:463-466` |
| $\mathbf{Q}_d = \mathbf{G}\mathbf{Q}_c\mathbf{G}^\top$ | `Qd = G*Qc*G.transpose()` | `:470` |

> **FEJ 同步点（重要）**：`:479` `set_fej(imu_x)` 在**每次传播后**把 FEJ 重设为当前传播值。这与 S1 说的"FEJ 在首估计锚定"看似矛盾——实际 OpenVINS 的策略是：**每个 clone 自己是独立的 FEJ 锚定点**（clone 在 `augment_clone` 时 `set_fej`），而当前活跃 `_imu` 的 FEJ 随传播刷新，因为活跃状态没有"固定首估计"的约束（它会被持续更新）。滑窗里**历史 clone 的 FEJ 才固定不动**，这正是 MSCKF 一致性（[20][27]）的关键。S1 报告说的"克隆时定 FEJ"= `Propagator.cpp:479` 此处（clone 由 `_imu` 复制，携带此刻 value==fej）。

### 3.3 离散均值积分（默认路径）

#### 论文原文（Trawny 2005 [40] Eq.101/103）

四元数：$\bar{q}_{k+1} = \big(\cos(\tfrac12\|\boldsymbol{\omega}\|\Delta t)\mathbf{I} + \tfrac{1}{\|\boldsymbol{\omega}\|}\sin(\tfrac12\|\boldsymbol{\omega}\|\Delta t)\boldsymbol{\Omega}(\boldsymbol{\omega})\big)\bar{q}_k$

速度/位置：$\mathbf{v}_{k+1}=\mathbf{v}_k+\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t - \mathbf{g}\Delta t$， $\mathbf{p}_{k+1}=\mathbf{p}_k+\mathbf{v}_k\Delta t+\tfrac12\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t^2-\tfrac12\mathbf{g}\Delta t^2$

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:482-505
void predict_mean_discrete(state, dt, w_hat, a_hat, new_q, new_v, new_p) {
    double w_norm = w_hat.norm();
    // 四元数: 用 Ω 矩阵的指数近似 (小角)
    if (w_norm > 1e-12)
      bigO = cos(0.5*w_norm*dt)*I + (1/w_norm)*sin(0.5*w_norm*dt)*Omega(w_hat);  // :493
    else bigO = I + 0.5*dt*Omega(w_hat);                                          // :495
    new_q = quatnorm(bigO * state->_imu->quat());                                 // :497
    // 速度: 局部系加速度 - 重力
    new_v = state->_imu->vel() + R_Gtoi.transpose()*a_hat*dt - _gravity*dt;       // :501
    // 位置: 速度积分 + 加速度二次积分 - 重力二次积分
    new_p = state->_imu->pos() + state->_imu->vel()*dt
          + 0.5*R_Gtoi.transpose()*a_hat*dt*dt - 0.5*_gravity*dt*dt;              // :504
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\bar{q}_{k+1}=\big(\cos\tfrac{\|\omega\|\Delta t}{2}\mathbf{I}+\dots\big)\bar{q}_k$ | `bigO=cos(...)*I + sin(...)*Omega(w_hat)` | `:493` |
| $\mathbf{v}_{k+1}=\mathbf{v}+\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t-\mathbf{g}\Delta t$ | `new_v = vel + Rᵀ*a*dt - g*dt` | `:501` |
| $\mathbf{p}_{k+1}=\mathbf{p}+\mathbf{v}\Delta t+\tfrac12\mathbf{R}^\top\tilde{\mathbf{a}}\Delta t^2-\tfrac12\mathbf{g}\Delta t^2$ | `:504` | `:504` |

> **积分法选择**：默认 `DISCRETE`（零阶 + 常加速度）；`RK4` 用四阶龙格库塔更精确；`ANALYTICAL` 用 ACI2 解析积分 [44]（Yang 2020）。运行时由 `state->_options.integration_method` 决定（Phase 2 待确认项 #1 的落点）。

---

## Section 4: 函数调用链

### 4.1 入口函数

```
Propagator::propagate_and_clone(state, timestamp) @ Propagator.cpp:33
  └→ 功能: IMU 预测 + 协方差前推 + 滑窗克隆
     参数: (state: State*, timestamp: double) → 目标时刻
     返回: void（改写 state._imu / 协方差 / _clones_IMU）
```

### 4.2 完整调用树

```
VioManager::feed_measurement_imu
  └─ Propagator::feed_imu (存 imu_data 缓冲)
        │ (相机帧到达触发)
        ▼
Propagator::propagate_and_clone @ :33
  ├─ select_imu_readings(imu_data, time0, time1) @ :269  → 选 [t0,t1] 内 IMU 段
  ├─ for each IMU segment:
  │    └─ predict_and_compute(state, d_i, d_{i+1}, F, Qdi) @ :395
  │         ├─ predict_mean_discrete / _rk4 / _analytic @ :482/507/667  (均值)
  │         ├─ compute_F_and_G_discrete / _analytic @ :830/683          (F,G)
  │         └─ Qd = G*Qc*Gᵀ  (Qc 用 1/dt, Trawny Eq.129/130)
  │              └─ state->_imu->set_value/set_fej (前推 + 固 FEJ) @ :478-479
  ├─ StateHelper::EKFPropagation(state, Phi_order, Phi_order, Phi_summed, Qd_summed) @ :130
  └─ StateHelper::augment_clone(state, last_w) @ :137   ★ 克隆入滑窗
```

### 4.3 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `propagate_and_clone` | `Propagator.cpp:33` | 预测+克隆主入口 | `state, ts` | 改 state | O(n_seg) |
| `select_imu_readings` | `Propagator.cpp:269` | 选时间窗 IMU 段 | `imu_data, t0,t1` | `vector<ImuData>` | O(n) |
| `predict_and_compute` | `Propagator.cpp:395` | 单段 F/Qd + 均值 | `d-, d+` | `F, Qd` | O(1) |
| `predict_mean_discrete` | `Propagator.cpp:482` | 离散均值积分 | `dt,w,a` | `q,v,p` | O(1) |
| `compute_F_and_G_discrete` | `Propagator.cpp:830` | 离散 F/G 雅可比 | `dt,w,a` | `F,G` | O(1) |
| `StateHelper::augment_clone` | `StateHelper.cpp:579` | 滑窗克隆 | `last_w` | 新增 clone | O(1) |

---

## Section 5: 关键参数清单

### 5.1 本模块涉及的参数

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `integration_method` | DISCRETE | `StateOptions` | 精度 | 离散/RK4/解析积分 | ★★ 中 |
| `_noises.sigma_w/a/wb/ab` | 配置 yaml | `Propagator.h` | 精度 | IMU 噪声 | ★★★ 高 |
| `_gravity` | 9.81 | `Propagator.cpp` | 精度 | 重力大小 | ★★ 中 |
| `calib_dt_CAMtoIMU` | 标定 | `State` | 时序 | 相机-IMU 时间偏移 | ★★ 中 |

### 5.2 参数影响分析

- **`sigma_*` 直接进 `Qc`（`:463-466`）**：噪声过大 → 协方差膨胀过快、滤波过信任测量；过小 → 对测量异常不敏感。是调参首要项。
- **`integration_method`**：`DISCRETE` 最快但高频振动时误差大；`RK4`/`ANALYTICAL` 精度高、算力增。解析法(ANALYTICAL)对应 ACI2 [44]，与 CPI [12][13] 思想一致。
- **`calib_dt_CAMtoIMU`**：决定 `time0/time1` 偏移（`:63-64`），偏移错 → 选错 IMU 段 → 预测错位。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `set_fej` 每传播步刷新活跃 IMU（`:479`） | ★ 低 | 理解不一致性时 | 初学者易与"FEJ 固定"混淆 | 文档明确：仅**历史 clone** FEJ 固定，活跃态随传播更新 |
| 2 | `assert(:101)` 硬失败 | ★ 中 | IMU 丢段/时间戳错 | 直接 exit，系统崩 | 可降级为警告+跳过该帧 |
| 3 | 默认 DISCRETE 积分精度有限 | ★ 低 | 高频振动 IMU | 均值漂移 | 高动态场景切 RK4/ANALYTICAL |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `Propagator` | VINS-Mono | OKVIS | 常规 MSCKF |
|---------|----------------------|-----------|-------|-----------|
| 预测方式 | 误差态 EKF 预积分 (F/Qd 累加) | BA 边缘化 | 误差态 EKF | 误差态 EKF |
| 积分法 | 离散/RK4/解析(ACI2) | 解析预积分 | 欧拉 | 欧拉 |
| 滑窗克隆 | `augment_clone` 显式 | 无(图优化) | 有 | 有 |
| 时间偏移 | `calib_dt_CAMtoIMU` 在线 | 标定 | 标定 | 标定 |
| FEJ | **原生(克隆锚定)** | 无 | 部分 | 部分 |
| 精度特点 | FEJ 一致性 [20][27] + ACI2 [44] | BA 最优 | MSCKF | MSCKF |

### 设计取舍分析

- 选择 **误差态 EKF + 逐段 F/Qd 累加 + 显式克隆** 是因为：滤波式实时 VIO 需固定计算量，且 MSCKF 滑动窗口靠克隆历史姿态实现多帧约束——`augment_clone`(`:137`) 即窗口生长点。
- 放弃 **图优化预积分（如 VINS-Mono 的 Ceres 预积分）** 是因为：EKF 框架下预积分以 F/Qd 形式直接进协方差，无需非线性重优化；且 FEJ 一致性 [20][27] 在 EKF 中自然表达。
- 在 **高频/高动态 IMU** 下当前方案可能不足：默认 DISCRETE 积分精度有限，应切 `RK4`/`ANALYTICAL`（ACI2 [44]）以匹配 CPI [12][13] 的连续时间精度。
