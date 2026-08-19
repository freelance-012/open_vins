# CpiV1 连续时间预积分 精读报告 (F3)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/cpi/CpiV1.cpp` / `CpiV1.h`（继承 `CpiBase`）
> **对应论文**: 连续时间 IMU 预积分（Continuous Pre-Integration, Dong-Si & Mourikis 2012）；与在线 Propagator 离散积分对照 [12][13][44]
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 F3（Phase 3 函数级路线，阶段 F 初始化，CpiV1 实现细节）

> 上游：F2（`DynamicInitializer` 内部构造并消费 `CpiV1`）。本篇聚焦 **`CpiV1`**（连续时间预积分器）——它把 IMU 积分做成一个**与状态解耦的相对测量**：以固定锚帧 $k$ 为参考，解析积分出相对位姿/速度/位置增量及偏置雅可比，供 Ceres 因子（`Factor_ImuCPIv1`）做批量 MLE。它与在线 `Propagator`（S5–S8）**思想同源、实现不同**（连续解析 vs 离散数值）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
DynamicInitializer::initialize  (F2, :243-307)
        │  make_shared<CpiV1>(sigma_w, sigma_wb, sigma_a, sigma_ab, true)
        │  setLinearizationPoints(bg, ba)
        │  feed_IMU(t0,t1,wm0,am0,wm1,am1) 逐段
        ▼
CpiV1 (连续时间预积分器)  ★ 本篇
   输出 (public 成员, CpiBase.h:118-144):
     DT, alpha_tau (位置增量), beta_tau (速度增量), q_k2tau/R_k2tau
     J_q,J_a,J_b,H_a,H_b (偏置雅可比), P_meas (15×15 协方差), Q_c (12×12 噪声)
        │
        ▼
Factor_ImuCPIv1 (ov_init/src/ceres/Factor_ImuCPIv1.cpp:28-75)
   残差 = 预积分测量 vs 当前状态 (含一阶 bias 重线性化)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 6 环节（初始化实现细节，被 F2 调用） |
| 输入 | IMU 段 `(t0,t1,wm0,am0,wm1,am1)`、`setLinearizationPoints(bg,ba)` |
| 输出 | 相对锚帧的预积分测量（均值+雅可比+协方差） |
| 调用方 | `DynamicInitializer`（F2） |
| 被调用方 | `quat_ops.h`（skew_x / rot_2_quat）、Ceres 因子 |

---

## Section 2: 类接口

### 2.1 继承与构造

```cpp
class CpiV1 : public CpiBase { ... };        // CpiV1.h:52
CpiV1(sigma_w, sigma_wb, sigma_a, sigma_ab, imu_avg_=false);  // CpiV1.h:63-64 → 转发 CpiBase
CpiBase(...)  // CpiBase.h:60-74: 填 Q_c (12×12 连续噪声), 单位向量 e_1/2/3 及其 skew
```

### 2.2 方法

| 方法 | 位置 | 说明 |
|------|------|------|
| `setLinearizationPoints(b_w_lin, b_a_lin, q_k_lin, grav)` | `CpiBase.h:88-95` | 设偏置线性化点（V1 不传 q/grav） |
| `feed_IMU(t0,t1,wm0,am0,wm1,am1)` | `CpiV1.cpp:33-337` | **唯一核心方法**，~300 行 |
| 结果读取 | `CpiBase.h:118-144` | 直接读 public 成员（无 get 函数） |

**用法三步**：构造 → `setLinearizationPoints` → 循环 `feed_IMU` → 读 public 变量。

---

## Section 3: feed_IMU 内部（连续解析积分 + RK4 协方差）

### 3.1 均值积分（解析闭式）(`:33-127`)

```cpp
:38    DT += delta_t
:46-47 w_hat = wm_0 - b_w_lin;  a_hat = am_0 - b_a_lin  (去偏)
:50-55 imu_avg 时两端点平均
:70    small_w = (mag_w < 0.008726646)  (0.5°/s, 双分支避免 1/mag_w^n 爆炸)
:88-89 R_tau2tau1 (闭式 Rodrigues, Eq.35); R_k2tau1 = R_tau2tau1 * R_k2tau (Eq.70)
:102-112 多项式系数 f_1..f_4 (连续时间核心): ∫R(τ)dτ, ∫∫R(τ)dτ² 的解析积分写成 w_x, w_x² 标量系数
          small_w 分支 → Taylor 展开 (dt³/3, dt⁴/8, dt²/2, dt³/6)
:116-122 alpha_arg=(dt²/2)I+f_1·w_x+f_2·w_x²; beta_arg=dt·I+f_3·w_x+f_4·w_x²
:126-127 alpha_tau += beta_tau*dt + H_al*a_hat;  beta_tau += H_be*a_hat;
:335-336 最后才提交 R_k2tau = R_k2tau1; q_k2tau = rot_2_quat(...)  (协方差用旧姿态)
```

**关键**：与在线 Propagator 的离散 RK4/ACI2 不同，CPI 把积分写成**闭式多项式系数** `f_1..f_4` 一次性算出（连续时间解析积分），无需数值步进。

### 3.2 偏置雅可比（解析）(`:135-238`)

```cpp
:135-137 右雅可比 J_r_tau1 (Eq.164)
:141     J_q = R_tau2tau1*J_q + J_r_tau1*delta_t  (Eq.81)
:145-147 H_a -= H_al; H_a += dt*H_b; H_b -= H_be  (Eq.49)
:151-153 d_R_bw_1/2/3 (Eq.84)
:172-220 df_i_dbw_j 12 个标量导数 (Eq.55-63)
:224-238 J_a/J_b 按列更新 (Eq.53/54/59)
```
- 给出预积分测量对 **bias 线性化点** 的一阶雅可比 `J_*`/`H_*`，使 Ceres 因子能**一阶重线性化**而不必重新积分。

### 3.3 协方差（15×15 RK4）(`:245-331`)

```cpp
:245-248 中点旋转 R_mid
:255-259 F_k1.. 状态序 [θ(0), b_w(3), v/beta(6), b_a(9), p/alpha(12)]; F.block(12,6)=I (α̇=β)
:269/290/301/324 Riccati: P_dot = F*P + P*F^T + G*Q_c*G^T
:330    P_meas += (dt/6)*(k1+2k2+2k3+k4)  (RK4 积分)
:331    强制对称化
```
- 协方差传播是**连续 Riccati** 的 RK4 数值解（与在线离散 `ΦPΦᵀ+Q_d` 对照）。

---

## Section 4: 输出（预积分测量）

`CpiBase.h:118-144`：

| 成员 | 行 | 含义 |
|------|------|------|
| `DT` | :119 | 总积分时长 |
| `alpha_tau` | :120 | 位置增量 |
| `beta_tau` | :121 | 速度增量 |
| `q_k2tau` / `R_k2tau` | :122-123 | 相对锚帧姿态 |
| `J_q` | :126 | θ wrt b_w |
| `J_a` / `J_b` | :127-128 | α/β wrt b_w |
| `H_a` / `H_b` | :129-130 | α/β wrt b_a |
| `P_meas` | :144 | 15×15 协方差 |
| `Q_c` | :141 | 12×12 连续噪声 |

---

## Section 5: Factor_ImuCPIv1 消费方式

`ov_init/src/ceres/Factor_ImuCPIv1.cpp:28-75`：

```cpp
:52   covariance.llt().solve(I) 求信息矩阵
:60-61 sqrtI_save = L^T (LLT)
:64-74 15 残差 + 10 参数块 [q,bg,v,ba,p]×2
:108-109 dbw = b_w1 - b_w_lin; dba = b_a1 - b_a_lin  (一阶偏置重线性化)
:130-136 残差:
   res(0:3)   = 2*q_res_plus.head(3)                      // 姿态 (含 q_b=[0.5*J_q*dbw;1])
   res(3:6)   = b_w2 - b_w1
   res(6:9)   = R_1*(v_2 - v_1 + g*dt) - J_b*dbw - H_b*dba - beta
   res(9:12)  = b_a2 - b_a1
   res(12:15) = R_1*(p_2 - p_1 - v_1*dt + 0.5*g*dt²) - J_a*dbw - H_a*dba - alpha
   res = sqrtI * res
```
- 残差即"预积分测量 vs 当前两帧状态预测"，偏置用一阶雅可比修正（避免重积分）。

---

## Section 6: 与在线 Propagator 对照

| 维度 | CpiV1（F3） | Propagator（S5–S8，在线） |
|------|------|------|
| 参考帧 | 固定锚帧 k（`R_k2tau` 累积），bias 冻结在 `b_*_lin` | 全局帧，绝对状态 |
| 均值积分 | **解析闭式**：`f_1..f_4` 多项式系数一次算出 | 离散 RK4 / ACI2 数值步进 |
| 协方差 | RK4 解连续 Riccati `Ṗ=FP+PFᵀ+GQGᵀ` | 离散 `Φ P Φᵀ + Q_d` |
| 偏置变化 | 一阶雅可比修正（`J_*`/`H_*`），不重积分 | 状态含 bias，直接更新 |
| 消费者 | Ceres 因子残差（批量 MLE） | EKF 状态预测 |
| 数值保护 | `small_w` 双分支 Taylor | 步长控制 |

**共同点**：都用 `skew_x`/`rot_2_quat`（`quat_ops.h`），同样的 `[θ,b_w,v,b_a,p]` 15 维误差排序，同样的 `F`/`G` 结构（`-w_x`、`-Rᵀa_x`、`F(12,6)=I`）。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **上游**：F2（`DynamicInitializer` 构造并调用 `feed_IMU`）。
- **下游**：`Factor_ImuCPIv1`（Ceres 残差）、`quat_ops.h`（S2 相关数学）。

### 待详细补充项

- **F3a — `f_1..f_4` 解析系数推导**：:102-112 与 Dong-Si & Mourikis 2012 Eq.35/70 的对应。
- **F3b — 偏置雅可比推导**：:135-238（Eq.49/53/54/59/81/84/164）。
- **F3c — RK4 协方差 F/G 矩阵构造**：:245-331 的 `F_k1..4`/`G_k1..4` 块填充。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"）。
