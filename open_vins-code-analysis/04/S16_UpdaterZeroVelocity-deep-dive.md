# UpdaterZeroVelocity 零速更新 精读报告 (S16)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/update/UpdaterZeroVelocity.cpp`
> **对应论文**: ZUPT (Zero Velocity Update) 静止约束；Maybeck 1982 Vol.1 Eq.(7-21a)-(7-21c)（白化/残差加权）；Trawny 2005 [40] (JPL 姿态)；Li & Mourikis 2013 [27]（时间偏移）
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S16（Phase 3 函数级路线，阶段 D 测量更新，ZUPT 零速约束）

> 上游：S5（`select_imu_readings` 选段）、S7（`predict_and_compute` 的 bias 校正框架）、S9（`EKFPropagation` / `EKFUpdate` 协方差更新）。本篇聚焦 **`UpdaterZeroVelocity::try_update` (:65)** —— 当系统检测到"载体静止"时，用 IMU 本体的**零速/零加速度约束**做一步 EKF 更新，把漂移的 velocity / bias 拉回。它不依赖任何视觉特征，是纯 IMU 的**自标定式**更新，常在手持静置、电梯停顿、地面机器人停驻时触发。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
VioManager (静止检测判定后)
        │
        ▼
UpdaterZeroVelocity::try_update(state, timestamp)  @ :65  ★ 本篇
   ├─ 1. 选段: Propagator::select_imu_readings(imu_data, time0, time1)  :97  → S5
   ├─ 2. 构造残差 + 雅可比 (零速约束, 白化)                          :136-180
   ├─ 3. measurement_compress_inplace 压缩过定系统                  :183  → S13c
   ├─ 4. R = _zupt_noise_multiplier * I                            :190
   ├─ 5. 偏置前向传播 Q_bias (dt*sigma_b^2)                        :192-196  → S9
   ├─ 6. χ² 门限 (预计算 chi_squared_table, 0.95)                  :198-216
   ├─ 7. 视差检查 disparity_passed (可选, 覆盖拒绝)                :218-237
   ├─ 8. 接受判定: 过 χ² 且 |v| < _zupt_max_velocity 且未过视差     :241
   └─ 9. 执行: EKFPropagation(bias) + EKFUpdate(state, Hx_order, H, res, R)  :263-277  → S9
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 4 环节（视觉测量更新之外，可选的纯 IMU 零速更新） |
| 输入 | `state`（含 `_imu`、`_clones_IMU`、`_calib_dt_CAMtoIMU`、`_calib_imu_*`）、`timestamp`（当前图像时间戳） |
| 输出 | 接受时改写 `state`（velocity 置零、bias 演化、协方差收缩）；返回 `true`/`false` 表示本次是否触发了 ZUPT |
| 调用方 | `VioManager`（静止检测逻辑判定应启用 ZUPT 时调用） |
| 被调用方 | `Propagator::select_imu_readings`（S5）、`UpdaterHelper::measurement_compress_inplace`（S13c）、`StateHelper::EKFPropagation` / `EKFUpdate`（S9） |

> 注：ZUPT 含残差更新数学（P0 级），但在 pipeline 中属"条件触发的可选更新"，整体模块归 P1（策略模块）。精读时并入 P0-2/3 的残差框架一起看（03 清单附录 A 备注）。

---

## Section 2: 函数骨架与调用关系

`try_update` 走一条"构造零速残差 → 白化 → 偏置演化 → χ² 门限 → 视差覆盖 → EKF 更新"的线：

```
try_update(state, timestamp)
│
├─ 0. 前置检查: imu_data 空 / state 已在该时刻 → false      :67-77
├─ 1. 时间偏移处理 (last_prop_time_offset / calib_dt)       :79-100
├─ 2. 选段: select_imu_readings(time0, time1)               :97  → S5
│       imu_recent.size()<2 → false                          :103
│
├─ 3. 组 Hx_order = [q, bg, ba] (+v 若 integrated)          :117-123
├─ 4. 大循环逐 IMU 段构造 H(6*N×9) 与 res(6*N)              :144-180
│       res = [-w_omega*w_hat; -w_accel*(a_hat - R*g)]       :161-166
│       H   = d(res)/d[q,bg,ba]                              :170-178
│
├─ 5. measurement_compress_inplace(H, res)                  :183  → S13c
├─ 6. R = _zupt_noise_multiplier * I                        :190
├─ 7. Q_bias (dt_summed * sigma_b^2)                        :192-196  → S9
├─ 8. χ² = res^T S^{-1} res, S = H P_marg H^T + R           :198-206
│       阈值 chi_squared_table[rows] (0.95)                  :208-216
│
├─ 9. 视差检查 disparity_passed (override_with_disparity)    :218-237
├─ 10. 接受判定: !disparity_passed && (χ²过阈 || |v|>max) → false  :241
│
└─ 11. 执行更新:
        model_time_varying_bias → EKFPropagation(bias, Q_bias) :263-273  → S9
        EKFUpdate(state, Hx_order, H, res, R)                  :276  → S9
        state->_timestamp = timestamp                          :277
```

---

## Section 3: 逐段精读

### 3.1 前置与时间偏移 (`:67-107`)

```cpp
if (imu_data.empty()) { ... return false; }                 // :67
if (state->_timestamp == timestamp) { ... return false; }   // :74 已在该时刻
```
- 两个快速退出：没有 IMU 数据、或状态已推进到目标时刻（无需更新）。

```cpp
double t_off_new = state->_calib_dt_CAMtoIMU->value()(0);   // :89 当前 IMU-CAM 偏移
double time0 = state->_timestamp + last_prop_time_offset;   // :93 上次用过的偏移
double time1 = timestamp + t_off_new;                       // :94 本次偏移
std::vector<ov_core::ImuData> imu_recent =
    Propagator::select_imu_readings(imu_data, time0, time1); // :97 → S5
```
- `time0` 用 `last_prop_time_offset`（上一次的值），`time1` 用 `t_off_new`（当前值），以兼容时间偏移的缓慢变化（[27]）。
- 选段逻辑复用 S5 的 `select_imu_readings`，保证与传播用的 IMU 区间一致。
- `imu_recent.size() < 2` 则无足够数据 → `false`（`:103`）。

### 3.2 残差与雅可比构造（零速约束）(`:116-180`)

#### 3.2.1 状态顺序 (`:117-123`)

```cpp
std::vector<std::shared_ptr<Type>> Hx_order;
Hx_order.push_back(state->_imu->q());   // :118 姿态
Hx_order.push_back(state->_imu->bg());  // :119 陀螺零偏
Hx_order.push_back(state->_imu->ba());  // :120 加速度零偏
// (integrated_accel_constraint 时再加 v)  —— 当前默认关闭
```
- 默认更新量：`[q_GtoI, bg, ba]`（9 维）；若开 `integrated_accel_constraint` 再加 `v`（12 维，代码标注 `untested`）。

#### 3.2.2 逐段残差 (`:144-166`)

```cpp
double dt = imu_recent.at(i+1).timestamp - imu_recent.at(i).timestamp;
Eigen::Vector3d a_hat = state->_calib_imu_ACCtoIMU->Rot() * Da * (imu_recent.at(i).am - state->_imu->bias_a());  // :148
Eigen::Vector3d w_hat = state->_calib_imu_GYROtoIMU->Rot() * Dw * (imu_recent.at(i).wm - state->_imu->bias_g() - Tg*a_hat);  // :149
```
- `a_hat` / `w_hat`：去外参旋转 + 标度 `Da/Dw` + 去零偏后的"真值估计"角速度与比力。注意 `w_hat` 还减了 `Tg*a_hat`（陀螺 g 灵敏度，[27] 时间/加速度耦合）。

```cpp
double w_omega = std::sqrt(dt) / _noises.sigma_w;   // :156 角速度白化权重
double w_accel = std::sqrt(dt) / _noises.sigma_a;   // :157 加速度白化权重
```
- **白化 (whitening)**：连续噪声 $\sigma_w^2$ 经 `dt` 离散化后标准差为 $\sigma_w/\sqrt{dt}$，其逆权重即 $\sqrt{dt}/\sigma_w$（见 5.2）。

```cpp
res.block(6*i+0, 0, 3, 1) = -w_omega * w_hat;                       // :161 真值角速度 = 0
res.block(6*i+3, 0, 3, 1) = -w_accel * (a_hat - state->_imu->Rot() * _gravity);  // :163 比力 - R*g = 0
```
- 零速约束的**核心物理假设**：静止时
  - 真值角速度 $\boldsymbol{\omega}_{true} = \mathbf{w}_m - \mathbf{b}_g - \mathbf{n}_w = \mathbf{0}$
  - 真值加速度 $\mathbf{a}_{true} = \mathbf{a}_m - \mathbf{b}_a - \mathbf{R}\mathbf{g} - \mathbf{n}_a = \mathbf{0}$
- 残差即"观测真值（应为 0）减去模型预测"，前面乘白化权重。

#### 3.2.3 雅可比 (`:168-178`)

```cpp
Eigen::Matrix3d R_GtoI_jacob = (state->_options.do_fej) ? state->_imu->Rot_fej() : state->_imu->Rot();  // :169 FEJ
H.block(6*i+0, 3, 3, 3) = -w_omega * Eigen::Matrix3d::Identity();   // :170  ∂(-w_omega*w_hat)/∂bg
H.block(6*i+3, 0, 3, 3) = -w_accel * skew_x(R_GtoI_jacob * _gravity);  // :172 ∂/∂q (via R*g)
H.block(6*i+3, 6, 3, 3) = -w_accel * Eigen::Matrix3d::Identity();      // :173 ∂/∂ba
```
- 对 `bg`：角速度残差直接对零偏求负（权重 `w_omega`）。
- 对 `q`：加速度残差里 `R*g` 对姿态的雅可比是 `-skew_x(R*g)`（罗德里格斯一阶，[40]）。
- 对 `ba`：加速度残差对加速度零偏求负。
- 全部乘对应白化权重。

### 3.3 压缩 + 噪声放大 (`:182-190`)

```cpp
UpdaterHelper::measurement_compress_inplace(H, res);   // :183 → S13c 过定系统行压缩
if (H.rows() < 1) return false;                        // :184
Eigen::MatrixXd R = _zupt_noise_multiplier * Eigen::MatrixXd::Identity(res.rows(), res.rows());  // :190
```
- `measurement_compress_inplace`（S13c）把 `6*(N-1)` 行压缩到最小行数（去掉线性相关行），得到满秩系统。
- **`_zupt_noise_multiplier`**：故意把测量噪声 `R` 放大（`>1`），使 ZUPT 不"过度自信"——IMU 零速假设本身有建模误差，放大噪声可防止滤波器对这一更新过于信任（注释 :188-189）。

### 3.4 偏置前向演化 Q_bias (`:192-196`)

```cpp
Eigen::MatrixXd Q_bias = Eigen::MatrixXd::Identity(6, 6);
Q_bias.block(0,0,3,3) *= dt_summed * _noises.sigma_wb_2;   // :195 陀螺随机游走
Q_bias.block(3,3,3,3) *= dt_summed * _noises.sigma_ab_2;   // :196 加速度随机游走
```
- 在 χ² 检查前，把 bias 的随机游走协方差 `dt_summed * σ_b²` 加到 `P_marg` 的 bias 块（`:202-204`），等价于"假设在接受更新前 bias 已按随机游走演化了 `dt_summed`"。
- 注释推导：`G*Qd*G^T = dt*(1/dt*Qc)*dt = dt*Qc`（连续→离散），即 `dt_summed * σ_b²`。

### 3.5 χ² 门限 (`:198-216`)

```cpp
Eigen::MatrixXd P_marg = StateHelper::get_marginal_covariance(state, Hx_order);  // :201
if (model_time_varying_bias) P_marg.block(3,3,6,6) += Q_bias;                    // :202-204
Eigen::MatrixXd S = H * P_marg * H.transpose() + R;                              // :205
double chi2 = res.dot(S.llt().solve(res));                                       // :206
double chi2_check = (res.rows() < 1000) ? chi_squared_table[res.rows()]
                                        : quantile(chi_squared(res.rows()), 0.95);  // :210-215
```
- `chi_squared_table` 在构造函数里预计算到 1000 自由度（0.95 置信，`:57-62`）。
- `chi2` 若超过 `_options.chi2_multipler * chi2_check` 说明"残差异常大"——即载体其实在动，零速假设不成立，应拒绝。

### 3.6 视差覆盖检查 (`:218-237`)

```cpp
bool disparity_passed = false;
if (override_with_disparity_check) {
    FeatureHelper::compute_disparity(_db, time0_cam, time1_cam, disp_avg, disp_var, num_features);  // :228
    disparity_passed = (disp_avg < _zupt_max_disparity && num_features > 20);  // :231
}
```
- 若图像平均视差 `disp_avg` 超过 `_zupt_max_disparity`（且特征数 > 20），说明载体在动，`disparity_passed=true` 会**强制拒绝** ZUPT（即使 χ² 勉强过），防止运动时被误判静止。

### 3.7 接受判定与执行更新 (`:241-277`)

```cpp
if (!disparity_passed && (chi2 > _options.chi2_multipler*chi2_check
                          || state->_imu->vel().norm() > _zupt_max_velocity)) {  // :241
    last_zupt_state_timestamp = 0.0; last_zupt_count = 0;
    return false;   // 拒绝：在动 / χ² 过大
}
```
- **三个拒绝条件**（任一即拒）：视差通过（在动）、χ² 超阈、速度模长 > `_zupt_max_velocity`。
- 接受路径（`:263-277`）：
  ```cpp
  if (model_time_varying_bias) {
      // 先把 bias 按随机游走演化 (Phi=I, Q=Q_bias)
      StateHelper::EKFPropagation(state, Phi_order, Phi_order, Phi_bias, Q_bias);  // :272 → S9
  }
  StateHelper::EKFUpdate(state, Hx_order, H, res, R);   // :276 → S9
  state->_timestamp = timestamp;                        // :277
  ```
- 即：先 `EKFPropagation` 演化 bias 协方差（若启用时变 bias），再 `EKFUpdate` 用零速残差收缩 `q/bg/ba`（及 `v`）的协方差与均值——把漂移的 velocity 拉向 0、零偏收敛。

> 另有 `explicitly_enforce_zero_motion=true` 分支（`:279-324`，代码标注为可选/非常规）：先 `propagate_and_clone` 前推，再用 `[ori, pos, vel]` 残差显式强制零运动并边缘化临时 clone。默认 `false`，本篇不展开。

---

## Section 4: 接口与数据结构

### 4.1 签名

```cpp
bool UpdaterZeroVelocity::try_update(std::shared_ptr<State> state, double timestamp);  // :65
```

| 参数 | 类型 | 方向 | 说明 |
|------|------|------|------|
| `state` | `std::shared_ptr<State>` | in/out | 读 `_imu`/`_clones_IMU`/标定型；写 velocity/bias 均值与协方差 |
| `timestamp` | `double` | in | 当前图像时间戳（ZUPT 要推进到的时刻） |
| 返回 | `bool` | out | `true`=本次触发并接受了 ZUPT；`false`=未触发/被拒绝 |

### 4.2 构造函数缓存

```cpp
// 构造函数 :42-63
_gravity << 0.0, 0.0, gravity_mag;                     // :49
_noises.sigma_w_2 / sigma_a_2 / sigma_wb_2 / sigma_ab_2 = pow(...,2);  // :52-55
for (int i=1; i<1000; i++)                             // :59-61
    chi_squared_table[i] = quantile(chi_squared(i), 0.95);
```
- `chi_squared_table` 预存 1–999 自由度、0.95 置信的 χ² 阈值，避免运行时反复计算。

### 4.3 关键成员变量

| 成员 | 含义 |
|------|------|
| `_zupt_max_velocity` | 速度模长上限，超过即认为在动（拒绝 ZUPT） |
| `_zupt_noise_multiplier` | 测量噪声放大系数（`R` 乘子，防止过度自信） |
| `_zupt_max_disparity` | 图像平均视差上限，`override_with_disparity_check` 时用 |
| `imu_data` | 外部喂入的 IMU 缓冲（由 `VioManager` 维护） |

---

## Section 5: 数学与公式对照

### 5.1 零速约束模型

静止时 IMU 测量满足：
$$
\boldsymbol{\omega}_{true} = \mathbf{w}_m - \mathbf{b}_g - \mathbf{n}_w = \mathbf{0}
$$
$$
\mathbf{a}_{true} = \mathbf{a}_m - \mathbf{b}_a - \mathbf{R}\,\mathbf{g} - \mathbf{n}_a = \mathbf{0}
$$
残差（白化后）：
$$
\mathbf{r}_{\omega,i} = -w_{\omega,i}\,\hat{\boldsymbol{\omega}}_i,\qquad
\mathbf{r}_{a,i} = -w_{a,i}\,(\hat{\mathbf{a}}_i - \mathbf{R}\,\mathbf{g})
$$
其中 $\hat{\boldsymbol{\omega}}_i,\hat{\mathbf{a}}_i$ 为去外参/标度/零偏后的"真值估计"（`:148-149`）。

### 5.2 白化权重（Maybeck Eq.7-21a-c, 注释 :153-155）

连续噪声功率谱密度 $\mathbf{Q}_c = \mathrm{diag}(\sigma_w^2,\sigma_a^2)$。离散化后单段方差 $\sigma^2/dt$，其逆（白化）权重为 $\sqrt{dt}/\sigma$：
$$
w_{\omega,i} = \frac{\sqrt{dt_i}}{\sigma_w},\qquad w_{a,i} = \frac{\sqrt{dt_i}}{\sigma_a}
$$
乘上权重等价于对 $\mathbf{R}_{meas}^{-1} = \mathbf{L}\mathbf{L}^\top$ 做 Cholesky 白化（与标准 EKF 更新等价）。

### 5.3 雅可比

对姿态（用 FEJ 姿态 `R_GtoI_jacob`，`:169`）：
$$
\frac{\partial (-\mathbf{R}\,\mathbf{g})}{\partial \boldsymbol{\theta}} = -\lfloor \mathbf{R}\,\mathbf{g} \rfloor_\times
$$
（罗德里格斯，JPL 小角扰动，见 [40]）。对零偏：
$$
\frac{\partial \mathbf{r}_{\omega}}{\partial \mathbf{b}_g} = -w_\omega \mathbf{I},\qquad
\frac{\partial \mathbf{r}_a}{\partial \mathbf{b}_a} = -w_a \mathbf{I}
$$

### 5.4 χ² 门限

$$
\mathbf{S} = \mathbf{H}\,\mathbf{P}_{marg}\,\mathbf{H}^\top + \mathbf{R},\qquad
\chi^2 = \mathbf{r}^\top \mathbf{S}^{-1} \mathbf{r}
$$
$\chi^2 > \chi^2_{0.95}(\text{dof})\times\text{multiplier}$ 则拒绝（疑似在动）。

### 5.5 偏置演化（注释 :193, :266）

$$
\mathbf{Q}_{bias} = dt_{sum}\,\mathrm{diag}(\sigma_{wb}^2\mathbf{I}, \sigma_{ab}^2\mathbf{I})
$$
离散随机游走：$G Q_d G^\top = dt\cdot Q_c = dt\cdot\sigma_b^2$。

---

## Section 6: 关键实现细节与易错点

1. **白化权重含 $\sqrt{dt}$**（`:156-157`）：不是简单的 `1/sigma`，必须乘段时长才能正确离散化连续噪声。
2. **`R = _zupt_noise_multiplier * I`**（`:190`）：放大噪声是刻意设计，避免 ZUPT 过度收敛；调参时应关注该乘子。
3. **FEJ 姿态**（`:169`）：雅可比用 `Rot_fej()`（若 `do_fej`），与全仓 FEJ 一致性约定一致（见 S1/S2）。
4. **`gyro g-sensitivity` `Tg*a_hat`**（`:149`）：`w_hat` 减去 `Tg*a_hat`，建模了加速度对陀螺读数的耦合（[27]），易被忽略。
5. **视差覆盖**（`override_with_disparity_check`，`:220`）：即使 χ² 勉强过，若图像视差大（在动）也会拒绝——这是 ZUPT 误触发的主要防线。
6. **默认 `explicitly_enforce_zero_motion=false`**（`:114`）：`:279-324` 的"显式强制零运动"分支是非默认/可选路径，常规运行走 `:263-277` 的 `EKFPropagation`+`EKFUpdate`。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **上游**：S5（`select_imu_readings` 选段）、S7（`predict_and_compute` 的 bias 校正框架，供对照 `a_hat/w_hat` 去偏逻辑）、S9（`EKFPropagation`/`EKFUpdate` 协方差更新，本篇 `:272`/`:276` 直接调用）。
- **对照 S13/S14/S15**：ZUPT 是**无特征**的纯 IMU 更新，残差来自运动学约束而非重投影；但它复用同一套 `measurement_compress_inplace`（S13c）与 `EKFUpdate`（S9）框架，可并入 D 阶段残差框架统一看。
- **下游调用方**：S17（`VioManager` 静止检测判定后触发）。

### 待详细补充项（核心复杂内部模块，后续需要时拆子文档）

- **S16a — `FeatureHelper::compute_disparity`**：图像间平均视差/方差的计算与 `num_features` 阈值逻辑（`UpdaterZeroVelocity.cpp:228` 调用，实现在 `ov_core/src/feat/FeatureHelper.*`）。
- **S16b — `explicitly_enforce_zero_motion` 分支**（`UpdaterZeroVelocity.cpp:279-324`）：`propagate_and_clone` 前推 + `[ori,pos,vel]` 残差强制零运动 + 临时 clone 边缘化，默认关闭，需深挖时补。
- **S16c — `EKFUpdate` 在 ZUPT 下的收缩行为**：见 S9。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"，主文档 + Sxxa/b/c 子文档 + 懒补充清单）。
