# StateHelper::EKFUpdate 深度精读 (S9b)

> **生成日期**：2026-08-19
> **前置文档**：`S9_StateHelper-deep-dive.md`（EKFPropagation 协方差传播）、`S11_StateHelper-marginalize-deep-dive.md`（状态边缘化）
> **核心代码**：`ov_msckf/src/state/StateHelper.cpp:116-197`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Trawny & Roumeliotis "Indirect Kalman Filter for 3D Attitude Estimation" 2005；Anderson & Moore "Optimal Filtering" 1979

---

## 1. 模块定位与职责

### 1.1 在系统中的位置

`EKFUpdate` 是 OpenVINS 中**最核心的更新函数**，被多个模块调用：

```
UpdaterMSCKF::update()          ─┐
UpdaterSLAM::update()            ├─→ EKFUpdate(state, H_order, H, res, R)  @ :116  ★ 本篇
UpdaterSLAM::delayed_init()     ─┘
UpdaterZeroVelocity::update()   ──→ EKFUpdate(...)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 所有 EKF 更新的最后一步 |
| 输入 | `state`（当前状态）、`H_order`（变量列表）、`H`（Jacobian）、`res`（残差）、`R`（噪声协方差） |
| 输出 | 更新后的 `state->_Cov` 和 `state->_variables`（均值） |
| 调用方 | `UpdaterMSCKF`、`UpdaterSLAM`、`UpdaterZeroVelocity` 等所有 updater |

### 1.2 一句话概括

`EKFUpdate` 执行标准 EKF 更新步：计算卡尔曼增益 $\mathbf{K}$，更新状态均值 $\delta\mathbf{x} = \mathbf{K}\mathbf{r}$，更新协方差 $\mathbf{P}' = \mathbf{P} - \mathbf{K}\mathbf{S}\mathbf{K}^\top$。

---

## 2. 数学模型

### 2.1 问题定义

给定线性化测量方程：

$$\mathbf{r} = \mathbf{H}_x\,\delta\mathbf{x} + \mathbf{n}, \quad \mathbf{n} \sim \mathcal{N}(\mathbf{0},\, \mathbf{R}) \tag{S9b-1}$$

其中：
- $\mathbf{r} \in \mathbb{R}^k$：残差向量
- $\mathbf{H}_x \in \mathbb{R}^{k \times m}$：状态 Jacobian
- $\delta\mathbf{x} \in \mathbb{R}^m$：状态误差
- $\mathbf{R} \in \mathbb{R}^{k \times k}$：测量噪声协方差

### 2.2 残差协方差 S

$$\mathbf{S} = \mathbf{H}_x\,\mathbf{P}\,\mathbf{H}_x^\top + \mathbf{R} \tag{S9b-2}$$

其中 $\mathbf{P}$ 是状态协方差矩阵（$m \times m$）。

**物理意义**：$\mathbf{S}$ 是残差的总不确定性，包含：
- $\mathbf{H}_x\mathbf{P}\mathbf{H}_x^\top$：状态不确定性通过测量模型传播到残差空间
- $\mathbf{R}$：传感器本身的测量噪声

### 2.3 卡尔曼增益 K

$$\mathbf{K} = \mathbf{P}\,\mathbf{H}_x^\top\,\mathbf{S}^{-1} \tag{S9b-3}$$

**物理意义**：卡尔曼增益决定了"残差中有多少应该归因于状态误差"。

### 2.4 状态更新

$$\delta\mathbf{x} = \mathbf{K}\,\mathbf{r} \tag{S9b-4}$$

### 2.5 协方差更新

$$\mathbf{P}' = \mathbf{P} - \mathbf{K}\,\mathbf{S}\,\mathbf{K}^\top \tag{S9b-5}$$

**等价形式**：

$$\mathbf{P}' = (\mathbf{I} - \mathbf{K}\mathbf{H}_x)\,\mathbf{P} \tag{S9b-6}$$

这两种形式在数学上等价，但 (S9b-5) 数值更稳定（对称性更好）。

---

## 3. 代码实现：逐行对照理论

### 3.1 代码全文

```cpp
// StateHelper.cpp:116-197
void StateHelper::EKFUpdate(state, H_order, H, res, R) {
    // 1. 计算 M = P * H^T
    Eigen::MatrixXd M_a = Eigen::MatrixXd::Zero(state->_Cov.rows(), res.rows());
    for (const auto &var : state->_variables) {
        Eigen::MatrixXd M_i = Eigen::MatrixXd::Zero(var->size(), res.rows());
        for (size_t i = 0; i < H_order.size(); i++) {
            std::shared_ptr<Type> meas_var = H_order[i];
            M_i += state->_Cov.block(var->id(), meas_var->id(), var->size(), meas_var->size())
                   * H.block(0, H_id[i], H.rows(), meas_var->size()).transpose();
        }
        M_a.block(var->id(), 0, var->size(), res.rows()) = M_i;
    }
    
    // 2. 计算 S = H * P_marg * H^T + R
    Eigen::MatrixXd P_small = get_marginal_covariance(state, H_order);
    Eigen::MatrixXd S(R.rows(), R.rows());
    S.triangularView<Eigen::Upper>() = H * P_small * H.transpose();
    S.triangularView<Eigen::Upper>() += R;
    
    // 3. 计算 K = M * S^{-1}
    Eigen::MatrixXd Sinv = Eigen::MatrixXd::Identity(R.rows(), R.rows());
    S.selfadjointView<Eigen::Upper>().llt().solveInPlace(Sinv);
    Eigen::MatrixXd K = M_a * Sinv.selfadjointView<Eigen::Upper>();
    
    // 4. 更新协方差 P' = P - K * S * K^T
    state->_Cov.triangularView<Eigen::Upper>() -= K * M_a.transpose();
    state->_Cov = state->_Cov.selfadjointView<Eigen::Upper>();
    
    // 5. 更新状态均值 δx = K * r
    Eigen::VectorXd dx = K * res;
    for (size_t i = 0; i < state->_variables.size(); i++) {
        state->_variables.at(i)->update(dx.block(state->_variables.at(i)->id(), 0,
                                                  state->_variables.at(i)->size(), 1));
    }
    
    // 6. 更新相机内参（如果在线标定）
    if (state->_options.do_calib_camera_intrinsics) {
        for (auto const &calib : state->_cam_intrinsics) {
            state->_cam_intrinsics_cameras.at(calib.first)->set_value(calib.second->value());
        }
    }
}
```

### 3.2 逐行代码 ↔ 理论对照

#### 3.2.1 计算 M = P * H^T

```cpp
// 分块累加：M_a[var_id, :] = Σ P[var_id, meas_var_id] * H[:, meas_var_col].T
for (const auto &var : state->_variables) {
    for (size_t i = 0; i < H_order.size(); i++) {
        M_i += state->_Cov.block(var->id(), meas_var->id(), var->size(), meas_var->size())
               * H.block(0, H_id[i], H.rows(), meas_var->size()).transpose();
    }
    M_a.block(var->id(), 0, var->size(), res.rows()) = M_i;
}
```

**理论对应**：公式 (S9b-3) 中的 $\mathbf{P}\mathbf{H}_x^\top$。

**为什么分块累加**：$\mathbf{P}$ 是 $n \times n$ 大矩阵（$n$ = 总状态维度），但 $\mathbf{H}_x$ 只涉及部分变量（`H_order`）。分块累加避免构造完整的 $\mathbf{H}_x$ 矩阵。

#### 3.2.2 计算 S = H * P_marg * H^T + R

```cpp
Eigen::MatrixXd P_small = get_marginal_covariance(state, H_order);
S.triangularView<Eigen::Upper>() = H * P_small * H.transpose();
S.triangularView<Eigen::Upper>() += R;
```

**理论对应**：公式 (S9b-2)。

**为什么用 `triangularView<Upper>()`**：$\mathbf{S}$ 是对称矩阵，只计算上三角部分，节省一半计算量。

#### 3.2.3 计算 K = M * S^{-1}

```cpp
Eigen::MatrixXd Sinv = Eigen::MatrixXd::Identity(R.rows(), R.rows());
S.selfadjointView<Eigen::Upper>().llt().solveInPlace(Sinv);
Eigen::MatrixXd K = M_a * Sinv.selfadjointView<Eigen::Upper>();
```

**理论对应**：公式 (S9b-3)。

**为什么用 LLT 分解**：$\mathbf{S}$ 对称正定，LLT（Cholesky 分解）比 LU 分解快一倍，且数值更稳定。

**`solveInPlace`**：直接修改 `Sinv`，避免额外内存分配。

#### 3.2.4 更新协方差 P' = P - K * S * K^T

```cpp
state->_Cov.triangularView<Eigen::Upper>() -= K * M_a.transpose();
state->_Cov = state->_Cov.selfadjointView<Eigen::Upper>();
```

**理论对应**：公式 (S9b-5)。

**注意**：代码用的是 $\mathbf{P}' = \mathbf{P} - \mathbf{K}\mathbf{M}^\top$（其中 $\mathbf{M} = \mathbf{P}\mathbf{H}_x^\top$），而非 $\mathbf{P} - \mathbf{K}\mathbf{S}\mathbf{K}^\top$。这两种形式等价：

$$\mathbf{K}\mathbf{M}^\top = \mathbf{K}(\mathbf{P}\mathbf{H}_x^\top)^\top = \mathbf{K}\mathbf{H}_x\mathbf{P} = \mathbf{K}\mathbf{S}\mathbf{K}^\top$$

**为什么用 `selfadjointView<Upper>()`**：强制对称化，消除数值误差导致的非对称性。

#### 3.2.5 更新状态均值 δx = K * r

```cpp
Eigen::VectorXd dx = K * res;
for (size_t i = 0; i < state->_variables.size(); i++) {
    state->_variables.at(i)->update(dx.block(state->_variables.at(i)->id(), 0,
                                              state->_variables.at(i)->size(), 1));
}
```

**理论对应**：公式 (S9b-4)。

**为什么循环更新**：每个 Type 变量有自己的 `update()` 方法（流形加法），例如：
- `Vec::update(dx)`：$\mathbf{x} \leftarrow \mathbf{x} + \delta\mathbf{x}$
- `JPLQuat::update(dx)`：$\bar{q} \leftarrow [\frac{1}{2}\delta\boldsymbol{\theta}^\top, 1]^\top \otimes \bar{q}$

#### 3.2.6 更新相机内参

```cpp
if (state->_options.do_calib_camera_intrinsics) {
    for (auto const &calib : state->_cam_intrinsics) {
        state->_cam_intrinsics_cameras.at(calib.first)->set_value(calib.second->value());
    }
}
```

**作用**：如果在线标定相机内参，更新后需要同步到相机模型对象。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{r}$ | 残差 | `res` | `:116` |
| $\mathbf{H}_x$ | 状态 Jacobian | `H` | `:116` |
| $\mathbf{R}$ | 测量噪声协方差 | `R` | `:116` |
| $\mathbf{P}$ | 状态协方差 | `state->_Cov` | `:116` |
| $\mathbf{M} = \mathbf{P}\mathbf{H}_x^\top$ | 中间矩阵 | `M_a` | `:119-128` |
| $\mathbf{S} = \mathbf{H}_x\mathbf{P}\mathbf{H}_x^\top + \mathbf{R}$ | 残差协方差 | `S` | `:131-134` |
| $\mathbf{K} = \mathbf{P}\mathbf{H}_x^\top\mathbf{S}^{-1}$ | 卡尔曼增益 | `K` | `:137-138` |
| $\delta\mathbf{x} = \mathbf{K}\mathbf{r}$ | 状态更新量 | `dx` | `:144` |

---

## 5. 关键设计决策

### 5.1 EKFUpdate vs EKFPropagation

| 特性 | EKFPropagation (S9) | EKFUpdate (S9b) |
|------|---------------------|-----------------|
| **阶段** | 预测步 | 更新步 |
| **输入** | $\boldsymbol{\Phi}, \mathbf{Q}_d$ | $\mathbf{H}_x, \mathbf{r}, \mathbf{R}$ |
| **输出** | $\mathbf{P}' = \boldsymbol{\Phi}\mathbf{P}\boldsymbol{\Phi}^\top + \mathbf{Q}_d$ | $\mathbf{P}' = \mathbf{P} - \mathbf{K}\mathbf{S}\mathbf{K}^\top$ |
| **状态均值** | 不变 | $\delta\mathbf{x} = \mathbf{K}\mathbf{r}$ |
| **计算量** | $O(n^2 m)$ | $O(n^2 k + k^3)$ |

### 5.2 为什么用 $\mathbf{P} - \mathbf{K}\mathbf{M}^\top$ 而非 $\mathbf{P} - \mathbf{K}\mathbf{S}\mathbf{K}^\top$

两种形式数学等价：

$$\mathbf{K}\mathbf{M}^\top = \mathbf{K}(\mathbf{P}\mathbf{H}_x^\top)^\top = \mathbf{K}\mathbf{H}_x\mathbf{P} = \mathbf{K}\mathbf{S}\mathbf{K}^\top$$

但 $\mathbf{P} - \mathbf{K}\mathbf{M}^\top$ 数值更稳定，因为：
1. $\mathbf{M} = \mathbf{P}\mathbf{H}_x^\top$ 已经计算过（用于计算 $\mathbf{K}$）
2. 避免重复计算 $\mathbf{K}\mathbf{S}\mathbf{K}^\top$（需要额外的矩阵乘法）

### 5.3 为什么用 LLT 而非 LU

$\mathbf{S}$ 是对称正定矩阵（$\mathbf{S} = \mathbf{H}_x\mathbf{P}\mathbf{H}_x^\top + \mathbf{R}$，$\mathbf{P} \succ 0, \mathbf{R} \succ 0$），LLT 分解：
- 速度：比 LU 快一倍（利用对称性）
- 稳定性：数值更稳定（不需要选主元）
- 存储：只需存储上三角或下三角

---

## 6. 代码位置速查

| 功能 | 代码行 |
|------|--------|
| 计算 $\mathbf{M} = \mathbf{P}\mathbf{H}_x^\top$ | `:119-128` |
| 计算 $\mathbf{S} = \mathbf{H}_x\mathbf{P}\mathbf{H}_x^\top + \mathbf{R}$ | `:131-134` |
| 计算 $\mathbf{K} = \mathbf{M}\mathbf{S}^{-1}$ | `:137-138` |
| 更新协方差 $\mathbf{P}' = \mathbf{P} - \mathbf{K}\mathbf{M}^\top$ | `:141-142` |
| 更新状态均值 $\delta\mathbf{x} = \mathbf{K}\mathbf{r}$ | `:144-148` |
| 更新相机内参 | `:151-155` |

---

## 7. 同类实现对比

| 特性 | OpenVINS | VINS-Mono | 标准 EKF |
|------|---------|-----------|---------|
| 协方差更新形式 | $\mathbf{P} - \mathbf{K}\mathbf{M}^\top$ | $\mathbf{P} - \mathbf{K}\mathbf{S}\mathbf{K}^\top$ | 两者皆可 |
| 对称化 | `selfadjointView<Upper>()` | 手动对称化 | 通常不需要 |
| 求逆方法 | LLT | LU | 视情况 |

### 设计取舍分析

- 选择 **$\mathbf{P} - \mathbf{K}\mathbf{M}^\top$** 是因为：复用已计算的 $\mathbf{M}$，减少一次矩阵乘法
- 选择 **LLT 分解** 是因为：$\mathbf{S}$ 对称正定，LLT 更快更稳定
- 选择 **`selfadjointView`** 是因为：强制对称化，消除数值误差

---

## 8. 完整流程图

### 8.1 EKFUpdate 内部流程

```
输入: state, H_order, H, res, R
        │
        ├─ 1. 计算 M = P * H^T
        │     └─ 分块累加：M_a[var_id, :] = Σ P[var_id, meas_var_id] * H[:, meas_var_col].T
        │
        ├─ 2. 计算 S = H * P_marg * H^T + R
        │     └─ 只计算上三角（对称矩阵）
        │
        ├─ 3. 计算 K = M * S^{-1}
        │     └─ LLT 分解求 S^{-1}
        │
        ├─ 4. 更新协方差 P' = P - K * M^T
        │     └─ 强制对称化
        │
        ├─ 5. 更新状态均值 δx = K * r
        │      每个 Type 变量调用自己的 update()（流形加法）
        │
        └─ 6. 更新相机内参（如果在线标定）
```

### 8.2 EKFUpdate 在系统中的调用关系

```
─────────────────────────────────────────────────────────────┐
│                      VioManager::update()                   │
└─────────────────────┬───────────────────────────────────────┘
                      │
        ┌─────────────┼─────────────┐
        ▼             ▼             ▼
┌───────────────┐ ┌──────────────┐ ┌──────────────────┐
│ UpdaterMSCKF  │ │ UpdaterSLAM  │ │ UpdaterZeroVel   │
│ ::update()    │ │ ::update()   │ │ ::update()       │
└───────┬───────┘ └──────┬───────┘ └─────────────────┘
        │                │                   │
        ────────────────┼───────────────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │  EKFUpdate()        │
              │  (StateHelper.cpp)  │
              └─────────────────────┘
```

### 8.3 各调用方的输入特点

| 调用方 | H 的构成 | res 的来源 | R 的特点 |
|--------|---------|-----------|---------|
| **UpdaterMSCKF** | 投影后的 $\tilde{\mathbf{H}}_{x2}$ | 投影后的残差 $\tilde{\mathbf{r}}_2$ | $\sigma_{\text{pix}}^2 \mathbf{I}$ |
| **UpdaterSLAM** | $[\mathbf{H}_x \mid \mathbf{H}_f]$ | 重投影残差 | 块对角（每特征独立） |
| **UpdaterZeroVelocity** | 速度/角速度选择矩阵 | 零速约束残差 | 小方差对角阵 |
