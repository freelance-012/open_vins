# CpiV1 偏置 Jacobian 推导深度解析 (F3b)

> **生成日期**：2026-08-25
> **前置文档**：`F3_CpiV1-continuous-preintegration-deep-dive.md`（CPI 主文档）、`F3a_CpiV1_analytic_coefficients-deep-dive.md`（解析系数）
> **核心代码**：`ov_core/src/cpi/CpiV1.cpp:128-238`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Eckenhoff et al. IJRR 2019 Eq.(49)-(63), Eq.(81), Eq.(84), Eq.(164)

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
CpiV1::feed_IMU()  (F3)
        │  均值积分完成
        ▼
偏置 Jacobian 计算  ★ 本篇
   ├─ 输入: Δt, ω̂, â, J_q, J_α, J_β, H_α, H_β, R_k2tau
   ├─ 计算: 右 Jacobian J_r, 偏置导数 df_i/dbw_j
   └─ 输出: 更新后的 J_q, J_α, J_β, H_α, H_β
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | CPI 预积分偏置 Jacobian 更新 |
| **输入** | 时间间隔、角速度、加速度、现有 Jacobian、旋转矩阵 |
| **输出** | 更新后的 bias Jacobian（5 个 3×3 矩阵） |
| **调用方** | `CpiV1::feed_IMU()`（F3） |
| **被调用方** | `Factor_ImuCPIv1`（F2b，用于一阶偏置修正） |

### 1.2 一句话概括

CpiV1 偏置 Jacobian 计算 CPI 预积分测量对陀螺零偏 $\mathbf{b}_\omega$ 和加表零偏 $\mathbf{b}_a$ 的一阶导数，使 Ceres 因子能在 bias 变化时进行一阶修正而不必重新积分。

---

## 2. 数学模型

### 2.1 偏置 Jacobian 定义

CPI 预积分测量对 bias 的 Jacobian：

$$\mathbf{J}_q = \frac{\partial\Delta\mathbf{R}}{\partial\mathbf{b}_\omega} \in \mathbb{R}^{3\times 3} \tag{F3b-1}$$

$$\mathbf{J}_\alpha = \frac{\partial\boldsymbol{\alpha}}{\partial\mathbf{b}_\omega} \in \mathbb{R}^{3\times 3}, \quad \mathbf{H}_\alpha = \frac{\partial\boldsymbol{\alpha}}{\partial\mathbf{b}_a} \in \mathbb{R}^{3\times 3} \tag{F3b-2}$$

$$\mathbf{J}_\beta = \frac{\partial\boldsymbol{\beta}}{\partial\mathbf{b}_\omega} \in \mathbb{R}^{3\times 3}, \quad \mathbf{H}_\beta = \frac{\partial\boldsymbol{\beta}}{\partial\mathbf{b}_a} \in \mathbb{R}^{3\times 3} \tag{F3b-3}$$

### 2.2 旋转对 $\mathbf{b}_\omega$ 的 Jacobian

由 $\Delta\mathbf{R} = \exp(-\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times \Delta t)$，对 $\mathbf{b}_\omega$ 求导（利用 $\hat{\boldsymbol{\omega}} = \boldsymbol{\omega}_m - \mathbf{b}_\omega$）：

$$\mathbf{J}_q^{(k+1)} = \Delta\mathbf{R}_{k\to k+1}\,\mathbf{J}_q^{(k)} + \mathbf{J}_r(-\hat{\boldsymbol{\omega}}\Delta t)\,\Delta t \tag{F3b-4}$$

其中 $\mathbf{J}_r(\boldsymbol{\phi})$ 是 $SO(3)$ 右 Jacobian（Eckenhoff IJRR 2019 Eq.164）：

$$\mathbf{J}_r(\boldsymbol{\phi}) = \mathbf{I} - \frac{1-\cos\|\boldsymbol{\phi}\|}{\|\boldsymbol{\phi}\|^2}\lfloor\boldsymbol{\phi}\rfloor_\times + \frac{\|\boldsymbol{\phi}\|-\sin\|\boldsymbol{\phi}\|}{\|\boldsymbol{\phi}\|^3}\lfloor\boldsymbol{\phi}\rfloor_\times^2 \tag{F3b-5}$$

小角度近似（$\|\boldsymbol{\phi}\| < 0.5°$）：

$$\mathbf{J}_r(\boldsymbol{\phi}) \approx \mathbf{I} - \frac{1}{2}\lfloor\boldsymbol{\phi}\rfloor_\times + \frac{1}{6}\lfloor\boldsymbol{\phi}\rfloor_\times^2 \tag{F3b-6}$$

### 2.3 增量对 $\mathbf{b}_a$ 的 Jacobian

由均值更新公式（F3 式 (F3-7)）：

$$\boldsymbol{\alpha} \mathrel{+}= \boldsymbol{\beta}\Delta t + \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_2\hat{\mathbf{a}}$$

$$\boldsymbol{\beta} \mathrel{+}= \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_1\hat{\mathbf{a}}$$

对 $\mathbf{b}_a$ 求导（注意 $\hat{\mathbf{a}} = \mathbf{a}_m - \mathbf{b}_a$，故 $\frac{\partial\hat{\mathbf{a}}}{\partial\mathbf{b}_a} = -\mathbf{I}$）：

$$\mathbf{H}_\alpha \mathrel{-}= \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_2, \quad \mathbf{H}_\alpha \mathrel{+}= \Delta t\,\mathbf{H}_\beta \tag{F3b-7}$$

$$\mathbf{H}_\beta \mathrel{-}= \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_1 \tag{F3b-8}$$

### 2.4 增量对 $\mathbf{b}_\omega$ 的 Jacobian

由于 $\boldsymbol{\Xi}_1, \boldsymbol{\Xi}_2$ 也依赖 $\hat{\boldsymbol{\omega}}$（通过 $f_1..f_4$），求导更复杂：

$$\mathbf{J}_\alpha \mathrel{+}= \mathbf{J}_\beta\Delta t + \frac{\partial(\mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_2\hat{\mathbf{a}})}{\partial\mathbf{b}_\omega} \tag{F3b-9}$$

利用乘积法则：

$$\frac{\partial(\mathbf{R}\boldsymbol{\Xi}_2\hat{\mathbf{a}})}{\partial\mathbf{b}_\omega} = \frac{\partial\mathbf{R}}{\partial\mathbf{b}_\omega}\boldsymbol{\Xi}_2\hat{\mathbf{a}} + \mathbf{R}\frac{\partial\boldsymbol{\Xi}_2}{\partial\hat{\boldsymbol{\omega}}}\frac{\partial\hat{\boldsymbol{\omega}}}{\partial\mathbf{b}_\omega}\hat{\mathbf{a}} + \mathbf{R}\boldsymbol{\Xi}_2\frac{\partial\hat{\mathbf{a}}}{\partial\mathbf{b}_\omega}$$

其中 $\frac{\partial\hat{\boldsymbol{\omega}}}{\partial\mathbf{b}_\omega} = -\mathbf{I}$，$\frac{\partial\hat{\mathbf{a}}}{\partial\mathbf{b}_\omega} = \mathbf{0}$（假设 $\hat{\mathbf{a}}$ 不依赖 $\mathbf{b}_\omega$）。

故：

$$\mathbf{J}_\alpha \mathrel{+}= \mathbf{J}_\beta\Delta t + \frac{\partial\mathbf{R}}{\partial\mathbf{b}_\omega}\boldsymbol{\Xi}_2\hat{\mathbf{a}} - \mathbf{R}\frac{\partial\boldsymbol{\Xi}_2}{\partial\hat{\boldsymbol{\omega}}}\hat{\mathbf{a}} \tag{F3b-10}$$

### 2.5 $\frac{\partial\mathbf{R}}{\partial\mathbf{b}_\omega}$ 的计算

由 $\Delta\mathbf{R} = \exp(-\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times \Delta t)$，对 $\mathbf{b}_\omega$ 的第 $j$ 个分量求导：

$$\frac{\partial\Delta\mathbf{R}}{\partial b_{\omega,j}} = -\Delta\mathbf{R}\,\lfloor\mathbf{J}_q\mathbf{e}_j\rfloor_\times \tag{F3b-11}$$

代码实现（[CpiV1.cpp:149-153](src/open_vins/ov_core/src/cpi/CpiV1.cpp#L149-L153)）：

```cpp
Eigen::MatrixXd d_R_bw_1 = -R_tau12k * skew_x(J_q * e_1);
Eigen::MatrixXd d_R_bw_2 = -R_tau12k * skew_x(J_q * e_2);
Eigen::MatrixXd d_R_bw_3 = -R_tau12k * skew_x(J_q * e_3);
```

### 2.6 $\frac{\partial f_i}{\partial\hat{\boldsymbol{\omega}}}$ 的计算

$f_1..f_4$ 是 $\omega = \|\hat{\boldsymbol{\omega}}\|$ 的函数，故：

$$\frac{\partial f_i}{\partial b_{\omega,j}} = \frac{\partial f_i}{\partial\omega}\frac{\partial\omega}{\partial b_{\omega,j}} = \frac{\partial f_i}{\partial\omega}\left(-\frac{\hat{\omega}_j}{\omega}\right) = \hat{\omega}_j \cdot \frac{df_i}{d\omega} \cdot \left(-\frac{1}{\omega}\right)$$

定义 $\frac{df_i}{d\omega_{\text{mag}}} = \frac{df_i}{d\omega} \cdot \left(-\frac{1}{\omega}\right)$，则：

$$\frac{\partial f_i}{\partial b_{\omega,j}} = \hat{\omega}_j \cdot \frac{df_i}{d\omega_{\text{mag}}} \tag{F3b-12}$$

#### 2.6.1 正常情况（$\omega \geq 0.5°/\text{s}$）

以 $f_3 = -\frac{1-\cos(\omega\Delta t)}{\omega^2}$ 为例：

$$\frac{df_3}{d\omega_{\text{mag}}} = \frac{2(\cos(\omega\Delta t) - 1) + \omega\Delta t\sin(\omega\Delta t)}{\omega^4} \tag{F3b-13}$$

代码实现（[CpiV1.cpp:209-210](src/open_vins/ov_core/src/cpi/CpiV1.cpp#L209-L210)）：

```cpp
double df_3_dw_mag = (2 * (cos_wt - 1) + w_dt * sin_wt) / (pow(mag_w, 4));
df_3_dbw_1 = w_1 * df_3_dw_mag;
df_3_dbw_2 = w_2 * df_3_dw_mag;
df_3_dbw_3 = w_3 * df_3_dw_mag;
```

#### 2.6.2 小角度近似（$\omega < 0.5°/\text{s}$）

以 $f_3 \approx -\frac{\Delta t^2}{2}$ 为例：

$$\frac{df_3}{d\omega_{\text{mag}}} \approx -\frac{\Delta t^4}{12} \tag{F3b-14}$$

代码实现（[CpiV1.cpp:186-187](src/open_vins/ov_core/src/cpi/CpiV1.cpp#L186-L187)）：

```cpp
double df_3_dw_mag = -(pow(delta_t, 4) / 12);
df_3_dbw_1 = w_1 * df_3_dw_mag;
```

### 2.7 完整 Jacobian 更新公式

$$\mathbf{J}_\alpha \mathrel{+}= \mathbf{J}_\beta\Delta t + \sum_{j=1}^3 \left(\frac{\partial\mathbf{R}}{\partial b_{\omega,j}}\boldsymbol{\Xi}_2 + \mathbf{R}\left(\frac{\partial f_1}{\partial b_{\omega,j}}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \frac{\partial f_2}{\partial b_{\omega,j}}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 - f_1\lfloor\mathbf{e}_j\rfloor_\times - f_2(\lfloor\mathbf{e}_j\rfloor_\times\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \lfloor\hat{\boldsymbol{\omega}}\rfloor_\times\lfloor\mathbf{e}_j\rfloor_\times)\right)\hat{\mathbf{a}}\right) \tag{F3b-15}$$

$$\mathbf{J}_\beta \mathrel{+}= \sum_{j=1}^3 \left(\frac{\partial\mathbf{R}}{\partial b_{\omega,j}}\boldsymbol{\Xi}_1 + \mathbf{R}\left(\frac{\partial f_3}{\partial b_{\omega,j}}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \frac{\partial f_4}{\partial b_{\omega,j}}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 - f_3\lfloor\mathbf{e}_j\rfloor_\times - f_4(\lfloor\mathbf{e}_j\rfloor_\times\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \lfloor\hat{\boldsymbol{\omega}}\rfloor_\times\lfloor\mathbf{e}_j\rfloor_\times)\right)\hat{\mathbf{a}}\right) \tag{F3b-16}$$

---

## 3. 代码实现：逐段对照

### 3.1 右 Jacobian 计算

```cpp
// CpiV1.cpp:133-137
Eigen::Matrix<double, 3, 3> J_r_tau1 =
    small_w ? eye3 - .5 * w_tx + (1.0 / 6.0) * w_tx * w_tx   // 小角度 (式 F3b-6)
            : eye3 - ((1 - cos_wt) / (pow((w_dt), 2.0))) * w_tx + ((w_dt - sin_wt) / (pow(w_dt, 3.0))) * w_tx * w_tx;  // 正常 (式 F3b-5)
```

### 3.2 $\mathbf{J}_q$ 更新

```cpp
// CpiV1.cpp:141
J_q = R_tau2tau1 * J_q + J_r_tau1 * delta_t;  // 式 F3b-4
```

### 3.3 $\mathbf{H}_\alpha, \mathbf{H}_\beta$ 更新

```cpp
// CpiV1.cpp:143-147
H_a -= H_al;           // -R*Ξ_2 (式 F3b-7)
H_a += delta_t * H_b;  // +Δt*H_β
H_b -= H_be;           // -R*Ξ_1 (式 F3b-8)
```

### 3.4 $\frac{\partial\mathbf{R}}{\partial\mathbf{b}_\omega}$ 计算

```cpp
// CpiV1.cpp:149-153
Eigen::MatrixXd d_R_bw_1 = -R_tau12k * skew_x(J_q * e_1);  // 式 F3b-11
Eigen::MatrixXd d_R_bw_2 = -R_tau12k * skew_x(J_q * e_2);
Eigen::MatrixXd d_R_bw_3 = -R_tau12k * skew_x(J_q * e_3);
```

### 3.5 $\frac{df_i}{d\omega_{\text{mag}}}$ 计算

```cpp
// CpiV1.cpp:172-220
if (small_w) {
    // 小角度 Taylor 展开 (式 F3b-14 等)
    double df_1_dw_mag = -(pow(delta_t, 5) / 15);
    double df_2_dw_mag = (pow(delta_t, 6) / 72);
    double df_3_dw_mag = -(pow(delta_t, 4) / 12);
    double df_4_dw_mag = (pow(delta_t, 5) / 60);
} else {
    // 正常解析公式 (式 F3b-13 等)
    double df_1_dw_mag = (pow(w_dt, 2) * sin_wt - 3 * sin_wt + 3 * w_dt * cos_wt) / pow(mag_w, 5);
    double df_2_dw_mag = (pow(w_dt, 2) - 4 * cos_wt - 4 * w_dt * sin_wt + pow(w_dt, 2) * cos_wt + 4) / (pow(mag_w, 6));
    double df_3_dw_mag = (2 * (cos_wt - 1) + w_dt * sin_wt) / (pow(mag_w, 4));
    double df_4_dw_mag = (2 * w_dt + w_dt * cos_wt - 3 * sin_wt) / (pow(mag_w, 5));
}

// 链式法则 (式 F3b-12)
df_1_dbw_1 = w_1 * df_1_dw_mag;
df_1_dbw_2 = w_2 * df_1_dw_mag;
df_1_dbw_3 = w_3 * df_1_dw_mag;
// ... 类似计算 df_2, df_3, df_4
```

### 3.6 $\mathbf{J}_\alpha, \mathbf{J}_\beta$ 完整更新

```cpp
// CpiV1.cpp:222-238
// J_α 更新 (式 F3b-15)
J_a += J_b * delta_t;
J_a.block(0, 0, 3, 1) +=
    (d_R_bw_1 * alpha_arg + R_tau12k * (df_1_dbw_1 * w_x - f_1 * e_1x + df_2_dbw_1 * w_x_2 - f_2 * (e_1x * w_x + w_x * e_1x))) * a_hat;
J_a.block(0, 1, 3, 1) +=
    (d_R_bw_2 * alpha_arg + R_tau12k * (df_1_dbw_2 * w_x - f_1 * e_2x + df_2_dbw_2 * w_x_2 - f_2 * (e_2x * w_x + w_x * e_2x))) * a_hat;
J_a.block(0, 2, 3, 1) +=
    (d_R_bw_3 * alpha_arg + R_tau12k * (df_1_dbw_3 * w_x - f_1 * e_3x + df_2_dbw_3 * w_x_2 - f_2 * (e_3x * w_x + w_x * e_3x))) * a_hat;

// J_β 更新 (式 F3b-16)
J_b.block(0, 0, 3, 1) +=
    (d_R_bw_1 * Beta_arg + R_tau12k * (df_3_dbw_1 * w_x - f_3 * e_1x + df_4_dbw_1 * w_x_2 - f_4 * (e_1x * w_x + w_x * e_1x))) * a_hat;
// ... 类似计算第 1, 2 列
```

---

## 4. 代码 - 公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{J}_r$ | 右 Jacobian | `J_r_tau1` | `:135-137` |
| $\mathbf{J}_q$ | 旋转对 $b_\omega$ 的 Jacobian | `J_q` | `:141` |
| $\mathbf{H}_\alpha$ | α 对 $b_a$ 的 Jacobian | `H_a` | `:145-147` |
| $\mathbf{H}_\beta$ | β 对 $b_a$ 的 Jacobian | `H_b` | `:147` |
| $\frac{\partial\mathbf{R}}{\partial b_{\omega,j}}$ | 旋转对 $b_{\omega,j}$ 的导数 | `d_R_bw_1/2/3` | `:151-153` |
| $\frac{df_i}{d\omega_{\text{mag}}}$ | 系数对 $\omega$ 的导数 | `df_i_dw_mag` | `:174-219` |
| $\frac{\partial f_i}{\partial b_{\omega,j}}$ | 系数对 $b_{\omega,j}$ 的导数 | `df_i_dbw_j` | `:175-219` |

---

## 5. 关键设计决策

### 5.1 为什么需要偏置 Jacobian

CPI 预积分是在固定线性化点 $(\mathbf{b}_\omega^{\text{lin}}, \mathbf{b}_a^{\text{lin}})$ 计算的。当 bias 变化时：

| 方法 | 优点 | 缺点 |
|------|------|------|
| **一阶修正**（代码采用） | 快速、只需矩阵乘法 | 线性近似误差 |
| 重新积分 | 精确 | 计算量大（需重新 feed_IMU） |

**选择理由**：Ceres 优化中 bias 变化小，一阶近似足够精确。

### 5.2 为什么分小角度分支

当 $\omega \to 0$ 时，$\frac{df_i}{d\omega}$ 的表达式涉及除以 $\omega^4, \omega^5, \omega^6$，数值不稳定。Taylor 展开避免此问题。

### 5.3 乘积法则的展开

式 (F3b-15)(F3b-16) 中的项：

$$\frac{\partial(\mathbf{R}\boldsymbol{\Xi}\hat{\mathbf{a}})}{\partial\mathbf{b}_\omega} = \frac{\partial\mathbf{R}}{\partial\mathbf{b}_\omega}\boldsymbol{\Xi}\hat{\mathbf{a}} + \mathbf{R}\frac{\partial\boldsymbol{\Xi}}{\partial\hat{\boldsymbol{\omega}}}\frac{\partial\hat{\boldsymbol{\omega}}}{\partial\mathbf{b}_\omega}\hat{\mathbf{a}}$$

代码中展开为：

```cpp
d_R_bw * alpha_arg  // ∂R/∂b_ω * Ξ * â
+ R * (df_dbw * w_x - f * e_x + df_dbw * w_x_2 - f * (e_x * w_x + w_x * e_x)) * a_hat
//   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//   ∂Ξ/∂ω̂ * ∂ω̂/∂b_ω * â (乘积法则展开)
```

---

## 6. 完整流程图

```
输入: Δt, ω, â, J_q, J_α, J_β, H_α, H_β, R_k2tau
        │
        ├─ 1. 计算右 Jacobian J_r
        │     └─ small_w ? Taylor : 解析公式
        │
        ├─ 2. 更新 J_q
        │     ─ J_q = ΔR * J_q + J_r * Δt
        │
        ├─ 3. 更新 H_α, H_β
        │     ├─ H_α -= R*Ξ_2
        │     ├─ H_α += Δt * H_β
        │     └─ H_β -= R*Ξ_1
        │
        ├─ 4. 计算 ∂R/∂b_ω
        │     ─ d_R_bw_j = -R * skew(J_q * e_j)
        │
        ├─ 5. 计算 df_i/dω_mag
        │     └─ small_w ? Taylor : 解析公式
        │
        ├─ 6. 计算 f_i/∂b_ω
        │     └─ df_i_dbw_j = w_j * df_i_dw_mag
        │
        └─ 7. 更新 J_α, J_β
              ─ 乘积法则展开 (式 F3b-15, F3b-16)
```

---

## 7. 与其他文档的关联

- **上游**：`F3a_CpiV1_analytic_coefficients-deep-dive.md`（提供 $f_1..f_4$ 系数）
- **下游**：`F2b_Factor_ImuCPIv1-deep-dive.md`（使用 Jacobian 进行一阶偏置修正）
- **相关**：Eckenhoff et al. IJRR 2019 Eq.(49)-(63), Eq.(81), Eq.(84), Eq.(164)

### 待详细补充项

- **F3b-1**：$\frac{\partial\boldsymbol{\Xi}}{\partial\hat{\boldsymbol{\omega}}}$ 的完整矩阵推导
- **F3b-2**：小角度 Taylor 展开的高阶项分析

> 以上子文档暂不展开，待需要逐项深挖时再补。
