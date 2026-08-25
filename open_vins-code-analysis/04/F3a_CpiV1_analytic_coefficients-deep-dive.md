# CpiV1 解析积分系数推导深度解析 (F3a)

> **生成日期**：2026-08-25
> **前置文档**：`F3_CpiV1-continuous-preintegration-deep-dive.md`（CPI 主文档）
> **核心代码**：`ov_core/src/cpi/CpiV1.cpp:70-127`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Eckenhoff et al. "Continuous Preintegration Theory for Graph-based Visual-Inertial Navigation" IJRR 2019 Eq.(35)-(37)

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
CpiV1::feed_IMU()  (F3)
        │  IMU 段积分
        ▼
解析积分系数 f_1..f_4  ★ 本篇
   ├─ 输入: δt, ω̂ (角速度), â (加速度)
   ├─ 计算: f_1, f_2, f_3, f_4 (标量系数)
   └─ 输出: α_arg, β_arg (3×3 矩阵)
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | CPI 预积分均值计算 |
| **输入** | 时间间隔 $\delta t$、角速度 $\hat{\boldsymbol{\omega}}$、加速度 $\hat{\mathbf{a}}$ |
| **输出** | $\boldsymbol{\Xi}_1, \boldsymbol{\Xi}_2$ 的标量系数 $f_1, f_2, f_3, f_4$ |
| **调用方** | `CpiV1::feed_IMU()`（F3） |

### 1.2 一句话概括

`f_1..f_4` 是连续时间预积分中 $\boldsymbol{\Xi}_1, \boldsymbol{\Xi}_2$ 矩阵的标量系数，将矩阵积分转化为标量多项式计算，避免数值积分。

---

## 2. 数学模型

### 2.1 连续时间积分定义

CPI 预积分需要计算两个矩阵积分（见 F3 式 (F3-5)(F3-6)）：

$$\boldsymbol{\Xi}_1(\Delta t) = \int_0^{\Delta t} \Delta\mathbf{R}(\tau)\,d\tau \tag{F3a-1}$$

$$\boldsymbol{\Xi}_2(\Delta t) = \int_0^{\Delta t}\int_0^s \Delta\mathbf{R}(\tau)\,d\tau\,ds \tag{F3a-2}$$

其中 $\Delta\mathbf{R}(\tau) = \exp(-\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times \tau)$ 是相对旋转矩阵。

### 2.2 Rodrigues 公式展开

利用 Rodrigues 旋转公式：

$$\exp(-\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times \tau) = \mathbf{I} - \frac{\sin(\omega\tau)}{\omega}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \frac{1-\cos(\omega\tau)}{\omega^2}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F3a-3}$$

其中 $\omega = \|\hat{\boldsymbol{\omega}}\|$。

### 2.3 积分计算

将式 (F3a-3) 代入式 (F3a-1)：

$$\boldsymbol{\Xi}_1(\Delta t) = \int_0^{\Delta t} \left[\mathbf{I} - \frac{\sin(\omega\tau)}{\omega}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \frac{1-\cos(\omega\tau)}{\omega^2}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2\right] d\tau$$

逐项积分：

$$\int_0^{\Delta t} \mathbf{I}\,d\tau = \Delta t\,\mathbf{I}$$

$$\int_0^{\Delta t} \frac{\sin(\omega\tau)}{\omega}\,d\tau = \frac{1-\cos(\omega\Delta t)}{\omega^2}$$

$$\int_0^{\Delta t} \frac{1-\cos(\omega\tau)}{\omega^2}\,d\tau = \frac{\omega\Delta t - \sin(\omega\Delta t)}{\omega^3}$$

故：

$$\boldsymbol{\Xi}_1(\Delta t) = \Delta t\,\mathbf{I} + f_3\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + f_4\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F3a-4}$$

其中：

$$f_3 = -\frac{1-\cos(\omega\Delta t)}{\omega^2}, \quad f_4 = \frac{\omega\Delta t - \sin(\omega\Delta t)}{\omega^3} \tag{F3a-5}$$

类似地，对式 (F3a-2) 双重积分：

$$\boldsymbol{\Xi}_2(\Delta t) = \frac{\Delta t^2}{2}\mathbf{I} + f_1\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + f_2\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F3a-6}$$

其中：

$$f_1 = \frac{\omega\Delta t\cos(\omega\Delta t) - \sin(\omega\Delta t)}{\omega^3}, \quad f_2 = \frac{(\omega\Delta t)^2 - 2\cos(\omega\Delta t) - 2\omega\Delta t\sin(\omega\Delta t) + 2}{2\omega^4} \tag{F3a-7}$$

### 2.4 小角度近似

当 $\omega < 0.5°/\text{s}$（代码中 `small_w = (mag_w < 0.008726646)`），直接计算会导致数值不稳定（除以小量）。使用 Taylor 展开：

$$\sin(x) \approx x - \frac{x^3}{6} + \frac{x^5}{120} - \cdots$$

$$\cos(x) \approx 1 - \frac{x^2}{2} + \frac{x^4}{24} - \cdots$$

代入式 (F3a-5)(F3a-7)，保留主导项：

$$f_1 \approx -\frac{\Delta t^3}{3}, \quad f_2 \approx \frac{\Delta t^4}{8}, \quad f_3 \approx -\frac{\Delta t^2}{2}, \quad f_4 \approx \frac{\Delta t^3}{6} \tag{F3a-8}$$

---

## 3. 代码实现：逐行对照

### 3.1 系数计算（正常情况）

```cpp
// CpiV1.cpp:102-112
if (small_w) {
    // 小角度 Taylor 展开 (式 F3a-8)
    f_1 = -(pow(delta_t, 3) / 3);
    f_2 = (pow(delta_t, 4) / 8);
    f_3 = -(pow(delta_t, 2) / 2);
    f_4 = (pow(delta_t, 3) / 6);
} else {
    // 正常解析公式 (式 F3a-5, F3a-7)
    f_1 = (w_dt * cos_wt - sin_wt) / (pow(mag_w, 3));
    f_2 = (pow(w_dt, 2) - 2 * cos_wt - 2 * w_dt * sin_wt + 2) / (2 * pow(mag_w, 4));
    f_3 = -(1 - cos_wt) / pow(mag_w, 2);
    f_4 = (w_dt - sin_wt) / pow(mag_w, 3);
}
```

**理论对应**：式 (F3a-5)(F3a-7)(F3a-8)。

### 3.2 矩阵构造

```cpp
// CpiV1.cpp:116-122
// Ξ_2 (式 F3a-6)
Eigen::Matrix<double, 3, 3> alpha_arg = ((dt_2 / 2.0) * eye3 + f_1 * w_x + f_2 * w_x_2);
// Ξ_1 (式 F3a-4)
Eigen::Matrix<double, 3, 3> Beta_arg = (delta_t * eye3 + f_3 * w_x + f_4 * w_x_2);

// 乘以 R^T 和 a_hat
Eigen::MatrixXd H_al = R_tau12k * alpha_arg;
Eigen::MatrixXd H_be = R_tau12k * Beta_arg;

// 更新均值 (式 F3-7)
alpha_tau += beta_tau * delta_t + H_al * a_hat;
beta_tau += H_be * a_hat;
```

**理论对应**：式 (F3-7)。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\boldsymbol{\Xi}_1$ | 速度增量积分矩阵 | `Beta_arg` | `:117` |
| $\boldsymbol{\Xi}_2$ | 位置增量积分矩阵 | `alpha_arg` | `:116` |
| $f_1$ | $\boldsymbol{\Xi}_2$ 的 $\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times$ 系数 | `f_1` | `:103-108` |
| $f_2$ | $\boldsymbol{\Xi}_2$ 的 $\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2$ 系数 | `f_2` | `:104-109` |
| $f_3$ | $\boldsymbol{\Xi}_1$ 的 $\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times$ 系数 | `f_3` | `:105-110` |
| $f_4$ | $\boldsymbol{\Xi}_1$ 的 $\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2$ 系数 | `f_4` | `:106-111` |

---

## 5. 关键设计决策

### 5.1 为什么用解析积分而非数值积分

| 方法 | 优点 | 缺点 |
|------|------|------|
| **解析积分**（代码采用） | 精确、快速、无累积误差 | 公式复杂 |
| 数值积分（RK4 等） | 简单、通用 | 计算量大、有累积误差 |

**选择理由**：
1. CPI 需要高精度（用于图优化残差）
2. 解析公式只需计算 4 个标量，极快
3. 无数值积分的累积误差

### 5.2 为什么需要小角度分支

当 $\omega \to 0$ 时：
- $f_1 \propto 1/\omega^3$ → 除以零
- $f_3 \propto 1/\omega^2$ → 除以零

Taylor 展开避免数值不稳定：

$$\lim_{\omega\to 0} f_3 = \lim_{\omega\to 0} -\frac{1-\cos(\omega\Delta t)}{\omega^2} = -\frac{\Delta t^2}{2} \tag{F3a-9}$$

### 5.3 阈值 0.5°/s 的选择

代码中 `small_w = (mag_w < 0.008726646)`，即 $0.5°/\text{s}$。

**选择理由**：
- 小于此值时，Taylor 展开误差 < 数值精度
- 大于此值时，解析公式数值稳定

---

## 6. 完整流程图

```
输入: δt, ω̂, â
        │
        ├─ 1. 计算 ω = ‖ω̂‖
        │
        ├─ 2. 判断小角度
        │     ─ small_w = (ω < 0.008726646)
        │
        ├─ 3. 计算系数 f_1..f_4
        │     ├─ if small_w: Taylor 展开 (式 F3a-8)
        │     └─ else: 解析公式 (式 F3a-5, F3a-7)
        │
        ├─ 4. 构造矩阵
        │     ├─ α_arg = (δt²/2)I + f₁[ω̂]× + f₂[ω]×²
        │     └─ β_arg = δt·I + f[ω̂]× + f₄[ω̂]×²
        │
        └─ 5. 更新均值
              ├─ H_al = R^T * α_arg
              ├─ H_be = R^T * β_arg
              ├─ α += β*δt + H_al*â
              └─ β += H_be*â
```

---

## 7. 与其他文档的关联

- **上游**：`F3_CpiV1-continuous-preintegration-deep-dive.md`（调用本函数）
- **下游**：`F2b_Factor_ImuCPIv1-deep-dive.md`（使用 CPI 测量）
- **相关**：Eckenhoff et al. IJRR 2019 Eq.(35)-(37)

### 待详细补充项

- **F3a-1**：Taylor 展开的高阶项分析（精度评估）
- **F3a-2**：小角度阈值的数值实验验证

> 以上子文档暂不展开，待需要逐项深挖时再补。
