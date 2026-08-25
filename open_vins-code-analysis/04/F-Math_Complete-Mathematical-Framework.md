# OpenVINS 初始化模块完整数学框架 (F-Math)

> **生成日期**：2026-08-25
> **前置文档**：`F1_InertialInitializer-deep-dive.md`, `F2_Static-Dynamic-Initializer-deep-dive.md`, `F3_CpiV1-continuous-preintegration-deep-dive.md`
> **子文档索引**：
> - F1a: `F1a_StaticInitializer-algorithm-deep-dive.md`
> - F1b: `F1b_DynamicInitializer-algorithm-deep-dive.md`
> - F2a: `F2a_compute_dongsi_coeff-deep-dive.md`
> - F2b: `F2b_Factor_ImuCPIv1-deep-dive.md`
> - F2c: `F2c_Factor_ImageReprojCalib-deep-dive.md`
> - F3a: `F3a_CpiV1_analytic_coefficients-deep-dive.md`
> - F3b: `F3b_CpiV1_bias_jacobians-deep-dive.md` (待补充)
> - F2d: `F2d_InitializerHelper-deep-dive.md` (待补充)
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Dong-Si & Mourikis IROS 2012; Eckenhoff et al. IJRR 2019

---

## 1. 初始化问题数学定义

### 1.1 待估计状态

给定初始化窗口 $[t_0, t_N]$ 内的 IMU 和相机观测，需要估计：

**IMU 状态**（每帧）：
$$\mathbf{x}_{I_k} = \begin{bmatrix} {}^{I_k}_G\bar{q}^\top & {}^G\mathbf{p}_{I_k}^\top & {}^G\mathbf{v}_{I_k}^\top & \mathbf{b}_\omega^\top & \mathbf{b}_a^\top \end{bmatrix}^\top \in \mathbb{R}^{16} \tag{F-Math-1}$$

**特征位置**：
$$\{\mathbf{p}_{F_j}^G\}_{j=1}^M, \quad \mathbf{p}_{F_j}^G \in \mathbb{R}^3 \tag{F-Math-2}$$

**重力方向**：
$$\mathbf{g}^G \in \mathbb{R}^3, \quad \|\mathbf{g}^G\| = 9.81\,\text{m/s}^2 \tag{F-Math-3}$$

### 1.2 观测模型

**IMU 预积分测量**（相邻帧 $k, k+1$）：
$$\mathbf{z}_{\text{IMU}}^{k,k+1} = \begin{bmatrix} \Delta\mathbf{R}_{k,k+1} \\ \boldsymbol{\alpha}_{k,k+1} \\ \boldsymbol{\beta}_{k,k+1} \end{bmatrix} + \mathbf{n}_{\text{IMU}} \tag{F-Math-4}$$

**视觉重投影**（特征 $j$ 在帧 $k$）：
$$\mathbf{z}_{\text{cam}}^{j,k} = \pi(\mathbf{p}_{F_j}^{C_k}) + \mathbf{n}_{\text{cam}} \tag{F-Math-5}$$

其中 $\pi(\cdot)$ 是相机投影模型（含畸变）。

---

## 2. CPI 预积分数学框架

### 2.1 连续时间 IMU 运动学

$$\dot{\mathbf{R}}_{G\to I}(t) = \mathbf{R}_{G\to I}(t)\lfloor\boldsymbol{\omega}(t)\rfloor_\times \tag{F-Math-6}$$

$${}^G\dot{\mathbf{p}}_I(t) = {}^G\mathbf{v}_I(t) \tag{F-Math-7}$$

$${}^G\dot{\mathbf{v}}_I(t) = \mathbf{R}_{G\to I}^\top(t)\mathbf{a}(t) + {}^G\mathbf{g} \tag{F-Math-8}$$

其中 $\boldsymbol{\omega}(t) = \boldsymbol{\omega}_m(t) - \mathbf{b}_\omega - \mathbf{n}_\omega$，$\mathbf{a}(t) = \mathbf{a}_m(t) - \mathbf{b}_a - \mathbf{n}_a$。

### 2.2 预积分测量定义

在 $I_k$ 系中定义相对运动量：

$$\Delta\mathbf{R}_{k\to k+1} = \mathbf{R}_{G\to I_k}^\top\mathbf{R}_{G\to I_{k+1}} \tag{F-Math-9}$$

$$\boldsymbol{\alpha}_{k\to k+1} = \int_{t_k}^{t_{k+1}}\int_{t_k}^s \Delta\mathbf{R}_{k\to\tau}\,\mathbf{a}(\tau)\,d\tau\,ds \tag{F-Math-10}$$

$$\boldsymbol{\beta}_{k\to k+1} = \int_{t_k}^{t_{k+1}} \Delta\mathbf{R}_{k\to\tau}\,\mathbf{a}(\tau)\,d\tau \tag{F-Math-11}$$

### 2.3 解析积分系数

利用 Rodrigues 公式展开 $\Delta\mathbf{R}(\tau) = \exp(-\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times \tau)$：

$$\Delta\mathbf{R}(\tau) = \mathbf{I} - \frac{\sin(\omega\tau)}{\omega}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + \frac{1-\cos(\omega\tau)}{\omega^2}\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F-Math-12}$$

代入积分得：

$$\boldsymbol{\Xi}_1(\Delta t) = \Delta t\,\mathbf{I} + f_3\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + f_4\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F-Math-13}$$

$$\boldsymbol{\Xi}_2(\Delta t) = \frac{\Delta t^2}{2}\mathbf{I} + f_1\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times + f_2\lfloor\hat{\boldsymbol{\omega}}\rfloor_\times^2 \tag{F-Math-14}$$

其中系数（见 F3a）：

$$f_1 = \frac{\omega\Delta t\cos(\omega\Delta t) - \sin(\omega\Delta t)}{\omega^3}, \quad f_2 = \frac{(\omega\Delta t)^2 - 2\cos(\omega\Delta t) - 2\omega\Delta t\sin(\omega\Delta t) + 2}{2\omega^4} \tag{F-Math-15}$$

$$f_3 = -\frac{1-\cos(\omega\Delta t)}{\omega^2}, \quad f_4 = \frac{\omega\Delta t - \sin(\omega\Delta t)}{\omega^3} \tag{F-Math-16}$$

小角度近似（$\omega < 0.5°/\text{s}$）：

$$f_1 \approx -\frac{\Delta t^3}{3}, \quad f_2 \approx \frac{\Delta t^4}{8}, \quad f_3 \approx -\frac{\Delta t^2}{2}, \quad f_4 \approx \frac{\Delta t^3}{6} \tag{F-Math-17}$$

### 2.4 均值更新公式

$$\boldsymbol{\alpha}_{k\to k+1} \mathrel{+}= \boldsymbol{\beta}_{k\to k+1}\Delta t + \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_2\hat{\mathbf{a}} \tag{F-Math-18}$$

$$\boldsymbol{\beta}_{k\to k+1} \mathrel{+}= \mathbf{R}_{k+1\to k}\boldsymbol{\Xi}_1\hat{\mathbf{a}} \tag{F-Math-19}$$

---

## 3. 线性系统构建

### 3.1 IMU 运动方程

在 $I_0$ 系中，利用 CPI 预积分：

$$\mathbf{p}_{I_k}^{I_0} = \mathbf{v}_{I_0}^{I_0}\Delta t_{0k} - \frac{1}{2}\mathbf{g}^{I_0}\Delta t_{0k}^2 + \boldsymbol{\alpha}_{0\to k} \tag{F-Math-20}$$

$$\mathbf{v}_{I_k}^{I_0} = \mathbf{v}_{I_0}^{I_0} - \mathbf{g}^{I_0}\Delta t_{0k} + \boldsymbol{\beta}_{0\to k} \tag{F-Math-21}$$

### 3.2 重投影约束

特征 $j$ 在帧 $k$ 的观测：

$$\begin{bmatrix} 1 & 0 & -u \\ 0 & 1 & -v \end{bmatrix} \mathbf{p}_{F_j}^{C_k} = \mathbf{0} \tag{F-Math-22}$$

坐标变换：

$$\mathbf{p}_{F_j}^{C_k} = \mathbf{R}_{I\to C}\left(\mathbf{R}_{I_0\to I_k}^\top(\mathbf{p}_{F_j}^{I_0} - \mathbf{p}_{I_k}^{I_0})\right) + \mathbf{p}_I^C \tag{F-Math-23}$$

### 3.3 线性方程

将式 (F-Math-20) 代入 (F-Math-23)，再代入 (F-Math-22)：

$$\mathbf{Y}\mathbf{p}_{F_j}^{I_0} - \Delta t_{0k}\mathbf{Y}\mathbf{v}_{I_0}^{I_0} + \frac{1}{2}\Delta t_{0k}^2\mathbf{Y}\mathbf{g}^{I_0} = \mathbf{Y}\boldsymbol{\alpha}_{0\to k} - \mathbf{H}_{\text{proj}}\mathbf{p}_I^C \tag{F-Math-24}$$

其中 $\mathbf{Y} = \mathbf{H}_{\text{proj}}\mathbf{R}_{I\to C}\mathbf{R}_{I_0\to I_k}^\top$。

堆叠所有观测得线性系统：

$$\mathbf{A}\mathbf{x} = \mathbf{b}, \quad \mathbf{x} = [\mathbf{p}_{F_1}^{I_0},\, \ldots,\, \mathbf{p}_{F_M}^{I_0},\, \mathbf{v}_{I_0}^{I_0},\, \mathbf{g}^{I_0}]^\top \tag{F-Math-25}$$

---

## 4. 重力约束求解

### 4.1 拉格朗日乘子法

约束 $\|\mathbf{g}\|^2 = g^2$，构造拉格朗日函数：

$$\mathcal{L}(\mathbf{g}, \lambda) = \|(\mathbf{D} - \lambda\mathbf{I})\mathbf{g} - \mathbf{d}\|^2 + \mu(\|\mathbf{g}\|^2 - g^2) \tag{F-Math-26}$$

最优性条件导出 6 次多项式（见 F2a）：

$$c_6\lambda^6 + c_5\lambda^5 + \cdots + c_1\lambda + c_0 = 0 \tag{F-Math-27}$$

### 4.2 伴随矩阵求根

构造伴随矩阵 $\mathbf{C}$，其特征值即为多项式的根。选择使 $\|\mathbf{g}\|$ 最接近 9.81 的实根。

---

## 5. Ceres MLE 优化

### 5.1 优化变量

$$\mathcal{X} = \{\mathbf{x}_{I_k}\}_{k=0}^N \cup \{\mathbf{p}_{F_j}^G\}_{j=1}^M \tag{F-Math-28}$$

### 5.2 残差因子

**IMU CPI 因子**（15 维，见 F2b）：

$$\mathbf{r}_{\text{IMU}}^{k,k+1} = \begin{bmatrix} 2\,\text{vec}(\mathbf{q}_{\text{res}}) \\ \mathbf{b}_{\omega 2} - \mathbf{b}_{\omega 1} \\ \mathbf{R}_1(\Delta\mathbf{v} - \boldsymbol{\beta}_{\text{corrected}}) \\ \mathbf{b}_{a 2} - \mathbf{b}_{a 1} \\ \mathbf{R}_1(\Delta\mathbf{p} - \boldsymbol{\alpha}_{\text{corrected}}) \end{bmatrix} \tag{F-Math-29}$$

**视觉重投影因子**（2 维，见 F2c）：

$$\mathbf{r}_{\text{cam}}^{j,k} = \mathbf{u}_{\text{dist}} - \mathbf{u}_{\text{meas}} \tag{F-Math-30}$$

### 5.3 目标函数

$$\min_{\mathcal{X}} \sum_{k} \|\mathbf{r}_{\text{IMU}}^{k,k+1}\|_{\mathbf{R}_{\text{IMU}}^{-1}}^2 + \sum_{j,k} \|\mathbf{r}_{\text{cam}}^{j,k}\|_{\sigma_{\text{pix}}^{-2}}^2 \tag{F-Math-31}$$

---

## 6. 完整算法流程

```
输入: IMU 数据, 特征观测
        │
        ├─ 1. 数据准备
        │     ├─ 窗口裁剪
        │     └─ 特征深拷贝
        │
        ├─ 2. 特征验证
        │     ├─ 测量数 ≥ min_num_meas
        │     ├─ 有效特征 ≥ min_valid_features
        │     ─ 角运动 ≥ init_dyn_min_deg
        │
        ├─ 3. CPI 预积分
        │     ├─ I0toIi (线性系统)
        │     └─ IitoIi1 (MLE)
        │
        ├─ 4. 线性系统求解
        │     ├─ 构造 A x = b
        │     ├─ compute_dongsi_coeff (重力约束)
        │     ├─ 伴随矩阵特征值分解
        │     └─ 回代：g, v, 特征位置
        │
        ├─ 5. 重力对齐
        │     ─ Gram-Schmidt → 全局坐标系
        │
        ├─ 6. Ceres MLE
        │     ├─ 添加变量 (位姿 + 速度 + bias + 特征)
        │     ├─ 添加因子 (IMU + 视觉 + 先验)
        │     └─ 求解 (DENSE_SCHUR + DOGLEG)
        │
        ├─ 7. 协方差恢复
        │     ├─ ceres::Covariance
        │     ├─ 膨胀 (10×/100×)
        │     └─ 对称化
        │
        └─ 8. 状态输出
              ├─ _imu (IMU 状态)
              ├─ _clones_IMU (位姿)
              └─ _features_SLAM (特征)
```

---

## 7. 公式索引

| 编号 | 公式 | 位置 |
|------|------|------|
| F-Math-1 | IMU 状态定义 | §1.1 |
| F-Math-2 | 特征位置定义 | §1.1 |
| F-Math-3 | 重力约束 | §1.1 |
| F-Math-4 | IMU 预积分测量 | §1.2 |
| F-Math-5 | 视觉重投影 | §1.2 |
| F-Math-6 | 旋转运动学 | §2.1 |
| F-Math-7 | 位置运动学 | §2.1 |
| F-Math-8 | 速度运动学 | §2.1 |
| F-Math-9 | 相对旋转 | §2.2 |
| F-Math-10 | 位置增量 | §2.2 |
| F-Math-11 | 速度增量 | §2.2 |
| F-Math-12 | Rodrigues 展开 | §2.3 |
| F-Math-13 | ₁ 公式 | §2.3 |
| F-Math-14 | Ξ₂ 公式 | §2.3 |
| F-Math-15 | 系数 f₁,f₂ | §2.3 |
| F-Math-16 | 系数 f₃,f₄ | §2.3 |
| F-Math-17 | 小角度近似 | §2.3 |
| F-Math-18 | α 更新 | §2.4 |
| F-Math-19 | β 更新 | §2.4 |
| F-Math-20 | 位置方程 | §3.1 |
| F-Math-21 | 速度方程 | §3.1 |
| F-Math-22 | 重投影约束 | §3.2 |
| F-Math-23 | 坐标变换 | §3.2 |
| F-Math-24 | 线性方程 | §3.3 |
| F-Math-25 | 线性系统 | §3.3 |
| F-Math-26 | 拉格朗日函数 | §4.1 |
| F-Math-27 | 6 次多项式 | §4.1 |
| F-Math-28 | 优化变量 | §5.1 |
| F-Math-29 | IMU 残差 | §5.2 |
| F-Math-30 | 视觉残差 | §5.2 |
| F-Math-31 | 目标函数 | §5.3 |

---

## 8. 子文档索引

| 文档 | 内容 |
|------|------|
| **F1a** | `F1a_StaticInitializer-algorithm-deep-dive.md` - 静态初始化算法 |
| **F1b** | `F1b_DynamicInitializer-algorithm-deep-dive.md` - 动态初始化算法 |
| **F2a** | `F2a_compute_dongsi_coeff-deep-dive.md` - 重力约束求解 |
| **F2b** | `F2b_Factor_ImuCPIv1-deep-dive.md` - IMU CPI 因子 |
| **F2c** | `F2c_Factor_ImageReprojCalib-deep-dive.md` - 视觉重投影因子 |
| **F2d** | `F2d_InitializerHelper-deep-dive.md` - 辅助函数 (gram_schmidt, select_imu_readings) |
| **F3a** | `F3a_CpiV1_analytic_coefficients-deep-dive.md` - CPI 解析系数 |
| **F3b** | `F3b_CpiV1_bias_jacobians-deep-dive.md` - CPI 偏置 Jacobian |

---

## 9. 代码位置索引

| 功能 | 文件 | 行号 |
|------|------|------|
| 静态初始化 | `ov_init/src/static/StaticInitializer.cpp` | 37-165 |
| 动态初始化 | `ov_init/src/dynamic/DynamicInitializer.cpp` | 44-1107 |
| CPI 预积分 | `ov_core/src/cpi/CpiV1.cpp` | 33-337 |
| 解析系数 | `ov_core/src/cpi/CpiV1.cpp` | 100-127 |
| 偏置 Jacobian | `ov_core/src/cpi/CpiV1.cpp` | 128-238 |
| 协方差 RK4 | `ov_core/src/cpi/CpiV1.cpp` | 240-331 |
| 重力约束 | `ov_init/src/utils/helper.h` | 183-319 |
| IMU 因子 | `ov_init/src/ceres/Factor_ImuCPIv1.cpp` | 28-281 |
| 视觉因子 | `ov_init/src/ceres/Factor_ImageReprojCalib.cpp` | 28-152 |

---

## 10. 待补充子文档

### F3b: CpiV1 偏置 Jacobian 推导

**内容**：
- $\mathbf{J}_q, \mathbf{J}_\alpha, \mathbf{J}_\beta$ 对 $\mathbf{b}_\omega$ 的 Jacobian 推导
- $\mathbf{H}_\alpha, \mathbf{H}_\beta$ 对 $\mathbf{b}_a$ 的 Jacobian 推导
- 小角度近似下的简化公式
- 代码对照（CpiV1.cpp:128-238）

### F2d: InitializerHelper 辅助函数

**内容**：
- `gram_schmidt()` 算法推导
- `select_imu_readings()` 插值策略
- `interpolate_data()` 线性插值

> 以上子文档待后续补充。
