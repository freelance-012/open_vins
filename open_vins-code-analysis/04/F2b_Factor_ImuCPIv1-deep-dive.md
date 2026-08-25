# Factor_ImuCPIv1 残差构造深度解析 (F2b)

> **生成日期**：2026-08-25
> **前置文档**：`F1b_DynamicInitializer-algorithm-deep-dive.md`（动态初始化算法）、`F3_CpiV1-continuous-preintegration-deep-dive.md`（CPI 预积分）
> **核心代码**：`ov_init/src/ceres/Factor_ImuCPIv1.cpp`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Dong-Si & Mourikis 2012 Eq.(91)-(107); Eckenhoff et al. IJRR 2019

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
Ceres MLE 优化 (F1b)
        │  添加因子
        ▼
Factor_ImuCPIv1  ★ 本篇
   ├─ 输入: CPI 测量 (α, β, ΔR, Δt) + 偏置线性化点 + Jacobian + 协方差
   ├─ Evaluate(): 计算残差 + Jacobian
   └─ 输出: 15 维残差向量 + 10 个参数块的 Jacobian
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | Ceres 图优化中的 IMU 预积分因子 |
| **输入** | CPI 测量、偏置线性化点、bias Jacobian、协方差 |
| **输出** | 15 维残差、10×15 Jacobian 矩阵 |
| **调用方** | `DynamicInitializer::initialize()`（F1b） |
| **被调用方** | Ceres 优化器 |

### 1.2 一句话概括

`Factor_ImuCPIv1` 计算 CPI 预积分测量与当前状态预测之间的残差（15 维），并计算对所有优化变量（两帧位姿+速度+bias+位置）的 Jacobian，用于 Ceres 图优化。

---

## 2. 数学模型

### 2.1 残差定义

CPI 预积分测量（见 F3）：
- $\boldsymbol{\alpha}$：位置增量
- $\boldsymbol{\beta}$：速度增量
- $\Delta\mathbf{R}$：相对旋转
- $\Delta t$：时间间隔
- $\mathbf{J}_q, \mathbf{J}_\alpha, \mathbf{J}_\beta, \mathbf{H}_\alpha, \mathbf{H}_\beta$：bias Jacobian

状态预测（从帧 1 到帧 2）：

$$\mathbf{R}_{1\to 2} = \mathbf{R}_2\mathbf{R}_1^\top \tag{F2b-1}$$

$$\Delta\mathbf{v} = \mathbf{v}_2 - \mathbf{v}_1 + \mathbf{g}\Delta t \tag{F2b-2}$$

$$\Delta\mathbf{p} = \mathbf{p}_2 - \mathbf{p}_1 - \mathbf{v}_1\Delta t + \frac{1}{2}\mathbf{g}\Delta t^2 \tag{F2b-3}$$

### 2.2 偏置重线性化

由于 CPI 是在线性化点 $(\mathbf{b}_\omega^{\text{lin}}, \mathbf{b}_a^{\text{lin}})$ 处计算的，当 bias 变化 $\delta\mathbf{b}_\omega = \mathbf{b}_\omega - \mathbf{b}_\omega^{\text{lin}}$ 时，用一阶近似修正：

$$\Delta\mathbf{R}_{\text{corrected}} \approx \Delta\mathbf{R} \otimes \delta\mathbf{q}_b, \quad \delta\mathbf{q}_b \approx \begin{bmatrix} \frac{1}{2}\mathbf{J}_q\delta\mathbf{b}_\omega \\ 1 \end{bmatrix} \tag{F2b-4}$$

$$\boldsymbol{\alpha}_{\text{corrected}} \approx \boldsymbol{\alpha} + \mathbf{J}_\alpha\delta\mathbf{b}_\omega + \mathbf{H}_\alpha\delta\mathbf{b}_a \tag{F2b-5}$$

$$\boldsymbol{\beta}_{\text{corrected}} \approx \boldsymbol{\beta} + \mathbf{J}_\beta\delta\mathbf{b}_\omega + \mathbf{H}_\beta\delta\mathbf{b}_a \tag{F2b-6}$$

### 2.3 残差向量（15 维）

$$\mathbf{r} = \begin{bmatrix} \mathbf{r}_q \\ \mathbf{r}_{b_\omega} \\ \mathbf{r}_v \\ \mathbf{r}_{b_a} \\ \mathbf{r}_p \end{bmatrix} = \begin{bmatrix} 2\,\text{vec}(\mathbf{q}_{\text{res}}) \\ \mathbf{b}_{\omega 2} - \mathbf{b}_{\omega 1} \\ \mathbf{R}_1(\Delta\mathbf{v} - \boldsymbol{\beta}_{\text{corrected}}) \\ \mathbf{b}_{a 2} - \mathbf{b}_{a 1} \\ \mathbf{R}_1(\Delta\mathbf{p} - \boldsymbol{\alpha}_{\text{corrected}}) \end{bmatrix} \tag{F2b-7}$$

其中：
- $\mathbf{q}_{\text{res}} = \mathbf{R}_{1\to 2} \otimes \Delta\mathbf{R}^{-1} \otimes \delta\mathbf{q}_b$（旋转残差四元数的向量部分）
- $\text{vec}(\mathbf{q})$ 取四元数的前 3 个分量（向量部分）

### 2.4 信息矩阵加权

残差乘以信息矩阵的平方根：

$$\mathbf{r}_{\text{weighted}} = \mathbf{L}^\top \mathbf{r}, \quad \mathbf{L}\mathbf{L}^\top = \mathbf{P}^{-1} \tag{F2b-8}$$

其中 $\mathbf{P}$ 是 CPI 测量的协方差矩阵（15×15）。

---

## 3. 代码实现：逐行对照

### 3.1 构造函数：保存测量

```cpp
// Factor_ImuCPIv1.cpp:28-75
Factor_ImuCPIv1::Factor_ImuCPIv1(double deltatime, Eigen::Vector3d &grav,
    Eigen::Vector3d &alpha, Eigen::Vector3d &beta, Eigen::Vector4d &q_KtoK1,
    Eigen::Vector3d &ba_lin, Eigen::Vector3d &bg_lin,
    Eigen::Matrix3d &J_q, Eigen::Matrix3d &J_beta, Eigen::Matrix3d &J_alpha,
    Eigen::Matrix3d &H_beta, Eigen::Matrix3d &H_alpha,
    Eigen::Matrix<double, 15, 15> &covariance) {
    
    // 保存 CPI 测量
    this->alpha = alpha;
    this->beta = beta;
    this->q_breve = q_KtoK1;
    this->dt = deltatime;
    this->grav_save = grav;
    
    // 保存线性化点
    this->b_a_lin_save = ba_lin;
    this->b_w_lin_save = bg_lin;
    
    // 保存 bias Jacobian
    this->J_q = J_q;
    this->J_a = J_alpha;
    this->J_b = J_beta;
    this->H_a = H_alpha;
    this->H_b = H_beta;
    
    // 计算信息矩阵的平方根
    Eigen::MatrixXd information = covariance.llt().solve(I);
    Eigen::LLT<Eigen::MatrixXd> lltOfI(information);
    sqrtI_save = lltOfI.matrixL().transpose();
    
    // 设置残差维度和参数块
    set_num_residuals(15);
    mutable_parameter_block_sizes()->push_back(4); // q_GtoI1
    mutable_parameter_block_sizes()->push_back(3); // bg_1
    // ... 共 10 个参数块
}
```

**理论对应**：式 (F2b-8) 的信息矩阵预处理。

### 3.2 Evaluate()：残差计算

```cpp
// Factor_ImuCPIv1.cpp:77-141
bool Factor_ImuCPIv1::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    
    // 获取状态估计
    Eigen::Vector4d q_1 = Eigen::Map<const Eigen::Vector4d>(parameters[0]);
    Eigen::Vector4d q_2 = Eigen::Map<const Eigen::Vector4d>(parameters[5]);
    Eigen::Vector3d b_w1 = Eigen::Map<const Eigen::Vector3d>(parameters[1]);
    Eigen::Vector3d b_w2 = Eigen::Map<const Eigen::Vector3d>(parameters[6]);
    Eigen::Vector3d b_a1 = Eigen::Map<const Eigen::Vector3d>(parameters[3]);
    Eigen::Vector3d b_a2 = Eigen::Map<const Eigen::Vector3d>(parameters[8]);
    Eigen::Vector3d v_1 = Eigen::Map<const Eigen::Vector3d>(parameters[2]);
    Eigen::Vector3d v_2 = Eigen::Map<const Eigen::Vector3d>(parameters[7]);
    Eigen::Vector3d p_1 = Eigen::Map<const Eigen::Vector3d>(parameters[4]);
    Eigen::Vector3d p_2 = Eigen::Map<const Eigen::Vector3d>(parameters[9]);
    
    // 偏置变化
    Eigen::Vector3d dbw = b_w1 - b_w_lin;
    Eigen::Vector3d dba = b_a1 - b_a_lin;
    
    // 偏置修正四元数 (式 F2b-4)
    Eigen::Vector4d q_b;
    q_b.block(0, 0, 3, 1) = 0.5 * J_q * dbw;
    q_b(3, 0) = 1.0;
    q_b = q_b / q_b.norm();
    
    // 相对旋转
    Eigen::Vector4d q_1_to_2 = ov_core::quat_multiply(q_2, ov_core::Inv(q_1));
    
    // 旋转残差
    Eigen::Vector4d q_res_minus = ov_core::quat_multiply(q_1_to_2, ov_core::Inv(q_breve));
    Eigen::Vector4d q_res_plus = ov_core::quat_multiply(q_res_minus, q_b);
    
    // 15 维残差 (式 F2b-7)
    Eigen::Matrix<double, 15, 1> res;
    res.block(0, 0, 3, 1) = 2 * q_res_plus.block(0, 0, 3, 1);          // 旋转
    res.block(3, 0, 3, 1) = b_w2 - b_w1;                                // 陀螺零偏
    res.block(6, 0, 3, 1) = R_1 * (v_2 - v_1 + gravity * dt) - J_b * dbw - H_b * dba - beta;  // 速度
    res.block(9, 0, 3, 1) = b_a2 - b_a1;                                // 加表零偏
    res.block(12, 0, 3, 1) = R_1 * (p_2 - p_1 - v_1 * dt + .5 * gravity * dt^2) - J_a * dbw - H_a * dba - alpha;  // 位置
    
    // 信息矩阵加权 (式 F2b-8)
    res = sqrtI * res;
}
```

### 3.3 Jacobian 计算

```cpp
// Factor_ImuCPIv1.cpp:147-278
if (jacobians) {
    Eigen::Matrix<double, 15, 30> Jacobian = Eigen::Matrix<double, 15, 30>::Zero();
    
    // 旋转对旋转 1 的 Jacobian (Dong-Si Eq.92)
    Jacobian.block(0, 0, 3, 3) = -((q_1_to_2(3) * I - skew(q_1_to_2.vec())) *
                                    (q_meas_plus(3) * I + skew(q_meas_plus.vec())) -
                                    q_1_to_2.vec() * q_meas_plus.vec().transpose());
    
    // 旋转对旋转 2 的 Jacobian (Dong-Si Eq.92)
    Jacobian.block(0, 15, 3, 3) = q_res_plus(3) * I + skew(q_res_plus.vec());
    
    // 速度对速度 1 的 Jacobian (Dong-Si Eq.97)
    Jacobian.block(6, 6, 3, 3) = -R_1;
    
    // 速度对速度 2 的 Jacobian (Dong-Si Eq.98)
    Jacobian.block(6, 21, 3, 3) = R_1;
    
    // 位置对位置 1 的 Jacobian (Dong-Si Eq.106)
    Jacobian.block(12, 12, 3, 3) = -R_1;
    
    // 位置对位置 2 的 Jacobian (Dong-Si Eq.107)
    Jacobian.block(12, 27, 3, 3) = R_1;
    
    // 应用信息矩阵加权
    Jacobian = sqrtI * Jacobian;
    
    // 存储到 Ceres 要求的格式
    if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 15, 4, Eigen::RowMajor>> J_th1(jacobians[0], 15, 4);
        J_th1.block(0, 0, 15, 3) = Jacobian.block(0, 0, 15, 3);
        J_th1.block(0, 3, 15, 1).setZero();  // 四元数第 4 分量无梯度
    }
    // ... 其他参数块类似
}
```

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\boldsymbol{\alpha}$ | CPI 位置增量 | `alpha` | `:33` |
| $\boldsymbol{\beta}$ | CPI 速度增量 | `beta` | `:34` |
| $\Delta\mathbf{R}$ | CPI 相对旋转 | `q_breve` | `:35` |
| $\mathbf{J}_q$ | 旋转对 $b_\omega$ 的 Jacobian | `J_q` | `:44` |
| $\mathbf{J}_\alpha, \mathbf{J}_\beta$ | 增量对 $b_\omega$ 的 Jacobian | `J_a`, `J_b` | `:45-46` |
| $\mathbf{H}_\alpha, \mathbf{H}_\beta$ | 增量对 $b_a$ 的 Jacobian | `H_a`, `H_b` | `:47-48` |
| $\mathbf{r}$ | 15 维残差 | `res` | `:130-136` |
| $\mathbf{J}$ | 15×30 Jacobian | `Jacobian` | `:151` |

---

## 5. 关键设计决策

### 5.1 为什么用一阶偏置修正

CPI 预积分是在固定线性化点计算的。当 bias 变化时，重新积分代价高。一阶近似：

$$\Delta\mathbf{R}(\mathbf{b}_\omega) \approx \Delta\mathbf{R}(\mathbf{b}_\omega^{\text{lin}}) \otimes \delta\mathbf{q}_b(\delta\mathbf{b}_\omega) \tag{F2b-9}$$

精度足够（bias 变化小），计算快（只需矩阵乘法）。

### 5.2 为什么残差是 15 维

| 分量 | 维度 | 物理意义 |
|------|------|---------|
| 旋转 | 3 | 姿态误差（四元数向量部分） |
| 陀螺零偏 | 3 | bias 随机游走 |
| 速度 | 3 | 速度增量误差 |
| 加表零偏 | 3 | bias 随机游走 |
| 位置 | 3 | 位置增量误差 |
| **总计** | **15** | — |

### 5.3 为什么用信息矩阵平方根加权

$$\mathbf{r}_{\text{weighted}} = \mathbf{L}^\top \mathbf{r}, \quad \mathbf{L}\mathbf{L}^\top = \mathbf{R}^{-1} \tag{F2b-10}$$

Ceres 默认最小化 $\|\mathbf{r}\|^2$。加权后等价于最小化 $\mathbf{r}^\top\mathbf{R}^{-1}\mathbf{r}$（马氏距离），考虑测量不确定性。

---

## 6. 完整流程图

```
输入: CPI 测量 (α, β, ΔR, Δt) + 偏置线性化点 + Jacobian + 协方差
        │
        ├─ 1. 获取状态估计
        │     └─ q1, q2, v1, v2, p1, p2, bw1, bw2, ba1, ba2
        │
        ├─ 2. 偏置修正
        │     ├─ δbw = bw1 - bw_lin
        │     ├─ δba = ba1 - ba_lin
        │     └─ δqb ≈ [0.5*J_q*δbw; 1] (一阶近似)
        │
        ├─ 3. 计算残差 (15 维)
        │     ├─ rq = 2*vec(q_res) (旋转)
        │     ├─ rbw = bw2 - bw1 (陀螺零偏)
        │     ├─ rv = R1*(Δv - β_corrected) (速度)
        │     ├─ rba = ba2 - ba1 (加表零偏)
        │     └─ rp = R1*(Δp - α_corrected) (位置)
        │
        ├─ 4. 信息矩阵加权
        │     └─ r_weighted = L^T * r
        │
        └─ 5. 计算 Jacobian (如果请求)
              ├─ 对 10 个参数块求导
              └─ 应用信息矩阵加权
```

---

## 7. 与其他文档的关联

- **上游**：`F3_CpiV1-continuous-preintegration-deep-dive.md`（提供 CPI 测量）
- **下游**：Ceres 优化器（F1b）
- **相关**：Dong-Si & Mourikis 2012 Eq.(91)-(107)（Jacobian 公式来源）

### 待详细补充项

- **F2b-1**：旋转残差 Jacobian 的详细推导（Dong-Si Eq.92）
- **F2b-2**：信息矩阵平方根的计算（LLT 分解）

> 以上子文档暂不展开，待需要逐项深挖时再补。
