# Factor_ImageReprojCalib 重投影因子深度解析 (F2c)

> **生成日期**：2026-08-25
> **前置文档**：`F1b_DynamicInitializer-algorithm-deep-dive.md`（动态初始化算法）
> **核心代码**：`ov_init/src/ceres/Factor_ImageReprojCalib.cpp`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Hartley & Zisserman "Multiple View Geometry" 2004; OpenVINS CPI 论文

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
Ceres MLE 优化 (F1b)
        │  添加因子
        ▼
Factor_ImageReprojCalib  ★ 本篇
   ├─ 输入: 观测像素 (u,v)、相机内参、位姿、特征位置
   ├─ Evaluate(): 计算重投影残差 + Jacobian
   └─ 输出: 2 维残差向量 + 6 个参数块的 Jacobian
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | Ceres 图优化中的视觉重投影因子 |
| **输入** | 观测像素、相机内参（8 维）、位姿（q, p）、特征位置、外参（q, p） |
| **输出** | 2 维残差、6×2 Jacobian 矩阵 |
| **调用方** | `DynamicInitializer::initialize()`（F1b） |
| **被调用方** | Ceres 优化器 |

### 1.2 一句话概括

`Factor_ImageReprojCalib` 计算特征重投影残差（观测像素 - 投影像素）及其对所有优化变量（位姿、特征、外参、内参）的 Jacobian，支持 radtan 和 fisheye 两种相机模型。

---

## 2. 数学模型

### 2.1 坐标变换链

特征从全局坐标系到相机坐标系的变换：

$$\mathbf{p}_F^{I_i} = \mathbf{R}_{G\to I_i}(\mathbf{p}_F^G - \mathbf{p}_{I_i}^G) \tag{F2c-1}$$

$$\mathbf{p}_F^{C_i} = \mathbf{R}_{I\to C}\mathbf{p}_F^{I_i} + \mathbf{p}_I^C \tag{F2c-2}$$

### 2.2 投影模型

归一化坐标：

$$\mathbf{u}_{\text{norm}} = \begin{bmatrix} x/z \\ y/z \end{bmatrix}, \quad \mathbf{p}_F^{C_i} = [x, y, z]^\top \tag{F2c-3}$$

畸变模型（radtan 或 fisheye）：

$$\mathbf{u}_{\text{dist}} = \text{distort}(\mathbf{u}_{\text{norm}}, \boldsymbol{\zeta}) \tag{F2c-4}$$

其中 $\boldsymbol{\zeta} = [f_x, f_y, c_x, c_y, k_1, k_2, k_3, k_4]^\top$ 是相机内参。

### 2.3 残差定义

$$\mathbf{r} = \mathbf{u}_{\text{dist}} - \mathbf{u}_{\text{meas}} \tag{F2c-5}$$

加权残差（信息矩阵平方根）：

$$\mathbf{r}_{\text{weighted}} = \mathbf{L}^\top \mathbf{r}, \quad \mathbf{L}\mathbf{L}^\top = \sigma_{\text{pix}}^{-2}\mathbf{I} \tag{F2c-6}$$

### 2.4 Jacobian 计算

对位姿 $\mathbf{q}_{G\to I_i}$ 的 Jacobian：

$$\frac{\partial\mathbf{r}}{\partial\delta\boldsymbol{\theta}} = \frac{\partial\mathbf{u}_{\text{dist}}}{\partial\mathbf{u}_{\text{norm}}} \cdot \frac{\partial\mathbf{u}_{\text{norm}}}{\partial\mathbf{p}_F^{C_i}} \cdot \mathbf{R}_{I\to C} \cdot [\mathbf{p}_F^{I_i}]_\times \tag{F2c-7}$$

对特征位置 $\mathbf{p}_F^G$ 的 Jacobian：

$$\frac{\partial\mathbf{r}}{\partial\mathbf{p}_F^G} = \frac{\partial\mathbf{u}_{\text{dist}}}{\partial\mathbf{u}_{\text{norm}}} \cdot \frac{\partial\mathbf{u}_{\text{norm}}}{\partial\mathbf{p}_F^{C_i}} \cdot \mathbf{R}_{I\to C}\mathbf{R}_{G\to I_i} \tag{F2c-8}$$

对外参 $\mathbf{q}_{I\to C}$ 的 Jacobian：

$$\frac{\partial\mathbf{r}}{\partial\delta\boldsymbol{\theta}_{I\to C}} = \frac{\partial\mathbf{u}_{\text{dist}}}{\partial\mathbf{u}_{\text{norm}}} \cdot \frac{\partial\mathbf{u}_{\text{norm}}}{\partial\mathbf{p}_F^{C_i}} \cdot [\mathbf{R}_{I\to C}\mathbf{p}_F^{I_i}]_\times \tag{F2c-9}$$

---

## 3. 代码实现：逐行对照

### 3.1 构造函数

```cpp
// Factor_ImageReprojCalib.cpp:28-44
Factor_ImageReprojCalib::Factor_ImageReprojCalib(
    const Eigen::Vector2d &uv_meas_, double pix_sigma_, bool is_fisheye_)
    : uv_meas(uv_meas_), pix_sigma(pix_sigma_), is_fisheye(is_fisheye_) {
    
    // 信息矩阵平方根 (式 F2c-6)
    sqrtQ = Eigen::Matrix<double, 2, 2>::Identity();
    sqrtQ(0, 0) *= 1.0 / pix_sigma;
    sqrtQ(1, 1) *= 1.0 / pix_sigma;
    
    // 设置残差维度和参数块
    set_num_residuals(2);
    mutable_parameter_block_sizes()->push_back(4); // q_GtoIi
    mutable_parameter_block_sizes()->push_back(3); // p_IiinG
    mutable_parameter_block_sizes()->push_back(3); // p_FinG
    mutable_parameter_block_sizes()->push_back(4); // q_ItoC
    mutable_parameter_block_sizes()->push_back(3); // p_IinC
    mutable_parameter_block_sizes()->push_back(8); // 内参
}
```

### 3.2 Evaluate()：残差计算

```cpp
// Factor_ImageReprojCalib.cpp:46-103
bool Factor_ImageReprojCalib::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    
    // 获取状态
    Eigen::Vector4d q_GtoIi = Eigen::Map<const Eigen::Vector4d>(parameters[0]);
    Eigen::Vector3d p_IiinG = Eigen::Map<const Eigen::Vector3d>(parameters[1]);
    Eigen::Vector3d p_FinG = Eigen::Map<const Eigen::Vector3d>(parameters[2]);
    Eigen::Vector4d q_ItoC = Eigen::Map<const Eigen::Vector4d>(parameters[3]);
    Eigen::Vector3d p_IinC = Eigen::Map<const Eigen::Vector3d>(parameters[4]);
    Eigen::Matrix<double, 8, 1> camera_vals = Eigen::Map<const Eigen::Matrix<double, 8, 1>>(parameters[5]);
    
    // 坐标变换 (式 F2c-1, F2c-2)
    Eigen::Vector3d p_FinIi = R_GtoIi * (p_FinG - p_IiinG);
    Eigen::Vector3d p_FinCi = R_ItoC * p_FinIi + p_IinC;
    
    // 归一化投影 (式 F2c-3)
    Eigen::Vector2d uv_norm;
    uv_norm << p_FinCi(0) / p_FinCi(2), p_FinCi(1) / p_FinCi(2);
    
    // 畸变 (式 F2c-4)
    Eigen::Vector2d uv_dist;
    if (is_fisheye) {
        ov_core::CamEqui cam(0, 0);
        cam.set_value(camera_vals);
        uv_dist = cam.distort_d(uv_norm);
    } else {
        ov_core::CamRadtan cam(0, 0);
        cam.set_value(camera_vals);
        uv_dist = cam.distort_d(uv_norm);
    }
    
    // 残差 (式 F2c-5)
    Eigen::Vector2d res = uv_dist - uv_meas;
    res = sqrtQ_gate * res;  // 信息矩阵加权 (式 F2c-6)
}
```

### 3.3 Jacobian 计算

```cpp
// Factor_ImageReprojCalib.cpp:106-150
if (jacobians) {
    // 投影 Jacobian (式 F2c-3 的导数)
    Eigen::MatrixXd H_dzn_dpfc = Eigen::MatrixXd::Zero(2, 3);
    H_dzn_dpfc << 1.0/p_FinCi(2), 0, -p_FinCi(0)/p_FinCi(2)^2,
                  0, 1.0/p_FinCi(2), -p_FinCi(1)/p_FinCi(2)^2;
    
    Eigen::MatrixXd H_dz_dpfc = H_dz_dzn * H_dzn_dpfc;
    
    // 对 q_GtoIi 的 Jacobian (式 F2c-7)
    if (jacobians[0]) {
        jacobian.block(0, 0, 2, 3) = H_dz_dpfc * R_ItoC * skew_x(p_FinIi);
        jacobian.block(0, 3, 2, 1).setZero();
    }
    
    // 对 p_IiinG 的 Jacobian
    if (jacobians[1]) {
        jacobian.block(0, 0, 2, 3) = -H_dz_dpfc * R_ItoC * R_GtoIi;
    }
    
    // 对 p_FinG 的 Jacobian (式 F2c-8)
    if (jacobians[2]) {
        jacobian.block(0, 0, 2, 3) = H_dz_dpfc * R_ItoC * R_GtoIi;
    }
    
    // 对 q_ItoC 的 Jacobian (式 F2c-9)
    if (jacobians[3]) {
        jacobian.block(0, 0, 2, 3) = H_dz_dpfc * skew_x(R_ItoC * p_FinIi);
        jacobian.block(0, 3, 2, 1).setZero();
    }
    
    // 对 p_IinC 的 Jacobian
    if (jacobians[4]) {
        jacobian.block(0, 0, 2, 3) = H_dz_dpfc;
    }
    
    // 对内参的 Jacobian
    if (jacobians[5]) {
        jacobian.block(0, 0, 2, 8) = H_dz_dzeta;
    }
}
```

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{u}_{\text{meas}}$ | 观测像素 | `uv_meas` | `:29` |
| $\sigma_{\text{pix}}$ | 像素噪声标准差 | `pix_sigma` | `:29` |
| $\mathbf{p}_F^G$ | 特征全局位置 | `p_FinG` | `:52` |
| $\mathbf{p}_F^{C_i}$ | 特征相机坐标 | `p_FinCi` | `:62` |
| $\mathbf{u}_{\text{norm}}$ | 归一化坐标 | `uv_norm` | `:65-66` |
| $\mathbf{u}_{\text{dist}}$ | 畸变后坐标 | `uv_dist` | `:73-93` |
| $\mathbf{r}$ | 残差 | `res` | `:100-101` |

---

## 5. 关键设计决策

### 5.1 为什么残差是 $\mathbf{u}_{\text{dist}} - \mathbf{u}_{\text{meas}}$

Ceres 默认最小化 $\|f(\mathbf{x})\|^2$。重投影约束是 $\mathbf{u}_{\text{meas}} = f(\mathbf{x}) + \mathbf{n}$，改写为 $0 = f(\mathbf{x}) + \mathbf{n} - \mathbf{u}_{\text{meas}}$，故残差为 $f(\mathbf{x}) - \mathbf{u}_{\text{meas}}$。

### 5.2 为什么支持两种相机模型

| 模型 | 适用场景 | 畸变参数 |
|------|---------|---------|
| **Radtan** | 普通镜头 | $k_1, k_2, k_3, k_4$（径向 + 切向） |
| **Fisheye** | 鱼眼镜头 | $k_1, k_2, k_3, k_4$（等距模型） |

### 5.3 为什么用信息矩阵平方根加权

$$\mathbf{r}_{\text{weighted}} = \mathbf{L}^\top\mathbf{r}, \quad \mathbf{L}\mathbf{L}^\top = \sigma_{\text{pix}}^{-2}\mathbf{I} \tag{F2c-10}$$

Ceres 最小化 $\|\mathbf{r}\|^2$。加权后等价于最小化 $\mathbf{r}^\top\mathbf{R}^{-1}\mathbf{r}$（马氏距离），考虑测量不确定性。

---

## 6. 完整流程图

```
输入: 观测像素 (u,v)、相机内参、位姿、特征位置、外参
        │
        ├─ 1. 坐标变换
        │     ├─ p_FinIi = R_GtoIi * (p_FinG - p_IiinG)
        │     └─ p_FinCi = R_ItoC * p_FinIi + p_IinC
        │
        ├─ 2. 投影
        │     └─ uv_norm = [x/z, y/z]
        │
        ├─ 3. 畸变
        │     ─ uv_dist = distort(uv_norm, 内参)
        │
        ├─ 4. 残差
        │     ├─ res = uv_dist - uv_meas
        │     └─ res_weighted = L^T * res
        │
        └─ 5. Jacobian (如果请求)
              ├─ 对 q_GtoIi: H_dz_dpfc * R_ItoC * skew(p_FinIi)
              ├─ 对 p_IiinG: -H_dz_dpfc * R_ItoC * R_GtoIi
              ├─ 对 p_FinG: H_dz_dpfc * R_ItoC * R_GtoIi
              ├─ 对 q_ItoC: H_dz_dpfc * skew(R_ItoC * p_FinIi)
              ├─ 对 p_IinC: H_dz_dpfc
              └─ 对内参: H_dz_dzeta
```

---

## 7. 与其他文档的关联

- **上游**：`F1b_DynamicInitializer-algorithm-deep-dive.md`（调用本因子）
- **下游**：Ceres 优化器
- **相关**：Hartley & Zisserman "Multiple View Geometry"（相机模型理论）

### 待详细补充项

- **F2c-1**：Radtan 畸变模型的 Jacobian 推导
- **F2c-2**：Fisheye 畸变模型的 Jacobian 推导

> 以上子文档暂不展开，待需要逐项深挖时再补。
