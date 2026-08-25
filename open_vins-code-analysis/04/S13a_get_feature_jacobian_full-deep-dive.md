# UpdaterHelper::get_feature_jacobian_full 深度精读 (S13a)

> **生成日期**：2026-08-19
> **前置文档**：`S13_UpdaterMSCKF-deep-dive.md`（主文档）、`S2_JPLQuat-deep-dive.md`（JPL 四元数）
> **核心代码**：`ov_msckf/src/update/UpdaterHelper.cpp:32-424`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
UpdaterMSCKF::update()
        │  （每个特征循环）
        ▼
UpdaterHelper::get_feature_jacobian_full(state, feat, H_f, H_x, res, x_order)  @ :192  ★ 本篇
        │
        ├─ get_feature_jacobian_representation()  @ :32  （特征表示的 ∂p_F^G/∂λ）
        │
        └─ 返回: H_f（特征 Jacobian）、H_x（状态 Jacobian）、res（残差）、x_order（变量顺序）
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | MSCKF/SLAM 更新中每个特征的 Jacobian 计算 |
| 输入 | `state`（当前状态）、`feat`（特征观测数据） |
| 输出 | `H_f`（$2N \times p$）、`H_x`（$2N \times m$）、`res`（$2N \times 1$）、`x_order`（变量列表） |
| 调用方 | `UpdaterMSCKF::update()`、`UpdaterSLAM::update()`、`UpdaterSLAM::delayed_init()` |

### 1.2 一句话概括

`get_feature_jacobian_full` 对单个特征的所有观测，计算重投影残差、特征 Jacobian $H_f$、状态 Jacobian $H_x$，并记录每个 Jacobian 列对应的状态变量（`x_order`）。

---

## 2. 数学模型

### 2.1 重投影残差

特征 $f$ 在 $N$ 个观测时刻被观测到，每个观测产生 2 维残差：

$$\mathbf{r}_c = \mathbf{u}_m - \pi(\mathbf{p}_F^{C_c}) \in \mathbb{R}^2 \tag{S13a-1}$$

其中 $\mathbf{u}_m$ 是观测像素坐标，$\pi(\cdot)$ 是相机投影模型（含畸变），$\mathbf{p}_F^{C_c}$ 是特征在相机 $c$ 坐标系中的位置。

### 2.2 坐标变换链

$$\mathbf{p}_F^{C_c} = \mathbf{R}_{I\to C}\,\mathbf{p}_F^{I_c} + \mathbf{p}_I^C \tag{S13a-2}$$

$$\mathbf{p}_F^{I_c} = \mathbf{R}_{G\to I_c}(\mathbf{p}_F^G - \mathbf{p}_{I_c}^G) \tag{S13a-3}$$

### 2.3 特征表示

OpenVINS 支持 5 种特征表示，每种对应不同的 $\mathbf{H}_f = \frac{\partial \mathbf{p}_F^G}{\partial \boldsymbol{\lambda}}$：

| 表示 | $\boldsymbol{\lambda}$ | $\mathbf{H}_f$ |
|------|------|------|
| GLOBAL_3D | $\mathbf{p}_F^G$ | $\mathbf{I}_3$ |
| GLOBAL_INVERSE_DEPTH | $(\theta, \phi, \rho)$ | 球坐标→笛卡尔 Jacobian |
| ANCHORED_3D | $\mathbf{p}_F^A$ | $\mathbf{R}_{C_{\text{anc}}}^G$ |
| ANCHORED_INVERSE_DEPTH | $(\theta, \phi, \rho)$ 在锚点系 | $\mathbf{R}_{C_{\text{anc}}}^G \cdot \frac{\partial \mathbf{p}_F^A}{\partial \boldsymbol{\lambda}}$ |
| ANCHORED_MSCKF_INVERSE_DEPTH | $(\alpha, \beta, \rho)$ | $\mathbf{R}_{C_{\text{anc}}}^G \cdot \begin{bmatrix} 1/\rho & 0 & -\alpha/\rho^2 \\ 0 & 1/\rho & -\beta/\rho^2 \\ 0 & 0 & -1/\rho^2 \end{bmatrix}$ |

其中 $\alpha = p_x/p_z$，$\beta = p_y/p_z$，$\rho = 1/p_z$。

### 2.4 Jacobian 的链式法则

对每个观测 $c$，Jacobian 通过链式法则计算：

$$\frac{\partial(u,v)}{\partial\mathbf{x}} = \underbrace{\frac{\partial(u,v)}{\partial(u_n,v_n)}}_{\text{畸变}} \cdot \underbrace{\frac{\partial(u_n,v_n)}{\partial\mathbf{p}_F^{C_c}}}_{\text{投影}} \cdot \frac{\partial\mathbf{p}_F^{C_c}}{\partial\mathbf{x}} \tag{S13a-4}$$

其中 $\frac{\partial(u_n,v_n)}{\partial\mathbf{p}_F^{C_c}} = \begin{bmatrix} 1/Z & 0 & -X/Z^2 \\ 0 & 1/Z & -Y/Z^2 \end{bmatrix}$（$X,Y,Z = \mathbf{p}_F^{C_c}$）。

---

## 3. 代码实现

### 3.1 变量顺序构建

代码（[UpdaterHelper.cpp:201-261](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L201-L261)）：

```cpp
// 1. 收集所有涉及的变量
for (每个观测的相机) {
    if (do_calib_camera_pose)  { x_order.push_back(calibration); }
    if (do_calib_camera_intrinsics) { x_order.push_back(distortion); }
    for (每个观测时刻) { x_order.push_back(clone_Ci); }
}
// 2. 如果是锚定表示，添加锚点 clone 和其标定
if (is_relative_representation) {
    x_order.push_back(anchor_clone);
    if (do_calib_camera_pose) { x_order.push_back(anchor_calibration); }
}
```

**`x_order` 的作用**：记录 `H_x` 的每一列对应哪个状态变量。后续拼接到全局 `Hx_big` 时，通过 `x_order[i]->id()` 找到该变量在 `_Cov` 中的位置。

### 3.2 逐观测计算 Jacobian

代码（[UpdaterHelper.cpp:292-423](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L292-L423)）：

```cpp
for (每个相机) {
    for (每个观测) {
        // 1. 计算 p_FinCi（特征在相机坐标系中的位置）
        p_FinIi = R_GtoIi * (p_FinG - p_IiinG);
        p_FinCi = R_ItoC * p_FinIi + p_IinC;
        
        // 2. 投影 + 畸变 → 残差
        uv_norm = p_FinCi.head(2) / p_FinCi(2);
        uv_dist = distort(uv_norm);
        res[2*c] = uv_m - uv_dist;
        
        // 3. FEJ：用第一估计值重新计算中间量
        if (do_fej) {
            R_GtoIi = clone->Rot_fej();
            p_IiinG = clone->pos_fej();
            // 重新计算 p_FinCi
        }
        
        // 4. Jacobian 链式法则
        dz_dzn, dz_dzeta = compute_distort_jacobian(uv_norm);
        dzn_dpfc = [1/Z, 0, -X/Z²; 0, 1/Z, -Y/Z²];
        dpfc_dclone = [R_ItoC * skew(p_FinIi), -R_ItoC * R_GtoIi];
        
        // 5. 组装 H_f 和 H_x
        H_f[2*c] = dz_dpfc * dpfg_dlambda;
        H_x[2*c, clone_col] = dz_dpfc * dpfc_dclone;
        H_x[2*c, anchor_col] += dz_dpfg * dpfg_dx;  // 锚点 Jacobian
        H_x[2*c, calib_col] += dz_dpfc * dpfc_dcalib;
    }
}
```

### 3.3 FEJ 的使用位置

代码（[UpdaterHelper.cpp:354-363](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L354-L363)）：

```cpp
if (state->_options.do_fej) {
    R_GtoIi = clone_Ii->Rot_fej();     // 用 FEJ 的旋转
    p_IiinG = clone_Ii->pos_fej();     // 用 FEJ 的位置
    // 注意：R_ItoC 和 p_IinC 不用 FEJ（注释掉的代码）
    p_FinIi = R_GtoIi * (p_FinG_fej - p_IiinG);
    p_FinCi = R_ItoC * p_FinIi + p_IinC;
}
```

**为什么只用 clone 的 FEJ**：clone 位姿是状态变量，会随 EKF 更新变化。用 FEJ 保证 Jacobian 在一致性分析中正确。外参标定变化慢，用当前值即可。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{r}_c$ | 单观测残差 | `res.block(2*c, 0, 2, 1)` | `:348` |
| $\mathbf{p}_F^{C_c}$ | 特征在相机系中位置 | `p_FinCi` | `:337` |
| $\mathbf{R}_{G\to I_c}$ | clone 旋转 | `R_GtoIi` | `:330` |
| $\mathbf{p}_{I_c}^G$ | clone 位置 | `p_IiinG` | `:331` |
| $\frac{\partial(u_n,v_n)}{\partial\mathbf{p}_F^{C_c}}$ | 投影 Jacobian | `dzn_dpfc` | `:370-371` |
| $\frac{\partial\mathbf{p}_F^{C_c}}{\partial\delta\mathbf{x}_c}$ | clone 位姿 Jacobian | `dpfc_dclone` | `:377-379` |
| $\frac{\partial\mathbf{p}_F^G}{\partial\boldsymbol{\lambda}}$ | 特征表示 Jacobian | `dpfg_dlambda` | `:304` |
| $\mathbf{H}_f$ | 特征 Jacobian | `H_f` | `:389` |
| $\mathbf{H}_x$ | 状态 Jacobian | `H_x` | `:392` |

---

## 5. 关键设计决策

### 5.1 为什么 `dpfg_dlambda` 只算一次

`dpfg_dlambda = ∂p_F^G/∂λ` 只依赖于特征表示参数，与观测无关。所有观测共用同一个 Jacobian，避免重复计算。

### 5.2 锚定表示的特殊处理

锚定表示（`ANCHORED_*`）的特征 Jacobian 需要额外的锚点 clone Jacobian：

$$\frac{\partial\mathbf{p}_F^G}{\partial\delta\mathbf{x}_{\text{anchor}}} = -\mathbf{R}_{C_{\text{anc}}}^G \cdot [\mathbf{R}_{I\to C}\lfloor\mathbf{p}_F^A - \mathbf{p}_I^C\rfloor_\times \;\; \mathbf{I}_3]$$

代码（[UpdaterHelper.cpp:99-106](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L99-L106)）：

```cpp
H_anc.block(0, 0, 3, 3) = -R_GtoI.transpose() * skew_x(R_ItoC.transpose() * (p_FinA - p_IinC));
H_anc.block(0, 3, 3, 3) = I;
x_order.push_back(anchor_clone);
H_x.push_back(H_anc);
```

### 5.3 MSCKF 逆深度的 Jacobian

$$\frac{\partial\mathbf{p}_F^A}{\partial(\alpha,\beta,\rho)} = \begin{bmatrix} 1/\rho & 0 & -\alpha/\rho^2 \\ 0 & 1/\rho & -\beta/\rho^2 \\ 0 & 0 & -1/\rho^2 \end{bmatrix}$$

代码（[UpdaterHelper.cpp:167-170](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L167-L170)）：

```cpp
d_pfinA_dpinv << (1.0/rho), 0, -(1.0/(rho*rho))*alpha,
                 0, (1.0/rho), -(1.0/(rho*rho))*beta,
                 0, 0, -(1.0/(rho*rho));
H_f = R_CtoG * d_pfinA_dpinv;
```

---

## 6. 代码位置速查

| 功能 | 代码行 |
|------|--------|
| `get_feature_jacobian_representation` 入口 | `:32-40` |
| 全局 XYZ | `:36-40` |
| 全局逆深度 | `:43-71` |
| 锚定表示公共部分 | `:77-115` |
| 锚定 XYZ | `:118-121` |
| 锚定逆深度 | `:124-150` |
| MSCKF 逆深度 | `:153-172` |
| 单逆深度 | `:175-186` |
| `get_feature_jacobian_full` 入口 | `:192-200` |
| 变量顺序构建 | `:201-261` |
| 逐观测 Jacobian 计算 | `:292-423` |
| FEJ 替换 | `:354-363` |
| 畸变 Jacobian | `:366-367` |
| 投影 Jacobian | `:370-371` |
| clone 位姿 Jacobian | `:377-379` |
| 链式法则组装 | `:384-398` |

---

## 7. 同类实现对比

| 特征维度 | OpenVINS | VINS-Mono | OKVIS |
|---------|---------|-----------|-------|
| 特征表示 | 5 种（全局/锚定 × XYZ/逆深度） | 逆深度 | 逆深度 |
| Jacobian 计算 | 链式法则 + FEJ | 解析 | 自动微分 |
| 畸变模型 | radtan/fisheye | radtan | radtan |
| 外参/内参 Jacobian | 支持 | 支持 | 支持 |

### 设计取舍分析

- 选择**5 种特征表示**是为了适应不同场景：远距离特征用逆深度（数值稳定），近距离用 XYZ（精度高）
- 选择**FEJ 只用 clone 位姿**是因为外参变化慢，用当前值不影响一致性
- 选择**链式法则而非自动微分**是因为解析 Jacobian 更快，且 OpenVINS 的数学结构清晰
