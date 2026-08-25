# StateHelper 新变量初始化 精读报告 (S12)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/StateHelper.cpp`
> **对应论文**: Solà 2017 [47] (Givens QR 初始化)、MSCKF-SLAM [32] (零空间投影 + χ² 检验)、Mourikis & Roumeliotis 2007 MSCKF [20]
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S12（Phase 3 函数级路线，阶段 C 新变量加入状态）

> 上游：S9/S11（协方差增删机制）。本篇聚焦 **`StateHelper::initialize` (:393)** 与 **`initialize_invertible` (:484)** —— SLAM 路标/地图点首次可见时的"入状态"路径：用 Givens QR 分离可逆系统与零空间投影系统，做 95% χ² Mahalanobis 检验拒绝异常，再把新变量协方差增广进 `_Cov`。是 S14 (UpdaterSLAM) 的后端支撑。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
UpdaterSLAM::delayed_init (S15) / UpdaterSLAM::update (S14)
        │  (新可见路标)
        ▼
StateHelper::initialize(state, new_var, H_order, H_R, H_L, R, res, mult)  @ :393  ★ 本篇
        ├─ Givens QR 分离 (H_L 含新变量列)
        ├─ χ² Mahalanobis 检验 (投影系统)
        ├─ initialize_invertible(state, new_var, H_order, Hxinit, H_finit, Rinit, resinit) @ :484 ★ 本篇
        │      └─ 增广 _Cov + new_var->update + set_local_id
        └─ EKFUpdate(state, H_order, Hup, resup, Rup)  (投影部分, 见 S13/S9)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 新变量加入状态（SLAM 路标首次可见） |
| 输入 | `new_variable`、`H_R`/`H_L`（QR 后的右/左系统）、`R`、`res`、`chi_2_mult` |
| 输出 | 增广 `state->_Cov` + `_variables`，新变量入状态 |
| 调用方 | `UpdaterSLAM` (S14/S15) |
| 被调用方 | `initialize_invertible` (:484)、`EKFUpdate`（S9/S13）、`Type::update`/`clone`/`set_local_id`（S1） |

### 1.2 一句话概括

`initialize` 用 Givens QR 把测量系统分离为"可逆初始化系统"（直接解新变量）和"零空间投影更新系统"（投影到旧变量更新），先做 χ² 检验拒绝异常，再分别用 `initialize_invertible` 增广新变量、用 `EKFUpdate` 更新旧变量。保证 SLAM 路标首次可见时的一致性与鲁棒性（[32][47]）。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp
bool StateHelper::initialize(state, new_var, H_order, H_R, H_L, R, res, chi_2_mult=1.0);  // :393
void StateHelper::initialize_invertible(state, new_var, H_order, H_R, H_L, R, res);       // :484
// H_L: 含新变量列的 Jacobian (new_var_size 列)
// H_R: 旧变量 Jacobian (H_order.size() 列)
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 初始化主流程 `initialize`（Givens QR + χ²）

#### 理论推导

**问题定义**（OpenVINS 论文 Eq.18-22）：设系统当前状态误差为 $\delta\mathbf{x} \in \mathbb{R}^m$，新变量误差为 $\delta\mathbf{x}_f \in \mathbb{R}^p$（$p$ 为新变量误差维度，如 3D 路标 $p=3$，单逆深度 $p=1$）。测量残差 $\mathbf{r} \in \mathbb{R}^k$ 的线性化模型为：

$$\mathbf{r} = \mathbf{H}_R\,\delta\mathbf{x} + \mathbf{H}_L\,\delta\mathbf{x}_f + \mathbf{n}, \quad \mathbf{n} \sim \mathcal{N}(\mathbf{0},\, \mathbf{R}) \tag{S12-1}$$

其中 $\mathbf{H}_R \in \mathbb{R}^{k \times m}$ 为旧变量 Jacobian，$\mathbf{H}_L \in \mathbb{R}^{k \times p}$ 为新变量 Jacobian，$\mathbf{R} = \sigma^2 \mathbf{I}_k$ 为各向同性噪声。

**Givens QR 分离**（Solà 2017 [47] / MSCKF-SLAM [32]）：对 $\mathbf{H}_L$ 做 QR 分解：

$$\mathbf{Q}^\top \mathbf{H}_L = \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix} \tag{S12-2}$$

其中 $\mathbf{T} \in \mathbb{R}^{p \times p}$ 上三角（满秩），$\mathbf{Q} = [\mathbf{Q}_1 \;\; \mathbf{Q}_2]$，$\mathbf{Q}_1 \in \mathbb{R}^{k \times p}$，$\mathbf{Q}_2 \in \mathbb{R}^{k \times (k-p)}$。将 $\mathbf{Q}^\top$ 左乘方程 (18)：

$$\begin{bmatrix} \mathbf{Q}_1^\top\mathbf{r} \\ \mathbf{Q}_2^\top\mathbf{r} \end{bmatrix} = \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{H}_R \\ \mathbf{Q}_2^\top\mathbf{H}_R \end{bmatrix}\delta\mathbf{x} + \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix}\delta\mathbf{x}_f + \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{n} \\ \mathbf{Q}_2^\top\mathbf{n} \end{bmatrix} \tag{S12-3}$$

分离为两个子系统：

- **可逆初始化系统**（$p$ 个方程）：$\mathbf{Q}_1^\top\mathbf{r} = \mathbf{Q}_1^\top\mathbf{H}_R\,\delta\mathbf{x} + \mathbf{T}\,\delta\mathbf{x}_f + \mathbf{Q}_1^\top\mathbf{n}$
- **零空间投影更新系统**（$k-p$ 个方程）：$\mathbf{Q}_2^\top\mathbf{r} = \mathbf{Q}_2^\top\mathbf{H}_R\,\delta\mathbf{x} + \mathbf{Q}_2^\top\mathbf{n}$

由于 $\mathbf{R} = \sigma^2\mathbf{I}$，正交变换保持噪声统计性质：$\mathbb{E}[\mathbf{Q}_i^\top\mathbf{n}\mathbf{n}^\top\mathbf{Q}_j] = \sigma^2\delta_{ij}\mathbf{I}$。

**χ² 检验**（对投影系统）：

$$\chi^2 = \mathbf{r}_{up}^\top \mathbf{S}^{-1} \mathbf{r}_{up} \sim \chi^2_{k-p}, \quad \mathbf{S} = \mathbf{H}_{up}\mathbf{P}_{\text{small}}\mathbf{H}_{up}^\top + \mathbf{R}_{up} \tag{S12-4}$$

**可逆初始化** `initialize_invertible`：从可逆系统解出新变量：

$$\delta\mathbf{x}_f = \mathbf{T}^{-1}(\mathbf{r}_1 - \mathbf{H}_{x1}\,\delta\mathbf{x}) \approx \mathbf{T}^{-1}\mathbf{r}_1 \tag{S12-5}$$

（三角化前已通过条件 Gauss-Newton 求解过初值，$\mathbf{r}_1$ 中已扣除 $\mathbf{H}_{x1}\delta\mathbf{x}$ 的贡献。）

**新变量协方差**：等效测量噪声为 $\mathbf{H}_{x1}\delta\mathbf{x} + \mathbf{n}_1$，其协方差 $\mathbf{M} = \mathbf{H}_{x1}\mathbf{P}_{\text{small}}\mathbf{H}_{x1}^\top + \mathbf{R}_1$。由误差传播：

$$\mathbf{P}_{ff} = \mathbf{T}^{-1}\mathbf{M}\mathbf{T}^{-\top} \tag{S12-6}$$

**交叉协方差**：定义 $\mathbf{M}_a = \mathbf{P}\mathbf{H}_R^\top \in \mathbb{R}^{n \times p}$（全状态对 $\mathbf{H}_R\delta\mathbf{x}$ 的交叉协方差），则：

$$\mathbf{P}_{\text{all},f} = -\mathbf{M}_a\mathbf{T}^{-\top} \tag{S12-7}$$

负号来源于 $\delta\mathbf{x}_f$ 对 $\delta\mathbf{x}$ 的负依赖（$\mathbf{r}_1$ 固定时，$\delta\mathbf{x}$ 增大则 $\delta\mathbf{x}_f$ 减小）。

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:393-482
bool StateHelper::initialize(state, new_variable, H_order, H_R, H_L, R, res, chi_2_mult=1.0) {
    // 1. 防御: new_variable 不能已在 _variables 中
    if (find(_variables.begin(),_variables.end(),new_variable)!=_variables.end()) { ...std::exit;} // :398-402
    // 2. 校验 R 各向同性对角 (R(0,0) 全等且非对角为0)
    for (r,c) { if(r==c && R(0,0)!=R(r,c)) exit; else if(r!=c && R(r,c)!=0) exit; }  // :408-420

    // 3. Givens QR: 把 H_L (含新变量列) 旋成上三角, 同步作用于 H_R/res/R
    Eigen::JacobiRotation<double> G;
    for (n=0; n<H_L.cols(); n++)                                            // :430
      for (m=H_L.rows()-1; m>n; m--) {                                      // :431
        G.makeGivens(H_L(m-1,n), H_L(m,n));                                 // :433
        H_L.block(m-1,n,2,H_L.cols()-n).applyOnTheLeft(0,1,G.adjoint());    // :437
        res.block(m-1,0,2,1).applyOnTheLeft(0,1,G.adjoint());               // :438
        H_R.block(m-1,0,2,H_R.cols()).applyOnTheLeft(0,1,G.adjoint());      // :439
      }

    // 4. 分离: 上 new_var_size 行=初始化系统, 余下=投影更新系统
    Hxinit=H_R.block(0,0,new_var_size,H_R.cols());  H_finit=H_L.block(0,0,new_var_size,new_var_size); // :445-446
    resinit=res.block(0,0,new_var_size,1);          Rinit=R.block(0,0,new_var_size,new_var_size);     // :447-448
    Hup=H_R.block(new_var_size,0,...);  resup=res.block(new_var_size,0,...);  Rup=R.block(...);        // :451-453

    // 5. χ² Mahalanobis 检验 (投影系统部分)
    P_up = get_marginal_covariance(state, H_order);                          // :459
    S = Hup*P_up*Hup.transpose()+Rup;                                       // :462
    chi2 = resup.dot(S.llt().solve(resup));                                  // :463
    boost::math::chi_squared chi_dist(res.rows());                           // :466
    double chi2_check = quantile(chi_dist, 0.95);                           // :467 95% 阈值
    if (chi2 > chi_2_mult * chi2_check) return false;                        // :468-469 拒绝

    // 6. 可逆初始化新变量 + 投影部分做 EKF 更新
    initialize_invertible(state, new_variable, H_order, Hxinit, H_finit, Rinit, resinit); // :475
    if (Hup.rows()>0) EKFUpdate(state, H_order, Hup, resup, Rup);            // :478-480
    return true;
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 | 公式编号 |
|--------|---------|------|---------|
| 测量模型 | `H_R, H_L, R, res` | `:393-395` | (18) |
| Givens QR 分离 | `G.makeGivens(...); applyOnTheLeft(...)` | `:430-439` | (19) |
| 分离为初始化/投影系统 | `Hxinit/H_finit/resinit` 与 `Hup/resup/Rup` | `:445-453` | (20) |
| $\mathbf{M} = \mathbf{P}\mathbf{H}_R^\top$ | `M_a.block(var.id,0,...) = M_i` | `:519-541` | (24) |
| $\mathbf{M} = \mathbf{H}_R\mathbf{P}_{\text{small}}\mathbf{H}_R^\top + \mathbf{R}$ | `M = H_R*P_small*H_R^T + R` | `:549-551` | (23) 的 $\mathbf{M}$ |
| $\mathbf{P}_{ff} = \mathbf{H}_L^{-1}\mathbf{M}\mathbf{H}_L^{-\top}$ | `P_LL = H_Linv*M*H_Linv^T` | `:557` | (23) |
| 交叉块 $\mathbf{P}_{iL} = -\mathbf{M}_a\mathbf{H}_L^{-\top}$ | `_Cov.block(0,oldSize) = -M_a*H_Linv^T` | `:562` | (24) |
| $\delta\mathbf{x}_f = \mathbf{H}_L^{-1}\mathbf{r}$ | `new_variable->update(H_Linv * res)` | `:568` | (22) |
| $\chi^2 = \mathbf{r}_{up}^\top\mathbf{S}^{-1}\mathbf{r}_{up}$ | `chi2 = resup.dot(S.llt().solve(resup))` | `:463` | (21) |
| 95% 阈值 | `quantile(chi_squared(res.rows()), 0.95)` | `:466-467` | — |

### 3.2 可逆初始化 `initialize_invertible`

#### 对应代码

```c++
// 文件: ov_msckf/src/state/StateHelper.cpp:484-577
void StateHelper::initialize_invertible(state, new_variable, H_order, H_R, H_L, R, res) {
    // 同款防御 + 各向同性 R 校验 (:488-511)
    // 1. M_a = P*H_Rᵀ  (与 EKFUpdate 同构, 但用 H_R)
    ... M_a.block(var.id,0,...) = M_i (含 H_R) ...                          // :519-541
    // 2. M = H_R*P_small*H_Rᵀ + R
    M.triangularView<Upper>()=H_R*P_small*H_R.transpose(); M+=R;            // :549-551
    // 3. 新变量协方差 P_LL = H_L⁻¹ * M * H_L⁻ᵀ
    H_Linv = H_L.inverse();                                                 // :556
    P_LL = H_Linv * M.selfadjointView<Upper>() * H_Linv.transpose();        // :557
    // 4. 扩协方差, 写交叉块 + 新变量自协方差
    oldSize=_Cov.rows();
    _Cov.conservativeResizeLike(Zero(oldSize+size, oldSize+size));          // :561
    _Cov.block(0,oldSize,oldSize,size).noalias() = -M_a*H_Linv.transpose(); // :562 交叉 (负)
    _Cov.block(oldSize,0,size,oldSize) = _Cov.block(0,oldSize,...).transpose(); // :563
    _Cov.block(oldSize,oldSize,size,size) = P_LL;                           // :564 自协方差
    // 5. 更新新变量均值 (应≈0, 因已用 CGN 解过初值)
    new_variable->update(H_Linv * res);                                     // :568
    new_variable->set_local_id(oldSize);                                    // :571 ★ 入末尾
    state->_variables.push_back(new_variable);                              // :572
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{P}_{LL}=\mathbf{H}_L^{-1}\mathbf{M}\mathbf{H}_L^{-\top}$ | `P_LL=H_Linv*M*H_Linvᵀ` | `:557` |
| 交叉块 $\mathbf{P}_{iL}=-\mathbf{M}_a\mathbf{H}_L^{-\top}$ | `_Cov.block(0,oldSize)=-M_a*H_Linvᵀ` | `:562` |
| 新变量入末尾 | `new_variable->set_local_id(oldSize)` | `:571` |

> **设计要点**：`initialize` 是 `EKFUpdate`（S9/S13）的"超集"——它先 QR 分离出可逆部分直接解新变量（`initialize_invertible`），再把零空间投影部分交给普通 `EKFUpdate`。`chi_2_mult`（默认 1.0）可调检验严格度。`:568` 的 `update` 再次调用 S1 boxplus。SLAM 路标、地图点首次可见都走这条路径。

> **R 各向同性限制**：`:408-420` 强制 `R` 为均匀对角，否则直接 `std::exit`——这是 Givens QR 分离的前提（非各向同性噪声需先白化，代码未实现，见 S11 §6 #1 同类技术债）。

---

## Section 4: 函数调用链

```
UpdaterSLAM (S14/S15)
  └─ StateHelper::initialize(state, new_var, H_order, H_R, H_L, R, res) @ :393  ★ 本篇
        ├─ Givens QR (Eigen JacobiRotation) @ :430
        ├─ get_marginal_covariance @ :459
        ├─ initialize_invertible @ :484  ★ 本篇
        │    └─ new_variable->update(H_Linv*res) @ :568  (S1 boxplus)
        │    └─ new_variable->set_local_id(oldSize) @ :571 (S1 连续 id)
        └─ EKFUpdate @ :478 (投影部分, → S9/S13)
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 |
|--------|-----------|------|
| `initialize` | `StateHelper.cpp:393` | 新变量(χ²+QR) |
| `initialize_invertible` | `StateHelper.cpp:484` | 可逆初始化 |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `chi_2_mult` | 1.0 | `initialize` 调用方 | 一致性 | χ² 检验倍率 | ★★ 中 |

> 调大→更容忍异常观测（漏检少但可能接纳 outlier）；调小→更严格（鲁棒但易拒绝正确初始化）。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `initialize` 要求 R 各向同性 | ★ 低 | 传入非对角/非均匀 R | `:408-420` 直接 exit | 放宽校验或支持一般 R |
| 2 | χ² 阈值固定 95% | ★ 低 | 高维残差 | `:467` quantile 0.95 | 可由 `chi_2_mult` 调 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `initialize` | VINS-Mono | OKVIS | MSCKF-SLAM [32] |
|---------|-----------------------|-----------|-------|-----------------|
| 新变量加入 | Givens QR + χ² ([47][32]) | 三角化 | 三角化 | QR + χ² |
| 零空间投影 | 是（Hup 部分） | 无(图优化) | 部分 | 是 |
| 一致性保障 | χ² 检验拒绝异常 | BA 残差 | 残差 | χ² |

### 设计取舍分析

- 选择 **Givens QR + χ²** 而非直接三角化，是因为：QR 把新变量系统分离为可逆+投影两部分，显式 χ² 检验可拒绝异常初始化（[32] MSCKF-SLAM 一致性保障）。
- 选择 **可逆部分直接解 + 投影部分交 EKFUpdate**，是因为：复用 S9 的卡尔曼更新逻辑，避免为初始化单独实现一套协方差更新。
