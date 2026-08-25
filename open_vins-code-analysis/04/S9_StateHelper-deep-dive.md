# StateHelper 协方差前推 精读报告 (S9)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/StateHelper.cpp`
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (误差态 EKF 协方差更新 Eq.103)、Mourikis & Roumeliotis 2007 MSCKF [20] (滑窗)、Li & Mourikis 2013 IJRR [27] (时间偏移)、Solà 2017 [47] (Givens QR)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S9（Phase 3 函数级路线，阶段 C 协方差前推）

> 上游：S5（`propagate_and_clone` `:130` 调用本函数）。本篇仅聚焦 **`StateHelper::EKFPropagation` (:36)** —— 把 S7 累加的 $\Phi/\mathbf{Q}_d$ 应用到全局协方差 $\mathbf{P}$，即误差态 EKF 预测步的协方差更新 $\mathbf{P}'=\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$。**克隆/增广 → S10**，**边缘化 → S11**，**新变量初始化 → S12** 已拆为独立文档。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
propagate_and_clone (S5 :130)
        │  (Phi_summed, Qd_summed, Phi_order)
        ▼
StateHelper::EKFPropagation(state, order_NEW, order_OLD, Phi, Q)  @ :36  ★ 本篇
        │  改写 state->_Cov (仅 IMU 相关块)
        ▼
augment_clone (S10 :579)  ← 后续滑窗克隆
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 预测第 3 步（协方差前推） |
| 输入 | `order_NEW`/`order_OLD`（被前推变量，均=IMU 块）、`Phi`（状态转移）、`Q`（离散噪声） |
| 输出 | 改写 `state->_Cov` |
| 调用方 | `Propagator::propagate_and_clone` (S5 :130) |
| 被调用方 | 无（纯 Eigen 块操作） |

### 1.2 一句话概括

`EKFPropagation` 是误差态 EKF 预测步的协方差更新：用 `order_OLD` 定位 $\Phi$ 列、逐块累加 $\mathbf{P}\Phi^\top$，再算 $\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$，写回 `_Cov` 的中心块与交叉块，并校验连续布局与正定性。它只动 IMU 相关行/列，clone/SLAM 路标的交叉协方差同步随 $\Phi$ 前推。

---

## Section 2: 核心数据结构与数学模型

### 2.1 函数签名

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:36-114
void StateHelper::EKFPropagation(state,
    const vector<Type> &order_NEW,   // 前推后变量 (S5 中 = IMU + 可选内参块)
    const vector<Type> &order_OLD,   // 前推前变量 (同上)
    const MatrixXd &Phi,             // 状态转移 (来自 S7 Phi_summed)
    const MatrixXd &Q);              // 离散噪声 (来自 S7 Qd_summed)
```

### 2.2 数学模型

**问题定义**（Trawny 2005 [40] Eq.103）：误差态 EKF 预测步的协方差更新为

$$\mathbf{P}' = \tilde{\boldsymbol{\Phi}}\,\mathbf{P}\,\tilde{\boldsymbol{\Phi}}^\top + \tilde{\mathbf{Q}}_d \tag{S9-1}$$

其中 $\tilde{\boldsymbol{\Phi}}$ 和 $\tilde{\mathbf{Q}}_d$ 是扩展到 $n\times n$ 的传播矩阵和噪声矩阵（非传播变量对应位置为单位阵和零）。

**分块结构**：将 $\mathbf{P}$ 按"传播变量"（下标 $p$）和"其他变量"（下标 $o$）分块：

$$\mathbf{P} = \begin{bmatrix} \mathbf{P}_{oo} & \mathbf{P}_{op} \\ \mathbf{P}_{po} & \mathbf{P}_{pp} \end{bmatrix} \tag{S9-2}$$

其中 $\mathbf{P}_{pp}$ 是 $m\times m$ 子块（$m$ = 传播变量总误差维度，通常 $m=15$ 对应 IMU）。传播只更新涉及 $p$ 的块：

$$\mathbf{P}_{pp}' = \boldsymbol{\Phi}\,\mathbf{P}_{pp}\,\boldsymbol{\Phi}^\top + \mathbf{Q}_d \tag{S9-3}$$

$$\mathbf{P}_{op}' = \mathbf{P}_{op}\,\boldsymbol{\Phi}^\top \tag{S9-4}$$

$$\mathbf{P}_{po}' = \boldsymbol{\Phi}\,\mathbf{P}_{po} = (\mathbf{P}_{op}')^\top \tag{S9-5}$$

$$\mathbf{P}_{oo}' = \mathbf{P}_{oo} \quad \text{（不变）} \tag{S9-6}$$

其中 $\boldsymbol{\Phi}$ 和 $\mathbf{Q}_d$ 是 $m\times m$ 的局部矩阵（代码中的 `Phi` 和 `Q`）。

**直接计算的代价**：直接构造 $n\times n$ 的 $\tilde{\boldsymbol{\Phi}}$ 并计算 $\tilde{\boldsymbol{\Phi}}\mathbf{P}\tilde{\boldsymbol{\Phi}}^\top$ 的复杂度为 $O(n^3)$。当 $n\approx 200$、$m\approx 15$ 时，绝大部分运算是在与单位阵/零做乘法。

### 2.3 分块计算推导

**核心思路**：定义中间矩阵

$$\mathbf{M} = \mathbf{P}_{\text{all},\,p}\,\boldsymbol{\Phi}^\top \in \mathbb{R}^{n\times m} \tag{S9-7}$$

即取 $\mathbf{P}$ 中传播变量对应的所有行、$m$ 列，右乘 $\boldsymbol{\Phi}^\top$。展开分块：

$$\mathbf{M} = \begin{bmatrix} \mathbf{P}_{op} \\ \mathbf{P}_{pp} \end{bmatrix} \boldsymbol{\Phi}^\top = \begin{bmatrix} \mathbf{P}_{op}\boldsymbol{\Phi}^\top \\ \mathbf{P}_{pp}\boldsymbol{\Phi}^\top \end{bmatrix} = \begin{bmatrix} \mathbf{P}_{op}' \\ \mathbf{P}_{pp}\boldsymbol{\Phi}^\top \end{bmatrix} \tag{S9-8}$$

$\mathbf{M}$ 的上半部分恰好就是更新后的交叉块 $\mathbf{P}_{op}'$（公式 (S9-4)）。再取 $\mathbf{M}$ 的下半部分左乘 $\boldsymbol{\Phi}$：

$$\boldsymbol{\Phi}(\mathbf{P}_{pp}\boldsymbol{\Phi}^\top) + \mathbf{Q}_d = \boldsymbol{\Phi}\mathbf{P}_{pp}\boldsymbol{\Phi}^\top + \mathbf{Q}_d = \mathbf{P}_{pp}' \tag{S9-9}$$

**三步算法**：

**Step 1**：计算 $\mathbf{M} = \mathbf{P}_{\text{all},\,p}\,\boldsymbol{\Phi}^\top$（$n\times m$）。设 `order_OLD` 中有 $K$ 个变量，第 $i$ 个变量的误差维度为 $s_i$（$\sum s_i = m$），在 $\mathbf{P}$ 中的起始列为 $\text{id}_i$，在 $\boldsymbol{\Phi}$ 中的起始列为 $\text{Phi\_id}_i$：

$$\mathbf{M} = \sum_{i=1}^{K} \mathbf{P}_{*,\,\text{id}_i:\,\text{id}_i+s_i} \cdot \boldsymbol{\Phi}_{\text{Phi\_id}_i:\,\text{Phi\_id}_i+s_i,\,*}^\top \tag{S9-10}$$

**Step 2**：计算 $\mathbf{P}_{pp}' = \boldsymbol{\Phi}\,\mathbf{M}_{p} + \mathbf{Q}_d$（$m\times m$）：

$$\mathbf{P}_{pp}' = \mathbf{Q}_d + \sum_{i=1}^{K} \boldsymbol{\Phi}_{*,\,\text{Phi\_id}_i:\,\text{Phi\_id}_i+s_i} \cdot \mathbf{M}_{\text{id}_i:\,\text{id}_i+s_i,\,*} \tag{S9-11}$$

**Step 3**：写回 $\mathbf{P}$。利用对称性 $\mathbf{P}_{po}' = (\mathbf{P}_{op}')^\top = \mathbf{M}^\top$：

$$\mathbf{P}' = \begin{bmatrix} \mathbf{P}_{oo} & \mathbf{M}_{\text{upper}} \\ \mathbf{M}_{\text{upper}}^\top & \mathbf{P}_{pp}' \end{bmatrix} \tag{S9-12}$$

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 协方差前推 $\mathbf{P}'=\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$

#### 理论推导

误差态 EKF 预测步的协方差更新（Trawny 2005 [40] Eq.103），即公式 (S9-1) 的局部形式：

$$\mathbf{P}' = \boldsymbol{\Phi}\,\mathbf{P}\,\boldsymbol{\Phi}^\top + \mathbf{Q}_d \tag{S9-13}$$

由于 $\boldsymbol{\Phi}$ 只作用于 IMU 相关的 $m$ 个变量（$m \ll n$），直接构造 $n\times n$ 扩展矩阵 $\tilde{\boldsymbol{\Phi}}$ 做全量乘法是 $O(n^3)$ 的浪费。

**分块计算**（Anderson & Moore 1979）：

定义 $\mathbf{M} = \mathbf{P}_{\text{all},\,p}\,\boldsymbol{\Phi}^\top \in \mathbb{R}^{n\times m}$（公式 (S9-7)），则：

- $\mathbf{M}$ 的上半部分（$n-m$ 行）= $\mathbf{P}_{op}\boldsymbol{\Phi}^\top = \mathbf{P}_{op}'$（新交叉块，公式 (S9-4)）
- $\mathbf{M}$ 的下半部分（$m$ 行）= $\mathbf{P}_{pp}\boldsymbol{\Phi}^\top$，左乘 $\boldsymbol{\Phi}$ 再加 $\mathbf{Q}_d$ 得 $\mathbf{P}_{pp}'$（公式 (S9-3)(S9-9)）

逐变量累加实现（公式 (S9-10)(S9-11)）：设 `order_OLD` 中第 $i$ 个变量的误差维度为 $s_i$，在 $\mathbf{P}$ 中起始列为 $\text{id}_i$，在 $\boldsymbol{\Phi}$ 中起始列为 $\text{Phi\_id}_i$。

**复杂度**：分块算法 $O(n\cdot m^2 + m^3)$ vs 朴素 $O(n^3)$。当 $n=200, m=15$ 时，分块减少约 99% 运算量。

**对称性利用**：$\mathbf{P}_{po}' = (\mathbf{P}_{op}')^\top = \mathbf{M}^\top$，只需计算一次 $\mathbf{M}$。

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:36-114
void StateHelper::EKFPropagation(state, order_NEW, order_OLD, Phi, Q) {
    // 1. 防御: 空数组直接崩溃
    if (order_NEW.empty() || order_OLD.empty()) { ... std::exit(...); }   // :41-44

    // 2. 校验 NEW 在内存中连续 (否则块写入错位)
    int size_order_NEW = order_NEW[0]->size();                            // :47
    for (i=0; i<order_NEW.size()-1; i++) {                                // :48
      if (order_NEW[i]->id()+order_NEW[i]->size() != order_NEW[i+1]->id())
        { ... std::exit(...); }                                           // :49-54 非连续→崩
      size_order_NEW += order_NEW[i+1]->size();                           // :55
    }
    int size_order_OLD = ...;                                             // :59-62
    assert(size_order_NEW==Phi.rows());                                   // :65
    assert(size_order_OLD==Phi.cols());                                   // :66
    assert(size_order_NEW==Q.cols() && ==Q.rows());                       // :67-68

    // 3. 给每个 OLD 变量定位其在 Phi 中的列起点 Phi_id
    int cur=0; std::vector<int> Phi_id;
    for (var : order_OLD) { Phi_id.push_back(cur); cur+=var->size(); }    // :72-76

    // 4. Cov_PhiT = P * Phiᵀ   (逐 OLD 变量块累加)
    Eigen::MatrixXd Cov_PhiT=Zero(_Cov.rows(), Phi.rows());               // :80
    for (var : order_OLD)                                                 // :81
      Cov_PhiT.noalias() += _Cov.block(0,var->id(),_Cov.rows(),var->size())
                            * Phi.block(0,Phi_id[i],Phi.rows(),var->size()).transpose();  // :83-84

    // 5. Φ*Cov*Φᵀ + Q = 逐 OLD 块累加
    Eigen::MatrixXd Phi_Cov_PhiT = Q.selfadjointView<Upper>();            // :88  (Q 对称)
    for (var : order_OLD)                                                 // :89
      Phi_Cov_PhiT.noalias() += Phi.block(0,Phi_id[i],Phi.rows(),var->size())
                                * Cov_PhiT.block(var->id(),0,var->size(),Phi.rows());  // :91

    // 6. 写回 _Cov 的三个区域 (start_id = order_NEW[0]->id())
    int start_id = order_NEW[0]->id();                                    // :95
    _Cov.block(start_id,0,   phi_size,total_size) = Cov_PhiT.transpose(); // :98  上块 (PΦᵀ)ᵀ
    _Cov.block(0,start_id,   total_size,phi_size) = Cov_PhiT;            // :99  左块 PΦᵀ
    _Cov.block(start_id,start_id, phi_size,phi_size) = Phi_Cov_PhiT;     // :100 中心块 ΦPΦᵀ+Q

    // 7. 对称正定检查: 对角线出现负值即崩溃
    Eigen::VectorXd diags=_Cov.diagonal();                                // :103
    for (i=0;i<diags.rows();i++) if (diags(i)<0) { ... std::exit; }       // :105-113
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 | 公式编号 |
|--------|---------|------|---------|
| $\mathbf{M} = \mathbf{P}\boldsymbol{\Phi}^\top$（Step 1） | `Cov_PhiT += _Cov.block(...)*Phi.block(...).transpose()` | `:83-84` | (S9-10) |
| $\mathbf{P}_{pp}' = \boldsymbol{\Phi}\mathbf{M}_p + \mathbf{Q}_d$（Step 2） | `Phi_Cov_PhiT += Phi*Cov_PhiT.block(...)`（起始于 `Q`） | `:88-91` | (S9-11) |
| 写回 $\mathbf{P}_{op}'$（交叉块） | `_Cov.block(0,start_id)=Cov_PhiT` | `:99` | (S9-4) |
| 写回 $\mathbf{P}_{po}'$（对称） | `_Cov.block(start_id,0)=Cov_PhiT.transpose()` | `:98` | (S9-5) |
| 写回 $\mathbf{P}_{pp}'$（对角块） | `_Cov.block(start_id,start_id)=Phi_Cov_PhiT` | `:100` | (S9-3) |

> **设计要点**：因为 `order_NEW == order_OLD`（都是 IMU 块），更新只覆盖协方差中对应 IMU 的**左上连续子块**及其与全局的交叉块。注意它**没有**重建整个 `_Cov`——只重算 IMU 相关行/列，其余变量（clone、SLAM 路标）与 IMU 的交叉协方差通过 `Cov_PhiT` 同步前推（`:98-99` 把 $\mathbf{P}\Phi^\top$ 转置写入，使交叉块也随 $\Phi$ 前推）。

> **连续布局 `assert`**：`:48-54` 强制 `order_NEW` 中各变量在 `_Cov` 中按 `_id` 连续无空洞。这是整个 StateHelper 的不变式（由 S11 的 `set_local_id` 重排维护）。

---

## Section 4: 函数调用链

```
Propagator::propagate_and_clone (S5 :130)
  └─ StateHelper::EKFPropagation(state, Phi_order, Phi_order, Phi_summed, Qd_summed) @ :36  ★ 本篇
        └─ 纯 Eigen 块操作, 无子调用
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 | 下沉 |
|--------|-----------|------|------|
| `EKFPropagation` | `StateHelper.cpp:36` | 预测前推 $\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$ | **S9 (本篇)** |
| `augment_clone` | `StateHelper.cpp:579` | 滑窗克隆 | **S10** |
| `marginalize` | `StateHelper.cpp:271` | 删单变量 | **S11** |
| `initialize` | `StateHelper.cpp:393` | 新变量初始化 | **S12** |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `order_NEW`/`order_OLD` | S5 构造 | S5 `:116-129` | 内部 | 被前推变量 | ★★ 中 |

> 本函数无独立外部参数，行为完全由 S5 传入的 `Phi_summed`/`Qd_summed` 决定（见 S7）。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 负值对角线硬崩溃 | ★ 中 | 数值病态 | `:105` 直接 `std::exit` | 可降级为重置协方差或告警 |
| 2 | 非连续布局硬崩溃 | ★ 中 | `set_local_id` 错乱 | `:49` 直接 `std::exit` | 理论上不应发生（由 S11 保证） |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `EKFPropagation` | VINS-Mono | OKVIS |
|---------|--------------------------|-----------|-------|
| 协方差更新 | $\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$ 块操作 | 边缘化先验 | $\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$ |
| 连续布局假设 | 强制（assert） | 无固定布局 | 无固定布局 |
| 只更新相关块 | 是（IMU 块） | 否（全量） | 是 |

### 设计取舍分析

- 选择 **按 `order` 块定位 + 只更新相关块** 是因为：EKF 框架中 $\Phi$ 仅作用于 IMU 块，全量重算 $\mathbf{P}$ 是浪费；块操作 O(N²·n_old) 远优于全量 O(N³)。
- 选择 **连续布局 + assert** 是因为：Eigen 块寻址要求变量在矩阵中连续；由 S11 边缘化时的 `set_local_id` 重排保证不变式，使 `assert` 既是校验也是文档。
