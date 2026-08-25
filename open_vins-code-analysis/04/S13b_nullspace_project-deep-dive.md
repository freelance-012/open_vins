# UpdaterHelper::nullspace_project_inplace 深度精读 (S13b)

> **生成日期**：2026-08-19
> **前置文档**：`S13_UpdaterMSCKF-deep-dive.md`（主文档）、`S13a_get_feature_jacobian_full-deep-dive.md`
> **核心代码**：`ov_msckf/src/update/UpdaterHelper.cpp:426-454`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Golub & Van Loan "Matrix Computations" 4th Ed. [48] Algorithm 5.2.4；MSCKF 论文 Section III-B

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
UpdaterMSCKF::update()
        │
        ├─ get_feature_jacobian_full() → H_f, H_x, res  （S13a）
        │
        ▼
nullspace_project_inplace(H_f, H_x, res)  @ :426  ★ 本篇
        │
        ├─ Givens QR 分解 H_f
        ├─ 消去 H_f 的下三角
        └─ 取 H_x, res 的下半部分（左零空间）
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | MSCKF 更新中每个特征的 Jacobian 后处理 |
| 输入 | `H_f`（$2N \times p$）、`H_x`（$2N \times m$）、`res`（$2N \times 1$） |
| 输出 | `H_x`（$(2N-p) \times m$）、`res`（$(2N-p) \times 1$），`H_f` 被消去 |
| 调用方 | `UpdaterMSCKF::update()`、`UpdaterSLAM::update()`、`UpdaterSLAM::delayed_init()` |

### 1.2 一句话概括

`nullspace_project_inplace` 用 Givens QR 分解对 $\mathbf{H}_f$ 做左零空间投影，消去特征位置变量 $\delta\mathbf{x}_f$，只保留关于状态 $\delta\mathbf{x}$ 的约束方程。

---

## 2. 数学模型

### 2.1 问题定义

给定线性化测量方程：

$$\mathbf{r} = \mathbf{H}_x\,\delta\mathbf{x} + \mathbf{H}_f\,\delta\mathbf{x}_f + \mathbf{n}, \quad \mathbf{n} \sim \mathcal{N}(\mathbf{0},\, \sigma^2\mathbf{I}) \tag{S13b-1}$$

其中 $\mathbf{H}_f \in \mathbb{R}^{2N \times p}$（$p=3$ 对于 3D 位置，$p=1$ 对于单逆深度），$\mathbf{H}_x \in \mathbb{R}^{2N \times m}$。

**目标**：消去 $\delta\mathbf{x}_f$，得到只关于 $\delta\mathbf{x}$ 的约束。

### 2.2 QR 分解

对 $\mathbf{H}_f$ 做 QR 分解（通过 Givens 旋转）：

$$\mathbf{Q}^\top \mathbf{H}_f = \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix} \tag{S13b-2}$$

其中 $\mathbf{T} \in \mathbb{R}^{p \times p}$ 上三角（满秩），$\mathbf{Q} = [\mathbf{Q}_1 \;\; \mathbf{Q}_2]$，$\mathbf{Q}_1 \in \mathbb{R}^{2N \times p}$，$\mathbf{Q}_2 \in \mathbb{R}^{2N \times (2N-p)}$。

### 2.3 投影后的系统

将 $\mathbf{Q}^\top$ 左乘方程 (S13b-1)：

$$\begin{bmatrix} \mathbf{Q}_1^\top\mathbf{r} \\ \mathbf{Q}_2^\top\mathbf{r} \end{bmatrix} = \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{H}_x \\ \mathbf{Q}_2^\top\mathbf{H}_x \end{bmatrix}\delta\mathbf{x} + \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix}\delta\mathbf{x}_f + \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{n} \\ \mathbf{Q}_2^\top\mathbf{n} \end{bmatrix} \tag{S13b-3}$$

**下半部分**（$2N-p$ 个方程）：

$$\mathbf{Q}_2^\top\mathbf{r} = \mathbf{Q}_2^\top\mathbf{H}_x\,\delta\mathbf{x} + \mathbf{Q}_2^\top\mathbf{n} \tag{S13b-4}$$

$\delta\mathbf{x}_f$ 被完全消去（$\mathbf{Q}_2^\top\mathbf{H}_f = \mathbf{0}$）。

### 2.4 噪声统计性质

由于 $\mathbf{R} = \sigma^2\mathbf{I}_{2N}$，正交变换保持噪声统计：

$$\mathbb{E}[\mathbf{Q}_2^\top\mathbf{n}\mathbf{n}^\top\mathbf{Q}_2] = \sigma^2\mathbf{Q}_2^\top\mathbf{Q}_2 = \sigma^2\mathbf{I}_{2N-p} \tag{S13b-5}$$

投影后的噪声仍是白噪声，方差不变。

### 2.5 从最小二乘角度理解

正交变换 $\mathbf{Q}^\top$ 不改变 2-范数：

$$\|\mathbf{r} - \mathbf{H}_x\delta\mathbf{x} - \mathbf{H}_f\delta\mathbf{x}_f\|^2 = \|\tilde{\mathbf{r}}_1 - \tilde{\mathbf{H}}_{x1}\delta\mathbf{x} - \mathbf{T}\delta\mathbf{x}_f\|^2 + \|\tilde{\mathbf{r}}_2 - \tilde{\mathbf{H}}_{x2}\delta\mathbf{x}\|^2 \tag{S13b-6}$$

前 $p$ 个自由度被 $\delta\mathbf{x}_f$ 吸收（$\mathbf{T}$ 可逆），后 $2N-p$ 个自由度才是关于 $\delta\mathbf{x}$ 的独立约束。

### 2.6 卡方检验：外点检测

零空间投影后，需要进行**卡方检验**来判断特征是否为外点（outlier）。

#### 2.6.1 残差的统计分布

投影后的系统为：

$$\tilde{\mathbf{r}}_2 = \tilde{\mathbf{H}}_{x2}\,\delta\mathbf{x} + \tilde{\mathbf{n}}_2 \tag{S13b-7}$$

如果模型正确且没有外点，残差 $\tilde{\mathbf{r}}_2$ 应该服从零均值高斯分布：

$$\tilde{\mathbf{r}}_2 \sim \mathcal{N}(\mathbf{0},\, \mathbf{S}) \tag{S13b-8}$$

其中 $\mathbf{S}$ 是残差的协方差矩阵。

#### 2.6.2 残差协方差 S 的计算

$$\mathbf{S} = \tilde{\mathbf{H}}_{x2}\, \mathbf{P}_{\text{marg}}\, \tilde{\mathbf{H}}_{x2}^\top + \mathbf{R} \tag{S13b-9}$$

其中：
- $\tilde{\mathbf{H}}_{x2}$：投影后的状态 Jacobian（$(2N-p) \times m$）
- $\mathbf{P}_{\text{marg}}$：相关状态的边缘协方差（$m \times m$）
- $\mathbf{R} = \sigma_{\text{pix}}^2 \mathbf{I}_{2N-p}$：测量噪声协方差

**物理意义**：

| 项 | 含义 |
|---|------|
| $\tilde{\mathbf{H}}_{x2} \mathbf{P}_{\text{marg}} \tilde{\mathbf{H}}_{x2}^\top$ | 状态不确定性通过测量模型传播到残差空间 |
| $\mathbf{R}$ | 传感器本身的测量噪声 |
| $\mathbf{S}$ | **残差的总不确定性**（状态不确定 + 测量噪声） |

#### 2.6.3 Mahalanobis 距离与卡方分布

**Mahalanobis 距离**定义为：

$$d^2 = \tilde{\mathbf{r}}_2^\top \mathbf{S}^{-1} \tilde{\mathbf{r}}_2 \tag{S13b-10}$$

如果模型正确，$d^2$ 服从自由度为 $2N-p$ 的 $\chi^2$ 分布：

$$d^2 \sim \chi^2(2N-p) \tag{S13b-11}$$

#### 2.6.4 外点检测原理

**正常情况**：如果特征投影正确，残差主要来自：
1. 状态估计误差（小的，因为 EKF 一直在更新）
2. 测量噪声（小的，$\sigma_{\text{pix}} \approx 1$ 像素）

此时 $d^2$ 应该较小，落在 $\chi^2(2N-p)$ 分布的 95% 分位数以内。

**外点情况**：如果特征是外点（误匹配、遮挡等），残差会异常大：

$$\|\tilde{\mathbf{r}}_2\| \gg \sqrt{\text{diag}(\mathbf{S})}$$

此时 $d^2$ 会远超 $\chi^2$ 阈值，被判定为外点。

#### 2.6.5 代码实现

```cpp
// UpdaterMSCKF.cpp:209-234
Eigen::MatrixXd P_marg = StateHelper::get_marginal_covariance(state, Hx_order);
Eigen::MatrixXd S = H_x * P_marg * H_x.transpose();
S.diagonal() += _options.sigma_pix_sq * Eigen::VectorXd::Ones(S.rows());
double chi2 = res.dot(S.llt().solve(res));

// 获取阈值（95% 分位数）
double chi2_check = chi_squared_table[res.rows()];

// 外点检测
if (chi2 > _options.chi2_multipler * chi2_check) {
    (*it2)->to_delete = true;  // 外点，删除
    it2 = feature_vec.erase(it2);
    continue;
}
```

#### 2.6.6 阈值的确定

- `chi_squared_table[n]`：$\chi^2(n)$ 分布的 95% 分位数（预计算表）
- `chi2_multipler`：可调节的倍率（默认 1.0），增大则更宽松

**示例**：
- $n = 4$（2 个观测，$p=0$）：$\chi^2_{0.95}(4) = 9.49$
- $n = 6$（3 个观测，$p=0$）：$\chi^2_{0.95}(6) = 12.59$
- $n = 2N-p$（一般情况）：查表或计算

#### 2.6.7 完整流程图

```
特征观测 → 计算残差 r
              ↓
         计算 S = H P H^T + R
              ↓
    计算 χ² = r^T S^{-1} r
              ↓
      χ² > 阈值 × 倍率？
         ↙          ↘
       是            否
       ↓             ↓
    外点，删除    正常，保留
```

---

## 3. 代码实现：逐行对照理论

### 3.1 代码全文

```cpp
// UpdaterHelper.cpp:426-454
void UpdaterHelper::nullspace_project_inplace(H_f, H_x, res) {
    Eigen::JacobiRotation<double> tempHo_GR;
    for (int n = 0; n < H_f.cols(); ++n) {                    // 外层：逐列
        for (int m = (int)H_f.rows() - 1; m > n; m--) {      // 内层：从下到上
            tempHo_GR.makeGivens(H_f(m-1,n), H_f(m,n));      // 构造 Givens 旋转
            (H_f.block(m-1,n,2,H_f.cols()-n)).applyOnTheLeft(0,1,tempHo_GR.adjoint());
            (H_x.block(m-1,0,2,H_x.cols())).applyOnTheLeft(0,1,tempHo_GR.adjoint());
            (res.block(m-1,0,2,1)).applyOnTheLeft(0,1,tempHo_GR.adjoint());
        }
    }
    H_x = H_x.block(H_f.cols(), 0, H_x.rows()-H_f.cols(), H_x.cols()).eval();
    res = res.block(H_f.cols(), 0, res.rows()-H_f.cols(), 1).eval();
    assert(H_x.rows() == res.rows());
}
```

### 3.2 逐行代码 ↔ 理论对照

#### 3.2.1 双循环结构：逐列消元

```
H_f 初始形态（以 p=3, 2N=8 为例）:        目标形态（上三角 T）:

[ ×  ×  × ]                              [ ×  ×  × ]
[ ×  ×  × ]                              [ 0  ×  × ]
[ ×  ×  × ]    → Givens QR →            [ 0  0  × ]
[ ×  ×  × ]                              [ 0  0  0 ]
[ ×  ×  × ]                              [ 0  0  0 ]
[ ×  ×  × ]                              [ 0  0  0 ]
[ ×  ×  × ]                              [ 0  0  0 ]
[ ×  ×  × ]                              [ 0  0  0 ]
```

**外层循环** `n = 0, 1, 2`（共 $p=3$ 列）：逐列处理。

**内层循环** `m = 7, 6, ..., n+1`（从最后一行到第 $n+1$ 行）：对当前列，从下到上逐个消去。

每次迭代用 Givens 旋转消去 `H_f(m, n)` 这个元素，使之为零。

#### 3.2.2 Givens 旋转的构造

```cpp
tempHo_GR.makeGivens(H_f(m-1,n), H_f(m,n));
```

**理论对应**：公式 (MT-5)(MT-6)。给定 $(a, b) = (H_f(m-1,n), H_f(m,n))$，构造 2×2 正交矩阵 $\mathbf{G}$ 使得：

$$\mathbf{G}^\top \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix}, \quad r = \sqrt{a^2+b^2}$$

`makeGivens` 内部计算 $c = a/r$, $s = b/r$，构造 $\mathbf{G} = \begin{bmatrix} c & s \\ -s & c \end{bmatrix}$。

#### 3.2.3 同步变换三个矩阵

```cpp
(H_f.block(m-1,n,2,H_f.cols()-n)).applyOnTheLeft(0,1, G.adjoint());
(H_x.block(m-1,0,2,H_x.cols())).applyOnTheLeft(0,1, G.adjoint());
(res.block(m-1,0,2,1)).applyOnTheLeft(0,1, G.adjoint());
```

**理论对应**：公式 (S13b-3)。$\mathbf{Q}^\top$ 左乘整个方程时，$\mathbf{H}_f$、$\mathbf{H}_x$、$\mathbf{r}$ 都要被同一个 $\mathbf{G}^\top$ 变换。

为什么同步变换？因为最终我们要的是 $\mathbf{Q}^\top\mathbf{H}_x$ 和 $\mathbf{Q}^\top\mathbf{r}$（公式 (S13b-3) 的下半部分）。不显式构造 $\mathbf{Q}$，而是将 $\mathbf{Q}^\top$ 的作用"累积"到 $\mathbf{H}_x$ 和 $\mathbf{r}$ 上。

**`applyOnTheLeft(0, 1, G.adjoint())`** 的含义：
- `block(m-1, ...)`：取第 $m-1$ 行和第 $m$ 行
- `applyOnTheLeft(0, 1, ...)`：用 $\mathbf{G}^\top$ 左乘这两行
- `G.adjoint()`：$\mathbf{G}$ 的共转置（实矩阵就是转置）

#### 3.2.4 为什么只作用于 `H_f.cols()-n` 列

```cpp
H_f.block(m-1, n, 2, H_f.cols() - n)  // 从第 n 列开始
```

**理论对应**：第 $n$ 列的前 $n$ 行已经被前面的迭代消成上三角了（第 $0$ 行到第 $n-1$ 行在第 $n$ 列上已经是零或非零的上三角元素）。从第 $n$ 列开始可以避免覆盖已经处理好的上三角部分。

这是一个**性能优化**，数学上作用于全部列结果一样（上三角部分会被 $\mathbf{G}^\top$ 变换但不影响最终结果）。

#### 3.2.5 取下半部分（左零空间）

```cpp
H_x = H_x.block(H_f.cols(), 0, H_x.rows()-H_f.cols(), H_x.cols()).eval();
res = res.block(H_f.cols(), 0, res.rows()-H_f.cols(), 1).eval();
```

**理论对应**：公式 (S13b-4)。QR 分解后：

$$\mathbf{Q}^\top \mathbf{H}_f = \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix} \quad \Rightarrow \quad \mathbf{Q} = [\mathbf{Q}_1 \;\; \mathbf{Q}_2]$$

$\mathbf{Q}_2^\top$ 对应 $\mathbf{Q}^\top$ 的下半部分（第 $p$ 行到第 $2N-1$ 行）。所以取 `H_x` 和 `res` 的第 $p$ 行到第 $2N-1$ 行，就是 $\mathbf{Q}_2^\top\mathbf{H}_x$ 和 $\mathbf{Q}_2^\top\mathbf{r}$。

**维度变化**：
- 输入：`H_x` 是 $2N \times m$，`res` 是 $2N \times 1$
- 输出：`H_x` 是 $(2N-p) \times m$，`res` 是 $(2N-p) \times 1$

#### 3.2.6 为什么用 `.eval()`

```cpp
H_x = H_x.block(...).eval();
```

**原因**：Eigen 的 `block()` 返回的是一个**视图**（view），不是新矩阵。直接赋值 `H_x = H_x.block(...)` 会导致**别名问题**（aliasing）——右边读取 `H_x` 的同时左边在写 `H_x`，数据会被覆盖。

`.eval()` 强制立即计算为一个临时矩阵，避免别名问题。

### 3.3 完整示例：$p=2, 2N=4$ 的消元过程

初始：

$$\mathbf{H}_f = \begin{bmatrix} a_{00} & a_{01} \\ a_{10} & a_{11} \\ a_{20} & a_{21} \\ a_{30} & a_{31} \end{bmatrix}$$

**第 1 列（n=0）**：

| 迭代 | 操作 | 效果 |
|------|------|------|
| m=3 | Givens 消去 $a_{30}$ | 第 3 行第 0 列 → 0 |
| m=2 | Givens 消去 $a_{20}$ | 第 2 行第 0 列 → 0 |
| m=1 | Givens 消去 $a_{10}$ | 第 1 行第 0 列 → 0 |

第 1 列完成后：

$$\mathbf{H}_f \to \begin{bmatrix} * & * \\ 0 & * \\ 0 & * \\ 0 & * \end{bmatrix}$$

**第 2 列（n=1）**：

| 迭代 | 操作 | 效果 |
|------|------|------|
| m=3 | Givens 消去 $a_{31}$ | 第 3 行第 1 列 → 0 |
| m=2 | Givens 消去 $a_{21}$ | 第 2 行第 1 列 → 0 |

第 2 列完成后：

$$\mathbf{H}_f \to \begin{bmatrix} * & * \\ 0 & * \\ 0 & 0 \\ 0 & 0 \end{bmatrix} = \begin{bmatrix} \mathbf{T} \\ \mathbf{0} \end{bmatrix}$$

取下半部分（第 2 行到第 3 行）：

$$\mathbf{Q}_2^\top\mathbf{H}_x = \begin{bmatrix} \text{行 2} \\ \text{行 3} \end{bmatrix}, \quad \mathbf{Q}_2^\top\mathbf{r} = \begin{bmatrix} r_2 \\ r_3 \end{bmatrix}$$

这就是公式 (S13b-4) 的最终结果。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{H}_f$ | 特征 Jacobian | `H_f` | `:426` |
| $\mathbf{H}_x$ | 状态 Jacobian | `H_x` | `:426` |
| $\mathbf{r}$ | 残差 | `res` | `:426` |
| $\mathbf{Q}$ | QR 正交矩阵（隐式） | Givens 旋转序列 | `:432` |
| $\mathbf{T}$ | $\mathbf{H}_f$ 的上三角部分 | `H_f.block(0,0,p,p)` | `:449` |
| $\mathbf{Q}_2^\top\mathbf{H}_x$ | 投影后状态 Jacobian | `H_x`（压缩后） | `:449` |
| $\mathbf{Q}_2^\top\mathbf{r}$ | 投影后残差 | `res`（压缩后） | `:450` |

---

## 5. 关键设计决策

### 5.1 Givens vs Householder

| 方法 | 适用场景 | 复杂度 |
|------|---------|--------|
| Givens | 稀疏矩阵、增量更新 | $O(m^2 n)$ |
| Householder | 稠密矩阵 | $O(m n^2)$ |

MSCKF 选择 Givens 是因为：
- 测量 Jacobian 通常是稀疏的（每个观测只涉及少量变量）
- 可以只作用于非零列范围（`H_f.cols()-n`）
- 与 `measurement_compress_inplace` 复用同一套逻辑

### 5.2 为什么"就地"（inplace）

函数名 `inplace` 表示直接修改输入参数，不创建新矩阵。这避免了额外的内存分配，对实时系统很重要。

---

## 6. 代码位置速查

| 功能 | 代码行 |
|------|--------|
| Givens QR 循环 | `:432-443` |
| 取左零空间 | `:449-450` |
| 断言检查 | `:453` |

---

## 7. 同类实现对比

| 特征维度 | OpenVINS | VINS-Mono | OKVIS | 常规 MSCKF |
|---------|---------|-----------|-------|-----------|
| 零空间投影 | Givens QR | 无 | 部分 | 解析 |
| 测量压缩 | Givens QR | 无 | Schur | 无 |

### 设计取舍分析

- 选择 **Givens QR 做投影 + 压缩** 是因为：数值稳定（Golub & Van Loan Algorithm 5.2.4），且投影与压缩可复用同一套旋转逻辑
- 选择 **inplace 修改** 是为了避免额外内存分配，适合实时系统
