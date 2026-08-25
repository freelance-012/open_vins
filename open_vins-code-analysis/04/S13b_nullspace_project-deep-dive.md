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

---

## 3. 代码实现

### 3.1 Givens QR 分解

代码（[UpdaterHelper.cpp:426-454](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L426-L454)）：

```cpp
void UpdaterHelper::nullspace_project_inplace(H_f, H_x, res) {
    Eigen::JacobiRotation<double> G;
    for (int n = 0; n < H_f.cols(); ++n) {                    // 逐列
        for (int m = H_f.rows()-1; m > n; m--) {              // 从下到上消去
            G.makeGivens(H_f(m-1,n), H_f(m,n));               // 计算 Givens 旋转
            H_f.block(m-1,n,2,H_f.cols()-n).applyOnTheLeft(0,1,G.adjoint());
            H_x.block(m-1,0,2,H_x.cols()).applyOnTheLeft(0,1,G.adjoint());
            res.block(m-1,0,2,1).applyOnTheLeft(0,1,G.adjoint());
        }
    }
    // 取下半部分（左零空间）
    H_x = H_x.block(H_f.cols(), 0, H_x.rows()-H_f.cols(), H_x.cols()).eval();
    res = res.block(H_f.cols(), 0, res.rows()-H_f.cols(), 1).eval();
}
```

### 3.2 Givens 旋转的工作原理

对矩阵 $\mathbf{A}$ 的第 $m-1$ 行和第 $m$ 行，构造 Givens 旋转 $\mathbf{G}$ 使得：

$$\begin{bmatrix} c & s \\ -s & c \end{bmatrix}^\top \begin{bmatrix} a_{m-1,n} \\ a_{m,n} \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix}$$

其中 $c = a_{m-1,n}/\sqrt{a_{m-1,n}^2 + a_{m,n}^2}$，$s = a_{m,n}/\sqrt{a_{m-1,n}^2 + a_{m,n}^2}$。

逐列从下到上消去，最终 $\mathbf{H}_f$ 变为上三角 $\mathbf{T}$。

### 3.3 为什么用 `eval()`

Eigen 的 lazy evaluation 可能导致 aliasing 问题（矩阵块引用自身）。`.eval()` 强制立即计算，避免读取被覆盖的数据。

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
