# UpdaterHelper::measurement_compress_inplace 深度精读 (S13c)

> **生成日期**：2026-08-19
> **前置文档**：`S13_UpdaterMSCKF-deep-dive.md`（主文档）、`S13b_nullspace_project-deep-dive.md`
> **核心代码**：`ov_msckf/src/update/UpdaterHelper.cpp:456-487`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Golub & Van Loan "Matrix Computations" 4th Ed. [48] Algorithm 5.2.4

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
UpdaterMSCKF::update()
        │
        ├─ 逐特征处理 → 拼成 Hx_big, res_big
        │
        ▼
measurement_compress_inplace(Hx_big, res_big)  @ :456  ★ 本篇
        │
        ├─ Givens QR 分解 Hx_big
        └─ 压缩到 min(行数, 列数) 行
        │
        ▼
EKFUpdate(state, Hx_compressed, res_compressed, R)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | MSCKF 更新中所有特征拼接后的压缩 |
| 输入 | `H_x`（$M \times N$）、`res`（$M \times 1$） |
| 输出 | `H_x`（$\min(M,N) \times N$）、`res`（$\min(M,N) \times 1$） |
| 调用方 | `UpdaterMSCKF::update()` |

### 1.2 一句话概括

`measurement_compress_inplace` 用 Givens QR 将超定系统 $H_x\,\delta x \approx res$（$M > N$）压缩为方阵系统（$N \times N$），不改变最小二乘解和 Hessian。

---

## 2. 数学模型

### 2.1 问题定义

多个特征拼接后，$H_x$ 是 $M \times N$ 矩阵（$M > N$，超定）：

$$\min_{\delta\mathbf{x}} \|H_x\,\delta\mathbf{x} - res\|^2 \tag{S13c-1}$$

### 2.2 QR 分解

对 $H_x$ 做 QR 分解：

$$H_x = Q \begin{bmatrix} R \\ 0 \end{bmatrix} \tag{S13c-2}$$

其中 $Q$ 是 $M \times M$ 正交矩阵，$R$ 是 $N \times N$ 上三角。

### 2.3 压缩后的等价问题

因为 $Q$ 正交（$Q^\top Q = I$），目标函数不变：

$$\|H_x\,\delta\mathbf{x} - res\|^2 = \|Q^\top(H_x\,\delta\mathbf{x} - res)\|^2 = \left\|\begin{bmatrix} R \\ 0 \end{bmatrix}\delta\mathbf{x} - Q^\top res\right\|^2 \tag{S13c-3}$$

令 $Q^\top res = \begin{bmatrix} r_1 \\ r_2 \end{bmatrix}$（$r_1$ 是 $N \times 1$，$r_2$ 是 $(M-N) \times 1$）：

$$\|H_x\,\delta\mathbf{x} - res\|^2 = \|R\,\delta\mathbf{x} - r_1\|^2 + \underbrace{\|r_2\|^2}_{\text{常数，与 } \delta\mathbf{x} \text{ 无关}} \tag{S13c-4}$$

所以 $\min \|H_x\,\delta\mathbf{x} - res\|^2$ **完全等价于** $\min \|R\,\delta\mathbf{x} - r_1\|^2$。

### 2.4 Hessian 不变性

$$H_x^\top H_x = \begin{bmatrix} R^\top & 0 \end{bmatrix} Q^\top Q \begin{bmatrix} R \\ 0 \end{bmatrix} = R^\top R \tag{S13c-5}$$

**压缩前后 Hessian 完全相同！** QR 分解只是把 $M$ 个方程的"信息"浓缩到了 $N$ 个方程中，没有丢失任何关于 $\delta\mathbf{x}$ 的信息。

### 2.5 为什么可以压缩

**核心原理**：超定系统中，多余方程不提供新信息。当 $M > N$ 时，$H_x$ 的行空间维度最多为 $N$，多余 $M-N$ 行是前 $N$ 行的线性组合（加噪声）。QR 分解把冗余信息吸收到了 $R$ 中。

---

## 3. 代码实现

### 3.1 提前返回

代码（[UpdaterHelper.cpp:456-460](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L456-L460)）：

```cpp
void UpdaterHelper::measurement_compress_inplace(H_x, res) {
    // 如果行数 ≤ 列数，无需压缩
    if (H_x.rows() <= H_x.cols()) return;
```

### 3.2 Givens QR 压缩

代码（[UpdaterHelper.cpp:462-477](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L462-L477)）：

```cpp
    Eigen::JacobiRotation<double> G;
    for (int n = 0; n < H_x.cols(); n++) {                    // 逐列
        for (int m = H_x.rows()-1; m > n; m--) {              // 从下到上消去
            G.makeGivens(H_x(m-1,n), H_x(m,n));               // 计算 Givens 旋转
            H_x.block(m-1,n,2,H_x.cols()-n).applyOnTheLeft(0,1,G.adjoint());
            res.block(m-1,0,2,1).applyOnTheLeft(0,1,G.adjoint());
        }
    }
```

### 3.3 取压缩后的矩阵

代码（[UpdaterHelper.cpp:479-486](src/open_vins/ov_msckf/src/update/UpdaterHelper.cpp#L479-L486)）：

```cpp
    int r = std::min(H_x.rows(), H_x.cols());  // = H_x.cols() (因为 rows > cols)
    H_x.conservativeResize(r, H_x.cols());      // 取前 r 行（上三角 R）
    res.conservativeResize(r, res.cols());      // 取前 r 行（r_1）
```

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $H_x$ | 原始 Jacobian（$M \times N$） | `H_x`（输入） | `:456` |
| $res$ | 原始残差（$M \times 1$） | `res`（输入） | `:456` |
| $R$ | 压缩后 Jacobian（$N \times N$） | `H_x`（输出） | `:485` |
| $r_1$ | 压缩后残差（$N \times 1$） | `res`（输出） | `:486` |
| $Q$ | QR 正交矩阵（隐式） | Givens 旋转序列 | `:464` |

---

## 5. 关键设计决策

### 5.1 压缩前后 Hessian 对比

| 阶段 | Hessian | 维度 |
|------|---------|------|
| 压缩前 | $H_x^\top H_x$ | $N \times N$ |
| 压缩后 | $R^\top R$ | $N \times N$ |

**完全相同**（公式 (S13c-5)）。

### 5.2 为什么在 EKFUpdate 之前压缩

EKFUpdate 中要计算 $S = H P H^\top + R$ 并求逆：
- 压缩前：$M \times M$ 矩阵求逆 → $O(M^3)$
- 压缩后：$N \times N$ 矩阵求逆 → $O(N^3)$

当 $M = 1000$、$N = 150$ 时，计算量从 $10^9$ 降到 $3.4 \times 10^6$，减少 **99.7%**。

### 5.3 与 nullspace_project 的关系

| 函数 | 输入 | 输出 | 目的 |
|------|------|------|------|
| `nullspace_project_inplace` | $H_f$（$2N \times p$） | 消去 $p$ 列 | 消去特征变量 |
| `measurement_compress_inplace` | $H_x$（$M \times N$，$M > N$） | 压缩到 $N$ 行 | 减少 EKFUpdate 计算量 |

两者都用 Givens QR，但目的不同：
- **nullspace_project**：消去列（特征变量）
- **measurement_compress**：消去行（冗余方程）

---

## 6. 代码位置速查

| 功能 | 代码行 |
|------|--------|
| 提前返回 | `:459-460` |
| Givens QR 循环 | `:464-476` |
| 取压缩后的矩阵 | `:481-486` |

---

## 7. 同类实现对比

| 特征维度 | OpenVINS | VINS-Mono | OKVIS |
|---------|---------|-----------|-------|
| 测量压缩 | Givens QR | 无 | Schur 补 |
| 压缩时机 | EKFUpdate 前 | — | 优化前 |

### 设计取舍分析

- 选择 **Givens QR 压缩** 是因为：与 `nullspace_project` 复用同一套逻辑，代码简洁
- 选择 **inplace 修改** 是为了避免额外内存分配
- 选择 **压缩后才构造 R** 是因为：压缩改变了残差维度，R 的维度要匹配
