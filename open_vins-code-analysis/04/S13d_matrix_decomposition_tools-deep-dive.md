# SLAM 中常用矩阵分解工具深度解析

> **生成日期**：2026-08-19
> **前置文档**：`S13_UpdaterMSCKF-deep-dive.md`、`S13b_nullspace_project-deep-dive.md`、`S13c_measurement_compress-deep-dive.md`、`S12_StateHelper-initialize-deep-dive.md`、`S11_StateHelper-marginalize-deep-dive.md`
> **核心代码**：`ov_msckf/src/update/UpdaterHelper.cpp`、`ov_msckf/src/state/StateHelper.cpp`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Golub & Van Loan "Matrix Computations" 4th Ed. [48]；Hartley & Zisserman "Multiple View Geometry" [49]

---

## 概述：同一把刀，切不同的菜

OpenVINS 和 VINS-Mono 等 SLAM 系统大量使用矩阵分解，但**同一个数学工具在不同场景下目的完全不同**。理解这一点，才能看透代码背后的数学意图。

| 工具 | 数学本质 | SLAM 用途 A | SLAM 用途 B | SLAM 用途 C |
|------|---------|------------|------------|------------|
| **SVD** | $A = U\Sigma V^\top$ | 三角化求深度 | PnP 求解位姿 | 可观性/零空间分析 |
| **QR** | $A = QR$ | 零空间投影消变量 | 测量压缩减方程 | 协方差初始化 QR 分离 |
| **Givens** | 2×2 旋转消元 | QR 的具体实现 | QR 的具体实现 | — |
| **Schur 补** | $A/A_{22} = A_{11} - A_{12}A_{22}^{-1}A_{21}$ | BA 边缘化地图点 | EKF 边缘化旧状态 | 信息矩阵融合 |

---

## 1. SVD 分解：$A = U\Sigma V^\top$

### 1.1 数学本质

任意矩阵 $A \in \mathbb{R}^{m \times n}$ 都可以分解为：

$$A = U\Sigma V^\top \tag{MT-1}$$

其中 $U \in \mathbb{R}^{m \times m}$、$V \in \mathbb{R}^{n \times n}$ 正交，$\Sigma = \text{diag}(\sigma_1, \ldots, \sigma_r, 0, \ldots, 0)$（$\sigma_1 \geq \cdots \geq \sigma_r > 0$）。

**核心性质**：
- $\sigma_i$ 是 $A^\top A$ 的特征值的平方根
- $V$ 的最后一列（对应最小奇异值）是 $\min\|Ax\|$ 的解（当 $Ax \approx 0$ 时）
- $V$ 的前 $r$ 列张成 $A$ 的列空间，后 $n-r$ 列张成 $A$ 的零空间

### 1.2 用途 A：三角化求特征深度

**场景**：已知特征在两个相机帧中的投影 $\mathbf{u}_1, \mathbf{u}_2$ 和两帧之间的相对位姿 $(R, t)$，求特征的 3D 位置 $\mathbf{p}_F$。

**数学**：每个观测给出方程 $\mathbf{u} \times \mathbf{p}_F^{C} = \mathbf{0}$，即：

$$\begin{bmatrix} 0 & -1 & v \\ 1 & 0 & -u \\ -v & u & 0 \end{bmatrix} \begin{bmatrix} R & t \end{bmatrix} \begin{bmatrix} \mathbf{p}_F \\ 1 \end{bmatrix} = \mathbf{0}$$

两个观测堆叠后得到 $A\mathbf{x} = \mathbf{0}$（$6 \times 4$ 矩阵），用 SVD 求最小奇异值对应的右奇异向量作为解。

**OpenVINS 代码**：`FeatureInitializer::single_triangulation()` 中用 `JacobiSVD` 求解：

```cpp
Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
Eigen::Vector4d x = svd.matrixV().col(3);  // 最小奇异值对应的右奇异向量
p_F = x.head(3) / x(3);  // 齐次坐标转 3D
```

### 1.3 用途 B：PnP 求解相机位姿

**场景**：已知 $N$ 个 3D-2D 点对应，求相机位姿 $(R, t)$。

**数学**：每个点对应给出 $\mathbf{u} \times (R\mathbf{p}_F + t) = \mathbf{0}$，线性化后得到 $A\mathbf{x} = \mathbf{b}$（$\mathbf{x}$ 包含 $R$ 的 9 个元素和 $t$ 的 3 个元素）。用 SVD 求最小二乘解，再从 $R$ 部分投影回 $SO(3)$。

**VINS-Mono 代码**：`initial_sfm.cpp` 中的 `solveRelativeRT()` 用 SVD 分解本质矩阵 $E = U\Sigma V^\top$，从 $U, V$ 恢复 $R, t$：

```cpp
Eigen::JacobiSVD<Eigen::MatrixXd> svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
Eigen::Matrix3d U = svd.matrixU();
Eigen::Matrix3d V = svd.matrixV();
// 从 U, V 恢复 R, t（四种可能解，选深度为正的那个）
```

### 1.4 用途 C：可观性/零空间分析

**场景**：分析 VIO 系统的可观性，确定哪些状态不可观。

**数学**：对可观性矩阵 $\mathcal{O}$ 做 SVD，零奇异值对应的右奇异向量就是不可观方向。

**OpenVINS 论文**：证明 VIO 有 4 个不可观自由度（yaw + 3D 全局位置），通过对 $\mathcal{O}$ 的零空间分析得出。

**VINS-Mono 代码**：`margin_factor.cpp` 中对先验信息矩阵做 SVD，检查零空间维度。

---

## 2. QR 分解：$A = QR$

### 2.1 数学本质

任意矩阵 $A \in \mathbb{R}^{m \times n}$（$m \geq n$）都可以分解为：

$$A = QR \tag{MT-2}$$

其中 $Q \in \mathbb{R}^{m \times m}$ 正交，$R \in \mathbb{R}^{m \times n}$ 上三角（下半部分为零）。

**核心性质**：
- $Q$ 正交 → $\|Qx\| = \|x\|$（保范数）
- $R$ 上三角 → 回代求解 $Rx = b$ 只需 $O(n^2)$
- $A^\top A = R^\top Q^\top Q R = R^\top R$（Cholesky 分解的替代）

### 2.2 用途 A：零空间投影（消去特征变量）

**场景**：MSCKF 中，特征位置 $\delta\mathbf{x}_f$ 不应加入状态，需要从测量方程中消去。

**数学**：对 $\mathbf{H}_f$ 做 QR 分解 $\mathbf{Q}^\top\mathbf{H}_f = [\mathbf{T}^\top, \mathbf{0}^\top]^\top$，将 $\mathbf{Q}^\top$ 左乘整个方程。下半部分（$\mathbf{Q}_2^\top$ 行）中 $\mathbf{H}_f$ 变为零，$\delta\mathbf{x}_f$ 被消去：

$$\mathbf{Q}_2^\top\mathbf{r} = \mathbf{Q}_2^\top\mathbf{H}_x\,\delta\mathbf{x} + \mathbf{Q}_2^\top\mathbf{n} \tag{MT-3}$$

**OpenVINS 代码**：`nullspace_project_inplace()`（S13b），用 Givens QR 实现。

```cpp
// 对 H_f 逐列 Givens 消元，同步变换 H_x 和 res
// 最后取 H_x 和 res 的下半部分（左零空间）
H_x = H_x.block(p, 0, 2N-p, m);
res = res.block(p, 0, 2N-p, 1);
```

### 2.3 用途 B：测量压缩（减少方程数）

**场景**：多个特征拼接后 $H_x$ 有 1000 行、150 列，EKFUpdate 中要对 $1000 \times 1000$ 矩阵求逆，计算量太大。

**数学**：对 $H_x$ 做 QR 分解，取前 $N$ 行（上三角 $R$），最小二乘解不变：

$$H_x^\top H_x = R^\top R, \quad H_x^\top r = R^\top r_1 \tag{MT-4}$$

**OpenVINS 代码**：`measurement_compress_inplace()`（S13c），同样用 Givens QR。

```cpp
// 对 H_x 逐列 Givens 消元，同步变换 res
// 最后取前 N 行（上三角 R 和对应的 r_1）
H_x.conservativeResize(N, N);
res.conservativeResize(N, 1);
```

### 2.4 用途 C：协方差初始化中的 QR 分离

**场景**：`StateHelper::initialize()` 中，新变量（SLAM 特征）首次加入状态时，用 QR 将测量系统分离为"可逆初始化系统"和"零空间投影更新系统"。

**数学**：对 $\mathbf{H}_L$（新变量 Jacobian）做 QR 分解，上半部分（$p$ 行）用于直接求解新变量，下半部分（$k-p$ 行）用于更新旧状态（S12）。

**OpenVINS 代码**：`StateHelper::initialize()`（S12），用 `JacobiRotation`（Givens）实现。

```cpp
// Givens QR 分离 H_L
// 上半部分：H_finit * δx_f = Hxinit * δx + resinit → 直接解新变量
// 下半部分：Hup * δx = resup → EKFUpdate 更新旧状态
```

---

## 3. Givens 旋转：2×2 平面旋转

### 3.1 数学本质

Givens 旋转是一个 2×2 正交矩阵：

$$G = \begin{bmatrix} c & s \\ -s & c \end{bmatrix}, \quad c^2 + s^2 = 1 \tag{MT-5}$$

给定 $(a, b)$，选择 $c = a/\sqrt{a^2+b^2}$，$s = b/\sqrt{a^2+b^2}$，使得：

$$G^\top \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix}, \quad r = \sqrt{a^2 + b^2} \tag{MT-6}$$

### 3.2 用途 A：QR 分解的实现

**场景**：上述 QR 分解的具体实现方式。

**与 Householder 的对比**：

| 方法 | 操作单位 | 适用场景 | 复杂度 |
|------|---------|---------|--------|
| Givens | 2×2 旋转 | 稀疏矩阵、只需消去少量元素 | $O(m^2 n)$ |
| Householder | 反射向量 | 稠密矩阵、需消去整列 | $O(m n^2)$ |

**OpenVINS 选择 Givens 的原因**：
1. 测量 Jacobian 是稀疏的（每个观测只涉及少量变量）
2. 可以只作用于非零列范围（`H_f.cols()-n`），避免无效计算
3. 与 `measurement_compress` 复用同一套逻辑

### 3.3 用途 B：增量 QR 更新

**场景**：当新观测到来时，只需在现有 QR 分解基础上追加一行，用 Givens 旋转消去新增元素，无需重新分解整个矩阵。

**VINS-Mono 中的应用**：滑动窗口优化中，新帧加入时增量更新 QR 分解。

**OpenVINS 中的体现**：虽然 OpenVINS 是 EKF 架构不做 BA，但 `nullspace_project_inplace` 和 `measurement_compress_inplace` 本质上都是对"当前批次的特征"做 QR，每帧独立处理。

---

## 4. Schur 补：$A/A_{22} = A_{11} - A_{12}A_{22}^{-1}A_{21}$

### 4.1 数学本质

对分块矩阵 $A = \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}$，若 $A_{22}$ 可逆，则 Schur 补定义为：

$$A/A_{22} = A_{11} - A_{12}A_{22}^{-1}A_{21} \tag{MT-7}$$

**核心性质**：
- $\det(A) = \det(A_{22}) \cdot \det(A/A_{22})$
- 若 $A$ 正定，则 $A/A_{22}$ 也正定
- 在信息矩阵（Hessian）中，Schur 补等价于"边缘化"（marginalize out）$\mathbf{x}_2$

### 4.2 用途 A：BA 中边缘化地图点加速求解

**场景**：Bundle Adjustment 中，状态向量 $\mathbf{x} = [\mathbf{x}_c^\top, \mathbf{x}_p^\top]^\top$（相机位姿 + 地图点），Hessian 矩阵为：

$$H = \begin{bmatrix} H_{cc} & H_{cp} \\ H_{pc} & H_{pp} \end{bmatrix} \tag{MT-8}$$

其中 $H_{pp}$ 是块对角的（每个地图点只与观测它的相机相连），求逆很快。

**Schur 补技巧**：先消去地图点：

$$H_{\text{reduced}} = H_{cc} - H_{cp}H_{pp}^{-1}H_{pc} \tag{MT-9}$$

$H_{\text{reduced}}$ 只包含相机位姿，维度远小于原 Hessian。求解 $\delta\mathbf{x}_c$ 后，回代求 $\delta\mathbf{x}_p$：

$$\delta\mathbf{x}_p = H_{pp}^{-1}(b_p - H_{pc}\,\delta\mathbf{x}_c) \tag{MT-10}$$

**VINS-Mono 代码**：`factor.cpp` 中的 `SchurComplement()` 函数。

**OpenVINS 中的体现**：EKF 架构不做 BA，但 `measurement_compress` 本质上也是"减少方程数"的类似思想。

### 4.3 用途 B：EKF 中边缘化旧状态

**场景**：滑动窗口 EKF 中，最老的 clone 位姿需要被边缘化，但希望保留它对剩余状态的约束信息。

**数学**：设状态 $\mathbf{x} = [\mathbf{x}_{\text{old}}^\top, \mathbf{x}_{\text{rest}}^\top]^\top$，信息矩阵 $\Lambda = P^{-1}$ 分块后，边缘化 $\mathbf{x}_{\text{old}}$：

$$\Lambda_{\text{rest}}' = \Lambda_{\text{rest}} - \Lambda_{\text{rest},\text{old}}\Lambda_{\text{old}}^{-1}\Lambda_{\text{old},\text{rest}} \tag{MT-11}$$

**MSCKF 的特殊性**：MSCKF 的特征约束已经通过零空间投影融入了 IMU 状态，被边缘化 clone 上的特征信息已经"转移"到了剩余状态中。因此 MSCKF 可以直接删除（不做 Schur 补），见 S11。

**VINS-Mono 的做法**：VINS-Mono 用 Schur 补做边缘化先验（marginalization prior），保留被删变量的约束信息。

### 4.4 用途 C：信息矩阵融合（因子图）

**场景**：因子图中，多个因子对同一变量的约束可以合并。

**数学**：若因子 1 给出信息矩阵 $\Lambda_1$，因子 2 给出 $\Lambda_2$，则总信息矩阵为 $\Lambda = \Lambda_1 + \Lambda_2$。这等价于对两个因子的残差做 Schur 补融合。

**OpenVINS 中的体现**：MSCKF 和 SLAM 特征的更新本质上是依次将信息添加到协方差矩阵中（通过 EKFUpdate 的 $P' = P - KSK^\top$）。

---

## 5. 四种工具的对比总结

### 5.1 数学关系

```
SVD:  A = UΣV^T          ← 最通用，但最慢
  ↓ 取 R = ΣV^T
QR:   A = QR             ← 保范数，用于消元
  ↓ 取 R 的上三角部分
Givens: 2×2 旋转消元      ← QR 的具体实现，适合稀疏/增量

Schur 补: A/A22 = A11 - A12 A22^{-1} A21   ← 用于边缘化/降维
  ↓ 在信息矩阵（Hessian/协方差逆）中操作
  ≈ QR 的"对偶"操作（QR 在测量矩阵上，Schur 在信息矩阵上）
```

### 5.2 在 OpenVINS 中的使用位置

| 工具 | 函数 | 文件 | 目的 |
|------|------|------|------|
| SVD | `single_triangulation()` | `FeatureInitializer.cpp` | 三角化求深度 |
| Givens QR | `nullspace_project_inplace()` | `UpdaterHelper.cpp:426` | 消去特征变量 |
| Givens QR | `measurement_compress_inplace()` | `UpdaterHelper.cpp:456` | 压缩方程数 |
| Givens QR | `initialize()` | `StateHelper.cpp:429` | QR 分离初始化/更新系统 |
| 直接删除 | `marginalize()` | `StateHelper.cpp:271` | 边缘化（MSCKF 特殊） |
| LLT | `EKFPropagation()` | `StateHelper.cpp:88` | 协方差对称正定求解 |
| LLT | `EKFUpdate()` | `StateHelper.cpp:161` | $S^{-1}$ 求解 |

### 5.3 在 VINS-Mono 中的对比

| 工具 | VINS-Mono 用途 | OpenVINS 用途 |
|------|---------------|--------------|
| SVD | 三角化、PnP、本质矩阵分解 | 三角化 |
| QR | 增量 BA 更新 | 零空间投影、测量压缩 |
| Schur 补 | BA 边缘化地图点、边缘化先验 | 不使用（MSCKF 直接删除） |
| Givens | — | QR 的具体实现 |

---

## 6. 设计取舍分析

### 6.1 为什么 MSCKF 不用 Schur 补做边缘化

MSCKF 的零空间投影已经将特征观测的约束**完全融入**了 IMU 状态（公式 (MT-3)）。被边缘化 clone 上的特征信息已经"转移"到了剩余 clone 的协方差中。因此直接删除不会丢失信息。

**对比**：VINS-Mono 的 BA 优化中，地图点的约束只通过因子连接，边缘化时必须用 Schur 补保留约束信息（否则丢失）。

### 6.2 为什么用 Givens 而非 Householder

| 因素 | Givens | Householder |
|------|--------|-------------|
| 稀疏性利用 | ✅ 只作用非零列 | ❌ 整列反射 |
| 增量更新 | ✅ 追加行只需新旋转 | ❌ 需重新分解 |
| 代码复用 | ✅ 投影和压缩同一套逻辑 | — |
| 数值稳定性 | ✅ 同样稳定 | ✅ 同样稳定 |
| 稠密矩阵速度 | 较慢 | 较快 |

OpenVINS 的测量 Jacobian 是稀疏的（每个观测只涉及 1-2 个 clone），Givens 更合适。

### 6.3 为什么压缩后才构造 R

测量压缩改变了残差维度（从 $M$ 行变为 $N$ 行）。噪声协方差 $\mathbf{R}$ 的维度必须与残差匹配。如果在压缩前构造 $\mathbf{R}$（$M \times M$），压缩后就不对了。

---

## 7. 一句话总结

| 工具 | 一句话 |
|------|--------|
| **SVD** | "找到矩阵的'骨架'——哪些方向有信息，哪些方向是零" |
| **QR** | "把方程'折叠'起来——多余方程不提供新信息" |
| **Givens** | "QR 的螺丝刀——一次拧一个螺丝（消一个元素）" |
| **Schur 补** | "在信息层面'压缩'状态——保留约束，去掉变量" |
