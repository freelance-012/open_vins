# compute_dongsi_coeff 重力约束求解深度解析 (F2a)

> **生成日期**：2026-08-25
> **前置文档**：`F1b_DynamicInitializer-algorithm-deep-dive.md`（动态初始化算法）
> **核心代码**：`ov_init/src/utils/helper.h:183-319`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Dong-Si & Mourikis "Estimator initialization in vision-aided inertial navigation with unknown camera-IMU calibration" IROS 2012; Gantmacher "The Theory of Matrices" (伴随矩阵理论)

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
DynamicInitializer::initialize()  (F1b)
        │  线性系统 Ax = b
        ▼
InitializerHelper::compute_dongsi_coeff(D, d, g)  ★ 本篇
   ├─ 输入: D (3×3 矩阵), d (3×1 向量), g (重力模长 9.81)
   ├─ 输出: coeff (7×1 多项式系数，最高次到常数项)
   ─ 用途: 构造伴随矩阵 → 特征值分解 → 求重力方向
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | 动态初始化线性系统求解后，重力约束处理 |
| **输入** | $\mathbf{D}$（3×3 对称矩阵）、$\mathbf{d}$（3×1 向量）、$g$（重力模长） |
| **输出** | $\text{coeff}$（7×1 多项式系数，$[\lambda^6, \lambda^5, \ldots, \lambda^0]^\top$） |
| **调用方** | `DynamicInitializer::initialize()`（F1b） |

### 1.2 一句话概括

`compute_dongsi_coeff` 计算带重力模长约束 $\|\mathbf{g}\|=g$ 的拉格朗日乘子法导出的 6 次多项式系数，用于后续伴随矩阵特征值分解求解重力方向。

---

## 2. 数学模型

### 2.1 问题来源

动态初始化线性系统求解后（见 F1b 式 (F1b-6)），得到关于重力 $\mathbf{g}^{I_0}$ 的约束方程：

$$(\mathbf{D} - \lambda\mathbf{I})\mathbf{g} = \mathbf{d} \tag{F2a-1}$$

其中 $\mathbf{D} \in \mathbb{R}^{3\times 3}$ 对称，$\mathbf{d} \in \mathbb{R}^3$，$\lambda$ 是拉格朗日乘子。

重力模长约束：

$$\|\mathbf{g}\|^2 = g^2 \tag{F2a-2}$$

### 2.2 拉格朗日乘子法

构造拉格朗日函数：

$$\mathcal{L}(\mathbf{g}, \lambda) = \|(\mathbf{D} - \lambda\mathbf{I})\mathbf{g} - \mathbf{d}\|^2 + \mu(\|\mathbf{g}\|^2 - g^2) \tag{F2a-3}$$

对 $\mathbf{g}$ 求导并令为零：

$$\frac{\partial \mathcal{L}}{\partial \mathbf{g}} = 2(\mathbf{D} - \lambda\mathbf{I})^\top[(\mathbf{D} - \lambda\mathbf{I})\mathbf{g} - \mathbf{d}] + 2\mu\mathbf{g} = \mathbf{0} \tag{F2a-4}$$

由于 $\mathbf{D}$ 对称，$(\mathbf{D} - \lambda\mathbf{I})^\top = (\mathbf{D} - \lambda\mathbf{I})$，故：

$$(\mathbf{D} - \lambda\mathbf{I})^2\mathbf{g} - (\mathbf{D} - \lambda\mathbf{I})\mathbf{d} + \mu\mathbf{g} = \mathbf{0} \tag{F2a-5}$$

整理得：

$$[(\mathbf{D} - \lambda\mathbf{I})^2 + \mu\mathbf{I}]\mathbf{g} = (\mathbf{D} - \lambda\mathbf{I})\mathbf{d} \tag{F2a-6}$$

### 2.3 消去 $\mu$

由式 (F2a-1)：$\mathbf{g} = (\mathbf{D} - \lambda\mathbf{I})^{-1}\mathbf{d}$（假设可逆）。

代入式 (F2a-2)：

$$\mathbf{d}^\top(\mathbf{D} - \lambda\mathbf{I})^{-2}\mathbf{d} = g^2 \tag{F2a-7}$$

这是关于 $\lambda$ 的有理方程。通分后得到多项式方程。

### 2.4 多项式系数推导

设 $\mathbf{D} - \lambda\mathbf{I} = \mathbf{M}(\lambda)$，则：

$$\mathbf{M}(\lambda) = \begin{bmatrix} D_{11}-\lambda & D_{12} & D_{13} \\ D_{21} & D_{22}-\lambda & D_{23} \\ D_{31} & D_{32} & D_{33}-\lambda \end{bmatrix} \tag{F2a-8}$$

$\det(\mathbf{M}(\lambda))$ 是 $\lambda$ 的 3 次多项式：

$$\det(\mathbf{M}(\lambda)) = -\lambda^3 + (D_{11}+D_{22}+D_{33})\lambda^2 - \cdots \tag{F2a-9}$$

$\mathbf{M}(\lambda)^{-1} = \frac{\text{adj}(\mathbf{M}(\lambda))}{\det(\mathbf{M}(\lambda))}$，其中 $\text{adj}(\mathbf{M})$ 是伴随矩阵（adjugate matrix），元素是 $\lambda$ 的 2 次多项式。

代入式 (F2a-7)：

$$\mathbf{d}^\top \frac{\text{adj}(\mathbf{M}(\lambda))^2}{\det(\mathbf{M}(\lambda))^2} \mathbf{d} = g^2 \tag{F2a-10}$$

整理得：

$$\mathbf{d}^\top \text{adj}(\mathbf{M}(\lambda))^2 \mathbf{d} = g^2 \det(\mathbf{M}(\lambda))^2 \tag{F2a-11}$$

左边是 $\lambda$ 的 4 次多项式（$\text{adj}$ 是 2 次，平方后 4 次），右边是 $\lambda$ 的 6 次多项式（$\det$ 是 3 次，平方后 6 次）。

移项后得到 6 次多项式：

$$c_6\lambda^6 + c_5\lambda^5 + c_4\lambda^4 + c_3\lambda^3 + c_2\lambda^2 + c_1\lambda + c_0 = 0 \tag{F2a-12}$$

### 2.5 代码实现

代码直接展开式 (F2a-12) 的系数，得到 7 个系数 $[c_6, c_5, c_4, c_3, c_2, c_1, c_0]$：

```cpp
// helper.h:183-319
static Eigen::Matrix<double, 7, 1> compute_dongsi_coeff(
    Eigen::MatrixXd &D, const Eigen::MatrixXd &d, double gravity_mag) {
    
    // 提取矩阵元素
    double D1_1 = D(0,0), D1_2 = D(0,1), D1_3 = D(0,2);
    double D2_1 = D(1,0), D2_2 = D(1,1), D2_3 = D(1,2);
    double D3_1 = D(2,0), D3_2 = D(2,1), D3_3 = D(2,2);
    double d1 = d(0), d2 = d(1), d3 = d(2);
    double g = gravity_mag;
    
    // 计算系数（直接展开）
    Eigen::Matrix<double, 7, 1> coeff;
    coeff(6) = ...;  // 常数项 c_0
    coeff(5) = ...;  // c_1
    coeff(4) = ...;  // c_2
    coeff(3) = ...;  // c_3
    coeff(2) = ...;  // c_4
    coeff(1) = ...;  // c_5
    coeff(0) = 1;    // c_6 = 1 (首一多项式)
    
    return coeff;
}
```

**系数表达式**（以 $c_0$ 为例，其余类似）：

$$c_0 = -\frac{1}{g^2}\big(-D_{11}^2D_{22}^2D_{33}^2g^2 + D_{11}^2D_{22}^2d_3^2 + \cdots + D_{23}^2D_{32}^2d_1^2\big) \tag{F2a-13}$$

完整表达式见代码 [helper.h:204-247](src/open_vins/ov_init/src/utils/helper.h#L204-L247)。

---

## 3. 代码实现：逐行对照

### 3.1 输入提取

```cpp
// helper.h:186-193
double D1_1 = D(0, 0), D1_2 = D(0, 1), D1_3 = D(0, 2);
double D2_1 = D(1, 0), D2_2 = D(1, 1), D2_3 = D(1, 2);
double D3_1 = D(2, 0), D3_2 = D(2, 1), D3_3 = D(2, 2);
double d1 = d(0, 0), d2 = d(1, 0), d3 = d(2, 0);
double g = gravity_mag;
```

**理论对应**：式 (F2a-8) 的矩阵元素提取。

### 3.2 平方项预计算

```cpp
// helper.h:196-200
double D1_1_sq = D1_1 * D1_1, D1_2_sq = D1_2 * D1_2, D1_3_sq = D1_3 * D1_3;
double D2_1_sq = D2_1 * D2_1, D2_2_sq = D2_2 * D2_2, D2_3_sq = D2_3 * D2_3;
double D3_1_sq = D3_1 * D3_1, D3_2_sq = D3_2 * D3_2, D3_3_sq = D3_3 * D3_3;
double d1_sq = d1 * d1, d2_sq = d2 * d2, d3_sq = d3 * d3;
double g_sq = g * g;
```

**目的**：避免重复计算平方项，提高计算效率。

### 3.3 系数计算

```cpp
// helper.h:202-315
Eigen::Matrix<double, 7, 1> coeff = Eigen::Matrix<double, 7, 1>::Zero();
coeff(6) = ...;  // 常数项（最复杂，~50 项）
coeff(5) = ...;  // 一次项
coeff(4) = ...;  // 二次项
coeff(3) = ...;  // 三次项
coeff(2) = ...;  // 四次项
coeff(1) = ...;  // 五次项
coeff(0) = 1;    // 六次项系数 = 1（首一多项式）
```

**理论对应**：式 (F2a-12) 的系数展开。

### 3.4 返回系数

```cpp
return coeff;
```

返回的 `coeff` 用于构造伴随矩阵（见 F1b）。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\mathbf{D}$ | 3×3 对称矩阵 | `D` | `:183` |
| $\mathbf{d}$ | 3×1 向量 | `d` | `:183` |
| $g$ | 重力模长 | `g` | `:193` |
| $c_6, \ldots, c_0$ | 6 次多项式系数 | `coeff(0)..coeff(6)` | `:203-315` |

---

## 5. 关键设计决策

### 5.1 为什么直接展开系数而非符号计算

| 方法 | 优点 | 缺点 |
|------|------|------|
| **直接展开**（代码采用） | 快速、无依赖、确定性 | 代码长、难维护 |
| 符号计算（SymPy 等） | 易维护、可读性好 | 需要外部依赖、运行时开销 |
| 数值微分 | 简单 | 精度低、不稳定 |

**选择理由**：
1. 初始化只执行一次，代码长度不是问题
2. 无外部依赖，部署简单
3. 解析表达式数值稳定

### 5.2 为什么是 6 次多项式

- $\det(\mathbf{M}(\lambda))$ 是 3 次（3×3 矩阵行列式）
- $\text{adj}(\mathbf{M}(\lambda))$ 元素是 2 次（余子式）
- 式 (F2a-11)：$\mathbf{d}^\top\text{adj}^2\mathbf{d}$ 是 4 次，$g^2\det^2$ 是 6 次
- 最高次项决定多项式次数 = 6

### 5.3 伴随矩阵 vs 逆矩阵

代码注释提到 "matlab constants"，暗示系数可能从 MATLAB 符号计算导出。

伴随矩阵 $\text{adj}(\mathbf{M})$ 的性质：
- $\mathbf{M}^{-1} = \frac{\text{adj}(\mathbf{M})}{\det(\mathbf{M})}$
- $\text{adj}(\mathbf{M})$ 元素是 $\lambda$ 的多项式（无分母）
- 避免分母为零的问题

---

## 6. 完整流程图

```
输入: D (3×3), d (3×1), g (标量)
        │
        ├─ 1. 提取矩阵元素
        │     ─ D1_1, D1_2, ..., D3_3, d1, d2, d3, g
        │
        ├─ 2. 预计算平方项
        │     └─ D1_1_sq, D2_2_sq, ..., d1_sq, d2_sq, d3_sq, g_sq
        │
        ├─ 3. 计算多项式系数
        │     ├─ coeff(6) = 常数项 (~50 项展开)
        │     ├─ coeff(5) = 一次项
        │     ├─ coeff(4) = 二次项
        │     ├─ coeff(3) = 三次项
        │     ├─ coeff(2) = 四次项
        │     ├─ coeff(1) = 五次项
        │     └─ coeff(0) = 1 (六次项)
        │
        └─ 4. 返回 coeff
              └─ 用于构造伴随矩阵 → 特征值分解 → 求 λ
```

---

## 7. 与其他文档的关联

- **上游**：`F1b_DynamicInitializer-algorithm-deep-dive.md`（调用本函数）
- **下游**：伴随矩阵构造 + 特征值分解（F1b 中描述）
- **相关**：Dong-Si & Mourikis 2012 论文（理论来源）

### 待详细补充项

- **F2a-1**：伴随矩阵理论（Gantmacher《矩阵论》）
- **F2a-2**：6 次多项式系数的手动推导过程（验证代码正确性）

> 以上子文档暂不展开，待需要逐项深挖时再补。
