# JPLQuat 四元数类型 精读报告 (S2)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/types/JPLQuat.h` + 依赖 `ov_core/src/utils/quat_ops.h`
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (Indirect Kalman Filter for 3D Attitude Estimation, Eq.9 / Eq.71)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S2（Phase 3 函数级路线，阶段 A 状态表示基础）

> 上游：S1 `Type` 基类（已读）。`JPLQuat` 是 `Type` 的第一个具体派生类，演示了"误差状态维度=3、存储=4"的流形 boxplus 如何落地。下游：S3 `PoseJPL` 用它组合 6-DOF 位姿。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
        Type(基类, S1)  ──继承──►  JPLQuat (本模块)
                                          │ 被组合
                                          ▼
        PoseJPL (S3)  ← IMU (S4)  ← State (ov_msckf)
                                          │ 旋转矩阵 _R / _Rfej
                                          ▼
        Propagator / Updater*  (所有旋转相关雅可比都从 _R / _Rfej 取)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 状态表示最底层（旋转变量的标准实现） |
| 输入 | 误差扰动 `dx`(3)（轴角）；构造初值 `q0`(4) |
| 输出 | 更新后的四元数 `_value` + 旋转矩阵 `_R`；FEJ 版本 `_fej`/`_Rfej` |
| 调用方 | `PoseJPL`/`IMU`/`Landmark` 组合它；`Propagator` 每帧 `update`/`set_fej` |
| 被调用方 | `quat_ops.h`: `quatnorm`(:496)、`quat_multiply`(:180)、`quat_2_Rot`(:152) |

### 1.2 一句话概括

`JPLQuat` 用 **JPL 约定四元数**（标量在末位 `q=[x,y,z,w]`、左乘误差态）表示旋转，重写 `Type::update` 实现**左乘小角扰动**的 boxplus，并缓存对应的旋转矩阵 `_R`/`_Rfej`，供全仓旋转雅可比统一取用。误差态维度 `_size=3`，但存储 `_value` 为 4 维——这是流形类型的通用约定。

---

## Section 2: 核心数据结构

### 2.1 类定义

```cpp
// 文件: ov_core/src/types/JPLQuat.h:92-188
class JPLQuat : public Type {
public:
    JPLQuat() : Type(3) {                     // ↑ 误差状态维度 = 3（轴角），非 4
        Eigen::Vector4d q0 = Eigen::Vector4d::Zero();
        q0(3) = 1.0;                          // ↑ 单位四元数 [0,0,0,1]
        set_value_internal(q0);
        set_fej_internal(q0);                  // ↑ 初值同时固化为 FEJ
    }

    void update(const Eigen::VectorXd &dx) override;   // boxplus 左乘扰动
    void set_value(const Eigen::MatrixXd &nv) override { set_value_internal(nv); }
    void set_fej(const Eigen::MatrixXd &nv)   override { set_fej_internal(nv); }
    std::shared_ptr<Type> clone() override;            // 深拷贝（滑窗克隆）

    Eigen::Matrix<double,3,3> Rot()    const { return _R; }     // 当前旋转矩阵
    Eigen::Matrix<double,3,3> Rot_fej() const { return _Rfej; } // FEJ 旋转矩阵

protected:
    Eigen::Matrix<double,3,3> _R;     // 由 _value 派生的旋转矩阵
    Eigen::Matrix<double,3,3> _Rfej;  // 由 _fej 派生的 FEJ 旋转矩阵

    void set_value_internal(const Eigen::MatrixXd &nv) {   // 写值 + 重算 _R
        _value = nv;
        _R = ov_core::quat_2_Rot(nv);
    }
    void set_fej_internal(const Eigen::MatrixXd &nv) {     // 写 FEJ + 重算 _Rfej
        _fej = nv;
        _Rfej = ov_core::quat_2_Rot(nv);
    }
};
```

### 2.2 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 (Trawny05 [40]) | 物理含义 | 类型 | 位置 |
|---------|--------------------------|---------|------|------|
| `_value` | $\bar{q}=[q_1,q_2,q_3,q_4]^\top$ | 当前 JPL 四元数（w 在末位） | `MatrixXd(4,1)` | `JPLQuat.h`(基类 `_value`) |
| `_fej` | $\check{\bar{q}}$ | 首估计四元数（线性化点） | `MatrixXd(4,1)` | 基类 `_fej` |
| `_R` | ${}^{I}_{G}\mathbf{R}$ | 当前旋转矩阵 | `Matrix3d` | `:154` |
| `_Rfej` | ${}^{I}_{G}\check{\mathbf{R}}$ | FEJ 旋转矩阵 | `Matrix3d` | `:157` |
| `dx` | $\delta\boldsymbol{\theta}$ | 轴角误差扰动（3 维） | `VectorXd(3)` | `update` 入参 |
| `q(3,0)` | $q_4$ | 标量部分（JPL 约定末位） | `double` | 各处 |

> **JPL vs Hamilton 约定**：JPL 把标量 $q_4$ 放**最后**（$\bar{q}=[q_1,q_2,q_3,q_4]^\top$），与常见 Hamilton（$\bar{q}=[w,x,y,z]$）相反。OpenVINS 全仓统一 JPL，且 `quat_multiply` 强制 $q_4>0$ 保证唯一性（`quat_ops.h:190`）。

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 左乘误差态更新（boxplus）

#### 论文原文（Trawny 2005 [40] Eq.71）

$$
{}^{I}_{G}\bar{q} \simeq \begin{bmatrix} \tfrac{1}{2}\delta\boldsymbol{\theta} \\ 1 \end{bmatrix} \otimes {}^{I}_{G}\hat{\bar{q}}
\qquad\text{(左乘，等价于 SO(3) 误差)}
$$

#### 推导说明

- **物理含义**：小角度轴角 $\delta\boldsymbol{\theta}$ 先转成"小扰动四元数" $[\tfrac12\delta\boldsymbol{\theta};1]$，再**左乘**当前估计，得到更新后的四元数。
- **假设条件**：$\delta\boldsymbol{\theta}$ 足够小（EKF 一阶）；左乘而非右乘（JPL 约定 + SO(3) 等价性，见 `JPLQuat.h:78-80`）。
- **适用范围**：所有旋转变量的标准更新；右乘需另写类型。

#### 对应代码

```cpp
// 文件: ov_core/src/types/JPLQuat.h:114-125
void update(const Eigen::VectorXd &dx) override {
    assert(dx.rows() == _size);                 // ↑ 校验 dx 维度 == 3（误差态）

    Eigen::Matrix<double, 4, 1> dq;
    dq << .5 * dx, 1.0;                         // ↑ 构造扰动四元数 [½δθ; 1]
    dq = ov_core::quatnorm(dq);                 // ↑ 单位化维持 |q|=1 约束

    set_value(ov_core::quat_multiply(dq, _value));  // ↑ 左乘：q ← dq ⊗ q_current
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 | 备注 |
|--------|---------|------|------|
| $\begin{bmatrix}\tfrac12\delta\boldsymbol{\theta}\\1\end{bmatrix}$ | `dq << .5*dx, 1.0` | `:120` | JPL 左乘扰动四元数 |
| 单位约束 $|\bar{q}|=1$ | `quatnorm(dq)` | `:121` | `quat_ops.h:496` |
| $\bar{q}\leftarrow dq\otimes\hat{\bar{q}}$ | `quat_multiply(dq, _value)` | `:124` | 左乘，调 `quat_ops.h:180` |

#### 差异记录

| # | 论文描述 | 代码实现 | 偏差原因 | 潜在影响 |
|---|---------|---------|---------|---------|
| 1 | 理论扰动四元数未显式单位化 | 代码 `quatnorm` 一步 | 数值累积需归一 | 无（标准做法） |

### 3.2 四元数乘法（JPL 左乘矩阵 $\mathcal{L}$）

#### 论文原文（Trawny 2005 [40] Eq.9）

$$
\bar{q}\otimes\bar{p}=
\begin{bmatrix}
q_4\mathbf{I}_3+\lfloor\mathbf{q}\times\rfloor & \mathbf{q} \\
-\mathbf{q}^\top & q_4
\end{bmatrix}
\begin{bmatrix}\mathbf{p}\\p_4\end{bmatrix}
$$

#### 对应代码

```cpp
// 文件: ov_core/src/utils/quat_ops.h:180-195
inline Eigen::Matrix<double,4,1> quat_multiply(const Eigen::Matrix<double,4,1> &q,
                                               const Eigen::Matrix<double,4,1> &p) {
    Eigen::Matrix<double,4,4> Qm;
    Qm.block(0,0,3,3) = q(3,0)*I3 - skew_x(q.block(0,0,3,1));  // ↑ q4*I - [q×]  = 左上 3x3
    Qm.block(0,3,3,1) = q.block(0,0,3,1);                      // ↑ q (右上 3x1)
    Qm.block(3,0,1,3) = -q.block(0,0,3,1).transpose();        // ↑ -q^T (左下 1x3)
    Qm(3,3) = q(3,0);                                          // ↑ q4 (右下)
    q_t = Qm * p;                                              // ↑ L(q) * p
    if (q_t(3,0) < 0) q_t *= -1;                               // ↑ 强制 q4>0 唯一性
    return q_t / q_t.norm();                                   // ↑ 归一化
}
```

#### 对照注释

| 公式块 | 对应代码 | 行号 |
|--------|---------|------|
| $q_4\mathbf{I}_3+\lfloor\mathbf{q}\times\rfloor$ | `q(3,0)*I3 - skew_x(q.head<3>())` | `:184` |
| $\mathbf{q}$ | `Qm.block(0,3,3,1)=q.head<3>()` | `:185` |
| $-\mathbf{q}^\top$ | `Qm.block(3,0,1,3)=-q.head<3>().transpose()` | `:186` |
| $q_4$ | `Qm(3,3)=q(3,0)` | `:187` |

### 3.3 四元数 → 旋转矩阵

#### 论文原文（Trawny 2005 [40]，JPL 约定）

$$
\mathbf{R}= (2q_4^2-1)\mathbf{I}_3 - 2q_4[\mathbf{q}\times] + 2\mathbf{q}\mathbf{q}^\top
$$

#### 对应代码

```cpp
// 文件: ov_core/src/utils/quat_ops.h:152-157
inline Eigen::Matrix<double,3,3> quat_2_Rot(const Eigen::Matrix<double,4,1> &q) {
    Eigen::Matrix3d q_x = skew_x(q.block(0,0,3,1));
    return (2*std::pow(q(3,0),2) - 1)*I3
         - 2*q(3,0)*q_x
         + 2*q.block(0,0,3,1)*q.block(0,0,3,1).transpose();
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $(2q_4^2-1)\mathbf{I}_3$ | `(2*pow(q(3),2)-1)*I3` | `:154` |
| $-2q_4[\mathbf{q}\times]$ | `-2*q(3,0)*q_x` | `:154` |
| $2\mathbf{q}\mathbf{q}^\top$ | `2*q.head<3>()*q.head<3>().transpose()` | `:155` |

> **缓存设计**：`set_value_internal`/`set_fej_internal` 每次写四元数时**顺带重算 `_R`/`_Rfej`**（`JPLQuat.h:171/186`）。好处是 Propagator/Updater 取旋转矩阵时零成本；代价是每次 update 多算一次 `quat_2_Rot`（O(1)，可忽略）。

---

## Section 4: 函数调用链

### 4.1 入口函数

```
JPLQuat::update(dx) @ JPLQuat.h:114
  └→ 功能: 轴角扰动左乘更新四元数
     参数: (dx: VectorXd(3)) 误差态增量
     返回: void（写 _value + 重算 _R）
```

### 4.2 完整调用树

```
JPLQuat::update(dx) @ JPLQuat.h:114
  ├── dq = [.5*dx; 1]
  ├── quatnorm(dq)                          @ quat_ops.h:496   (归一)
  └── set_value(quat_multiply(dq, _value))
        ├── quat_multiply(dq, _value)       @ quat_ops.h:180  (L(q)*p)
        └── set_value_internal(nv)          @ JPLQuat.h:163
              ├── _value = nv
              └── _R = quat_2_Rot(nv)       @ quat_ops.h:152  (重算旋转矩阵)

JPLQuat::set_fej(nv) @ JPLQuat.h:137
  └── set_fej_internal(nv) @ JPLQuat.h:178
        ├── _fej = nv
        └── _Rfej = quat_2_Rot(nv)
```

### 4.3 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `JPLQuat::update` | `JPLQuat.h:114` | 左乘扰动 boxplus | `dx(3)` | 写 `_value`+`_R` | O(1) |
| `quat_multiply` | `quat_ops.h:180` | JPL 四元数乘法 (Eq.9) | `q,p(4)` | `q⊗p(4)` | O(1) |
| `quat_2_Rot` | `quat_ops.h:152` | 四元数→旋转矩阵 | `q(4)` | `R(3x3)` | O(1) |
| `quatnorm` | `quat_ops.h:496` | 归一化 | `q(4)` | 单位四元数 | O(1) |
| `Rot_fej` | `JPLQuat.h:150` | 取 FEJ 旋转矩阵 | — | `_Rfej` | O(1) |

---

## Section 5: 关键参数清单

### 5.1 本模块涉及的参数

| 参数名 | 默认值 | 定义位置 | 所属类别 | 含义 | 敏感度 |
|--------|--------|---------|---------|------|--------|
| `_size` | `3`（`Type(3)`） | `JPLQuat.h:95` | 结构 | 误差态维度（轴角），**≠存储 4** | ★★ 中 |
| 单位四元数初值 | `[0,0,0,1]` | `JPLQuat.h:96-99` | 初始化 | 恒等旋转 | ★ 低 |

### 5.2 参数影响分析

- **`_size=3` 是流形核心约定**：误差态用最小表示（轴角 3 维），但状态存储用 4 维四元数。`update` 里 `assert(dx.rows()==_size)` 保证传入 3 维；这要求**所有调用方传入的 `dx` 块必须是 3 维**，否则崩。S4 `IMU`/`PoseJPL` 组合多个 `JPLQuat`/`Vec` 时，总误差维度是各子变量 `_size` 之和——这是 StateHelper 拼协方差块的依据。
- **唯一性 `q4>0`**：`quat_multiply` 强制标量非负（`quat_ops.h:190`）。若某处手算四元数未走 `quat_multiply`，可能出现 $q_4<0$ 的双覆盖问题（同一旋转两种表示），导致雅可比符号错。统一走库函数可避免。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `set_value`/`set_fej` 只更新四元数+旋转矩阵，不校验单位长度 | ★ 低 | 外部直接写非单位四元数 | `_R` 基于非单位四元数算，旋转失真 | 可在 `set_value_internal` 加 `quatnorm` |
| 2 | `_R`/`_Rfej` 缓存与 `_value`/`_fej` 强耦合，需成对维护 | ★ 低 | 新增派生类型忘重写 `set_*_internal` | 旋转矩阵与四元数不一致 | S1 已提：基类加断言/默认实现 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `JPLQuat` | Eigen `Quaterniond` | OKVIS/VI-SLAM | VINS-Mono |
|---------|-------------------|---------------------|---------------|-----------|
| 四元数约定 | **JPL**（标量末位） | Hamilton（标量首位） | JPL/Hamilton 视实现 | Eigen(Hamilton) |
| 误差态 | 左乘轴角(3) | 无原生误差态 | 取决于类型 | 无 |
| 旋转矩阵缓存 | **预存 `_R`/`_Rfej`** | 按需 `toRotationMatrix()` | 按需 | 按需 |
| FEJ 支持 | **原生 `_Rfej`** | 无 | 部分 | 无 |
| 唯一性 | q4>0 强制 | 不强制 | 视实现 | 不强制 |

### 设计取舍分析

- 选择 **JPL 约定 + 预存 `_Rfej`** 是因为：与 Trawny 2005 [40] 推导一致，误差态直接是 SO(3) 等价（左乘），且 FEJ 旋转矩阵缓存让雅可比计算零成本——这对滤波式实时 VIO 很关键。
- 放弃 **Eigen 原生 `Quaterniond`** 是因为：Eigen 用 Hamilton 约定、无误差态/无 FEJ，无法直接承载 OpenVINS 的滤波数学。
- 在 **与 Hamilton 约定的外部库（如 Eigen/ROS tf）互操作** 时需小心标量位置转换——这是全仓统一 JPL 的代价，但换来内部一致性。
