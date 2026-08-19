# PoseJPL 复合类型 精读报告 (S3)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/types/PoseJPL.h`（依赖 `JPLQuat.h` + `Vec.h`）
> **对应论文**: Trawny & Roumeliotis 2005 TR [40] (JPL 姿态)、Huang 2009 [20]/Li&Mourikis 2013 [27] (FEJ 一致性)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S3（Phase 3 函数级路线，阶段 A 状态表示基础）

> 上游：S1 `Type`（基类）+ S2 `JPLQuat`（旋转）。本篇聚焦**复合类型（composite）的第一层**——`PoseJPL` 如何把 `JPLQuat`(旋转) 与 `Vec(3)`(位置) 组合成 6-DOF 位姿，并通过 `set_local_id` + `check_if_subvariable` 把协方差分块映射到连续内存。下游：S4 `IMU` 再组合 `PoseJPL` 成 15-DOF。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
   S2 JPLQuat ─┐
   S4 Vec ─────┼─► S3 PoseJPL(=JPLQuat+Vec, 6DOF) ─┐
               │                                    ├─► S4 IMU(=PoseJPL+3×Vec, 15DOF)
   (Vec 单独用)┘                                    ┘            │
                                                          State::_variables[0] = IMU
                                                          │ 通过 set_local_id 分配协方差块
                                                          ▼
                              StateHelper::EKFPropagation / EKFUpdate / marginalize
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 状态表示复合层（把基础流形拼成 6-DOF 位姿） |
| 输入 | 误差扰动 `dx`（6 维：前 3 旋转、后 3 位置） |
| 输出 | 更新后的位姿估计；子变量 `_q`/`_p` 的 `_id` 连续排布于协方差 |
| 调用方 | `IMU`（组合它）、`State`（IMU clone 内含 PoseJPL）、`Propagator` |
| 被调用方 | `JPLQuat::update` / `Vec::update` / `quat_ops` 工具 |

### 1.2 一句话概括

`PoseJPL` 把 `JPLQuat`(旋转) 与 `Vec(3)`(位置) 组合成 6-DOF 位姿，演示了 type-based 系统里**复合类型的标准范式**：`update` 按子变量顺序分别施加各自的 boxplus（旋转左乘、位置加法），`set_local_id` 把子变量 `_id` 连续排布，`check_if_subvariable` 让上层按子类型定位协方差块。

---

## Section 2: 核心数据结构

### 2.1 类定义

```cpp
// 文件: ov_core/src/types/PoseJPL.h:37-187
class PoseJPL : public Type {
public:
    PoseJPL() : Type(6) {                          // ↑ 误差态 = 6 (3 旋转 + 3 位置)
        _q = make_shared<JPLQuat>();               // 子变量: 旋转
        _p = make_shared<Vec>(3);                  // 子变量: 位置
        Eigen::Matrix<double,7,1> pose0; pose0.setZero(); pose0(3)=1.0;  // [0,0,0,1, 0,0,0]
        set_value_internal(pose0); set_fej_internal(pose0);
    }

    void set_local_id(int new_id) override {       // ↓ 关键: 子变量 _id 连续排布
        _id = new_id;
        _q->set_local_id(new_id);
        _p->set_local_id(new_id + ((new_id!=-1)? _q->size() : 0));  // 位置块接旋转块后
    }

    void update(const Eigen::VectorXd &dx) override {   // boxplus: 旋转用 JPLQuat, 位置用加法
        assert(dx.rows()==_size);                        // dx(0:3)=δθ, dx(3:6)=δp
        Eigen::Matrix<double,7,1> newX = _value;
        Eigen::Matrix<double,4,1> dq; dq << .5*dx.block(0,0,3,1), 1.0; dq=quatnorm(dq);
        newX.block(0,0,4,1) = quat_multiply(dq, quat());  // 旋转: 左乘 (S2)
        newX.block(4,0,3,1) += dx.block(3,0,3,1);          // 位置: 欧式加法
        set_value(newX);
    }

    std::shared_ptr<Type> check_if_subvariable(const std::shared_ptr<Type> c) override {
        if (c==_q) return _q; else if (c==_p) return _p; return nullptr;  // 子变量定位
    }

    Eigen::Matrix3d Rot()/Rot_fej()  { return _q->Rot()/_q->Rot_fej(); }  // 转发旋转
    Eigen::Vector4d quat()/quat_fej(){ return _q->value()/_q->fej(); }
    Eigen::Vector3d pos()/pos_fej()  { return _p->value()/_p->fej(); }

    std::shared_ptr<JPLQuat> q() { return _q; }   // 子类型访问器
    std::shared_ptr<Vec>     p() { return _p; }

protected:
    std::shared_ptr<JPLQuat> _q;   // 旋转子变量
    std::shared_ptr<Vec>     _p;   // 位置子变量

    void set_value_internal(const MatrixXd &nv) {   // 拆 7 维 → 子变量
        _q->set_value(nv.block(0,0,4,1));
        _p->set_value(nv.block(4,0,3,1));
        _value = nv;
    }
    void set_fej_internal(const MatrixXd &nv) {
        _q->set_fej(nv.block(0,0,4,1));
        _p->set_fej(nv.block(4,0,3,1));
        _fej = nv;
    }
};
```

### 2.2 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 物理含义 | 类型 | 位置 |
|---------|---------|---------|------|------|
| `_q` (JPLQuat) | ${}^{I}_{G}\bar{q}$ | 旋转（JPL 四元数） | `JPLQuat` | `PoseJPL.h:147` |
| `_p` (Vec) | ${}^{I}_{G}\mathbf{p}$ | 位置 | `Vec(3)` | `PoseJPL.h:150` |
| `dx` 顺序 | $[\delta\boldsymbol{\theta},\delta\mathbf{p}]$ | 6 维误差态 | `VectorXd(6)` | `PoseJPL.h:72` |

> **6-DOF 最小表示**：旋转用 3 维误差（轴角）+ 位置 3 维 = 6，而非 7（四元数 4 + 位置 3）。这是误差状态 EKF 的特征——存储 `_value` 是 7 维，但误差态 `_size=6`。

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 复合 boxplus（位姿分块更新）

#### 论文原文（误差状态 EKF，MSCKF [32]）

位姿误差态 $\delta\boldsymbol{\xi}_{pose}=[\delta\boldsymbol{\theta}_I,\delta\mathbf{p}_I]^\top$，子块独立施加：
- 旋转：$\bar{q}\leftarrow[\tfrac12\delta\boldsymbol{\theta};1]\otimes\bar{q}$（左乘，S2）
- 位置：$\mathbf{p}\leftarrow\mathbf{p}+\delta\mathbf{p}$（欧式加法）

#### 对应代码

```cpp
// 文件: ov_core/src/types/PoseJPL.h:74-91
void update(const Eigen::VectorXd &dx) override {
    assert(dx.rows() == _size);                  // ↑ 校验 6 维
    Eigen::Matrix<double,7,1> newX = _value;     // 存储 7 维(4+3)

    Eigen::Matrix<double,4,1> dq;
    dq << .5*dx.block(0,0,3,1), 1.0;            // ↑ δθ → 扰动四元数
    dq = ov_core::quatnorm(dq);
    newX.block(0,0,4,1) = quat_multiply(dq, quat());  // ↑ 旋转: 左乘 (S2)
    newX.block(4,0,3,1) += dx.block(3,0,3,1);          // ↑ 位置: 加法
    set_value(newX);
}
```

#### 对照注释

| 子块 | 公式 | 代码行 | 偏移 |
|------|------|--------|------|
| 旋转 δθ(0:3) | 左乘四元数 | `:81,:85` | val(0:4) |
| 位置 δp(3:6) | 加法 | `:88` | val(4:7) |

> **顺序契约**：`dx` 块顺序 = `[δθ, δp]`，`set_value_internal` 拆 7 维也按 `[4,3]`。上层构造 `dx` 必须守此顺序（见 S4 IMU 的 15 维扩展）。

### 3.2 协方差块映射：`set_local_id` 连续排布

#### 推导说明

EKF 协方差 $\mathbf{P}$ 是单个大矩阵。每个状态变量需在 $\mathbf{P}$ 中占一个**连续子块**，由左上角索引 `_id` 与大小 `_size` 决定。复合类型必须让子变量 `_id` **首尾相接**，这样 `dx.block(id,0,size,1)` 才能正确切出增量、`P.block(id_i,id_j,size_i,size_j)` 取出子协方差。

#### 对应代码

```cpp
// 文件: ov_core/src/types/PoseJPL.h:63-67
void set_local_id(int new_id) override {
    _id = new_id;
    _q->set_local_id(new_id);                                        // 旋转块从 new_id 起
    _p->set_local_id(new_id + ((new_id!=-1)? _q->size() : 0));       // 位置块接旋转后
}
```

#### 对照注释

| 子变量 | 起始 `_id` | 大小 | 协方差块 |
|--------|----------|------|---------|
| `_q` | `new_id` | 3 | `P[id..id+3, *]` |
| `_p` | `id+3` | 3 | `P[id+3..id+6, *]` |

> **`new_id==-1` 守卫**：当 `_id=-1`（变量不在协方差中，如临时 clone），偏移 `+size` 被跳过，避免错误索引。

### 3.3 子变量定位：`check_if_subvariable`

```cpp
// 文件: ov_core/src/types/PoseJPL.h:112-119
std::shared_ptr<Type> check_if_subvariable(const std::shared_ptr<Type> check) override {
    if (check == _q) return _q;
    else if (check == _p) return _p;
    return nullptr;
}
```

> 给定某子类型指针，找回它所属的这个 `PoseJPL`，并可直接取 `_q`/`_p` 的 `_id` 写协方差。这是复合类型能"分层访问"的关键（呼应 S1 Section 6 易错点 #2：复合类型须正确重写此方法）。

---

## Section 4: 函数调用链

### 4.1 入口函数

```
PoseJPL::update(dx) @ PoseJPL.h:74   ──复合 boxplus──► JPLQuat::update + Vec::update(位置)
PoseJPL::set_local_id(id) @ PoseJPL.h:63  ──协方差映射──► 排布 _q/_p 的 _id
```

### 4.2 完整调用树

```
PoseJPL() @ PoseJPL.h:40
  ├── _q = JPLQuat()      (_q._size=3)
  ├── _p = Vec(3)         (_p._size=3)
  └── set_value_internal(pose0)  → 拆 7 维写 _q/_p + _value

State/IMU 构造时:
  pose->set_local_id(id) @ PoseJPL.h:63
    ├── _q.set_local_id(id)         → 旋转块 id..id+3
    └── _p.set_local_id(id+3)       → 位置块 id+3..id+6

EKFUpdate 时 (StateHelper::EKFUpdate @ StateHelper.cpp:186):
  v->update(dx.block(v->id(), 0, v->size(), 1))  → PoseJPL::update → JPLQuat::update + Vec::update
```

### 4.3 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `PoseJPL::update` | `PoseJPL.h:74` | 6DOF boxplus | `dx(6)` | 写 `_value` | O(1) |
| `PoseJPL::set_local_id` | `PoseJPL.h:63` | 子变量连续排布 | `id` | 设 `_q`/`_p` 的 `_id` | O(1) |
| `PoseJPL::check_if_subvariable` | `PoseJPL.h:112` | 定位 JPLQuat/Vec | `Type*` | 子变量或 nullptr | O(1) |

---

## Section 5: 关键参数清单

### 5.1 本模块涉及的参数

| 参数名 | 默认值 | 定义位置 | 含义 | 敏感度 |
|--------|--------|---------|------|--------|
| `PoseJPL._size` | `6` | `PoseJPL.h:40` | 误差态维度(3旋转+3位置) | ★★ 中 |
| `dx` 块顺序 | `[δθ, δp]` | `PoseJPL.h:72` | 误差增量布局契约 | ★★★ 高 |

### 5.2 参数影响分析

- **`dx` 顺序是全局契约**：`update` 与 `set_value_internal` 都按 `[4,3]` 切块。上层（Propagator/Updater）拼 `dx` 必须遵守，否则旋转/位置错位。
- **`_size` 求和 = 协方差块大小**：`StateHelper::EKFPropagation` 用 `order[i]->id()+size()==order[i+1]->id()` 断言块连续（`StateHelper.cpp:49`）——若 `set_local_id` 写错，这里直接 assert 失败。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 复合类型须手动维护 `set_local_id`/`check_if_subvariable`/`set_*_internal` 三处一致 | ★★ 中 | 新增复合类型漏写某处 | 协方差块错位或子变量找不到 | 抽基类模板自动递归 |
| 2 | `check_if_subvariable` 依赖指针 `==` 比较 | ★ 低 | 同一子变量被多个复合共享 | 返回首个匹配，可能歧义 | 当前无共享，安全 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `PoseJPL` | Eigen `Isometry3d` | OKVIS/ceres | g2o `VertexSE3` |
|---------|-------------------|---------------------|-------------|-----------------|
| 姿态表示 | JPL 四元数 + 位置 | 4x4 变换矩阵 | 四元数/SO3 | 四元数 |
| 误差态 | 最小(3+3) | 无原生 | 参数块 | 无 |
| 协方差映射 | `set_local_id` 连续 | 无 | problem 管理 | optimizer 管理 |
| 复合定位 | `check_if_subvariable` | 无 | 无 | 无 |
| FEJ | 原生(`_fej`/`_Rfej`) | 无 | 无 | 无 |

### 设计取舍分析

- 选择 **复合 Type + 最小表示 + set_local_id 连续映射** 是因为：误差状态 EKF 需要固定维度误差向量 `dx` 与协方差 $\mathbf{P}$ 严格对齐；复合类型让位姿在数学上分开推导雅可比，又在内存上连续排布——`StateHelper` 因此能用统一循环处理任意状态（S1 提到的统一 update 在此落地）。
- 在 **新增自定义位姿变量** 时当前方案可能不足：必须同步实现三件套，否则破坏契约——这是 type-based 系统的扩展性代价。
