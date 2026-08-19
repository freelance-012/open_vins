# IMU / Vec 状态类型 精读报告 (S4)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/types/{IMU.h, Vec.h}`（IMU 依赖 `PoseJPL.h`）
> **对应论文**: MSCKF Mourikis 2007 [32] (15-DOF IMU 状态)、Trawny 2005 [40] (JPL 姿态)、Huang 2009 [20]/Li&Mourikis 2013 [27] (FEJ 一致性)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S4（Phase 3 函数级路线，阶段 A 状态表示基础）

> 上游：S1 `Type` + S2 `JPLQuat` + S3 `PoseJPL`。本篇聚焦**复合类型的第二层**——`IMU` 把 `PoseJPL`(6DOF) 与三个 `Vec(3)`(速度/陀螺bias/加速度bias) 组合成 15-DOF IMU 状态，并演示 `check_if_subvariable` 的**递归**写法；以及最简派生 `Vec`（欧式加法）。`IMU` 即 `State` 中的 `_variables[0]`，是整条滤波链的状态根。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
   S3 PoseJPL(6DOF) ─┐
   S4 Vec(3) ×3 ────┼─► S4 IMU(=PoseJPL+Vec+Vec+Vec, 15DOF)
                     │        │
                     │        ▼
                     │   State::_variables[0] = IMU
                     │        │ 通过 set_local_id 分配协方差块
                     ▼        ▼
              StateHelper::EKFPropagation / EKFUpdate / marginalize
              （按 _id / _size 切协方差块；IMU 是状态根，clone 增广出滑窗）
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 状态表示复合顶层（MSCKF 的 15-DOF IMU 状态根） |
| 输入 | 误差扰动 `dx`（15 维，按约定顺序） |
| 输出 | 更新后的 IMU 估计；各子变量 `_id` 连续排布于协方差 |
| 调用方 | `State`（构造为 `_variables[0]`）、`Propagator`/`StateHelper` |
| 被调用方 | `PoseJPL::update` / `Vec::update` / `quat_ops` 工具 |

### 1.2 一句话概括

`IMU` 是 OpenVINS 滤波状态的**根变量**：15-DOF = `PoseJPL`(6) + `Vec`速度(3) + `Vec`陀螺bias(3) + `Vec`加速度bias(3)。它递归复用 S3 的复合范式，`update` 按 `[q,p,v,bg,ba]` 顺序施加各子块 boxplus，`set_local_id` 递归连续排布，`check_if_subvariable` 递归下钻到 `JPLQuat`/`Vec`。`Vec` 是最简派生——纯欧式加法，无流形。

---

## Section 2: 核心数据结构

### 2.1 IMU 类（复合 15-DOF）

```cpp
// 文件: ov_core/src/types/IMU.h:37-232
class IMU : public Type {
public:
    IMU() : Type(15) {                          // ↑ 误差态 = 15
        _pose = make_shared<PoseJPL>();         // 子: 6-DOF 位姿
        _v = make_shared<Vec>(3);               // 子: 速度
        _bg = make_shared<Vec>(3);              // 子: 陀螺 bias
        _ba = make_shared<Vec>(3);              // 子: 加速度 bias
        Eigen::VectorXd imu0=VectorXd::Zero(16,1); imu0(3)=1.0;
        set_value_internal(imu0); set_fej_internal(imu0);
    }

    void set_local_id(int new_id) override {    // ↓ 递归连续排布 (pose→v→bg→ba)
        _id = new_id;
        _pose->set_local_id(new_id);
        _v ->set_local_id(_pose->id() + ((new_id!=-1)? _pose->size() : 0));
        _bg->set_local_id(_v->id()    + ((new_id!=-1)? _v->size()    : 0));
        _ba->set_local_id(_bg->id()   + ((new_id!=-1)? _bg->size()   : 0));
    }

    void update(const Eigen::VectorXd &dx) override {  // 顺序: q(0:3) p(3:6) v(6:9) bg(9:12) ba(12:15)
        assert(dx.rows()==_size);
        Eigen::Matrix<double,16,1> newX=_value;
        Eigen::Matrix<double,4,1> dq; dq << .5*dx.block(0,0,3,1),1.0; dq=quatnorm(dq);
        newX.block(0,0,4,1) = quat_multiply(dq, quat());   // 旋转
        newX.block(4,0,3,1) += dx.block(3,0,3,1);           // 位置
        newX.block(7,0,3,1) += dx.block(6,0,3,1);           // 速度
        newX.block(10,0,3,1)+= dx.block(9,0,3,1);           // 陀螺 bias
        newX.block(13,0,3,1)+= dx.block(12,0,3,1);          // 加速度 bias
        set_value(newX);
    }

    std::shared_ptr<Type> check_if_subvariable(const std::shared_ptr<Type> c) override {
        if (c==_pose) return _pose;
        else if (c==_pose->check_if_subvariable(c)) return _pose->check_if_subvariable(c); // 递归进 PoseJPL
        else if (c==_v)  return _v;
        else if (c==_bg) return _bg;
        else if (c==_ba) return _ba;
        return nullptr;
    }

    // 转发访问器: Rot/quat/pos/vel/bias_g/bias_a 及其 *_fej 版本
protected:
    std::shared_ptr<PoseJPL> _pose;
    std::shared_ptr<Vec> _v, _bg, _ba;
    void set_value_internal(const MatrixXd &nv) {   // 拆 16 维 → 4 子变量
        _pose->set_value(nv.block(0,0,7,1));
        _v->set_value(nv.block(7,0,3,1));
        _bg->set_value(nv.block(10,0,3,1));
        _ba->set_value(nv.block(13,0,3,1));
        _value=nv;
    }
    void set_fej_internal(const MatrixXd &nv) { /* 同上, 写 _fej + 各子 fej */ }
};
```

### 2.2 Vec 类（欧式加法，最简派生）

```cpp
// 文件: ov_core/src/types/Vec.h:32-69
class Vec : public Type {
public:
    Vec(int dim) : Type(dim) { _value=VectorXd::Zero(dim); _fej=VectorXd::Zero(dim); }
    void update(const Eigen::VectorXd &dx) override {
        assert(dx.rows()==_size);
        set_value(_value + dx);              // ↑ 欧式加法, 无流形 (对比 JPLQuat 左乘)
    }
    std::shared_ptr<Type> clone() override { ... }
    // check_if_subvariable 用基类默认 (返回 nullptr) —— Vec 是叶子, 无子变量
};
```

### 2.3 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 物理含义 | 类型 | 位置 |
|---------|---------|---------|------|------|
| `IMU::_pose`(PoseJPL) | ${}^{I}_{G}\bar{q},{}^{I}_{G}\mathbf{p}$ | 旋转+位置 (6DOF) | `PoseJPL` | `IMU.h:188` |
| `IMU::_v` | ${}^{I}\mathbf{v}$ | 速度 | `Vec(3)` | `IMU.h:191` |
| `IMU::_bg` | $\mathbf{b}_g$ | 陀螺零偏 | `Vec(3)` | `IMU.h:194` |
| `IMU::_ba` | $\mathbf{b}_a$ | 加速度零偏 | `Vec(3)` | `IMU.h:197` |
| `dx` 顺序 | $[\delta\boldsymbol{\theta},\delta\mathbf{p},\delta\mathbf{v},\delta\mathbf{b}_g,\delta\mathbf{b}_a]$ | 15 维误差态 | `VectorXd(15)` | `IMU.h:76` |

> **MSCKF 15-DOF 状态**（Mourikis 2007 [32]）：位姿(6, 最小表示) + 速度(3) + 陀螺bias(3) + 加速度bias(3) = 15。注意 OpenVINS 用**最小表示**（四元数 3 误差 + 向量），而非 16/18 维过参数化——这是误差状态 EKF 的特征。

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 复合 boxplus（IMU 分块更新）

#### 论文原文（误差状态 EKF，MSCKF [32]）

误差态按变量顺序拼接：$\delta\boldsymbol{\xi}=[\delta\boldsymbol{\theta}_I,\delta\mathbf{p}_I,\delta\mathbf{v}_I,\delta\mathbf{b}_g,\delta\mathbf{b}_a]^\top$，各子块独立施加：
- 旋转：$\bar{q}\leftarrow[\tfrac12\delta\boldsymbol{\theta};1]\otimes\bar{q}$（左乘，S2）
- 位置/速度/bias：$\mathbf{x}\leftarrow\mathbf{x}+\delta\mathbf{x}$（欧式加法，见 `Vec`）

#### 对应代码（IMU::update）

```cpp
// 文件: ov_core/src/types/IMU.h:78-96
void update(const Eigen::VectorXd &dx) override {
    assert(dx.rows() == _size);                  // ↑ 校验 15 维
    Eigen::Matrix<double,16,1> newX = _value;    // 存储 16 维(4+3+3+3+3)

    Eigen::Matrix<double,4,1> dq;
    dq << .5*dx.block(0,0,3,1), 1.0;             // ↑ δθ → 扰动四元数
    dq = ov_core::quatnorm(dq);
    newX.block(0,0,4,1) = quat_multiply(dq, quat());  // ↑ 旋转: 左乘 (S2)
    newX.block(4,0,3,1) += dx.block(3,0,3,1);    // ↑ 位置: 加法
    newX.block(7,0,3,1) += dx.block(6,0,3,1);    // ↑ 速度: 加法
    newX.block(10,0,3,1)+= dx.block(9,0,3,1);    // ↑ 陀螺 bias: 加法
    newX.block(13,0,3,1)+= dx.block(12,0,3,1);   // ↑ 加速度 bias: 加法
    set_value(newX);
}
```

#### 对照注释

| 子块 | 公式 | 代码行 | 偏移 |
|------|------|--------|------|
| 旋转 δθ(0:3) | 左乘四元数 | `:85,:88` | val(0:4) |
| 位置 δp(3:6) | 加法 | `:89` | val(4:7) |
| 速度 δv(6:9) | 加法 | `:91` | val(7:10) |
| 陀螺bias δbg(9:12) | 加法 | `:92` | val(10:13) |
| 加速度bias δba(12:15) | 加法 | `:93` | val(13:16) |

> **顺序契约是全局契约**：`dx` 块顺序 = `[q,p,v,bg,ba]`，且 `set_value_internal` 拆 16 维时也按 `[4,3,3,3,3]` 切。Propagator/Updater 构造 `dx` 时必须遵守此顺序，否则偏差。

### 3.2 协方差块映射：`set_local_id` 连续排布

#### 对应代码（递归连续排布）

```cpp
// 文件: ov_core/src/types/IMU.h:64-70
void set_local_id(int new_id) override {
    _id = new_id;
    _pose->set_local_id(new_id);                                       // 位姿块从 new_id 起
    _v ->set_local_id(_pose->id() + ((new_id!=-1)? _pose->size() : 0));// 速度块接位姿后
    _bg->set_local_id(_v->id()    + ((new_id!=-1)? _v->size()    : 0));// 陀螺块接速度后
    _ba->set_local_id(_bg->id()   + ((new_id!=-1)? _bg->size()   : 0));// 加速块接陀螺后
}
```

#### 对照注释

| 子变量 | 起始 `_id` | 大小 | 协方差块 |
|--------|----------|------|---------|
| `_pose` | `new_id` | 6 | `P[id..id+6, *]` |
| `_v` | `id+6` | 3 | `P[id+6..id+9, *]` |
| `_bg` | `id+9` | 3 | `P[id+9..id+12, *]` |
| `_ba` | `id+12` | 3 | `P[id+12..id+15, *]` |

> **递归特性**：`_pose->set_local_id` 内部又调用 `_q`/`_p` 的 `set_local_id`，因此一次 `imu->set_local_id(0)` 就把整棵 15-DOF 子树的 `_id` 全部正确排布。`StateHelper::augment_clone` 分配新 `_id` 时调此函数，保证滑窗增广后协方差块连续。

### 3.3 子变量定位：`check_if_subvariable`（递归下钻）

#### 对应代码

```cpp
// 文件: ov_core/src/types/IMU.h:117-130
std::shared_ptr<Type> check_if_subvariable(const std::shared_ptr<Type> check) override {
    if (check == _pose) return _pose;
    else if (check == _pose->check_if_subvariable(check))   // ↓ 递归进 PoseJPL
        return _pose->check_if_subvariable(check);          //   (找 JPLQuat / Vec)
    else if (check == _v)  return _v;
    else if (check == _bg) return _bg;
    else if (check == _ba) return _ba;
    return nullptr;
}
```

#### 用处（StateHelper 实证）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:186-187  (EKFUpdate 应用 dx)
for (size_t i = 0; i < state->_variables.size(); i++) {
    state->_variables.at(i)->update(
        dx.block(state->_variables.at(i)->id(), 0,
                 state->_variables.at(i)->size(), 1));   // ↑ 按 _id/_size 切本变量增量
}
```

> `check_if_subvariable` 用于"给定某子类型指针，找回它所属的那个复合变量"——例如 Updater 只对某子块残差用到 `JPLQuat` 时，能递归定位回对应 `IMU`/clone 并取其 `_id` 写协方差。这是复合类型能"分层访问"的关键（呼应 S1 Section 6 易错点 #2）。

---

## Section 4: 函数调用链

### 4.1 入口函数

```
IMU::update(dx) @ IMU.h:78   ──复合 boxplus──► 内部按序调用 JPLQuat::update / Vec::update ×3
IMU::set_local_id(id) @ IMU.h:64  ──协方差映射──► 递归排布 _pose/_v/_bg/_ba 的 _id
```

### 4.2 完整调用树（IMU 构造 → 协方差映射）

```
IMU() @ IMU.h:40
  ├── _pose = PoseJPL()          → JPLQuat()+Vec(3)   (_q._size=3, _p._size=3)
  ├── _v/_bg/_ba = Vec(3)
  └── set_value_internal(imu0)   → 拆 16 维写各子变量 + _value

State 构造时:
  imu_state->set_local_id(0) @ IMU.h:64
    ├── _pose.set_local_id(0)        → _q.id=0, _p.id=0+3=3
    ├── _v.set_local_id(3+3=6)
    ├── _bg.set_local_id(6+3=9)
    └── _ba.set_local_id(9+3=12)      → IMU 占 P[0..15]

EKFUpdate 时:
  StateHelper::EKFUpdate @ StateHelper.cpp:116
    └── for v in _variables: v->update(dx.block(v->id(), 0, v->size(), 1))  @ :186
          └── IMU::update → JPLQuat::update + Vec::update×3
```

### 4.3 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `IMU::update` | `IMU.h:78` | 15DOF boxplus | `dx(15)` | 写 `_value` | O(1) |
| `IMU::set_local_id` | `IMU.h:64` | 递归连续排布 | `id` | 设 4 子变量 `_id` | O(1) |
| `IMU::check_if_subvariable` | `IMU.h:117` | 递归定位子变量 | `Type*` | 子变量或 nullptr | O(1) |
| `Vec::update` | `Vec.h:55` | 欧式加法 | `dx(n)` | 写 `_value` | O(1) |

---

## Section 5: 关键参数清单

### 5.1 本模块涉及的参数

| 参数名 | 默认值 | 定义位置 | 含义 | 敏感度 |
|--------|--------|---------|------|--------|
| `IMU._size` | `15` | `IMU.h:40` | 误差态维度(6+3+3+3) | ★★ 中 |
| `Vec._size` | 构造传入 `dim` | `Vec.h:39` | 向量维度 | ★ 低 |
| `dx` 块顺序 | `[q,p,v,bg,ba]` | `IMU.h:76` | 误差增量布局契约 | ★★★ 高 |

### 5.2 参数影响分析

- **`dx` 顺序是全局契约**：`IMU::update` 与各 `set_*_internal` 的切块都按 `[4,3,3,3,3]`。任何生成 `dx` 的上层（Propagator 的 F 矩阵、Updater 的 H 矩阵）必须严格按此顺序拼装，否则旋转/位置/bias 增量错位 → 估计发散。这是调试 VIO 收敛问题的首选排查点。
- **`_size` 求和 = 协方差块大小**：`StateHelper::EKFPropagation` 用 `order_NEW[i]->id()+size() == order_NEW[i+1]->id()` 断言块连续（`StateHelper.cpp:49`）——若复合类型 `set_local_id` 写错，这里直接 assert 失败。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 复合类型必须手动维护 `set_local_id`/`check_if_subvariable`/`set_*_internal` 三处一致性 | ★★ 中 | 新增复合类型漏写某处 | 协方差块错位或子变量找不到 | 可抽基类模板自动递归 |
| 2 | `check_if_subvariable` 递归依赖指针 `==` 比较 | ★ 低 | 同一子变量被多个复合共享 | 返回首个匹配，可能歧义 | 当前场景无共享，安全 |
| 3 | `dx` 顺序硬编码于注释 + 切块 | ★★ 中 | 上游拼错顺序 | 估计静默发散 | 集中定义为枚举/常量，避免魔法偏移 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS (type-based) | OKVIS (ceres) | VINS-Mono (ceres) | g2o |
|---------|----------------------|---------------|-------------------|-----|
| 状态组织 | 复合 `Type`(IMU=15DOF) | `ceres::Vector` 拼块 | `ceres::Vector` 拼块 | `VertexSE3`/`VertexSBAPoint` |
| 误差表示 | **最小表示(四元数3+向量)** | 四元数/SO3 参数块 | 四元数/SO3 | 四元数 |
| 协方差映射 | `set_local_id` 自动连续排布 | 由 problem 管理 | 由 problem 管理 | 由 optimizer 管理 |
| 复合定位 | `check_if_subvariable` 递归 | 无（扁平参数块） | 无 | 无 |
| FEJ | **原生(_fej/_Rfej)** | 无 | 无 | 无 |
| 精度特点 | 误差态 EKF + FEJ 一致性 [20][27] | BA 最优 | BA 最优 | BA 最优 |

### 设计取舍分析

- 选择 **复合 Type + 最小表示 + set_local_id 连续映射** 是因为：误差状态 EKF 需要固定维度的误差向量 `dx` 与协方差 $\mathbf{P}$ 严格对齐；复合类型让"位姿/速度/bias"既能在数学上分开推导雅可比，又能在内存上连续排布——`StateHelper` 因此能用统一循环 `v->update(dx.block(v->id(),...))` 处理任意状态（S1 提到的统一 update 流程在此落地）。
- 放弃 **ceres/g2o 的扁平参数块** 是因为：图优化无 FEJ 概念，且滑动窗口边缘化在 EKF 框架下表达更自然（citelist [20][27] 的核心贡献依赖 FEJ）。
- 在 **新增自定义状态变量（如轮速、磁力计）** 时当前方案可能不足：必须同步实现 `set_local_id`/`check_if_subvariable`/`set_*_internal` 三件套，否则破坏契约——这是 type-based 系统的扩展性代价。
