# StateHelper 滑窗克隆 精读报告 (S10)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/StateHelper.cpp`
> **对应论文**: Mourikis & Roumeliotis 2007 MSCKF [20] (随机克隆/滑窗)、Li & Mourikis 2013 IJRR [27] (时间偏移雅可比)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S10（Phase 3 函数级路线，阶段 C 状态增广/滑窗克隆）

> 上游：S5（`propagate_and_clone` `:137` 调用 `augment_clone`）、S9（协方差前推后滑窗生长）。本篇聚焦 **`StateHelper::augment_clone` (:579)** 与底层 **`clone` (:341)** —— MSCKF 滑窗的"生长点"：把当前 IMU pose 的深拷贝追加到协方差末尾，使历史姿态成为可观测约束，并可选注入相机-IMU 时间偏移雅可比（[27]）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
propagate_and_clone (S5 :137)
        │  (last_w)
        ▼
augment_clone(state, last_w)  @ :579  ★ 本篇
        ├─ clone(state, _imu->pose())  @ :589 → clone :341
        ├─ _clones_IMU[_timestamp] = pose   (入滑窗)
        └─ (可选) 时间偏移雅可比注入  (do_calib_camera_timeoffset)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 预测第 4 步（滑窗生长） |
| 输入 | `last_w`（末段角速度，来自 S5 `:105-113`） |
| 输出 | 新增 IMU clone 入 `_clones_IMU`，协方差扩维 |
| 调用方 | `Propagator::propagate_and_clone` (S5 :137) |
| 被调用方 | `clone` (:341)、`Type::clone`/`set_local_id`（S1） |

### 1.2 一句话概括

`augment_clone` 是 MSCKF 随机克隆 [20] 的落地：调用 `clone` 把当前 IMU pose 深拷贝到协方差末尾（携带此刻 `value==fej`，即**克隆 FEJ 锚定**），记入 `_clones_IMU[timestamp]` 滑窗；开启时间偏移标定时，注入克隆 pose 对 $t_d$ 的雅可比（[27]）。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp
std::shared_ptr<Type> StateHelper::clone(state, variable_to_clone);  // :341 底层
void StateHelper::augment_clone(state, Eigen::Vector3d last_w);     // :579 滑窗入口
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 滑窗克隆 `augment_clone`

#### 推导说明

随机克隆 [20]：在时刻 $t_k$ 把当前 IMU pose 的**拷贝**追加到状态末尾，协方差为原变量协方差块的复制（克隆态与母态初始完全相关，后续经测量解耦）。若开启时间偏移在线标定（[27]），克隆姿态是偏移 $t_d$ 的函数，需注入雅可比 $\partial\mathbf{n}_c/\partial t_d=[\boldsymbol{\omega};\mathbf{v}]$。

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:579-616
void StateHelper::augment_clone(state, last_w) {
    // 防御: 同时间戳不可再克隆
    if (_clones_IMU.find(_timestamp)!=_clones_IMU.end()) { ... std::exit; } // :582-585
    // 克隆当前 IMU pose 到协方差末尾 (调底层 clone)
    Type* posetemp = clone(state, state->_imu->pose());                    // :589 → :341
    PoseJPL* pose = dynamic_pointer_cast<PoseJPL>(posetemp);               // :592
    if (pose==nullptr) { ... std::exit; }                                  // :593-596
    _clones_IMU[_timestamp] = pose;                                       // :599 ★ 入滑窗

    // 时间偏移在线标定: 克隆 pose 对 dt 的雅可比 ([27])
    if (state->_options.do_calib_camera_timeoffset) {                     // :604
      Eigen::Matrix<double,6,1> dnc_dt=Zero(6,1);
      dnc_dt.block(0,0,3,1) = last_w;                                     // :607 角速度→姿态漂移
      dnc_dt.block(3,0,3,1) = state->_imu->vel();                         // :608 速度→位置漂移
      _Cov.block(0,pose->id(),_Cov.rows(),6)   += _Cov.block(0,dt_id,_,1)*dnc_dt.transpose(); // :611
      _Cov.block(pose->id(),0,6,_Cov.rows())   += dnc_dt*_Cov.block(dt_id,0,1,_Cov.rows());   // :613
    }
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 克隆入滑窗 | `_clones_IMU[_timestamp]=pose` | `:599` |
| $\partial\mathbf{n}_c/\partial t_d=[\boldsymbol{\omega};\mathbf{v}]$ | `dnc_dt=[last_w; vel]` | `:607-608` |
| 协方差增广偏移雅可比 | `_Cov.block(0,pose.id)+= _Cov.block(0,dt_id)*dnc_dtᵀ` | `:611-613` |

### 3.2 底层克隆 `clone`（通用变量深拷贝）

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:341-391
std::shared_ptr<Type> StateHelper::clone(state, variable_to_clone) {
    int total_size = variable_to_clone->size();                             // :344
    int old_size = _Cov.rows();                                             // :345
    int new_loc  = _Cov.rows();                                             // :346 新位置=原末尾
    // 1. 协方差末尾扩 total_size 行/列
    _Cov.conservativeResizeLike(Zero(old_size+total_size, old_size+total_size)); // :349
    // 2. 定位要克隆的变量 (支持子变量 check_if_subvariable, S4)
    std::shared_ptr<Type> new_clone = nullptr;
    for (k=0; k<_variables.size(); k++) {                                   // :356
      Type* type_check = _variables[k]->check_if_subvariable(variable_to_clone); // :360
      if (_variables[k]==variable_to_clone) type_check=_variables[k];      // :361-362
      else if (type_check != variable_to_clone) continue;                   // :363-364
      int old_loc = type_check->id();                                       // :368
      // 3. 复制协方差三块
      _Cov.block(new_loc,new_loc, total_size,total_size) = _Cov.block(old_loc,old_loc,total_size,total_size); // :371 自协方差
      _Cov.block(0,new_loc, old_size,total_size)       = _Cov.block(0,old_loc,old_size,total_size);          // :372 上交叉
      _Cov.block(new_loc,0, total_size,old_size)       = _Cov.block(old_loc,0,total_size,old_size);         // :373 左交叉
      // 4. 生成新克隆 (深拷贝 value/fej), 设到新位置
      new_clone = type_check->clone();                                      // :376 S1 虚 clone
      new_clone->set_local_id(new_loc);                                     // :377 ★ id=新末尾
      break;
    }
    if (new_clone==nullptr) { ... std::exit; }                              // :382-386
    state->_variables.push_back(new_clone);                                 // :389
    return new_clone;
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 协方差扩维 | `conservativeResizeLike(Zero(old+total, old+total))` | `:349` |
| 克隆自协方差块 | `_Cov.block(new_loc,new_loc)=_Cov.block(old_loc,old_loc)` | `:371` |
| 克隆交叉块 | `_Cov.block(0,new_loc)=_Cov.block(0,old_loc)` (左/上) | `:372-373` |
| 新克隆 id | `new_clone->set_local_id(new_loc)` | `:377` |

> **设计要点**：`clone` 是通用机制（**任何 `Type` 含子变量都能克隆**，靠 `check_if_subvariable` S4 递归定位）。克隆体通过 `type_check->clone()`（S1 虚函数）深拷贝，**携带此刻 `value==fej`**——即 S7/FEJ 讨论中的"克隆 FEJ 锚定"：历史 clone 的 FEJ 固定不动，是 MSCKF 一致性 [20][27] 的关键。

> **时间偏移雅可比技术债**：`:610` 代码注释自承此处本可复用 `EKFPropagation`（S9）但暂未重构——属已知待清理项。

---

## Section 4: 函数调用链

```
Propagator::propagate_and_clone (S5 :137)
  └─ StateHelper::augment_clone(state, last_w) @ :579  ★ 本篇
        └─ StateHelper::clone(state, _imu->pose()) @ :341
              ├─ _variables[k]->check_if_subvariable(variable)  (S4 递归)
              ├─ type_check->clone()  (S1 虚 clone 深拷贝)
              └─ new_clone->set_local_id(new_loc)  (S1 连续 id)
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 |
|--------|-----------|------|
| `augment_clone` | `StateHelper.cpp:579` | 滑窗克隆入口 |
| `clone` | `StateHelper.cpp:341` | 通用变量深拷贝到末尾 |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 |
|--------|--------|---------|------|------|
| `do_calib_camera_timeoffset` | false | `StateOptions` | 时序 | 时间偏移在线标定开关 |

> 开启时 `augment_clone` 注入偏移雅可比（`:611-613`，对应 [27]），否则跳过。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 偏移雅可比未复用 EKFPropagation | ★ 低 | `do_calib_camera_timeoffset` | 注释 `:610` 自承技术债 | 重构为统一前推 |
| 2 | 同时间戳克隆崩溃 | ★ 中 | 重复调用 | `:582` 直接 `std::exit` | 理论不应发生（每帧一次） |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `augment_clone` | VINS-Mono | OKVIS | 常规 MSCKF [20] |
|---------|--------------------------|-----------|-------|-----------------|
| 滑窗增广 | `clone`+`augment_clone` 显式 | 无(图优化) | 有 | 有 |
| 时间偏移 | `augment_clone` 雅可比 [27] | 标定 | 标定 | 标定 |
| FEJ 锚定 | 克隆携带 value==fej | 无 | 部分 | 部分 |

### 设计取舍分析

- 选择 **显式 `clone` 深拷贝 + 末尾追加** 是因为：MSCKF 滑窗需要历史姿态作为独立可观测变量；通用 `clone`（支持子变量）避免为 IMU/SLAM/clone 各写一份。
- 选择 **克隆即 FEJ 锚定** 是因为：历史 clone 的雅可比固定在其生成时刻的首估计，保证多帧约束的一致性（[20][27]），这是 OpenVINS FEJ 机制的核心（区别于活跃 `_imu` 每步刷新 FEJ，见 S7）。
