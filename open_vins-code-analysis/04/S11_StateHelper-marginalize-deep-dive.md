# StateHelper 边缘化 精读报告 (S11)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/StateHelper.cpp`
> **对应论文**: Mourikis & Roumeliotis 2007 MSCKF [20] (滑窗边缘化/删旧 clone)、Li & Mourikis 2013 IJRR [27] (时间偏移)、常规 Schur 补边缘化
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S11（Phase 3 函数级路线，阶段 C 状态边缘化）

> 上游：S9（协方差前推）、S10（滑窗增广）。本篇聚焦 **`StateHelper::marginalize` (:271)** 与两个策略封装 **`marginalize_old_clone` (:618)** / **`marginalize_slam` (:631)** —— MSCKF 滑窗的"收缩点"：删除最旧 clone（维持 `max_clone_size`）或已丢失的 SLAM 路标，压缩协方差并用 `set_local_id` 重排保持连续布局（S1 不变式）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
滑窗超限 / SLAM 路标丢失
        │
        ▼
marginalize_old_clone (S11 :618)  ─┐
marginalize_slam       (S11 :631)  ─┤→ marginalize(state, marg) @ :271 ★ 本篇核心
                                       │
                                       ▼ 压缩 _Cov + set_local_id 重排
                                  维持连续布局 (S1 不变式)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 预测/更新后的滑窗收缩 |
| 输入 | `marg`（要删的变量）；`marginalize_old_clone`/`marginalize_slam` 无参（从 state 取） |
| 输出 | 压缩 `state->_Cov`，重排 `_variables`，更新 `_clones_IMU`/`_features_SLAM` |
| 调用方 | `VioManager`（触发）、UpdaterSLAM（S14/S12 后续） |
| 被调用方 | `Type::set_local_id`（S1）、`State::margtimestep` |

### 1.2 一句话概括

`marginalize` 把单变量从 `_Cov` 中物理删除（保留其余变量对应块，重排 id 保持连续）；`marginalize_old_clone` 维持滑窗长度 `max_clone_size`（MSCKF 核心 [20]）；`marginalize_slam` 删除应边缘化的 SLAM 路标（保留 aruco 标签路标）。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp
void StateHelper::marginalize(state, std::shared_ptr<Type> marg);          // :271
void StateHelper::marginalize_old_clone(state);                            // :618
void StateHelper::marginalize_slam(state);                                 // :631
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 状态边缘化 `marginalize` —— 删除单变量并压缩协方差

#### 理论推导

**两种边缘化策略**（Anderson & Moore 1979）：

**(a) Schur 补边缘化**（图优化常用）：设 $\mathbf{x} = [\mathbf{x}_a^\top, \mathbf{x}_b^\top]^\top$，要边缘化 $\mathbf{x}_b$。信息矩阵 $\mathbf{\Lambda} = \mathbf{P}^{-1}$ 分块为：

$$\mathbf{\Lambda} = \begin{bmatrix} \mathbf{\Lambda}_{aa} & \mathbf{\Lambda}_{ab} \\ \mathbf{\Lambda}_{ba} & \mathbf{\Lambda}_{bb} \end{bmatrix} \tag{S11-1}$$

Schur 补边缘化后，$\mathbf{x}_a$ 的信息矩阵为：

$$\mathbf{\Lambda}_a' = \mathbf{\Lambda}_{aa} - \mathbf{\Lambda}_{ab}\mathbf{\Lambda}_{bb}^{-1}\mathbf{\Lambda}_{ba} \tag{S11-2}$$

特点：保留 $\mathbf{x}_b$ 对 $\mathbf{x}_a$ 的约束，但信息矩阵变稠密（fill-in）。

**(b) 直接删除**（MSCKF 采用）：直接保留 $\mathbf{x}_a$ 对应的协方差子块：

$$\mathbf{P}_a' = \mathbf{P}_{aa} \tag{S11-3}$$

特点：不保留 $\mathbf{x}_b$ 的约束，协方差保持稀疏。

**MSCKF 为什么可以直接删除**：MSCKF 的零空间投影（S13）已将特征观测的约束完全融入 IMU 状态。设 $\mathbf{x} = [\mathbf{x}_{\text{clone}}^\top, \mathbf{x}_{\text{rest}}^\top]^\top$，MSCKF 更新后的协方差为：

$$\mathbf{P}' = \begin{bmatrix} \mathbf{P}_{cc} & \mathbf{P}_{cr} \\ \mathbf{P}_{rc} & \mathbf{P}_{rr} \end{bmatrix} \tag{S11-4}$$

直接删除 $\mathbf{x}_{\text{clone}}$ 后得 $\mathbf{P}_{\text{rest}}' = \mathbf{P}_{rr}$，而 Schur 补结果为 $\mathbf{P}_{\text{rest}}^{\text{Schur}} = \mathbf{P}_{rr} - \mathbf{P}_{rc}\mathbf{P}_{cc}^{-1}\mathbf{P}_{cr}$。差异项 $\mathbf{P}_{rc}\mathbf{P}_{cc}^{-1}\mathbf{P}_{cr}$ 是 $\mathbf{x}_{\text{clone}}$ 对 $\mathbf{x}_{\text{rest}}$ 的"额外约束"——在 MSCKF 中，这部分约束已在零空间投影时被考虑（通过 $\mathbf{Q}_2^\top\mathbf{H}_x$），因此可直接删除（MSCKF 论文 Section III-B）。

**代码中的分块结构**：设被删变量 $\mathbf{x}_m$ 将 $\mathbf{P}$ 分为三区 $\mathbf{x} = [\mathbf{x}_1, \mathbf{x}_m, \mathbf{x}_2]$：

$$\mathbf{P} = \begin{bmatrix} \mathbf{P}_{11} & \mathbf{P}_{1m} & \mathbf{P}_{12} \\ \mathbf{P}_{m1} & \mathbf{P}_{mm} & \mathbf{P}_{m2} \\ \mathbf{P}_{21} & \mathbf{P}_{2m} & \mathbf{P}_{22} \end{bmatrix} \xrightarrow{\text{删除 }\mathbf{x}_m} \mathbf{P}' = \begin{bmatrix} \mathbf{P}_{11} & \mathbf{P}_{12} \\ \mathbf{P}_{21} & \mathbf{P}_{22} \end{bmatrix} \tag{S11-5}$$

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:271-339
void StateHelper::marginalize(state, marg) {
    // 1. 防御: marg 必须在 _variables 内 (不支持子变量边缘化!)
    if (find(_variables.begin(),_variables.end(),marg)==_variables.end())    // :274
      { ... std::exit; }                                                     // :275-278

    // 2. 计算三区尺寸
    int marg_size=marg->size(), marg_id=marg->id();                          // :293-294
    int x2_size = _Cov.rows() - marg_id - marg_size;                         // :295

    // 3. 构造新协方差 (去掉 marg 行/列)
    Eigen::MatrixXd Cov_new(_Cov.rows()-marg_size, _Cov.rows()-marg_size);   // :297
    Cov_new.block(0,0,        marg_id, marg_id) = _Cov.block(0,0,marg_id,marg_id);                 // :300 P11
    Cov_new.block(0,marg_id,  marg_id, x2_size) = _Cov.block(0,marg_id+marg_size,marg_id,x2_size); // :303 P12
    Cov_new.block(marg_id,0,  x2_size, marg_id) = Cov_new.block(0,marg_id,marg_id,x2_size).transpose(); // :306 P21
    Cov_new.block(marg_id,marg_id, x2_size,x2_size)=_Cov.block(marg_id+marg_size,marg_id+marg_size,x2_size,x2_size); // :309 P22

    state->_Cov = Cov_new;                                                   // :313
    assert(_Cov.rows()==Cov_new.rows());                                     // :315

    // 4. 重排剩余变量 id: 排在 marg 之后的变量向前挪 marg_size
    std::vector<Type> remaining;
    for (var : _variables) {                                                 // :320
      if (var != marg) {
        if (var->id() > marg_id)
          var->set_local_id(var->id() - marg_size);                          // :325 ★ 连续化重排
        remaining.push_back(var);
      }
    }
    marg->set_local_id(-1);                                                  // :335 标记已删
    state->_variables = remaining;                                           // :338
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{P}_{11}$ | `Cov_new.block(0,0,marg_id,marg_id)` | `:300` |
| $\mathbf{P}_{12}$ | `Cov_new.block(0,marg_id,...)` | `:303` |
| $\mathbf{P}_{21}=\mathbf{P}_{12}^\top$ | `Cov_new.block(marg_id,0,...).transpose()` | `:306` |
| $\mathbf{P}_{22}$ | `Cov_new.block(marg_id,marg_id,...)` | `:309` |
| id 重排（保持连续） | `var->set_local_id(var->id()-marg_size)` | `:325` |

> **设计要点（连续布局维护）**：边缘化会留下"空洞"（被删变量原来占的行/列）。`:325` 用 `set_local_id` 把所有排在后面的变量偏移 `-marg_size`，使协方差在**逻辑上重新连续**。这保证了 S9 的连续 `assert`（`:48`）始终成立。注释 `:318` 明确说**尚不支持子变量边缘化**——只能删顶层 `Type`（如整个 clone、整个 SLAM 路标）。

### 3.2 滑窗/SLAM 专用边缘化

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:618-645
void StateHelper::marginalize_old_clone(state) {
    if ((int)_clones_IMU.size() > _options.max_clone_size) {               // :619 滑窗容量超限
      double marginal_time = state->margtimestep();                        // :620 最旧 clone 时刻
      std::lock_guard<std::mutex> lock(state->_mutex_state);              // :622 线程安全
      marginalize(state, _clones_IMU.at(marginal_time));                   // :624 调 3.1
      _clones_IMU.erase(marginal_time);                                    // :627 删指针
    }
}

void StateHelper::marginalize_slam(state) {
    auto it=_features_SLAM.begin();
    while (it!=_features_SLAM.end()) {                                      // :636
      if (it->second->should_marg && (int)it->first > 4*max_aruco_features) { // :637 标记+非aruco
        marginalize(state, it->second);                                    // :638 调 3.1
        it = _features_SLAM.erase(it);                                      // :639 删指针
      } else it++;
    }
}
```

#### 对照注释

| 调用 | 触发条件 | 行号 |
|------|---------|------|
| `marginalize_old_clone` | `_clones_IMU.size() > max_clone_size` | `:619` |
| `marginalize_slam` | `should_marg==true` 且非 aruco 路标 | `:637` |

> **设计要点**：这两个是 `marginalize`（§3.1）的**策略封装**——决定"何时删谁"。`marginalize_old_clone` 维持滑窗长度 `max_clone_size`（MSCKF 核心：用固定数量历史 clone 提供多帧约束 [20]）；`marginalize_slam` 删除已丢失/应边缘化的 SLAM 路标，但**保留 aruco 标签路标**（`< 4*max_aruco_features` 的 id 段，见 `:637` 的保护条件）。

---

## Section 4: 函数调用链

```
VioManager / UpdaterSLAM
  ├─ marginalize_old_clone(state) @ :618  ★ S11
  │     └─ marginalize(state, oldest_clone) @ :624 → :271
  └─ marginalize_slam(state)      @ :631  ★ S11
        └─ marginalize(state, slam_feat) @ :638 → :271

StateHelper::marginalize(state, marg) @ :271  ★ 核心
  └─ var->set_local_id(var->id()-marg_size)  (S1 连续 id 重排)
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 |
|--------|-----------|------|
| `marginalize` | `StateHelper.cpp:271` | 删单变量+压缩协方差 |
| `marginalize_old_clone` | `StateHelper.cpp:618` | 删最旧 clone（滑窗） |
| `marginalize_slam` | `StateHelper.cpp:631` | 删 SLAM 路标 |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `max_clone_size` | yaml | `StateOptions` | 滑窗 | 最大历史 clone 数 | ★★★ 高 |
| `max_aruco_features` | yaml | `StateOptions` | SLAM | aruco 路标保护范围 | ★ 低 |

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `marginalize` 不支持子变量 | ★ 中 | 试图删 PoseJPL 内单元素 | 注释 `:318` 未实现 | 靠"删整个 clone/路标"规避 |
| 2 | 线程锁仅保护 SLAM 分支 | ★ 低 | 并发访问 `_clones_IMU` | `:622` 仅 `marginalize_old_clone` 加锁 | clone 访问侧需自保 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `marginalize` | VINS-Mono | OKVIS | 常规 MSCKF [20] |
|---------|------------------------|-----------|-------|-----------------|
| 边缘化方式 | 物理删块 + id 重排 | Schur 补 | Schur 补 | 删 clone |
| 滑窗维持 | `max_clone_size` | 无(图优化) | 有 | 有 |
| aruco 保护 | `:637` 条件保留 | 无 | 无 | 无 |

### 设计取舍分析

- 选择 **物理删块 + `set_local_id` 重排** 是因为：EKF 框架中边缘化后变量不再需要，直接压缩协方差保持连续布局（S1 不变式），避免矩阵维度无限增长。
- 选择 **`max_clone_size` 固定滑窗** 是因为：MSCKF [20] 用固定数量历史 clone 提供多帧约束，长度决定精度/算力权衡。
- 选择 **保留 aruco 路标** 是因为：标签路标是全局已知基准，不应被边缘化（`:637` 的 `> 4*max_aruco_features` 保护）。
