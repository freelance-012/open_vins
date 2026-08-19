# InertialInitializer 初始化调度 精读报告 (F1)

> **仓库**: open_vins
> **模块路径**: `ov_init/src/init/InertialInitializer.cpp` / `.h`
> **对应论文**: 初始化策略调度（静态 [静止] / 动态 [任意运动]，Dong-Si & Mourikis 2012 背景）；对齐 [11]
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 F1（Phase 3 函数级路线，阶段 F 初始化，调度层）

> 上游：S17（`try_to_initialize` 调用本层）。本篇聚焦 **`InertialInitializer`**——它**不实现**具体初始化算法，而是**调度**两种策略：`StaticInitializer`（静止启动，纯 IMU）与 `DynamicInitializer`（任意运动，IMU+视觉+CPI+Ceres）。判定依据是**视差检测**（窗口前后两半的运动状态）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
VioManager::try_to_initialize  (S17)
        │
        ▼
InertialInitializer::initialize(timestamp, covariance, order, t_imu, wait_for_jerk)  @ :73  ★ 本篇
   ├─ 0. 求最新相机时间 + 裁剪窗口                            :77-96
   ├─ 1. 视差检测: 窗口分两半, compute_disparity            :98-125
   ├─ 2. 策略判定树 (has_jerk / is_still)                    :130-146
   │      ├─ 静态 → init_static->initialize(...)             :132-134
   │      ├─ 动态 → init_dynamic->initialize(...,_clones,_features_SLAM)  :135-139
   │      └─ 其余 → return false (等待/纯惯导不初始化)
   └─ (下游) StaticInitializer (F-stat) / DynamicInitializer (F2)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 6 环节（初始化调度，F 阶段入口） |
| 输入 | `imu_data`（IMU 缓冲）、`_db`（特征库）、`oldest_time`（滑窗裁剪） |
| 输出 | 引用出参 `timestamp`/`covariance`/`order`/`t_imu`（静态）；动态额外出 `_clones_IMU`/`_features_SLAM` |
| 调用方 | `VioManager::try_to_initialize`（S17） |
| 被调用方 | `StaticInitializer::initialize`（F-static）、`DynamicInitializer::initialize`（F2）、`FeatureHelper::compute_disparity` |

---

## Section 2: 函数骨架

```
InertialInitializer (构造)  :38-47
├─ 建共享 imu_data 缓冲、_db
└─ 构造 init_static (StaticInitializer) + init_dynamic (DynamicInitializer)，三者共享 imu_data/_db

feed_imu(message, oldest_time=-1)  :49
├─ imu_data.emplace_back(message)                       :52
└─ oldest_time!=-1 → 删 timestamp<oldest_time 旧数据    :61-70  (滑窗裁剪)

initialize(timestamp, covariance, order, t_imu, wait_for_jerk)  :73
├─ 最新相机时间 newest_cam_time + oldest_time 裁剪      :77-96
├─ 视差检测 (init_max_disparity>0 时):                  :98-125
│     窗口中点分两半 → compute_disparity 得 avg_disp0/avg_disp1
│     feat<15 → return false (特征不足)
│     moving_1to0 / moving_2to1 (两半是否超阈值)
├─ 策略树:                                             :130-146
│     has_jerk = !moving_1to0 && moving_2to1
│     is_still = !moving_1to0 && !moving_2to1
│     ├─ 静态: (has_jerk&&wait_for_jerk)||(is_still&&!wait_for_jerk) && init_imu_thresh>0
│     │        → init_static->initialize(...)
│     ├─ 动态: init_dyn_use && !is_still
│     │        → init_dynamic->initialize(...,_clones_IMU,_features_SLAM)
│     └─ 其余 → return false
```

---

## Section 3: 逐段精读

### 3.1 构造函数 (`:38-47`)

```cpp
imu_data = std::make_shared<std::vector<ov_core::ImuData>>();  // :42 共享缓冲
init_static = std::make_shared<StaticInitializer>(params, imu_data, _db);   // :45
init_dynamic = std::make_shared<DynamicInitializer>(params, imu_data, _db); // :46
```
- 三个初始化器**共享同一份 `imu_data` 与 `_db`**，避免数据拷贝。

### 3.2 feed_imu (`:49-70`)

```cpp
imu_data->emplace_back(message);                       // :52 追加
if (oldest_time != -1) {                               // :61 滑窗裁剪
    while (!imu_data->empty() && (*imu_data)[0].timestamp < oldest_time)
        imu_data->erase(imu_data->begin());
}
```
- 仅做缓冲追加 + 过期裁剪，无算法。

### 3.3 视差检测 (`:98-125`)

```cpp
double mid_time = oldest_time + (newest_cam_time - oldest_time) * 0.5;  // :106 窗口中点
FeatureHelper::compute_disparity(_db, oldest_time, mid_time, avg_disp0, var0, n0);   // :111 旧半窗
FeatureHelper::compute_disparity(_db, mid_time, newest_cam_time, avg_disp1, var1, n1); // :112 新半窗
if (n0 < 15 || n1 < 15) return false;                  // :116-119 特征不足
moving_1to0 = (avg_disp0 > init_max_disparity);        // :123
moving_2to1 = (avg_disp1 > init_max_disparity);        // :124
```
- **关键思想**：把初始化窗口分成前后两半，分别算平均视差。若前静后动 → "急动"（jerk）；若两半皆静 → 静止；若前动后动 → 一直在动。
- 最少特征数硬编码 `feat_thresh=15`（`:115`）。

### 3.4 策略判定树 (`:130-146`)

```cpp
bool has_jerk = (!moving_1to0 && moving_2to1);   // :130 旧静新动
bool is_still = (!moving_1to0 && !moving_2to1);   // :131 两半皆静

if (((has_jerk && wait_for_jerk) || (is_still && !wait_for_jerk)) && params.init_imu_thresh > 0.0) {
    return init_static->initialize(timestamp, covariance, order, t_imu, wait_for_jerk);  // :132-134
}
if (params.init_dyn_use && !is_still) {
    // 建局部 _clones_IMU / _features_SLAM 传出
    return init_dynamic->initialize(timestamp, covariance, order, t_imu, _clones_IMU, _features_SLAM);  // :135-139
}
// 其余: 打印原因, return false
return false;   // :140-146
```
- **静态分支**：刚经历 jerk 后静止（且 `wait_for_jerk`），或整体静止（不要求 jerk）。走 `StaticInitializer`。
- **动态分支**：启用动态（`init_dyn_use`）且非静止。走 `DynamicInitializer`，额外传出 clones 与 SLAM 特征（用于初始化后直接进入 SLAM 状态）。
- **其余**：既不静止、又没开动态 → 不初始化，等条件满足。

### 3.5 默认阈值 (`InertialInitializerOptions.h`)

| 参数 | 默认 | 含义 |
|------|------|------|
| `init_window_time` | 1.0 | 初始化窗口时长（秒） |
| `init_imu_thresh` | 1.0 | IMU 静态判定阈值（>0 才允许静态初始化） |
| `init_max_disparity` | 1.0 | 视差检测阈值 |
| `init_max_features` | 50 | 最大特征数 |
| `init_dyn_use` | false | 是否启用动态初始化 |
| `init_dyn_num_pose` | 5 | 动态初始化使用的相机帧数 |
| `init_dyn_min_deg` | 45.0 | 最小角运动（度） |

---

## Section 4: 接口与数据结构

### 4.1 签名

```cpp
bool InertialInitializer::initialize(double &timestamp, Eigen::MatrixXd &covariance,
                                     std::vector<std::shared_ptr<Type>> &order,
                                     std::shared_ptr<ov_type::IMU> t_imu, bool wait_for_jerk = true);  // :73 (.h:97)
```

| 参数 | 方向 | 说明 |
|------|------|------|
| `timestamp` | out | 初始化时刻 |
| `covariance` | out | 初始协方差（IMU 块） |
| `order` | out | 协方差块顺序（含 `t_imu`） |
| `t_imu` | out | 初始 IMU 状态（q/p/v/bg/ba） |
| `wait_for_jerk` | in | 是否等"jerk 后静止" |

- **无 `InitResults` 结构体**：结果全部经引用出参回传（动态分支额外经 `_clones_IMU`/`_features_SLAM` 成员传出）。

---

## Section 5: 数学与公式对照

本篇为**调度层**，无新数学。核心判定是启发式（视差阈值），不涉及公式推导。
- 静态初始化算法 → 见 StaticInitializer（F-stat，在 F2 文档中覆盖）
- 动态初始化算法（CPI + Ceres MLE）→ 见 F2 / F3

---

## Section 6: 关键实现细节与易错点

1. **共享缓冲**：`imu_data`/`_db` 三器共享，注意 `feed_imu` 的裁剪会影响下游读取的窗口（`:61-70`）。
2. **视差分两半**：用窗口中点 `mid_time` 切分（`:106`），不是简单前后帧。
3. **`wait_for_jerk` 语义**：静态初始化可要求"先动后静"（jerk）以确认系统已稳定，也可宽松为"整体静止即初始化"（关 `wait_for_jerk`）。
4. **动态默认关闭**（`init_dyn_use=false`）：默认只走静态路径。开启动态需设 `init_dyn_use=true`。
5. **返回 false 的歧义**：本层 `return false` 可能意味着"尚未满足条件"或"不支持该运动模式"，上层（S17）需区分是"继续等"还是"初始化失败"。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **上游**：S17（`try_to_initialize` 调用 `initialize`）。
- **下游**：`StaticInitializer`（F-stat，见 F2）、`DynamicInitializer`（F2）、`CpiV1`（F3）、Ceres 因子（F2）。

### 待详细补充项

- **F1a — 静态初始化内部算法**（重力方向 gram_schmidt、bias 均值估计、速度置零、协方差构造）：见 F2 文档中 StaticInitializer 部分，或后续拆 F1a。
- **F1b — 动态初始化内部算法**（线性系统求 v/g + Ceres MLE）：见 F2。
- **F1c — `compute_disparity` 实现**：见 S16a（UpdaterZeroVelocity 已引用）。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"）。
