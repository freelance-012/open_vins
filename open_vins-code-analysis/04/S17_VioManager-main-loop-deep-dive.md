# VioManager 主循环 精读报告 (S17)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/core/VioManager.cpp` / `VioManager.h`
> **对应论文**: 全系统编排（串联 B→C→D 阶段）；MSCKF [32]、ZUPT、初始化 [11]
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S17（Phase 3 函数级路线，阶段 E 主循环编排）

> 上游：A–D 全阶段（S1–S16）。本篇聚焦 **`VioManager`** 三大入口——`feed_measurement_imu` (:166)、`track_image_and_update` (:256)、`do_feature_propagate_update` (:323)——它把 B（预测）+ C（增广/边缘化）+ D（MSCKF/SLAM/ZUPT 更新）串成一条完整流水线，并负责**特征分组**（feats_lost / feats_marg / feats_slam）与滑动窗口调度。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
外部回调 (ROS/Sim) → feed_measurement_imu / feed_measurement_camera
        │
        ▼
┌──────────────── VioManager (中枢) ────────────────┐
│ feed_measurement_imu (:166)                        │
│   ├─ propagator->feed_imu (缓存)                   │
│   ├─ initializer->feed_imu (未初始化时)            │
│   └─ updaterZUPT->feed_imu (ZUPT 时)               │
│                                                     │
│ track_image_and_update (:256)  [相机回调]           │
│   ├─ trackFEATS->feed_new_camera (KLT 跟踪)        │
│   ├─ trackARUCO->feed_new_camera (Aruco)          │
│   ├─ updaterZUPT->try_update (:294) 静止检测       │
│   ├─ try_to_initialize (:311) 未初始化时           │
│   └─ do_feature_propagate_update (:320)            │
│                                                     │
│ do_feature_propagate_update (:323)  [核心编排]      │
│   ├─ propagator->propagate_and_clone (:340) → S5   │
│   ├─ 特征分组 feats_lost/marg/slam/maxtracks       │
│   ├─ updaterMSCKF->update (:525) → S13             │
│   ├─ updaterSLAM->update (:533) → S14             │
│   ├─ updaterSLAM->delayed_init (:547) → S15       │
│   ├─ updaterSLAM->change_anchors (:585) → S14     │
│   └─ StateHelper::marginalize_old_clone (:596)→S11 │
└─────────────────────────────────────────────────────┘
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 5 环节（主循环编排，P1-1，枢纽） |
| 输入 | IMU 消息（`feed_measurement_imu`）、相机消息（`track_image_and_update`） |
| 输出 | 改写 `state`（克隆增删、协方差、路标）、特征库清理、触发各下游更新 |
| 调用方 | ROS/Sim 回调层（G1） |
| 被调用方 | `Propagator`（S5）、`UpdaterMSCKF`（S13）、`UpdaterSLAM`（S14/S15）、`UpdaterZeroVelocity`（S16）、`StateHelper`（S9/S10/S11/S12）、`InertialInitializer`（F1） |

---

## Section 2: 三大函数骨架

```
feed_measurement_imu(message)  @ :166
├─ oldest_time = state->margtimestep() (未初始化时改窗口)     :170-176
├─ propagator->feed_imu(message, oldest_time)                :177
├─ 未初始化 → initializer->feed_imu(...)                     :180-182
└─ 已初始化+ZUPT → updaterZUPT->feed_imu(...)                :186-188

track_image_and_update(message)  @ :256
├─ trackFEATS->feed_new_camera (KLT)                          :281
├─ trackARUCO->feed_new_camera (若有)                        :286-288
├─ updaterZUPT->try_update(state, ts) (:294) → 命中则清 IMU+return
├─ 未初始化 → try_to_initialize (:311) → 失败 return
└─ do_feature_propagate_update(message)                      :320

do_feature_propagate_update(message)  @ :323
├─ 乱序 → return                                            :330-334
├─ propagator->propagate_and_clone(state, ts)                 :340  → S5
├─ 克隆数 < min(max_clone_size,5) → return (等待可三角化)     :348-352
├─ 特征分组 (见 Section 3)                                    :368-500
├─ updaterMSCKF->update (:525) / propagate invalidate_cache  :526
├─ updaterSLAM->update (分批)                                 :533-544  → S14
├─ updaterSLAM->delayed_init (新特征)                         :547  → S15
├─ retriangulate_active_tracks (base cam 0)                   :555-558
├─ 标记 to_delete + 双库 cleanup                              :563-582
├─ updaterSLAM->change_anchors                               :585  → S14
├─ cleanup_measurements(margtimestep)                         :588-593
└─ StateHelper::marginalize_old_clone(state)                  :596  → S11
```

---

## Section 3: 特征分组逻辑（核心）

`do_feature_propagate_update` 内部把特征库里的特征分成几组，决定各自走哪条更新路径：

| 组 | 来源 | 行号 | 去向 |
|----|------|------|------|
| `feats_lost` | `features_not_containing_newer(state->_timestamp, false, true)` —— 当前帧丢失的 KLT 特征 | :368-369 | MSCKF 更新（已无新观测，可边缘化） |
| `feats_marg` | `features_containing(margtimestep,...)`，仅当 `clones > max_clone_size \|\| clones > 5` | :372-377 | MSCKF 更新 |
| `feats_maxtracks` | 遍历 `feats_marg`，若某相机 `timestamps.size() > max_clone_size`（:417）即移入 —— 轨迹达上限，候选升级 SLAM | :411-429 | 部分进 SLAM |
| `feats_slam` | `max_slam_features>0` 且超时 `dt_slam_delay` 且未达上限时，从 `feats_maxtracks` 尾部取 `valid_amount` 个（:450）；以及 `_features_SLAM` 中仍在跟踪的轨迹（:464/:468） | :442-476 | SLAM 更新 / delayed_init |
| `feats_slam_UPDATE` | `feats_slam` 中已在 `state->_features_SLAM` 的 | :484-495 | `UpdaterSLAM::update`（S14） |
| `feats_slam_DELAYED` | `feats_slam` 中的新特征 | :484-495 | `UpdaterSLAM::delayed_init`（S15） |
| `featsup_MSCKF` | `feats_lost + feats_marg + feats_maxtracks`（非 SLAM 部分） | :498-500 | `UpdaterMSCKF::update`（S13） |

**SLAM 路标修剪**：遍历 `state->_features_SLAM`，若轨迹丢失（`feat2==nullptr`）或 `update_fail_count>1`，标记 `landmark->should_marg`（:473/:475），随后 `StateHelper::marginalize_slam(state)`（:481，S11）剔除非 Aruco 路标。

---

## Section 4: 接口与数据流

### 4.1 三个入口签名

```cpp
void feed_measurement_imu(const ov_core::ImuData &message);                       // :166 (VioManager.h:80)
void feed_measurement_camera(const ov_core::CameraData &message);                 // .h:81 → 内部调 track_image_and_update
void track_image_and_update(const ov_core::CameraData &message_const);            // :256
void do_feature_propagate_update(const ov_core::CameraData &message);             // :323
```

### 4.2 关键状态变量

| 变量 | 改写方 | 说明 |
|------|--------|------|
| `state` | propagate_and_clone / 各 updater / marginalize | 克隆增删、误差态、路标 |
| `state->_features_SLAM` | delayed_init（增）、marginalize_slam（删） | 持久路标集合 |
| `state->_clones_IMU` | augment_clone（增）、marginalize_old_clone（删） | 滑窗 |
| `feature database` | cleanup / cleanup_measurements / to_delete | 测量生命周期 |
| `has_moved_since_zupt` / `good_features_MSCKF` | do_feature_propagate_update | ZUPT 与统计 |

---

## Section 5: 数学与公式对照

本篇为**编排层**，无新数学；所有公式见上游：
- 预测/克隆 → S5/S10
- 协方差更新 → S9
- MSCKF 残差/零空间 → S13
- SLAM 残差/锚点 → S14/S15
- ZUPT 残差 → S16
- 初始化 → F1–F3

核心"不变式"：每帧图像 → 先 `propagate_and_clone` 把状态时间推到该帧并长出新克隆 → 再把"已无新观测"的特征用 MSCKF 方式边缘化掉 → 把"轨迹够长"的特征升格 SLAM → SLAM 路标实时更新 → 滑窗最旧克隆边缘化。这构成 MSCKF+SLAM 滑动窗口的闭环。

---

## Section 6: 关键实现细节与易错点

1. **ZUPT 短路**（`:294-306`）：若 `try_update` 命中（判定静止），直接清理 IMU 并返回，**不跑后续 MSCKF/SLAM 更新**——因为静止时本帧不动，无需视觉更新。
2. **未初始化分支**（`:310-317`）：先 `try_to_initialize`，失败则 return，不进入主循环。初始化成功后才正式跑 `do_feature_propagate_update`。
3. **克隆数门槛**（`:348-352`）：`clones < min(max_clone_size, 5)` 时 return，等滑窗攒够帧数才能三角化。
4. **特征分组顺序敏感**：`feats_lost` 要先剔除不含当前相机、且已被 `feats_marg` 包含的（:382-408），避免重复更新。
5. **SLAM 分批更新**（`:533-544`）：`max_slam_in_update` 限制单次 SLAM 更新特征数，防止单帧计算爆炸。
6. **base 相机 0 才重三角化**（`:555-558`）：`retriangulate_active_tracks` 只在 base 相机时调用。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **下游**：S5（propagate_and_clone）、S9（EKFPropagation/EKFUpdate）、S10/S11（augment/marginalize）、S12（initialize）、S13（MSCKF update）、S14（SLAM update + change_anchors）、S15（delayed_init）、S16（ZUPT try_update）。
- **上游**：G1（ROS/Sim 回调把消息喂入 `feed_measurement_imu`/`feed_measurement_camera`）、F1（try_to_initialize 触发初始化）。

### 待详细补充项

- **S17a — `try_to_initialize`**：如何调 `InertialInitializer::initialize`（F1）并把结果（IMU 状态、clones、SLAM 特征）装回 `state`（VioManager.cpp:311 附近，实现在 `VioManagerHelper.cpp`）。
- **S17b — `retriangulate_active_tracks`**：base 相机重三角化逻辑（:555-558）。
- **S17c — `cleanup_measurements` / 双库 cleanup**：测量回收时机（:579-593）。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"）。
