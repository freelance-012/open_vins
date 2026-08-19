# ROS I/O 工程支撑 精读报告 (G1)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/ros/{ROS1Visualizer,ROS2Visualizer,ROS1Subscriber,ROS2Subscriber}.cpp` + `ROSVisualizerHelper.cpp`
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 G1（Phase 3 函数级路线，阶段 G 工程支撑，P2-1）

> 本篇为**工程支撑层**（只需了解话题名/参数/回调入口，无需逐行推导）。它是最外层 I/O 适配：把 `sensor_msgs` 转成 `ov_core::ImuData`/`CameraData`，喂给 `VioManager`（S17）；并把 `VioManager` 的结果可视化发布。

---

## Section 1: 模块定位

ROS1 与 ROS2 两套实现，**话题名完全一致**。核心逻辑在 `ROSVisualizerHelper.cpp` 共用，各版本包一层 publisher/subscriber。

## Section 2: Publish 话题（两版同名）

| 话题名 | 类型 | 发布函数 |
|--------|------|----------|
| `poseimu` | PoseWithCovarianceStamped | `publish_state()` |
| `odomimu` | Odometry | `visualize_odometry()` |
| `pathimu` | Path | `publish_state()` |
| `points_msckf` / `points_slam` / `points_aruco` / `points_sim` | PointCloud2 | `publish_features()` |
| `trackhist` | Image (image_transport) | `publish_images()` |
| `posegt` / `pathgt` | PoseStamped / Path | `publish_groundtruth()` |
| `loop_pose` / `loop_extrinsic` | Odometry | `publish_loopclosure_information()` |
| `loop_feats` | PointCloud | 同上 |
| `loop_intrinsics` | CameraInfo | 同上 |
| `loop_depth` / `loop_depth_colored` | Image | 同上 |

## Section 3: Subscribe 话题

| 话题 | 默认 | 回调 |
|------|------|------|
| IMU | `/imu0` (`topic_imu`) | `callback_inertial` → `_app->feed_measurement_imu()` |
| 图像（每相机） | `/cam{i}/image_raw` (`topic_camera{i}`) | `callback_stereo`（双目，message_filters 同步）/ `callback_monocular`（单目）→ `_app->feed_measurement_camera()` |

- 相机标定**不在 ROS 订阅**，而是从 YAML（`relative_config_imucam`）解析（见 G2）。

## Section 4: 在 pipeline 中的位置

```
sensor_msgs → ROS{1,2}Subscriber::callback_* → VioManager::feed_measurement_imu / _camera (S17)
VioManager 结果 → ROS{1,2}Visualizer::publish_* → rviz 话题
```

## Section 5: 接口要点

- 三个主要外部回调：`callback_inertial`、`callback_stereo`、`callback_monocular`。
- 构造函数里 `advertise`/`create_publisher` 注册上述话题。
- 注意：ROS 层**不实现算法**，只做消息转换与可视化。

## Section 6: 待详细补充项

- **G1a — `ROSVisualizerHelper` 共用逻辑**：位姿/特征/回环信息的实际拼装。
- **G1b — 回环话题**（`loop_*`）：若启用 VIO 回环模块，相关发布逻辑。

> 工程支撑层，按需查。遵循 SKILL.md "04/ 文档拆分规则"，子文档待需要时补。
