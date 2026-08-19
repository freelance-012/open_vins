# 相机模型 工程支撑 精读报告 (G2)

> **仓库**: open_vins
> **模块路径**: `ov_core/src/cam/{CamBase,CamRadtan,CamEqui}.h/.cpp`
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 G2（Phase 3 函数级路线，阶段 G 工程支撑，P3-1）

> 本篇为**工程支撑层**（了解投影/反投影接口即可）。相机模型是 MSCKF/SLAM 特征重投影残差的基础（被 `UpdaterHelper`、`State` 调用），位于"状态表示 + 测量更新"环节。

---

## Section 1: 基类接口（CamBase）

`ov_core/src/cam/CamBase.h` 定义纯虚接口（注意命名：**`distort`/`undistort`**，非 project/unproject）：

```cpp
Eigen::Vector2f distort_f(const Eigen::Vector2f& uv_norm);    // 归一化坐标 → 像素 (投影)
Eigen::Vector2f undistort_f(const Eigen::Vector2f& uv_dist);  // 像素 → 归一化坐标 (反投影)
void compute_distort_jacobian(...);
// 另有 distort_d / undistort_d (double 封装)
```

## Section 2: 两个实现

| 模型 | 文件 | 畸变类型 | 说明 |
|------|------|----------|------|
| `CamRadtan` | `CamRadtan.h` | radtan / Brown-Conrady（pinhole-radtan） | `undistort_f` 调 OpenCV `cv::undistortPoints`；`distort_f` 手算 `r²/r⁴` 径向 + 切向 `[k1,k2,p1,p2]` |
| `CamEqui` | `CamEqui.h` | fisheye / equidistant（pinhole-equi） | `undistort_f` 调 OpenCV `cv::fisheye::undistortPoints`；`distort_f` 用 `θ=atan(r)` + `θ_d=θ+k1θ³+...+k4θ⁹` |

## Section 3: 工厂选择

`VioManagerOptions.h` 读 YAML 字段 `distortion_model`（默认 `"radtan"`）：
- `"equidistant"` → `make_shared<CamEqui>`（第 273 行）
- 否则 → `make_shared<CamRadtan>`（第 276 行）
- 存入 `params.camera_intrinsics`（map<int, shared_ptr<CamBase>>，按相机 id 索引）。

## Section 4: 在 pipeline 中的位置

```
State (_cam_intrinsics_cameras) + UpdaterHelper::distort_d (:343) → 重投影残差 (S13/S14)
VioManagerHelper::undistort_cv → 跟踪/特征初始化
```
- 是 MSCKF/SLAM 特征投影残差的几何基础（见 S13a `get_feature_jacobian_full` 内部投影调用）。

## Section 5: 待详细补充项

- **G2a — `CamRadtan::distort_f` 畸变公式**：手算 `r²/r⁴` + 切向项推导。
- **G2b — `CamEqui::distort_f` 等距模型**：`θ=atan(r)` 反函数实现。

> 工程支撑层，按需查。遵循 SKILL.md "04/ 文档拆分规则"，子文档待需要时补。
