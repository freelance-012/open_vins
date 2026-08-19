# CPI 在线侧对照 工程支撑 精读报告 (G3)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/Propagator.cpp`（在线，未 include cpi）
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 G3（Phase 3 函数级路线，阶段 G 工程支撑，P3-2 对照）

> 本篇为**对照文档**（不新增代码，只澄清关系）。`ov_core/src/cpi/CpiV1.*`（F3）具有**双重身份**：
> ① 在线 `Propagator` **未 include cpi**，仅借鉴其连续时间预积分思想自实现（见 S5–S8 的离散 RK4/ACI2）；
> ② 动态初始化 `DynamicInitializer` **直接调用 `CpiV1`** 并接入 Ceres MLE（属 F2/F3 核心实现，已上移计入 P0-6）。

---

## Section 1: 在线 Propagator 与 CPI 的关系

- **在线 `Propagator`（`ov_msckf/src/state/Propagator.cpp`）**：实现 EKF 预测（S5–S8），采用**离散**积分（RK4 / ACI2 / 离散 F/G），状态含 bias 直接更新，**不** include `CpiV1`。
- **`CpiV1`（F3）**：连续时间预积分，输出相对锚帧的测量 + 偏置雅可比，供**批量 MLE**（Ceres）使用。

## Section 2: 思想对照（同源异实现）

| 维度 | 在线 Propagator (S5–S8) | CpiV1 (F3) |
|------|------|------|
| 积分方式 | 离散 RK4 / ACI2 / 离散 F-G | 连续解析积分（`f_1..f_4` 多项式系数） |
| 协方差 | 离散 `ΦPΦᵀ+Q_d` | 连续 Riccati 的 RK4 数值解 |
| 参考帧 | 全局绝对状态 | 固定锚帧相对测量 |
| 偏置 | 状态变量，直接更新 | 一阶雅可比修正，不重积分 |
| 用途 | EKF 逐帧预测 | 初始化批量优化残差 |

- 两者共享同一套 IMU 运动学（`skew_x`、`rot_2_quat`、`[θ,b_w,v,b_a,p]` 15 维误差排序、相同的 `F`/`G` 结构）。

## Section 3: 何时回查 CPI

- 精读 **P0-1（Propagator, S5–S8）** 时若遇预积分公式分歧，可回查 **F3（CpiV1）** 比对连续时间思想来源。
- 精读 **P0-6（初始化, F1–F3）** 时，CPI 是核心实现，必须读 F3。

## Section 4: 结论

在线滤波器与离线初始化**共用 IMU 预积分的数学本质**，但为实现目标（EKF 流式更新 vs 批量 MLE）选择了不同的数值路径。本篇仅作对照索引，无新代码分析。

> 遵循 SKILL.md "04/ 文档拆分规则"，对照类文档不另拆子文档。
