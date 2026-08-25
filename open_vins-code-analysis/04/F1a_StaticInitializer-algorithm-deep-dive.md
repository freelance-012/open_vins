# 静态初始化 StaticInitializer 算法深度解析 (F1a)

> **生成日期**：2026-08-25
> **前置文档**：`F1_InertialInitializer-deep-dive.md`（调度层）、`F2_Static-Dynamic-Initializer-deep-dive.md`（算法层概览）
> **核心代码**：`ov_init/src/static/StaticInitializer.cpp`
> **代码版本**：`17b73cfe4b870ade0a65f9eb217d8aab58deae19`
> **理论依据**：Trawny & Roumeliotis "Indirect Kalman Filter for 3D Attitude Estimation" 2005 Section II-A

---

## 1. 模块定位与职责

### 1.1 在数据流中的位置

```
InertialInitializer::initialize()  (F1)
        │  策略判定: has_jerk / is_still
        ▼
StaticInitializer::initialize(timestamp, covariance, order, t_imu, wait_for_jerk)  ★ 本篇
   ├─ 1. 窗口划分与统计量计算
   ├─ 2. 静止判定（加速度方差 + jerk 检测）
   ├─ 3. Gram-Schmidt 重力对齐求姿态
   ├─ 4. 零偏估计（陀螺均值 + 加表减重力）
   ├─ 5. 设置 IMU 状态（速度=0，位置=0）
   └─ 6. 构造初始协方差
```

| 属性 | 内容 |
|------|------|
| **数据流位置** | 静态初始化算法实现 |
| **输入** | `imu_data`（IMU 缓冲）、`wait_for_jerk`（是否等急动） |
| **输出** | `timestamp`（初始化时刻）、`covariance`（初始协方差）、`order`（变量顺序）、`t_imu`（IMU 状态） |
| **调用方** | `InertialInitializer`（F1） |

### 1.2 一句话概括

`StaticInitializer` 假设载体**静止启动**，利用静止期间加速度计测量对齐重力方向求姿态，用陀螺仪均值估计零偏，速度/位置设为零（静止假设），构造对角初始协方差。

---

## 2. 数学模型

### 2.1 静止假设下的 IMU 模型

静止时，载体加速度为零，角速度为零。IMU 测量方程简化为：

$$\mathbf{a}_m = \mathbf{R}_{G\to I}\,\mathbf{g}^G + \mathbf{b}_a + \mathbf{n}_a \tag{F1a-1}$$

$$\boldsymbol{\omega}_m = \mathbf{b}_\omega + \mathbf{n}_\omega \tag{F1a-2}$$

其中 $\mathbf{g}^G = [0,\,0,\,g]^\top$，$g = 9.81\,\text{m/s}^2$。

### 2.2 姿态求解

对式 (F1a-1) 取时间平均（消除噪声）：

$$\bar{\mathbf{a}} \approx \mathbf{R}_{G\to I}\,\mathbf{g}^G + \mathbf{b}_a \tag{F1a-3}$$

由于静止时 $\mathbf{b}_a$ 是常值偏差，而 $\mathbf{R}_{G\to I}\,\mathbf{g}^G$ 是重力在机体系中的投影，两者方向不同（重力沿 $-z_G$，bias 是传感器固有偏差）。但在**短时间窗口内**，bias 变化可忽略，故：

$$\bar{\mathbf{a}} - \mathbf{b}_a \approx \mathbf{R}_{G\to I}\,\mathbf{g}^G \tag{F1a-4}$$

**关键洞察**：$\mathbf{R}_{G\to I}\,\mathbf{g}^G$ 的方向就是重力在 IMU 系中的方向。设 $\mathbf{z}_{\text{axis}} = \frac{\bar{\mathbf{a}} - \mathbf{b}_a}{\|\bar{\mathbf{a}} - \mathbf{b}_a\|}$，则 $\mathbf{z}_{\text{axis}}$ 是 $R_{G\to I}$ 的第三列（z 轴）。

但 $\mathbf{b}_a$ 未知！代码中**先忽略 bias**，直接用 $\bar{\mathbf{a}}$ 的方向作为重力方向：

$$\mathbf{z}_{\text{axis}} \approx \frac{\bar{\mathbf{a}}}{\|\bar{\mathbf{a}}\|} \tag{F1a-5}$$

这是因为初始化窗口很短（~1s），bias 相对于重力（9.81）很小（通常 <0.1），方向误差可忽略。

### 2.3 Gram-Schmidt 正交化

已知 $\mathbf{z}_{\text{axis}}$，需要构造完整的旋转矩阵 $\mathbf{R}_{G\to I} = [\mathbf{x}_{\text{axis}}\;\; \mathbf{y}_{\text{axis}}\;\; \mathbf{z}_{\text{axis}}]$。

**问题**：只有 z 轴方向确定，x/y 轴（即 yaw 角）是**不可观**的（VIO 4-DoF 不可观性之一）。

**Gram-Schmidt 方法**（[helper.h:138-171](src/open_vins/ov_init/src/utils/helper.h#L138-L171)）：选择与 $\mathbf{z}_{\text{axis}}$ 夹角最大的标准基向量 $\mathbf{e}_i$（$\mathbf{e}_1=[1,0,0]^\top$ 或 $\mathbf{e}_2=[0,1,0]^\top$），构造正交基：

$$\mathbf{x}_{\text{axis}} = \frac{\mathbf{z}_{\text{axis}} \times \mathbf{e}_i}{\|\mathbf{z}_{\text{axis}} \times \mathbf{e}_i\|} \tag{F1a-6}$$

$$\mathbf{y}_{\text{axis}} = \mathbf{z}_{\text{axis}} \times \mathbf{x}_{\text{axis}} \tag{F1a-7}$$

**选择 $\mathbf{e}_i$ 的策略**：选与 $\mathbf{z}_{\text{axis}}$ 夹角最大的（即 $|\mathbf{e}_i \cdot \mathbf{z}_{\text{axis}}|$ 最小的），避免叉积结果接近零向量（数值不稳定）。

$$i = \arg\min_{j \in \{1,2\}} |\mathbf{e}_j \cdot \mathbf{z}_{\text{axis}}| \tag{F1a-8}$$

代码实现：

```cpp
double inner1 = e_1.dot(z_axis) / z_axis.norm();
double inner2 = e_2.dot(z_axis) / z_axis.norm();
if (fabs(inner1) < fabs(inner2)) {
    x_axis = z_axis.cross(e_1);  // e_1 与 z 轴夹角更大
} else {
    x_axis = z_axis.cross(e_2);  // e_2 与 z 轴夹角更大
}
```

### 2.4 零偏估计

**陀螺零偏**：静止时角速度为零，测量值即为 bias：

$$\mathbf{b}_\omega = \bar{\boldsymbol{\omega}} = \frac{1}{N}\sum_{k=1}^N \boldsymbol{\omega}_m^{(k)} \tag{F1a-9}$$

**加表零偏**：由式 (F1a-4) 解出：

$$\mathbf{b}_a = \bar{\mathbf{a}} - \mathbf{R}_{G\to I}\,\mathbf{g}^G \tag{F1a-10}$$

其中 $\mathbf{R}_{G\to I}$ 由 Gram-Schmidt 构造（式 (F1a-6)(F1a-7)）。

### 2.5 速度与位置

静止假设下：

$$\mathbf{v}_{I_0}^{G} = \mathbf{0} \tag{F1a-11}$$

$$\mathbf{p}_{I_0}^{G} = \mathbf{0} \tag{F1a-12}$$

位置设为原点是任意的（全局位置不可观），速度设为零是静止假设的直接推论。

### 2.6 初始协方差

对角协方差矩阵：

$$\mathbf{P}_0 = \text{diag}(\sigma_q^2\mathbf{I}_3,\; \sigma_p^2\mathbf{I}_3,\; \sigma_v^2\mathbf{I}_3,\; \sigma_{bg}^2\mathbf{I}_3,\; \sigma_{ba}^2\mathbf{I}_3) \tag{F1a-13}$$

默认值：
- $\sigma_q = 0.02\,\text{rad}$（姿态不确定度 ~1°）
- $\sigma_p = 0.05\,\text{m}$（位置固定为 0，不确定度小）
- $\sigma_v = 0.01\,\text{m/s}$（静止假设下速度不确定度极小）
- $\sigma_{bg}$、$\sigma_{ba}$：由 bias 估计方法决定（均值估计的不确定度）

---

## 3. 代码实现：逐行对照

### 3.1 窗口划分与统计量

```cpp
// StaticInitializer.cpp:57-64
// 切出两半窗
std::vector<ImuData> window_1to0, window_2to1;
for (const ImuData &data : *imu_data) {
    if (data.timestamp > newesttime - 0.5 * params.init_window_time && data.timestamp <= newesttime) {
        window_1to0.push_back(data);  // 最近半窗
    }
    if (data.timestamp > newesttime - params.init_window_time && data.timestamp <= newesttime - 0.5 * params.init_window_time) {
        window_2to1.push_back(data);  // 前一半窗
    }
}
```

**理论对应**：将初始化窗口 $[t_{\text{oldest}}, t_{\text{newest}}]$ 分为前后两半，分别计算统计量。

### 3.2 静止判定

```cpp
// StaticInitializer.cpp:73-97
// 计算加速度均值和方差
Eigen::Vector3d a_avg_1to0 = Eigen::Vector3d::Zero();
for (const ImuData &data : window_1to0) {
    a_avg_1to0 += data.am;
}
a_avg_1to0 /= (int)window_1to0.size();

double a_var_1to0 = 0;
for (const ImuData &data : window_1to0) {
    a_var_1to0 += (data.am - a_avg_1to0).dot(data.am - a_avg_1to0);
}
a_var_1to0 = std::sqrt(a_var_1to0 / ((int)window_1to0.size() - 1));  // 样本标准差
```

**判定逻辑**（[StaticInitializer.cpp:101-119](src/open_vins/ov_init/src/static/StaticInitializer.cpp#L101-L119)）：

```cpp
// wait_for_jerk = true 时：需要"先静后动"
if (a_var_1to0 < params.init_imu_thresh && wait_for_jerk) {
    return false;  // 还在静止，没检测到运动
}
if (a_var_2to1 > params.init_imu_thresh && wait_for_jerk) {
    return false;  // 之前就在动，不是"静止启动"
}

// wait_for_jerk = false 时：整体静止即可
if ((a_var_1to0 > params.init_imu_thresh || a_var_2to1 > params.init_imu_thresh) && !wait_for_jerk) {
    return false;  // 有运动，不满足静止假设
}
```

### 3.3 Gram-Schmidt 姿态求解

```cpp
// StaticInitializer.cpp:122-125
Eigen::Vector3d z_axis = a_avg_2to1 / a_avg_2to1.norm();  // 重力方向
Eigen::Matrix3d Ro;
InitializerHelper::gram_schmidt(z_axis, Ro);              // Gram-Schmidt 正交化
Eigen::Vector4d q_GtoI = rot_2_quat(Ro);                  // 转四元数
```

**理论对应**：式 (F1a-5)(F1a-6)(F1a-7)。

### 3.4 零偏估计

```cpp
// StaticInitializer.cpp:128-131
Eigen::Vector3d gravity_inG;
gravity_inG << 0.0, 0.0, params.gravity_mag;  // g^G = [0,0,9.81]^T
Eigen::Vector3d bg = w_avg_2to1;              // b_ω = ω̄ (式 F1a-9)
Eigen::Vector3d ba = a_avg_2to1 - quat_2_Rot(q_GtoI) * gravity_inG;  // b_a = ā - R·g (式 F1a-10)
```

### 3.5 IMU 状态设置

```cpp
// StaticInitializer.cpp:135-141
timestamp = window_2to1.at(window_2to1.size() - 1).timestamp;
Eigen::VectorXd imu_state = Eigen::VectorXd::Zero(16);  // 16 维：q(4)+p(3)+v(3)+bg(3)+ba(3)
imu_state.block(0, 0, 4, 1) = q_GtoI;     // 姿态
// 位置 = 0, 速度 = 0（静止假设，式 F1a-11/F1a-12）
imu_state.block(10, 0, 3, 1) = bg;        // 陀螺零偏
imu_state.block(13, 0, 3, 1) = ba;        // 加表零偏
```

### 3.6 初始协方差构造

```cpp
// StaticInitializer.cpp:144-146
// 对角协方差
Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(15, 15);
covariance.block(0, 0, 3, 3) *= std::pow(0.02, 2);   // q: 0.02² rad²
covariance.block(3, 3, 3, 3) *= std::pow(0.05, 2);   // p: 0.05² m²
covariance.block(6, 6, 3, 3) *= std::pow(0.01, 2);   // v: 0.01² (m/s)²
```

**理论对应**：式 (F1a-13)。

---

## 4. 代码-公式变量映射表

| 公式符号 | 含义 | 代码变量 | 代码位置 |
|---------|------|---------|---------|
| $\bar{\mathbf{a}}$ | 加速度均值 | `a_avg_1to0`/`a_avg_2to1` | `:73-97` |
| $\sigma_a$ | 加速度标准差 | `a_var_1to0`/`a_var_2to1` | `:85-97` |
| $\mathbf{z}_{\text{axis}}$ | 重力方向 | `z_axis` | `:122` |
| $\mathbf{R}_{G\to I}$ | 姿态矩阵 | `Ro` | `:123` |
| $\mathbf{b}_\omega$ | 陀螺零偏 | `bg` | `:130` |
| $\mathbf{b}_a$ | 加表零偏 | `ba` | `:131` |
| $\mathbf{P}_0$ | 初始协方差 | `covariance` | `:144-146` |

---

## 5. 关键设计决策

### 5.1 为什么用 Gram-Schmidt 而非其他方法

| 方法 | 优点 | 缺点 |
|------|------|------|
| **Gram-Schmidt** | 简单、快速、数值稳定 | yaw 任意（但 yaw 不可观，无所谓） |
| SVD 分解 | 理论上更优 | 过度设计，没必要 |
| 优化方法 | 可融合更多信息 | 计算量大，初始化不需要 |

### 5.2 为什么忽略 bias 求重力方向

式 (F1a-5) 忽略了 $\mathbf{b}_a$，直接用 $\bar{\mathbf{a}}$ 的方向。原因：
1. 初始化窗口短（~1s），bias 变化可忽略
2. bias 量级小（~0.1）相对于重力（9.81），方向误差 <1°
3. 后续 MLE 会精化 bias 和姿态

### 5.3 为什么速度/位置设为零

- **速度**：静止假设的直接推论（式 F1a-11）
- **位置**：全局位置不可观（VIO 4-DoF 不可观性），设为原点是任意选择

---

## 6. 完整流程图

```
输入: imu_data, wait_for_jerk
        │
        ├─ 1. 窗口划分
        │     ─ window_1to0 (最近半窗), window_2to1 (前一半窗)
        │
        ├─ 2. 统计量计算
        │     ─ a_avg, a_var, w_avg (两窗分别计算)
        │
        ├─ 3. 静止判定
        │     ├─ wait_for_jerk=true: 需"先静后动"
        │     ─ wait_for_jerk=false: 整体静止即可
        │
        ├─ 4. Gram-Schmidt 姿态
        │     ├─ z_axis = ā/‖ā‖
        │     ├─ 选 e_i (与 z 轴夹角最大)
        │     ├─ x_axis = z × e_i / ‖z × e_i
        │     └─ y_axis = z × x
        │
        ├─ 5. 零偏估计
        │     ├─ bg = ω̄
        │     └─ ba = ā - R·g
        │
        ├─ 6. IMU 状态
        │     ├─ q = rot_2_quat(R)
        │     ├─ p = 0, v = 0
        │     └─ bg, ba (上一步结果)
        │
        └─ 7. 初始协方差
              └─ P = diag(0.02², 0.05², 0.01², σ_bg², σ_ba²)
```

---

## 7. 与其他文档的关联

- **上游**：`F1_InertialInitializer-deep-dive.md`（调度层，策略判定）
- **下游**：`StateHelper::set_initial_covariance()`（S12，协方差设置）
- **相关**：`InitializerHelper::gram_schmidt()`（F2d，Gram-Schmidt 实现）
