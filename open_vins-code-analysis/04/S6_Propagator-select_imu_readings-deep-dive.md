# Propagator 选段 精读报告 (S6)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/Propagator.cpp`
> **对应论文**: 无专门论文（IMU 数据缓冲管理，工程实现）
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S6（Phase 3 函数级路线，阶段 B 选段时间窗）

> 上游：S5（预测编排入口，`:68` 调用本函数）。本篇聚焦 **`Propagator::select_imu_readings` (:269)** —— 从异步到达的 IMU 缓冲中，**截取出落在 $[t_0,t_1]$ 时间窗内的测量段**，并在窗口边界做线性插值，使每段都有非零 `dt`。这是预测准确性的前提：选段错了，后续 F/Qd 累加（S7/S8）全错。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
propagate_and_clone (S5 :68)
        │
        ▼
select_imu_readings(imu_data, time0, time1)  @ :269  ★ 本篇
        │ 返回 prop_data: vector<ImuData>，首尾可能插值生成
        ▼
predict_and_compute (S7) 逐段 (prop_data[i], prop_data[i+1]) 调用
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 预测第 1 步（选段） |
| 输入 | 全量 IMU 缓冲 `imu_data`、`time0`（IMU 起点）、`time1`（IMU 终点）、`warn`（告警开关） |
| 输出 | `prop_data`：落在 $[t_0,t_1]$ 内的 IMU 段（边界插值） |
| 调用方 | `propagate_and_clone` (S5 :68) |
| 被调用方 | `interpolate_data`（线性插值，本文件内） |

### 1.2 一句话概括

`select_imu_readings` 遍历 IMU 缓冲，按"**起点插值 / 中段整段保留 / 终点插值截断**"三种情况（CASE 1/2/3.x）截取 $[t_0,t_1]$ 内的测量，并把末尾不足窗长的段"拉伸"补齐，最后剔除 `dt<1e-12` 的零间隔段，保证输出每段都有有效时间间隔供 S7 积分。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:269-393
std::vector<ov_core::ImuData> Propagator::select_imu_readings(
    const std::vector<ov_core::ImuData> &imu_data,  // 全量缓冲
    double time0, double time1,                     // 目标窗 [t0,t1] (已是 IMU 时间系)
    bool warn = true);                              // 缺失时是否告警
// 内部依赖: Propagator::interpolate_data(a, b, t) —— 在 t 时刻线性插值 a/b 两测量
```

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 三种截取情况

#### 推导说明

IMU 缓冲是离散采样的，而传播窗 $[t_0,t_1]$ 的边界一般**不落在**某个采样点上。需要：
- **起点**：若边界落在两采样之间，在 $t_0$ 处插值生成一段（CASE 1）。
- **中段**：整段都落在窗内的采样，原样保留（CASE 2）。
- **终点**：边界后的第一段，在 $t_1$ 处插值截断（CASE 3.x），并补一段恰在 $t_1$ 的端点。

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:284-342
for (size_t i = 0; i < imu_data.size() - 1; i++) {
  // CASE 1: 起点落在 (imu[i], imu[i+1]) 之间 → 在 time0 插值
  if (imu_data[i+1].timestamp > time0 && imu_data[i].timestamp < time0) {     // :290
    auto data = interpolate_data(imu_data[i], imu_data[i+1], time0);          // :291
    prop_data.push_back(data);                                                // :292
    continue;
  }
  // CASE 2: 整段落在窗内 → 原样保留
  if (imu_data[i].timestamp >= time0 && imu_data[i+1].timestamp <= time1) {   // :301
    prop_data.push_back(imu_data[i]);                                         // :302
    continue;
  }
  // CASE 3: 下一段越过 time1 → 在 time1 截断
  if (imu_data[i+1].timestamp > time1) {                                      // :312
    if (imu_data[i].timestamp > time1 && i == 0) { break; }                   // :317 全在窗后, 无数据
    else if (imu_data[i].timestamp > time1) {                                 // :322
      auto data = interpolate_data(imu_data[i-1], imu_data[i], time1);        // :323
      prop_data.push_back(data);
    } else {
      prop_data.push_back(imu_data[i]);                                       // :328 CASE 3.2
    }
    // 若末尾段终点 != time1, 补一段恰在 time1 的插值
    if (prop_data.back().timestamp != time1) {                                // :334
      auto data = interpolate_data(imu_data[i], imu_data[i+1], time1);        // :335
      prop_data.push_back(data);                                              // :336 CASE 3.3
    }
    break;                                                                    // :340 已到窗末
  }
}
```

#### 对照注释

| 情况 | 触发条件 | 代码行 | 动作 |
|------|---------|-------|------|
| CASE 1 起点插值 | `imu[i+1]>t0 && imu[i]<t0` | `:290-295` | 在 `t0` 插值 |
| CASE 2 中段保留 | `imu[i]>=t0 && imu[i+1]<=t1` | `:301-305` | 原样 push |
| CASE 3.1 全越界 | `imu[i]>t1 && i==0` | `:317-321` | break（无数据） |
| CASE 3.2 末段在窗内 | `imu[i]<=t1` | `:327-331` | push `imu[i]` |
| CASE 3.3 补终点 | 末段终点 `!=t1` | `:334-339` | 在 `t1` 插值 |

### 3.2 末尾拉伸与零 dt 清理

```cpp
// 文件: ov_msckf/src/state/Propagator.cpp:354-389
// 若最后一段仍未到 time1 (缓冲不够) → 拉伸最后一段到 time1
if (prop_data.back().timestamp != time1) {                                    // :358
  auto data = interpolate_data(imu_data[size-2], imu_data[size-1], time1);    // :362
  prop_data.push_back(data);                                                  // :363
}
// 剔除零 dt 段 (否则 Qc=σ²/dt → 无穷大)
for (size_t i = 0; i < prop_data.size()-1; i++) {                            // :371
  if (abs(prop_data[i+1].timestamp - prop_data[i].timestamp) < 1e-12) {       // :372
    prop_data.erase(prop_data.begin()+i);  i--;                              // :376-377
  }
}
// 至少需 2 段才能积分 (dt = 段间差)
if (prop_data.size() < 2) { ... return prop_data; }                          // :382-389
```

#### 对照注释

| 步骤 | 对应代码 | 行号 |
|------|---------|------|
| 末尾拉伸补齐 | `interpolate_data(..., time1)` | `:362-363` |
| 零 dt 剔除 | `erase` if `dt<1e-12` | `:372-377` |
| 最少 2 段校验 | `prop_data.size()<2` | `:382` |

> **设计要点**：
> - **零 dt 清理（`:372`）是数值安全的必须**——S7 的离散噪声 $\mathbf{Q}_c=\sigma^2/\Delta t$（Trawny Eq.129/130，见 S7），$\Delta t\to 0$ 会使噪声协方差趋于无穷。
> - **末尾拉伸（`:358`）** 是当 IMU 缓冲跟不上相机频率时的兜底（代码注释 `:357` 自承"逻辑不精确，应修"），避免传播时间窗不足。
> - 本函数**不修改** `imu_data` 缓冲（只读 + 加锁由 S5 的 `lock_guard` 负责，见 S5 `:67`），返回的是新 vector。

---

## Section 4: 函数调用链

```
Propagator::propagate_and_clone (S5 :68)
  └─ select_imu_readings(imu_data, time0, time1) @ :269  ★ 本篇
        └─ interpolate_data(a, b, t)  (同文件内线性插值)
              ▼ 返回 prop_data
        └─ predict_and_compute (S7 :87) 逐段调用
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `select_imu_readings` | `Propagator.cpp:269` | 选时间窗 IMU 段 | `imu_data,t0,t1` | `vector<ImuData>` | O(n) |
| `interpolate_data` | `Propagator.cpp` (内部) | 线性插值单点 | `a,b,t` | `ImuData` | O(1) |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `time0`/`time1` | 由 S5 算 | `Propagator.cpp:63-64` | 时序 | IMU 系窗边界 | ★★★ 高 |
| `warn` | true | 本函数签名 | 调试 | 缺失告警 | ★ 低 |

> `time0/time1` 由 S5 `:63-64` 经 `calib_dt_CAMtoIMU` 偏移得到；偏移错 → 选错段 → 后续全错。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 末尾拉伸逻辑不精确 | ★ 低 | IMU 缓冲滞后 | 注释 `:357` 自承 TODO | 重写使上界精确等于 time1 |
| 2 | 无数据直接返回空 | ★ 中 | IMU 丢包严重 | S5 后续 `prop_data.size()<2` 跳过积分 | 可加重试/告警升级 |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `select_imu_readings` | VINS-Mono | OKVIS |
|---------|-------------------------------|-----------|-------|
| 选段方式 | 边界线性插值 + 零 dt 清理 | 预积分自带时间对齐 | 按帧选段 |
| 缓冲管理 | 加锁 `imu_data_mtx` | 队列弹出 | 队列弹出 |
| 缺失处理 | 拉伸/跳过 | 等待足够 IMU | 等待足够 IMU |

### 设计取舍分析

- 选择 **边界插值 + 零 dt 剔除** 是因为：EKF 预积分要求每段有精确 `dt`，且 $\mathbf{Q}_c\propto 1/\Delta t$（S7）对零间隔极敏感；插值保证窗边界对齐相机时刻，使克隆（S10）时间戳精确。
- 放弃 **"整段丢弃只留整采样"** 是因为：会引入时间对齐误差，累积到多帧后漂移明显；插值代价低且精度足够。
