# 第三方解读 — OpenVINS ZUPT 零速更新（知乎）

> 来源: https://zhuanlan.zhihu.com/p/1959364758944842374 （知乎专栏「OpenVINS星图志」）
> 缓存日期: 2026-08-19 | 语言: 中文 | 质量: ★★★☆

## 零速检测如何工作

零速检测先于 MSCKF 与 SLAM 更新；若检测到零速则仅做零速更新。使用 **IMU 判据** 与 **视觉判据**，可单独或交叉验证（视觉视差判据非必须）。

### IMU 判据
- 取两帧图像间的 IMU 数据，遍历计算残差与雅可比（STEP3）→ 测量压缩（STEP4）→ IMU 判据（STEP5：卡方检验判断 IMU 噪声水平）。
- 状态量更新：姿态、偏置、速度（可选）。
- 配置开关：
  - `integrated_accel_constraint`：false=简单加速度约束（残差=加速度计测量-重力）；true=速度积分约束（残差=预测速度-0）。
  - `model_time_varying_bias`：true 在协方差中考虑偏置随机游走噪声 Q_bias。
  - `override_with_disparity_check`：true 时即使 IMU 静止、视觉视差过大也判为运动。
  - `explicitly_enforce_zero_motion`：false 通过 IMU 残差间接约束；true 直接强制位姿速度不变。
- 残差定义：
  - 残差 A（角速度）：`r_ω = w_m - b_g`，静止时 ω_true=0。
  - 残差 B1（直接加速度）：`r_a = a_hat - R_IG * g`。
  - 残差 B2（加速度积分）：`r_v = v_k - g*dt + R_IG^T * a_hat * dt`。
- 雅可比：A 对 b_g 为 -I；B1 对 R(左扰动) 为 -(R_IG*g)^∧、对 b_a 为 -I；B2 对 R 为 R_IG^T*a_hat^∧*dt、对 b_a 为 -R_IG^T*dt、对 v 为 I。
- 卡方检验（STEP5）：H0 静止，残差 r~N(0,S)，S=H*P_marg*H^T+R；χ²=r^T*S^{-1}*r，α=0.05。超阈值则拒绝 H0（在运动）。`model_time_varying_bias` 为真时 P_marg 偏置块加 Q_bias。

### 视觉判据（STEP6）
- 前后帧特征视差均值：平均视差 < 阈值(1) 且特征对 > 阈值(20) 则通过。
- **已知 BUG**：`disparities.size() < 2` 时仅赋值 `disp_mean=-1` 但未 `return`，后续均值有除零风险。

### 零速判断（STEP7）
```
if (视差太大未通过 && IMU未通过(卡方不通过 || 状态速度>0.1))
    两者均不通过  // 一个判据通过即算通过
```

## 零速更新方程

含完整传播、更新、边缘化。两种模式：
- **模式 A（非强制）**：IMU 测量多少按多少更新，不一定严格静止。
- **模式 B（强制）**：残差强制约束姿态/位置变化量为 0，速度严格为 0。

| 场景 | 模式A | 模式B |
|------|-------|-------|
| 实时性高 | ✅ | ❌ |
| 长时间静止 | ❌ | ✅ |
| IMU 质量好 | ✅ | ❌ |
| 计算有限 | ✅ | ❌ |
| 精度极高 | ❌ | ✅ |
| 频繁切换 | ✅ | ❌ |

## 代码文件
- 相关实现位于 ZUPT 更新器（OpenVINS 通常为 `UpdaterZeroVelocity` 类，`ov_msckf/src/update/`）。
- 代码片段：遍历 `imu_recent` 计算 `res` 与 `H`；`StateHelper::get_marginal_covariance(state, Hx_order)`；卡方检验用 `Eigen::LLT` 与 `boost::math::chi_squared`。

## 注意事项
1. 视觉判据 BUG（除零风险）。
2. 配置参数影响更新约束范围。
3. 卡方阈值：维度 ≥1000 时用 boost 动态分位数并告警。
4. 模式选择权衡；多数情况 A 平衡性好。
5. 视觉判据非必须，但交叉验证更稳。

> 注：本缓存为 AI 摘要，原文含代码截图与公式，精读时建议结合 `UpdaterZeroVelocity` 源码与官方文档 `update-zerovelocity.html`。
