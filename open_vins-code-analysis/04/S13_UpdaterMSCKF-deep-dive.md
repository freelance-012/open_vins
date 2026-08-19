# UpdaterMSCKF 多帧投影更新 精读报告 (S13)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/update/UpdaterMSCKF.cpp` + `UpdaterHelper.cpp`
> **对应论文**: Mourikis & Roumeliotis 2007 MSCKF [20] (多帧投影残差 / 零空间投影边缘化特征)、Solà 2017 [47] (Givens/QR 零空间投影)、Li & Mourikis 2013 IJRR [27] (时间偏移)、MSCKF 2.0 [32] (逆深度 / 特征表示)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S13（Phase 3 函数级路线，阶段 D 测量更新，MSCKF 部分）

> 上游：S9（`EKFPropagation` 预测后协方差）、S10（滑窗 clone 已入 `_clones_IMU`）、S11（边缘化删旧 clone）。本篇聚焦 **`UpdaterMSCKF::update` (:58)** —— MSCKF 的"测量更新"核心：对所有可见特征构建**多帧投影残差**，经 **nullspace 投影**（把特征位置变量消掉，只留状态约束）与 **χ² 检验** 后，合并成一个大线性系统交给 `EKFUpdate`（S9）一次性更新。这是 MSCKF "不把特征放进状态" 的关键——特征通过零空间投影被边缘化掉 [20]。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
VioManager::track_image_and_update / do_feature_propagate_update
        │  (feature_vec: 本帧所有特征)
        ▼
UpdaterMSCKF::update(state, feature_vec)  @ :58  ★ 本篇
   ├─ 1. 清洗特征 (删 <2 测量)                       :75-94
   ├─ 2. 构造各相机在 clone 时刻的位姿 clones_cam     :97-115
   ├─ 3. 三角化 + GN 精修每个特征                     :117-143
   ├─ 4. 逐特征: get_feature_jacobian_full → nullspace_project → χ² 检验 → 拼大 Hx/res  :169-256
   ├─ 5. measurement_compress_inplace (测量压缩)      :275
   └─ 6. EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big)  :285 → S9
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 4 环节（视觉测量更新，MSCKF 部分） |
| 输入 | `state`（含滑窗 `_clones_IMU`）、`feature_vec`（本帧特征列表） |
| 输出 | 经 `EKFUpdate` 改写 `state->_Cov` 与 `_variables` 均值；特征标记 `to_delete` |
| 调用方 | `VioManager`（每帧图像触发） |
| 被调用方 | `FeatureInitializer`(三角化)、`UpdaterHelper`(Jacobian/QR/压缩)、`StateHelper::EKFUpdate`(S9) |

### 1.2 一句话概括

`UpdaterMSCKF::update` 对每个特征先三角化其 3D 位置，再算"测量−投影"残差和关于**状态(clone/标定/内参)与特征位置**的联合雅可比；用 **nullspace 投影**（Givens QR，[47]）消去特征位置变量，得到只关于状态的约束；逐特征做 χ² 检验剔除异常，拼成全局大系统后经**测量压缩**交给 `EKFUpdate`（S9）——特征本身始终不进状态（MSCKF 精髓 [20]）。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/update/UpdaterMSCKF.cpp:42-56
UpdaterMSCKF::UpdaterMSCKF(UpdaterOptions &options, FeatureInitializerOptions &feat_init_options)
    : _options(options) {
  _options.sigma_pix_sq = pow(_options.sigma_pix, 2);                    // :45 像素噪声平方
  initializer_feat = make_shared<FeatureInitializer>(feat_init_options); // :48 三角化器
  // 预计算 χ² 表 (0.95 置信, 1..499 维)
  for (i=1;i<500;i++) chi_squared_table[i] = quantile(chi_squared(i),0.95); // :52-55
}
// 主入口
void UpdaterMSCKF::update(state, vector<Feature> &feature_vec);          // :58
```

### 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 含义 | 位置 |
|---------|---------|------|------|
| `res` (per-feat) | $\mathbf{r}_f$ | 单特征投影残差 | `:199`→`:348` |
| `H_f` | $\mathbf{H}_f$ | 关于**特征位置**雅可比 | `:197`→`:296` |
| `H_x` | $\mathbf{H}_x$ | 关于**状态**雅可比 | `:198`→`:297` |
| `Hx_big`/`res_big` | $\mathbf{H}$,$\mathbf{r}$ | 全局大系统 | `:161-162` |
| `R_big` | $\mathbf{R}$ | 各向同性像素噪声 | `:282` |
| `chi2`/`chi2_check` | $\chi^2$ | Mahalanobis 检验 | `:212`/`:217` |

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 特征清洗与滑窗位姿构造

#### 对应代码

```cpp
// 文件: ov_msckf/src/update/UpdaterMSCKF.cpp:58-115
void UpdaterMSCKF::update(state, feature_vec) {
    if (feature_vec.empty()) return;                                     // :61-62
    // 0. 取所有 clone 时间戳 (有效测量时刻)
    vector<double> clonetimes;
    for (clone : _clones_IMU) clonetimes.push_back(clone.first);         // :69-72
    // 1. 清洗: 删测量不足 2 个的特征 (无法三角化/投影)
    auto it0=feature_vec.begin();
    while (it0!=feature_vec.end()) {
      (*it0)->clean_old_measurements(clonetimes);                        // :79 只留滑窗内测量
      int ct_meas=0; for (pair : (*it0)->timestamps) ct_meas+=pair.second.size(); // :83-85
      if (ct_meas<2) { (*it0)->to_delete=true; it0=erase(it0); }         // :88-93 不足2→删
      else it0++;
    }
    // 2. 构造每个相机在 clone 时刻的位姿 (IMU clone → 相机 clone)
    unordered_map<cam_id, map<timestamp, ClonePose>> clones_cam;         // :98
    for (calib : _calib_IMUtoCAM) {                                      // :99
      for (clone : _clones_IMU) {                                        // :103
        R_GtoCi = clone_calib.Rot() * clone_imu.Rot();                   // :106 IMU→CAM 旋转
        p_CioinG = clone_imu.pos() - R_GtoCi.transpose()*clone_calib.pos(); // :107 IMU→CAM 平移
        clones_cami[clone.first] = ClonePose(R_GtoCi, p_CioinG);         // :110
      }
      clones_cam[calib.first] = clones_cami;                             // :114
    }
```

#### 对照注释

| 步骤 | 对应代码 | 行号 |
|------|---------|------|
| 取滑窗时刻 | `clonetimes` from `_clones_IMU` | `:69-72` |
| 测量<2 删除 | `ct_meas<2 → erase` | `:88-93` |
| 相机 clone 位姿 $\mathbf{R}_{G\to C_i},\mathbf{p}_{C_i}$ | `R_GtoCi = R_calib*R_imu; p = p_imu - Rᵀ*p_calib` | `:106-107` |

> **要点**：MSCKF 只用**落在滑窗 `_clones_IMU` 内**的测量（`:79` `clean_old_measurements`）。`:106-107` 把 IMU clone 位姿转成相机位姿（外参 `R_calib_IMUtoCAM`/`p_calib_IMUtoCAM` 来自 S4/S10 的标定变量）。

### 3.2 三角化与 GN 精修

```cpp
// 文件: ov_msckf/src/update/UpdaterMSCKF.cpp:117-143
auto it1=feature_vec.begin();
while (it1!=feature_vec.end()) {
    bool success_tri = true;
    if (initializer_feat->config().triangulate_1d)                      // :123 单逆深度三角化
      success_tri = initializer_feat->single_triangulation_1d(*it1, clones_cam);  // :124
    else success_tri = initializer_feat->single_triangulation(*it1, clones_cam);  // :126
    bool success_refine = true;
    if (initializer_feat->config().refine_features)                     // :131 GN 精修
      success_refine = initializer_feat->single_gaussnewton(*it1, clones_cam);     // :132
    if (!success_tri || !success_refine) { (*it1)->to_delete=true; it1=erase(it1); continue; } // :136-139
    it1++;
}
```

> **要点**：每个特征先三角化出 3D 位置（或逆深度表示），再用高斯牛顿精修。失败（退化/远处）即删除。三角化器 `FeatureInitializer` 用多帧 clone 位姿做几何三角化（见 S4 特征表示 / [32]）。

### 3.3 逐特征：Jacobian → 零空间投影 → χ²

```cpp
// 文件: ov_msckf/src/update/UpdaterMSCKF.cpp:160-256
Eigen::VectorXd res_big=Zero(max_meas_size);                            // :161 全局残差
Eigen::MatrixXd Hx_big=Zero(max_meas_size, max_hx_size);                // :162 全局状态雅可比
unordered_map<Type,size_t> Hx_mapping; vector<Type> Hx_order_big;       // :163-164 变量→列映射

auto it2=feature_vec.begin();
while (it2!=feature_vec.end()) {
    UpdaterHelper::UpdaterHelperFeature feat = {...};                    // :173-194 组装特征(含锚定表示)
    Eigen::MatrixXd H_f, H_x; Eigen::VectorXd res; vector<Type> Hx_order;
    UpdaterHelper::get_feature_jacobian_full(state, feat, H_f, H_x, res, Hx_order); // :203 ★ S14 下沉
    UpdaterHelper::nullspace_project_inplace(H_f, H_x, res);            // :206 ★ 零空间投影
    // χ² 检验 (投影后系统)
    Eigen::MatrixXd P_marg = get_marginal_covariance(state, Hx_order);  // :209
    Eigen::MatrixXd S = H_x*P_marg*H_x.transpose();                      // :210
    S.diagonal() += _options.sigma_pix_sq*Ones(S.rows());                // :211 加 R (各向同性)
    double chi2 = res.dot(S.llt().solve(res));                           // :212
    double chi2_check = (res.rows()<500) ? chi_squared_table[res.rows()] // :216-217 预计算表
                                      : quantile(chi_squared(res.rows()),0.95); // :219-220
    if (chi2 > _options.chi2_multipler * chi2_check) {                  // :225 异常→删
      (*it2)->to_delete=true; it2=erase(it2); continue;
    }
    // 拼入全局大系统 (按 Hx_mapping 去重合并变量)
    size_t ct_hx=0;
    for (var : Hx_order) {                                               // :238
      if (Hx_mapping.find(var)==end) { Hx_mapping[var]=ct_jacob; Hx_order_big.push_back(var); ct_jacob+=var->size(); } // :241-245
      Hx_big.block(ct_meas, Hx_mapping[var], H_x.rows(), var->size()) = H_x.block(0,ct_hx,...,var->size()); // :248
      ct_hx += var->size();
    }
    res_big.block(ct_meas,0,res.rows(),1)=res; ct_meas+=res.rows();      // :253-254
    it2++;
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{r}=\mathbf{z}_{meas}-\pi(\mathbf{p}_f)$ | `res = uv_m - uv_dist` | `:348` (在 `get_feature_jacobian_full`) |
| 零空间投影 $\tilde{\mathbf{H}}_x=\mathbf{N}^\top\mathbf{H}_x$ | `nullspace_project_inplace(H_f,H_x,res)` | `:206` |
| $\chi^2=\mathbf{r}^\top\mathbf{S}^{-1}\mathbf{r}$ | `chi2=res.dot(S.llt().solve(res))` | `:212` |
| 95% 阈值 | `chi_squared_table[res.rows()]` | `:217` |

> **零空间投影（MSCKF 核心 [20]）**：`get_feature_jacobian_full` 返回关于**特征位置**的 $\mathbf{H}_f$（3 列或逆深度维）和关于**状态**的 $\mathbf{H}_x$。`nullspace_project_inplace`（S14 下沉，Givens QR，[47]）把 $\mathbf{H}_f$ 旋成上三角，取左零空间部分——即**消去特征变量、只保留状态约束**。这样特征从不被加入状态，却用多帧约束更新了 clone/标定/内参。

> **χ² 检验**：投影后系统仍可算 Mahalanobis 距离（`:210-212`），超 95% 阈值（可调 `chi2_multipler`）即判该特征异常，删除不更新。预计算 `chi_squared_table`（`:52-55`）避免每次 `quantile` 开销。

### 3.4 测量压缩与全局更新

```cpp
// 文件: ov_msckf/src/update/UpdaterMSCKF.cpp:259-285
for (f=0; f<feature_vec.size(); f++) feature_vec[f]->to_delete=true;    // :261 用完即标记删
if (ct_meas<1) return;                                                  // :266 无有效特征
res_big.conservativeResize(ct_meas,1);                                  // :271 裁剪
Hx_big.conservativeResize(ct_meas, ct_jacob);                           // :272
// 5. 测量压缩 (Givens, 减行不减信息)
UpdaterHelper::measurement_compress_inplace(Hx_big, res_big);           // :275
if (Hx_big.rows()<1) return;                                            // :276
// 各向同性噪声 (压缩后才构造, 因为压缩后维度变了)
Eigen::MatrixXd R_big = _options.sigma_pix_sq * Identity(res_big.rows(), res_big.rows()); // :282
// 6. 一次性更新所有特征
StateHelper::EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big);    // :285 → S9 :116
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 测量压缩 (QR 减行) | `measurement_compress_inplace(Hx_big, res_big)` | `:275` |
| 各向同性 $\mathbf{R}=\sigma_{pix}^2\mathbf{I}$ | `R_big = sigma_pix_sq*I` | `:282` |
| 卡尔曼更新 | `StateHelper::EKFUpdate(...)` | `:285` |

> **测量压缩（`measurement_compress_inplace`，S14 下沉）**：把所有特征拼成的大系统 $\mathbf{H}\in\mathbb{R}^{M\times N}$（$M\gg N$ 时）用 Givens QR 压缩到 $\min(M,N)$ 行——等价于保留全部信息但减少 `EKFUpdate` 的计算量（避免超大 $M\times M$ 的 $\mathbf{S}$ 矩阵）。压缩后才构造各向同性 $\mathbf{R}$（`:282`），因为压缩改变了残差维度。

---

## Section 4: 函数调用链

```
VioManager (每帧图像)
  └─ UpdaterMSCKF::update(state, feature_vec) @ :58  ★ 本篇
        ├─ FeatureInitializer::single_triangulation[_1d]  (三角化)   → S4 特征表示
        ├─ FeatureInitializer::single_gaussnewton         (GN 精修)
        ├─ UpdaterHelper::get_feature_jacobian_full       @ :203 → 下沉 S14
        ├─ UpdaterHelper::nullspace_project_inplace       @ :206 → 下沉 S14 (Givens QR, [47])
        ├─ StateHelper::get_marginal_covariance            (χ² 用, S9 工具)
        ├─ UpdaterHelper::measurement_compress_inplace    @ :275 → 下沉 S14 (Givens QR)
        └─ StateHelper::EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big) @ :285 → S9 :116
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 |
|--------|-----------|------|
| `UpdaterMSCKF::update` | `UpdaterMSCKF.cpp:58` | 多帧投影更新主入口 |
| `get_feature_jacobian_full` | `UpdaterHelper.cpp:192` | 特征 Jacobian（状态+特征） |
| `nullspace_project_inplace` | `UpdaterHelper.cpp:426` | 零空间投影（消去特征） |
| `measurement_compress_inplace` | `UpdaterHelper.cpp:456` | 测量压缩（减行） |
| `StateHelper::EKFUpdate` | `StateHelper.cpp:116` | 卡尔曼更新（S9） |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `sigma_pix` | yaml | `UpdaterOptions` | 精度 | 像素重投影噪声 | ★★★ 高 |
| `chi2_multipler` | 1.0 | `UpdaterOptions` | 鲁棒性 | χ² 检验倍率 | ★★ 中 |
| `triangulate_1d` | 配置 | `FeatureInitOptions` | 精度 | 单逆深度三角化 | ★ 低 |
| `refine_features` | true | `FeatureInitOptions` | 精度 | GN 精修 | ★ 低 |
| `feat_rep_msckf` | 配置 | `StateOptions` | 表示 | 特征参数化(逆深度等) | ★★ 中 |

> `sigma_pix` 直接进 `R_big`（`:282`）——过大→滤波过信任测量；过小→对异常像素敏感。`chi2_multipler` 调大更容忍异常（但可能接纳 outlier）。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 特征不足 2 测量即删 | ★ 低 | 短暂遮挡 | `:88` 直接删 | 可改"暂存"而非删 |
| 2 | 全局大系统 `Hx_big` 可能极大 | ★ 中 | 多相机多特征 | `:162` 预分配 max 尺寸 | 可分批 `EKFUpdate` |
| 3 | `R` 各向同性假设 | ★ 低 | 非均匀像素噪声 | `:282` 均匀对角 | 支持各向异性 R |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `UpdaterMSCKF` | VINS-Mono | OKVIS | 常规 MSCKF [20] |
|---------|------------------------|-----------|-------|-----------------|
| 特征进状态? | **否**（零空间投影消去） | 否(图优化) | 否 | 否 |
| 多帧约束 | 全部滑窗 clone | 滑窗关键帧 | 滑窗 | 滑窗 |
| 零空间投影 | Givens QR ([47]) | 无 | 部分 | 解析 |
| 测量压缩 | `measurement_compress` | 无 | Schur | 无 |
| 特征表示 | 逆深度/锚定 ([32]) | 逆深度 | 逆深度 | XYZ |
| 时间偏移 | clone 雅可比注入 (S10, [27]) | 标定 | 标定 | 标定 |

### 设计取舍分析

- 选择 **零空间投影消去特征 + 不进状态** 是因为：MSCKF [20] 的核心优势——特征不被边缘化进协方差，避免 O(n²) 膨胀；多帧约束通过投影保留，精度接近 SLAM 但算力恒定。
- 选择 **Givens QR 做投影+压缩** 是因为：数值稳定（[47] 算法 5.2.4），且投影与压缩可复用同一套旋转逻辑（`nullspace_project_inplace` 与 `measurement_compress_inplace` 同源）。
- 选择 **三角化+GN 精修前置** 是因为：投影残差需要精确特征初值；精修后用其雅可比做线性更新，比纯线性化更准（见 S4 特征表示 / [32]）。

> 与 **S14 UpdaterSLAM** 的对比：MSCKF 路径（本篇）把特征完全边缘化掉；SLAM 路径（S14）则把**部分**特征作为持久 Landmark 加入状态（经 S12 `initialize`），用 `feat_rep_msckf`/`ANCHORED_*` 表示区分两类特征的处理。
