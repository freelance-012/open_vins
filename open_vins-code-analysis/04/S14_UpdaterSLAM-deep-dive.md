# UpdaterSLAM SLAM特征更新 精读报告 (S14)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/update/UpdaterSLAM.cpp`
> **对应论文**: MSCKF 2.0 / MSCKF-SLAM [32] (SLAM 特征持久化进状态 + 锚定表示)、Li & Mourikis 2013 IJRR [27] (时间偏移)、Solà 2017 [47] (Givens QR 零空间投影)、Genova et al. 2020 [21] (锚点切换一致性)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S14（Phase 3 函数级路线，阶段 D 测量更新，SLAM 部分）

> 上游：S9（`EKFUpdate` 预测后协方差）、S12（`initialize` 把 Landmark 加进状态）、S13（MSCKF 路径，本篇对照对象）。本篇聚焦 **`UpdaterSLAM::update` (:253)** —— 与 S13 不同，SLAM 路径把**一部分特征作为持久 `Landmark` 变量留在状态里**，残差同时含"特征位置雅可比 $\mathbf{H}_f$"和"状态雅可比 $\mathbf{H}_x$"，两者**直接拼接进一个大雅可比**（不再零空间投影消去特征，除非是单逆深度表示）。还包括内部子函数 **`change_anchors` (:481)** / **`perform_anchor_change` (:505)** —— 当最旧 clone 被边缘化时，把锚定在它上面的 Landmark **重锚定**到当前帧（保持表示一致性 [21]）。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
VioManager::track_image_and_update / do_feature_propagate_update
        │  (feature_vec: 本帧 SLAM 特征)
        ▼
UpdaterSLAM::update(state, feature_vec)  @ :253  ★ 本篇
   ├─ 1. 清洗特征 (SLAM 特征已进状态, 至少 1 测量)            :270-297
   ├─ 4. 逐特征: get_feature_jacobian_full → 拼 H_x+H_f → χ² → 拼大系统  :321-450
   ├─ 5. EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big)  :470 → S9
   └─ (外部触发) change_anchors(state) @ :481 → perform_anchor_change @ :505

与 S13 的对照:
   S13 MSCKF: 特征不进状态 → 零空间投影消去 H_f → 只更新状态
   S14 SLAM:  特征=Landmark 进状态 → H_f 直接拼入 H_xf → 同时更新状态+特征
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 4 环节（视觉测量更新，SLAM 部分） |
| 输入 | `state`（含 `_features_SLAM`）、`feature_vec`（本帧 SLAM 特征） |
| 输出 | 经 `EKFUpdate` 改写 `state->_Cov` + `_variables`（含 Landmark 均值） |
| 调用方 | `VioManager`（每帧图像触发，与 S13 并行处理两类特征） |
| 被调用方 | `UpdaterHelper::get_feature_jacobian_full`（S13a）、`StateHelper::EKFUpdate`（S9）、`StateHelper::EKFPropagation`（S9，锚点切换用） |

### 1.2 一句话概括

`UpdaterSLAM::update` 与 S13 结构类似（清洗→逐特征 Jacobian→χ²→拼大系统→`EKFUpdate`），**关键差异**是 SLAM 特征已是状态变量：特征位置雅可比 $\mathbf{H}_f$ 直接拼到状态雅可比后面（或单逆深度时先零空间投影 bearing 部分），使卡尔曼更新**同时修正 IMU 状态与 Landmark 位置**。`change_anchors`/`perform_anchor_change` 在滑窗收缩时把 Landmark 的锚点从待删 clone 重锚到当前帧，保持锚定表示的有效性 [21]。

---

## Section 2: 核心数据结构

```cpp
// 文件: ov_msckf/src/update/UpdaterSLAM.cpp:43-59
UpdaterSLAM::UpdaterSLAM(UpdaterOptions &options_slam, UpdaterOptions &options_aruco, FeatureInitializerOptions &feat_init)
    : _options_slam(options_slam), _options_aruco(options_aruco) {
  _options_slam.sigma_pix_sq  = pow(_options_slam.sigma_pix, 2);   // :47 SLAM 特征噪声
  _options_aruco.sigma_pix_sq = pow(_options_aruco.sigma_pix, 2);  // :48 aruco 标签噪声 (独立参数!)
  initializer_feat = make_shared<FeatureInitializer>(feat_init);   // :51 三角化器
  for (i=1;i<500;i++) chi_squared_table[i]=quantile(chi_squared(i),0.95); // :55-58 χ² 表
}
void UpdaterSLAM::update(state, vector<Feature> &feature_vec);     // :253 ★ 本篇
void UpdaterSLAM::change_anchors(state);                           // :481 锚点切换入口
void UpdaterSLAM::perform_anchor_change(state, landmark, new_ts, new_cam_id); // :505 单路标重锚
```

### 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 含义 | 位置 |
|---------|---------|------|------|
| `H_f` | $\mathbf{H}_f$ | 特征位置雅可比 | `:357`→`:296`(S13a) |
| `H_x` | $\mathbf{H}_x$ | 状态雅可比 | `:358` |
| `H_xf` | $[\mathbf{H}_x\,|\,\mathbf{H}_f]$ | 拼接后总雅可比 | `:365` |
| `landmark` | $\mathbf{x}_f$ | 状态中的 Landmark 变量 | `:285`/`:329` |
| `chi2`/`chi2_check` | $\chi^2$ | Mahalanobis 检验 | `:395`/`:400` |
| `_options_aruco`/`_options_slam` | — | aruco/SLAM 独立噪声+χ²阈值 | `:393`/`:409` |

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 特征清洗（与 S13 的差异）

#### 对应代码

```cpp
// 文件: ov_msckf/src/update/UpdaterSLAM.cpp:270-297
auto it0=feature_vec.begin();
while (it0!=feature_vec.end()) {
    (*it0)->clean_old_measurements(clonetimes);                     // :274 只留滑窗内测量
    int ct_meas=0; for (pair:(*it0)->timestamps) ct_meas+=pair.second.size(); // :278-280
    std::shared_ptr<Landmark> landmark = state->_features_SLAM.at((*it0)->featid); // :285 取状态中的路标
    // 单逆深度表示需要 ≥2 测量 (因为要零空间投影 bearing)
    int required_meas = (landmark->_feat_representation==ANCHORED_INVERSE_DEPTH_SINGLE) ? 2 : 1; // :286
    if (ct_meas<1) { (*it0)->to_delete=true; it0=erase(it0); }      // :289-291 不足1→删
    else if (ct_meas<required_meas) { it0=erase(it0); }             // :292-293 不足required→删(不标to_delete)
    else it0++;
}
```

#### 对照注释

| 步骤 | 对应代码 | 行号 |
|------|---------|------|
| 取状态中的 Landmark | `_features_SLAM.at(featid)` | `:285` |
| 单逆深度需≥2测量 | `required_meas = (SINGLE)? 2 : 1` | `:286` |
| 清洗 | `clean_old_measurements` + 计数删 | `:274-293` |

> **与 S13 的关键差异**：① SLAM 特征已是状态变量（`:285` 从 `_features_SLAM` 取），MSCKF 路径（S13）特征不在状态里；② 单逆深度表示因需零空间投影 bearing，要求 ≥2 测量（注释 `:284` 明确说明）。

### 3.2 逐特征：Jacobian 拼接（SLAM 核心差异）

#### 推导说明

MSCKF 路径（S13）把 $\mathbf{H}_f$ **零空间投影消去**（特征不进状态）。SLAM 路径相反：特征 $\mathbf{x}_f$ 已在状态向量中，故把 $\mathbf{H}_f$ **直接拼到** $\mathbf{H}_x$ 右侧得到 $\mathbf{H}_{xf}=[\mathbf{H}_x\,|\,\mathbf{H}_f]$，残差同时约束状态与特征。只有**单逆深度表示**例外——其 bearing（方位）部分已被边缘化，需先对 bearing 做零空间投影，再把剩余的深度雅可比拼入。

#### 对应代码

```cpp
// 文件: ov_msckf/src/update/UpdaterSLAM.cpp:361-383
UpdaterHelper::get_feature_jacobian_full(state, feat, H_f, H_x, res, Hx_order); // :362 → S13a
Eigen::MatrixXd H_xf = H_x;                                                    // :365 初始 = H_x
if (landmark->_feat_representation == ANCHORED_INVERSE_DEPTH_SINGLE) {         // :366 单逆深度特例
    // 把深度雅可比列追加到 H_xf
    H_xf.conservativeResize(H_x.rows(), H_x.cols()+1);                          // :369
    H_xf.block(0,H_x.cols(),H_x.rows(),1) = H_f.block(0,H_f.cols()-1,H_f.rows(),1); // :370 取深度列
    H_f.conservativeResize(H_f.rows(), H_f.cols()-1);                           // :371 去掉深度列
    // 对 bearing 部分零空间投影 (特征已被边缘化, 不取 bearing 为真值 → 保一致性 [21])
    UpdaterHelper::nullspace_project_inplace(H_f, H_xf, res);                   // :376 → S13b
} else {
    // 全特征在状态中 → 直接拼接 H_f 所有列
    H_xf.conservativeResize(H_x.rows(), H_x.cols()+H_f.cols());                 // :381
    H_xf.block(0,H_x.cols(),H_x.rows(),H_f.cols()) = H_f;                       // :382
}
// 雅可比顺序: 状态变量 + landmark 自身
std::vector<Type> Hxf_order = Hx_order; Hxf_order.push_back(landmark);          // :386-387
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{H}_{xf}=[\mathbf{H}_x\,|\,\mathbf{H}_f]$ (全特征) | `H_xf.block(0,H_x.cols())=H_f` | `:382` |
| 单逆深度: 仅深度列拼入 + bearing 零空间投影 | `:369-376` | `:369-376` |

> **要点**：`:386-387` 把 `landmark` 自身加入 `Hxf_order`——这是 SLAM 与 MSCKF 最根本的区别：卡尔曼更新（S9）会**同时修正 Landmark 的均值**（经 S1 `update` boxplus），使路标位置被多帧观测持续精炼。

### 3.3 χ² 检验与全局更新（aruco/SLAM 双参数）

```cpp
// 文件: ov_msckf/src/update/UpdaterSLAM.cpp:389-470
Eigen::MatrixXd P_marg = get_marginal_covariance(state, Hxf_order);            // :390
Eigen::MatrixXd S = H_xf*P_marg*H_xf.transpose();                              // :391
// aruco 标签与 SLAM 特征用**独立**噪声/χ² 阈值
double sigma_pix_sq = ((int)feat.featid < max_aruco_features) ? _options_aruco.sigma_pix_sq : _options_slam.sigma_pix_sq; // :392-393
S.diagonal() += sigma_pix_sq*Ones(S.rows());                                   // :394
double chi2 = res.dot(S.llt().solve(res));                                     // :395
double chi2_check = (res.rows()<500) ? chi_squared_table[res.rows()] : quantile(chi_squared(res.rows()),0.95); // :399-404
double chi2_multipler = ((int)feat.featid < max_aruco_features) ? _options_aruco.chi2_multipler : _options_slam.chi2_multipler; // :408-409
if (chi2 > chi2_multipler*chi2_check) {                                        // :410 异常
    if ((int)feat.featid < max_aruco_features) { /* aruco: 告警并删 */ }        // :411-414
    else { landmark->update_fail_count++; }                                     // :415 SLAM: 失败计数(不直接删)
    (*it2)->to_delete=true; it2=erase(it2); continue;                           // :417-419
}
// 拼入大系统 (同 S13 的 Hx_mapping 去重)
... Hx_big.block(ct_meas, Hx_mapping[var], ...) = H_xf.block(...); ...          // :429-441
R_big.block(ct_meas,ct_meas,res.rows(),res.rows()) *= sigma_pix_sq;            // :444 各向异性 R (aruco/slam 不同)
res_big.block(ct_meas,0,res.rows(),1)=res; ct_meas+=res.rows();                // :447-448
...
StateHelper::EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big);           // :470 → S9 :116
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\chi^2=\mathbf{r}^\top\mathbf{S}^{-1}\mathbf{r}$ | `chi2=res.dot(S.llt().solve(res))` | `:395` |
| aruco/SLAM 双噪声 | `sigma_pix_sq = (featid<max_aruco)? aruco : slam` | `:393` |
| 各向异性 R | `R_big.block(...)*=sigma_pix_sq` | `:444` |
| 卡尔曼更新 | `StateHelper::EKFUpdate(...)` | `:470` |

> **要点**：① aruco 标签与 SLAM 路标用**独立** `sigma_pix` 与 `chi2_multipler`（`:393`/`:409`），因为标签是已知基准、噪声模型不同；② `R_big` 是**各向异性块对角**（每个特征块用各自的 `sigma_pix_sq`，`:444`），而非 S13 的统一各向同性；③ SLAM 特征 χ² 失败时只 `update_fail_count++`（`:415`），不立即删——多次失败后才在 S11 `marginalize_slam` 中边缘化。

### 3.4 锚点切换 `change_anchors` / `perform_anchor_change`（内部子函数）

#### 理论推导

**问题定义**（Geneva MSCKF 2.0）：锚定表示把特征位置表示为"相对某个 anchor clone 的位姿"。设旧锚点 $A_{\text{old}}$ 和新锚点 $A_{\text{new}}$ 在全局系中的位姿分别为 $(\mathbf{R}_{G\to A_{\text{old}}}, \mathbf{p}_{A_{\text{old}}}^G)$ 和 $(\mathbf{R}_{G\to A_{\text{new}}}, \mathbf{p}_{A_{\text{new}}}^G)$。

**刚体变换**：旧→新的旋转和平移为：

$$\mathbf{R}_{A_{\text{old}}\to A_{\text{new}}} = \mathbf{R}_{G\to A_{\text{new}}}\mathbf{R}_{G\to A_{\text{old}}}^\top \tag{S14-1}$$

$$\mathbf{p}_{A_{\text{old}}}^{A_{\text{new}}} = \mathbf{R}_{G\to A_{\text{new}}}(\mathbf{p}_{A_{\text{old}}}^G - \mathbf{p}_{A_{\text{new}}}^G) \tag{S14-2}$$

特征在新锚点系中的坐标：

$$\mathbf{p}_F^{A_{\text{new}}} = \mathbf{R}_{A_{\text{old}}\to A_{\text{new}}}\,\mathbf{p}_F^{A_{\text{old}}} + \mathbf{p}_{A_{\text{old}}}^{A_{\text{new}}} \tag{S14-3}$$

**锚点切换 Jacobian** $\boldsymbol{\Phi}$：对 (S14-3) 做一阶展开，设旧锚点系中特征坐标为 $\mathbf{p}_F^{A_{\text{old}}}$，新锚点系中为 $\mathbf{p}_F^{A_{\text{new}}}$：

$$\delta\mathbf{p}_F^{A_{\text{new}}} \approx \mathbf{R}_{A_{\text{old}}\to A_{\text{new}}}\,\delta\mathbf{p}_F^{A_{\text{old}}} + \mathbf{H}_{x,\text{old}}\,\delta\mathbf{x}_{\text{old}} - \mathbf{H}_{x,\text{new}}\,\delta\mathbf{x}_{\text{new}} \tag{S14-4}$$

其中 $\mathbf{H}_{x,\text{old}}$ 和 $\mathbf{H}_{x,\text{new}}$ 分别是 $\mathbf{p}_F^{A_{\text{new}}}$ 对旧/新锚点位姿误差的 Jacobian（$\mathbf{H}_{x,\text{new}}$ 带负号，因为 $\mathbf{p}_{A_{\text{old}}}^{A_{\text{new}}}$ 依赖于新锚点位姿）。

将 $\mathbf{p}_F^{A_{\text{old}}}$ 用特征表示的 Jacobian 展开：$\delta\mathbf{p}_F^{A_{\text{old}}} = \mathbf{H}_{f,\text{old}}\,\delta\boldsymbol{\lambda}_{\text{old}}$。同时新锚点表示下 $\delta\mathbf{p}_F^{A_{\text{new}}} = \mathbf{H}_{f,\text{new}}\,\delta\boldsymbol{\lambda}_{\text{new}}$。联立解出：

$$\delta\boldsymbol{\lambda}_{\text{new}} = \mathbf{H}_{f,\text{new}}^{-1}\left(\mathbf{R}_{A_{\text{old}}\to A_{\text{new}}}\mathbf{H}_{f,\text{old}}\,\delta\boldsymbol{\lambda}_{\text{old}} + \mathbf{H}_{x,\text{old}}\,\delta\mathbf{x}_{\text{old}} - \mathbf{H}_{x,\text{new}}\,\delta\mathbf{x}_{\text{new}}\right) \tag{S14-5}$$

定义锚点切换 Jacobian $\boldsymbol{\Phi}$：

$$\delta\boldsymbol{\lambda}_{\text{new}} = \boldsymbol{\Phi}\begin{bmatrix} \delta\mathbf{x}_{\text{old}} \\ \delta\boldsymbol{\lambda}_{\text{old}} \\ \delta\mathbf{x}_{\text{new}} \end{bmatrix}, \quad \boldsymbol{\Phi} = \mathbf{H}_{f,\text{new}}^{-1}\begin{bmatrix} \mathbf{H}_{x,\text{old}} & \mathbf{R}_{A_{\text{old}}\to A_{\text{new}}}\mathbf{H}_{f,\text{old}} & -\mathbf{H}_{x,\text{new}} \end{bmatrix} \tag{S14-6}$$

**伪逆处理**：
- `phisize=1`（单逆深度）：$\mathbf{H}_{f,\text{new}}$ 是 $2N\times 1$ 向量，伪逆为 $\mathbf{H}_{f,\text{new}}^+ = \frac{\mathbf{H}_{f,\text{new}}^\top}{\|\mathbf{H}_{f,\text{new}}\|^2}$
- `phisize=3`（3D 位置）：$\mathbf{H}_{f,\text{new}}$ 是 $2N\times 3$ 矩阵，用列主元 QR 求伪逆

**协方差传播**：锚点切换是确定性变换（$\mathbf{Q}=\mathbf{0}$），复用 `EKFPropagation`（S9）：

$$\mathbf{P}_{\lambda\lambda}' = \boldsymbol{\Phi}\,\mathbf{P}_{\text{old}}\,\boldsymbol{\Phi}^\top \tag{S14-7}$$

**FEJ 一致性**：代码（`:558-572`）同步计算 FEJ 版本的变换。因为锚点切换 Jacobian $\boldsymbol{\Phi}$ 是在当前估计值处线性化的，若后续 EKF 更新使用 FEJ 原则，则锚点切换后的特征坐标也必须用 FEJ 值变换，以保持一致性。

#### 对应代码

```cpp
// 文件: ov_msckf/src/update/UpdaterSLAM.cpp:481-647
void UpdaterSLAM::change_anchors(state) {
    if ((int)_clones_IMU.size() <= max_clone_size) return;                     // :484 滑窗未超限→不需切换
    double marg_timestep = state->margtimestep();                              // :491 最旧 clone 时刻
    for (auto &f : state->_features_SLAM) {                                    // :492
      if (f.second 是 GLOBAL 表示) continue;                                   // :494-496 全局表示不锚定
      assert(marg_timestep <= f.second->_anchor_clone_timestamp);              // :498
      if (f.second->_anchor_clone_timestamp == marg_timestep)                  // :499 锚在待删 clone 上
        perform_anchor_change(state, f.second, state->_timestamp, f.second->_anchor_cam_id); // :500 重锚到当前帧
    }
}

void UpdaterSLAM::perform_anchor_change(state, landmark, new_anchor_timestamp, new_cam_id) {
    // 1. 算旧/新 anchor 的相机位姿 (含 FEJ)
    R_GtoOLD/R_GtoNEW = calib_IMUtoCAM.Rot()*clone.Rot();  p_OLDinG/p_NEWinG = ...;  // :538-547 (含 _fej :558-567)
    // 2. 旧→新 anchor 变换
    R_OLDtoNEW = R_GtoNEW*R_GtoOLD.transpose();  p_OLDinNEW = R_GtoNEW*(p_OLDinG-p_NEWinG); // :550-551 (FEJ版 :570-571)
    new_feat.p_FinA = R_OLDtoNEW * landmark->get_xyz(false) + p_OLDinNEW;       // :552 新锚下坐标
    new_feat.p_FinA_fej = R_OLDtoNEW_fej * landmark->get_xyz(true) + p_OLDinNEW_fej; // :572
    // 3. 算旧/新表示的雅可比 (get_feature_jacobian_representation, S13a)
    get_feature_jacobian_representation(state, old_feat, H_f_old, H_x_old, x_order_old); // :525
    get_feature_jacobian_representation(state, new_feat, H_f_new, H_x_new, x_order_new); // :578
    // 4. 锚点切换雅可比 Phi = H_f_new^{-1} * (H_f_old*Δx_f + H_x_old*Δx_old - H_x_new*Δx_new)
    H_f_new_inv = (phisize==1)? (1/||H_f_new||²)*H_f_newᵀ : H_f_new.colPivHouseholderQr().solve(I); // :616-621
    Phi.block(...)= H_f_new_inv*H_x_old[i];   // :625 旧状态块
    Phi.block(...)= H_f_new_inv*H_f_old;      // :629 旧特征块
    Phi.block(...)= -H_f_new_inv*H_x_new[i]; // :633 新锚状态块(负)
    // 5. 协方差传播 (只动 landmark 块)
    StateHelper::EKFPropagation(state, phi_order_NEW={landmark}, phi_order_OLD, Phi, Q=0); // :637 → S9 :36
    // 6. 写回新锚定
    landmark->_anchor_clone_timestamp = new_anchor_timestamp;                   // :643
    landmark->set_from_xyz(new_feat.p_FinA, false);                             // :644
    landmark->set_from_xyz(new_feat.p_FinA_fej, true);                          // :645
    landmark->has_had_anchor_change = true;                                     // :646
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 旧→新 anchor 刚体变换 $\mathbf{R}_{O\to N},\mathbf{p}_{O\in N}$ | `:550-551` (FEJ `:570-571`) | `:550-551` |
| 新锚坐标 $\mathbf{p}_{F\in A_{new}}=\mathbf{R}_{O\to N}\mathbf{p}_{F\in A_{old}}+\mathbf{p}_{O\in N}$ | `:552` | `:552` |
| 锚点切换雅可比 $\Phi=\mathbf{H}_{f,new}^{-1}(\mathbf{H}_{f,old}\,\Delta\mathbf{x}_f+\mathbf{H}_{x,old}\,\Delta\mathbf{x}_{old}-\mathbf{H}_{x,new}\,\Delta\mathbf{x}_{new})$ | `:616-633` | `:616-633` |
| 协方差传播 | `EKFPropagation(phi_order_NEW={landmark}, phi_order_OLD, Phi, Q=0)` | `:637` |

> **设计要点**：`perform_anchor_change` 本质是**状态变量的线性重参数化**（[21]），用 `EKFPropagation`（S9）以 $\Phi$ 传播协方差，保持 Landmark 的均值与协方差在新旧锚定下一致。FEJ 版本（`:558-572`）同步变换，保证一致性。这是 SLAM 路径比纯 MSCKF（S13）多出的维持成本——但换来持久地图点。

---

## Section 4: 函数调用链

```
VioManager (每帧图像, 与 S13 并行)
  └─ UpdaterSLAM::update(state, feature_vec) @ :253  ★ 本篇
        ├─ state->_features_SLAM.at(featid)            (取状态中的 Landmark)  :285/:329
        ├─ UpdaterHelper::get_feature_jacobian_full   @ :362 → S13a
        ├─ UpdaterHelper::nullspace_project_inplace   @ :376 → S13b (仅单逆深度)
        ├─ StateHelper::get_marginal_covariance        (χ² 用, S9 工具)
        └─ StateHelper::EKFUpdate(state, Hx_order_big, Hx_big, res_big, R_big) @ :470 → S9 :116

滑窗收缩触发 (VioManager / S11 marginalize_old_clone 后)
  └─ UpdaterSLAM::change_anchors(state) @ :481  ★ 内部子函数
        └─ perform_anchor_change(state, landmark, ts, cam) @ :505
              ├─ get_feature_jacobian_representation (old/new) @ :525/:578 → S13a
              └─ StateHelper::EKFPropagation(state, {landmark}, OLD, Phi, 0) @ :637 → S9 :36
```

### 关键函数速查

| 函数名 | 文件:行号 | 功能 |
|--------|-----------|------|
| `UpdaterSLAM::update` | `UpdaterSLAM.cpp:253` | SLAM 特征更新主入口 |
| `UpdaterSLAM::change_anchors` | `UpdaterSLAM.cpp:481` | 锚点切换入口（滑窗收缩时） |
| `UpdaterSLAM::perform_anchor_change` | `UpdaterSLAM.cpp:505` | 单路标重锚定 + 协方差传播 |
| `get_feature_jacobian_full` | `UpdaterHelper.cpp:192` | 特征 Jacobian（S13a） |
| `nullspace_project_inplace` | `UpdaterHelper.cpp:426` | 零空间投影（S13b，仅单逆深度） |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `_options_slam.sigma_pix` | yaml | `UpdaterOptions` | 精度 | SLAM 特征像素噪声 | ★★★ 高 |
| `_options_aruco.sigma_pix` | yaml | `UpdaterOptions` | 精度 | aruco 标签像素噪声 | ★★★ 高 |
| `_options_slam.chi2_multipler` | 1.0 | `UpdaterOptions` | 鲁棒性 | SLAM χ² 倍率 | ★★ 中 |
| `_options_aruco.chi2_multipler` | 1.0 | `UpdaterOptions` | 鲁棒性 | aruco χ² 倍率 | ★★ 中 |
| `max_aruco_features` | yaml | `StateOptions` | 调度 | aruco/SLAM 分界 id | ★ 低 |
| `feat_rep_msckf` / landmark 表示 | 配置 | `StateOptions` | 表示 | 特征参数化 | ★★ 中 |

> aruco 与 SLAM 用**独立**噪声/χ²（`:393`/`:409`）——标签是已知基准，噪声模型不同。SLAM 特征 χ² 失败仅计数不删（`:415`），多次失败后才经 S11 `marginalize_slam` 边缘化。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | 锚点切换只在 clone 被边缘化时触发 | ★ 低 | 滑窗边界 | `change_anchors` 仅 `:484` 超限才跑 | 可加周期性重锚以均衡误差 |
| 2 | 单逆深度需≥2测量 | ★ 低 | 短暂单测量 | `:292` 直接删 | 可暂存待下一帧 |
| 3 | 各向异性 R 构造较隐式 | ★ 低 | aruco/slam 混合 | `:444` 块乘 sigma | 可显式构造块对角 R |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `UpdaterSLAM` (S14) | OpenVINS `UpdaterMSCKF` (S13) | VINS-Mono | OKVIS |
|---------|------------------------------|-------------------------------|-----------|-------|
| 特征进状态? | **是**（Landmark 持久） | **否**（零空间投影消去） | 否(图优化) | 否 |
| 雅可比拼接 | $\mathbf{H}_{xf}=[\mathbf{H}_x\,\|\,\mathbf{H}_f]$ 直接拼 | 零空间投影消 $\mathbf{H}_f$ | — | — |
| 锚定表示 | 是（[32]），可重锚（[21]） | 无（特征不存） | 逆深度 | 逆深度 |
| 噪声模型 | aruco/SLAM 双参数、各向异性 R | 统一各向同性 R | 统一 | 统一 |
| 锚点切换 | `perform_anchor_change` (协方差传播) | 无 | 无 | 无 |
| 时间偏移 | 同 S13 (clone 雅可比, S10, [27]) | 同 | 标定 | 标定 |

### 设计取舍分析

- 选择 **SLAM 特征持久化进状态 + 雅可比直接拼接** 是因为：持久地图点提供**跨会话/回环**的全局约束，精度优于纯 MSCKF；代价是协方差维度随路标数增长（O(n²)），故只让"部分"特征走 SLAM 路径（其余走 S13 MSCKF）。
- 选择 **锚定表示 + 重锚定（[21]）** 是因为：相对锚点的参数化比全局 XYZ 数值更稳定（尤其远特征），且滑窗收缩时通过 `perform_anchor_change`（协方差传播，S9）保持一致性，避免引用已删 clone。
- 选择 **aruco/SLAM 双参数 + 各向异性 R** 是因为：标签是已知基准（噪声小/阈值松），普通路标是估计值（噪声大/阈值严），分开调参更鲁棒。

> 与 **S13 MSCKF** 的分工：S13 把特征完全边缘化（算力恒定、精度接近），S14 把部分特征持久化（精度上限更高、算力随路标增长）。两者共用 `get_feature_jacobian_full`（S13a）、`EKFUpdate`（S9）、χ² 框架，仅"特征是否进状态 + 雅可比是否消去"不同。
