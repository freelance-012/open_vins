# UpdaterSLAM delayed_init 延迟初始化 精读报告 (S15)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/update/UpdaterSLAM.cpp`
> **对应论文**: MSCKF-SLAM [32] (SLAM 特征持久化)、Genova et al. 2020 [21] (锚定表示)、Trawny 2005 [40] (JPL 姿态)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S15（Phase 3 函数级路线，阶段 D 测量更新，SLAM 延迟初始化）

> 上游：S12（`initialize` 把 Landmark 加进状态并做 χ² 门限）、S13（MSCKF 路径，特征不进状态）、S14（`UpdaterSLAM::update` 实时更新已进状态的 SLAM 特征）。本篇聚焦 **`UpdaterSLAM::delayed_init` (:61)** —— 与 S14 不同，它是一个**批量延迟初始化**入口：当 MSCKF 滑动窗口里某个特征已经积累了足够多帧观测、三角化收敛后，把它**升格 (promote)** 为一个持久的 `Landmark` 变量插进 `state->_features_SLAM`。这一步是"MSCKF 短期特征 → SLAM 长期路标"的桥接点。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
VioManager::do_feature_propagate_update
        │  (feats_slam / 已 lost 但满足延迟初始化条件的特征)
        ▼
UpdaterSLAM::delayed_init(state, feature_vec)  @ :61  ★ 本篇
   ├─ 0. 收集所有 clone 时间戳 (滑窗有效观测时刻)               :71-75
   ├─ 1. 清洗特征测量，ct_meas<2 的直接删除                    :78-97
   ├─ 2. 构造每相机、每 clone 时刻的相机位姿 clones_cam          :100-118
   ├─ 3. 逐特征三角化 (+GN 精修)，失败删除                     :120-145
   ├─ 4. 逐特征: get_feature_jacobian_full → (单逆深度时零空间投影 bearing) → 建 Landmark → StateHelper::initialize → 成功则插 _features_SLAM  :148-241
   └─ 内部复用: get_feature_jacobian_full (S13a)、nullspace_project_inplace (S13b)、StateHelper::initialize (S12)

与 S13 / S14 的对照:
   S13 MSCKF:          特征从不进状态 → 零空间投影消去 H_f → 只更新状态 (短期)
   S14 UpdaterSLAM::update: 特征=Landmark 已进状态 → H_f 直接拼入 → 实时联合更新 (长期)
   S15 delayed_init:   特征从"未进状态"→ 升格为 Landmark 进状态 (MSCKF→SLAM 的晋升门槛)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | 第 4 环节（视觉测量更新，SLAM 延迟初始化段） |
| 输入 | `state`（含 `_clones_IMU` / `_calib_IMUtoCAM`）、`feature_vec`（待升格的候选特征，通常来自 `feats_lost`/`feats_slam`） |
| 输出 | 成功者写入 `state->_features_SLAM`，并从 `feature_vec` 移除（`to_delete=true`）；失败者从 `feature_vec` 删除 |
| 调用方 | `VioManager`（在特征分组 `feats_lost`/`feats_slam` 处理阶段调用） |
| 被调用方 | `FeatureInitializer::single_triangulation`/`_1d`/`single_gaussnewton`、`UpdaterHelper::get_feature_jacobian_full`（S13a）、`UpdaterHelper::nullspace_project_inplace`（S13b）、`StateHelper::initialize`（S12） |

---

## Section 2: 函数骨架与调用关系

`delayed_init` 是一个**四阶段流水线**，三个阶段带计时 (`boost::posix_time`) 用于调试：

```
delayed_init(state, feature_vec)
│
├─ 阶段 0  : clonetimes = keys of state->_clones_IMU          :71
│
├─ 阶段 1  : 对每个 feature:
│            - clean_old_measurements(clonetimes)             :82
│            - 统计 ct_meas (各相机各时间戳的 uvs 数)           :85-88
│            - ct_meas < 2 → to_delete + erase                :91-93
│
├─ 阶段 2  : 构造 clones_cam[cam_id][clone_timestamp]         :100-118
│            R_GtoCi = calib.Rot() * clone.Rot()
│            p_CioinG = clone.pos() - R_GtoCi^T * calib.pos()
│
├─ 阶段 3  : 对每个 feature:
│            - single_triangulation(_1d)  (可配)              :126-130
│            - single_gaussnewton 精修 (可配 refine_features)  :134-135
│            - 任一失败 → to_delete + erase                   :139-142
│
└─ 阶段 4  : 对每个 feature (核心初始化):                     :148
             - 转 UpdaterHelperFeature + 选 feat_rep          :152-165
             - get_feature_jacobian_full → H_f, H_x, res      :185  (S13a)
             - 若 ANCHORED_INVERSE_DEPTH_SINGLE:
                 追加 depth 列到 H_xf + 零空间投影 bearing      :190-206  (S13b)
             - 建 Landmark (size=1 或 3) + 设锚点/表示/xyz     :210-223
             - R = sigma_pix_sq * I (aruco/slam 双 sigma)      :226-228
             - StateHelper::initialize(...) → 成功插 _features_SLAM :233-234
```

**关键点**：`delayed_init` 与 `UpdaterSLAM::update`（S14）**共用** `get_feature_jacobian_full` 与 `nullspace_project_inplace`（S13a/S13b），但目的不同——S14 是"已经进状态的 Landmark 做实时残差更新"，S15 是"还没进状态的特征，先建好 Landmark 并一次性塞进状态"。

---

## Section 3: 逐段精读

### 3.1 阶段 0–1：滑窗时刻收集 + 测量清洗 (`:71-97`)

```cpp
std::vector<double> clonetimes;
for (const auto &clone_imu : state->_clones_IMU)
    clonetimes.emplace_back(clone_imu.first);          // :71-75
```
- 滑窗里所有 IMU clone 的时间戳，就是"有效观测时刻"。后面 `clean_old_measurements` 会丢掉不落在这些时刻上的测量。

```cpp
(*it0)->clean_old_measurements(clonetimes);            // :82
int ct_meas = 0;
for (const auto &pair : (*it0)->timestamps)
    ct_meas += (*it0)->timestamps[pair.first].size();  // :85-88
if (ct_meas < 2) {                                     // :91
    (*it0)->to_delete = true;
    it0 = feature_vec.erase(it0);                      // 删除：观测不足，三角化不可靠
}
```
- `ct_meas` 统计的是**所有相机 × 所有时间戳**下的归一化坐标 `uvs_norm` 数量。少于 2 个观测点时无法三角化，直接剔除。

### 3.2 阶段 2：相机位姿表 clones_cam (`:100-118`)

```cpp
for (const auto &clone_calib : state->_calib_IMUtoCAM) {
  for (const auto &clone_imu : state->_clones_IMU) {
    Eigen::Matrix3d R_GtoCi = clone_calib.second->Rot() * clone_imu.second->Rot();
    Eigen::Vector3d p_CioinG = clone_imu.second->pos()
                              - R_GtoCi.transpose() * clone_calib.second->pos();
    clones_cami.insert({clone_imu.first, FeatureInitializer::ClonePose(R_GtoCi, p_CioinG)});
  }
  clones_cam.insert({clone_calib.first, clones_cami});
}
```
- 对每个**相机外参** `IMUtoCAM` 和每个 **IMU clone**，预计算该时刻的相机在世界系的位姿 `(R_GtoCi, p_CioinG)`。
- 三角化时，特征在多个相机时刻的观测都投影到这个统一的 `clones_cam` 表上做线性/非线性求解（见 S13a / FeatureInitializer）。

### 3.3 阶段 3：三角化 + GN 精修 (`:120-145`)

```cpp
bool success_tri = true;
if (initializer_feat->config().triangulate_1d)        // 1D 逆深度三角化（FAST 模式）
    success_tri = initializer_feat->single_triangulation_1d(*it1, clones_cam);
else
    success_tri = initializer_feat->single_triangulation(*it1, clones_cam);

bool success_refine = true;
if (initializer_feat->config().refine_features)       // GN 精修（可选）
    success_refine = initializer_feat->single_gaussnewton(*it1, clones_cam);

if (!success_tri || !success_refine) {                // 三角化失败 → 丢弃
    (*it1)->to_delete = true;
    it1 = feature_vec.erase(it1);
    continue;
}
```
- **三角化**把多帧 2D 观测反投影成一个 3D 点（相对锚点的逆深度表示，见 [32][21]）。
- **GN 精修**（`single_gaussnewton`）对该 3D 点再做若干次高斯-牛顿迭代，最小化重投影误差，提升初始化精度。
- 任一失败即丢弃——延迟初始化宁缺毋滥，避免坏路标污染状态。

### 3.4 阶段 4：特征雅可比 + 单逆深度零空间投影 (`:148-206`)

#### 3.4.1 表示选择 (`:159-165`)

```cpp
auto feat_rep = ((int)feat.featid < state->_options.max_aruco_features)
                    ? state->_options.feat_rep_aruco   // Aruco 用 aruco 表示
                    : state->_options.feat_rep_slam;   // 普通 SLAM 用 slam 表示
if (feat_rep == ANCHORED_INVERSE_DEPTH_SINGLE)
    feat.feat_representation = ANCHORED_MSCKF_INVERSE_DEPTH;  // 单逆深度等价于 MSCKF 逆深度
```
- 用 `featid` 是否小于 `max_aruco_features` 区分 Aruco 板角点与普通特征，分别取不同表示。
- 若配置为 `ANCHORED_INVERSE_DEPTH_SINGLE`（单逆深度 = 用 bearing + 1 个逆深度），则**临时降级为** `ANCHORED_MSCKF_INVERSE_DEPTH`（bearing + 逆深度，2 自由度），因为延迟初始化要把它直接建为深度的 Landmark（见 3.4.3）。

#### 3.4.2 取特征雅可比 (`:178-185`)

```cpp
Eigen::MatrixXd H_f, H_x; Eigen::VectorXd res;
std::vector<std::shared_ptr<Type>> Hx_order;
UpdaterHelper::get_feature_jacobian_full(state, feat, H_f, H_x, res, Hx_order);  // :185 → S13a
```
- `H_f`：残差对**特征参数**（bearing ± 逆深度）的雅可比。
- `H_x`：残差对**状态变量**（pose/bias/landmark 等）的雅可比。
- `Hx_order`：状态变量的顺序（对应 `H_x` 的列）。

#### 3.4.3 单逆深度的 bearing 零空间投影 (`:190-206`)

**数学模型**：对于 `ANCHORED_INVERSE_DEPTH_SINGLE` 表示，特征参数为 $\boldsymbol{\lambda} = [\mathbf{b}^\top, \rho]^\top$（bearing $\mathbf{b} = [\alpha, \beta]^\top$ + 逆深度 $\rho$）。系统为：

$$\mathbf{r} = \mathbf{H}_x\,\delta\mathbf{x} + \mathbf{H}_f\begin{bmatrix} \delta\mathbf{b} \\ \delta\rho \end{bmatrix} + \mathbf{n} \tag{S15-1}$$

其中 $\mathbf{H}_f = [\mathbf{H}_{f,b} \;\; \mathbf{h}_{f,\rho}]$，$\mathbf{H}_{f,b} \in \mathbb{R}^{2N\times 2}$（bearing 的 2 列），$\mathbf{h}_{f,\rho} \in \mathbb{R}^{2N\times 1}$（depth 的 1 列）。

**步骤**（对应代码 `:190-206`）：

**(a) 将 depth 列并入状态侧**：构造 $\mathbf{H}_{xf} = [\mathbf{H}_x \;\; \mathbf{h}_{f,\rho}]$，令 $\tilde{\mathbf{H}}_f = \mathbf{H}_{f,b}$（仅 bearing）。方程重写为：

$$\mathbf{r} = \mathbf{H}_{xf}\begin{bmatrix} \delta\mathbf{x} \\ \delta\rho \end{bmatrix} + \tilde{\mathbf{H}}_f\,\delta\mathbf{b} + \mathbf{n} \tag{S15-2}$$

**(b) 对 bearing 做零空间投影**（Givens QR，S13）：对 $\tilde{\mathbf{H}}_f$ 做 QR 分解 $\mathbf{Q}^\top\tilde{\mathbf{H}}_f = [\mathbf{T}_b^\top, \mathbf{0}^\top]^\top$，左乘 $\mathbf{Q}^\top$：

$$\begin{bmatrix} \mathbf{Q}_1^\top\mathbf{r} \\ \mathbf{Q}_2^\top\mathbf{r} \end{bmatrix} = \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{H}_{xf} \\ \mathbf{Q}_2^\top\mathbf{H}_{xf} \end{bmatrix}\begin{bmatrix} \delta\mathbf{x} \\ \delta\rho \end{bmatrix} + \begin{bmatrix} \mathbf{T}_b \\ \mathbf{0} \end{bmatrix}\delta\mathbf{b} + \begin{bmatrix} \mathbf{Q}_1^\top\mathbf{n} \\ \mathbf{Q}_2^\top\mathbf{n} \end{bmatrix} \tag{S15-3}$$

下半部分（$2N-2$ 个方程）中 bearing 误差 $\delta\mathbf{b}$ 被完全消去：

$$\mathbf{Q}_2^\top\mathbf{r} = \mathbf{Q}_2^\top\mathbf{H}_{xf}\begin{bmatrix} \delta\mathbf{x} \\ \delta\rho \end{bmatrix} + \mathbf{Q}_2^\top\mathbf{n} \tag{S15-4}$$

**(c) 拆分回 $H_x$ 和 $H_f$**：投影后的 $\mathbf{Q}_2^\top\mathbf{H}_{xf}$ 拆回状态部分 $\tilde{\mathbf{H}}_x$ 和深度部分 $\tilde{\mathbf{h}}_\rho$。

**为什么需要投影**：三角化得到的 bearing $(\alpha_0, \beta_0)$ 有误差，若直接作为真值初始化，Jacobian 中不包含 bearing 的误差项，滤波器会低估不确定性（过于乐观）。零空间投影将 bearing 的误差"吸收"到投影算子 $\mathbf{Q}_2^\top$ 中，投影后的 Jacobian 隐式地包含了 bearing 不确定性对深度和状态的影响。数学上等价于将 bearing 视为随机变量并边缘化：

$$p(\delta\mathbf{x}, \delta\rho \mid \mathbf{r}) = \int p(\delta\mathbf{x}, \delta\rho, \delta\mathbf{b} \mid \mathbf{r})\,d\delta\mathbf{b} \tag{S15-5}$$

投影后的有效方程数为 $2N-2$（减去 2 个 bearing 维度）。最终 `landmark_size = 1`（只有 depth 进状态）。

### 3.5 阶段 4（续）：建 Landmark + 调 initialize (`:208-241`)

```cpp
int landmark_size = (feat_rep == ANCHORED_INVERSE_DEPTH_SINGLE) ? 1 : 3;  // :210
auto landmark = std::make_shared<Landmark>(landmark_size);
landmark->_featid = feat.featid;
landmark->_feat_representation = feat_rep;
landmark->_unique_camera_id = (*it2)->anchor_cam_id;
if (LandmarkRepresentation::is_relative_representation(feat.feat_representation)) {
    landmark->_anchor_cam_id = feat.anchor_cam_id;
    landmark->_anchor_clone_timestamp = feat.anchor_clone_timestamp;
    landmark->set_from_xyz(feat.p_FinA, false);    // 值
    landmark->set_from_xyz(feat.p_FinA_fej, true); // FEJ 值
} else {
    landmark->set_from_xyz(feat.p_FinG, false);
    landmark->set_from_xyz(feat.p_FinG_fej, true);
}
```
- `landmark_size`：单逆深度 → 1（只有深度），否则 → 3（xyz 全局或相对坐标）。
- 锚点信息（`anchor_cam_id` / `anchor_clone_timestamp`）和 3D 点（值 + FEJ 值）都从三角化结果填上。

```cpp
double sigma_pix_sq =
    ((int)feat.featid < state->_options.max_aruco_features)
        ? _options_aruco.sigma_pix_sq : _options_slam.sigma_pix_sq;   // :226-227
Eigen::MatrixXd R = sigma_pix_sq * Eigen::MatrixXd::Identity(res.rows(), res.rows());  // :228

double chi2_multipler =
    ((int)feat.featid < state->_options.max_aruco_features)
        ? _options_aruco.chi2_multipler : _options_slam.chi2_multipler;  // :231-232

if (StateHelper::initialize(state, landmark, Hx_order, H_x, H_f, R, res, chi2_multipler)) {  // :233 → S12
    state->_features_SLAM.insert({(*it2)->featid, landmark});  // 成功：升格进 SLAM 状态
    (*it2)->to_delete = true;
    it2++;
} else {
    (*it2)->to_delete = true;          // 失败：χ² 门限不过，丢弃
    it2 = feature_vec.erase(it2);
}
```
- **测量噪声 `R`** 用 `sigma_pix_sq * I`，Aruco 与普通特征各用各的 sigma（双参数，与 S14 一致）。
- **`StateHelper::initialize`（S12）**：用 `H_x`（状态雅可比）、`H_f`（特征雅可比）、`R`、`res` 做"新变量加入状态"的 EKF 更新（含 Givens QR 求逆 + χ² 门限）。返回 `true` 表示该特征初始化成功、可信任。
- 成功 → 以 `featid` 为键插入 `state->_features_SLAM`（之后就归 S14 的 `update` 实时维护）；失败 → 删除候选。

---

## Section 4: 接口与数据结构

### 4.1 签名

```cpp
void UpdaterSLAM::delayed_init(std::shared_ptr<State> state,
                               std::vector<std::shared_ptr<Feature>> &feature_vec);  // :61
```

| 参数 | 类型 | 方向 | 说明 |
|------|------|------|------|
| `state` | `std::shared_ptr<State>` | in/out | 读取 `_clones_IMU`/`_calib_IMUtoCAM`/`_options`；写出 `_features_SLAM` 与协方差 |
| `feature_vec` | `std::vector<std::shared_ptr<Feature>> &` | in/out | 候选特征；成功/失败者均标记 `to_delete`，部分被 `erase` |

### 4.2 成员依赖

| 成员 | 作用 |
|------|------|
| `initializer_feat` | `FeatureInitializer` 实例，负责三角化与 GN 精修（阶段 3） |
| `_options_aruco` / `_options_slam` | Aruco / SLAM 各自的 `sigma_pix_sq`、`chi2_multipler` |

---

## Section 5: 数学与公式对照

### 5.1 三角化（阶段 3）

多帧观测 $\mathbf{z}_{c,i} = \pi(\mathbf{R}_{Gc_i}, \mathbf{p}_{C_i}, \mathbf{p}_F)$（相机投影），反解逆深度参数 $\boldsymbol{\rho}$。线性 `single_triangulation` 用 DLT / 逆深度最小二乘；`single_gaussnewton` 在此基础上迭代最小化
$$
\min_{\boldsymbol{\rho}} \sum_i \| \mathbf{z}_{c,i} - \pi(\mathbf{T}_{Gc_i}, \boldsymbol{\rho}) \|^2_{\boldsymbol{\Sigma}_{pix}}
$$
（详见 S13a / FeatureInitializer 文档；本篇不重复展开，列入"待详细补充项"）。

### 5.2 单逆深度的零空间投影（阶段 4.3, :190-206）

设特征参数 $\boldsymbol{\lambda} = [\mathbf{b}^\top, \rho]^\top$（bearing $\mathbf{b}$ + 逆深度 $\rho$）。系统为
$$
\mathbf{r} = \mathbf{H}_f \begin{bmatrix} \mathbf{b} \\ \rho \end{bmatrix} + \mathbf{H}_x \delta \mathbf{x}
$$
把 $\rho$ 列并入状态侧得 $\mathbf{H}_{xf} = [\mathbf{H}_x \mid \mathbf{H}_{f,\rho}]$，并令 $\tilde{\mathbf{H}}_f = \mathbf{H}_{f,\mathbf{b}}$（仅 bearing）。执行（S13b）
$$
\mathbf{N} = \mathrm{null}(\tilde{\mathbf{H}}_f^\top), \quad
\begin{bmatrix} \mathbf{H}'_{xf} \\ \mathbf{r}' \end{bmatrix} = \mathbf{N}^\top \begin{bmatrix} \mathbf{H}_{xf} \\ \mathbf{r} \end{bmatrix}
$$
投影后 $\tilde{\mathbf{H}}_f$ 被消去，特征只剩 $\rho$ 一维，且与状态联合满足零空间约束——**保证初始化时不把 bearing 当作已知量，维持 FEJ 一致性**（对应 [21][32] 的锚定逆深度一致性要求）。

### 5.3 新变量加入状态（阶段 4.5, `StateHelper::initialize` → S12）

对"状态 + 新 Landmark"联合系统做标准 EKF 更新：
$$
\mathbf{H}_{xf} = \begin{bmatrix} \mathbf{H}_x & \mathbf{H}_f \end{bmatrix},\quad
\mathbf{S} = \mathbf{H}_{xf}\mathbf{P}\mathbf{H}_{xf}^\top + \mathbf{R},\quad
\chi^2 = \mathbf{r}^\top \mathbf{S}^{-1} \mathbf{r}
$$
$\chi^2$ 超过 `chi2_multipler * χ²_table[自由度]` 则拒绝。通过后把 Landmark 的均值/协方差拼入 `state`。

---

## Section 6: 关键实现细节与易错点

1. **`ct_meas < 2` 的硬门槛**（`:91`）：少于 2 个观测点无法三角化，直接丢。注意 `ct_meas` 统计的是全部相机×时间戳的 `uvs_norm` 总数，不是帧数。
2. **单逆深度的"先降后升"**（`:163-165` + `:210`）：配置 `ANCHORED_INVERSE_DEPTH_SINGLE` 时，先把特征表示降级为 MSCKF 逆深度（2 自由度 bearing+depth），再在 4.3 用零空间投影去掉 bearing，最终 Landmark 是**单深度 1 维**。这是延迟初始化与 S14 实时更新在单逆深度上的关键差异。
3. **Aruco / SLAM 双参数**（`:226-232`）：`sigma_pix_sq` 和 `chi2_multipler` 按 `featid < max_aruco_features` 分流，与 S14 的 `update` 一致。
4. **`to_delete` 的两种后果**：成功 → 标记 `to_delete=true` 但保留在 vector（由上层统一回收）；失败 → 直接 `erase`。两者最终都会让该特征退出 MSCKF 候选池。
5. **`get_feature_jacobian_full` / `nullspace_project_inplace` 复用**：与 S13 共用（S13a/S13b），但此处是为了"建新变量"而非"消去旧变量"。

---

## Section 7: 与其他文档的关联 & 待详细补充项

- **上游**：S12（`StateHelper::initialize` 的 Givens QR + χ² 细节）、S13（`UpdaterMSCKF::update` 的零空间投影框架）、S14（`UpdaterSLAM::update` 实时维护已进状态的 Landmark）。
- **下游调用方**：S17（`VioManager::do_feature_propagate_update` 在 `feats_lost`/`feats_slam` 分组时触发本函数）。
- **对照**：S15 是"MSCKF 特征 → SLAM Landmark"的**晋升门槛**；S14 是"已晋升 Landmark 的实时更新"；S13 是"始终不晋升的纯 MSCKF 特征更新"。

### 待详细补充项（核心复杂内部模块，后续需要时拆子文档）

- **S15a — `FeatureInitializer::single_triangulation` / `single_triangulation_1d`**：DLT / 逆深度线性三角化的具体构造与数值稳定性处理（`UpdaterSLAM.cpp:126-130` 调用，实现在 `ov_core/src/FeatureInitializer.*`）。
- **S15b — `FeatureInitializer::single_gaussnewton`**：GN 精修的迭代步长、收敛判据、雅可比对锚点相机位姿的求导（`UpdaterSLAM.cpp:135` 调用）。
- **S15c — `get_feature_jacobian_full` 相对/全局表示分支**：本篇只用到其输出 `H_f/H_x/res`，内部如何按 `ANCHORED_*` 表示分支构造雅可比，见 S13a。
- **S15d — `nullspace_project_inplace` 的 Givens QR 实现**：见 S13b。

> 以上子文档暂不展开，待需要逐项深挖时再补（遵循 SKILL.md "04/ 文档拆分规则"，主文档 + Sxxa/b/c 子文档 + 懒补充清单）。
