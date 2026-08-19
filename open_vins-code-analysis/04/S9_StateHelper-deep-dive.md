# StateHelper 协方差代数 精读报告 (S9)

> **仓库**: open_vins
> **模块路径**: `ov_msckf/src/state/StateHelper.{h,cpp}`
> **对应论文**: Mourikis & Roumeliotis 2007 MSCKF [20] (滑窗边缘化/克隆)、Li & Mourikis 2013 IJRR [27] (时间偏移雅可比)、Trawny 2005 [40] (误差态 EKF 协方差更新)、Solà 2017 [47] (Givens QR 初始化)
> **生成日期**: 2026-08-19
> **分析者**: Gavin + AI
> **精读锚点**: §7 S9（Phase 3 函数级路线，阶段 C 状态增广与边缘化）

> 上游：S5（预测+克隆入口）。本篇深入**协方差矩阵本身的所有线性代数操作**——`StateHelper` 是把理论公式翻译成 `Eigen` 块操作的无状态工具集。它不持有状态，只对 `State` 内的 `_Cov` 和 `_variables` 列表做"预测前推 / 卡尔曼更新 / 增广 / 边缘化"。理解本文件，等于掌握了 MSCKF 滑窗滤波器每一帧协方差变形的底层机制。

---

## Section 1: 模块定位与职责

### 1.1 在数据流中的位置

```
S5 Propagator::propagate_and_clone
        │
        ├─► StateHelper::EKFPropagation(...)      【预测】ΦPΦᵀ+Qd 写回 _Cov
        ├─► StateHelper::augment_clone(state,w)   【增广】克隆 _imu->pose() 入滑窗
        │         └─ clone(state, pose)
        │
VioManager::feed_measurement_monocular/stereo
        │
        └─► UpdaterMSCKF::update(...) / UpdaterSLAM::update(...)
                  ├─ StateHelper::EKFUpdate(...)          【更新】K·res 写回 _Cov + _variables
                  ├─ StateHelper::initialize(...)         【新变量】SLAM 路标首次可见
                  │         ├─ initialize_invertible(...)
                  │         └─ EKFUpdate(...) (零空间投影部分)
                  └─ StateHelper::marginalize(...)        【边缘化】删旧 clone / SLAM 路标
                            ├─ marginalize_old_clone(...)
                            └─ marginalize_slam(...)
```

| 属性 | 内容 |
|------|------|
| 数据流位置 | **第 3 环节（协方差代数中枢）**：被预测器、更新器、初始化器共同调用的"纯函数式"工具集 |
| 输入 | `std::shared_ptr<State>`（持有 `_Cov`/`_variables`）、各族雅可比 `H`/`Phi`、残差 `res`、噪声 `R` |
| 输出 | 原地改写 `state->_Cov` 与 `state->_variables`（增广/边缘化还改 `_clones_IMU`/`_features_SLAM`） |
| 调用方 | `Propagator`（预测+克隆）、`UpdaterMSCKF`/`UpdaterSLAM`（更新+初始化+边缘化）、`VioManager`（触发边缘化） |
| 被调用方 | 各 `Type` 的 `update`/`clone`/`set_local_id`/`check_if_subvariable`（S1–S4） |

### 1.2 一句话概括

`StateHelper` 是 MSCKF 滤波器在 `_Cov`（全局协方差）上执行 **预测前推、卡尔曼更新、状态增广、状态边缘化** 四类块操作的统一实现。它把"误差态 EKF 协方差更新公式"翻译成"按变量 `_id` 寻址的 Eigen 块拷贝/乘法"，并用 `set_local_id` 重排机制维护滑窗增减时所有变量的连续内存布局。

---

## Section 2: 核心数据结构

### 2.1 关键接口清单

```cpp
// 文件: ov_msckf/src/state/StateHelper.h (节选)
class StateHelper {
public:
    // ===== 协方差更新 =====
    static void EKFPropagation(state, order_NEW, order_OLD, Phi, Q);  // :36  P'=ΦPΦᵀ+Q
    static void EKFUpdate(state, H_order, H, res, R);                 // :116 P-=K(PHᵀ)ᵀ ; x+=K·res

    // ===== 协方差查询/设置 =====
    static void set_initial_covariance(state, cov, order);            // :199 启动协方差注入
    static Eigen::MatrixXd get_marginal_covariance(state, small_vars);// :226 取子块协方差
    static Eigen::MatrixXd get_full_covariance(state);                // :256 全量副本

    // ===== 状态增广 / 边缘化 =====
    static void marginalize(state, marg);                            // :271 删单个变量
    static std::shared_ptr<Type> clone(state, variable);             // :341 复制变量到协方差末尾
    static bool initialize(state, new_var, H_order, H_R, H_L, R, res, mult=1.0); // :393 新变量(含χ²检验)
    static void initialize_invertible(state, new_var, H_order, H_R, H_L, R, res);// :484 可逆初始化

    // ===== 滑窗/SLAM 专用 =====
    static void augment_clone(state, last_w);                        // :579 克隆 IMU pose 入滑窗
    static void marginalize_old_clone(state);                        // :618 删最旧 clone
    static void marginalize_slam(state);                             // :631 删应边缘化 SLAM 路标
};
```

### 2.2 变量 ↔ 论文符号映射表

| 代码变量 | 论文符号 | 物理含义 | 位置 |
|---------|---------|---------|------|
| `state->_Cov` | $\mathbf{P}$ | 全局误差态协方差（活动变量） | `State.h` |
| `Phi` | $\Phi$ | 状态转移（来自 S5 `Phi_summed`） | `EKFPropagation:37` |
| `Q` | $\mathbf{Q}_d$ | 离散过程噪声（来自 S5 `Qd_summed`） | `EKFPropagation:38` |
| `H` | $\mathbf{H}$ | 测量雅可比（按 `H_order` 排列） | `EKFUpdate:116` |
| `res` | $\mathbf{r}$ | 残差 | `EKFUpdate:117` |
| `R` | $\mathbf{R}$ | 测量噪声 | `EKFUpdate:117` |
| `M_a` | $\mathbf{M}=\mathbf{P}\mathbf{H}^\top$ | 卡尔曼增益左半 | `EKFUpdate:124` |
| `S` | $\mathbf{S}=\mathbf{H}\mathbf{P}\mathbf{H}^\top+\mathbf{R}$ | 残差协方差 | `EKFUpdate:154` |
| `K` | $\mathbf{K}=\mathbf{M}\mathbf{S}^{-1}$ | 卡尔曼增益 | `EKFUpdate:162` |
| `Cov_PhiT` | $\mathbf{P}\Phi^\top$ | 前推中间量 | `EKFPropagation:80` |
| `H_L`/`H_R` | — | Givens QR 后的左/右（新变量/旧变量）系统 | `initialize:394` |

> **关键约定**：OpenVINS 的 `_Cov` **只存活动变量**（不含被边缘化掉的），且保证变量按 `_id` 在矩阵中**连续无空洞**（见 §3.4 `set_local_id` 重排）。所有块操作都依赖这一"连续布局"假设——`EKFPropagation` 入口甚至会 `assert` 非连续即崩溃（`:48-54`）。

---

## Section 3: 理论推导 → 代码逐行对照 ⭐

### 3.1 预测前推 `EKFPropagation` —— $\mathbf{P}'=\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$

#### 推导说明

误差态 EKF 预测步协方差更新（Trawny 2005 [40] Eq.103）：

$$\mathbf{P}' = \Phi\,\mathbf{P}\,\Phi^\top + \mathbf{Q}_d$$

其中 $\Phi$ 来自 S5 累加的 `Phi_summed`，仅作用于 `order_NEW==order_OLD`（即 IMU 块自身，滑窗 clone 不参与预测）。

#### 对应代码（EKFPropagation 主体）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:36-114
void StateHelper::EKFPropagation(state, order_NEW, order_OLD, Phi, Q) {
    // 1. 防御: 空数组直接崩溃
    if (order_NEW.empty() || order_OLD.empty()) { ... std::exit(...); }   // :41-44

    // 2. 校验 NEW 在内存中连续 (否则块写入错位)
    int size_order_NEW = order_NEW[0]->size();                            // :47
    for (i=0; i<order_NEW.size()-1; i++) {                                // :48
      if (order_NEW[i]->id()+order_NEW[i]->size() != order_NEW[i+1]->id())
        { ... std::exit(...); }                                           // :49-54 非连续→崩
      size_order_NEW += order_NEW[i+1]->size();                           // :55
    }
    // OLD 尺寸 (用于 Phi 列数校验)
    int size_order_OLD = ...;                                             // :59-62
    assert(size_order_NEW==Phi.rows());                                   // :65
    assert(size_order_OLD==Phi.cols());                                   // :66
    assert(size_order_NEW==Q.cols() && ==Q.rows());                       // :67-68

    // 3. 给每个 OLD 变量定位其在 Phi 中的列起点 Phi_id
    int cur=0; std::vector<int> Phi_id;
    for (var : order_OLD) { Phi_id.push_back(cur); cur+=var->size(); }    // :72-76

    // 4. Cov_PhiT = P * Phiᵀ   (逐 OLD 变量块累加)
    Eigen::MatrixXd Cov_PhiT=Zero(_Cov.rows(), Phi.rows());               // :80
    for (var : order_OLD)                                                 // :81
      Cov_PhiT.noalias() += _Cov.block(0,var->id(),_Cov.rows(),var->size())
                            * Phi.block(0,Phi_id[i],Phi.rows(),var->size()).transpose();  // :83-84

    // 5. Φ*Cov*Φᵀ + Q  =  逐 OLD 块 (Φ*Cov_PhiT_block) 累加
    Eigen::MatrixXd Phi_Cov_PhiT = Q.selfadjointView<Upper>();            // :88  (Q 对称)
    for (var : order_OLD)                                                 // :89
      Phi_Cov_PhiT.noalias() += Phi.block(0,Phi_id[i],Phi.rows(),var->size())
                                * Cov_PhiT.block(var->id(),0,var->size(),Phi.rows());  // :91

    // 6. 写回 _Cov 的三个区域 (start_id = order_NEW[0]->id())
    int start_id = order_NEW[0]->id();                                    // :95
    _Cov.block(start_id,0,   phi_size,total_size) = Cov_PhiT.transpose(); // :98  上块 (PΦᵀ)ᵀ
    _Cov.block(0,start_id,   total_size,phi_size) = Cov_PhiT;            // :99  左块 PΦᵀ
    _Cov.block(start_id,start_id, phi_size,phi_size) = Phi_Cov_PhiT;     // :100 中心块 ΦPΦᵀ+Q

    // 7. 对称正定检查: 对角线出现负值即崩溃
    Eigen::VectorXd diags=_Cov.diagonal();                                // :103
    for (i=0;i<diags.rows();i++) if (diags(i)<0) { ... std::exit; }       // :105-113
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{P}\Phi^\top$ | `Cov_PhiT += _Cov.block(...)*Phi.block(...).transpose()` | `:83-84` |
| $\Phi\mathbf{P}\Phi^\top+\mathbf{Q}_d$ | `Phi_Cov_PhiT += Phi*Cov_PhiT.block(...)`（`Phi_Cov_PhiT` 起始于 `Q`） | `:88-91` |
| 写回中心块 | `_Cov.block(start_id,start_id)=Phi_Cov_PhiT` | `:100` |
| 写回交叉块 | `_Cov.block(start_id,0)` 与 `_Cov.block(0,start_id)` | `:98-99` |

> **设计要点**：因为 `order_NEW == order_OLD`（都是 IMU 块），更新只覆盖协方差中对应 IMU 的**左上连续子块**及其与全局的交叉块。注意它**没有**重建整个 `_Cov`——只重算 IMU 相关行/列，其余变量（clone、SLAM 路标）与 IMU 的交叉协方差通过 `Cov_PhiT` 同步前推（`:98-99` 把 $\mathbf{P}\Phi^\top$ 转置写入，使交叉块也随 $\Phi$ 前推）。

### 3.2 卡尔曼更新 `EKFUpdate` —— $\mathbf{P}\mathrel{-}=\mathbf{K}(\mathbf{P}\mathbf{H}^\top)^\top$

#### 推导说明

标准 EKF 测量更新（[40]）：

$$\mathbf{S} = \mathbf{H}\mathbf{P}\mathbf{H}^\top + \mathbf{R},\quad \mathbf{K} = \mathbf{P}\mathbf{H}^\top\mathbf{S}^{-1},\quad \mathbf{P} \mathrel{-}= \mathbf{K}\mathbf{S}\mathbf{K}^\top = \mathbf{K}(\mathbf{P}\mathbf{H}^\top)^\top,\quad \delta\mathbf{x} = \mathbf{K}\mathbf{r}$$

#### 对应代码（EKFUpdate 主体）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:116-197
void StateHelper::EKFUpdate(state, H_order, H, res, R) {
    assert(res.rows()==R.rows()); assert(H.rows()==res.rows());          // :122-123
    // 1. M_a = P * Hᵀ   (逐活动变量 var 累加)
    Eigen::MatrixXd M_a = Zero(_Cov.rows(), res.rows());                  // :124
    std::vector<int> H_id;                                               // :128-132 给 H_order 定位列起点
    for (var : state->_variables) {                                      // :137  遍历所有活动变量
      Eigen::MatrixXd M_i = Zero(var->size(), res.rows());               // :139
      for (meas_var : H_order)                                          // :140
        M_i.noalias() += _Cov.block(var->id(),meas_var->id(),var->size(),meas_var->size())
                         * H.block(0,H_id[i],H.rows(),meas_var->size()).transpose();  // :142-143
      M_a.block(var->id(),0,var->size(),res.rows()) = M_i;               // :145
    }

    // 2. S = H*P_small*Hᵀ + R  (P_small = H_order 对应边际协方差)
    Eigen::MatrixXd P_small = get_marginal_covariance(state, H_order);   // :151
    Eigen::MatrixXd S(R.rows(),R.rows());                                // :154
    S.triangularView<Upper>() = H*P_small*H.transpose();                 // :155
    S.triangularView<Upper>() += R;                                      // :156

    // 3. S 稳定求逆 (LLT) → K = M_a * S⁻¹
    Eigen::MatrixXd Sinv = Identity(R.rows(),R.rows());                  // :160
    S.selfadjointView<Upper>().llt().solveInPlace(Sinv);                 // :161
    Eigen::MatrixXd K = M_a * Sinv.selfadjointView<Upper>();             // :162

    // 4. P -= K * M_aᵀ   (即 K*(PHᵀ)ᵀ = KS Kᵀ), 再对称化
    _Cov.triangularView<Upper>() -= K * M_a.transpose();                 // :166
    _Cov = _Cov.selfadjointView<Upper>();                                // :167

    // 5. 对称正定检查
    ... (同 3.1 的负值崩溃) ...                                          // :172-182

    // 6. 状态均值更新: δx = K*res, 逐变量 boxplus
    Eigen::VectorXd dx = K * res;                                        // :185
    for (i=0;i<state->_variables.size();i++)                            // :186
      state->_variables[i]->update(dx.block(var->id(),0,var->size(),1)); // :187 S1 boxplus

    // 7. 在线内参标定: 把 cam_intrinsics 的 value 同步到相机对象
    if (state->_options.do_calib_camera_intrinsics)                     // :192
      for (calib : state->_cam_intrinsics)
        state->_cam_intrinsics_cameras[calib.first]->set_value(calib.second->value());  // :194
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{M}=\mathbf{P}\mathbf{H}^\top$ | `M_i += _Cov.block(var.id, meas.id)*H.block(...).transpose()` | `:142-145` |
| $\mathbf{S}=\mathbf{H}\mathbf{P}\mathbf{H}^\top+\mathbf{R}$ | `S = H*P_small*Hᵀ + R` | `:155-156` |
| $\mathbf{K}=\mathbf{M}\mathbf{S}^{-1}$ | `K = M_a * Sinv` | `:162` |
| $\mathbf{P}\mathrel{-}=\mathbf{K}\mathbf{S}\mathbf{K}^\top$ | `_Cov -= K * M_a.transpose()` | `:166` |
| $\delta\mathbf{x}=\mathbf{K}\mathbf{r}$ | `dx = K*res; var->update(dx.block(...))` | `:185-187` |

> **设计要点**：`M_a` 是**全变量** $\times$ **残差维度** 的矩阵，第 `var->id()` 行块就是该变量对应的 $\mathbf{P}_{i\cdot}\mathbf{H}^\top$。这样卡尔曼增益 `K = M_a·S⁻¹` 一次算出，**所有变量共用一个 K**。`:187` 的 `update` 调用 S1 的 boxplus，使均值按各 Type 的几何约束（如四元数左乘）更新——这是"类型系统"在第 2 次出场（第 1 次是 S5 `set_value`）。

### 3.3 状态边缘化 `marginalize` —— 删除单变量并压缩协方差

#### 推导说明

将变量 $\mathbf{x}_m$ 从 $\mathbf{x}=[\mathbf{x}_1,\mathbf{x}_m,\mathbf{x}_2]$ 中删去，协方差变为（Schur 补/条件边缘化在 MSCKF 中已通过测量完成，此处只是**物理删除**：

$$
\mathbf{P}' =
\begin{bmatrix}
\mathbf{P}_{11} & \mathbf{P}_{12} \\
\mathbf{P}_{21} & \mathbf{P}_{22}
\end{bmatrix}
$$

删 $\mathbf{x}_m$ 后保留 $[\mathbf{x}_1, \mathbf{x}_2]$ 的对应块（代码注释 `:280-289` 图示）。

#### 对应代码（marginalize 主体）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:271-339
void StateHelper::marginalize(state, marg) {
    // 1. 防御: marg 必须在 _variables 内 (不支持子变量边缘化!)
    if (find(_variables.begin(),_variables.end(),marg)==_variables.end())    // :274
      { ... std::exit; }                                                     // :275-278

    // 2. 计算三区尺寸
    int marg_size=marg->size(), marg_id=marg->id();                          // :293-294
    int x2_size = _Cov.rows() - marg_id - marg_size;                         // :295

    // 3. 构造新协方差 (去掉 marg 行/列)
    Eigen::MatrixXd Cov_new(_Cov.rows()-marg_size, _Cov.rows()-marg_size);   // :297
    Cov_new.block(0,0,        marg_id, marg_id) = _Cov.block(0,0,marg_id,marg_id);                 // :300 P11
    Cov_new.block(0,marg_id,  marg_id, x2_size) = _Cov.block(0,marg_id+marg_size,marg_id,x2_size); // :303 P12
    Cov_new.block(marg_id,0,  x2_size, marg_id) = Cov_new.block(0,marg_id,marg_id,x2_size).transpose(); // :306 P21
    Cov_new.block(marg_id,marg_id, x2_size,x2_size)=_Cov.block(marg_id+marg_size,marg_id+marg_size,x2_size,x2_size); // :309 P22

    state->_Cov = Cov_new;                                                   // :313
    assert(_Cov.rows()==Cov_new.rows());                                     // :315

    // 4. 重排剩余变量 id: 排在 marg 之后的变量向前挪 marg_size
    std::vector<Type> remaining;
    for (var : _variables) {                                                 // :320
      if (var != marg) {
        if (var->id() > marg_id)
          var->set_local_id(var->id() - marg_size);                          // :325 ★ 连续化重排
        remaining.push_back(var);
      }
    }
    marg->set_local_id(-1);                                                  // :335 标记已删
    state->_variables = remaining;                                           // :338
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| $\mathbf{P}_{11}$ | `Cov_new.block(0,0,marg_id,marg_id)` | `:300` |
| $\mathbf{P}_{12}$ | `Cov_new.block(0,marg_id,...)` | `:303` |
| $\mathbf{P}_{21}=\mathbf{P}_{12}^\top$ | `Cov_new.block(marg_id,0,...).transpose()` | `:306` |
| $\mathbf{P}_{22}$ | `Cov_new.block(marg_id,marg_id,...)` | `:309` |
| id 重排（保持连续） | `var->set_local_id(var->id()-marg_size)` | `:325` |

> **设计要点（连续布局维护）**：边缘化会留下"空洞"（被删变量原来占的行/列）。`:325` 用 `set_local_id` 把所有排在后面的变量偏移 `-marg_size`，使协方差在**逻辑上重新连续**。这保证了后续 `EKFPropagation` 的连续 `assert`（`:48`）始终成立。注意代码注释 `:318` 明确说**尚不支持子变量边缘化**——只能删顶层 `Type`（如整个 clone、整个 SLAM 路标）。

### 3.4 状态克隆 `clone` —— 把变量（含协方差）复制到矩阵末尾

#### 推导说明

MSCKF 随机克隆 [20]：在时刻 $t_k$ 把当前 IMU pose 的**拷贝**追加到状态末尾，协方差为原变量协方差块的复制（克隆态与母态初始完全相关，后续经测量解耦）。

#### 对应代码（clone 主体）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:341-391
std::shared_ptr<Type> StateHelper::clone(state, variable_to_clone) {
    int total_size = variable_to_clone->size();                             // :344
    int old_size = _Cov.rows();                                             // :345
    int new_loc  = _Cov.rows();                                             // :346 新位置=原末尾

    // 1. 协方差末尾扩 total_size 行/列
    _Cov.conservativeResizeLike(Zero(old_size+total_size, old_size+total_size)); // :349

    // 2. 在 _variables 中定位要克隆的变量 (支持子变量 check_if_subvariable)
    std::shared_ptr<Type> new_clone = nullptr;
    for (k=0; k<_variables.size(); k++) {                                   // :356
      Type* type_check = _variables[k]->check_if_subvariable(variable_to_clone); // :360 S4 递归
      if (_variables[k]==variable_to_clone) type_check=_variables[k];      // :361-362
      else if (type_check != variable_to_clone) continue;                   // :363-364

      int old_loc = type_check->id();                                       // :368
      // 3. 复制协方差三块 (中心自协方差 + 与全局的交叉)
      _Cov.block(new_loc,new_loc, total_size,total_size) = _Cov.block(old_loc,old_loc,total_size,total_size); // :371
      _Cov.block(0,new_loc, old_size,total_size)       = _Cov.block(0,old_loc,old_size,total_size);          // :372
      _Cov.block(new_loc,0, total_size,old_size)       = _Cov.block(old_loc,0,total_size,old_size);         // :373

      // 4. 生成新克隆 (深拷贝 value/fej), 设到新位置
      new_clone = type_check->clone();                                      // :376 S1 虚 clone
      new_clone->set_local_id(new_loc);                                     // :377 ★ id=新末尾
      break;
    }
    if (new_clone==nullptr) { ... std::exit; }                              // :382-386
    state->_variables.push_back(new_clone);                                 // :389
    return new_clone;
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 协方差扩维 | `conservativeResizeLike(Zero(old+total, old+total))` | `:349` |
| 克隆自协方差块 | `_Cov.block(new_loc,new_loc)=_Cov.block(old_loc,old_loc)` | `:371` |
| 克隆交叉块 | `_Cov.block(0,new_loc)=_Cov.block(0,old_loc)` (左/上) | `:372-373` |
| 新克隆 id | `new_clone->set_local_id(new_loc)` | `:377` |

> **设计要点**：`clone` 是 `augment_clone`（S5 `:137` 调用）的底层。它**不区分 IMU/SLAM/clone**——任何 `Type`（含子变量，靠 `check_if_subvariable` S4 递归定位）都能克隆。克隆体通过 `type_check->clone()`（S1 虚函数）深拷贝，携带此刻 `value==fej`（即 S5 提到的**克隆 FEJ 锚定**）。

### 3.5 滑窗克隆 `augment_clone` 与时间偏移雅可比

#### 对应代码

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:579-616
void StateHelper::augment_clone(state, last_w) {
    // 防御: 同时间戳不可再克隆
    if (_clones_IMU.find(_timestamp)!=_clones_IMU.end()) { ... std::exit; } // :582-585
    // 克隆当前 IMU pose 到协方差末尾
    Type* posetemp = clone(state, state->_imu->pose());                    // :589 调 3.4
    PoseJPL* pose = dynamic_pointer_cast<PoseJPL>(posetemp);              // :592
    if (pose==nullptr) { ... std::exit; }                                 // :593-596
    _clones_IMU[_timestamp] = pose;                                       // :599 ★ 入滑窗

    // 时间偏移在线标定: 克隆 pose 对 dt 的雅可比 (Li & Mourikis 2013 [27])
    if (state->_options.do_calib_camera_timeoffset) {                     // :604
      Eigen::Matrix<double,6,1> dnc_dt=Zero(6,1);
      dnc_dt.block(0,0,3,1) = last_w;                                     // :607 角速度→姿态漂移
      dnc_dt.block(3,0,3,1) = state->_imu->vel();                         // :608 速度→位置漂移
      _Cov.block(0,pose->id(),_Cov.rows(),6)   += _Cov.block(0,dt_id,_,1)*dnc_dt.transpose(); // :611
      _Cov.block(pose->id(),0,6,_Cov.rows())   += dnc_dt*_Cov.block(dt_id,0,1,_).transpose(); // :613
    }
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| 克隆入滑窗 | `_clones_IMU[_timestamp]=pose` | `:599` |
| $\partial \mathbf{n}_c / \partial t_d = [\boldsymbol{\omega}; \mathbf{v}]$ | `dnc_dt=[last_w; vel]` | `:607-608` |
| 协方差增广偏移雅可比 | `_Cov.block(0,pose.id) += _Cov.block(0,dt_id)*dnc_dtᵀ` | `:611-613` |

> **设计要点**：`augment_clone` 是 S5 `propagate_and_clone` 在每帧末尾调用的滑窗**生长点**。当开启相机-IMU 时间偏移在线标定（`do_calib_camera_timeoffset`，对应 [27]）时，克隆姿态是时间偏移 $t_d$ 的函数，需把该雅可比注入协方差交叉块（`:611-613`）。代码注释 `:610` 自承此处本可复用 `EKFPropagation` 但暂未重构——属已知技术债。

### 3.6 新变量初始化 `initialize` + `initialize_invertible`（Givens QR）

#### 推导说明

SLAM 路标首次被观测到时，需把新变量加入状态。为避免可观性缺失，用 **Givens 旋转将测量系统 QR 分解**（Solà 2017 [47] / MSCKF-SLAM [32]）：把含新变量的行旋到左上，分离成
- **可逆初始化系统** $\mathbf{H}_f \delta\mathbf{x}_f = \mathbf{H}_L \delta\mathbf{x}_L + \mathbf{r}_L$（直接求新变量）
- **零空间投影更新系统** $\mathbf{H}_{up}\delta\mathbf{x}_{old}=\mathbf{r}_{up}$（投影到旧变量更新）

先做 Mahalanobis $\chi^2$ 检验（异常则拒绝初始化），再分别处理两部分。

#### 对应代码（initialize 主体）

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:393-482
bool StateHelper::initialize(state, new_variable, H_order, H_R, H_L, R, res, chi_2_mult=1.0) {
    // 1. 防御: new_variable 不能已在 _variables 中
    if (find(_variables.begin(),_variables.end(),new_variable)!=_variables.end()) { ...std::exit;} // :398-402
    // 2. 校验 R 各向同性对角 (R(0,0) 全等且非对角为0)
    for (r,c) { if(r==c && R(0,0)!=R(r,c)) exit; else if(r!=c && R(r,c)!=0) exit; }  // :408-420

    // 3. Givens QR: 把 H_L (含新变量列) 旋成上三角, 同步作用于 H_R/res/R
    Eigen::JacobiRotation<double> G;
    for (n=0; n<H_L.cols(); n++)                                            // :430
      for (m=H_L.rows()-1; m>n; m--) {                                      // :431
        G.makeGivens(H_L(m-1,n), H_L(m,n));                                 // :433
        H_L.block(m-1,n,2,H_L.cols()-n).applyOnTheLeft(0,1,G.adjoint());    // :437
        res.block(m-1,0,2,1).applyOnTheLeft(0,1,G.adjoint());               // :438
        H_R.block(m-1,0,2,H_R.cols()).applyOnTheLeft(0,1,G.adjoint());      // :439
      }

    // 4. 分离: 上 new_var_size 行=初始化系统, 余下=投影更新系统
    Hxinit=H_R.block(0,0,new_var_size,H_R.cols());  H_finit=H_L.block(0,0,new_var_size,new_var_size); // :445-446
    resinit=res.block(0,0,new_var_size,1);          Rinit=R.block(0,0,new_var_size,new_var_size);     // :447-448
    Hup=H_R.block(new_var_size,0,...);  resup=res.block(new_var_size,0,...);  Rup=R.block(...);        // :451-453

    // 5. χ² Mahalanobis 检验 (投影系统部分)
    P_up = get_marginal_covariance(state, H_order);                          // :459
    S = Hup*P_up*Hup.transpose()+Rup;                                       // :462
    chi2 = resup.dot(S.llt().solve(resup));                                  // :463
    boost::math::chi_squared chi_dist(res.rows());                           // :466
    double chi2_check = quantile(chi_dist, 0.95);                           // :467 95% 阈值
    if (chi2 > chi_2_mult * chi2_check) return false;                        // :468-469 拒绝

    // 6. 可逆初始化新变量 + 投影部分做 EKF 更新
    initialize_invertible(state, new_variable, H_order, Hxinit, H_finit, Rinit, resinit); // :475
    if (Hup.rows()>0) EKFUpdate(state, H_order, Hup, resup, Rup);            // :478-480
    return true;
}
```

#### `initialize_invertible` 核心（`:484-577`）

```cpp
void StateHelper::initialize_invertible(state, new_variable, H_order, H_R, H_L, R, res) {
    // 同款防御 + 各向同性 R 校验 (:488-511)
    // 1. M_a = P*H_Rᵀ  (与 EKFUpdate 同构, 但用 H_R)
    ... M_a.block(var.id,0,...) = M_i (含 H_R) ...                          // :519-541
    // 2. M = H_R*P_small*H_Rᵀ + R
    M.triangularView<Upper>()=H_R*P_small*H_R.transpose(); M+=R;            // :549-551
    // 3. 新变量协方差 P_LL = H_L⁻¹ * M * H_L⁻ᵀ
    H_Linv = H_L.inverse();                                                 // :556
    P_LL = H_Linv * M.selfadjointView<Upper>() * H_Linv.transpose();        // :557
    // 4. 扩协方差, 写交叉块 + 新变量自协方差
    oldSize=_Cov.rows();
    _Cov.conservativeResizeLike(Zero(oldSize+size, oldSize+size));          // :561
    _Cov.block(0,oldSize,oldSize,size).noalias() = -M_a*H_Linv.transpose(); // :562 交叉 (负)
    _Cov.block(oldSize,0,size,oldSize) = _Cov.block(0,oldSize,...).transpose(); // :563
    _Cov.block(oldSize,oldSize,size,size) = P_LL;                           // :564 自协方差
    // 5. 更新新变量均值 (应≈0, 因已用 CGN 解过初值)
    new_variable->update(H_Linv * res);                                     // :568
    new_variable->set_local_id(oldSize);                                    // :571 ★ 入末尾
    state->_variables.push_back(new_variable);                              // :572
}
```

#### 对照注释

| 公式项 | 对应代码 | 行号 |
|--------|---------|------|
| Givens QR 分离 | `G.makeGivens(...); applyOnTheLeft(...)` | `:430-439` |
| $\chi^2=\mathbf{r}_{up}^\top\mathbf{S}^{-1}\mathbf{r}_{up}$ | `chi2=resup.dot(S.llt().solve(resup))` | `:463` |
| 95% 阈值 | `quantile(chi_squared(res.rows()),0.95)` | `:466-467` |
| $\mathbf{P}_{LL}=\mathbf{H}_L^{-1}\mathbf{M}\mathbf{H}_L^{-\top}$ | `P_LL=H_Linv*M*H_Linvᵀ` | `:557` |
| 交叉块 $\mathbf{P}_{iL}=-\mathbf{M}_a\mathbf{H}_L^{-\top}$ | `_Cov.block(0,oldSize)=-M_a*H_Linvᵀ` | `:562` |

> **设计要点**：`initialize` 是 `EKFUpdate` 的"超集"——它先 QR 分离出可逆部分直接解新变量（`initialize_invertible`），再把零空间投影部分交给普通 `EKFUpdate`。`chi_2_mult`（默认 1.0）可调检验严格度。`:568` 的 `update` 再次调用 S1 boxplus。SLAM 路标、地图点首次可见都走这条路径。

### 3.7 滑窗/SLAM 专用边缘化

```cpp
// 文件: ov_msckf/src/state/StateHelper.cpp:618-645
void StateHelper::marginalize_old_clone(state) {
    if ((int)_clones_IMU.size() > _options.max_clone_size) {               // :619 滑窗容量超限
      double marginal_time = state->margtimestep();                        // :620 最旧 clone 时刻
      std::lock_guard<std::mutex> lock(state->_mutex_state);              // :622 线程安全
      marginalize(state, _clones_IMU.at(marginal_time));                   // :624 调 3.3
      _clones_IMU.erase(marginal_time);                                    // :627 删指针
    }
}

void StateHelper::marginalize_slam(state) {
    auto it=_features_SLAM.begin();
    while (it!=_features_SLAM.end()) {                                      // :636
      if (it->second->should_marg && (int)it->first > 4*max_aruco_features) { // :637 标记+非aruco
        marginalize(state, it->second);                                    // :638 调 3.3
        it = _features_SLAM.erase(it);                                      // :639 删指针
      } else it++;
    }
}
```

#### 对照注释

| 调用 | 触发条件 | 行号 |
|------|---------|------|
| `marginalize_old_clone` | `_clones_IMU.size() > max_clone_size` | `:619` |
| `marginalize_slam` | `should_marg==true` 且非 aruco 路标 | `:637` |

> **设计要点**：这两个是 `marginalize`（§3.3）的**策略封装**——决定"何时删谁"。`marginalize_old_clone` 维持滑窗长度 `max_clone_size`（MSCKF 核心：用固定数量历史 clone 提供多帧约束 [20]）；`marginalize_slam` 删除已丢失/应边缘化的 SLAM 路标，但**保留 aruco 标签路标**（`< 4*max_aruco_features` 的 id 段，见 `:637` 的保护条件）。

---

## Section 4: 函数调用链

### 4.1 入口分类

```
【预测类】  Propagator::propagate_and_clone(:33)
   └─ EKFPropagation(:36)                          [ΦPΦᵀ+Qd]
   └─ augment_clone(:579) → clone(:341)           [滑窗生长]

【更新类】  UpdaterMSCKF::update / UpdaterSLAM::update
   └─ EKFUpdate(:116)                             [K·res]
   └─ initialize(:393) → initialize_invertible(:484) + EKFUpdate(:116)
   └─ marginalize(:271) ← marginalize_old_clone(:618) / marginalize_slam(:631)

【查询类】  get_marginal_covariance(:226) / get_full_covariance(:256) / set_initial_covariance(:199)
```

### 4.2 完整调用树（含 S1–S5 依赖）

```
StateHelper::EKFUpdate @ :116
  ├─ get_marginal_covariance(state, H_order) @ :226     → 取 H_order 子块协方差
  ├─ _Cov.block(...) 读/写                              → 卡尔曼增益 + 协方差更新
  └─ var->update(dx.block(...)) @ :187                  → S1 Type::update (boxplus)
        └─ JPLQuat::update / Vec::update / PoseJPL / IMU (S2–S4)

StateHelper::clone @ :341
  ├─ _variables[k]->check_if_subvariable(variable) @ :360  → S4 递归定位
  ├─ type_check->clone() @ :376                           → S1 Type::clone (深拷贝)
  └─ new_clone->set_local_id(new_loc) @ :377              → S1 连续 id 设置

StateHelper::initialize @ :393
  ├─ Givens QR (Eigen JacobiRotation) @ :430
  ├─ get_marginal_covariance @ :459
  ├─ initialize_invertible @ :484
  │    └─ new_variable->update(H_Linv*res) @ :568         → S1 boxplus
  └─ EKFUpdate @ :478 (投影部分)

StateHelper::augment_clone @ :579
  ├─ clone(state, _imu->pose()) @ :589                   → 3.4
  └─ (可选) 时间偏移雅可比注入 @ :611                     → [27]
```

### 4.3 关键函数速查

| 函数名 | 文件:行号 | 功能 | 输入 | 输出 | 复杂度 |
|--------|-----------|------|------|------|--------|
| `EKFPropagation` | `StateHelper.cpp:36` | 预测前推 | `Phi, Q, order` | 改 `_Cov` | O(N²·n_old) |
| `EKFUpdate` | `StateHelper.cpp:116` | 卡尔曼更新 | `H, res, R, H_order` | 改 `_Cov`+`_variables` | O(N²·n_meas) |
| `marginalize` | `StateHelper.cpp:271` | 删单变量 | `marg` | 压缩 `_Cov`+重排 id | O(N²) |
| `clone` | `StateHelper.cpp:341` | 复制变量到末尾 | `variable` | 新 `Type*` | O(N·size) |
| `initialize` | `StateHelper.cpp:393` | 新变量(χ²+QR) | `H_R,H_L,R,res` | 增广 `_Cov` | O(N²+N·qr) |
| `initialize_invertible` | `StateHelper.cpp:484` | 可逆初始化 | `H_R,H_L,R,res` | 增广 `_Cov` | O(N²) |
| `augment_clone` | `StateHelper.cpp:579` | 滑窗克隆 | `last_w` | 新 clone | O(N·6) |
| `marginalize_old_clone` | `StateHelper.cpp:618` | 删最旧 clone | — | 压缩 | 调 `marginalize` |
| `marginalize_slam` | `StateHelper.cpp:631` | 删 SLAM 路标 | — | 压缩 | 调 `marginalize` |

---

## Section 5: 关键参数清单

| 参数名 | 默认值 | 定义位置 | 类别 | 含义 | 敏感度 |
|--------|--------|---------|------|------|--------|
| `max_clone_size` | 配置 yaml | `StateOptions` | 滑窗 | 最大历史 clone 数 | ★★★ 高 |
| `max_aruco_features` | 配置 yaml | `StateOptions` | SLAM | aruco 路标保护范围 | ★ 低 |
| `do_calib_camera_timeoffset` | false | `StateOptions` | 时序 | 时间偏移在线标定开关 | ★★ 中 |
| `do_calib_camera_intrinsics` | false | `StateOptions` | 标定 | 内参在线标定开关 | ★★ 中 |
| `chi_2_mult` | 1.0 | `initialize` 调用方 | 一致性 | χ² 检验倍率 | ★★ 中 |

### 参数影响分析

- **`max_clone_size`**：直接决定 MSCKF 滑窗长度（`:619`）。过大→计算量线性增、但多帧约束更紧；过小→约束不足、精度降。是 MSCKF 与 EKF-SLAM 的分水岭参数 [20]。
- **`chi_2_mult`**：缩放 95% χ² 阈值（`:467-469`）。调大→更容忍异常观测（漏检少但可能接纳 outlier）；调小→更严格（鲁棒但易拒绝正确初始化）。
- **`do_calib_camera_timeoffset`**：开启时 `augment_clone` 注入偏移雅可比（`:611-613`，对应 [27]），否则跳过——省算力但失去时间同步自适应。

---

## Section 6: 问题与改进方向

| # | 问题 | 严重度 | 触发条件 | 当前表现 | 改进建议 |
|---|------|-------|---------|---------|---------|
| 1 | `marginalize` 不支持子变量 | ★ 中 | 试图删 PoseJPL 内单元素 | `assert`/注释表明未实现，直接 exit | 当前靠"删整个 clone/路标"规避；如需部分边缘化需扩展 |
| 2 | `augment_clone` 偏移雅可比未复用 `EKFPropagation` | ★ 低 | `do_calib_camera_timeoffset` | 代码注释 `:610` 自承技术债，手动块加法 | 重构为调用统一前推，降低维护风险 |
| 3 | 负值对角线硬崩溃 | ★ 中 | 数值病态（协方差非正定） | `EKFPropagation:105` / `EKFUpdate:172` 直接 `std::exit` | 可降级为重置协方差或告警，避免整系统崩溃 |
| 4 | `initialize` 要求 R 各向同性 | ★ 低 | 传入非对角/非均匀 R | `:408-420` 直接 exit | 放宽校验或支持一般 R（需改 QR 推导） |

---

## Section 7: 同类实现对比

| 特征维度 | OpenVINS `StateHelper` | VINS-Mono | OKVIS | 常规 MSCKF [20] |
|---------|----------------------|-----------|-------|-----------------|
| 协方差表示 | 单 `_Cov` 活动变量块 | 边缘化先验(Hessian) | 误差态 `_Cov` | `_Cov` |
| 滑窗增广 | `clone`+`augment_clone` 显式 | 无(图优化) | 有 | 有 |
| 边缘化 | `marginalize`+`set_local_id` 重排 | Schur 补 | Schur 补 | 删 clone |
| 新变量 | Givens QR + χ² (`initialize`) | 三角化 | 三角化 | 无(纯 MSCKF) |
| 时间偏移 | `augment_clone` 雅可比 [27] | 标定 | 标定 | 标定 |
| 类型系统耦合 | `var->update`/`clone` 虚函数 (S1–S4) | 直接 vector | 直接 vector | 直接 vector |

### 设计取舍分析

- 选择 **单 `_Cov` + 连续 id + `set_local_id` 重排** 是因为：EKF 框架需要 O(N²) 原地块操作，连续布局让 Eigen 块寻址 O(1)，且 `EKFPropagation` 可 `assert` 连续性防错。
- 选择 **`clone` 通用化（支持子变量 `check_if_subvariable`）** 是因为：IMU、SLAM 路标、clone 都是 `Type` 子类，统一克隆逻辑避免重复代码；但代价是 `marginalize` 必须限制为顶层变量（见 §6 #1）。
- 在 **SLAM 路标首次可见** 场景用 **Givens QR + χ²** 而非直接三角化，是因为：QR 把新变量系统分离为可逆+投影两部分，显式 χ² 检验可拒绝异常初始化（[32] MSCKF-SLAM 一致性保障）。
