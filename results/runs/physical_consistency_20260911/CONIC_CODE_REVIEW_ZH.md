# 完整联合锥通路只读代码复核

结论：未发现不等式方向、缩放、offset、端点固定或结果分类的阻断错误。新通路保留全部独立振幅方向，使用原JointProblem映射；它不是已有振幅的hull。需要限定一处理论/元数据表述：删除自由高能R对应的Gram后，严格散射内域有Schur完成，但闭边界仍有range相容条件，不能无条件称为原完整闭可行集的精确投影。最终原式audit保留完整Gram检查，因此当前代码不会据此误报原问题可行。

本复核只读代码与给定设计，没有运行新计算、求解、测试或修改src/tests。两组compile pack检查的通过来自主线；完整suite及真实M50/L10运行仍由主线负责。

## 各族映射

`conic.py:25–28`令原变量z=offset+D*x，并存储A_sys=−A_map*D、b_sys=b_map+A_map*offset。因此Clarabel的s=b_sys−A_sys*x恰还原b_map+A_map*z，常数与负号一致。

| 族 | 实际cone slack | 核对 |
|---|---|---|
| 散射 | (1+y, √2 x, √(2t)y, 1−y) | 平方差2(2y−x²−ty²)，SOC头非负，不引入负支；保留t项 |
| χ | (1,Cg*z/ε) | 同一分组Euclidean球 |
| L4局部 | (1+w_i,2r_i,w_i−1) | 等价w_i≥r_i²≥0 |
| L4整体 | (1,w) | 与局部合起来等价sum r_i⁴≤1 |
| Gram | g+J*z | 原svec顺序(11,√2·12,22,√2·13,√2·23,33)，与设计给定PSD上三角按列规范一致 |
| FF | (1,fc+F*z) | 逐行实/虚部重排与原归一化FF圆盘一致 |
| FESR | (1−moment,1+moment) | 由(-W,+W)、常数(1+target,1−target)恢复，同一四个绝对误差 |
| 领先高能 | (1+y,√2 x,1−y) | 等价2y−x²≥0，使用已配置的同五行 |

L4的r直接取现有actual-density坐标/B；ρ2非对角权重已由JointProblem坐标承担，无额外π或2。B/nd^(1/4)仅是正列缩放；辅助w保留为独立变量。完整L4球没有改成逐项箱界。Gram的正对角congruence仍来自同一P，没有对PSD作额外放松或τ平移。

## 端点、变量与恢复

`compile_joint_cones`只删除endpoint_pivots，把其固定值写入offset；明确禁止删除density列。T0固定坐标为0，两个FF通道各固定(1,0)，由既有精确完成解释依赖显示值。其它活动C、独立ImF和矩有界R都保留，辅助密度变量位于原活动变量之后。

返回时先按完整columns/scale恢复extended向量，再取前len(P.active)个原联合坐标，正确丢弃辅助w；没有把求解器压缩向量的前段误当原坐标。`P.raw`继续实施C0精确零与FF依赖变量完成。实际组装使用浮点系数，仍须按设计经过原式及后续解析源误差核验。

**自由高能R的闭边界限制：** 对块矩阵[[A,b],[bᵀ,R]]且R无上界，有限PSD完成需要A≥0且b属于A的值域；A严格正定时range条件自动满足，Schur下界有限。编译器仅保留P.keep的Gram，故在自由R所对应的闭散射边界上可能包含不满足range条件的约化点。不能据此声称整个闭约化cone无条件等于原完整问题的投影，也不应无额外论证把联合问题称为该投影的闭包。

现有 `joint_audit`在delta>0时提升自由高能谱，随后检查所有原Gram；无法完成的点不会接受。delta=0而相容的点也可能需要伪逆Schur提升，当前严格delta恢复分支未必能取到它，因此一次恢复失败不证明该闭边界点或整个问题不可行。建议将conic_model的“algebraically equivalent cones”限定为保留各族的代数提升，并明确自由谱消元的严格内域/边界相容范围；不需要为此扩大当前生产求解目标。

## current、reference与派发

`gauge.initialize`在新三行派发之前已经通过support_data/current_data校验H/current路径与分辨率，核对UV配置和endpoint/tail旗标，并检查输入C与joint.npz完整一致。

conic初始化从同一个当前H/current构造JointProblem。若使用旧Newton Phase I的reference_point，先要求两套preparation路径一致，只读取reference作可逆坐标与正congruence，不读取旧z、tau或传播值作为锥解。无需继承旧χ/B可行性；这些参数取当前args且重新进入所有锥。原reference仅需满足构造所需的正对角条件，否则明确失败，不是物理不可行结论。`warm_start_used=False`与实际调用一致，Clarabel没有接收旧点作为热启动。

## 零目标与状态

求解器P矩阵和目标向量都为0；代码保留native status及候选，但不把Solved作为接受判据。只有原 `joint_audit`检查后的点进入joint_result。

method明确phase_I=True、support_search=False、infeasibility_certified=False、representative_center_converged=False。`resolution_report`因此强制支持方向[0,0]；joint_result不会据零目标gap产生support_optimality，初始化器还显式将其置False，所以不会生成代表selection。原生PrimalInfeasible状态没有转成严格排除；非有限候选抛错、有限但原式失败则保留joint_candidate。新的phase_I_model为full-conic-feasibility-v1，与Newton缓存模型分开，不能伪装成已保存的Newton中心。

建议主线在首次真实返回后继续按现有顺序区分native状态、原式有限可行性与解析源核验，不预宣称此算法会优于此前路径。

## 输入指纹

复核时间：2026-09-11T09:00:42.246548+00:00。

| 文件 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/conic.py` | `087b0d3d41b64b9e4466cefaa72f926e2ccd9d92591056e75137e02571109667` |
| `src/smatrix_bootstrap/gauge.py` | `6c0e1e626744973a07e1735740ce56722ef273680a14c2769a6ac5cce3d36b5d` |
| `src/smatrix_bootstrap/operators.py` | `5850eb2add0c8a8638b7b92bf4868b4d0a9015e74be3737a0b9f4006ae1a1c0c` |
| `src/smatrix_bootstrap/endpoints.py` | `6b4dacd438e750ca61f285744cc0a152d16e3090581e22c77962f23b867005db` |
| `src/smatrix_bootstrap/certificates.py` | `726151e3225eda13f82260d589b4a503c1ebd3afdbdddd05831e440a7efc66a4` |
| `src/smatrix_bootstrap/__init__.py` | `9e65f9de80eaceb49fe41feb6d29c8335cafa24f7fa588f46d1b63cd66ff6cbb` |
| `tests/test_mainline.py` | `fc003ed7721e477bedd8fb9a173034d1557a6b435a04cdbc963b54ff63773c4d` |
| `results/runs/physical_consistency_20260911/CONIC_DESIGN_ZH.md` | `38d5de36ef7787858a89f060e1ee6fd93793540a7827fed09f3df79ed61b24e3` |
