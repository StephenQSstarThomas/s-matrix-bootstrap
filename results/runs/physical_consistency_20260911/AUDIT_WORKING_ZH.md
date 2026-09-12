# He–Kruczenski 2309.12402v3：原型机、ρ机制与物理一致性审计

2026-09-11，工作稿。运行中的问题以 [PROGRESS.json](PROGRESS.json) 为准。本稿不能作为全部计算已结束的记录。

## 当前结论与验收对象

**3582盘按几何规则预选的IR对照与UV三点链已完成，但正确且稳健的ρ复现仍有Major差异。** 中心、支持、解析和冻结求值均已交付，同一185能量×3波的5MeV网格全部通过、0严格违例。max η=1包含阈值，不是continuum证书；质量分散、强非弹性、S0偏差及原基线L依赖仍须分别判断。正式Watsonian首轮因线性残差检查停止，非不可行结论；修补已通过106项测试并重启，同一固定目标已取得多个新标准中心；支持间隙仍待收敛，尚无已完成的 Watson 物理结果。

历史3123盘模型的非空证据是 [interlaced_phase_I_07](interlaced_phase_I_07/report.json)的τ=−.001098097655及[解析核验](interlaced_witness_analytic/report.json)。其IR/UV冻结点和[旧原生相位](interlaced_native_phases/phases.json)、独立5MeV违例全部保留；这些数值不移用于3582盘。

3582盘的[完整联合见证](fine_witness_analytic/report.json)、[IR中心/支持](fine_IR_ref_02/report.json)及[解析核验](fine_IR_ref_analytic/source_audit.json)、UV [tip](fine_tip_support_02/report.json)/[ref](fine_UV_ref_support_02/report.json)/[mid](fine_UV_mid_support_02/report.json)的中心、支持和解析均通过。[fine_IR_selected](fine_IR_selected/selection.json)与[fine_UV_selected](fine_UV_selected/selection.json)在相位检查前固定完整C；fine_C1/C3及原生/网格结果已经交付。先前[旧见证6项S0违例](fine_witness_replay/source_audit.json)仍保留为本模型新见证之前的失败。

原文§4.3把P1相移穿过π/2的能量作为这里的共振读数，Fig.9–10比较预选代表，Appendix A/Fig.11比较分辨率。原文不是“所有可行解必有同一个ρ极点”的数学定理。本文不把唯一极点、无限自旋或连续证书强加为原有限原型的新前置条件；但对于已经发现的具体违例，也不能回避或用画图平滑消除。

Fig.7只有IR输入且未出现合适ρ，是预期对照。正确的主线是：同一散射架构加入手征条件改善S波；再通过电流Gram、FF归一化、QCD谱矩及高能FF界约束P波；最后比较完整预选振幅、谱与稳定性。看到Fig.7无ρ本身不能判断失败。

[八点中心重验](CENTER_RECHECK_RESULTS.md)已完成：fine三点和原F五配置均零迭代通过新检查点标准，完整C、ImF与全部R不变。当前数值问题是同一Watson目标的支持收敛速度，不是把已交付几何点再误称为未收敛。

## 逐项对齐和差异等级

原生1500盘、历史3123盘和当前3582盘分别标识。下文保留旧数值，并单列3582盘最新结果，不把它们混成一条收敛或代表轨迹。

| 对象 | 已完成的证据 | 当前差异 |
|---|---|---|
| A1–A3：底层振幅 | 完整3876系数的共同解析族、交叉/同位旋/投影、解析行误差与独立角积分 | 有效的显式有限完成；未恢复作者未公开的全部插值细节 |
| B1：中心和支持 | 3582盘IR/UV三点全链通过，完整C相位前冻结 | 中心属于实际约化障碍的数值中心，不证明完整QCD振幅唯一 |
| B2–B3 / Fig.3–4 | 纯散射强收缩方向见证、手征右端支持界 | **Major**：完整新区域未交付；手征范围偏大 |
| C1–C3 / Fig.5–7 | fine_IR及fine_C1/C3交付，旧三色保留 | fine P1末端8.07855°、S0零点括区(.4,.5)；IR对照不替代全部新区域 |
| D1–D3：IR/UV耦合 | 3582盘联合见证及全部几何代表通过；原UV独立排除fine IR C | upper=−.1484093366072121，不借新FF端点；只排除此完整C |
| E1 / Fig.8 | 原基线同一截面的收缩支持证书 | **Major**：只见轻微不对称，未复现图中的强不对称 |
| E2–E3 / Fig.9–10 | 3582盘三点原生/同5MeV网格交付，网格零严格违例 | **Major**：质量分散、min ηP1约.20–.31、S0偏差；有限网格通过非连续认证 |
| F1–F3 / Fig.11 | 原生基线五组完整支持、核验、代表及求值比较 | **Major**：L跨度显著大于原图；M跨度与原图同量级，须分别判断 |

[原生五组完整比较](../sequential_reproduction_20260910/F_five_configs_verified_comparison/resolution.json)已包括(50,8)、(50,10)、(50,12)、(45,10)、(60,10)。固定M50的P1读数L跨度为46.47MeV，原图约12.54MeV；固定L10的M跨度为54.06MeV，原图约54.35MeV。后者同量级，不能笼统称“M不收敛”；它也不构成M→∞收敛证明。这组结果属于原生基线，不能与加强模型的新见证拼成同一组分辨率比较。

旧1500盘原生基线三代表的90°读数约803.49、701.26、696.63MeV（tip、mid、ref）；同C在ρ窗口加点后约802.01、697.13、688.38MeV，最大η约1.093、1.119、1.141。旧数值及独立角积分违例全部保留，不与以下3123盘结果混用。

[3123盘模型的原生求值](interlaced_native_phases/phases.json)中，IR对照P1末端为8.0989°。UV三点的完整C先冻结，再得到：

| UV代表 | 原生P1 90°读数（MeV） | min ηP1 | 所报原生能区末端δS0（°） |
|---|---:|---:|---:|
| tip | 751.9317 | .304105 | 190.22 |
| mid | 680.2284 | .214759 | 202.06 |
| ref | 678.5991 | .196052 | 209.74 |

不能只选tip较接近770MeV就宣称成功。质量分散、强非弹性和S0偏差仍是Major；90°数据也只是已声明原生网格上的描述性读数，不是极点证书。

[独立5MeV验证](interlaced_fine_profiles/direct_profiles.json)覆盖.28–1.2GeV，每点555个能量/主波样本：

| UV代表 | 严格违例样本数 | 最大η |
|---|---:|---:|
| tip | 98 | 1.0066646 |
| mid | 57 | 1.0047289 |
| ref | 59 | 1.0044671 |

这些历史违例确认3123盘没有覆盖相应中间能量，促成后来新增459行。旧谱、相位和C均保留；当前3582盘同网格通过情况另列，不抹去旧失败。

现有证据没有把剩余Major问题证明成原论文错误。未唯一给定的范数、正则化配方、离散插值、端点及代表规则都会影响有限结果；这些选择必须明示，不能从输出曲线反推并追称独立复现。

## 3582盘几何预选交付与同网格核验

[fine_C1](fine_C1_profiles/profiles.json)及[fine_C3](fine_C3_profiles/profiles.json)给出IR P1末端8.07855°、S0零点括区s∈(.4,.5)，不声称唯一。UV三点的[原生相位](fine_native_phases/phases.json)与[同一5MeV网格](fine_grid_profiles/direct_profiles.json)使用冻结后的完整C：

| 3582盘UV代表 | 原生P1 90°（MeV） | 5MeV网格读数（MeV） | 原生min ηP1 | 原生末端δS0（°） |
|---|---:|---:|---:|---:|
| tip | 752.363 | 747.621 | .311947 | 190.145 |
| mid | 681.596 | 680.882 | .215344 | 201.450 |
| ref | 677.890 | 679.351 | .200378 | 211.310 |

fine_tip_grid、fine_mid_grid、fine_ref_grid在.28–1.2GeV同一185能量×3波网格上均sampled_passed、555项/点、0严格违例。max η=1包含阈值；该网格通过不能扩大为点间、全能量或无限自旋认证。原生与网格90°读数都完整报告，不择较好看的数；质量分散、强非弹性和S0定量问题仍是Major。

## 为什么IR与UV能约束ρ，但并不自动证明正确ρ

原散射归一化是S=1+iκf、κ=π√(1−4/s)，因此

\[
\Delta=1-|S|^2=2\kappa\operatorname{Im}f-\kappa^2|f|^2\ge0.
\]

它来自完整S矩阵在两pion子空间上的压缩，只要求弹性块收缩；没有擅加η=1。只检查Im f≥0是不充分的。

电流Gram的Schur补在|S|<1时严格等价于

\[
R\ge |\mathcal F|^2+
\frac{|\mathcal F-S\mathcal F^*|^2}{1-|S|^2}.
\]

FF相位和S的失配消耗额外正谱；F(0)=1防止任意同质缩小FF，FESR限制可用谱预算，高能FF界进一步限制解析形状。原Fig.7绿色振幅的全部C固定后，其电流扩展已有包含解析误差的负Farkas上界−9.1096171369。它精确排除这个完整IR振幅，证明UV条件确实产生了非平凡的新约束。

对历史3123盘IR固定C，[同模型IR→UV分离记录](../../evidence/SAME_MODEL_IR_UV_EXCLUSION_20260911.md)又给出两份不同条件下的证书。含FF二阶端点的加强模型负界约−.185；随后在**不增加FF端点**的原有限UV条件下独立重算，[解析上界](original_UV_fixed_new_IR_analytic/report.json)仍为−.1484359517403。后者不是前者的转移，原生行前缀及误差包络另经核对；两份对偶的负数大小也不度量条件强弱或ρ误差。它仍只排除这个新IR C，不是联合域或所有无ρ振幅的不可行证书。

同时，Gram一个零特征值不等于秩一、两pion饱和或弹性。两个谱矩只给加权平均能标，不能决定峰数和宽度。边界FF模长也不在缺少内部零点信息时唯一决定相位。完整推导、适用前提和明确反例见[有限证明](../../evidence/ANALYTIC_CARDINAL_FINITE_PROOFS_20260911.md)与[根因说明](ROOT_CAUSE_ZH.md)。因此“有ρ式上穿”是数值现象；“正确且稳定的ρ复现”还须实际完成图形比较与验证。

[fine IR在原UV下的独立排除](FINE_IR_UV_EXCLUSION.json)再次给出解析upper=−.1484093366072121；固定的是fine_IR完整C，不借新增FF端点或新增高能散射条件，不沿用3123盘固定C的负界，也不声称整个联合域不可行。

[Gram圆盘与非弹性证明](../../evidence/GRAM_DISK_AND_INELASTICITY_20260911.md)及[3123盘原数据诊断](GRAM_DISK_DIAGNOSTIC.json)说明，谱份额r=|ℱ|²/R约.6时，PSD允许η约.2；r不是η²。其定量例子属于3123盘，不能改写成3582盘新测量。它解释原条件的许可范围，没有修改谱、相位或参数，也没有消除ρ定量问题；新的同网格通过同样不是完整连续QCD证明。

[有限电流矩投影判据](../../evidence/FINITE_CURRENT_MOMENT_PROJECTION_20260911.md)补全了给定S/FF、同一正权重支持集上两条FESR的谱存在性判据；弹性边界仍需range相容条件。保存的浮点程序必须使用实际权重比W₁/W₀，不能直接替换为未经逐位核对的节点s。端点稀疏补谱是存在性构造，不是预测的ρ峰，也不建立连续谱版本的前提。

[闭式IR族及UV排除证明](../../evidence/GLOBAL_IR_SEED_AND_UV_OBSTRUCTION_20260911.md)给出同一解析族内、全能区和全部允许分波满足弹性块收缩且没有P1 90°上穿的弱耦合例子，并证明原M50有限UV输入排除其中的充分弱耦合部分。它不构造完整多通道幺正S矩阵，不替代Fig.7边界代表，也不证明所有无ρ振幅都被UV排除。

[峰存在的弱严格结论](../../evidence/RHO_SIGNAL_PEAK_EXISTENCE_20260911.md)对冻结的旧1500盘原生基线UV三代表验证了P1强度和向量FF平方的内点峰存在。所有原生局部极大候选均保留；结论不确定精确峰位、唯一性或极点，不能转给离散电流R的折线峰，也不自动转给新模型的C。旧振幅的幺正违例仍成立，故有解析内点峰不等于正确物理ρ已复现。

## 已实施的修补及其限度

1. **中心保存。** 六个既有UV终端点的活动变量逐位相同，仅发生允许的自由高能谱Schur提升；116项输入哈希复核支持这一范围。[旧三色IR身份审计](IR_CENTER_IDENTITY_REVIEW_ZH.md)则确认每份3876个C都在中心求解后被乘1−10⁻⁸并舍入，中心标签不能转移到交付点。它们的原式可行性、有限支持及解析primal证据仍对应交付C，保留有效范围。新路径以完整C身份检查区分数值中心与恢复后的可行点。另须更正旧联合中心的收敛范围：旧实现以步前的小Newton decrement标记步后点；尚未按拟修正的“保存刚通过检查的点”规则重验，不能追称旧中心已按新标准通过。其独立primal、support和完整C身份核对仍有效。
2. **同一函数上的加点。** 第一轮window变体保留原生1500盘及全部4076原变量，加入ρ窗口93盘；其联合见证τ<0，并通过1593盘及所有原电流条件的解析审计。
3. **精确端点。** 本有限族h_j(∞)=0，全能量幺正必需T0=0。FF的1/s幂次要求z=−1处双零点。五个仿射条件以依赖变量的精确代数定义实现，不能用“小浮点残差”代替。第一份仅加端点的联合见证通过了当时约束，但该中间点仍有领先高能违例。
4. **领先高能条件。** 从完整核的受控极限推得两个线性、三个凸二次条件；在1593盘加端点的中间变体中已实现并得到另一完整联合见证，五式严格正。该结论给固定有限波集合的最终高能幺正性存在，不给有限起始能量或整个中间能区。
5. **全波间隙。** 在查看新相移前固定交错网格，发现上一点668项严格违例，再实际加入1530盘。早期Phase I停滞与失败记录保留；第七轮取得的3123盘完整见证之后，IR对照及UV三点均已完成中心、支持、解析核验和冻结。这些有限证据不覆盖所有中间能量或无限自旋。
6. **新IR切片。** [M3/L2控制与解析核验](IR_slice_control_analytic/report.json)、[固定截面重新居中](IR_slice_section_recenter/report.json)、[选择](IR_slice_section_selected/selection.json)和[求值](IR_slice_section_evaluate/evaluation.json)作为早期真实小控制保留。其后M50/L10的[3123盘IR支持](interlaced_IR_ref_02/report.json)、[解析核验](interlaced_IR_ref_analytic/source_audit.json)及冻结求值已经交付；这一个新IR对照不等于全部新IR区域或新的三色Fig.5–7均已完成。
7. **更细主波网格。** 3123盘的98/57/59项旧违例及3582首次重放的6项S0违例均保留。新增459行后，3582盘见证、IR及UV三点中心/支持/解析、几何冻结和求值已完成；同一185×3的5MeV网格现全部通过、零严格违例。它修复这组已知网格问题，不证明连续区间或无限自旋。

新增约束均是明确命名的修补变体，物理质量、s0、QCD输入、χ容差/范数、矩误差、FF界和B(M)保持固定。原生基线及失败记录保留，未将新约束追认为作者程序本来已有。

## 正则化与求解能证明什么

实际双密度L4界只直接约束双密度。结合原生S0/S2虚部的自由列恒等式可控制两组单密度，再用实部控制T0，才得到有限散射可行集的紧性；不是凭L4或不可能成立的“全H列满秩”宣告紧性。原生M50/L10的H是3011×3876，不可能全列满秩。

加入电流后，自由高能R无上界，完整联合集不紧。Newton中心可在严格内域消去这些谱后作Schur提升；全谱锥通路则可保留全部R和Gram，不能把消元说成唯一算法或把闭边界range条件省略。严格可行截面上的有限μ中心与μ→0的唯一极值振幅也不同。中心的支持间隙不控制相移误差，原式浮点H的支持证书与包含解析源误差的primal证书分开归档。

[正则化目标值证明](../../evidence/REGULARIZATION_VALUE_PROOF_20260911.md)在固定其它有限输入、非空可行域上证明有限B的V(B)有限、非减且凹，并用已有对偶及见证重放得到原M50/L10下V(B)−V(B/2)≥.0054931289277。它严格说明该B规范会影响有限投影边界；没有改变生产B、扫描相移或证明B→∞有有限极限，也没有证明B单独造成质量、Fig.8形状或全部L差异。

早期[缓存及实际Newton方向审计](INTERLACED_NEWTON_DIRECTION_AUDIT.md)未发现足以解释当时停滞的缓存/所测导数失配，且没有据此宣布不可行。3123与3582盘后来都取得完整可行见证；3582盘几何主线也已交付。数值中心、有限支持和同网格通过仍不保证正确ρ或唯一物理振幅，当前正式Watson诊断另行记录。

## 后续研究的作用与不能替代的工作

作者2403/2505及正式PTEP版有公开离散实现和额外迭代选点信息，可以核对来源选择；这些版本含新的谱表示、输入或步骤，不能整套移植后宣称原2309原参数独立复现。正式版Appendix D仍使用有限节点，没有提供当前欠缺的连续幺正证明。逐版本与commit核查见[来源报告](FOLLOWUP_SOURCE_CHECK_ZH.md)。

3582盘几何主线已按当时报告交付，后续为已声明的[Watsonian诊断](../../evidence/WATSONIAN_DIAGNOSTIC_SCOPE_20260911.md)。[fine_watson_01_support_01](fine_watson_01_support_01/report.json)因linear residual约1.08257e−5停止，报告inconclusive，不是不可行证书。本次同步时，允许真实下降的不精确方向、保持原中心门槛及保存刚通过检查点的修补仍在106项测试中，尚未重启，尚无已完成的 Watson 物理结果。它不是原2309算法，不覆盖原几何三点的已测差异。区域、加强模型分辨率及ρ定量问题仍待处理；原五组F始终保留原基线身份。

## 证据入口

- [当前状态与运行名](PROGRESS.json)、[逐步骤台账](../sequential_reproduction_20260910/CLAIM_LEDGER.json)。
- [联合中心保存审计](../sequential_reproduction_20260910/CENTER_REPRESENTATIVE_ROUNDING_AUDIT.md)、[三色IR中心身份审计](IR_CENTER_IDENTITY_REVIEW_ZH.md)、[缓存审计](INTERLACED_STATE_AUDIT.md)、[密度曲率补充](DENSITY_CURVATURE_SUPPLEMENT.md)。
- [有限解析族和正则化证明](../../evidence/ANALYTIC_CARDINAL_FINITE_PROOFS_20260911.md)、[高能证明](../../evidence/ANALYTIC_TAIL_UNITARITY_20260911.md)、[阈值证明](../../evidence/ANALYTIC_THRESHOLD_UNITARITY_20260911.md)。阈值条件尚未加入生产优化。
- [正则化目标值依赖](../../evidence/REGULARIZATION_VALUE_PROOF_20260911.md)、[全能区IR族及有限UV排除](../../evidence/GLOBAL_IR_SEED_AND_UV_OBSTRUCTION_20260911.md)、[给定S/FF的有限谱充要判据](../../evidence/FINITE_CURRENT_MOMENT_PROJECTION_20260911.md)、[冻结UV振幅的内点峰存在](../../evidence/RHO_SIGNAL_PEAK_EXISTENCE_20260911.md)。各自的严格范围不互相替代。
- [新IR在原UV条件下的独立排除](../../evidence/SAME_MODEL_IR_UV_EXCLUSION_20260911.md)、[Gram圆盘与非弹性](../../evidence/GRAM_DISK_AND_INELASTICITY_20260911.md)及[原数据诊断](GRAM_DISK_DIAGNOSTIC.json)、[3123盘原生相位](interlaced_native_phases/phases.json)、[独立细网格](interlaced_fine_profiles/direct_profiles.json)。
- [原论文](https://arxiv.org/abs/2309.12402v3)、[作者公开仓库](https://github.com/hyfysics/gauge-theory-bootstrap#paper-snapshots)、[正式后续版本](https://academic.oup.com/ptep/article/2026/6/063B05/8686484)。

全部科学计算仍用 `python -m smatrix_bootstrap.run`。[此前Watson接口验证](verification_Watson_adapter_final.json)记录106测试、21核心模块、3测试文件；负R、partial及projection的[实际重放](WATSON_ADAPTER_REAL_REPLAY.json)与[代码审阅](WATSON_ADAPTER_CODE_REVIEW.md)已经归档。[预算控制更正](WATSON_BUDGET_CONTROL_CORRECTION.json)明确goal_00是**零预算的一步控制，返回原可行点且不授予新中心资格**：solver_seconds=0但时间在步后检查，实际iterations=1，中心和支持标记均False。这些验证只证明接口控制范围，不是正式Watson物理结果，也不代表正在修补的新代码已通过测试。本稿继续作为工作记录，不宣称全部复现完成。