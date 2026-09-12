# 复现状态：He–Kruczenski 2309.12402v3

更新：2026-09-07。 **进展审计后已暂停新增长任务；B阶段未验收通过，当前没有可靠完成ETA。** 唯一科学主线是 [Bootstrapping gauge theories](references/2309.12402v3.pdf) 的 Fig.3–11；Snowmass 只作背景，标量基准和有限秩证明只作校准与诊断。**论文主要区域及相移尚未复现，尚无同时通过支持最优性与独立物理审核的边界点。**

阶段 A 已完成条件定义及算子，**B1–B3仍未整体验收**。B1已在多个明确有限问题上得到有效支持区间；当前按作者公开的Legendre-Q/PV-midpoint离散化推进M50/L10原型，保留全部3876变量和1500个原生幺正约束。大regulator的手征右端约.1059，与Fig.4约.0826存在明显差异；source核换成作者collocation后差异仍在，不能靠换solver或按图调参掩盖。按原作者方法完成预定regulator序列后，五个手征+x支持区间均已闭合，但相邻增量均超过1%，未找到共同稳定平台；不能冻结一个最贴图的B。后续作者代码B377500的条件区域比较仍在推进。

## B1 方法审计后的执行口径

当前首要目标是完整支持论文claims的原型机。[共同regulator诊断](results/runs/stage_B1_source_method_audit_20260907/COMMON_REGULATOR_DIAGNOSIS_ZH.md)及[Fig.4采样审计](results/runs/stage_B1_source_method_audit_20260907/FIG4_SAMPLING_ZH.md)说明为什么不能把差异简单归为solver或显示角度。已完成[原文审计](results/runs/stage_B1_source_method_audit_20260907/ORIGINAL_PAPER_ZH.md)、[我方方法review](results/runs/stage_B1_source_method_audit_20260907/OUR_METHOD_REVIEW_ZH.md)及[离散基后续研究核对](results/runs/stage_B1_source_method_audit_20260907/DISCRETE_BASIS_ZH.md)。论文的理论幺正性要求全域成立，但其公开计算采用M50/L10有限节点，Fig.11给五组分辨率的相移稳定性；没有报告排除全部违例的连续域证明，亦不能因此认定作者振幅有违例。

原型推进与加强审核分别记录：**有限求解已经闭合过一个M50支持问题；与论文设置对应的B1原型验收仍未完成，B2/B3区域未完成。** T0zero、五尾、global FG及10062行是本项目的加强分支，不能默认为作者设置或让全域认证阻塞全部D阶段。有限原型／加强分支、作者解析分波缩放、PV节点算子与处方身份已接通；当前优先完成B，FF/current连接保持已有代码状态；仍保留全部已知违例，关键相移主张要接受相关能区的数值稳定性检查。

违例应按物理问题判断：新8221点的全域最坏eta≈1.139在22.72GeV；论文E≤1.2GeV窗口内三主波最坏eta≈1.000090，共45个已列违例。旧eta≈170属于140GeV，旧ROI主波最坏约1.017。它们既不能混写，也不能因窗口内较小就事后改判通过。所有数值及来源见方法review。

## 全流程六阶段、十八步骤；A 已完成，剩余五阶段十五步

A1–A3 已按“明确条件分支”的验收口径完成。B1–F3 共十五项仍待完成；已有诊断和代码不折算为论文区域／相移的完成比例。

| ID | 阶段／具体交付 | 原文对应 | 当前状态与验收缺口 | 依赖／实现位置 |
|---|---|---|---|---|
| A1 | 明确散射有限函数族及有界性处方 | (2.7)、(3.58–62) | **条件定义完成**：全部 sine-cardinal 密度、未减除坐标、ordinary l4 上限 B(M)=100×双密度维数。A 原始算子保留 T0；free 分支只是采样松弛，当前 `boundary` 默认 zero，是本族全能量幺正性的必要非充分条件 | `kernels/SineCardinalBasis`、`convert_subtraction`；零条件有本族推导，未认定为作者原始处方 |
| A2 | 锁定手征、电流、FESR 全部单位和容差 | (2.33)、(2.50–57)、(3.64)、(3.72–75) | **条件定义完成**：独立四矩 raw ε=.002、printed targets、hard-midpoint、平均 m_q；归一化误差／partial-cell／RMS 等具名比较分支；OPE 围道积分独立核对通过 | `model/UVConfig`、`fesr_data`；后续 D 实际施加约束 |
| A3 | 同一振幅的通用 M/L、阈下与物理求值 | (2.9–12)、§3.1 | **完成**：通用求值、精度自适应、并行能量算子和精确减除换元；最新密集算子为M50/L14、832-bit/order40 | `operators/run`；算子准备与初始化均不等于支持结果 |
| B1 | 首个非平凡支持点及有效内外界 | §3.1 的凸最优化 | **未完成验收**：PV M50/L10五个预定B的手征+x gap均<1e−4；明确排除该序列中的1%平台。B377500 pure有效区间[2.2348894,2.2362318]。原图处方差异仍未解决，B2/B3只能作具名条件比较 | A1–A3；全H可行下界与有效对偶上界共同形成间隙 |
| B2 | 纯 S 区域 | Fig.3 | 未完成；PV、B377500条件pure四初始方向已齐，形成有效内外包；最新10个不同法向、12份有效支持记录；几何距离约.11635>.01。新增派发已暂停 | B1；旁路基线，不阻塞 C/D |
| B3 | 六种手征容差的区域 | Fig.4 | 未完成；六种ε的+x支持系列正在计算；六种ε右端均已闭合，六个正半平面外包已闭合；最新几何距离约.52–.65>.01，0/6通过。新增派发已暂停，已启动5项正常收尾 | B1；`analysis` |
| C1 | 阈下 S0/S2/P1 与线性近似 | Fig.5 | 未完成；检查整个 0<s<4 的形状和 S0 零点，不能只检查四个约束点 | B3+A3；`operators/analysis` |
| C2 | 手征区域上的代表点与完整振幅 | Fig.6 | 未完成；固定靠近参考黑点的选点规则，保存系数，记录选点歧义 | B3；`analysis/io` |
| C3 | 仅手征信息的三种相移 | Fig.7 | 未完成；同一解恢复 δ、η，重现 S0/S2 较好而 P1 不足的物理对照 | C2+A3；`analysis` |
| D1 | 形状因子及两个电流 Gram 算子 | (2.33–41)、(3.65–71) | 部分基础：独立电流历史 witness；已有D1 affine算子及独立代数控制；100个3阶实PSD／4矩／14FF cap已接线，缺M50生产准备和联合求解 | A2；可与散射线并行，`model/operators` |
| D2 | 四项 FESR 与高能 FF 上界 | (2.56)、(3.72–75) | A2 已提供通用四矩、单位／权重、容差及 FF caps；仍需与 D1 的当前电流变量连接并在求解中施加 | A2+D1；`operators/linear` |
| D3 | 当前有界散射与电流联合可行点 | (3.70) | 未完成；共享同一 S0/P1，逐项检查散射、手征、PSD、四矩和 FF | A3+D2及有限散射初值；不等待B的全域认证，`model/analysis` |
| E1 | 联合 QCD 区域及上下边界变化 | Fig.8 | 未完成；先验证联合支持上下界，再扫青色区域，并与同设定手征区域比较 | B3+D3；`analysis` |
| E2 | 尖端及两个邻近 QCD 代表点 | Fig.8 的红／粉点 | 未完成；保存完整散射及电流解，使用不依赖实验曲线的选点规则 | E1；`analysis/io` |
| E3 | QCD 的三种相移和 ρ 共振 | Fig.9–10 | 未完成；连续选支，比较三点、90° 交点及论文约 6% 偏移；不以实验调未定参数 | E2+A3；`analysis` |
| F1 | 五组分辨率复算 | Fig.11／附录 A | 未完成；(M,L)=(50,8),(50,10),(50,12),(45,10),(60,10)，固定选点规则 | E3；`run/analysis` |
| F2 | 独立可行性、支持及未采样约束审计 | 加强论文的数值检验 | 未完成；离节点／高自旋检查、有效支持残差与误差预算；采样检查不得称连续认证 | B–F1；`linear/analysis` |
| F3 | 可重跑的逐图比较与最终来源报告 | Fig.3–11 全链 | 未完成；数据、系数、配置、源哈希、误差及处方差异齐全；缺作者设定时只可报告条件版本 | 所有交付；`io/run` |

关键路径仍是 **A → B1/B3 → C 与 D3 → E → F**。当前主算子为[PV M50/L10](results/runs/stage_B_pv_M50_L10_prepare_20260907/report.json)：约7秒完成3011×3876矩阵准备，角投影采用解析Legendre-Q；原生物理节点及阈下求值有明确定义，未定义的物理离节点求值被明确拒绝。仅允许sampled/free，不套用sine的FG／五尾或解析重采样证明。[精确核映射及差式](results/runs/stage_B1_source_method_audit_20260907/COLLOCATION_KERNEL_MAP_ZH.md)解释两种处方的区别；所有历史证明和失败继续分开保存。

[Regulator协议](results/runs/stage_B_regulator_protocol_20260907/protocol.json)先固定对数序列100、1000、10000、100000、377500，[五点结果](results/runs/stage_B_regulator_protocol_20260907/RESULTS_ZH.md)已用有效区间排除预定序列中的1%平台；不能外推序列以外，也不得根据PDF距离选上限。[源方法审查](results/runs/stage_B1_source_method_audit_20260907/REGULATOR_METHOD_ZH.md)还说明怎样安全复用同题dual。此前增强的10062行长任务已停止并保存stop_reason／已有状态；不会将未收敛结果计作完成。

数值求解使用相同物理问题的barrier和Clarabel路线。监督器已固定BLAS／OMP线程预算为8；原生锥采用有依据的正缩放、可行参考仿射平移及辅助变量自然单位。它们均保持原约束和全部幅度变量，效果由实际M50结果判断。

[B3条件支持及图形](results/runs/stage_B3_PV_B377500_region_20260907/manifest.json)保留全部负x系数，按原Fig.4的x≥0窗口截取有效内外凸包；.01门槛未改。

[B2条件扫描批次](results/runs/stage_B2_PV_B377500_region_20260907/manifest.json)固定后续作者代码所对应的B377500作实际区域比较，不代表共同regulator已获接受。B3同时计算全部六种ε的右端，不能把这些单点称为六个区域。平台的1%、支持请求1e−4及区域Hausdorff .01都是本项目预先定义的数值门槛，不是2309报告的定理或原始参数。

## 已实现的增强与论文交付分开计数

| 类型 | 已有实现／结构 | 对应主线与限制 |
|---|---|---|
| 方法来源整理 | 原文22项证据、官方后续18项、我方13项差异及独立离散基审计 | 约束来源、数值等价、加强条件和未恢复处方分别标注；不计区域完成 |
| 基础物理一致性 | `model`统一κ、S、幺正裕量及Gram；`kernels/operators`同一解析振幅支持物理／阈下求值 | A1–A3的条件版本；全部M50幅度3876个系数保留 |
| 数值等价优化 | `quotient/linear`精确换元、SOC密度消元、QR／迭代精化；`run`按能量并行生成 | 不改变所声明可行集；性能以实际完整运行判断 |
| 有限支持与恢复 | `analysis/imaginary`全H复核、真实系数凸恢复、工作集续跑及独立对偶上界 | B1已有有限支持闭合证据，不能自动替代物理审核 |
| 加强物理诊断 | 五尾、三个global FG多项式；M50六个25阶Gram；10062行扩展问题 | 本项目加强分支，不倒写成2309公开算法；完整连续认证尚无 |
| 区域后处理 | `io`保存完整系数的目标包络、内外多边形误差、下一法向、共同设置下的嵌套检查 | B2/B3工具已实现；实际pure及六ε区域仍未交付 |
| 电流／UV连接 | `operators.current_operators`与`prepare-current`：4076变量、100个3阶PSD、四矩、14个FF cap | D1/D2算子及独立代数控制，尚无M50生产准备或联合解 |

代码保持十核心模块、三测试文件；有限／加强分支和原生接口更新后，完整97项测试通过（正半平面几何更新后14.48秒）。CVXOPT显式primalstart及9组独立数学控制已接通；M3原H间隙约2.36e−5。现成ldl2使用密集约化KKT，未称为稀疏后端；M50对照在1200秒总预算内未完成，未得到可接受的新支持点；没有声称提速。数学测试和增强数量不折算成论文完成百分比。

## 当前结果与运行

| 证据 | 已知事实及限制 |
|---|---|
| [PV手征M50/L10有限支持](results/runs/stage_B_pv_chiral_M50_barrier_20260907/report.json) | B377500、ε.002；lower=.105885098499、upper=.105940175784，gap≈5.51e−5；原1500盘通过。它是明确有限程序的支持结果，仍不等于Fig.4已复现或连续物理解 |
| [PV pure M50/L10](results/runs/stage_B_pv_pure_M50_barrier_20260907/report.json) | 同B：lower=2.234889398548、upper=2.236231782347，gap≈.00134238；有效粗区间，未达该次请求的1e−4。其右端与Fig.3标记约2.2329接近，但整个区域仍未计算 |
| [原8221行手征支持](results/runs/stage_B_chiral_round3_recovered_support_M50_20260907/report.json)／[独立审核](results/runs/stage_B_chiral_round3_recovered_support_audit_frozen_20260907/report.json) | lower=.090157025369、upper=.090249617483、gap≈9.25921e−5；原H、FG、5尾通过。832-bit/order40的52962对中1841违例，最坏s≈26336.1861、I=2、ell=2，margin≈−.297128；故仍未通过B1物理验收 |
| [10062行恢复初值及算子](results/runs/stage_B_chiral_expanded_recover_M50_20260907/report.json) | 1841新物理行已加入，原2209个附加行保留前缀；恢复目标(.080624150415,−.005281374935)，原H/FG/5尾通过。约150秒完成并行追加算子及恢复；仅为续跑初值 |
| [M50有效初值v2](results/runs/stage_B_FG_repaired_M50_v2_20260907/report.json)／[独立审计](results/runs/stage_B_FG_repaired_M50_audit_20260907/report.json) | 目标(.080550644094,−.005276824472)；实际832-bit/order40、707×42=29694检查全过；三个p_I全局下界约.4740803、.1563960、.1616019。仅为初值及必要条件验证，非连续幺正证明 |
| [M50 pure恢复点](results/runs/stage_B_recover_pure_M50_raw_20260907/report.json)／[独立审核](results/runs/stage_B_recover_pure_M50_raw_audit_20260907/report.json) | 目标(.831310061640,.012168707442)；832-bit/order40、1261×42=52962检查全过，global FG与5尾通过，实际ρ L4≈625.912。保留[.825886初值及交接](results/runs/stage_B_strong_pure_M50_round2_20260907/handoff.json)；均不计边界或连续认证 |
| [重新验证的pure松弛对偶](results/runs/stage_B_reused_relaxation_dual_M50_20260907/report.json) | 当前完整8221行上的区间[.831310061640,2.559916601968]。旧无效primal未转移；旧6012个dual权重补零后重新按当前原H/5尾/FG条件计算。仅为有效粗区间 |
| [M50手征首轮CG原式支持界](results/runs/stage_B_barrier_CG_M50_20260907/report.json) | lower=.080550643289、upper=.092990467638；完整8221行、5尾及global FG通过，gap≈.0124398，仍非已验收边界。保存2280个物理索引供续跑；历史LP上界21478已被改进 |
| [M50手征第二轮CG](results/runs/stage_B_barrier_CG_M50_round2_20260907/report.json)／[第三轮运行](results/runs/stage_B_barrier_CG_M50_round3_20260907/worker.json) | 第二轮完整支持区间[.080550643289,.090916316050]；第三轮所得有限支持及密集审核结果见上方；不计B1完成 |
| [M50 pure第二轮CG](results/runs/stage_B_barrier_CG_pure_M50_round2_20260907/report.json) | 完整支持区间[.831310053327,2.218540783076]；工作候选2.201830仍有1509个未加入的幺正违例，不能作下界。保存3090个物理索引及有效最佳对偶继续计算 |
| [M3/L2直接裕量CG](results/runs/stage_B_barrier_CG_M3_control_20260907/report.json) | 6轮，物理工作集18→26／完整60行，FG切点21→29；lower=.04298711321828、upper=.04306848660773、gap≈8.14e−5。只验证方法，不替代M50 |
| [原生坐标／密度锥控制](results/runs/stage_B_equilibration_controls_20260907/report.json) | 放宽equilibration未改善，已恢复默认；同历史输入subtracted+SOC再得约1.60e−7的M3 gap。M50原生尝试仍未闭合 |
| 已结束的完整M50：[手征](results/runs/stage_B_native_subtracted_M50_round1_20260907/report.json)、[pure](results/runs/stage_B_native_subtracted_pure_M50_round1_20260907/report.json) | 原生分别MaxTime（13／12迭代），并未收敛。恢复后手征lower=.04174070348279、upper=289961049.6036；pure lower=.04779318187414、upper=2525930.456843。报告保留有限可行点，均未通过支持最优性或物理接受 |
| M3减除坐标控制：[qdldl](results/runs/stage_B_native_subtracted_qdldl_M3_20260907/report.json)、[faer](results/runs/stage_B_native_subtracted_faer_M3_20260907/report.json) | 两者Solved，原坐标Arb支持间隙约1.1155e−7、1.6043e−7。只验证数值方法，不转移到M50 |
| [完整密集算子](results/runs/stage_B_dense_M50_prepare_20260907/report.json) | M50/L14、832-bit/order40；767×42=32214个散射行，矩阵64439×3876，约385秒生成；这是准备产物，不是支持计算或通过记录 |

`--solver clarabel`把FG多项式在[−1,1]上的非负性精确提升为PSD：M50有六个25阶Gram矩阵。方法依据本地[Roh–Vandenberghe 2006，§5.1](references/roh_vandenberghe_2006.pdf)的标准加权平方和表示。Gram变量没有删掉密度方向，也没有加入逐分量谱密度正性。独立外界验证moment Toeplitz矩阵及原C_flat残差中的+νD/B，不依赖solver的成功标志。

`--native-coordinates subtracted`只改变数值坐标：保留自由T0_sub，用一条等式施加原T0_unsub=0；绝不另设T0_sub=0。双密度、B、κ、χ容差及原始物理公式不变。输出恢复原坐标，保留较好的已验证incumbent；真实系数收缩或凸混合的权重、目标损失和最终原式检查均记录，不剪裁S或p。

`absorptive`还可用原生S0/S2的Im f作可逆坐标，全部ρ保留；该换元的数学一致性已经检查，数值效果未解决M50障碍。Arb零中心区间的平方改为乘法，修复了理想参考点被NaN误拒绝的问题。推导见[SCIENCE](SCIENCE.md#equivalent-coordinates-and-constraint-generation-b1)。

## 判据、失败与证据保留

完整S†S=1投影到两pion弹性块给|S|≤1；eta=1只表示该通道无其它跃迁，并未强制。有限采样、固定分波的高能必要条件、固定能量的大自旋必要多项式须分开验收；即使p_I全区间认证为正，也不等于完整连续幺正性。公式与理想mode1全局散射witness证明见[SCIENCE](SCIENCE.md)，该witness不替代支持边界或gauge结果。

| 保留的失败 | 科学结论 |
|---|---|
| 原生未减除完整M50：[手征](results/runs/stage_B_native_M50_round1_20260907/report.json)、[pure](results/runs/stage_B_native_pure_M50_round1_20260907/report.json) | 第1次迭代NumericalError、返回零向量；零解不计进展，已有非平凡初值保留 |
| 旧M50径向点：[pure](results/runs/stage_B_radial_pure_M50_20260906/report.json)、[chiral](results/runs/stage_B_radial_chiral_M50_20260906/report.json) | 高spin多项式有严格负证据，含E≤1.2 GeV，不接受为物理解 |
| [原50点支持](results/runs/stage_B_equilibrated_M50_20260906/report.json)／[物理探测](results/runs/stage_B_physical_probe_20260906/evaluation.json) | 小采样支持间隙仍伴随eta≈170，不能认作物理边界 |
| [起点精度失败](results/runs/stage_B_seed_representation_20260906/progress.jsonl) | 384-bit阈值投影消减已修复，失败保留；提精度不自动给出积分包络 |
| 旧PV失败：[原始](results/runs/pv_support_377500_x_20260906/report.json)、[缩放](results/runs/pv_support_377500_x_scaled_20260906/report.json)、[中止](results/runs/pv_support_377500_x_square_20260906/report.json) | 与一致sine族分开；数值停止不是不可行证明 |

## 主线结构与历史入口

- 活跃计算为`prepare/prepare-current/boundary/evaluate/dual/embed`，`regions`仅汇总合格支持结果，`recover`仅作原系数可行恢复；rank、quotient、integer、control及PV执行已退休。旧核心和测试完整保留于[pre-global-SDP快照](results/evidence/pre_global_sdp_core_20260907.tar.gz)，证明输入与失败未删除。
- `quotient`现负责精确坐标换元、全局多项式Gram／moment锥及normal dual，`imaginary`负责高spin多项式及全局下界；其余十文件职责见[README](README.md)。当前全套97项测试通过，含精确换元、稀疏Krylov、约束续跑和独立原式支持检查；它们只校验实现，不计为科学完成量。
- 历史[finite-sine秩](results/analytic_full_source_rank3011_certificate.json)、[整数审计](results/remaining_rank_exact_integer_audit.json)、[PV完整秩](results/runs/pv_remaining_rank_20260906/consequence.json)、[小控制](results/pv_midpoint_small_rank_controls.json)和[M50重构](results/runs/core_control_20260906/report.json)保留各自有限处方的结论，不转移到作者实现或连续问题。
- [标量基准](results/scalar_benchmark_necessary_N8_M400_L40.json)、[旧PV起点](results/runs/pv_regularized_seed_20260906/report.json)、[独立电流witness](results/source_M50_exact_current_fiber_256.json)及[M10初始化序列](results/runs/stage_B_FG_initializer_M10_sequence_20260906/report.json)只作其声明范围的基础。整理的[执行记录](results/runs/repo_reorganization_20260906/report.json)、[清理记录](results/runs/repo_reorganization_20260906/cleanup.json)、[恢复索引](results/runs/repo_reorganization_20260906/recovery.json)保留。

按既定验收仍剩**五阶段、十五步骤**。先完成当前M50支持点与有效间隙，再交付区域、代表点和相移，最后做分辨率及覆盖审计。原处方未恢复时，最终结果仍只能称条件复现。

## 25小时节点的进展审计

B1–B3仍为0/3正式验收通过。B1已有完整M50有限支持证据，未证实理论不可行；原2309密度/手征范数等设置差异与共同平台仍未解决。B2条件pure的最新距离为.1163525233（10个不同法向），见[pure审计](results/runs/stage_B2_PV_B377500_region_20260907/progress_audit.json)。B3六个条件窗口的距离为.520132/.530110/.555605/.579364/.647629/.617327，见[B3审计](results/runs/stage_B3_PV_B377500_region_20260907/PROGRESS_AUDIT_ZH.md)。这些比值不能转换为完成百分比或剩余小时数。

停止新增批量扫描的原因是：当前源设置尚不能支持所要求的论文claims，部分数值路线也未形成可靠的支持间隙和几何收敛成本；继续按时间追加任务不能自动解决它们。原有失败和已启动任务正常保留/收尾。已复核并修复Newton中心停止顺序及前向/dual缩放dtype一致性，但最窄条带的实际宽dual尚未解决；新增barrier减除坐标适配仅完成数学review和相关实现检查，未启动M50生产对照，不计已证实的提速或科学结果。
