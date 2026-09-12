# He–Kruczenski 2309.12402v3 复现状态

更新：2026-09-09。**核心claims尚未全部完成；本轮固定最有来源依据的combined-l2 + hard-midpoint，从B/C重新对齐后推进D–F。** [原Fig.5残差指纹](../../../results/runs/mainline_alignment_20260909/CHIRAL_SOURCE_RESIDUALS_ZH.md)与作者后续stack8一致；[截止规则](../../../results/runs/mainline_alignment_20260909/CUTOFF_DECISION_ZH.md)对应原文统一π/M。输入及选点在新相位前固定于[PAPER_MAINLINE](../../../results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json)。原H、pure区域及电流准备可复用；旧手征/联合可行性须新范数重新验回。

## 本轮按研究依赖推进

| 环节 | 已完成/正在执行 | 结论范围 |
|---|---|---|
| A2 / χ定义 | 原Fig.5三色合并L2/ε约1.000007/1.000020/1.000146，与作者后续代码一致 | 有明确源图指纹的执行推断；不是根据复现相位选范数 |
| A2/D2 / 权重与单位 | hard统一权重；[双谱归一化](../../../results/runs/mainline_alignment_20260909/NORMALIZATION_TRANSFER_ZH.md)确认B无需另乘π或2 | 不把clipped端点修正追认为原文规定，不改变质量/匹配能标/容差 |
| B3 / 六ε | 166条旧查询按新球径向恢复primal、按新支持函数重验dual；然后补必要几何缺口 | [重验清单](../../../results/runs/mainline_alignment_20260909/B_combined/manifest.json)；pure不受χ分组影响 |
| C1–C3 | 三色局部上侧证书及C1→C2→C3已完成；132项新物理检查通过 | 阈下线性改善及IR的P1不足得到支持；[新C结果](../../../results/runs/mainline_alignment_20260909/C_RESULT_ZH.md)，不替代B全区域 |
| D1–D3 | [D3_paper_hull](../../../results/runs/mainline_alignment_20260909/D3_paper_hull/report.json)完整4076变量通过新联合约束 | current算子不施加χ范数，D1_hard_M50保持；见[当前D](../../../results/runs/mainline_alignment_20260909/D_RESULT_ZH.md) |
| E / 支持 | combined-hard tip/ref/mid原式支持均完成，gap=9.70e−5/7.68e−4/7.71e−4 | 三点已在相位前固定；新相移正在评估，Fig8下侧正在求解 |
| F | 旧五组结果保留 | 新输入下的相位/稳定性仍须重算，不能继承旧F验收 |

数值教训：两项步长/权重试探没有闭合旧mid，已恢复既有QR与线搜索策略。可行混合不是障碍中心；使用完整中心路径并保存中心间S/F变化。旧hard-separate ref已经原式达标、P1中心间变化约.157°，但它在新合并球下不能直接当作可行解。

以下审计与阶段表记录此前separate/clipped条件模型；其来源、向量及失败继续保留。最新判定以本节和PAPER_MAINLINE为准。

## A–F本轮重新审计

| 对齐步骤 | 已执行工作 | 结论 |
|---|---|---|
| R1 → 全阶段 | [原文](../../../results/runs/claims_alignment_20260909/SOURCE_REAUDIT_ZH.md)、[定义](../../../results/runs/claims_alignment_20260909/DEFINITION_REAUDIT_ZH.md)、[实际选点/producer](../../../results/runs/claims_alignment_20260909/SELECTION_REAUDIT_ZH.md)独立重审 | 未发现应直接改χ、Gram、rho2、FF平方界或FESR单位的重大错误；未公开数值设置逐项列出 |
| R2 → C/E/F | 新Fig.5原270点；48组bootstrap/phenomenology对照，36组同无量纲原生s比较；新增统一`compare`入口 | C的IR-only对照得到支持；E非tip P1与高能S0的主要差异仍在，不能由约0.3%显示质量差消除 |
| R3 → E2/E3 | 固定原mid的x、全部输入、500y目标，继续真实checkpoint；[Newton接线复核](../../../results/runs/claims_alignment_20260909/NEWTON_REAUDIT_ZH.md) | 有界诊断完成，精度目标未达；借界dy≤3.2712e−5，P1最大变化23.0°、仍无90°上穿。[结果](../../../results/runs/claims_alignment_20260909/MID_PRECISION_ZH.md)不替代原支持完成 |
| R4 → 报告/代码 | 修正prepare的χ范数误标；同步本表、SCIENCE及路线，保留旧向量与结果 | 仍为15个核心模块、3个测试文件；120项测试通过。报告按实际claim判定，不按文件存在判定 |

尚需依次处理三个科学关口：E代表振幅身份/精度 → E三点物理机制与形态 → 对最终E流程验证分辨率。它们不是“再跑三次必然完成”的承诺；只在有明确证据时修正设置，不做贴图扫描。以下保留各阶段实际交付与限制。

## F逐步验收

| 步骤 | 已完成结果 | 核心claim与范围 |
|---|---|---|
| F1 / Fig.11 | 五组(50,8)/(50,10)/(50,12)/(45,10)/(60,10)完整UV +x支持；最大gap=9.0232e−5 | 同一预选规则、全部幅度方向、固定物理输入及B(M)离散配方；不按相移选点 |
| F2 / 原式与分辨率 | 7650幺正盘、510 Gram、10 χ球、5 L4界、20 FESR、74 FF全部通过；669项新原生主波检查通过 | 五组P1均有90°上穿，读数840.6/812.4/798.3/752.1/773.8MeV；S2较稳，低能S0较接近 |
| F2 / 保留差异 | 固定M的P1读数跨度42.3MeV，原图约12.5；固定L跨度60.2，原图约54.3 | L依赖强度未完全重现；中高能S0绝对相移仍明显偏离原图。显示插值与支持gap不是连续或相移误差证明 |
| F3 / 图谱与实现 | [逐图来源](../../../results/runs/stage_F_mainline_20260909/FIGURE_ATLAS_ZH.md)、[重放](../../../results/runs/stage_F_mainline_20260909/REPLAY_ZH.md)、[120项测试与结构](../../../results/runs/stage_F_mainline_20260909/verification.json) | 15个必要核心模块、3个测试；未采用的F投影池初始化已移除，失败及生产源码快照保留 |

本轮只新增`resolution.py`这一必要科学模块。当前Phase I保持所有散射、χ和实际密度约束，只松弛电流子系统寻找起点；原式通过后完成无松弛全空间支持。五组均已闭合，所有本轮计算已收尾。E的既有定量差异继续保留。

## 本轮E逐步复现与审计

| 步骤 | 已完成结果 | 核心claim与范围 |
|---|---|---|
| A–D基础复核 | [独立Eq.(2.7)投影、Gram、单位审计](../../../results/runs/stage_E_closure_audit_20260909/OPERATOR_AUDIT_ZH.md)；M50代表行覆盖全部3876系数 | 所审计的κ、rho2、共轭、FESR幂和FF平方界未发现重大错误；不宣称排除了所有程序错误或连续违例 |
| E1 / Fig.8 | [7份联合点/支持及Arb参考截面](../../../results/runs/stage_E_closure_20260909/E1_regions_certified/regions.json)；[图](../../../results/runs/stage_E_closure_20260909/E1_regions_certified/regions.pdf) | UV右端严格小于IR；xref上侧变化更大已证明，比值[1.026389,1.167740]。原图显示读数约5.35的强非对称未复现；UV仍为稀疏内外包络 |
| E2 / 三点与后续方法 | 原E2三向量保留；2403三主波Watsonian一次固定目标全部求完。tip/mid/ref原式gap为.00062164/.00066370/.00092079 | 各4076变量及全部原约束通过。目标来自各自旧FF/S；投影等式按作者后续代码释放，[新旧坐标](../../../results/runs/stage_E_closure_20260909/E3_W_phases/watson_diagnostics.pdf)同时报告，不冒充仍是旧边界点 |
| E3 / Fig.9–10 | [P1](../../../results/runs/stage_E_closure_20260909/E3_W_phases/fig9.pdf)、[S0/S2](../../../results/runs/stage_E_closure_20260909/E3_W_phases/fig10.pdf)、[η](../../../results/runs/stage_E_closure_20260909/E3_W_phases/inelasticities.pdf)；396项原生检查通过 | 后续方法下三点均出现90°上穿，读数807.57/728.39/697.60MeV；tip/mid约5%偏差，ref约−9.4%。三点读数跨约110MeV，原图三曲线接近程度未复现 |
| E3内层精度 | [ref同目标补核](../../../results/runs/stage_E_closure_20260909/precision_comparison.json)完成当前μ中心，gap=.00064681 | P1相移最大变化.08271°、读数变化.029MeV、最低η .84699→.84949；不足以解释该点的主要差异。原冻结三解未替换，不把这一检查当作M/L收敛 |

[来源审计](../../../results/runs/stage_E_closure_audit_20260909/SOURCE_AUDIT_ZH.md)说明：Watsonian是具名后续方法，不是2309已经公开的原始算法；它改变选中振幅，不改变Fig.8可行域。主线仍是separate-l2、B377500及D的打印/raw FESR与FF设置，没有按图调参。combined-l2只完成数学接线、新M50联合见证及IR候选核验，不转移旧证书或另开六ε/C扫描。

本轮修正了“η接近1就足以保证Watson相位”的理解：原Gram给出 |S−F/F*|²≤(1−|S|²)(1/r−1)，r=|mathcalF|²/rho_current。r很小时，相位关系仍可能偏离；详见报告。FF引导的模π lift单独作为条件诊断，主图保留nearest-S分支；不以逐点完全饱和新增原论文的硬门槛。

[112项测试及结构检查](../../../results/runs/stage_E_closure_audit_20260909/verification.json)通过：14个必要核心模块、3个测试文件。两次未改善M50的完整native conic试验已退出活跃代码，源/候选及压缩测试保留。新Watsonian默认μ=.001，同目标续算自动读取checkpoint；已有原式界足以验收时及时停止，不再强求无必要的后续中心。所有本轮计算已收尾。

## B逐步验收与结论范围

固定PV/Legendre-Q、M50/L10、3876变量、1500个原生幺正盘、自由T0、actual-ρ组合普通L4上限377500；手征为两个四维L2球。六种ε只改变容差，其余条件相同。来源未指定的范数/正则化细节继续具名记录；不靠图形调整参数。

| 步骤 | 已完成的条件交付 | 保留的科学限制 |
|---|---|---|
| B1 | **条件有限计算器已就绪**：M50非平凡点、多项独立原H支持gap<1e−4；六ε的+x均已验界 | 最窄条带停止错误已修复并通过原H支持验界；v3原图设置身份未恢复，单独保留 |
| B2 | **条件pure区域几何已完成**：38份有效支持、36个法向，Arb几何距离.00967577335<=.01 | 原图显示包的局部差异约.016单独报告，不宣称设置身份或连续域证明 |
| B3 | **六ε条件区域几何已完成**：166份有效支持；距离依次.009813/.009730/.009850/.009567/.009587/.009074，全部≤.01 | Fig.4方法复现已完成；范数／正则化明确声明，数值差异保留，不作为无限前置任务 |

B2采用完整pure区域；B3只采用原Fig4展示的x>=0闭包，所有负x解仍保留。已验证宽支持界可以参与外包；最终按整体误差验收，不把每份旧记录都强求为严格最优点。**条件原型完成、原图数值吻合、连续幺正证明是三种不同结论。**

[论文主张对应表](../../../results/runs/stage_B_delivery_20260907/B_CLAIMS_ZH.md)、[参考点判定](../../../results/runs/stage_B_delivery_20260907/REFERENCE_POINTS_ZH.md)和[全轮廓比较](../../../results/runs/stage_B_delivery_20260907/DISPLAY_HULL_COMPARISON_ZH.md)分别给出证据。ε=.002及.001均包含两种参考点，.0006/.0002排除两者；C1现已完成Fig.5阈下形态比较，较小容差更接近线性；不宣称.002唯一或最优。

## C逐步验收

| 步骤 | 已完成结果 | 科学解释 |
|---|---|---|
| C1 / Fig.5 | 同物理xref的三个预选上边界振幅，81个阈下点；[结果](../../../results/runs/stage_C_mainline_20260908/C1_RESULT_ZH.md) | ε缩小后三波偏离下降；蓝色网格未见S0零点，橙/绿零点在.30–.35/.40–.45，趋近Weinberg的.5 |
| C2 / Fig.6 | 原样保存C1绿色的[3876系数](../../../results/runs/stage_C_mainline_20260908/C2_fig6/coefficients.json)；[结果](../../../results/runs/stage_C_mainline_20260908/C2_RESULT_ZH.md) | 原图magenta复用绿色上边界点；选点先于相移，不按相移重新筛选 |
| C3 / Fig.7 | 同一振幅的三相移与η，43原生节点加阈值；132项物理检查通过；[结果](../../../results/runs/stage_C_mainline_20260908/C3_RESULT_ZH.md) | S0/S2低能形态较接近，P1缺正确ρ；保留S0中能/高能等差异，不作逐点一致或全域排除主张 |

代表点采用已验证内包边的凸组合，并再检查原1500盘、χ和密度约束；是有区域误差说明的边界近似，不冒称原作者唯一精确极值。98项测试及十文件／三测试结构检查通过。B阶段教训已写入[AGENTS](../../../AGENTS.md)：先明确物理决策、复用数据、限制计算范围和并发，不能以假设来源不完全或逐点差异无限阻断主线。

## 原E投影基线（2026-09-08，保留）

[原报告](../../../results/runs/stage_E_mainline_20260908/E_RESULT_ZH.md)保留6份联合支持、原三套E2向量及未加Watsonian时的相移。原tip/ref读数约812.36/699.96MeV，mid在既定展开中未见90°上穿。新结果没有覆盖这些证据，也不能用投影移动后的相移追认旧上边界三点的原始claim。E1的旧xref符号未决区间已由本轮下边界补算和Arb截面更新。

## D逐步验收

| 步骤 | 已完成结果 | 结论范围 |
|---|---|---|
| D1 | [同一PV/M50散射算子的电流准备](../../../results/runs/stage_D_mainline_20260908/D1_current_M50/report.json)：100个Gram、四矩、14个FF锥；独立复／实Gram接线测试 | F0(0)≈1、F1(0)=1、mathcalF=kF，电流共用原S0/P1，未强制η=1 |
| D2 | [四矩、截止与FF输入冻结](../../../results/runs/stage_D_mainline_20260908/D1_D2_RESULT_ZH.md) | 求解前选定clipped-phi、打印矩、raw误差.002及epsilonFF=6e−5；常谱诊断不作为未知谱的误差界 |
| D3 | [完整4076变量联合见证](../../../results/runs/stage_D_mainline_20260908/D3_hull_009/joint.npz)；[768-bit原式复核](../../../results/runs/stage_D_mainline_20260908/D3_original_audit_768/joint_audit.json)全部通过 | 原1500盘、两个χ球、actual-ρ L4、100×7个Gram主子式、四矩及14个FF全部检查；不宣称极值／区域完成 |

该点投影约(.06491413702,−.004112130549)，两项χ范数约.001999815307/.000627365577，L4=370688.3651<377500。[重新生成768-bit核](../../../results/runs/stage_D_mainline_20260908/D3_fresh_native_768/evaluation.json)的1500个原生分波也通过，χ保持在容差内；新核与原H收缩的最大复f差约5.10e−11，仅是此点的舍入诊断。

D3复用B同设置的204个完整可行振幅，并显式联立χ、电流及UV。390维凸包搜索产生的是原M50全4076变量的合法见证；最终按原式核验，不按求解器AlmostSolved标签验收。求解约.263秒，含输出与原式核验共5.779秒。它只回答联合可行性，不能代替E1的全振幅空间支持；所有[先前失败](../../../results/runs/stage_D_mainline_20260908/ATTEMPTS.json)保留。未奏效的QP/L2初始化分支、旧工作集派发和独立recover入口已退出活跃代码，必要数学检查与生产源保留。[全部101项测试及十文件／三测试结构检查](../../../results/runs/stage_D_mainline_20260908/verification.json)通过。

## 科学来源及已知差异

- [原文/我方方法审计](../../../results/runs/stage_B1_source_method_audit_20260907/OUR_METHOD_REVIEW_ZH.md)：原文采用有限M/L计算，未给连续全域认证；T0zero、五尾、FG和10062行是本项目加强分支，当前不再扩展。
- [PV精确核映射](../../../results/runs/stage_B1_source_method_audit_20260907/COLLOCATION_KERNEL_MAP_ZH.md)：当前节点核对应后续作者解析投影。PV与finite-sine作用范围、秩及可行性结论不混用。
- [regulator结果](../../../results/runs/stage_B_regulator_protocol_20260907/RESULTS_ZH.md)：手征预设窗口未满足我们后来增加的1%/两decade平台。失败保留，该附加认证移出B原型交付的前置条件；六ε仍共用同一个固定B。
- [Fig4采样审计](../../../results/runs/stage_B1_source_method_audit_20260907/FIG4_SAMPLING_ZH.md)：ε=.002本模型右端约.1059，原图最右显示点约.0826；ε=.0002分别约.0545和.0224。maxsample不是严格支持上界，但没有证据把全部差距归于裁切、读数或径向角。不得宣称数值已吻合。
- [共同regulator诊断](../../../results/runs/stage_B1_source_method_audit_20260907/COMMON_REGULATOR_DIAGNOSIS_ZH.md)及[尚需原始设置](../../../results/runs/stage_B1_source_method_audit_20260907/ORIGINAL_SETTINGS_REQUEST.md)：缺失项明确记录，不通过扫描范数拟合原图；未向外发送消息。

## 执行与证据

[主算子](../../../results/runs/stage_B_pv_M50_L10_prepare_20260907/report.json)为3011×3876，约7秒准备；仅定义原生物理节点、阈下和阈值极限，不冒充任意物理离节点的连续振幅。

[最终七域数据](../../../results/runs/stage_B_delivery_20260907/regions/regions.json)同时给出B2/B3_geometry_ready=true、204有效记录、0拒绝、无嵌套矛盾。[B3最终派发记录](../../../results/runs/stage_B3_PV_B377500_region_20260907/manifest.json)已收尾，无运行或待派发作业。证书适用于保存的float64有限算子，不包含离散化或连续域误差。

Newton提前停止错误已修复：小减量须同时通过梯度支持宽度检查，否则执行已算出的可行Newton末步；M50同题原H区间为[.7842651767,.7847169067]，gap=.00045173，约122秒。该修复后全套97项通过。随后修复近重复法向的几何排序，使用半平面和Arb叉积；相关2项测试及实际M50重汇总通过。算法与范围详见[最终报告](../../../results/runs/stage_B_delivery_20260907/B_RESULT_ZH.md)。已删除未改善M50的CVXOPT活跃分支；[结构检查](../../../results/runs/stage_B_delivery_20260907/structure.json)确认10核心文件（各≤350行/24KiB）、3测试文件（各≤350行）。

## 全主线

| 阶段 | 对应核心步骤 | 状态 |
|---|---|---|
| A | A1有限合同；A2手征/UV单位；A3同族算子 | 条件版本已完成；sine与PV范围分别声明 |
| B | B1有效有限支持；B2 pure全域；B3六手征正窗口 | 方法与有限计算完成；原图定量差异保留 |
| C | C1 Fig.5阈下形态；C2 Fig.6代表振幅；C3 Fig.7 IR相移 | 条件三步完成；原图曲线直接比较已补齐，选点近似保留 |
| D | D1 FF/电流Gram；D2 FESR/FF caps；D3共同M50联合解 | 三步完成；完整原式验收及新核原生幺正检查通过 |
| E | E1 Fig.8区域；E2原三点及具名后续步骤；E3 Fig.9–10/η | 核心claims部分：原mid精度需独立判定，原三点稳健性/中高能S0及强非对称未闭合；后续三点不能代替原三点 |
| F | F1 Fig.11五组M/L；F2独立物理审核；F3逐图来源报告 | 五组原UV+x计算完成；ρ保留/S2稳定得到支持，L依赖及中高能S0仍不支持原文全部稳定性claim |

原计划A–F都有计算交付，五组既定F任务已收尾；这不等于核心claims全部完成。最新[逐项审计](../../../results/runs/claims_alignment_20260909/CLAIMS_ZH.md)收紧了此前“只剩定量差异”的口径，并给出三个尚未闭合关口。现有F检验原UV+x流程，没有检验后续Watsonian三代表点的收敛。不得用额外扫描、重选相位或新增验收门槛改写差异。

## 历史证据保留

[此前完整STATUS及失败索引](../../../results/runs/stage_B_plan_audit_20260907/STATUS_BEFORE.md)保留所有旧运行入口。[源核心快照](../../../results/evidence/pre_global_sdp_core_20260907.tar.gz)、finite-sine/PV秩证明、初值、10062行状态及违例均未删除，也不默认重扫。典型已知违例：8221行候选在论文窗口主三波有45项偏离、最坏eta≈1.000090；全域最坏eta≈1.139在22.72GeV。它们属于相应sine候选，不转移给PV或原论文。
