# s-matrix-bootstrap 求解器与模型决策史复盘（2026-09-05 → 09-12）

审计范围：只读 `results/runs/*` 入口文档、`results/evidence/*.md`、`SCIENCE.md`、`STATUS.md`、`AGENTS.md`、`src/smatrix_bootstrap/{conic,linear}.py` 及约 20 个 run 的 `report.json` 头部。未运行计算，未改仓库。下文路径均相对 `/home/shiqiu/s-matrix-bootstrap/`。

## 核心判断（先结论）

1. **仓库最终走上"自研 barrier+Newton+QR/CG+Phase I+支撑函数证书"路线，不是因为论文问题本身需要它，而是 09-06/07 两天里三件事叠加**：(a) 先把论文的有限原型"加强"成 finite-sine+T0zero+五尾+global-FG+8221/10062 行的更难问题；(b) 把这个病态缩放（行尺度跨 74 个量级）的问题交给 cvxpy→Clarabel，并用 300–1200 s 的时限判它"失败"；(c) 自研的 dense-QR Newton 在同一台机器上先跑出了结果，于是被 `SOLVER_DECISION_ZH.md` 以"已有可验证产出"锁定。
2. **现成求解器从未在公平条件下被试过**：`SOLVER_DECISION_ZH.md` 自认没有"同一冻结模型、同一目标、同一精度、同一硬件"的对照；我核实的 `D3_joint_001` 显示 Clarabel 在被 300 s 时限杀掉时（26 步）primal 残差已 2.5e−7、正在收敛。09-11 的全锥重试给 Clarabel/SCS 各 600 s、47M 非零的稀疏表示，同样"无首个迭代"即停。失败形式全部是 MaxTime/时限或第 1 步 NumericalError，指向表示与预算，不指向不可解。
3. **中途大量分支是自造的**：T0zero、五尾、global-FG、8221→10062 行、3123→3582 盘、FF 双零点、五条领先高能条件、B(M) 配方、hard-midpoint、障碍权重 750、μ→1e−10 的"数值中心"代表点规则、Arb 768-bit 证书——论文都没要求。真正来自论文/作者代码的只有 PV 核、1500 盘、四矩、Gram、B=100Mρ 与 combined-8 χ；而 combined-8 χ 又在 09-10 被以"输出辅助"为由退回 separate-L2。
4. 结论：**论文的有限 SDP（4076 变量、1500 个 SOC(3)、100 个 PSD(3)）是作者用 MATLAB CVX+MOSEK 直解的小问题；仓库把它做成了 450 个 run、17 GB、23 模块的自研内点法工程。这是"没苦硬吃"为主，叠加了一部分真实的来源缺口（χ 范数/正则化/FESR 打包 2309 未公开）。**

---

## 1. 时间线表

| 日期 | run 名 | 当时模型 | 当时求解器 | 结果 | 放弃/转向的自述理由 | 我的评价 |
|---|---|---|---|---|---|---|
| 09-06 | `stage_A_M50_L10_20260906` | finite-sine-cardinal，M50/L10，3876 变量，T0 保留但 `boundary` 默认 zero | — | 算子准备 37 s | — | 函数族是自选的（`OUR_METHOD_REVIEW_ZH.md` 表 2 行："未恢复假设／我们的明确函数族"） |
| 09-06 | `stage_B_first_20260906` | 同上 + 两个 4 维 L2 χ 球 + L4 B=377500 | cvxpy→Clarabel，`solver_seconds=300` | `user_limit`，320 s 无解 | report 未写理由；`scattering_scales` 从 1 到 2.4e−74 | 74 个量级的行尺度直接喂给双精度 IPM，300 s 判死：**问题表示病态 + 预算过短**，不是不可解 |
| 09-06 | `stage_B_barrier_small`、`stage_B_equilibrated_M50` | 同上 | 自研"analytic density barrier with full dense QR" | 小例 0.05 s；M50 56 步 173 s，outer.upper=.1059 | — | 第一次"自研能出数"。但该点后来在 `stage_B_physical_probe` 被查出 140 GeV 处 η≈170（`STATUS_BEFORE.md` 失败表） |
| 09-06 | `stage_B_qr_power`、`stage_B_native_M10` | 同上 | Clarabel（cvxpy） | 188 s "Solver CLARABEL failed"；M10 inconclusive | — | 同上：缩放未处理 |
| 09-06/07 | `stage_B_FG_initializer_M10_round1..6/near1..3`、`stage_B_tail_*`、`stage_B_unitarity*` | **加强分支**：T0zero + 五尾 + global-FG（6 个 25 阶 Gram PSD）+ 扩展采样 | barrier / Clarabel | 反复初始化 | "本族全能量幺正性的必要条件" | 把论文没要求的连续幺正必要条件前置为主线，自己把问题做难（09-07 审计已承认） |
| 09-07 | `stage_B_native_M50_round1`、`stage_B_native_subtracted_M50_round1` | 加强 sine：8221 盘、50106 行、**1.06e8 非零**、1950 个 Gram 变量 | Clarabel（qdldl/faer），1000–1100 s | 第 1 步 NumericalError；subtracted 版 13 步 MaxTime（1256 s，primal res .033） | `STATUS_BEFORE.md`："原生分别 MaxTime，并未收敛" | 这是**加强后**问题的失败，不能拿来评价论文问题（`FOLLOWUP_METHODS_ZH.md` 末段自己也这么说） |
| 09-07 | `stage_B_barrier_CG_M50_round1..3`、`chiral_round3_recovered_support_*` | 加强 sine，8221→10062 行，working-set | 自研 barrier+CG | gap 9.26e−5 闭合，但 832-bit 审核 52962 对中 1841 违例，最坏 Δ=−.297（22.7 GeV） | — | 支撑闭合≠物理通过；说明自研路线"出数"但出的是坏点 |
| 09-07 12:20 | `stage_B1_source_method_audit/OUR_METHOD_REVIEW_ZH.md`、`FOLLOWUP_METHODS_ZH.md` | 审计 | — | 结论：B1 是"我们选定的 sine 族＋后续代码正则化＋我们加强的必要条件"，非论文原型 | **明确建议**"使用成熟凸求解器…不必自行重写全部锥求解算法"、"先核对到可运行的成熟 conic 基线，再讨论求解器强化"、"不能把自研 Newton 的耗时当成作者原型机固有复杂度" | 正确诊断；但随后并未执行 |
| 09-07 13:13 | `stage_B_prototype_settings/SCAN_PLAN_ZH.md`、`stage_B_pv_M50_L10_prepare` | **切换到 PV-midpoint/Legendre-Q**（2403 代码核），sampled/free，1500 盘，B=377500 | barrier start_mu 1e−5，solver 600 s；"对照后端 Clarabel/faer" | PV 手征 M50 gap 5.5e−5 | 采用作者后续代码的核 | 正确方向；但把时限 600 s 写进计划 |
| 09-07 | `stage_B_regulator_protocol` | PV，B∈{100,1e3,1e4,1e5,377500} | barrier | 相邻增量均 >1%，无"两 decade 平台"；ε=.002 右端 .1059 vs 论文 .0826 | "不能冻结一个最贴图的 B" | 1%/两 decade 平台是自设门槛（`AGENTS.md` 后来也写"不要把发明的 regulator 平台变成前提"）。**28% 的区域差异从此一直未解释** |
| 09-07 | `stage_B_pv_cvxopt_M3_control`、`stage_B_pv_chiral_M50_reg10000_cvxopt` | PV | CVXOPT（dense ldl2，primalstart） | M3 gap 2.4e−5；M50 1200 s 内未完成 | `B_RESULT_ZH.md`："已删除未改善 M50 的 CVXOPT 活跃分支和依赖" | CVXOPT 本就慢，1200 s 也是预算判死；删除分支使后续无法对照 |
| 09-07 晚 | `stage_B_delivery_20260907` | PV 1500 盘，separate-L2 | 自研 barrier Newton，204 份支撑 | 7 区域几何 ≤.01 通过；pure 接近论文，χ 区域差 28% | "已修复 Newton 近中心提前停止"；"停止新的 solver、范数、regulator…扩展" | 交付的是"本方几何验收"，论文 Fig.4 未复现。此处冻结了自研路线 |
| 09-08 | `stage_D_mainline/ATTEMPTS.json`、`D_RESULT_ZH.md` | PV + 电流：4076 变量、100 个 Gram(3)、4 矩、14 FF SOC，clipped-phi | Clarabel（faer）全空间，`solver_seconds=300` | 001/002：MaxTime（26/25 步；001 primal res **2.5e−7**、dual .019）；003 第 1 步 NumericalError；Phase-I 004–007 NumericalError/MaxTime；**hull_009**：390 维凸包，Clarabel 23 步 0.26 s AlmostSolved→Arb 验收 | "最终成功来自收束为 D3 存在性问题" | **001 明显在收敛，被 300 s 杀掉**。正确动作是把时限改成 1 小时或改 dense KKT，而不是退到 390 维凸包 |
| 09-08/09 | `stage_E_mainline/PLAN_E_ZH.md`、`stage_E_closure_20260909` | 同上，E1 全变量支撑、Watsonian 具名目标 | 自研 centered Newton | mid 在 μ=1e−6 停 53 步不收敛（`SOLVER_DIAGNOSIS_ZH.md`）；"未奏效的完整 Clarabel 编译…已从活跃源码移除" | Gram 块占曲率 99.6%，"慢中心化，非 bug" | 自研 IPM 的中心化在 PSD 块上遇到经典困难；成熟 IPM 的预测-校正/自对偶嵌入正是为此设计 |
| 09-09 | `stage_E_closure_audit/PLAN_ZH.md` 末段 | 同上 | 两次原生 Clarabel 锥模型（unsubtracted/auto；faer/absorptive） | 第 1 步 NumericalError；15 步大残差 MaxTime | "该试验后端已从活跃核心直接删除" | 仍是预算内判死，随后删代码 |
| 09-09 | `mainline_alignment/SOLVER_DECISION_ZH.md`、`B_SOLVER_FIX_ZH.md` | — | — | 结论"保留 Newton"；发现复用 QR 时 CG 真实残差高达 .04 仍被当方向使用 | "Newton 在本仓库声明问题上已有原式可行点／支撑界证据" | 自认无公平对照；线性解质量门槛缺失是自研实现的典型风险 |
| 09-09 | `mainline_alignment/CHIRAL_SOURCE_RESIDUALS_ZH.md`、`CUTOFF_DECISION_ZH.md`、`PAPER_MAINLINE.json` | 切 **combined-L2 + hard-midpoint**；B/C/E/F 重算 | Newton | E 三点 789/700/696 MeV；F 四组，L12 Phase-I 失败 | Fig.5 残差指纹 1.000007/1.000020/1.000146 + 作者代码 | 这是最有来源依据的 χ 读法 |
| 09-10 | `foundations_review/REVIEW_ZH.md`、`FOUNDATIONS_PROOFS` | 证明 PV 混合配点非单一解析函数（Res f_R(s_i)=−w_i Im f_PV(s_i)）；把 combined-L2、Fig.8 标记选点标为"输出辅助" | — | `select-gauge` 拒绝标记；χ 范数须显式 | "不能改写成自始独立的推导" | 留数障碍是真发现；但把 combined-L2 判为"输出辅助"并退回 separate-L2 是教条（见 §4） |
| 09-10 | `sequential_reproduction_20260910` | **analytic-cardinal** 族（sine 换名回归），1500 盘，separate-L2，hard，B(M) | Newton + Phase I | 全链 C/E/F；E5 离节点 η 达 1.14；T0≠0 ⇒ 无穷远障碍；`HELD_INTER_NODE_ISSUE` | 用户要求反例单列不改基线 | 自选族的必要条件再次成为议题 |
| 09-11 | `physical_consistency/interlaced_*`、`fine_*` | analytic-cardinal + T0=0 + FF 双零点 + 5 领先高能条件 + ρ 窗交错行 = **3123 盘**；再 +459 行 = **3582 盘** | Phase I ×7 轮；两次全锥（Clarabel faer 47M nnz；SCS indirect）各 600 s **无首个迭代** | 3582 盘三点 752/682/678 MeV，min η .20–.31 | `CONIC_DESIGN_ZH.md`："多轮完整 Newton Phase I 尚未得到联合见证…准备用标准锥算法区分求解路径问题与约束相容性" | 唯一一次想用现成求解器做"判别"，仍给 600 s；稀疏 LDL 遇到密集块必然如此 |
| 09-11 | `physical_consistency/PV_*` | **退回 PV 主线**：1500 盘，separate-L2，hard，B=377500，有界坐标 | Newton（bounded coordinates，`INEXACT_NEWTON`、`LOCAL_BARRIER_DESCENT_BOUND`） | tip/mid/ref 803/701/697 MeV；PV 五组 F；Fig.8 比 [1.003,1.104] | "作者原文 PV 有限程序的同输入对照优先" | 数值修补占大头 |
| 09-12 | `major_claims_20260912` | 同上，支撑收紧 100×（tip gap 7.8e−7），μ 到 9.5e−11 | Newton + Arb 局部自协调下降证明 | 证书 Re S0>.881、Im S1<−.5055 ⇒ 与论图符号不符；L 跨度 46 MeV vs 12.5 | "当前明确原型的近极值区域与原文不兼容；须恢复作者数值合同" | 证明的是**自建模型**排除论文曲线，不是论文难 |

## 2. 现成求解器的真实战绩

| 求解器 | 用在什么问题 | 规模/表示 | 失败形式 | 是否做过缩放/预处理/换表示/调容差 | 出处 |
|---|---|---|---|---|---|
| Clarabel（经 cvxpy，09-06） | finite-sine M50 χ 支撑，3876 变量、1500 SOC | 行尺度 1…2.4e−74 | 300 s `user_limit`；188 s "failed" | 仅"正缩放"，无作者 Λ_ℓ 阈值重缩放 | `stage_B_first_20260906/report.json`、`stage_B_qr_power` |
| Clarabel 原生（09-07） | **加强 sine** 8221 盘 + 6 个 25 阶 SOS-Gram | 50106 行、1.06e8 nnz | NumericalError@1；MaxTime@13（1256 s） | 试过 subtracted 坐标、faer/qdldl；M3 控制 Solved | `stage_B_native_*_M50_round1_20260907/report.json` |
| CVXOPT（09-07） | PV M3/M50 χ 支撑 | dense ldl2 | M3 gap 2.4e−5 通过；M50 1200 s 未完 | primalstart | `STATUS_BEFORE.md`、`B_RESULT_ZH.md` |
| Clarabel 原生（09-08 D3） | PV 联合 4076 变量、86 个 PSD(3)、3775 密度 SOC | faer，unsubtracted，300 s | MaxTime@26（primal 2.5e−7）、@25；NumericalError@1 | "有依据的正缩放、可行参考仿射平移、辅助变量自然单位"；无时限放宽 | `stage_D_mainline_20260908/D3_joint_00{1,2,3}/report.json` |
| Clarabel 原生（09-08 D3_hull_009） | 同上但限制在 204 父点凸包，390 维 | — | **AlmostSolved 0.26 s**，Arb 原式验收通过 | — | `D_RESULT_ZH.md` |
| Clarabel 原生（09-09 E） | E ref 目标全空间 | auto/unsub；faer/absorptive | NumericalError@1；15 步大残差 MaxTime | 换坐标一次 | `stage_E_closure_audit/PLAN_ZH.md` |
| Clarabel 原生（09-11） | 3123 盘 analytic 联合可行性（零目标） | 28184×7832，**47,087,255 nnz**，faer，600 s | 674 s 无首个迭代 | 密度列缩放 B/nd^¼；无 Λ_ℓ 重缩放 | `physical_consistency/interlaced_full_conic_01/{conic_model,bounded_stop}.json` |
| SCS indirect（09-11） | 同上保留全谱 | 28268×7846，47.2M nnz，eps 1e−8，warm start | 755 s 无迭代返回，stdout 0 字节 | 暖启动 | `interlaced_full_conic_scs_01/bounded_stop.json` |
| highspy | 只做固定法向的对偶 LP（`highs-primal/ds`） | — | 正常 | — | `stage_B_native_M50_round1/report.json`（normal_dual） |
| MOSEK | **从未安装/尝试** | — | — | — | `SOLVER_DECISION_ZH.md` §2 |

核实：`SOLVER_DECISION_ZH.md` 第 2 节明言"本次所读记录没有同一个冻结模型…下 MOSEK 与本 Newton 实现的公平性能对照"。我进一步发现：**同一冻结问题上 Newton vs Clarabel 的对照也没有**——每次 Clarabel 跑的问题都比 Newton 跑的问题更大（加强分支、3123 盘）或时限更短（300–600 s vs Newton 的多轮 900 s 续算）。

**判断：现成求解器的失败是表示与预算造成的，不是问题不可解。理由三条：**

1. **密集块喂给稀疏 LDL。** 散射算子 H 是密集 3011×3876（`STATUS_BEFORE.md`：约 7 s 准备）；锥嵌入后 A 有 47M 非零，其中散射块 12492×7832 几乎全密。Clarabel/QDLDL/faer 的稀疏 KKT 对 m+n≈36k 的含密集块矩阵会产生 GB 级填充，600 s 连符号分析都做不完。仓库自研 Newton 用"full dense QR"（`stage_B_barrier_small/report.json`）恰好回避了这一点——**Newton 的"成功"来自密集线性代数，不来自障碍法优于内点法**。等价的正确做法是给成熟 IPM 一个 n×n 密集正规方程（7832² 的 Cholesky 约 1 s），或用 MOSEK（自动处理密集列）、CVXOPT 自定义 KKT。
2. **未采用作者的阈值/自旋重缩放。** `FOLLOWUP_METHODS_ZH.md` 明确记录 2403 用 Λ_ℓ(s)=((√s−2)/(√s+2))^{ℓ/2} 重缩放幺正锥并对 Gram 换基"消除已知运动学小因子"，并评价"不是把接近 1 的原始 S 直接塞给机器精度 PSD 程序"。仓库 09-06 的 2e−74 行尺度正是没做这一步；后续选择用 long double、Arb 768-bit、"有界坐标"（`PV_CORE_CLAIM_AUDIT_ZH.md` §3）去硬扛病态，而不是重参数化。
3. **预算把正在收敛的求解判死。** `D3_joint_001` 26 步 primal 2.5e−7；典型 IPM 需 30–80 步；再给 10 分钟大概率 Solved。09-11 两次全锥同样 600 s。`AGENTS.md`"bounded runtime"被执行成了"5–10 分钟"。作者的 CVX+MOSEK 就是解这一规模；`OUR_METHOD_REVIEW_ZH.md` 引 2511.11513 报告 M49/nmax37/lmax11 约 1 分钟。

至于 PV 核的"病态性"：它没有阻止自研 Newton 出结果，也没有阻止 Clarabel 在 390 维凸包上 0.26 s 收敛，故不是现成求解器失败的原因。17 GB 只是历史结果堆积，与单次求解无关。

## 3. 每个分支的功过

- **finite-sine（sine-cardinal 族）**：想给双谱积分一个明确插值以得到全能量可求值的解析振幅；论文未规定（`OUR_METHOD_REVIEW` 表）。后果：引出 T0zero/五尾/FG 一串必要条件，两天内把问题做成 8221–10062 行。09-07 退休，09-10 以 analytic-cardinal 之名回归。09-12 不在主线。
- **T0zero**：本族 f00→5T0/2 的全能量幺正必要条件；自造。2403 代码 `sig` 保持自由。改变可行集；09-07 退休，09-11 在 3123/3582 盘分支再次强加。PV 主线不用。
- **五尾**：S 波单密度首项半空间＋三个非零 spin 尾 SOC；自造。09-07 退休；09-11 以"五条领先高能条件"复活。PV 主线不用。
- **global-FG**：固定 s、ℓ→∞ 多项式在 [−1,1] 非负，SOS 提升成 6 个 25 阶 Gram PSD（Roh–Vandenberghe）；自造。使 Clarabel 问题膨胀到 1e8 非零并成为"Clarabel 失败"的主要样本。09-09 退休。
- **PV（pv-midpoint/Legendre-Q）**：来自 2403 代码（`COLLOCATION_KERNEL_MAP_ZH.md`）；09-07 起主线。功：把 1500 盘的论文原型接上。过：09-10 证明它不是单一解析函数（留数障碍），离节点值"明确拒绝"，于是既不能做 5 MeV 网格检查，也逼出了 analytic 分支。至今主线。
- **B 正则化（377500=100×3775）**：来自 2403 `norm(rho/Mrho,4)<=100`；2309 未给。功：使有限问题紧。过：`REGULARIZATION_VALUE_PROOF` 证明 B/2 使 x 上界降 ≥.0055，它是真实的模型输入而非预条件；09-07 的"1%/两 decade 平台"是自设门槛。至今主线，且很可能是 χ 区域 .1059 vs .0826 的主因之一，从未定量归因。
- **combined/separate χ**：论文"some norm"。separate 是 09-06 自设；combined-8 有 Fig.5 指纹（偏离 ≤.06%）+ 作者代码双重支持，09-09 采用；09-10 以"输出辅助"退回 separate；09-12 证明 combined 严格排除当前 tip（范数 .00211>.002）。后果：B3/C/E 两次重算，最终"严格证书"建立在来源依据更弱的那个读法上。
- **hard-midpoint vs clipped**：FESR 截止求积；论文 Eq.3.72 未给 mask。clipped 09-08 用；hard 09-09 以"更字面"采用，尽管 `SCIENCE.md` 自记常谱下 hard 相对误差 14.6%/28.9% vs clipped 0.9%/2.9%——比 raw .002 容差（相对 1.4%/1.8%）大一个量级。字面主义压过数值常识。至今主线。
- **analytic-cardinal**：09-10 为解决 PV 非解析而建的一致族；交付了完整 C/E/F 与 Fig.11 五组；代价是暴露离节点 η>1、T0≠0，随后 3123/3582 盘。09-12 标为"历史"。
- **1500/3123/3582 盘**：1500 是论文原型；3123=+ρ 窗交错行+端点等式+5 领先条件；3582=+459 行 5 MeV 网格。全部为修自选族的自造违例；Phase I 跑 7 轮才得 τ<0。不在 09-12 主线。
- **Arb 高精度证书**：384/512/768-bit 原式 primal/dual 复核，功在把"求解器状态≠证书"落实；过在把精力引向 1e−14 中心门槛、μ≈1e−10 的对数值抵消（`LOCAL_DESCENT_PROOF_ZH.md`），而物理差异是 46 MeV 与 28% 区域。`SOLVER_DECISION` 自己写"不能用测试数或高 bits 抵消差异"，但工作流仍在这样做。
- **数值中心作代表点**（未列入题目但关键）：χ 障碍权重 750=n_disks/2，代表点=固定 μ 的解析中心。`PV_CORE_CLAIM_AUDIT_ZH.md` 承认"正的障碍权重会改变有限 μ 中心；若极值不唯一，极限选取可能依赖权重"。论文用 CVX 直接取极值点。这是影响相移读数的自造规则。

## 4. AGENTS.md / 工作流层面的问题

对照 `AGENTS.md` 与 09-10 前的版本 `expert_response_20260910/AGENTS_BEFORE_MODULE_RELAXATION.md`：

1. **"one scientific calculation at a time with a bounded runtime… avoid broad parameter, solver or boundary scans"** → 被执行为 300–900 s 硬时限与"不做 solver 对照"。这条直接造成 §2 的全部 MaxTime，并让 `SOLVER_DECISION` 只能以"没有对照"收场。它本意防扫参，实际封死了换求解器的路。
2. **"Keep exactly ten active Python modules… ≤350 lines/24KiB… at most three test files"**（旧版；新版放宽到"ten is not a hard limit"）→ `conic.py`/`linear.py` 全是分号压缩的单行代码（如 `nv=len(P.active);nd=P.p-P.free;total=nv+nd`）。字节上限逼出不可读代码，CG 残差 .04 被当方向用（`B_SOLVER_FIX`）这类错误正是在这种代码里滋生的。同时 6 次模块整理（15→10→12→13→17→21→23）本身消耗了 09-09/10 的大量精力。
3. **"Use published output curves for comparison, not for choosing new physical inputs… Label historical output-assisted choices"** → 把 Fig.5 残差指纹+作者代码共同支持的 combined-8 χ 判为"输出辅助"而退回 separate-L2。这是把"来源鉴定"误当"拟合"。结果是 09-12 的"严格排除证书"针对的是来源依据更弱的模型，其"论文不兼容"结论因此打折。
4. **"Do not turn invented regulator platforms or continuum certificates into prerequisites"** —— 这条是 09-07 从 B 阶段吸取的正确教训，但 09-10/11 的 T0=0、FF 双零点、五领先条件、3582 盘又重犯，说明规则没有拦住"论文没规定=可以发明"的惯性。`AGENTS.md` 缺一条对称规则："论文没规定的，优先取作者公开代码的选择，其次取最简单、最接近字面的选择，不得自造更强条件。"
5. **"Preserve historical proof inputs and failures… Keep the removed legacy tree removed"** → 450 个 run 目录、17 GB、每份 STATUS 都有"历史原样保留"的长尾，阅读成本本身阻碍决策；`SOLVER_DECISION` 自认"移除长期没有收益的具体分支以减少维护与上下文负担"——即为了上下文而删掉现成求解器代码。
6. **以测试数/bits 替代复现**：STATUS 每节报"97/101/112/120/79/86/91/99/106/107 项测试通过"、"768-bit"、"1e−14 中心"，同时 P1 三点跨 107 MeV、L 跨 46 MeV、χ 区域差 28%。文档在措辞上反复否认"测试≠复现"，行动上却把最后三天投入 Newton 精度而非模型识别。
7. **"Never relax fundamental unitarity"** 与 "Sampled feasibility is not continuum certification" 两条相互拉扯，是 09-06 加强分支、09-11 3582 盘的直接思想来源：既不许放松，又不承认有限采样，只能不断加约束。论文本身只做 1500 盘采样。

## 5. 结论

**仓库是"没苦硬吃"为主。** 论文的核心计算是一个 4076 变量、1500 SOC(3)+100 PSD(3)+少量线性/L4 约束的小型 SDP，作者用 CVX+MOSEK 直解，作者也明确 2309 无需额外幺正化步骤（`MAJOR_SOURCE_NORM_AUDIT/ITERATION_SCOPE_ADDENDUM_ZH.md`）。真实困难只有一层：2309 未公开 χ 范数打包、密度正则化、FESR 误差打包和代表点规则；但 2403 代码把前两项给了，仓库采了 PV 核和 B 却拒了 combined χ。

**三个最大误判：**
1. 09-06/07 先把论文原型"加强"（T0zero/五尾/FG/8221 行），把病态缩放的加强问题在 300–1200 s 内交给 Clarabel/CVXOPT，据此认定现成求解器不行，转而锁定自研 dense-QR Newton；而 09-07 12:27 的自家审计已写明应"先核对到成熟 conic 基线"。
2. 09-08 `D3_joint_001` Clarabel 在收敛中被 300 s 杀掉，却退到 390 维凸包，再用自研 Newton 做 E/F 中心化，把 09-09→09-12 全部耗在 Newton 病理（QR 复用残差、慢中心化、不精确方向、1e−14 抵消、有界坐标）上。
3. 09-10 以"输出辅助"退回 separate-L2，并以障碍权重 750、μ→1e−10 的数值中心作代表点，然后在 09-12 证明这个自建模型"严格排除"论文曲线——把自造模型的失配当作对论文的判断。

**如果重来，应停下改走 SDP 直解的节点：09-07 13:13（`SCAN_PLAN_ZH.md`）。** 那时 PV 核、1500 盘、B=100Mρ 已齐，只需再取 2403 的 Λ_ℓ 重缩放、Gram 换基、combined-8 χ，在 cvxpy 里用 Clarabel/SCS（给密集正规方程或 ≥1 小时预算）或申请 MOSEK 学术许可，先解 max f00(3) 与三个代表点，再直接比 Fig.9–11。最晚的补救点是 09-08 `D3_joint_001`：把 300 s 改成 3600 s 即可知道答案。
