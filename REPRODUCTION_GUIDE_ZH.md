# 2309.12402v3：复现路线与运行接口

> ## ⚠ 读之前：这个仓库有两条平行主线
>
> **(1) Newton 主线**（`results/runs/mainline_alignment_20260909/`、`major_claims_20260912/`）——
> 本文件描述的就是它。它已被 `REVISED_CLAIMS_ZH.md` 的只读审计判定为
> **"解的是一个自建的近似模型，不是论文的问题"**（论文没有的密度正则化 B、比字面更紧的 L2 手征球、
> 以障碍中心而非极值点作代表、自研 barrier–Newton 求解器）。它的 P0–P5 计划已不再执行。
>
> **(2) SDP 直解路线**（`src/smatrix_bootstrap/sdp/`）——当前在做的。求解器已改为 **SDPB-only**
> （任意精度），MOSEK/MATLAB/CVX 都不在链路上。
>
> **当前进度、运行方法、精度判断与待补缺口的唯一入口是
> [`HANDOFF_PHYSICS_AUDIT_ZH.md`](HANDOFF_PHYSICS_AUDIT_ZH.md)。**
> 本文件保留为历史记录，其中与下列事实冲突的陈述以 handoff 为准：
>
> - 当前SDP路线同样通过统一入口：`python -m smatrix_bootstrap.run sdp ...`。
>   下文旧Newton命令及 `scripts/sdp/pmp_run.py` 不作为当前生产入口。
>   Mathematica用于独立公式检查，SDPB用于优化；自包含运行、精度理由和真实资源范围见handoff。
>

## 2026-09-12当前入口与下一阶段

[完整严格审计](results/runs/major_claims_20260912/final_report/REPORT_ZH.pdf)、[六项未闭合工作、A1–F3台账及P0–P5计划](results/runs/major_claims_20260912/CORE_CLAIM_LEDGER_AND_PLAN_ZH.md)。当前定性IR→UV→共振链及五组有限计算已完成；稳健定量ρ、S0、L比较及原文区域/数值身份仍按台账保留，不以代码/图齐备宣称物理成功。

计算仍只经 `python -m smatrix_bootstrap.run`。本轮新入口 `support-probe` 在原可行集上优化线性组合，并给原支撑附近的观测量上界，不添加物理约束，也不生成代表selection。`--probe-audit-only` 可重验保存的完整点/对偶；失败路径的 `joint_path_candidate.npz` 只在目标/模型身份匹配并重新原式验收后使用。

`--probe-anchor`、`--probe-node`、`--probe-wave`、`--probe-component`、`--probe-sign` 明确目标；`--probe-xref` 配合两个 `--support-runs` 在audit-only模式下由完整父点证明精确截面的可行内界。所有原H/current、χ/B/UV与 `--start-mu` 仍须显式给定。辅助组合目标未达最优性时保持false，严格负的条件上界可独立完成排除证明。不要将probe输出送进论文代表选择。

23个模块均≤350行/24KiB，三个测试文件均≤350行，107项通过。新增Arb下降核验的范围与完整推导见[数值附录](results/runs/major_claims_20260912/LOCAL_DESCENT_PROOF_ZH.md)；物理primal、原式支撑和数值中心资格分开。聚合入口现在记录显式支持来源，选择/区域重验与原输出逐项一致。

先完成计划P0/P1的来源/算子问题，再决定是否有依据改变合同；合同不变时，本轮近tip和精确xref的相位排除继续成立。不要通过改变质量、范数、截止或后续泛函来追图。

[4页执行计划PDF](results/runs/major_claims_20260912/final_report/PLAN_ZH.pdf)。

## 2026-09-11及以前的记录（原样保留）

[正式中文报告（18页）](results/runs/physical_consistency_20260911/PV_mainline_report_final/REPORT_ZH.pdf)。

## 当前入口：PV有限链与同模型五组F均完成（2026-09-11）

按原文约束的mixed-PV有限实现已完成独立IR→UV→共振信号链；稳定合理的ρ和三分波定量比较仍有major差异，不能称为全复现成功。先读[STATUS](STATUS.md)、[SCIENCE](SCIENCE.md)、[PV代表规则](results/runs/physical_consistency_20260911/PV_REPRESENTATIVE_RULE.json)和[F预登记规则](results/runs/physical_consistency_20260911/PV_F_REPRESENTATIVE_RULE.json)。所有科学计算仍经 `python -m smatrix_bootstrap.run`，使用新的 `results/runs` 输出目录。

| 项目 | 当前PV值 |
|---|---|
| 散射准备 | `results/runs/stage_B_pv_M50_L10_prepare_20260907` |
| 电流准备 | `results/runs/mainline_alignment_20260909/D1_hard_M50` |
| 处方／分辨率 | `--prescription pv-midpoint --nodes 50 --waves 10`，1500盘 |
| IR／正则化 | `--chiral-norm separate-l2 --chiral-tolerance .002 --density-limit 377500` |
| UV／无穷远 | `--moment-source printed --sr-error raw-absolute --fesr-cutoff hard-midpoint --mq-rule arithmetic-mean --infinity free` |
| 中心选点 | `--require-center --chiral-barrier-weight 750`；+x tip、物理xref和中间截面；相位前冻结完整C/ImF/R |

已完成链为[PV_IR_ref_05](results/runs/physical_consistency_20260911/PV_IR_ref_05/report.json)→[PV_IR_selected](results/runs/physical_consistency_20260911/PV_IR_selected/report.json)→[PV_IR_native](results/runs/physical_consistency_20260911/PV_IR_native/evaluation.json)，以及[PV_UV_selected](results/runs/physical_consistency_20260911/PV_UV_selected/report.json)冻结的tip/mid/ref→各自原生求值→[PV_native_phases](results/runs/physical_consistency_20260911/PV_native_phases/phases.json)。IR支持gap=.00159036，UV tip/ref/mid的gap分别为9.22923e−5/.001106419/.001131500，四份幅度的528项主波检查均sampled_passed。

相位JSON保存IR P1末值8.56389°、UV tip/mid/ref的90°读数803.386/701.114/696.561 MeV、min η(P1)=.775985/.434728/.243206及严格邻点峰证据，原生P1强度主峰节点为792.136/680.414/680.414 MeV。同文件的ref/mid支持平面以正余量.09777174521/.08328539559排除该冻结IR幅度的UV扩展。90°插值、离散峰和投影排除各有独立范围，不证明极点、连续散射完成或全部IR幅度的性质。

F已复用[四组H/current与完整初值](results/runs/physical_consistency_20260911/PV_F_REUSE_INVENTORY.md)，保持每组配套目录。流程为旧点在新separate-L2下原式重验、必要时联合初始化、新+x中心支持、同配置 `resolution --require-center` 选择、原生 `evaluate --primary-waves`，最后五组 `resolution --profile-runs ... --require-center` 汇总。B(M)=100[M²+M(M+1)/2]；χ障碍权重显式取原生盘数/2。旧combined点只作初值，L12的非正谱初始化缺陷修补后已有完整见证和最终中心支持，旧失败不构成不可行证明；不能继承旧中心或支持。

[五组完整F](results/runs/physical_consistency_20260911/PV_F_five_complete/resolution.json)已通过全部中心、原式支持与669项原生检查。[Fig.8截面](results/runs/physical_consistency_20260911/PV_E1_regions_complete/regions.json)两侧严格收缩，上/下比[1.00317959428,1.10388598635]。固定M的L稳定性、S0形状和三代表ρ定量仍有major差异，详见[最终审计](results/runs/physical_consistency_20260911/PV_CORE_CLAIM_AUDIT_ZH.md)。

mixed-PV是本仓库的具名有限实现。作者的散射有限插值、χ范数及B等离散约定未唯一恢复，不能称为原代码。保持已保存H/current，不添加解析族端点／高能条件或额外采样；PV只在原生物理节点及已定义阈下范围求值，连续证书不是该有限交付的前置条件。

已有[完整F五配置](results/runs/sequential_reproduction_20260910/F_five_configs_verified_comparison/resolution.json)属于analytic-cardinal，不能代替本轮PV F。解析加强模型W1诊断已完成，W2未启动。最新主线为107项测试、21核心模块和3测试文件；严格原生峰与投影排除已在同CLI交付，测试通过不代表论文全部物理claim完成。

后续源码的误差球、输入及节点变化只记录在[一手审计](results/evidence/PV_PRIMARY_SOURCE_AUDIT_20260911/AUDIT_ZH.md)，不作为当前CLI默认替换。

本次只更新顶部摘要；[此次修改前全文/SHA](results/evidence/pv_final_summary_before_20260911/manifest.json)及[修改前全文／SHA](results/evidence/pv_chain_summary_before_20260911/manifest.json)保留为非可执行证据。以下旧接口、默认范数和“当前”等文字均按历史记录阅读，不覆盖本摘要。

**下文保留历史说明。** 其中PV、十模块、旧Fig.8标记选点、以及“解析后端尚未实现”等描述属于此前分支；当前入口、模块规模、选点和可行性范围以上述摘要及STATUS顶部为准。历史命令、失败和原始结果没有被改写或追认为新解析分支的证据。


先读[STATUS](STATUS.md)的当前18步骤，再读[SCIENCE](SCIENCE.md)的公式和来源。既有有限结果的冻结记录是[PAPER_MAINLINE](results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json)：PV原生节点、合并八维手征L2、打印四矩/raw误差、hard-midpoint及固定物理输入。旧separate/clipped运行保留作对照，不继承其完成状态。

## 1. 论文真正要回答什么

目标是把高能处的夸克／胶子信息与低能处的 pion 信息结合，约束中间能区的 pion 散射。论文用 Nc=3、Nf=2 的例子检验这种方法。它的关键结果是一组连续的比较：基本散射条件允许很多振幅；手征关系把区域压成很窄的一条，并能较好解释 S0/S2；P1 仍不合理；加入电流、形状因子和 QCD 求和规则后，P1 才出现较合理的 rho 共振。这种“加入哪类物理信息，改善哪一部分结果”的因果比较，比单独画出一条像实验的曲线更重要。

S0、S2 的角动量都为0，同位旋分别为0、2；P1 的角动量和同位旋都为1。分波 S=eta exp(2i delta)，delta 是相移，eta 表示留在弹性通道中的程度；模型约束 eta≤1，未强制 eta=1。允许非弹性过程不违反总概率守恒。

幺正性在这里有明确的可计算标准。论文第9页(2.9–12)给出 \(S=1+i\kappa f\)、\(\kappa=\pi\sqrt{1-4/s}\)，所以应逐波检查

\[
1-|S|^2=\kappa(2\operatorname{Im}f-\kappa|f|^2)\ge0.
\]

它等价于第24页(3.71)的散射2×2矩阵 \(\begin{pmatrix}1&S\\S^*&1\end{pmatrix}\) 半正定；只检查 Im f≥0 会漏掉真实违例。精确阈值s=4处，若f有限则S=1、裕量为0；阈下的振幅用于手征约束，不用上述物理幺正disk检验。共用`unitarity_margin/scattering_gram`及弹性、非弹性、违例解析测试固定了这一归一化。

判据正确与覆盖完整是两件事。计算器可以在所列能量和分波上报告通过、违例或无法判定；50点或150点都没有检查所有能量，增加几个高分波也没有覆盖无限尾部。sine的Arb角积分尚未包含严格积分余项包络；PV则采用声明的有限节点程序。物理主张须结合对应能区／分波的偏离量级、稳定性和已知违例审查，有限程序通过不等于连续认证。

区域图横轴 x=f00(3)、纵轴 y=f11(3)。s 是以 pion 质量为单位的能量平方；物理阈值在4，所以3是解析延拓的坐标，并非一次真实碰撞。图上的一点仅提供两个数，不能代替完整振幅；不同振幅可能投到同一点，即使在边界也没有唯一性的证明。相移必须由保存的完整解恢复。

[Snowmass 白皮书](references/2203.02421v2.pdf)提供方法背景，不是与本论文并列的第二套交付。标量值2.6613是方法校准参考，既不约束 QCD 电流，也不代表已经复现了 pion 区域的一部分。

## 研究依赖与顺序

A统一振幅、投影和输入；B用同一有限模型产生合法支持/区域；C锁定B的完整绿色幅度，依次复现阈下、代表点和IR-only相移。D独立建立电流/FF/UV算子后与同一散射幅度联立；E比较IR与UV区域，预选三套完整幅度后检查相移与弹性；F只改变原文的五组分辨率。

B2纯S全区域是基线旁路。C代表点只需对应截面的已证明上侧间隙，不必等待B3其它方向的几何细化；全区域误差仍单独报告。D见证只证明共同可行性，不能替代E的全幅度支持优化。支持完成也不自动证明物理曲线稳定。

## 固定的计算接口

所有生产计算使用`python -m smatrix_bootstrap.run`，数据与运行参数写到新的`results/runs`目录。当前环境：

```bash
export PYTHONPATH=/tmp/collocation_arb:src
export SMATRIX_BLAS_THREADS=2
/home/shiqiu/miniconda3/bin/python -m smatrix_bootstrap.run --help
```

| 步骤 | CLI | 关键输入／输出 |
|---|---|---|
| A3 | `prepare` | 明确`--prescription pv-midpoint --nodes 50 --waves 10`；输出原H |
| B1/B3 | `boundary --mode chiral` | 原H、ε、方向、B；输出完整系数及支持区间 |
| B2 | `boundary --mode pure` | 同一模型，关闭χ球；输出pure支持 |
| B汇总 | `regions --support-runs ...` | 同模型合法支持形成内外包络 |
| C1/C2 | `select`，`evaluate`，`profiles` | 几何选点后固定全部系数；阈下与物理求值分开 |
| D1/D2 | `prepare-current` | 必须绑定同一H；输出Gram/FESR/FF算子 |
| D3 | `boundary --mode gauge --joint-feasibility` | 完整联合可行见证 |
| E1/E2 | `boundary --mode gauge --solver centered` | 全幅度支持；tip为+x，另两点固定预定x最大化y |
| E1汇总 | `gauge-regions` | 同模型IR区域及UV原式支持 |
| E3 | `evaluate`，`gauge-phases`，`compare` | 已冻结tip/mid/ref完整向量、原图CSV |
| F | `resolution` | 既定(50,8)/(50,10)/(50,12)/(45,10)/(60,10) |
| 原式复核 | `dual` | 保存的完整点及可重验的同算子对偶 |

PV支持必须显式指定`--unitarity-scope sampled --infinity free`。物理优化、区域汇总及新分辨率求解必须显式指定χ范数；不再隐含默认combined-l2。生产UV固定：

```text
--chiral-norm combined-l2 --chiral-tolerance .002
--density-limit 377500 --moment-source printed --sr-error raw-absolute
--fesr-cutoff hard-midpoint --mq-rule arithmetic-mean
```

B(M)=100[M²+M(M+1)/2]是统一声明的离散正则化配方；跨M必须按此更新。物理参数、FESR容差与FF界不按相位调整。原H路径为`results/runs/stage_B_pv_M50_L10_prepare_20260907`，当前M50电流路径为`results/runs/mainline_alignment_20260909/D1_hard_M50`。

同目标续算读取实际`barrier_state.npz`的μ；不是把任意可行混合当作小μ中心。新目标或新模型须重建并原式验回。可逆换元、QR预条件及正障碍权重属于数值实现；删掉物理方向或放宽幺正性不属于允许的加速。

## 选点、相位与比较

C使用同一物理xref的蓝/橙/绿上侧代表点；C2/C3原样复用绿色。E的tip是完整+x支持，另两点的x按原Fig.8标记与黑点的比例固定，完整规则记录于PAPER_MAINLINE。所有选点先于相位，无实验相位筛选。

PV物理求值使用原生节点；不支持任意物理离节点延拓。相移取从S(4)=1展开的arg(S)/2，η=|S|。90°线性读数是描述性共振指标，不是解析极点定位。原图bootstrap、phenomenology和实验比较分列；原图坐标质量的映射仅用于同原生s对照，不能改变物理输入。

2403的Watsonian是具名后续方法，当前原2309投影主线与F不自动加入它。使用时必须保存旧F构造的固定目标、完成该目标的原式支持、报告移动后的x/y，并把目标最优与新S/F固定点分开；不能用它改写原始代表点的身份。

## 科学验收与执行规模

幺正性、χ、actual-ρ、同一S的电流Gram、FESR和FF必须共同满足。独立数学检查保留在三个测试文件中；测试及原式有限证书不能替代区域、物理曲线和分辨率。全能量/无限高spin的连续认证不作为论文有限原型的额外前置任务，也不冒称已完成。

默认一个有明确物理决策的有界计算；完全独立的任务最多两个生产进程。运行中不改生产源码。同一问题有实际中心/gap进展才续算；一个科学比较已有有效界时，不为额外小数位无限推进。不得通过参数扫描贴图或把数值失败称为理论不可行。

旧接口说明和全部历史设置见[整理前指南](results/runs/mainline_alignment_20260909/GUIDE_BEFORE_CONSOLIDATION.md)；失败、原系数与证明输入留在原位置，继续当前主线无需默认扫描。

## 当前代码与重放兼容

2026-09-10已归并为指定十模块，逐文件限额及79项检查通过。当前仅保留PV生产入口；sine/FG/五尾/旧全空间conic、embed/regulators执行分支已退休，其证明和原producer保持独立历史范围。[接口与退休清单](results/runs/mainline_alignment_20260909/core10_review/REVIEW_ZH.md)说明从旧report重放时应剔除的无效CLI参数；现有PV系数、联合点、两类缓存、selection和report数据仍可读取。

`resolution`对少于五个已验配置只能输出明确partial；固定L/M子比较的complete只表示配置齐全，不代表原图稳定性自动通过。原Fig.3–11及η图仍生成，额外差值/FF-lift/Watson调试图已退出自动绘制，对应数值JSON/CSV及诊断保持。当前F已完成四组，L12联合可行性仍未闭合，详见STATUS。

## 专家反馈后的机制审计接口（2026-09-10）

模块数量按用户新指示不再固定为十个。新增 `spectra` 与 `certificates` 分别承载原生谱诊断和原式联合证书，现为12模块；仍用同一CLI、三个测试文件。数学与结果见[新中文报告](results/runs/expert_response_20260910/RHO_MECHANISM_AUDIT_ZH.md)。

`gauge-phases` 增加 `--baseline-profile`：显式传入同模型C原生求值目录，输入不一致会报错，避免静默省略IR对照。输出增加 `rho_mechanism.pdf/png`，并在 `phases.json` 保存全部离散峰、Schur谱预算和矩分解。不会重新选点或平滑数据。

固定振幅检验使用 `boundary --mode gauge --joint-feasibility --fixed-amplitude --direction 0 0`。继续提供原preparation/current-preparation、全部原约束参数与 `--coefficients`；region-summary只提供同模型签名。它固定整个C向量，只寻找电流扩展，不求全幅度支持。`dual` 接受同一组固定振幅参数以重验保存的对偶。负零目标上界只排除该振幅；`fixed_amplitude_inconclusive` 不得解释为不可行。

本轮完整命令在 [C_fixed_current_fiber_03](results/runs/expert_response_20260910/C_fixed_current_fiber_03/report.json) 与 [768-bit replay](results/runs/expert_response_20260910/C_fixed_current_certificate_768/report.json)。原式约束没有变化；本轮固定部分的精确代入与λ=1消去修正的是数值构造。全部历史失败及当时producer仍可重查。

## 基础复审后的接口与解释（2026-09-10）

当前新增 `analytic-audit`，使用独立推导的完整解析比较族和Arb严格角积分。`--profile-runs` 传冻结系数目录；`--comparison-manifest` 可复用已完成角积分，补全代数核差界、留数与输入审计。结果不是旧PV点的可行性转移，不能送回旧支持汇总冒充同一模型。

新 `select-gauge` 禁止 `--selection-reference` 作为选点输入，采用xref和已计算tip的几何规则。旧Fig.8标记选点仍能通过其保存记录审查，不会被重新标成独立结果。底层问题、全部A1–F3证明范围与下一次计算协议见 [REVIEW](results/runs/foundations_review_20260910/REVIEW_ZH.md) 和 [INDEPENDENT_PROTOCOL](results/runs/foundations_review_20260910/INDEPENDENT_PROTOCOL.json)。

现有PV求解器仍只交付声明的有限近似。新 `analytic` 尚不是全套B–F优化后端；共同解析性障碍不能通过继续画旧曲线或提升bits消除。全论文完成状态以新STATUS/claim台账为准。

本轮原文主线判定见 [逐图复审](results/runs/followup_research_20260910/MAINLINE_ALIGNMENT_ZH.md)。analytic-audit 的 comparison-manifest 后处理现附带固定系数的积分分辨率充分界和局部代数反例。该功能不改变生产H、系数或物理约束；全能全自旋认证和唯一极点不作为原Fig.3–11的新增前置条件。
