# SDPB-only 主线交接：2309.12402v3

> ## ⚠ 本文件是历史记录（2026-09-13 起）
>
> 现行入口只有 [`PLAN_SDPB_2309_ZH.md`](PLAN_SDPB_2309_ZH.md) 与 [`TRIAGE_ZH.md`](TRIAGE_ZH.md)。
> 本文件中以下说法已被证伪（证据与正确版本见 TRIAGE §1）：
> "2309 没有正则化、密度界 B 应删除"；"去掉 B 后 M=50 才认证、失败是精度或截断问题"；
> "M=50 不收敛、双精度结果只是下界"；"mixed-PV 需换成 sine-cardinal"；
> "chi-c、SR-a 是 2309 字面合同"。当前主线：mixed-pv、chi-b、ε^χ=0.002、
> 正则化 |ρ_ij|≤10²（物理判据定），SDPB 192 bit；门槛已通过，端点与论文差 0.22%。
> 本文件其余内容保留作证据，其中的数值、结论与"当前"等措辞不再代表项目状态。
>


当前主线未完成论文复现。唯一有效入口为 `python -m smatrix_bootstrap.run sdp`。Newton、MOSEK、MATLAB/CVX 优化均不在当前计算链路；历史结果不能跨解析族移植。

[有限离散化解释更正](results/runs/sdpb_mainline_20260912/FINITE_DISCRETIZATION_SCOPE_CORRECTION_ZH.md)：不以“有限M下不是同一个精确函数”单独否定混合配点近似；sine-only是当前已验证路线的身份约束，非论文唯一要求。既有数值和物理失败保持原记录。

本次完整源文研究见 [中文报告](results/runs/sdpb_mainline_20260912/deepresearch/REPORT_ZH.md) 和 [PDF](results/runs/sdpb_mainline_20260912/deepresearch/REPORT_ZH.pdf)。它核对了原文、2403/2505后续方法、物理推导、精度和代码复用。修改前本文件及相关代码的原字节与哈希保存在 [归档清单](results/runs/sdpb_mainline_20260912/deepresearch/before_sha256.json) 和同目录 `before_source.json.gz`；旧交接不再作为当前指令。

## 已证实的断裂与修复

1. 旧 double PMP 在SDPB内收敛，但原式严重违例。源算子、投影、PMP及完整y回读必须贯通高精度。已实现40位源目标和30位输出；最初192bit，当前按实测故障保留320bit SDPB，不加B掩盖误差。
2. 高精度mixed-PV运行确实通过其全套有限原生约束，但P1首次90°约562MeV，S0/S2仍偏离；它不是精度不足的证据。
3. mixed-PV的cut/off-cut有限规则不能同时解释为所写有限有理函数的精确边界与延拓；这不单独否定其作为同一色散表示的配点近似。归档明确写midpoint approximation，求积误差与M/L稳定性仍待验证。当前sine-cardinal在所有腿使用同一已声明函数族，有M8/MMA控制，但不是论文唯一指定的插值。
4. 原默认chi-b是合并八维球，SR-b是10%相对误差，frozen把P1最后高能节点的平方界放宽约240倍，不能称为2309字面合同。新默认为chi-c/SR-a/node：两组L2仍是声明范数；raw误差0.002和逐节点FF直接对应原文。见 [源文决定](results/runs/sdpb_mainline_20260912/deepresearch/SOURCE_READING_DECISION.json)。
5. 代表规则tip/ref/mid是预登记的截面代理，不是已证明的最近物理边界点。相移前冻结完整振幅；不按相移、rho位置调整输入或选点。
6. 已修复验收缺项可通过、prepared旧约束配新算子目标、IR复用被无效UV字段阻止等问题。测试通过不表示论文claim完成。

原 [预登记](results/runs/sdp_reproduction_20260912/preregistration.json) SHA256仍为 `3fe62f034bc11e6908a600d0c5e09f0a98504c6befe20f7c86eeb19fb1acfa4a`，原字节不改。其“FESR量纲不闭合”等解释不成立；旧九组合扫描、B扫描不再执行。新源文合同与旧预登记分别记录，验收阈值不为结果放宽。

## 当前计算与已完成的同基底对照

优化停止误差对照已完成：`IR_M50_L12_eps00021_gap1e8_warm/report.json`保持同一个804,546,355字节PMP、basis、mask、模型和320bit，只将gap1e−6收紧到1e−8；6步主迭代、总423.48秒（7.06分钟）。原式数值验收通过，新目标仍在旧数值支持区间内。E≤1.2GeV最大相位变化S0/S2/P1为.706173°/.112199°/.056124°，旧新C4_source均FAIL；旧源checkpoint原字节不变。本次低能停止误差已小于预先设定的1°诊断预算，不以更严gap替代截断定位。5.944GeV处η变化可达.83–.86，不能声称整幅度稳定或唯一。证据见 `optimization_gap_control_20260913/REPORT_ZH.md`。短时与L13同时各8ranks，现只剩L13在跑。

名义ε=.002的C3部分源图对照也已由既有131点评价直接完成，90个源横坐标全部命中保存网格，零新求值。S0/S2/P1的RMS/f00(3)为66.17%/20.49%/3.06%；前两波超过原8%预算，P1通过。绿色源图的描述性f00(3)≈.072453671，与本地直接值差1.18%，仍有样本符号/形状差异。作者高亮点身份未认证；完整三ε C3仍not run，未把无采样过零称作全局无零。见 `nominal_eps002_C3_partial_comparison_20260913/REPORT_ZH.md`。

当前主要截断优化为 `IR_M50_L13_fixed_e541_eps00021_upper/ref/report.json`，driver PID2598007、2026-09-13 13:16:48 UTC启动。它保持M50、χ=.00021、IR-only、xref上端目标及e541原basis，只把L12增至L13，增加I0/I2 ell24与I1 ell25的150个原生盘。属于一次截断辨识，不是Fig.11第六配置或连续证书前置门槛。320bit、源40/PMP30位、8ranks、原gap1e−6/原可行性1e−10及7200秒预算不变；同M50 MMA资料已实际校验复用，18000秒上限保留。计划、相位读取规则和启动记录见 `fixed_basis_L13_IR_diagnostic_20260913/`。

同M跨L基底运输已实现，仅放宽坐标矩阵的L相等项，新增source_L/target_L/cross_L记录；目标索引、mask、PMP和原式验收必须重新建立。prepared旧PMP仍不能跨L复用。13项定向测试、8组实际小PMP控制和完整222项回归（77.18秒）通过，见 `same_M_cross_L_basis_implementation_20260913/`。实际M50新旧PMP已按(I,ell,node)重排核对：1804个公共块、目标和normalization全部相同，basis哈希相同，仅150条新盘。原L12幅度在两条真实新PMP块上有严格负行列式；这只排除此旧完整幅度，不代替新最优值。见 `fixed_basis_L13_IR_diagnostic_20260913/ACTUAL_INPUT_COMPARISON.json` 与 `OLD_POINT_ON_NEW_BLOCKS.json`。已完成的17-case便携包保持冻结，它是此接口更新前的可验证版本。

χ=.00021 的同basis UV tip/ref/near已全部终态通过，完整选定IR＋三UV对比见 `FPI_CALIBRATED_FULL_CONTRAST/`：C4、C6、C7均FAIL。tip：x=.205068159213、y=−.0136248072904、63.7分钟，P1首次90°≈911.779MeV、minη=.469667；ref：x=xref、y=−.00443854642136、64.3分钟，无原生低能90°、minη=.140574；near：x=.0740546044436、y=−.00448739123093、12.49分钟，也无90°、minη=.144387。各波图及1116行支持CSV保持保存相位，论文η留空。实际basis/完整y/源合同均通过核对，只有IR ref与UV ref匹配目标和截面；不同代表的差不能全部归因于UV。

near截面已在新UV相移读取前由 `FPI_CALIBRATED_UV_NEAR_PLAN/report.json` 固定，随后 `FPI_CALIBRATED_UV_SELECTION/report.json` 核验三完整幅度。内部mid槽明确映射near，旧算术中点登记仍为未完成。近点没有因曲线差异再选择。保持M50/L12、χ=.00021、chi-c、无B、sine-cardinal、SR-a raw四矩、逐节点FF、mqmean和原QCD参数；此为fπ校准实验，不是恢复作者ε。

**C5上下截面已完成；当前L13对照独立于该已交付结果。** 同PMP续算 `UV_M50_L12_fpi_calibrated_eps00021_lower_resumed/report.json` 实际加载原终态检查点，7步主迭代、额外466.57秒（7.78分钟），其中SDPB264.32秒；gap2.35333e−8、primalError7.60864e−11，全部原式数值验收通过。原lower_warm的1200秒预算超时（总1463.45秒、未接受）完整保留；只延长timeout，没有提高位数或放宽阈值。

同basis、同xref、同χ=.00021下，UV lower y=−.00556223619319。上侧降低1.29561e−5，下侧抬升7.95028e−5，上／下收缩比.162964，数值支持误差区间[.150202,.174294]；论文Fig.8源图在同xref的分支插值约6.30322。上下均收缩，但不对称方向相反，不能以末位误差解释，也不能称完成全Fig.8区域。下点相位未用于选择或比较。图、CSV与来源见 `FPI_CALIBRATED_C5_SECTION/`、`FIRST_FPI_CALIBRATED_C5_LOWER_RESUMED.json`。

前一组χ=.00021 IR两端均已终态通过：上端2686.45秒（44.8分钟）、y=−.00442559034316，下端2753.12秒（45.9分钟）、y=−.00564173901187。`FPI_CALIBRATED_IR_GEOMETRIC_SELECTION/report.json`在新相移读取前选择upper；顺序证据为 `FPI_CALIBRATED_IR_PHASE_READOUT_ORDER.json`。原生盘分别1725/1723项严格通过、75/77项未定、0失败，均满足原数值预算。与旧χ=.002同xref upper相比，1804个PMP块只有两个χ块改变；1800原生盘、两截面块、目标、normalization和basis完全相同，见 `CALIBRATED_IR_ONLY_CHI_CHANGE.json`。

选定IR的S0阈下RMS/xref从.6335降到.08942，S2从.2135降到.02849，P1从.04145降到.00594；S0在[.45,.5319955]由Arb端点异号确认至少一个零，不声称唯一。C4源图比较仍失败：S0/S2/P1原生相位RMS约48.53°/3.94°/11.87°。逐点模180°的最小相位RMS为17.10°/3.94°/11.87°，最高公开能量处仍差34.35°/17.41°/67.52°，所以换支不能使完整比较通过。预先固定的五个独立角积分均吻合（27.98秒），首遗漏I0/I2 ell24、node42的∣S∣仍为6.68358/8.06104。结果、CSV及统一CLI图表见 `FIRST_FPI_CALIBRATED_IR_PAIR.json`、`FPI_CALIBRATED_IR_SELECTED_ANALYSIS/`、`FPI_CALIBRATED_IR_FIGURES_FINAL/`。

复用零幅、两端和.21倍旧χ=.001 tip构造内凸包，可在固定basis/原生数值scope下夹住最近边界距离：[.0004597107016150954,.0004633990329218793]，相对宽度.8023%。上端使用原out.txt的数值支撑上界，不把候选误作严格边界。没有全空间SVD等价、严格对偶或幅度/相移唯一性结论；也不代表作者未公开的距离度量。此精度已足以停止额外昂贵几何微调，见 `calibrated_IR_nearest_distance_bound_20260913/REPORT_ZH.md`。

这一容差由已验收IR upper、χ=.001 tip及零幅的凸组合取得物理黑点来确定，未用论文相移调参。高精度组合χ两组范数为.000194724286619/.000192844443310；保守三角上界.000203320960745，在χ=.00021下投影见证域有半径约1.0337e−5的黑点邻域。实际χ也小于论文网格中的.0002；保留.00021的保守余量，不为5%差异重跑。此为fπ校准实验，不是恢复作者值或全空间最小容差；数值来源下的凸复用论证不是严格全凸可行证书。见 `chiral_convex_witness_20260913/` 与 `chiral_tolerance_calibration_design_20260913/`。组合幅度未作为代表，也未读phase。

此前两项任务均已终态验收：名义 `UV_M50_L10_source_tip_20260913/tip/report.json` 耗时4001.15秒（66.7分钟）、x=.6061065134、y=−.04134905754、gap9.74675e−7，ρ穿越545.174MeV、P1 minη=.932437；仍未复现论文。独立角积分确认两个穿越相邻节点，首遗漏I0/ell20/node42的∣S∣=55.2958019，19.43秒。原700个Gram项严格通过；较松5.27秒电流诊断698项通过、2项符号未定、0失败，不推翻原验收。P1两pion谱矩比例仅20.62%/10.17%。名义旧ref与新tip的1624物理块、normalization和basis相同，已保存两点相位/η/CSV比较，mid未完成。

`IR_M50_L12_chi_001_tip/tip/report.json` 耗时4157.89秒（69.3分钟）、x=.4037901852、y=−.02710509066、gap7.10539e−7，原1800盘数值通过。它本身只证明+x范围，完整黑点准入证据来自上述同basis凸组合。两新解的完整y哈希均已实际核验，并导出精确Mathematica系数。结果见 `FIRST_NOMINAL_SOURCE_UV_TIP.json`、`FIRST_ADJACENT_CHI_001_IR_TIP.json`、`nominal_UV_tip_ref_comparison_20260913/`。

`IR_M50_L12_shared_basis/ref/report.json` 已完成：2664.193秒（44.4分钟）、gap8.50961e−7，x=0.07332139057293464，y=−0.002100421085515907。它与UV_L12使用逐字节相同的3876×1989基底（SHA256 `e541340318a1eb8c7deeb8c7c9725fae8c270998240a2371a5ff7808a03c8269`）。原生1800盘数值通过，1724项严格通过、76项符号未定、0项严格失败。IR/UV同截面上端点的数值收缩区间约[3.5473e−5,3.6849e−5]，不是完整Fig.8或严格对偶证书。

`IR_M50_L12_section_lo/report.json` 的同xref下端点也已完成：2459.740秒（41.0分钟）、y=−0.010036119335789418。读取下点相移前，`IR_L12_GEOMETRIC_SELECTION/report.json` 按到Weinberg黑点的距离选择了**上端点**，两端距离误差区间分离。规则只比较这两个截面端点，未证明全边界最近点。下点保持未选定，不为曲线更接近而替换。统一CLI为 `sdp ir-select --source-report 上端点leaf --compare-report 下端点leaf --out 新目录`。

`SOURCE_L12_PAPER_STYLE_REFERENCE/` 已交付相移、η、数据及来源报告：一个几何选定IR幅度与UV ref比较，实际Fig.7 IR曲线单独显示。基底和目标均匹配，可以量化这组有限问题的UV影响；三UV代表尚不齐全。IR的Fig.7 S0/S2/P1相移RMS差约49.3°/28.4°/15.4°，物理比较失败。

`UV_M50_L12_source_ref_p320/report.json` 已通过数值及原式验收，共3702.896秒（61.7分钟），gap=5.25175e−7、目标回读误差0。原1800盘有1725项严格通过、75项符号未定、0项严格失败；均满足既定数值预算。700个Gram子式、4个FESR和14个FF界严格通过。完整y的保存哈希已实际复核。见 `FIRST_SOURCE_UV_L12.json`、`L12_PRECISION_REPAIR_RESULT.json`。

原L12任务在192bit、主迭代81、gap≈8.15e−6时发生Cholesky(Q)数值正定性丢失；观测条件数≈4.26e67，不能将其判作物理不可行。旧out.txt/y.txt来自计时阶段，未获验收，不是失败时主迭代的物理振幅。320bit已经越过该失败并完成求解；最大观测条件数≈1.11e69，剩余十进制位估计约27。没有测试256bit，也未声称320是严格最小精度。

320bit重跑保持原L12 PMP、保存基底及盘掩码的确切字节，仍为M50/L12、chi-c、SR-a、逐节点FF、无B、sine-cardinal、40/30位源/输出、8ranks、gap1e−6和同一xref截面。仅将两阶段算术提高至320bit，并每300秒写检查点。原失败任务没有可用解检查点，因此此次是相同输入的新求解，未声称断点续算。当时仅实际写出8个非空rank检查点；后续near已实际验证加载，见下述warm-start记录。见 `L12_PRECISION_FAILURE_AND_RESTART.json`、`L12_P320_CHECKPOINT_OBSERVED.json`。L12来自论文附录，不是新增物理模型。

M50/L10的IR与UV参考均已完成原式数值验收（1798秒、2363秒），原生相移均未复现论文。IR在加UV前已有S0/S2明显偏差；UV无原生网格90°读数、P1 min eta≈0.2295。见 FIRST_SOURCE_IR.json、FIRST_SOURCE_UV.json 和 SOURCE_REFERENCE_CONTRAST/。这些离散读数不构成连续谱中没有共振的排除。

L10首个未施加波在低能原生节点有|S|=52.9057（I0,l20,node42）与1.84553（I1,l21,node38），独立积分确认。L12已施加这些波和l22/l23，但新遗漏波仍有|S|=7.85590（I0,l24,node42）及5.76876（I2,l24,node42）；I1,l25,node38为0.99217。三点的独立积分对照均通过，18.37秒。两项大违例证明该冻结L12幅度的截断仍未控制，不能靠提高算术精度修复。见 `source_UV_L12_first_omitted_waves/report.json`。

UV_L12的P1仍无原生低能90°读数，min eta=0.275914；S0在0.9/1.196GeV的相移约135.09°/163.80°，物理验收未通过。电流诊断5.54秒，S0两pion占两矩约92.4%/92.9%，P1约61.8%/63.6%。见 `source_UV_L12_ref_mechanism/report.json`。L10/L12低能复S和相移差异仍大，但两组保留基不同；`SOURCE_UV_L10_L12_DRIFT.json`只给描述性漂移，不把全部差异归因于L。暂不扩展更高L扫描；共享基底IR对照、名义tip与相邻ε已完成，当前已转入同容差的完整UV tip/ref。

另已测得旧IR/UV分别计算的SVD保留子空间不完全相同，最大主角正弦≈8.73e-5，见 IR_UV_BASIS_COMPARISON.json。已有对照已附 BASIS_SCOPE_ADDENDUM.json，不把差异全部归因于UV。现已统一SVD未缩放源行顺序，并实现 `--basis-source-report` 显式复用确切保存基底；比较分组核验实际基底身份，无法核验时不允许宣称隔离了UV因素。IR_L12已使用已验收UV_L12的同一基底完成重新装配、求解和独立验收。此修复不证明SVD截断对最优值无影响。

旧输入sine控制rho577.72MeV属于chi-b/SR-b/frozen，不混入新基线。两个超时/中断目录保留为工程失败，均未产生SDPB科学结果。MMA上限为18000秒，同M/同脚本结果可经校验复用。

完整范围和状态：EXECUTION_FROM_RESEARCH.json；下一进程定位：CONTINUATION_STATE.json；首要缺陷：LIMITATIONS_FIRST_ZH.md。

## 自包含运行方式

在仓库根目录使用已安装项目的Python环境；首次安装为 `python -m pip install -e '.[sdpb,test]'`。本机已验证的解释器为 `/home/shiqiu/miniconda3/bin/python`。Docker需可调用，镜像为 `bootstrapcollaboration/sdpb:3.1.0` 和 `wolframresearch/wolframengine:15.0.0`；本机实际镜像ID和包版本保存在 `RUNTIME_ENVIRONMENT.json`。Mathematica许可证目录默认 `/home/shiqiu/Licensing`，也可用 `MMA_LICENSE_DIR` 指定含mathpass的目录。运行材料不包含许可证内容。

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 \
python -m smatrix_bootstrap.run sdp solve \
  --workdir /playpen1/shiqiu/sdpb_mainline_20260912/新的基线目录 \
  --M 50 --L 10 --chiral --chi chi-c \
  --uv --uv-parts gram fesr ff --sr SR-a --ff-factor node \
  --scattering-prescription sine-cardinal --reduce-basis \
  --operator-dps 40 --digits 30 --precision 320 \
  --nproc 8 --duality-gap 1e-6 --timeout 7200
```

默认通过 `scripts/mma/wolfram.sh` 调用Mathematica独立公式检查，SDPB通过 `scripts/sdpb/sdpb.sh` 调用固定3.1.0镜像。不要打印许可证内容。`--points tip ref mid`执行既有截面规则；它不自动解决作者最近边界点的几何身份问题。`sdp support --source-report 完整leaf/report.json --point ref --out 新目录`仅可复用相同源算子和合同，改变SR/FF/chi或源算子后须新建约束。

Python装配可使用4线程；SDPB启动器现单独固定每个MPI进程的BLAS/OMP/MKL为1线程，并显式传入容器，遵循[SDPB 3.1.0 Usage](https://github.com/davidsd/sdpb/blob/3.1.0/docs/Usage.md)。宿主原生分支和实际容器均已验证三个值为1。这不更改已运行进程。本机旧容器每rank曾建立42线程，但短时采样仅主线程消耗CPU，不能把线程数量直接当作持续算力消耗。见 `runtime_thread_fix/`、[资源核验](results/runs/sdpb_mainline_20260912/runtime_resource_review/REPORT_ZH.md)。

同PMP的8rank预处理测量受Docker创建延迟影响，在预算边界才取得CID，随后被所属进程清理；966.845秒、返回124，没有完成矩阵转换，因此没有可报告的加速比。此前单rank预处理外部570.678秒、程序内部105.198秒，差额也不能当作MPI计算时间。当前保留原预处理进程数，不重复该测试、不虚报多节点性能。打印检查遗留的一个从未启动容器也已按确切身份删除，科学容器未受影响。

192bit曾完成部分L10，但L10/L12的观测条件数均可达约1e69，L12已发生192bit失败；当前新主线保留实测成功的320bit，不追求更高精度。本次重跑入口为 `sdp refine --source-report 终态leaf/report.json --out 新目录 --precision 320`，它保留确切PMP输入，不重新选择基底。不要对仍存活的计算调用重跑。非零求解退出现在会保留但隔离可能过时的部分输出，不再用计时阶段的y作物理读数。

源式复验使用 `sdp audit-run --source-report 完整leaf/report.json --out 新目录 --bits 384`。它同时输出原生Gram机制诊断，不添加Watson或幺正化约束。独立对偶和连续幺正证书不是当前有限论文复现的新前置门槛。

本机只确认一个80逻辑CPU、约503GB内存主机；实际8 MPI ranks。没有真实远端调度入口或已分配hostfile，不能声称验证100节点。真实集群须加载兼容SDPB/MPI模块，使用共享目录，再设置 `SDPB_NATIVE=1 SDPB_HOSTFILE=/实际/hostfile`；`--nproc`是总进程数。默认一项科学优化；只为已记录的独立主线问题允许有限并行。已完成的UV双任务合计16ranks，不代表16节点；C5续算使用单个8ranks任务，现已完成。

原生启动器使用OpenMPI语法；应在已分配的计算作业内运行，并使用与SDPB构建匹配的MPI模块和远端可见路径。现在显式用 `-x` 传递三个单线程环境变量。本机OpenMPI4.1.6的两rank localhost hostfile调用已通过；这验证入口和环境传递，不代表测试过远端SDPB。见 `runtime_thread_fix/explicit_mpi_exports_validation.json` 和 [OpenMPI环境导出说明](https://www.open-mpi.org/doc/v4.0/man1/mpirun.1.php)。不要把其中的localhost测试文件当作集群分配。

## 验收与剩余范围

保持C1–C8和Fig.3–11完整范围：纯区域、手征ε几何/线性检查、近物理IR幅度、加入UV后的截面非对称收缩、三代表相移/eta、五组M/L变化。当前未完成；一份可行UV幅度不替代区域或稳定性。

先判断原式数值是否通过，再判断是否解了同一个具名模型，再比较物理。Gram PSD不自动意味着弹性和两pion谱密度均饱和；`accuracy.current_mechanism`分开记录这些量。不要把2505 Watson迭代或2403的B正则化当2309必需步骤。

按现有M50实测，冷启动单个支持约45–70分钟；同合同邻点warm-start实测12.49分钟。固定M50/L12的一轮IR两端＋UV三代表，按两组各8ranks并行及已验收邻点复用，约2–3小时（理想实测工段合计约2小时）；冷启动串行或环境重试约3–5小时。新名义tip的转换外部1688秒、内部75.443秒，SDPB外部2089秒；环境启动等待明显，失败重试会延长。完整Fig.3–11另需区域和分辨率计算，不能承诺同样工期。不要只为缩小已远低于物理容差的数值误差而提高精度。

最新完整工程测试为 [续算回归日志](results/runs/sdpb_mainline_20260912/checkpoint_resume_control_20260913/regression.log)：217项通过，75.19秒。精确输入/失败输出、共享基底、任意方向复用、几何选点、完整阈下求值及参考覆盖保护均有独立控制材料。早期证据不删除；工程通过不能称为物理完成。

## 原文验收修正、参数查找与完整振幅复用

原C4把S0在0.9GeV限制于85–110°，却会拒绝论文Fig.7自己的约82.57°。新的 `C4_SOURCE_ACCEPTANCE_AMENDMENT.json` 在读取新下点相移前冻结：与真实Fig.7 IR曲线比较，各波RMS和最后公开能量点差均不超过10°。这是项目的比较容差，不是论文统计误差。旧预登记及旧C4仍保留；新规则下当前IR仍失败。

旧C7的1.196GeV略超Fig.10实际最后点约1.19574GeV，导致参考曲线自身也被拒绝。`C7_source` 使用公开数据域内末点，保持原S0/S2区间、10°RMS和三代表完整性要求。若待比较幅度没有覆盖公开末点，返回未完成，不能缩短比较窗获得通过。

C3的8%是相同ε的Fig.5彩色曲线复现误差，并非相对Weinberg直线的理论容差。完整阈下求值已覆盖全部Fig.5横坐标（共131个s点），禁止插值器默默钳位到旧3.95终点。`sdp subthreshold --source-report leaf/report.json --out 新目录`只重评已保存幅度；`sdp figures --out 科学运行根目录 --subthreshold-replay 求值/report.json`绑定身份后绘图，不改原叶。参考纵坐标不进入优化。详见 `paper_acceptance_self_audit_20260913/`、`subthreshold_completion_20260913/`。

作者公开仓库README明确2309当时未发布代码；已检查版本源包、仓库历史、作者页面、会议材料和后续代码，尚未找到2309 producer。见 [搜索报告](results/runs/sdpb_mainline_20260912/author_implementation_search_20260913/REPORT_ZH.md)。后续2403的手征代码使用bare f、无κ或Λ隐藏权重，但这不证明2309用了同一范数；不能把后续B或Watson回填。当前仍有截断和近物理选点问题，不能宣称只剩作者缺参。相邻ε试验是有论文来源的有限尝试，不是按相移拟合输入。

完整y与保存binary64基底现可通过精确整数点积恢复全部3876个原密度系数：`sdp amplitude-export --source-report leaf/report.json --out 新目录`。产物包括纯数据 `amplitude.wl`、JSON分块、完整y和哈希，可用Mathematica `Get`读取；拒绝未终态、超时残留和缺少求解时y哈希的旧结果。相对保存输入的转换误差严格为0，不等于提高物理解精度。已完成L12 IR/UV两份导出，约7秒/份，见 `canonical_amplitude_export_20260913/README_ZH.md`。Mathematica15已实际Get完整UV导出，精确整数/有理数结构、ρ2对称与两项独立角积分均通过（56秒）；见 `canonical_amplitude_export_20260913/mma_readin_20260913/report.json`。这些材料仅本地整理，未向作者外发。

旧 `report.py` / `scripts/sdp/finalize.py` 对空输入也会生成已被否定的B、精确SVD及M50失败结论。已将两入口关闭，在任何写入前明确拒绝；原始producer字节、SHA及空输入反例保存在 `legacy_report_retirement_20260913/`，7项控制通过。当前figures/contrast链不受影响，不能再调用旧finalize汇总SDPB结果。


## 最新执行优先级

先完成底层推导与实现自洽检查，再讨论作者未公开约定。MMA核验和sine独立积分核验上限统一18000秒。若最终只有未公开参数阻碍精确对应，按最有来源依据的选择完成有限模型计算，打包已复现结果与代码，把缺陷放在最前面，供课题组向作者请教；不擅自外发。

此前`UV_M50_source_ref_reuse_mma`在交互中断时终止于PMP写出，残缺输入和interruption.json已保留，未启动SDPB。已改为独立进程运行同一新基线，避免交互中断导致重复装配。当前物理输入未因该工程问题改变。

## 底层独立验证的新证据

[归一化与手征工作流复核](results/runs/sdpb_mainline_20260912/foundations_normalization_review/REPORT_ZH.md)（[三页PDF](results/runs/sdpb_mainline_20260912/foundations_normalization_review/pdf/REPORT_ZH.pdf)）直接对照v3 PDF页码并独立重算冻结IR幅度，未发现整体符号、1/4投影、8π²或κ因子错误。正确的Weinberg黑点为(0.07332139057,−0.00488809270)；当前IR ref的第二坐标却为−0.00169984366。其阈值f00约0.3701，Weinberg值约0.10265，且S2阈值符号也不同；单个归一化因子无法修复这些形状差异。独立原始角积分也重现S0在0.40235GeV的71.0055°。

论文先检查阈下线性，再使用近黑点边界幅度讨论IR相移。当前固定xref上端点是声明代理，这一步的近线性与最近点身份尚未成立。现有相移失败不能当作论文最近边界点的反例，也不证明该代理绝不是最近点。L12截断检查和两截面端点的几何选择已完成；仍不按相移另挑代表，全边界最近点及插值稳定性未被证明。

`foundations_independent_control/report.json`：六个独立mpmath原始角积分/打包控制通过。`foundations_M50_angle_check/report.json`：完整M50系数的S0 node42、P1 node33和I1 ell19 node0三点通过，16.34秒；是独立数值积分对照，不是连续证书。入口：`sdp angle-audit --source-report 完成leaf/report.json --source-snapshot 对应启动源码.json.gz --out 新目录`。

`currents`按完整y重建，仅核验两主波和电流：`sdp currents --source-report 完成leaf/report.json --source-snapshot 对应启动源码.json.gz --out 新目录`。已完成旧输入sine控制的100节点诊断，耗时5.83秒；P1两pion谱部分仅占两FESR矩约7.5%和5.1%，虽然弹性接近饱和。这是该振幅中UV传递未闭合的直接量化证据，不能自动归为精度不足或作者参数缺失。较松诊断积分使6个极小Gram子式符号未定，不推翻原700项更紧核验的通过。见SINE_CONTROL_MECHANISM.json。

尚需：SVD对主线支持/相位稳定性的检查；physical-off-node虚部尾界已修复并通过独立解析反例；新源文基线与完整IR/UV、C1–C8链。不要把这批独立检查说成全目标完成。

## 当前研究报告与便携材料

[综合报告](results/runs/sdpb_mainline_20260912/reproduction_review_20260913/REPORT_ZH.md)及[PDF](results/runs/sdpb_mainline_20260912/reproduction_review_20260913/pdf/REPORT_ZH.pdf)纳入名义tip、相邻χ、凸组合校准、两IR端点和三UV代表，缺陷优先。`reproduction_package_20260913/`持续更新完成案例；以 `PACKAGE_INDEX.json` 的实际归档哈希与验证为准，不能把staging新增文件说成旧压缩包已有。完整PMP/basis单独压缩，mask、MMA复用资料和完整y保留；搬移后实际PMP目标重建、MMA资料复用及幅度重导出已通过。许可证和账号凭据不入包。

## 选定IR与近黑点UV代表的统一交付接口

`python -m smatrix_bootstrap.run sdp ir-figures --source-report 选定IR叶 --ir-selection-report 几何收据 --out 新目录`会重新核对几何身份，输出相移、η、阈下图及完整支持CSV。它不借用其它容差的UV幅度，也不把单点线性图叫完整C3。实际CLI已经验收，数据逐项与保存值相同，论文η列为空。

旧mid=(xref+xtip)/2在当前名义tip很远时不符合“黑点附近”。旧规则和登记保留；新 `uv-plan --source-report 完成UVtip叶 --out 新计划` 冻结 x_near=xref+min(.01*xref,(xtip−xref)/2)。1%只是数值横坐标邻点上限，完整二维距离单列，不是作者坐标、物理误差条或近点合格半径。随后用通用support分别求ref及该固定near截面，再调用 `uv-select --source-report 计划/report.json --ref-report 完成ref叶 --near-report 完成near叶 --out 新收据`。

`contrast --uv-report 收据/report.json`自动识别并重新核验新收据；结合既有 `--ir-selection-report` 比较同一选定IR。内部mid槽显式映射near，图例写UV near，旧登记完成标记为false。它核验完整ModelSpec、实际basis、有效mask、完整y、真实PMP目标/截面/物理块及producer身份。36项独立控制和28项CLI/交付检查通过，实际名义tip计划已核验；真实三UV收据和完整contrast已通过实际CLI核验，物理C4/C6/C7仍FAIL。证据在 `uv_selection_validation_20260913/`。

prepared复用现在按ModelSpec核对disk_mask.npy的完整内容：无显式mask时先全true；有UV Gram时两主波对应盘置false，由Gram承担其幺正条件。不能把UV mask补成全true；M50/L10 UV是1500形状、1400个true。错误mask提前拒绝，正常源投影/目标字符串不变。26项独立控制、17项定向检查和217项完整回归通过；历史错误合成fixture保留，新版本36+28项选择/交付控制均通过。见 `prepared_mask_identity_fix_20260913/` 和 `uv_selection_validation_20260913/mask_guard_compatible/`。

便携依赖审查确认旧11-case快照曾漏带mask，已补全部13个已完成科学叶的原mask及11组M50 MMA report/independent资料。实际搬移后Pmp.from_saved和uv-plan通过，目标与原PMP头逐项相同，MMA结果可校验复用；不含许可证。新UV选择接口是后续版本升级，不追认为旧快照已有。见 `reproduction_package_20260913/PACKAGE_DEPENDENCY_VALIDATION_ZH.md`。


## 已实测的SDPB检查点复用

`python -m smatrix_bootstrap.run sdp support --source-report 完成UVref叶 --direction 0 1 --fix-f00 0.07405460447866398 --warm-start --timeout 1200 --out 新目录`可复现本次near初始化方式。source报告必须已接受；完整转换块结构/采样、basis、normalization、SDPB3.1.0、320bit、8ranks和分块分配必须匹配。只允许目标/截面改变，原式验收重新执行。检查点独立复制到warm_start_input，输入与输出sdp.ck分离，原checkpoint不改。

near的SDPB日志实际记录加载binary checkpoint，主迭代14步，SDPB520.36秒、转换121.61秒、总749.37秒，最终primalError4.72e−32、dualError3.51e−11，原式数值通过。它说明邻近问题可减少耗时；没有同PMP冷启动对照，不宣称严格加速倍数。当前仅验证单宿主Docker；多节点须保持兼容MPI及分块身份，不能把本地加载说成百节点验证。原始初始化检查点及哈希留在运行目录，便携归档是否含其二进制以PACKAGE_INDEX为准。


同输入续算入口为 `python -m smatrix_bootstrap.run sdp refine --source-report 终态未验收叶 --precision 320 --resume --timeout 7200 --out 新目录`。source必须有终态进程记录和有效检查点；PMP、目标、截面、ModelSpec及除timeout外的Settings完全相同，源进程不能仍活跃。不同精度须使用普通refine新求解，不能二进制续算。30项checkpoint控制、9项CLI控制及217项回归通过；实际加载由新日志单列标记，不把复制完文件当成恢复成功。C5续算已实际通过加载与最终原式验收。

[剩余主线决定](results/runs/sdpb_mainline_20260912/remaining_mainline_decision_20260913/REPORT_ZH.md)：已存χ=.001的大x见证，在ε=.002下满足三种常见未加权范数；扩大包含旧U的SVD空间仍会保留该允许点。因此不能靠这些改法缩回过宽区域。当前物理fπ校准已足以交付本轮条件性失败，minε不是论文要求的新前置门槛。C1/C2/C3/C8的完整缺项继续保留；下一求解先明确截断辨识作用域，不按相移继续调χ。
