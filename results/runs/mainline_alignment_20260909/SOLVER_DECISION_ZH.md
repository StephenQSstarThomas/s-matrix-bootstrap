# 求解器选择、实测依据与有限原型的严谨性

2026-09-09。只读核对本地2309v3源、官方后续代码和既有实验记录；未增加求解、web研究或archive扫描，未改src。

**结论：不能把现状描述成“舍弃了2309指定的求解器，因为我们开发了更优方法”。2309v3正文没有指定具体solver；后续程序明确CVX/MOSEK。仓库保留Newton的依据是本问题上的可验证产出和具体数值修复，不是对MOSEK的全面优越性证明。**

## 1. 2309与后续代码究竟指定了什么

- v3 TeX931描述凸集合及最大化线性泛函，并举径向最大化t的例子；998–1017给PSD约束。数值方法正文没有指定Newton、CVX、MOSEK、内点参数或solver容差。
- 全源中CVX条目位于TeX1576–1579，但该书目块以`%%`注释，不能据此认定实际指定了CVX，更不能据此认定底层求解器。
- 作者README81–100明确：2309是原始提案，当时没有公开代码；2403是首次公开程序。该声明不证明作者未使用某个solver，只说明没有2309原producer可直接照搬的公开证据。
- 后续2403 `GTB_numerics.m:51–52`、2505 `optimize_core.m:28–29` 明确写 `cvx_begin sdp` 和 `cvx_solver mosek`。

**CVX是MATLAB凸优化建模系统，负责表达和转换约束；MOSEK是这些脚本指定的底层数值求解器。** “用CVX”与“用哪个solver”不是同一层问题。我们采用的障碍函数、Newton步、QR/CG及迭代精化也属于凸优化数值实施，不是新增物理原理，不能仅凭算法名称宣称比成熟conic solver更先进。

官方文档进一步说明：[MOSEK conic求解](https://docs.mosek.com/latest/pythonapi/solving-conic.html)采用homogeneous/self-dual interior-point框架，终止同时考察primal/dual residual及gap；[Clarabel](https://clarabel.org/stable/)也明确属于interior-point solver。Newton步骤与成熟锥优化并不互斥。但本仓库的可行中心／显式Phase I路线不是MOSEK的同质自对偶实现，不能把名称相近当成算法等同，更不能当成效率优越证据。

来源：[v3本地TeX](../../../references/2309.12402v3-source/prd_submission_2.tex)、[作者README](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/README.md#L81-L100)、[2403源码](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L51-L65)、[2505源码](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m#L27-L38)。

后续程序同时改变了中心ν0、UV匹配尺度、电流数、谱表示、FESR和部分范数／选点方法；把整个后续输入包搬来并不自动复现2309。后续核的100位十进制预计算也不是优化器以100位运行的证据。准确版本及边界见 [FOLLOWUP_METHODS_ZH.md](../stage_B1_source_method_audit_20260907/FOLLOWUP_METHODS_ZH.md)。

## 2. 已有实验能支持什么选择，不能支持什么排名

| 已有事实 | 可以得出的结论 | 不能得出的结论 |
|---|---|---|
| D3_joint_001/002在同一PV/M50/L10、B377500、χ=.002、clipped-current模型上，直接全空间Clarabel约316.90/320.36秒后MaxTime；原式均无有效lower。003在第1步NumericalError | 当时这些具体矩阵／数值表示没有交出联合可行支持解 | Clarabel普遍不适用；原问题不可行；MOSEK也会失败 |
| D3_hull_009仍用Clarabel，23步约.263秒，在全部204个同模型父幅度凸包内联立电流／UV，得到完整4076变量原式通过的见证 | 任务收束为存在性并复用数据有效；成熟solver在合适子问题上仍有实际价值 | 390维受限搜索比全空间优化快，就证明任何solver全面更优；该凸包可替代E1全空间外界 |
| B最窄条带的Newton停止修复后，原式区间为[.7842651766943,.7847169067406]，gap≈4.52e−4；含监督／核验130.69秒 | 该同题实例中的具体错误修复有可验证效果 | 所有方向有相同速度；Newton不存在其它误差 |
| B_SOLVER_FIX发现复用QR时CG真实相对残差最高约.04027，旧逻辑仍使用方向；修正为残差不达标就刷新，fresh QR仍失败则退出 | 缺少线性解质量门槛是我们实现中的问题，不能用“步数用完”代替解准 | 新门槛本身就是原式最优性或相移精度证明 |
| NEWTON_REAUDIT的另一MID运行仍在慢中心化，但13步的原／变换方向减量差≤2.47e−10、曲率差≤6.74e−9、线性相对残差≤1.35e−7 | 没有该方向的明确梯度／Hessian错配证据；α接近2本身不是过冲证明 | 慢中心化已解决，或所有小特征方向都足够准确 |

具体记录：[D全部尝试](../stage_D_mainline_20260908/ATTEMPTS.json)、[D3结果](../stage_D_mainline_20260908/D_RESULT_ZH.md)、[B同题原式结果](../stage_B3_pv_M50_eps0002_central_fixed_20260907/report.json)、[B_SOLVER_FIX](B_SOLVER_FIX_ZH.md)、[NEWTON_REAUDIT](../claims_alignment_20260909/NEWTON_REAUDIT_ZH.md)。

因此，当前选择的合理表述是：**继续维护已经在本仓库声明问题上产出原式可行点／支持界的Newton路径，移除长期没有收益的具体分支以减少维护与上下文负担，同时保留Clarabel在成功初始化和数学控制中的用途。** `linear.py:237–249` 也明确区分联合hull初始化与全空间centered/barrier支持。B交付记录还说明未改善M50的CVXOPT分支已移除；CVXOPT同样不是MOSEK。

本次所读记录没有“同一个冻结模型、同一个目标、同一个验收精度、同一硬件”下MOSEK与本Newton实现的公平性能对照。因此不能声称作者solver被实验证伪、Newton全面更快、更准，或当初的实现路线是最优决策。此前大量未收敛尝试和后来查出的Newton问题都应如实承认，而不能全部归咎于论文或成熟solver。

历史 [OUR_METHOD_REVIEW](../stage_B1_source_method_audit_20260907/OUR_METHOD_REVIEW_ZH.md) 还指出，早期曾把finite-sine、T0zero、五尾、global-FG及扩展采样作为唯一主线。这既改变问题也改变计算负担；这些更强问题的性能不能拿来比较作者公开的有限设置。该报告是历史偏离与修正记录，不是当前已完成D/E/F仍未实现的状态说明。

## 3. 等价数值变换与改变科学问题的分界

| 类别 | 例子 | 所需结论范围 |
|---|---|---|
| 同一有限问题的数值表示 | 可逆subtracted/absorptive坐标、正单位缩放、PSD正对角合同／实化、准确消元辅助变量、QR/CG预条件 | 全部原物理方向与约束须保留，输出逆变换后仍需原式检查；浮点实现可能有误差 |
| 求解路径选择 | 正的障碍权重、可行线搜索、中心续算、足够准确的Newton方向 | 不改变最终可行集合，但不同路径／目标面可能选择不同完整幅度，需保存代表点规则 |
| 仅用于产生候选的子问题 | D3父幅度凸包、Phase I暂时松弛UV | 只有原问题重新通过才是原见证；受限最优值不是全空间支持外界，τ>0不是原联合可行 |
| 改变科学或离散模型 | χ两球改八维球、改变B或ρ度量、改FESR误差尺度／截止、PV改finite-sine、增加T0zero/FG/新采样点 | 必须具名、重建相应身份并复核；不能称作单纯换solver |

例如H=q+b换元的T0、σ必须一起变换；设T0_sub=0不是这项恒等式。密度SOC辅助变量消元保留所有actual双密度方向；不能把实际减小M或截掉方向也称“消元”。这些界线详见 [定义重审](../claims_alignment_20260909/DEFINITION_REAUDIT_ZH.md)，以及 [截止规则决策](CUTOFF_DECISION_ZH.md)。

## 4. 严谨性来自分层验收，不能由solver标签保证

1. **定义及接线。** 明确源版本、H/K/k、C_flat/actual-ρ、χ范数、FESR与FF单位。独立从论文重建角投影、复Gram和矩关系，避免仅比较两个共享错误的实现。该层仍有未公开原始设置，不能假称身份完全恢复。
2. **原式primal。** 将保存候选恢复为完整原变量，检查声明的全部有限幺正盘、χ、actual-L4、current Gram全部主子式、四矩及FF。NumericalError/AlmostSolved/小内部残差均不自动决定通过；D3即以原式验界接受AlmostSolved返回的可行点。
3. **原问题dual。** PSD乘子须合法，回传原系数后消去无界自由残差，对剩余密度／电流残差给保守支持界；把可行目标作为lower、有效支持界作为upper。间隙约束的是这个具名有限问题的指定线性目标。内部Newton减量和CG残差不能取代它。这层原式验证也可接成熟solver的候选／对偶，不是自研Newton独有的严谨性优势。
4. **观测量及离散化。** 支持gap小不等于每个系数、相移或ρ位置有同样误差；投影最优也不保证完整幅度唯一。代表点先锁定，随后比较S、η、相位分支和规定M/L配置。保存float64算子的Arb证书不含未知谱插值误差、连续域误差或未指定输入的不确定性。
5. **失败与主张。** 原型计算链跑通、有限约束/目标获证、原图定量吻合、连续幺正证明是不同结论。失败输入及生产源保留，未通过的科学claims继续明确报告，不能用测试数或高bits抵消差异。

证据入口：[独立算子审计](../stage_E_closure_audit_20260909/OPERATOR_AUDIT_ZH.md)、[当前STATUS](../../../STATUS.md)、[D原式核验](../stage_D_mainline_20260908/D3_original_audit_768/joint_audit.json)。这些提高结果的可核查性，但不是“保证没有任何程序错误”的承诺。

## 可直接回答用户的简洁结论

2309v3正文没有指定具体solver，明确采用CVX/MOSEK的是后续公开代码；CVX是建模层，MOSEK才是底层求解器。我们没有同题MOSEK对照，因此不能宣称自研Newton更优。保留Newton是因为它在当前声明模型上已有独立原式可行性和支持界证据；失败的具体native分支已退出，成功的Clarabel初始化仍保留。严谨性依靠原变量约束、合法对偶与误差界，而非solver名称；这些有限证书也不意味着所有原图定量结论或连续物理主张已通过。若此前表述成“2309指定的solver被证明不如新方法”，应纠正。

可附来源：[2309v3数值方法](https://arxiv.org/pdf/2309.12402v3#page=21)、[作者README说明](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/README.md#L81-L100)、[官方CVX/MOSEK调用](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L51-L65)、[本仓库尝试记录](../stage_D_mainline_20260908/ATTEMPTS.json)。
