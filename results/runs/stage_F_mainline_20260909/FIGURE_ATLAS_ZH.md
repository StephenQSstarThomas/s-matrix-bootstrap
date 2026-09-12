# Fig.3–11：主线复现图谱与claim判定

本图谱将2309.12402v3的九张结果图对应到A–F数据、完整振幅及已成立/仍有差异的结论。**F1五组完整支持、F2原式/新原生审核及F3图谱已交付；原图定量差异仍保留。** 原论文见[本地PDF](../../../references/2309.12402v3.pdf)，科学定义见[SCIENCE](../../../SCIENCE.md)。

## 共用科学底座：A与D

A固定同一完整幅度、crossing、分波归一化和单位；当前图示主线采用PV/Legendre-Q、M50/L10、sampled/free、3876幅度变量及actual-ρ组合L4上限377500。L按每个isospin分波数计数，原生幺正盘共1500个。历史finite-sine族与当前PV不能交换可行性/秩/离节点结论。[生产算子](../stage_B_pv_M50_L10_prepare_20260907/report.json)与[独立投影/Gram/单位审核](../stage_E_closure_audit_20260909/OPERATOR_AUDIT_ZH.md)提供底层入口；后者在所审查项目中未发现重大归一化错误，不是排除所有程序错误的证明。

D为Fig.8–10补入共同S0/P1电流Gram、F0(0)≈1/F1(0)=1、四FESR及高能FF平方界。输入在求解前固定为clipped-phi、打印矩、四raw绝对误差±.002、epsilonFF=6e−5、1.2GeV匹配尺度。[D报告](../stage_D_mainline_20260908/D_RESULT_ZH.md)、[电流算子](../stage_D_mainline_20260908/D1_current_M50/report.json)及[4076变量见证](../stage_D_mainline_20260908/D3_hull_009/joint.npz)可直接追踪。D见证只证明联立可行，不代替E的全幅度支持或相移结果。

## Fig.3 / B2：基本散射条件允许多大的区域

**物理问题：** 只用解析性、crossing及幺正性，投影到x=f00(3)、y=f11(3)后允许哪些值；这一基线并不能识别QCD。

**实际入口：** [共用区域PDF左面板](../stage_B_delivery_20260907/regions/regions.pdf)、[完整数据与支持索引](../stage_B_delivery_20260907/regions/regions.json)、[B报告](../stage_B_delivery_20260907/B_RESULT_ZH.md)。

**已成立：** 完整pure平面内外包已交付，38份有效支持、36个法向；本项目规定度量下的几何距离上界.0096758≤.01。

**差异与范围：** 原图显示轮廓仍有局部约.016的差别；本方几何误差不是与论文的误差，也不是连续能量/无限分波证明。论文§4.1将pure图作为独立基线，不再用它直接生成后续QCD结果。

## Fig.4 / B3：手征信息怎样把区域压成窄带

**物理问题：** 六个χ容差怎样限制区域、形成Weinberg关系y=−x/15附近的薄带，并保留物理fπ参考点。

**实际入口：** [区域PDF中/右面板](../stage_B_delivery_20260907/regions/regions.pdf)、[六域数据](../stage_B_delivery_20260907/regions/regions.json)、[Fig.4方法审计](../stage_C_mainline_20260908/FIG4_AUDIT_ZH.md)、[参考点判定](../stage_B_delivery_20260907/REFERENCE_POINTS_ZH.md)。

**已成立：** ε=.006/.004/.002/.001/.0006/.0002共用同一H、B与两个四维L2球；六个x≥0窗口均满足既定.01几何标准。嵌套与窄带关系成立；.002及.001包含物理参考点，.0006/.0002排除。

**差异与范围：** .002的本方右端约.1059，原图最右显示点约.0826；更小ε也保留明显差异。不能把作者显示点当作其真实支持上界，也不能称.002已被证明唯一最优。两套χ范数组合方式的诊断没有追认成此图的新结果。

## Fig.5 / C1：四个χ点是否带来合理的阈下形状

**物理问题：** 区域靠近Weinberg线是否也对应0<s<4内近线性的S0/S2/P1，尤其S0的手征零点是否保留。

**实际入口：** [三波图](../stage_C_mainline_20260908/C1_fig5/profiles.pdf)、[曲线数据](../stage_C_mainline_20260908/C1_fig5/profiles.json)、[C1报告](../stage_C_mainline_20260908/C1_RESULT_ZH.md)。

**已成立：** 同一预定物理xref、ε=.006/.004/.002的三代表振幅在81个阈下点比较；缩小ε后偏离下降。蓝色曲线所用网格未见S0零点，橙/绿零点分别位于.30–.35/.40–.45，向Weinberg的.5靠近。

**差异与范围：** 这是预选三个完整振幅的形状比较，不是所有可行振幅的线性定理，也不把四点χ约束当成整个区间精确线性。零点括区是当前网格读数。

## Fig.6 / C2：从几何点锁定同一完整振幅

**物理问题：** 选取物理手征参考点附近的上边界振幅，为下一张物理相移图固定输入。

**实际入口：** [代表点图](../stage_C_mainline_20260908/C2_fig6_display/profiles.pdf)、[全部3876系数](../stage_C_mainline_20260908/C2_fig6/coefficients.json)、[选点记录](../stage_C_mainline_20260908/C2_fig6/selection.json)、[C2报告](../stage_C_mainline_20260908/C2_RESULT_ZH.md)。

**已成立：** 原样复用Fig.5绿色振幅，先选点再看相移；论文图示magenta/绿色的对应已审查。

**差异与范围：** 本方代表振幅是已验证内包边上的凸组合，原盘、χ及L4都验回；不是作者唯一精确极值。凸混合可增加非弹性，二维区域精度不能转换为相移精度或极值唯一性。

## Fig.7 / C3：仅靠IR能说明哪些相移

**物理问题：** 用同一个χ-only振幅，比较S0/S2的低能形态与P1的rho不足，建立加入UV前的对照。

**实际入口：** [相移图](../stage_C_mainline_20260908/C3_fig7_final/profiles.pdf)、[δ/η数据](../stage_C_mainline_20260908/C3_fig7_final/profiles.json)、[原生求值](../stage_C_mainline_20260908/C3_evaluate/evaluation.json)、[C3报告](../stage_C_mainline_20260908/C3_RESULT_ZH.md)。

**已成立：** 三波43个原生点加阈值的132项所列检查通过；S0/S2低能形态较接近参考，P1在约.792GeV仅约5.29°，没有正确rho形态。

**差异与范围：** S0中/高能仍有可见差异；此结果不排除整个χ可行集中存在rho。主相位为unwrap(arg S)/2，η另存；没有强制η=1或使用arg(f)替代。灰线/数据来自PDF读数，误差条和协方差没有完整导入，不是拟合数据集。

## Fig.8 / E1：QCD约束是否收紧区域及上边界

**物理问题：** 除共同散射/χ条件外，增加电流PSD、FESR、FF高能抑制后区域怎样变化。

**实际入口：** [本轮区域图](../stage_E_closure_20260909/E1_regions_certified/regions.pdf)、[7份联合点/支持及Arb截面](../stage_E_closure_20260909/E1_regions_certified/regions.json)、[E总报告](../stage_E_closure_20260909/E_RESULT_ZH.md)。

**已成立：** 全4076变量联合+x支持位于[.09924806,.09925920]，严格小于同设置IR右端约.1059。xref处上侧变化严格大于下侧，比值在[1.026389,1.167740]。

**差异与范围：** 原图标记凸包的同位置比值约5.35，强非对称程度未复现。当前UV是稀疏内外包，局部Arb截面结论不等于密集全轮廓收敛。改变后续Watsonian目标不会改变这张图的数学可行域。

## Fig.9 / E2–E3：UV能否产生稳健的P1/rho

**物理问题：** 比较tip及两个预选附近完整解的P1相移，检查rho的90°上穿和三点稳健性。

**实际入口：** [最新P1图](../stage_E_closure_20260909/E3_W_phases/fig9.pdf)、[全部曲线/η](../stage_E_closure_20260909/E3_W_phases/phases.json)、[原三点身份](../stage_E_mainline_20260908/E2_representatives/selection.json)。后续完整解分别为[tip](../stage_E_closure_20260909/literal_W_tip_path_03/report.json)、[mid](../stage_E_closure_20260909/literal_W_mid_path_03/report.json)、[ref](../stage_E_closure_20260909/literal_W_ref_path_02/report.json)。

**已成立：** 三个预定角色各完成一次2403三主波Watsonian固定目标求解，全部原约束保留且通过，主nearest-S相位均出现90°上穿；tip/mid/ref线性读数约807.57/728.39/697.60MeV。

**差异与范围：** 续算是具名后续方法，不是2309明确记录的步骤；它释放原投影等式，因此新曲线不能冒称仍属于旧边界点。三点读数仍跨约110MeV，原图三条曲线的接近程度未复现；ref相对770MeV约−9.4%。90°读数是原生节点间插值，不是复平面极点质量。旧未续算结果（tip约812.36、ref约699.96、mid未见上穿）仍保留在[原相移数据](../stage_E_mainline_20260908/E3_phases/phases.json)。

## Fig.10 / E3：同样三套解的S0/S2及非弹性

**物理问题：** UV加入后S0/S2是否仍主要由低能约束决定，以及三代表点在哪些能区开始分散。

**实际入口：** [S0/S2图](../stage_E_closure_20260909/E3_W_phases/fig10.pdf)、[η图](../stage_E_closure_20260909/E3_W_phases/inelasticities.pdf)、[旧/新投影及谱诊断](../stage_E_closure_20260909/E3_W_phases/watson_diagnostics.pdf)，完整解与数据同Fig.9。

**已成立：** 全程使用Fig.9相同三套4076变量解，没有为S0/S2换点；396项三解原生主波检查通过。低能较接近、高能更分散的趋势保留。

**差异与范围：** 全能窗三曲线范围约S0 57.8°、S2 26.9°，定量一致性未全部复现。一次旧Q目标收敛不是新S/新F的Watson固定点证明。η接近1但两π谱份额很小时，仍可存在明显Watson相位偏离；[FF引导模π图](../stage_E_closure_20260909/E3_W_phases/ff_guided_phases.pdf)只作预先定义的诊断，未用于替换主图或补出PV散射节点间绕行。

## Fig.11 / F1：五组分辨率下结果是否稳定

**物理问题：** 对完整UV流程比较(50,8)、(50,10)、(50,12)、(45,10)、(60,10)；作者认为S0/S2较稳，P1保留rho但位置随M略移。

**实际入口与当前状态：** [F报告](F_RESULT_ZH.md)、[Fig.11](fig11/fig11.pdf)、[η图](fig11/fig11_eta.pdf)、[完整数据与审核索引](fig11/resolution.json)。五组全空间UV支持均达到原式gap≤1e−4，669项新原生主波检查通过。

**已有来源材料：** [附录A/Fig.11审计](SOURCE_AUDIT_ZH.md)、[S0原图CSV](../../../references/figure11_s0_phases.csv)、[S2原图CSV](../../../references/figure11_s2_phases.csv)、[P1原图CSV](../../../references/figure11_p1_phases.csv)及各自metadata。每份CSV218点，是五条印刷图曲线的矢量标记读数，不是作者全部振幅系数。

**已成立与差异：** 五组均有P1的90°上穿，读数约840.6/812.4/798.3/752.1/773.8MeV；S2较稳、低能S0接近。固定M的读数跨度42.3MeV，原图约12.5；固定L跨度60.2MeV，原图约54.3。L依赖强度未完全重现，中高能S0绝对相移仍偏离原图。原附录没有给唯一代表身份、M变化的regulator缩放或百分比门槛，本方保持预先声明的+x规则和B(M)配方，未按相移挑点。图例颜色只对应M/L，不提供Fig.8点身份。

## 统一解读

已证明的有限支持界、预选振幅的相移现象、原图显示读数及连续物理命题分别记载。通过本方计算验收不等于论文全部claims通过；原图定量差异也不能仅凭图形比较升级成论文理论错误。

当前A–F条件计算链和逐图证据均已交付。F没有抹平E的定量差异，也没有新增物理参数、事后择点、额外M/L扫描或连续认证前置门槛。
