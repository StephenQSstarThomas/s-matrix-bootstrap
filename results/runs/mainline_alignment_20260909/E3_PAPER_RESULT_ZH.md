# 本轮E3：combined-L2与hard-midpoint的三代表振幅

**三套已冻结完整解均通过768-bit重新生成的原生主波求值，三点P1均出现90度上穿；原论文三点rho位置接近及高能S0形态仍未全部复现。** 本次只有选点记录、求值、绘图和原图比较，没有新增优化、Watsonian步骤、物理参数变化或相位择优。

## 输入与执行

设置由[PAPER_MAINLINE](PAPER_MAINLINE.json)在新相位之前冻结：PV/M50/L10、combined-l2 χ=.002、B377500、free T0、原打印/raw四FESR误差.002、hard-midpoint、FF容差6e−5及算术平均mq。三套均保存全部4076变量，原盘/Gram/χ/L4/FESR/FF在支持报告中通过。

| 角色 | 完整支持报告与向量目录 | 选点记录 |
|---|---|---|
| tip | [paper_tip_path_02](paper_tip_path_02/report.json) | [selection](paper_tip_path_02/selection.json) |
| ref | [paper_ref_path_02](paper_ref_path_02/report.json) | [selection](paper_ref_path_02/selection.json) |
| mid | [paper_mid_path_03](paper_mid_path_03/report.json) | [selection](paper_mid_path_03/selection.json) |

mid选点记录在读取任何新相位之前写入：请求x=.07745128711688606，保存点原H的x约.07745128714200201、y约−.00462508726922234。有效支持方向为(36.32888889008706,500)，原式区间[.5011755703648473,.5019460976605346]，gap=.0007705273≤.001；对应保存点实际x的上边距离≤1.5411e−6。微小x数值偏差已记录，未据相位重新选点。

通过统一CLI依次执行[tip求值](E3_paper_tip_evaluate/evaluation.json)、[ref求值](E3_paper_ref_evaluate/evaluation.json)、[mid求值](E3_paper_mid_evaluate/evaluation.json)，每套使用与历史M50/F相同的43个物理原生节点加精确阈值、768 bits、S0/S2/P1三主波。三套各132项、共396项所列幺正检查全部passed。随后单进程完成gauge-phases及9组compare；原式源/输入未发生变化。 主线另对冻结tip执行了独立的 [768-bit全原生重建](E_tip_fresh_all_native/report.json)：50个节点、全部30个列定分波、共1500个盘全部sampled_passed，最小幺正裕量约3.6082e−85；范围仍为已声明原生采样，不是连续能量/无限自旋证书。

## 三点P1及eta

[Fig.9 P1](E3_paper_phases/fig9.pdf)、[Fig.10 S0/S2](E3_paper_phases/fig10.pdf)、[eta](E3_paper_phases/inelasticities.pdf)、[完整相位数据](E3_paper_phases/phases.json)均已生成。主相位固定为threshold-anchored unwrap(arg S)/2，eta=abs(S)，未加pi改善曲线。

| 角色 | P1原生90度括区/MeV | 线性读数/MeV | 相对770MeV | P1最小eta | S0最小eta |
|---|---|---:|---:|---:|---:|
| tip | [731.675, 792.136] | 788.621 | +2.418% | 0.630996 | 0.898599 |
| mid | [680.414, 731.675] | 700.403 | -9.039% | 0.387090 | 0.206946 |
| ref | [680.414, 731.675] | 696.267 | -9.576% | 0.237792 | 0.572231 |

三点读数跨度约92.35MeV，原图三条曲线的对应读数约812.5–826.6MeV。tip较接近真实rho尺度，但mid/ref早约9%–10%；不能只报告三条都有上穿就称稳健性claim已闭合。所有读数只是已列节点间的线性相位读数，不是复平面极点质量或连续绕行证明。

P1三点最大展开相位差约100.53度，最大模pi差约84.70度；末节点展开差降至7.63度，因此“中间差异更大、两端较接近”的定性分离趋势存在，但强度明显大于原图。S0与S2最大三点范围约71.09度与25.25度，均在高端；mid的高能S0 eta低至.20695，是实际允许的明显非弹性现象，不是已违反盘约束。

## 同规则原图数值比较

[比较manifest](E3_paper_comparison_manifest.json)、[比较JSON](E3_paper_comparison/comparison.json)、[summary CSV](E3_paper_comparison/summary.csv)及[差值图](E3_paper_comparison/E_paper_differences.pdf)使用现有统一compare入口。九组均识别为原M50节点前缀；下表按相同无量纲s节点比较，以免把原图约0.3%的显示能标差当作相位差。没有调整生产mπ或振幅。

| 角色 | S0 RMS/度 | S2 RMS/度 | P1 RMS/度 | P1 RMS模180/度 |
|---|---:|---:|---:|---:|
| tip | 22.755 | 1.299 | 7.119 | 7.119 |
| mid | 11.523 | 2.098 | 25.093 | 20.099 |
| ref | 9.328 | 1.521 | 24.931 | 19.989 |

原图输入明确为 [Fig.9 P1](../../../references/figure9_p1_phases.csv)、[Fig.10 S0](../../../references/figure10_s0_phases.csv)、[Fig.10 S2](../../../references/figure10_s2_phases.csv) 的本地矢量marker读数及同名metadata；三点rho接近、S波高端与P波中段分离的原文在 [2309v3 TeX](../../../references/2309.12402v3-source/prd_submission_2.tex)1169–1172行，选点说明在1145–1149行。

这些是原图marker的描述性差值，不是统计chi-squared或误差条；原图角色与本方的对应是预声明的几何角色，不意味着完整作者系数身份相同。P1中段、tip高能S0仍是主要差异，不能由显示能标差消除。

绘图中灰实验/现象学曲线来自既有原图参考CSV，棕色细线来自原v3 Fig.9–10三套marker。旧separate-l2的C控制没有混入当前combined模型图；本轮新的C对照由主线单独交付。自动输出的FF引导lift/Watson诊断只是同一完整解的辅助观察，未用于改变主图或宣布连续极点。

本任务到此完成，无仍运行的生产进程。后续F必须基于明确的同一最终代表解生成流程；本报告不提前承诺分辨率变化会消除上述物理差异。

## 已保存中心路径的物理精度诊断

无需新求解，直接检查各运行最后两份已收敛且原式通过的完整中心。P1最大模π相位变化为tip .22895°、ref .18150°、mid .09404°；S0相应为.10418°/.39511°/.89488°。这些变化远小于当前三点P1和部分S0的主要曲线差异，不能继续把全部差异笼统归因于未完成中心化。数据见[E_CENTER_STABILITY](E_CENTER_STABILITY.json)。这是有限路径稳定性诊断，不是到精确极值的统一相位误差界；最终预选点可能是中心之后已验回的Newton步，原三点不替换。

另见[新tip全部1500原生幺正盘](E_TIP_FULL_NATIVE_ZH.md)：768-bit重新生成原生分波后全部通过，不仅复查主三波或旧float64 H。
