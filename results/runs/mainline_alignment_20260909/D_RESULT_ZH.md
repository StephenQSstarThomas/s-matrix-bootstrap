# 当前D：combined-l2 + 原文统一权重

来源见[原Fig5残差指纹](CHIRAL_SOURCE_RESIDUALS_ZH.md)、[截止决策](CUTOFF_DECISION_ZH.md)及[归一化转移](NORMALIZATION_TRANSFER_ZH.md)。当前合同为[PAPER_MAINLINE](PAPER_MAINLINE.json)：一个合并8维L2手征球，hard-midpoint统一权重，其他物理输入和actual-ρ上限不变。

[D3_paper_hull](D3_paper_hull/report.json)已给出完整4076变量联合见证，投影约(.07143138818,−.00455776509)，原1500盘、合并χ球、L4、100 Gram、四FESR及14FF全部通过。它使用新球明确重加约束，并没有继承旧separate的可行性；原B完整幅度仅作为初始化列。

D1_hard_M50的电流算子独立于χ范数，因此保持原文件；prepare阶段并不施加χ球。E仍需完成全空间支持和相位比较。以下保留此前只改截止规则的诊断。

# 此前hard-separate诊断

D1/D2已重建：M50同一散射H，100个原电流Gram、四个raw FESR与14个FF上界；截止采用Eq3.72统一π/M权重配合s_i≤s0。其余物理输入、χ及双密度界不变。[来源和包含关系](CUTOFF_DECISION_ZH.md)明确该选择与旧clipped变体的区别。

D3得到完整4076变量联合见证，坐标约(0.06068381549, -0.00380921130)，原1500盘、两个χ球、actual-ρ L4、100 Gram、四FESR及14FF全部验回。生成过程仍复用同设置B的完整幅度凸包及新电流SDP；这是存在性见证，不是E支持最优点。求解器标签为AlmostSolved，结论依据原式审计，而非该标签。

[准备](D1_hard_M50/report.json)、[完整见证与原式结果](D3_hard_hull/report.json)、[联合向量](D3_hard_hull/joint.npz)。原有clipped证书未转移给新模型。

E的两near-black位置在新相位求值前已经固定：[HARD_MAINLINE](HARD_MAINLINE.json)。使用原Fig8横坐标相对黑点的比例及正文140/92参考，得到ref=.07298305468121799、mid=.07745128711688606；tip仍是完整+x支持。此规则匹配几何角色，不假定恢复了作者完整向量。
