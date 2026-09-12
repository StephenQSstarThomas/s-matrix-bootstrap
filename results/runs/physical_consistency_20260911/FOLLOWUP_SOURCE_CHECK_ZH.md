# 原论文与后续公开实现的来源核查

核查日期：2026-09-11。对应 A1/A3（参数化与求值）、D（联合约束）、F2/F3（收敛与来源）。本次仅查原始论文和作者公开代码，未执行作者程序、未优化、未更改生产源码或物理输入。

**结论：原 2309 与后续资料都明确使用有限节点约束。后续资料补充了离散核、正则化、谱参数化和 Watsonian 选点，但没有提供把当前振幅的节点可行性提升为离节点／无穷远可行性的依据。** 因此后续实现可帮助解释方法，不能追认为原 Fig.3–11 的完整生产代码，也不证明当前失败的振幅不存在其他解析完成。

## 1. 固定版本，避免混用

| 一手来源 | 本次核查版本／日期 | 定位 |
|---|---|---|
| *Bootstrapping gauge theories* | arXiv:2309.12402v3，2024-12-04；本仓库目标 | [版本记录](https://arxiv.org/abs/2309.12402v3)，[PDF](https://arxiv.org/pdf/2309.12402v3)，本地原 TeX `references/2309.12402v3-source/prd_submission_2.tex` |
| *Gauge Theory Bootstrap: Pion amplitudes and low energy parameters* | arXiv:2403.10772v4，2025-08-08；首次提交 2024-03-16 | [版本及附件入口](https://arxiv.org/abs/2403.10772v4)，[全文](https://arxiv.org/html/2403.10772v4) |
| *The Gauge Theory Bootstrap: Predicting pion dynamics from QCD* | arXiv:2505.19332v3，2026-02-04 | [版本记录](https://arxiv.org/abs/2505.19332v3)，[全文](https://arxiv.org/html/2505.19332v3) |
| 同题正式发表版 | PTEP **2026 (6), 063B05**，2026-05-18 在线发表，2026-06-30 corrected/typeset；DOI **10.1093/ptep/ptag087** | [期刊正文及日期](https://academic.oup.com/ptep/article/2026/6/063B05/8686484) |

日期取版本履历；HTML 内自动排版日期不当成 arXiv 修订日期。正式发表版新增 Appendix D，不能拿 arXiv v3 的方程号直接引用正式版。

## 2. 2309 已明示与未明示的内容

- **节点和约束已明示。** pp.21–25，Eqs.(3.58–62) 定义圆盘映射、`φ_i=π(i−1/2)/M` 及离散谱变量；Eqs.(3.65–67) 明确 FF 节点实部的离散 Hilbert 核；Eqs.(3.70–71) 为带节点指标的 Gram／幺正约束。不能说原文完全没有给离散核。
- **有限高能条件已明示。** p.25，Eq.(3.75) 对有限集合 `s_i>s_0` 约束重标度 FF；它没有单独列 `F(∞)=0` 等端点等式。理论上的高能衰减与有限节点实现须分开。p.24，Eqs.(3.72–74) 给均匀 `π/M` 权重、Jacobian 和绝对残差形式；但打印求和没有完整说明截止穿过网格单元时的处理。不能因此把 raw 绝对误差说成毫无原文依据。
- **仍未恢复的实现选择。** §3 未明确给出所有单／双谱密度的连续插值族、其端点／尾部条件、散射在切上与离切线求值的统一有限函数构造，以及实际正则化范数的打包与缩放。Eq.(3.64) 仅称使用某种范数；本次不根据原图反推它。
- **收敛检查的范围。** pp.36–37、Appendix A/Fig.11 比较 `M=45,50,60` 与每同位旋 `L=8,10,12`；作者明确指出改变 M 同时改变变量、幺正条件和 SVZ 条件。这是有限分辨率敏感性检查，未报告固定振幅的连续幺正误差界。原文也没有 Watsonian 迭代作为 Fig.3–11 的步骤。

上述定位已逐项与本地 TeX 和 [2309v3 PDF](https://arxiv.org/pdf/2309.12402v3) 对读。作者仓库说明原 2309 当时未公开代码。[作者 README](https://github.com/hyfysics/gauge-theory-bootstrap#paper-snapshots)

## 3. 后续论文真正增加了什么

**2403v4** 的 §8.3–8.4（pp.57–68）明确采用谱变量配点／电流谱 B-spline 混合方案；切外积分求和与切上离散核分别见 Eqs.(8.31–35)，幺正在同一 `s_j` 集合上施加。映射中心改为 `ν_0=−20`；Eq.(8.56) 给双谱 L4 正则化。Eq.(8.72) 明确把 FESR 求和与归一化截止取到最后一个低于匹配能标的节点；Eq.(8.74) 仍是高能节点 FF 界。Appendix B.1/Fig.17 单项测试 `M:50→60`、最大自旋 `19→23`、`N_B:40→46`、`M_reg:100→500`。这没有提供独立加密幺正网格。[§8.3–8.4](https://arxiv.org/html/2403.10772v4#S8.SS3)，[Appendix B.1](https://arxiv.org/html/2403.10772v4#A2.SS1)

该文 §2.5（pp.19–20，Eqs.(2.27–29)）新增普通弹性梯度与 FF 相位驱动的 Watsonian 目标，可重复至收敛。它在原不等式内重新优化，并不把基本幺正性放宽；也不是对一个冻结振幅的画线修正。[§2.5](https://arxiv.org/html/2403.10772v4#S2.SS5)

**2505v3** 的 §2.2、§3（pp.7–13）进一步报告从不同极值／可行起点迭代得到相近结果。其 Eq.(2.8) 用旧 `|h| F/|F|` 构造目标；Appendix B（pp.31–32，Eqs.(B.20–22)）改用满足阈值幂且在无穷远有限的电流谱函数，并限制系数。这是额外电流谱假设和选点程序；它不等于散射或 FF 已满足无穷远端点条件。[§2.2–3](https://arxiv.org/html/2505.19332v3#S2.SS2)，[Appendix B](https://arxiv.org/html/2505.19332v3#A2)

**正式 PTEP 版** Appendix D、Eqs.(D.37–38) 更直接写明：幺正与 S0/P1/D0 正半定约束施加于 `j=1,…,M` 的离散点，采用 `M=50, l_max=21, n_max=50`。它仍未给固定变量数下的独立加密网格或无穷远认证。旧振幅相位替换现为 Eq.(2.10)，零模目标为 Eq.(2.18)；§3 的收敛是作者数值结果，不能从单步凸性推出全局唯一固定点。[正式版 §2.2、§3、Appendix D](https://academic.oup.com/ptep/article/2026/6/063B05/8686484)

这里关于“尚不能推出”的判断是本次逻辑范围审查，不是作者给出的反证。2505 两种目标一般不等价、Gram 一个零模不等于秩一的问题已有具名[历史代数记录](../followup_research_20260910/HELD_ISSUES_ZH.md)；正式版仍保留相关等同表述，不能当作已修复。

## 4. 作者代码的只读核对

公开来源确实存在：2403 arXiv 附件 `GTB_numerics.m/.nb`；作者 GitHub 提供[2403 冻结 release](https://github.com/hyfysics/gauge-theory-bootstrap/releases/tag/arxiv-2403.10772)（2025-09-30，commit `c8409331d5df911f592e53e5a43628418a39b213`）和[2505 release](https://github.com/hyfysics/gauge-theory-bootstrap/releases/tag/arxiv-2505.19332)（2025-10-08，commit `801684d9ece3de2918a20145178f569459081098`）。本次查询时 main 也指向后一 commit；这是读取到的身份，不推测尚未公开的生产设置。

| 代码定位 | 直接观察 | 不可转移的结论 |
|---|---|---|
| [2403 脚本 L7–15](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L7) | `M=50,l=10`；同一 midpoint 映射；`n0=sum(vnu<=s0)`；高能 FF 从 `n0+1` 开始 | “all energies”注释对应有限矩阵行，不是无穷远检查 |
| [2403 L59–65、L119–120](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L59) | `norm(rho/Mrho,4)<=100`；固定行幺正；初始优化后一次 Watsonian 重优化 | 不能只抄论文中“100”而忽略变量数及缩放，也不是原 2309 的代码身份 |
| [2505 驱动 L26–48、L91–101](https://github.com/hyfysics/gauge-theory-bootstrap/blob/801684d9ece3de2918a20145178f569459081098/theories/qcd/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m#L26) | 默认 `M=50,l=10,Na=49,nWatson=5`；目标使用旧 `abs(ht).*FF./abs(FF)` | `l=10` 每同位旋至最大自旋19，不能直接等同正式版 D.38 的21；谱系数计数也须另映射 |
| [2505 optimize_core L37–98](https://github.com/hyfysics/gauge-theory-bootstrap/blob/801684d9ece3de2918a20145178f569459081098/theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m#L37) | 幺正矩阵固定；电流 Gram 到 `n0`；FF 界从 `n0` 开始，包含最后一个低于匹配点的节点 | 与2403的 `n0+1` 边界选择有实际差别；无显式 `F(∞)=0` 或独立离节点约束 |

只读检查没有找到自动增加约束网格的循环；本次未查看输出曲线选择范数或参数。后续 chiral、FF、FESR 和谱定义本就有变化，不能把整套后续默认值搬回2309以获得相似曲线。

## 5. 对当前主线的用途

最小、已获来源支持的决定是**继续分别记录有限节点复现与连续物理诊断**。原 Fig.11 要求的 M/L 比较仍有必要；它既不会自动补救冻结振幅的离节点违例，也不能被另造的连续认证前置条件取代。若以后要加入新网格或端点等式，应先固定同一解析族及约束定义、报告它是新有限模型，再重新求完整联合见证；本次没有实施这种改动。

下载文件按 `.txt` 保存代码供阅读；未执行任何作者代码。来源 URL、字节数、SHA256 在 [downloads.json](source_check/downloads.json)，release/ref/tree 元数据及期刊 HTML 在 `source_check/`。2403 两文件的 Git blob SHA1 与其冻结 tag 相同；2505 脚本及核心与该 release 的冻结目录相同。原论文引用均对应仓库提供的 PDF／TeX；未恢复旧执行树。
