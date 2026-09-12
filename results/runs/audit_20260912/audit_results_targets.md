# 复现验收目标表与逐图判定 —— He–Kruczenski arXiv:2309.12402v3 §4 / Fig.3–11 / Appendix A

审计日期 2026-09-12。只读审计；论文页码按 PDF 页码（§4 起 p.25，Fig.3 p.27 … Fig.11 p.37）。仓库数值全部取自任务指定入口文件；"本审计计算"指我用 python 对数字化 CSV 与仓库 JSON 做的节点配对差值。

## 0. 核心缺陷判断摘要（≤300字）

论文核心主张（摘要 p.1、§5 p.35）：仅用 Nc、Nf、mq、ΛQCD、mπ、fπ 与凝聚，bootstrap 得到的 S0/P1/S2 相移"in good agreement with experimental results"；ρ 由 UV 信息决定且"resonance shifted by roughly 6%"、"results seem quite robust"；S0/S2"largely independent of the high energy data"。仓库结果否定或未支持其中三条定量主张：(1) 三代表 ρ 读数 803/701/696 MeV，跨度 107 MeV（论文 812.5–826.6 MeV，跨度 14 MeV），且 mid/ref 在共振处 η 低至 0.43/0.24，非弹性"共振"；(2) tip/mid 的 S0 在 1.2 GeV 达 178°/190°（论文 87–106°），即加入 UV 后 S0 剧变，直接与"S0 与高能无关"矛盾；(3) Fig.8 强非对称与 Fig.11 固定 M 的 L 稳定性（46 vs 12.5 MeV）未复现，且仓库可行区域比论文宽约 2 倍、+x 端远 28%。这些属于物理模型层面差异（约束集合/约定与作者不同），是核心缺陷：审稿人会认定论文主结果尚未被复现，但也不能据此判定论文错误。

## 1. 论文验收目标表

数字化精度（各 `.metadata.json`）：PyMuPDF 矢量路径中心提取，工作分辨率 0.01 pt（Fig.3 折合 f00 6e-5、f11 1.6e-5；Fig.11 刻度拟合残差 <0.002 pt）；能量节点与 M=50 中点公式吻合到 3e-6 GeV，说明作者作图用 mπ=139.57 MeV 而非正文 140（`figure9_p1_phases.metadata.json: independent_energy_node_check`）。metadata 明示"published graphical marker coordinates, not original optimization output"。

| 图 | 坐标/物理量 | 论文定量陈述（原文） | 数字化 CSV 关键数值（本审计计算） |
|---|---|---|---|
| Fig.3 (p.27) | 横 f00(s=3)、纵 f11(s=3)，mπ=1 单位；纯 S-matrix bootstrap，M=50、10 分波/同位旋 | "Inside the shape are all the possible values … under analyticity, crossing and unitarity"; "the shape in fig. 3 will play no further role" (p.26) | 498 点；f00∈[−2.902, 2.233]，f11∈[−0.734, 0.079]（`figure3_boundary.csv`） |
| Fig.4 (p.28) | 同平面，仅 f00>0；六个 εχ=6e-3…2e-4 | "shrinks to a thin region around … f11(3)=−f00(3)/15 (4.79)"; "εχ=0.002 (green points) … allows the physical value of fπ"; "the shape terminates … lower bound on fπ" | 无独立 CSV；ε=0.002 绿形与 Fig.8 绿组相同：f00 max=0.0826，xref 处上/下边界 −0.004368/−0.005128，宽 0.00076（`figure8_boundary.csv` green） |
| Fig.5 (p.29) | 0<s<4 的 f00、f20、f11 实值分波，蓝/橙/绿=ε 6e-3/4e-3/2e-3 | "most notable deviation from linearity is in the S0 wave since the chiral zero … disappears for the blue points"; "green points match better the linear prediction" | S0 零点：绿 s=0.425，橙 s=0.305，蓝无（f>0 全程）；f00(3)=0.0725/0.0758/0.0797；f11(3)=−0.0043/−0.0037/−0.0032（`figure5_subthreshold.csv`） |
| Fig.6 (p.30) | ε=0.002 绿形 + 黑点(fπ=92) + 洋红选点 | "choose a point closest to the black dot"; "only points at the boundary (green points) have partial waves" | 黑点(0.07131, −0.004754)（与 mπ=139.57、fπ=93 一致，metadata `independent_chiral_black_check`）；仓库黑点用 140/92 → (0.07332, −0.004888) |
| Fig.7 (p.31) | δ00、δ20、δ11(°) vs E(GeV)，仅手征约束（洋红）vs 实验 [38,39] 与拟合 [54] | "remarkably compatible with experimental results for the S0 and S2 waves but not the P1"; "P1 … requires further input from QCD" | 无仓库级 CSV；图读：S0 在 0.9 GeV 约 100° 后回落至 ~75°；P1 在 1.2 GeV ≲20° |
| Fig.8 (p.32) | 同 Fig.4 平面：绿=仅手征，青=加 FESR(εSR=2e-3, εFF=6e-5)，红/粉/浅粉三选点 | "The upper boundary shrinks notable but the lower not so much"; "the upper curve is modified much more than the lower one" | xref=0.0733 处：上边界收缩 2.32e-4，下边界上移 3.7e-5，比 6.3（x=0.02–0.08 比 1.3→7.6）；青 f00 max 0.0811（绿 0.0826）；选点 red(0.0811,−0.00525)、pink(0.0753,−0.00476)、light_pink(0.0710,−0.00442) |
| Fig.9 (p.33) | δ11 vs E，三选点 | "resonance … crosses π/2 is slightly shifted from … 770 MeV by roughly 6%"; "shape of the resonance seems quite good"; "results seem quite robust"; "coincide both at low and high energy but they spread in the middle" | 90° 过点：light_pink 812.5、pink 824.3、red 826.6 MeV（+5.5/+7.1/+7.4%），三点跨度 14 MeV；0.7294 GeV 节点 δ11≈31–33°；1.2 GeV 末值 169–173° |
| Fig.10 (p.34) | δ00、δ20 vs E，三选点 | "three bootstrap curves … coincide at low energy but they spread at high energy" (S0); "similar but not equal phase shifts"; "agreement with experiment is remarkable" | S0：0.51 GeV 43–51°，0.79 GeV 76–99°，1.06 GeV 96–109°，1.196 GeV 99.6/87.3/106.4°（红/粉/浅粉），末端跨度 19°；S2：1.196 GeV −36.4/−26.4/−18.3° |
| Fig.11 (p.37, App.A p.36) | 五组 (M,L) 的三波相移 | "For M=50 … L=8,10,12 … reasonably good convergence"; "S0,S2 … very good convergence, while the P1 … peak slightly shifted depending on M" | P1 90°：(50,8) 820.1、(50,10) 832.7、(50,12) 831.6、(45,10) 786.6、(60,10) 778.3 MeV；固定 M50 L 跨度 12.5 MeV，固定 L10 M 跨度 54.3 MeV；S0 五组在 1.0 GeV 87–103°（跨 16°），1.19 GeV 89–106° |

与实验对比的陈述：§4.2 "S0 and S2 waves are correct but the P1 evidently not"（Fig.7 图注）；§4.3 "shifted … 770 MeV by roughly 6%. This might be improved by choosing other value of s0"；"we consider the results to be compatible with the real world QCD"（p.33）；§4.3 末 "the shape of the S0 wave is a consequence of pion low energy dynamics whereas the ρ is a consequence of QCD … agrees with our findings"（p.35）。

## 2. 逐图复现判定

仓库最终模型（`refined_phases/phases.json: model_signature`）：M50/L10、1500 盘、pv-midpoint、separate-l2 χ=0.002、B=377500、printed 四矩 raw 绝对 0.002、hard-midpoint、εFF=6e-5、mπ=0.14、s0=1.2 GeV、αs=0.4。

| 图 | 仓库最终数值（文件+字段） | 与论文/数字化差距 | 仓库判定 | 独立判定 / 差距层面 |
|---|---|---|---|---|
| Fig.3 | STATUS "B2 沿用 38 份合法支持；整体几何误差 .00967577"；+x 见证 x≈1.6066（STATUS B2/B3 段）；`CLAIM_LEDGER.steps[B2]`="历史有限几何已有通过记录" | 入口文件无与 `figure3_boundary.csv` 的直接比较；论文 +x 极值 2.233，仓库仅给可行见证 1.607（下界，差 −28%）；"几何误差 .0097"是仓库内外壳自洽度，不是对论文的偏差 | "通过记录"（历史） | (b) 部分复现/未验证。层面：无法判定；至少 +x 支持未证明达到论文值 |
| Fig.4 | `refined_regions/regions.json: IR_outer_vertices` x max 0.10594，`reference_section_certificate.IR_upper/IR_lower`=−0.004117/−0.005524（xref 处宽 0.00141） | 论文 ε=0.002 绿形 f00 max 0.0826（仓库 +28%）；xref 处宽 0.00076（仓库 1.85 倍）；六 ε 区域 STATUS 称".01 标准尚未全部达到" | B3 "全区域尚未达标"，另称"印刷图不是数学排除证书" | (c) 未复现（定量）。层面：物理模型（χ 范数/打包、B 未恢复导致可行集为论文的超集） |
| Fig.5 | `C_comparison/comparison.json: C1_0.00x_{S0,S2,P1}.direct_coordinate.rms/max_abs`：S0 rms 0.0040/0.0016/0.0015，max 0.0062/0.0025/0.0049；S2 rms ≤0.0013；P1 rms ≤0.001；零点：橙 [0.3,0.4]、绿 [0.4,0.5] 已证，蓝采样未见（STATUS C1） | 幅值差 ≤ 8% of f00(3)；零点位置与数字化（0.305/0.425/无）一致 | C1 "采样线性改善；零点已证" | (a) 已复现（半定量）。注意：该链属 20260909 combined 范数支线，且范数识别参考过 Fig.5 输出（STATUS "输出辅助"），非当前 separate 主线独立结果 → 选点/约定层面保留 |
| Fig.6 | C2 "原样锁定 C1 绿色 3876 系数，选点先于相位"；黑点 (0.07332,−0.004888) | 黑点约定差 (−0.0020, +0.00013)（139.57/93 vs 140/92） | 已交 | (a) 已复现（选点/约定层面差异，可接受） |
| Fig.7 | `C_comparison: C3_paper_bootstrap_S0/S2/P1`：rms 5.5°/0.55°/2.8°，max 12.2°/2.5°/12.0°，bias +3.8°/+0.2°/−1.2°；当前主线 IR 基线 `phases.json: chiral_only_baseline`：S0 max 91.7°@0.864，末 85.8°；P1 末 8.56° | IR-only 下 S0/S2 接近、P1 无 ρ，定性与论文一致；P1 在 E>0.75 偏低平均 7°，max 12° | C3 "新完整 IR 对照通过" | (b) 部分复现（数值精度层面，≤12°）。同样注意比较数据来自 combined 支线，当前 PV 基线无对论文的数值比较 |
| Fig.8 | `refined_regions/regions.json: reference_section_certificate`：upper_shrink [1.959e-4,1.991e-4]，lower_rise [1.820e-4,1.949e-4]，ratio [1.0049,1.0938]；UV_upper/UV_lower −0.004314/−0.005335（宽 0.00102）；UV +x tip 0.09910 | 论文上收缩 2.32e-4（仓库上收缩量吻合，差 15%），但论文下上移 3.7e-5，仓库 1.9e-4（5 倍）；论文比 6.3 vs 仓库 ≤1.09；仓库 UV 宽度为论文青形 2.1 倍；UV tip x 0.0991 vs 0.0811（+22%） | E1 "只复现轻微不对称"，"明显非对称未复现" | (c) 未复现（非对称是该图唯一定量陈述）。层面：物理模型（下边界几乎不动是论文 UV 约束的特征，仓库不具备） |
| Fig.9 | `phases.json: profiles[tip/mid/ref].P1_first_upward_90.linear_in_energy_gev`=0.80339/0.70103/0.69649；`min η(P1)`=0.777/0.435/0.242（`REFINEMENT_SUMMARY.json`）；节点配对最大差 21°/110°/115°（本审计）；`three_curve_separation.P1` 中段跨度 >100° | 论文 812.5–826.6 MeV、跨度 14 MeV、全部 +6%；仓库 +4.3/−9.0/−9.5%，跨度 107 MeV；仓库 mid/ref 共振宽仅一个节点（0.68→0.73 GeV 相移跳 84°）且强非弹性 | E3 "ρ稳健性未通过，major" | (c) 未复现。层面：物理模型（可行集不同 → 极值点不同）+ 选点（仓库 tip/mid/ref x=0.0991/0.0862/0.0733 均比论文 0.0811/0.0753/0.0710 偏 +x） |
| Fig.10 S0 | `phases.json: profiles[].phase_degrees[:,0]`：1.196 GeV tip 178.2°、mid 189.6°、ref 112.4°；`three_curve_separation.S0.range_at_last_node`=77.1°；`E3_paper_comparison: E_paper_tip_S0.rms/max_abs`=23.1°/78.7°（20260909 点） | 论文末端 99.6/87.3/106.4°，跨度 19°；仓库 tip 高 79°、mid 高 102°、ref 高 6°；E≥0.9 平均偏高 76°/53°/6°（本审计） | E3 "S0 系统偏高，major" | (c) 未复现（tip/mid）、(b) 部分（ref）。层面：物理模型；且与论文"S0 与 UV 无关"主张直接冲突（IR 基线 S0 末 86° → 加 UV 后 178°） |
| Fig.10 S2 | 同上 [:,1]：1.196 GeV −39.1/−24.3/−12.8°；节点配对 rms 1.5°/0.6°/1.5°（tip↔red、mid↔pink、ref↔light_pink，本审计） | 论文 −36.4/−26.4/−18.3°；最大差 ≤5.5° | "S2 仅有比较差异" | (a) 已复现（数值精度层面） |
| Fig.11 | `FINAL_SCIENCE_SUMMARY.json: F_results[].native_90_mev`=835.5/803.4/789.3/749.3/773.5 MeV；`fixed_M_L_spread_mev`=46.21，`fixed_L_M_spread_mev`=54.06；`PV_F_five_complete/resolution.json: profiles[].phase_degrees` S0(1.0 GeV)=162–169° 五组；`M45_M60_reused_not_recomputed: true` | P1：逐组差 +15/−29/−42/−37/−5 MeV；L 跨度 3.7 倍；L 序反转（论文 L8 最低，仓库 L8 最高）；M 跨度同量级但 M45/M60 序反转；S0 五组均比论文高 ~75°（论文 87–103°）；S2 差 ≤7° | F2 "L 未对齐；M 不能笼统判失败"；"中高能 S0 差异仍在" | (c) 未复现（S0 与 L 稳定性）、(b) 部分（P1 位置在 0.6–5% 内、S2）。层面：物理模型；"M 跨度同量级"是巧合级指标 |


### 2.1 逐图差距解读（层面归因）

- **Fig.3**：仓库把纯 bootstrap 区域视为"已有通过记录"，但入口文件中的证据只有内外壳自洽误差与一个 +x 可行见证（1.607）。论文 Fig.3 的边界是线性泛函极大值的轨迹，其 +x 极值 2.233 是可直接验收的数；仓库若已有 38 份支撑，应能给出各方向支撑值与数字化极值的逐向比较，目前没有。`results/figures/reproduction_progress.png` 中的对照点是早期低阶（N=4/6/8）候选，全部位于论文边界内侧很远，不能当最终结果。故本图只能定为"未验证"。
- **Fig.4**：论文明确写出 ε 序列与"shape terminates"，即区域在 +x 方向有终点。仓库 ε=0.002 区域终点 0.1059 比论文 0.0826 远 28%，且在物理 fπ 截面上宽 1.85 倍。若 χ 范数为 4 维 L2 且容差为 0.002，区域尺度对范数归一化敏感（仓库自己在 CORE_CLAIM_LEDGER §6 指出 αN(v)≤ε 等价于 N(v)≤ε/α），所以这是约定层面进入物理模型层面的典型差距：约定未恢复导致约束集合不同。
- **Fig.5/6/7**：这三张图是 IR 链，仓库的定性结论（S0 手征零点随 ε 变小而出现、绿色最接近线性、IR-only 无 ρ、S0/S2 接近实验）全部成立，数值差在幅值 ≤8%、相移 ≤12°。问题在于链条身份：这些比较来自 combined 范数支线，而最终 UV 链用 separate 范数；两条链的 IR 代表振幅不是同一份。论文的逻辑是"同一 ε=0.002 区域内选点 → Fig.7 → 加 UV → Fig.8–10"，仓库尚未在同一合同下把 Fig.5→7 与 Fig.8→10 接通。
- **Fig.8**：论文唯一定量陈述是上/下边界改变量的悬殊。数字化显示上收缩 2.3e-4、下上移 3.7e-5，比值随 x 从 1.3 升到 7.6，即越靠近 tip 非对称越强。仓库的上收缩量 1.96e-4 与论文相近，说明 FESR 对上边界的约束作用被复现；差距全部来自下边界：论文几乎不动，仓库上移 1.9e-4。下边界对应 f11(3) 更负（fπ 更小、耦合更强）的振幅，论文的 UV 约束对它们几乎无效，仓库的却显著有效——这是约束集合不同的直接证据，属于物理模型层面。
- **Fig.9**：论文三点 ρ 都在 812–827 MeV，相互差 14 MeV，且共振区跨越 0.79→0.86 两节点相移升 45–49°，形状平滑。仓库 tip 的形状与此相似（803 MeV，跨两节点升 51°，min η 0.78），单看 tip 可称"部分复现"；但 mid/ref 在 0.68→0.73 GeV 一个节点内跳 84°、η 降到 0.24–0.43，这不是论文式的 ρ。仓库三点坐标（0.0991/0.0862/0.0733）与论文三点（0.0811/0.0753/0.0710）不重合——因为仓库区域更大，同样的"tip/近黑点"规则落到了不同位置。差距是模型层面（区域）叠加选点层面（规则相同但区域不同）。
- **Fig.10**：S2 三点全部在 5.5° 内，且颜色对应（tip↔red 最负、ref↔light_pink 最浅）与论文一致，是本次复现最干净的成功项。S0 则分裂：ref 与论文 light_pink 的 rms 9.6°，可接受；tip/mid 在 0.9 GeV 以上急升到 160–190°，论文三条都停在 87–109°。论文正文对 S0 的定量描述是"coincide at low energy but spread at high energy"，仓库低能确实重合（<0.6 GeV 平均差 3–6°），高能跨度 77° 而论文 19°。
- **Fig.11**：论文用此图证明数值参数不敏感。仓库五组之间 S0 相互一致（1.0 GeV 跨 7°）但整体偏高 75°，说明"偏高"不是离散化噪声而是模型系统偏差；P1 的 M 依赖同量级但 L 依赖大 3.7 倍且方向相反（论文 L8 最低、仓库 L8 最高）。仓库对 L12 做了自协调下降修补，修补后读数 789 MeV 与修补前 789.08 几乎不变，说明 L 跨度确为模型结果而非收敛问题（仓库亦如此表述）。

## 3. 核心缺陷判断

论文核心主张（摘要 p.1；§5 p.35）："a numerical implementation of the method gives the phase shifts of the S0, P1 and S2 waves in good agreement with experimental results"；"the S0 and S2 waves are largely independent of the high energy data, whereas the P1 is determined by the UV information"；"This information seems enough to determine the P1 wave"；§4.3 "the results seem quite robust"。

仓库已支持的部分：(i) 定性机制"UV → ρ"成立：IR 基线 P1 末值 8.6°，加 FESR/FF 后三点均上穿 90°，且冻结 IR 点被 UV 支撑平面严格排除（`phases.json: chiral_only_baseline.UV_projection_exclusion.excluded=true`）；(ii) S2 三点与论文差 ≤5.5°；(iii) Fig.5–7 的 IR 定性图像与论文一致。

直接否定/未支持核心主张的差距：
1. ρ 不稳健：803/701/696 MeV，跨度 107 MeV（论文 14 MeV）；两点低于 770 MeV，与"shifted by roughly 6%"（向上）方向相反；mid/ref 的 min η(P1)=0.43/0.24 表明"共振"由强非弹性驱动，与实验 ρ（弹性）不符。论文图无 η，但"agreement with experiment"隐含近弹性。
2. S0 对 UV 敏感：tip/mid 的 S0 在 1.2 GeV 达 178°/190°，比 IR 基线（86°）高 90°+，直接违反"S0 … largely independent of the high energy data"。仓库自证 Re S0(1.199)>0.881 在近 tip 区域全域成立——即在其模型内这一冲突不可靠精度修复。
3. 区域几何为论文的超集：ε=0.002 区域 +x 端 0.1059 vs 0.0826，xref 处宽 1.85 倍；UV 区域宽 2.1 倍、下边界上移 5 倍于论文。这说明约束集合（χ 范数/打包、B、FESR 误差口径）与作者不同，随之极值代表点全部偏 +x，是 1、2 两项的上游原因。
4. Fig.11：L 稳定性 46 vs 12.5 MeV；S0 五组一致偏高 75°。

判定：以上属核心缺陷。审稿人拿到这份复现会认为论文主结果（Fig.9–10 的三点稳健相移、6% ρ 偏移、S0 由 IR 决定）尚未被复现；但因作者数值合同未唯一恢复（仓库 `CLAIM_LEDGER.json: full_original_parameter_impossibility_proved=false`），也不能反推论文有误。区分：Fig.10 S2、Fig.5/6/7 属数值精度或约定层面；Fig.4/8/9/10-S0/11 属物理模型层面。


**审稿人视角。** 若把这份复现作为独立验证提交，审稿人会问三个问题：(1) 是否得到与论文相同的可行区域？——否，区域是论文的超集，宽约 2 倍；(2) 是否在区域上得到相同的三条相移？——S2 是，S0 仅一条接近，P1 仅一条形状相似且 ρ 位置分散 107 MeV；(3) 是否验证了"S0 由 IR 决定、ρ 由 UV 决定"的物理结论？——后半句定性成立，前半句在仓库模型内被自己的证书否定。因此结论只能是"论文主结果未被复现"，而非"论文主结果被证伪"。仓库对此的自我判定（`CLAIM_LEDGER.json: quantitative_paper_success=false`；STATUS "论文稳健定量核心仍未复现成功"）与本审计一致。

**复现成功所需的最小条件。** 按论文文本可验收的量：Fig.8 上/下改变量比 ≥3（论文 6.3）；三点 ρ 90° 读数在 770 MeV 的 +3%…+9% 内且互差 ≤30 MeV；三点 S0 末端 ≤120° 且互差 ≤30°；S2 差 ≤10°；固定 M50 的 L 跨度 ≤25 MeV。这些阈值是本审计据数字化 CSV 与论文措辞（"roughly 6%"、"not so much"、"reasonably good convergence"）提出的建议，论文本身未给统一阈值。当前仓库满足其中 S2 一项。

## 4. 仓库报告的诚实度与自洽性

**"严格排除"的适用范围。** `tip_S0_certificate_final/probe.json: scope`="All complete feasible amplitudes in this declared finite near-support set; no new physical constraints, pole or continuum claim"；`ref_P1_certificate_final/probe.json` 同。MAJOR_CLAIMS §6 明写"不是这样的无条件结论……不能说这些原物理输入在整个可行域中分别禁止正确符号"。因此 Re S0(1.1994)>0.88103、Im S1(0.7317)<−0.5055 是仓库有限模型 D 在其自选近支撑层内的对偶证书，不是对论文的反证；仓库没有把它误写成对论文的反证。但 STATUS 顶部把"严格差异……排除原图相应相移符号"列为首条要点，措辞易被读成"论文曲线被排除"；且该证书的实质只是"我们的 D 与论文的可行集不同"的另一种表述（D 本身比论文宽 2 倍），其科学增量有限。

**"通过"的含义。** 18 步表与 `CLAIM_LEDGER.steps[].completed_scope` 中的"通过"多为内部自洽（primal/支撑/中心/107 项测试），STATUS 也声明"测试通过不等于物理复现成功"。但两处"通过"易被过度解读：B2 "Fig.3 历史有限几何已有通过记录"——入口文件中没有与数字化 Fig.3 的比较，"几何误差 .0097"是内外壳自洽；C3 "新完整 IR 对照通过"——新 PV 基线只有 132 项内部检查，对 Fig.7 的数值比较（rms 5.5°等）来自 20260909 combined 支线。仓库自己在 CORE_CLAIM_LEDGER §5 承认"不能追认证明本轮 separate 合同下的独立三色链"，与此一致，但顶层表格未标明这一差别。

**"做了但差很远"的定级。** Fig.9/10-S0 定为 major、Fig.8 定为"只复现轻微不对称"、Fig.11 定为"L 未对齐"——定级恰当。偏乐观之处：F2 反复强调"M 跨度 54.06 vs 54.35 同量级，不能笼统判失败"，但逐组差达 37 MeV、M45/M60 序反转、五组 S0 均偏高 75°，"同量级"是单一巧合指标；Fig.8 汇总只报比值区间，未指出仓库上收缩量其实与论文吻合（1.96e-4 vs 2.32e-4）而下边界上移是论文 5 倍、区域宽度 2 倍——这一分解更能定位模型差异，却未出现在顶层。

**自洽性核对。** MAJOR_CLAIMS §4/§5 引用的论文节点值（S0 末 106.41/87.26/99.57°；P1 0.7294 GeV 约 31°）与我从 CSV 独立读出一致；`FINAL_SCIENCE_SUMMARY.fixed_M_L_spread_mev`=46.2096 与 F_results 一致；`M45_M60_reused_not_recomputed=true`、`quantitative_paper_success=false` 均如实标注；"没做"未被包装成"做了"的情况未发现，反而普遍使用大量限定语。整体结论：仓库报告诚实、过度谨慎；主要风险是把内部证书的严密性（Arb 区间、对偶界）与对论文的复现程度并列陈述，使读者可能高估已完成的物理对齐。

**具体引文核对。** STATUS 顶部："在 x≥.09909604549 的近 tip 区域，Re S0(1.1994 GeV)>.88103……它们排除原图相应相移符号，不能靠同区域另一完整振幅或 π 分支修复；不是对所有原参数实现的无条件排除。"——前半句"排除原图相应相移符号"在语法上以论文相移为宾语，但证书对象是仓库的 D；建议改为"排除本模型近 tip 区域出现原图相移符号"。MAJOR_CLAIMS §13："本轮可修补的实际数值问题已完成修补；当前明确有限原型的两个近支撑相位差异已获严格证书。论文稳健定量核心仍未完整复现。"——三句层次清楚，属诚实表述。STATUS 2026-09-11 段："mixed-PV 是本仓库按原文约束建立的具名有限实现……不能称为作者原代码。"——明确承认模型身份差异。CORE_CLAIM_LEDGER §5："历史 separate/PV 区域中 pure 有 38 份支持……这是存档所声明有限绘图域/metric 下的完成状态"——如实限定了 Fig.3 "通过"的含义，但 STATUS 18 步表中的 B2 行没有带这一限定。

**未发现的问题。** 没有发现用相移反选输入、放宽容差、改质量或截止以贴图的迹象（`phases.json: selection_uses_phase_data=false, experimental_phase_fitting=false`；F 五组 `criterion="Same +f00(3) UV representative; no phase ranking or Watsonian step"`）；失败记录（L12 修补前、B/2、combined 重验）均保留。
