# He–Kruczenski 2309.12402v3 数值规定审计：论文到底规定了什么，仓库又自造了什么

审计范围：只读 `references/2309.12402v3.txt`（全文）、`references/2309.12402v3-source/prd_submission_2.tex`（下文 "TeX L"）、PDF 图 3–11、作者公开代码库 `references/upstream-gauge-theory-bootstrap`（README + `papers/arxiv-2403.10772/GTB_numerics.m`）、以及仓库 `SCIENCE.md/STATUS.md/README.md` 与 `results/evidence/PV_PRIMARY_SOURCE_AUDIT_20260911/AUDIT_ZH.md`。未运行任何生产计算，未改动仓库文件。

---

## 0. 核心裁决

**结论：论文把数值问题规定到了"标准 SDP/SOCP 可直接实现"的程度，但有两处物理上要紧的量没有唯一定义：FESR 容差 ε^SR 的归一化口径，以及 χ 约束 (3.64) 的范数打包。除此之外，仓库所列"未恢复"项绝大多数要么是论文已明确规定（M/L、节点、Hilbert 核、FF 界、xref），要么是任何实现者都会做的常规离散化（PV 核、hard 截止、径向扫描），要么是仓库自己换了表示/求解器后才产生的自造问题（analytic-cardinal、finite-sine、T0=0、infinity=zero、3123/3582 盘、Newton-barrier 支持求解器、Arb 高精度"证书"）。**

责任划分（置信度约 80%）：
- **论文的锅（约 30%）**：(a) TeX L1030–1033 的 ε^SR 只写 `||...−QCD value|| ≤ ε^SR`，而 (2.56) 的"QCD value"是按 1/s0^{n+2} 归一化的 ~10⁻⁷–10⁻⁵ 量级；ε^SR=2×10⁻³ 若按字面套在归一化矩上是 10³–10⁴ 倍的空约束，若套在 raw 矩上则 S0 n=0 盒宽是目标的 84%、P1 n=−1 是 47%（本审计独立算得），两者都不可能是作者实际用的口径。而 P1 两矩的比值直接固定谱密度的平均能量（≈822 MeV），即 ρ 质量。这是与"仓库 ρ 质量 803/701/696 MeV 对论文 ≈815 MeV"直接相关的**唯一**论文级缺口。(b) (3.64) "with some norm"（TeX L951）。(c) 作者 README 明言 2309 "No public code was provided at the time"，所以"原 producer 不可得"是事实，但这不等于"问题不可复现"。
- **仓库的锅（约 70%）**：论文要求的是一个约 4000 实变量、1500 个二阶锥 + 100 个 3×3 Hermitian PSD 块 + 4 条线性 FESR + 14 条 FF 上界 + 1–2 个范数球、线性目标的凸锥规划（作者后续 2403 代码正是 `cvx_begin sdp; cvx_solver mosek` 一段 60 行完成）。仓库把它换成自写的 Newton/barrier"支持"求解器、Arb 384–768 bit 认证、"analytic-cardinal/finite-sine" 离节点解析族、5 MeV 细网格附加幺正盘（3582 盘）、T0=0 / 高能渐近条件等一整套论文里没有的对象，然后再去"恢复"这些对象在原论文中的对应物——这些对应物本来就不存在。"B 密度正则化"是从作者 2403 代码借来的，2309 正文一字未提。

最关键 3 条证据：(1) TeX L926/L993/L999–1011：变量集 (3.62)/(3.69) 与约束 (3.70)/(3.71) 已完整、封闭；(2) TeX L1030–1033 + L876–880：ε^SR 的口径与 (2.56) 归一化互相矛盾；(3) 作者 2403 快照 `GTB_numerics.m` L51–111：同类问题在 CVX+MOSEK 下一段即解，且 χ 用 8 维合并 L2（L102）、FESR 用逐通道归一化 L2 球 ε=10⁻⁷/6×10⁻⁶（L39, L91）——与仓库当前"separate-L2 + raw 绝对 .002 盒"都不同。

---

## A. 论文规定的数值问题清单

### A1 变量与参数化
- **振幅**：Mandelstam 表示 (2.7)，TeX L554–558：`A = T0 + (1/π)∫σ1/(x−s) + (1/π)∫σ2[1/(x−t)+1/(x−u)] + (1/π²)∬ρ1/(x−s)[1/(y−t)+1/(y−u)] + (1/π²)∬ρ2/((x−t)(y−u))`，ρ2 对称。变量 (2.8)：`T0, σ_{1,2}(x), ρ_{1,2}(x,y)`（L561）。**T0 是自由变量**。
- **同位旋/分波**：(2.4)–(2.6) L538–550；分波 (2.9) `f = ¼∫P_ℓ T^I dμ`，t = −(s−4)(1−μ)/2 (L566–571)；S = 1 + iπ√((s−4)/s) f (2.11, L574)。
- **共形映射**：(3.58) `z = (2−√(4−ν))/(2+√(4−ν))`（L910），即 1708.06765 (rhoVar) 取 s0=0 的特例（1708 源 L815 明写同式）；(3.60) ν(φ)=8/(1+cos φ)。
- **网格**：(3.61) `φ_i = (π/M)(i−½), i=1..M`（L922）。这是**整个上半圆**的 M 个中点节点，s_i 从 4.001 到 1.6×10⁴（M=50 时 43 个节点 ≤ s0=73.47）。
- **离散变量** (3.62)：`{T0, σ_{α,i}, ρ_{α,ij}}`，ρ2 对称（L926）→ 1 + 2M + M² + M(M+1)/2 = 3876 实数（M=50）。
- **M、L 取值**：§4.1 L1065 "M=50 … unitarity for 10 partial waves per isospin"；Appendix A L1206：`M=50, L=10`；测试 (50,8),(50,10),(50,12),(45,10),(60,10)。"10 partial waves per isospin"按作者 2403 代码 L8/L68–69 的索引方式是每同位旋 10 条（I=0,2：ℓ=0,2,…,18；I=1：ℓ=1,…,19），故幺正盘数 = 3·L·M = 1500。
- **形状因子**：(2.37) 一次减除色散（L723）；离散 (3.65)–(3.67) L967–985：`ReF_i = 1 + K_ij ImF_j`，`K_ij = K̃_{i+j−2M−1} − K̃_{i−j}`，`K̃_m = (1−(−1)^m)/(2M)·cot(mπ/2M)`（引 [9] He:2018uxa）。谱密度 ρ_{ℓ,i} 为节点变量 (3.68)。总变量 (3.69) L993。

### A2 约束
- **解析性/交叉**：由 (2.7) 自动满足（L931 "The parameterization used already satisfies crossing and analyticity"）。
- **幺正性**：(2.12) |S|≤1，离散实现为 2×2 PSD (3.71) L1009——即每个 (I,ℓ,i) 一个圆盘/二阶锥；其它波只用 (3.71)，S0、P1 用 3×3 (3.70) L999–1007：`[[1,S,𝓕],[S*,1,𝓕*],[𝓕*,𝓕,ρ]] ⪰ 0`，𝓕 由 (2.33) 的运动学因子换算（L701–705）。论文没说 PSD 只施加到 s0 以下还是全部 M 点（2403 代码只到 n0）。
- **手征对称破缺**：§2.2 线性 Weinberg 振幅 (2.17)，分波 (2.18) L604；§3.2 用比值 R^χ_{01}=3(2s−1)/(s−4)，R^χ_{21}=3(2−s)/(s−4) (3.63) L939；约束 (3.64) L945–950：`||f00(s_j) − R01 f11(s_j)|| ≤ ε^χ`，`||f20(s_j) − R21 f11(s_j)|| ≤ ε^χ`，**原文 "with some norm and tolerance ε^χ"**（L951）；s_j = 1/2, 1, 3/2, 2（L951）。ε^χ 扫描值 6e−3,4e−3,2e−3,1e−3,6e−4,2e−4（Fig.4 caption L1083），选定 **ε^χ = 0.002**（L1079："The value ε^χ=0.002 (green points), that we now choose"）。注意 (3.64) 写成两行、各带 ‖·‖，字面上更像两个独立范数；但作者 2403 代码 L102 是把 8 个残差堆成一个向量取 L2。
- **形状因子 bootstrap**：(2.41) Gram 矩阵 L749–770；(2.42) |𝓕|² ≤ ρ；流：jS=m_q(ūu+d̄d)、jV（2.24）；只做 S0、P1（脚注 15，L962）。
- **SVZ/FESR**：(2.47) L803–807 含 αs 修正与两个凝聚项；FESR (2.48) L819；矩公式 (2.50) L827–835；**每波两条**：S0 取 n=0,1，P1 取 n=−1,0（L836, L1027）。离散 (3.72)–(3.73) L1021–1027：`∫₄^{s0} ρ x^n dx → (π/M)Σ_i (ds/dφ)_i s_i^n ρ_i`，ds/dφ = 8 sinφ/(1+cosφ)²——注意求和写的是 i=1..M，而积分上限是 s0，截止如何处理没说。容差 (3.74) L1030：`||(π/M)Σ … − QCD value|| ≤ ε^SR`，"The QCD values are those from (2.56)"（L1033），**ε^SR = 2×10⁻³**（L1141）。(2.56) L876–880 的数值是 `(1/s0^{n+2})∫ρ x^n`：S0 = 3.09e−8{27.38/(n+2)+0.61δn}，P1 = −4.34e−6{−13.26/(n+2)+0.41δn}。本审计独立复算：这些系数与 N_f m_q² → m_u²+m_d² 一致（8.40e−7 vs 印刷 8.46e−7），用算术平均 m_q 会差 9%。
- **FF 渐近**：(2.51)–(2.53) L842–850；实现 (3.75) L1038：`||𝓕0(s_i)||² ≲ 2m_q² ε^FF, ||𝓕1(s_i)||² ≲ ½ε^FF, s_i > s0`；ε^FF ≃ O(10⁻⁵) (L1051–1053)，取 **6×10⁻⁵**（L1141）。M=50 时 s_i>s0 的节点 7 个 → 14 条界。

### A3 目标泛函与图的生成
- Fig.3/4/8：投影到 (f00(3), f11(3)) 平面；径向扫描 L931："f00(3)=a+t cosα, f11(3)=b+t sinα … maximize t … Sweeping α∈[0,2π] gives the shape"。这是标准支持函数/凸包扫描；(a,b) 任取内点即可。
- Fig.9/10 代表点（L1148）："Previous experience with the bootstrap would suggest looking at the tip of the shape (here the red dot). However the results seem quite robust so we choose two other points (pink, light pink) near the chiral point (black dot)." 无算法；Fig.8 图上红点在 +x 尖端（≈0.085,−0.0053），两个粉点在黑点（0.0733,−0.00489）上侧。Fig.10 caption 自认 "similar but not equal phase shifts"。
- Fig.6/7：ε^χ=0.002 下"a point closest to the black dot"（L1111）。

### A4 求解器/算法
正文**没有任何**求解器、Newton、barrier、interior point、CVX、MOSEK、SDPB、principal value、regularization、midpoint 字样（grep 全文仅命中 [45] 标题 "regularization and dual convex problem" 与一般词 "feasib/toleran/semidefinite/convex"）。唯一算法性句子：L931 "Since the constraints define a convex space such shape is convex and can be mapped out by maximizing linear functionals"；L1033/L1049 容差"chosen such that … not too strictly to make the problem infeasible"。**被注释掉的 TeX**：L1337 `%%\bibitem{Tikhonov}`、L1576–1579 `%%\bibitem{cvx}`（Grant–Boyd CVX、Nemirovski conic programming）——是作者旧 bib 残留，未被正文引用；L853–857 注释掉一版 FF 容差写法 `|F0| ≲ 2m_q² ε^FF`（未平方）；L882–887 注释掉一版旧 (2.56) 数值（2.2e−6, 3.2e−4）；L897 注释掉的共振位置估计 E_P1≈820 MeV；L1231–1256 注释掉的 "pacman" 共形映射附录。

### A5 QCD 输入（§2.4）
s0=(1.2 GeV)², αs=0.4, m_u=4 MeV, m_d=7.3 MeV (2.54, L866)；⟨αs G²/π⟩≈0.023 GeV⁴，⟨jS⟩≈−(0.1 GeV)⁴ (2.55, L870)；μ=√s0；mπ=140 MeV（(2.2) 单位）；fπ=92 MeV 仅用于标记黑点（脚注 14，L888）。凝聚项"small…included for completeness"。无误差处理——都是中心值。

### A6 Appendix A（L1204–1208）
"M=50 interpolation points and L=10 partial waves per isospin … compromise between precision and speed." M=50 下 L=8,10,12 "reasonably good convergence"；L=10 下 M=45,50,60 "slightly more complicated … modifying M changes both the number of bootstrap variables as well as the unitarity and SVZ sum rules constraints … S0,S2 very good convergence, while the P1 … ρ resonance with the peak slightly shifted depending on M"。无数字；未说明 Fig.11 用哪个代表点。

---

## B. 论文确实没有唯一给出的东西

| # | 未定项 | 出处 | 判断 |
|---|---|---|---|
| B1 | 在割线上求分波时 (2.7) 的 1/(x−s) 奇异积分怎么离散（PV 核） | L929 只说 "evaluated using (2.9)" | (i) 常规。论文已为 FF 给出同一 Hilbert 核 (3.67)，沿用即可；作者 2403 notebook 亦用同一 cot 核 + Legendre-Q 角积分 |
| B2 | FESR 求和的 s0 截止（节点 ≤ s0 硬截止 vs 端点修正） | L1021 Σ_{i=1}^{M} 与 ∫₄^{s0} 不一致 | (i)。2403 代码 L14/L89 硬截止 `n0=sum(vnu<=s0)` |
| B3 | **ε^SR 的口径**（归一化/raw/相对） | L1030–1033 vs L876 | **(ii)**。字面归一化读法空约束 10³–10⁴ 倍；raw 读法四盒占目标 84%/1.8%/47%/1.4%（本审计算）；作者 2403 用归一化逐通道 L2 球，S0 ε=1e−7≈目标 20%，P1 ε=6e−6≈10–20%。这直接控制 P1 谱矩比即 ρ 质量 |
| B4 | **χ 范数打包**（两个 4 维球 vs 一个 8 维球） | L951 "some norm" | (i)–(ii) 之间：差别至多 √2 的有效容差，论文自称"small changes in this tolerance"不影响结果（L1079）；后续代码与仓库自己的 Fig.5 残差指纹（L2/ε≈1.000007…）都指向 8 维合并 L2。仓库当前却锁定 separate-L2，是自己选了更不像作者的口径 |
| B5 | PSD (3.70) 施加范围（全部 M 点还是 ≤s0） | L999 | (i)。s>s0 处 ρ_i 无 FESR 约束，PSD 平凡可满足 |
| B6 | FF 界中的 m_q 取值（算术平均/均方根） | L1038 | (i)。量级界，9% 无关紧要 |
| B7 | 代表点选择（tip + "near the chiral point"） | L1148 | (ii) 但被论文自身的三点散布界定（Fig.9/10 caption）；无 tie-break |
| B8 | 径向扫描内点 (a,b)、α 步长 | L931 | (i) |
| B9 | 是否用任何正则化/双谱密度上界 | 全文无 | (i)/(ii) 取决于是否活跃：2309 未提；2403 加了 ‖ρ/Mρ‖₄ ≤ 100 这种极松的截断（M=50 时 ‖ρ‖₄≤377500）。若在解处不活跃则无物理影响；仓库应报告其活跃性 |
| B10 | 求解器、精度、可行性容差 | 全文无 | (i)。凸锥规划，任何内点法给同一最优值到 1e−8 |
| B11 | Fig.11 各配置用哪个代表点、B 随 M 如何变 | L1206 | (ii) 小：论文只做定性比较 |

---

## C. 仓库术语逐项对照

| 术语 | 论文有无 / 出处 | 判定 |
|---|---|---|
| **M50/L10** | 有：L1065, L1206 | 论文规定 |
| **1500盘** | 有（推论）：3 同位旋 × 10 波 × 50 节点，(3.61)+(3.71) | 论文规定的直接后果；与 2403 代码 L65 一致 |
| **xref=.07332** | 有：(2.18) 在 s=3、fπ=92/140 → f00(3)=0.073321（本审计复算 0.0733214），即 Fig.4/6/8 黑点 | 论文规定 |
| **FF平方界 6e−5** | 有：(3.75) + L1141 ε^FF=6×10⁻⁵ | 论文规定 |
| **separate-L2 χ=.002** | ε=.002 有（L1079）；"separate" 无——(3.64) 两行各带 ‖·‖ 但注明 "some norm" | 仓库选择；与作者后续代码（8 维合并 L2）及仓库自己的 Fig.5 指纹相反 |
| **combined-l2** | 同上 | 更可能是作者口径；README 与 STATUS 对此互相矛盾（README 说当前 combined，STATUS/SCIENCE 说名义合同 separate） |
| **printed 四raw矩盒 .002** | "printed" 指用 (2.56) 印刷数字（合理）；"raw 绝对 .002 独立盒"是对 (3.72)–(3.74) 的一种字面读法 | 论文口径不唯一（B3）。这种读法使 S0 n=0 与 P1 n=−1 近乎空约束，与作者 2403 的相对 10–20% 球差异巨大，且直接影响 ρ 质量。是仓库最该优先重审的一项，而非"锁定不动" |
| **hard-midpoint** | 半有：(3.61) 中点节点 + (3.72) π/M 权是论文的；"masked to s_i≤s0"是常规截止 | 常规离散化，且与 2403 代码 L14 一致 |
| **clipped** | 无 | 仓库自创的端点修正变体 |
| **PV H / current** | 半有：FF 的 PV-Hilbert 核 (3.67) 是论文明写的；对振幅 (2.7) 在割线上必然需要同样的 PV 处理，论文未展开 | 论文表示的必然离散化，不是仓库额外引入；"H"只是把 (2.7)+(2.9) 在节点上写成线性映射的名字 |
| **B=377500 / B(M)=100[M²+M(M+1)/2]**（密度正则化） | **无**。2309 全文无 regularization。来源是作者 2403 代码 L10/L60 `norm(rho/Mrho,4)<=100`（Mρ=3775） | 从后续代码借来，非论文规定。论文的表示并不"必须"有它：目标与所有观测量都是节点上的线性泛函，被幺正盘有界；ρ_ij 的零空间方向只影响求解器数值稳定，不影响最优值。若仓库的解处该界活跃，则它改变了物理 |
| **infinity=free** | 有：T0 是 (2.8)/(3.62) 的自由变量 | 论文规定（free 正确） |
| **T0=0 / infinity=zero / 高能渐近条件** | 无 | 仓库自创（依附于 analytic-cardinal 族） |
| **analytic-cardinal / finite-sine** | 无。论文只在节点上工作，从未定义离节点的连续振幅 | 仓库自创的插值族；随之而来的"离节点幺正违例""残差恒等式""秩证明"都是自造问题 |
| **3123盘 / 3582盘** | 无 | 仓库自创：在 1500 原生盘之外加"ρ 窗口、全波交错、5 MeV 主波网格"附加幺正样本（STATUS "历史中间模型"段） |
| **Newton/barrier 支持求解、Arb 证书、analytic center、support gap** | 无 | 仓库自选的求解路线；论文/作者代码是 CVX+MOSEK 一次求解 |

---

## D. 忠实实现应当长什么样，以及责任划分

**问题规模**：实变量 ≈ 3876（振幅）+ 2×50（ImF）+ 2×50（ρ_current）≈ 4076；约束：1500 个二阶锥 |1+iπβf|≤1（或 |h|²≤2 Im h）、2×43（或 2×50）个 3×3 Hermitian PSD 块、4 条带容差的线性 FESR、14 条 FF 二次上界、1–2 个 χ 范数球、线性目标 t（径向）或直接 f11(3)。**这就是一个标准的 SDP/SOCP**，CVX/MOSEK、Clarabel、SCS 均可在秒–分钟级求解；作者 2403 快照 `GTB_numerics.m` L51–111 正是这样一段 60 行代码（`cvx_begin sdp; cvx_solver mosek; … norms([Reht,Imht],2,2)<=sqrt(2*Imhh); Bmat(:,:,i)>=0; norm(wS0)<=eS0; … maximize(F1); cvx_end`）。Fig.3/4/8 = 对 α 扫描若干十次求解；Fig.9/10 = 三次求解后读节点相位。整个论文的计算量是"一台笔记本一个下午"。

**论文的锅**：ε^SR 口径（B3）与 χ 打包（B4）确实不唯一，且前者物理上要紧；作者也承认 2309 当时无公开代码。仓库"作者 2309 完整 producer 尚未唯一恢复"这句话在**字面上成立**，但它把两个可以用作者后续代码 + 论文图指纹合理判定的口径问题，扩大成了对整套数值模型的不可知论。

**仓库的锅**：把一个凸锥规划换成自研的高精度 barrier/Newton 支持求解器、引入论文没有的解析插值族（analytic-cardinal/finite-sine）、离节点附加幺正样本（3582 盘）、T0=0 与高能渐近条件、以及从 2403 借来却当作"名义合同"锁死的 B 正则化与 raw-绝对 FESR 盒；随后用"不按相移选参数"的规则拒绝重审 B3——而 B3 恰是论文文本本身就不自洽、必须靠外部证据（作者后续代码给出归一化相对量级）来定的项。这些自造对象对应的"论文原物"本来就不存在，所谓"恢复"无从谈起。

**置信度约 80%**。最关键三条证据：(1) TeX L926、L993、L999–1011（变量集与约束集完整封闭）；(2) TeX L1030–1033 对照 L876–880（ε^SR=2e−3 与 (2.56) 归一化矛盾，本审计复算四矩 raw 值 2.39e−3/0.112/4.23e−3/0.146）；(3) 作者 2403 `GTB_numerics.m` L39、L60、L91、L102（归一化 FESR 球、8 维合并 χ、L4 正则化仅为后续截断、CVX+MOSEK 一次求解）。

建议的下一步（供上游决策，不在本审计范围内执行）：以 (2.56) 归一化矩为基准、按作者 2403 的相对量级（10–20%）重设四条 FESR 容差并用 combined-L2 χ 球，用现成锥求解器（MOSEK/Clarabel）在原生 1500 盘 + 100 Gram 模型上直接重做 Fig.8–10，再看 ρ 质量与上/下边界不对称是否复现；同时报告 B=377500 在各代表解处是否活跃。
