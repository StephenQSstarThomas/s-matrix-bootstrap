# 上游代码审计：He–Kruczenski 是怎样把 GTB 变成数值优化的

2026-09-12。只读审计；路径缩写：`UP=/home/shiqiu/s-matrix-bootstrap/references/upstream-gauge-theory-bootstrap/theories/qcd`，`NUM=UP/papers/arxiv-2403.10772/GTB_numerics.m`，`OPT=UP/papers/arxiv-2505.19332/src/matlab/optimize_core.m`，`DRV=UP/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m`，`T=/home/shiqiu/s-matrix-bootstrap/references/2309.12402v3.txt`，`SRC=/home/shiqiu/s-matrix-bootstrap/src/smatrix_bootstrap`。Mathematica notebook 已用脚本抽出 Input 单元到本目录 `nb2403_input.txt`、`nb2505_gen.txt`、`nb2505_plot.txt`，引用时标 cell 号。仓库 git：`main` 与 tag `arxiv-2505.19332` 的两个 .m 逐字节相同（diff 为空）；tag `arxiv-2403.10772` 只含一个 219 行的 `GTB_numerics.m` 与生成核的 `.nb`。

## §3 判断（先给结论）

1. **按上游写法，这个问题就是一个标准的 SOCP+小块 SDP，MOSEK/Clarabel/SCS 可以直接解。** 锥类型只有：1500 个三维旋转二阶锥（幺正性）、117 个 3×3 Hermitian PSD（电流 Gram）、一个 3775 维 L4 范数球、两三个 8–12 维二阶锥（手征）、若干线性/绝对值约束（FESR、FF 渐近）。`NUM:51-111`、`OPT:28-124` 全部用 CVX 建模、`cvx_solver mosek`，没有任何 `cvx_precision`/`cvx_solver_settings`，即完全交给默认内点法。本仓库自己的 `SRC/conic.py:12-58` 也已把同一联合问题编译成 Clarabel 的 SOC/PSD/L4 锥，说明规模和锥类型上没有障碍。
2. **自研 barrier+Newton+QR/CG+Phase I+中心证书路线不是论文要求，属于"等价替代但自找麻烦"。** 2309 §3 只说"maximize linear functionals"；2403/2505 用现成 MOSEK；第三方 2511.11513 同样 CVX+MOSEK 且明说"用 SDP 求解器即可"。Newton 本身与 MOSEK 的内点法在数学上同源，但本仓库的实现引入了论文没有的自由度：不精确方向门槛 1e-6（`SRC/linear.py:111`）、中心判据 1e-14（`linear.py:11`）、支撑 gap 1e-4–1e-3 才停（STATUS 记录 gap .00113/9.2e-5）、以及"中心/代表点"规则（`linear.py:201-206`）。这些都可能在不改物理约束的情况下改变最终代表幅度，进而改变 90° 读数。**ρ 质量 46 MeV 的 L 分散不能直接归咎于求解器**，因为模型本身（χ 范数、B、FESR 误差、s0）也与上游不同；但求解精度松于 MOSEK 默认 1e-8 量级这一点是可核查的额外系统误差来源。
3. **最省事且最可信的路线**：以 `NUM` 的约束组装为模板（变量、锥、目标逐行对应），用 cvxpy→MOSEK（有 license）或 Clarabel（tol 1e-9）改写，把参数换回 2309：ν0=0、s0=1.2 GeV、αs=0.4、两个电流 S0/P1、(3.70) 三阶 PSD 只作用于 S0/P1、其余波用 |S|≤1、手征两条 separate 4 维范数、FESR n=0,1/−1,0 共 4 矩、FF 用 (3.75) 平方界；目标改为 §3.1 的径向最大化 t。核函数最好直接移植 2403 nb 的构造（cell 11–40，参数 nu0 可设 0），而不是复用本仓库另一套 PV 离散化。Newton 代码保留为独立校核器，不做主线。

## 1. 作者数值管线的精确描述

### 1.1 变量与参数化
- 振幅变量：`sig(Msigma=2M+1=101)` = (T0, σ1,i, σ2,i)，`rho(Mrho=M²+M(M+1)/2=3775)` = (ρ1,ij 全矩阵, ρ2,ij 上三角)；`NUM:10,53-54`、`OPT:31-32`。生成器对 ρ2 用 `2*symtr` 再把对角列除 2、取上三角（nb2403 cell 32），所以向量存的是 actual ρ2 上三角。
- 网格：M=50 个上半圆点 φ_i=(i−½)π/M，映射 s_j=ν0+(8−2ν0)/(1+cos φ)，ν0=−20（`NUM:7,11-12`；`DRV:26,33-35`）。每同位旋 l=10 个分波（L=20，nb2403 cell 2-3）。s0=(2000/139.57)²≈205.3，n0=39 个节点 ≤ s0（`NUM:13-14`；本机复算）。
- 分波矩阵：`Imht=B2*rho+B1*sig; Reht=A2*rho+A1*sig; Imhh=Bhat2*rho+Bhat1*sig`（`NUM:63-64`），B/A 为 1500×3775、1500×101 的稠密预计算核；核用 Mathematica 100 位精度（cell 4 `p=100`）计算，包含直接道的离散 cot 型 PV 核（cell 11，与 2309 式 (3.67) 相同）、交叉道 W/W00 有限和（cell 12-14）与 Legendre-Q 角投影（cell 15-16）。行序为 I=0,2,1 分块（cell 38-39）。
- 电流部分：`ImFF(3M=150)`；`ReFF=KFF*ImFF+1`（`NUM:78`；KFF3 = K + K0 的减除修正，cell 59-60）。谱密度 2403 用 40 个三次 BSpline 系数 `specC(120)`，`specC>=0`、首系数==1（`NUM:74-75`；cell 61-63）；2505 改为 49 个解析系数 `a(147)`，Π(z)=(1−z)^{k}(a0+(1−z)Σ a_n z^{n−1})（nb2505_gen cell 76），约束 `abs(a(2:Na))<=aXb`、`a(1)∈[1,10]`（`OPT:65-67`）。
- 辅助变量 `Bmat(3,3,3*n0) hermitian`（`NUM:57`）。

### 1.2 幺正性
`norms([Reht,Imht],2,2) <= sqrt(2*Imhh)`（`NUM:65`；`OPT:45`），逐行即逐（I,ℓ,s_j）一条，共 3·10·50=1500 条；等价 |h̃|² ≤ 2 Im ĥ，CVX 编译为旋转二阶锥。h̃=h/Λ、Im ĥ=Im h/Λ²，Λ_ℓ(s)=((√s−2)/(√s+2))^{ℓ/2}（cell 54-56），即先解析重缩放再交给求解器。没有 2×2 PSD 写法。

### 1.3 手征约束
`fchi=fsigmachi*sig+frhochi*rho` 给出 0<s<4 处 4 个点 s*=½,1,3/2,2 的分波（`NUM:25,97`；cell 7、52）。2403：`norm([S0−R01·P1; S2−R21·P1]) <= echi`（8 维**合并** L2）另加 `norm([D0;D2;F1]) <= echi`（12 维）；`echi=2e-3`（`NUM:38,102-103`）。2505 只保留 8 维合并一条（`OPT:109`，`DRV:48`）。2309 (3.64) 写的是两条**分开**的范数（`T:952-954`）。

### 1.4 电流 Gram PSD
对 S0、P1、D0 三条波、每个 s_j≤s0 的节点（3·39=117 块）：`Bmat(:,:,i) == [1, h̃, F; conj(h̃), 2Imĥ, 2ImF̂; conj(F), 2ImF̂, spec]; Bmat>=0`（`NUM:83-86`；`OPT:78-83`），ImF̂=ImF/Λ3（`NUM:80`）。2309 (3.70) 只对 S0、P1 两条波施加（`T:1015-1027`）。

### 1.5 FESR/SVZ
2403：每通道 3 矩 `fesr=[intS0*spec;intP1*spec;intD0*spec]`，`norm(w_X) <= e_X`（3 矩组内 L2、绝对误差），eS0=1e-7、eP1=6e-6、eD0=5e-6（`NUM:39,89-91`）；pQCD 值由 nb cell 67-69 给出（αs(2GeV)=0.3146）。2505：6 矩、pinched 权重、先把矩缩放到 O(1)（SRscale 1e7/1e5/1e5，nb2505_gen cell 70-73），再 `abs(w)<=e` 逐矩盒式，e=0.05/0.1/0.1（`DRV:45`；`OPT:86-93`）。2309：S0 用 n=0,1、P1 用 n=−1,0 共 4 矩，`||·||≤ε_SR`（`T:1035-1053`），s0=1.2 GeV、αs=0.4。
FF 渐近：2403 `abs(FF_S0(s>s0))<=0.05`、`|s·F_P1|<=2·6.87`、`|s·F_D0|<=6·20.6`（`NUM:40-41,94`，11 个高能点）；2505 S0 改为 (s/log s)|F|≤8、P1/D0 用 Brodsky–Lepage 常数×2.2/8（`DRV:21-23,46`；`OPT:96-98`）。2309 (3.75) 是 |F|²≲ε_FF~1e-5 的平方界（`T:1068-1087`）。

### 1.6 目标、扫描、代表点
- 2309 §3.1：在 (f00(3), f11(3)) 平面上取内点 (a,b)，加约束 f00=a+t cosα、f11=b+t sinα，**最大化 t**，扫 α 得边界（`T:918-928`）；代表点为边界上靠近手征点的红/粉/浅粉三点（`T:1500-1503`）。
- 2403：`F0=2(f00+5f02)|_{s=3}`、`F1=2(3f11+7f13)|_{s=3}`；`F0==F0fix; maximize(F1)`，F0fix=ni/50·0.163108，ni=44/46/48 对应蓝/绿/红（`NUM:43-45,106-109`）。之后**一次** Watsonian 再解：目标 `maximize(v1+v2)`，v = 上一解方向·h − ΣImĥ，S0/P1/D0 的方向取自形状因子相位（`NUM:117-120,178-180`）。
- 2505：ni=0 时 `minimize(0)` 纯可行性，然后 5 次 Watsonian 迭代，方向 `|h_prev|·F_prev/|F_prev|`（`DRV:50,91-101`；`OPT:116-122`）；每步输出 λ=3πA(4/3,4/3)/4（`OPT:112-113`）。

### 1.7 求解器、规模、耗时
`cvx_begin sdp` + `cvx_solver mosek`（`NUM:51-52,122-123`；`OPT:28-29` 加 `quiet`）。无精度设置；README 要求 CVX + MOSEK 学术许可（`UP/../README.md` Requirements 节）。规模：原变量 3876+150+120(或 147)，辅助 Hermitian 117×9 实自由度，L4 锥 CVX 内部再展开约 2×3775 个三维锥；约束 1500 SOC + 117 PSD + 3 SOC(3) + 33 abs + 8/12 维 SOC + 120 非负。核矩阵 B2/A2/Bhat2 各为 1500×3775 稠密阵，是内存与每次内点迭代的主要开销；第三方改用离散基（nmax=37，变量数 1+5·nmax(3nmax+1)/2≈1.0×10⁴，`cordoba.txt:810-812`；其脚注 `:779-781` 也明说 ρ(x,y) 的范数正则化"is necessary"并引 2103.11484；L4 锥只作用于两个 nmax² 维向量，`f_00f_11ExperimentsSVZcircle.m:361-362`）后报告 30 分钟→1 分钟，但论文未剖析提速的具体来源，这里只记录事实。README/注释没有给耗时；第三方 2511.11513 称 2403 的实现"roughly 30 minutes per run"，其自己的离散基版本约 1 分钟（`scratchpad/cordoba.txt:1944-1946,2347-2349`；[2]=2403 见 `:2361-2362`）。

### 1.8 PV、正则化、无自研求解器
- PV：只在核生成阶段以离散 cot 核处理直接道 Hilbert 变换（cell 11 "kernel for principal integral eq. (8.36)"），优化阶段全是线性代数。
- 正则化：`norm(rho/Mrho,4) <= 1e2`（`NUM:60`；`OPT:38`），即 ‖ρ‖₄ ≤ 377500。**2309 正文没有任何此类约束**（`T` 中 "regulari" 只命中参考文献 [45] 标题 `T:1915-1916`，[45] 在 `T:447` 只用于一个对偶界）。
- 无 Newton/barrier/CG 自研代码：两个 .m 全文只有 CVX 块、索引整理、绘图与 JSON 输出（`NUM:1-219`；`OPT:1-173`；`DRV:1-162`）。

### 1.9 后处理
2403：`deltaS0=angle(h̃)·deg`、`deltaD2=½angle(1+i h̃)·deg`（`NUM:190-192`；注意 `NUM:192` 的 deltaP1 漏乘 deg，是原代码小 bug），只画相移，不算 η、不找 ρ 质量。2505 绘图 nb：S=1+iΛh̃，η=|S|，δ=arg S/2（S0,S2,D2,F1）或 arg h（D0,P1）加 +180° 分支修正（nb2505_plot cell 3-4）。**上游没有任何"找 90°"或"找极点"代码**；2309 文字说 ρ 位置取"phase shift crosses π/2"（`T:1508-1510`）。

### 1.10 2403/2505 相对 2309 的改动
ν0：0→−20；s0：1.2→2 GeV；电流：2→3 条（加 D0）；谱表示：逐点 ρ_ℓ,i→BSpline(40)→解析系数(49)；FESR：4 矩 L2→9 矩组内 L2→18 矩逐矩盒；手征：两条分开 4 维→合并 8 维(+12 维)→合并 8 维；FF：平方界→绝对值×渐近常数；目标：径向 t→固定 F0 最大化 F1+一次 Watson→可行性+5 次 Watson；新增 L4 正则化。因此 2403/2505 代码只能当**约束组装与求解器接口的模板**，不能当 2309 的参数表。

## 2. 对 2309 的推断
2309 §3 的变量集 (3.69) 与 2403 完全同构（T0、σ、ρ、ImF、ρ_ℓ；M=50 时 3876+100+100=4076 个实变量），约束全部是 (3.70)/(3.71) 的 3×3/2×2 PSD、(3.64)/(3.74)/(3.75) 的范数球和线性等式，目标是线性泛函 t 的最大化（`T:898-928,1004-1090`）。Appendix A 说 M=50、L=10 是"precision and speed 的折中"并做了 M=45/50/60、L=8/10/12 的变参（`T:1651-1670`）。2403 是同一作者、同一 Mathematica→MATLAB 管线的首次公开，TeX 源里还有被注释的 CVX 书目（`2309.12402v3-source/prd_submission_2.tex:1576-1577`）。综合判断：2309 的数值问题几乎必然是**一个 ~4×10³ 变量、~1.5×10³ 个 SOC 加 ~10² 个 3×3 PSD 的锥规划，用 CVX+MOSEK（或同类内点法）求解，单次求解从分钟到数十分钟量级**（Cordoba 对 2403 的 30 分钟报告是唯一外部计时；2309 无 L4 锥、PSD 更少，应不慢于此）。理由补充：(a) 2309 的 §3.1 要"sweeping α∈[0,2π]"画出 Fig.3/4/8 的整条边界，每个 α 一次求解，图上边界点数以十计到百计，若单次求解需数小时则整篇文章的多张变参图不可能完成，反过来说明单次求解应在分钟量级；(b) 所有约束都是范数球或 2×2/3×3 PSD，没有任何非凸或需要自定义内点法的结构，2309 又反复强调"convex space … maximizing linear functionals"，这正是 CVX/MOSEK 这类工具的设计场景；(c) 同作者 2021 年的 2103.11484 已在 3+1 维 S-matrix bootstrap 中使用同类锥规划与对偶方法，2309 在 `T:447` 直接引用其对偶结果，工具链连续。按 2309 参数复算：ν0=0 的映射下 s0=(1200/140)² 对应 n0=43 个节点，(3.70) 的 3×3 PSD 块为 2×43=86 个，(3.71) 的 |S|≤1 共 1500 条，FF 高能界 2×7=14 条，FESR 4 条，手征 2 条——这比 2403 的问题还小。不确定性：2309 是否已用 L4 正则化、χ 范数是否真为两条分开、FF 界的确切缩放、ρ_ℓ,i 是否加了非负或光滑性约束，代码未公开，无法证实；上述"分钟量级"也只是从边界扫描的可行性反推，不是实测。

## 3. 对本仓库路线的详细判断
- **可直接用现成锥求解器**：`SRC/conic.py:38-54` 已把散射 SOC、手征 SOC、密度 L4（SOC3+SOC(nd+1) 精确提升）、Gram PSDTriangle(3)、FF SOC3、矩盒 NonnegativeCone(8) 全部编译给 Clarabel/SCS；`conic.py:109-111` 设 tol 1e-11/1e-10。仓库记录 D3_joint_001/002 用 Clarabel 约 320 s 达 MaxTime 无解（`results/runs/mainline_alignment_20260909/SOLVER_DECISION_ZH.md:26`），但同文 :36 承认没有同题 MOSEK 对照。上游在交给 MOSEK 前做了 Λ 重缩放（cell 54-57）和电流 Gram 换基（`NUM:80`），仓库的 `operators.py:52-57` 也做了对角合同缩放；Clarabel 失败更像 conditioning/超时设置问题，不是锥类型或规模超界。
- **必要 / 等价 / 自找麻烦**：不是必要；数学上等价（同一凸集、同一线性目标）；工程上自找麻烦。`linear.py:100-111` 的 QR 复用 + ≤64 步 CG、`:111` 的 1e-6 不精确方向容忍、`:132-143` 的自制线搜索、`:181,205` 的 μ 减半续算，替代了 MOSEK 的同质自对偶内点法，且没有 MOSEK 那样的 primal-dual gap 终止；仓库用外部 `joint_audit` 支撑界补证（`linear.py:73,172-190`），gap 只做到 1e-4–1e-3（STATUS 表格）。单次求解 108–353 s（`results/runs/physical_consistency_20260911/PV_UV_mid_support_02/report.json`、`PV_tip_support_01/report.json` 的 solver.seconds），并不比上游快。
- **系统偏差风险**：(i) 代表点规则不同——上游代表点是边界上的线性泛函极大点（2309 径向 t；2403 固定 F0 极大 F1），仓库则把"收敛中心"作代表（`linear.py:201-206` `require_center`），中心随 μ 路径和 barrier 权重（`linear.py:49-52` chiral_barrier_weight）变化，这是论文没有的自由度；(ii) 松 gap 下的幅度未必是该目标的极大点，P1 相位穿越 90° 的位置对幅度小变化敏感，可能贡献 ρ 读数分散；(iii) 但 STATUS 的 46 MeV L 分散也可能来自 B=377500（借自 2403，2309 未提）、χ 范数选择、FESR 误差口径等模型差异，不能仅凭本审计归因于求解器。
- **为什么说"中心"与论文的代表点不是一回事**：内点法在 μ→0 极限才逼近目标的极大点；仓库若在 μ 仍有限时以"收敛中心"作代表，取到的是被 barrier 权重推离边界的内点，其 f00(3)、f11(3) 与论文的边界极大点系统性不同。2309 图 8 三个代表点全在边界上（`T:1500-1503`），2403 的三个代表点也是"固定 F0、极大化 F1"的边界解（`NUM:44-45,109`）。仓库的 `require_center` 开关（`linear.py:29,202,206`）把两种规则并存，报告时必须分清哪一种在用，否则与论文相移曲线的比较没有意义。
- **Watsonian 步骤的定位**：`NUM:114-182` 的第二次求解与 `DRV:91-101` 的 5 次迭代是 2403/2505 新增的幺正化选点手段，目的是从边界解中挑出接近饱和幺正性、且 S0/P1/D0 相位与形状因子相位一致（Watson 定理）的幅度；2309 正文没有这一步。忠实复现 2309 时它只能作为可选后处理，不能当作必做约束；仓库 STATUS 已注意到这一点，本审计确认其正确。
- **对 Clarabel 超时的看法**：上游把幺正锥按 Λ_ℓ(s) 重缩放（cell 54-56）、把 Im F 除以 Λ3 后再进 PSD（`NUM:80`），并把谱密度写成非负 BSpline 系数（`NUM:74-75`），这些都是在建模层面改善条件数的做法；仓库 `conic.py:19` 只对密度列做了 `bound/nd^0.25` 缩放，散射行靠 `operators.py:37-39` 的 t 缩放。在没有同题 MOSEK 对照、也没有按上游方式缩放后再试 Clarabel 的情况下，"现成求解器不行"这个结论并未成立。
- **建议路线**：见 §3 判断第 3 条。具体做法是照 `NUM:51-111` 逐行翻译成 cvxpy（变量→`cp.Variable`，`norms(...)<=sqrt(2*Imhh)`→`cp.SOC`/`cp.quad_over_lin`，`Bmat>=0`→`cp.PSD` 或 3×3 Hermitian 实化，`norm(rho,4)`→`cp.pnorm`），核矩阵优先移植 2403 nb 的构造并把 nu0 设 0、s0 设 (1200/140)²、L=20；用 MOSEK 或 Clarabel 求解，先复现 2403 自身（有公开参数可对拍），再切 2309 参数做径向扫描；Newton 代码只留作对同一冻结问题的独立复核。

## 4. 第三方旁证（2511.11513）
论文 §4："such coefficients can be found, for instance, using semi-definite programming solvers [29–32]"（`cordoba.txt:785-786`；[29]=MOSEK 手册、[32]=CVX，`:2420,2426`）。代码 `cordoba-discrete-gtb-2511.11513/f_00f_11ExperimentsSVZcircle.m:302-304` `cvx_begin quiet; cvx_solver mosek`；`:355` `maximize(f2S0)`；`:361-362` `norm(x1,4)<=Mreg*10^power`（=1000）；`:366-367` 二次锥幺正性；`:371` 合并 8 维手征范数 2e-3；`:395-400` `hermitian_semidefinite(3)`；`:407-409` SVZ 3 矩 L2；参数 `:94-114`。即社区两个独立实现都用 CVX+MOSEK，没有人自写内点法。
