# 仓库实现 vs He–Kruczenski 2309.12402v3：只读审计

审计对象：`/home/shiqiu/s-matrix-bootstrap/src/smatrix_bootstrap/*.py`（23 模块 6033 行）、`tests/`（3 文件、62 个函数、`pytest --collect-only` 收集 107 项）、`SCIENCE.md`/`REPRODUCTION_GUIDE_ZH.md`/`STATUS.md`、论文 `references/2309.12402v3-source/prd_submission_2.tex`（下文 TeX 行号均指此文件）。未修改任何仓库文件，未运行生产计算。

## 0. 结论先行

1. **仓库名义主线（"PV-midpoint M50/L10, 1500 盘"）求解的是论文 §3 的同一个凸可行集的一个具名有限离散版本**：变量集 (3.62)+(3.69) 完全一致（3876+200=4076），幺正盘 (3.71)、S0/P1 三阶 Gram (3.70)、四点手征残差 (3.64)、四个 FESR 矩 (3.72–3.74)、高能 FF 平方上界 (3.75) 逐项对应。散射色散积分的离散规则（中点 φ 网格＋节点 PV/跳跃＋解析 Legendre-Q 角投影）是论文只给了形状因子版本 (3.66–3.67) 的自然推广，属合理离散化。
2. **但有四项决定可行集大小/代表点身份的东西不是论文规定的**：(a) 双密度 L4 正则化 B=377500 —— 从作者 2024 年后续代码移植，2309 全文无此约束，仓库自己的扫描显示结果随 B 单调变化、无平台；(b) 手征范数取 L2 球（separate/combined）—— 论文 (3.64) 字面是 8 条标量不等式（L∞ 盒），L2 球严格更小；(c) FESR 误差取"raw 绝对 .002"—— 在归一化单位下四个盒的相对宽度分别约 84%/1.8%/47%/1.4%，这个口径直接决定 ρ 的能标；(d) 代表点取"加权障碍（χ 权重 750）在小 μ 处的解析中心"而不是极值点 —— 论文的点是求解器返回的边界极值。
3. **Newton 路线不是论文要求的，也不是必要的。** 论文 §3.1（TeX 931）只说"凸集合、最大化线性泛函"，后续作者代码用 CVX/MOSEK。仓库问题是标准 SOCP+小块 SDP（≈4076 变量、1500 个 SOC(4)、100 个 PSD(3)、2 个 SOC(5)、14 个 SOC(3)、8 条线性、1 个 L4 球），现成内点求解器可解。仓库放弃 Clarabel/SCS 的依据是 300–600 s 的时限内未完成、且两次全空间运行被 SIGTERM 杀掉（无任何迭代记录）。之后为自研 barrier-Newton 又补了 Phase I、QR/CG 残差门槛修复、近中心对数抵消的 Arb 下降证书、径向居中等 ≈1800 行代码。这是典型的"没苦硬吃"——不过它附带的 Arb 对偶证书本身有价值，且可以直接接在锥求解器的对偶上（仓库已在 hull 路径这么做）。
4. 3123/3582 盘、analytic-cardinal、T0=0、五条渐近不等式是仓库自建的**更强诊断模型**，不是论文问题；仓库文档已如实标注为历史/诊断分支。

---

## 1. 仓库实际求解的问题（数学表述）

### 1.1 决策变量：3876 + 200 = 4076

`kernels.density_labels` (kernels.py:99-103)：`T0`(1) + `sigma1_i`(M) + `sigma2_i`(M) + `rho1_ij`(M²，不对称) + `rho2_ij, i≤j`(M(M+1)/2)。M=50 → 1+50+50+2500+1275 = **3876**，即论文 (3.62) 的 {T0, σ_{α,i}, ρ_{α,ij}}，ρ₂ 对称（TeX 926-929）。电流侧 `model.current_operators` (model.py:112, 155-156)：ImF0(50)+ImF1(50)+ρ_cur0(50)+ρ_cur1(50) = 200，即论文 (3.69) 的 {ImF_{ℓ,i}, ρ_{ℓ,i}}。**判定：论文规定。** 存储约定 `C_flat_ij = 2ρ2_ij (i≠j)`（kernels.py:148-151）只是打包。

### 1.2 参数化与"基"

论文 (2.7) 是 Mandelstam 双谱色散表示，(3.58)-(3.61) 把切割映到单位圆、在上半圆取 M 个中点 φ_i。仓库 `operators.midpoint_grid` (operators.py:207-214)：x_j = 4/cos²(φ_j/2)，w_j = x'(φ_j)/M（含 (2.7) 的 1/π）。**这里没有"多项式基"**：主线 `PVSourceRows` 的未知量就是节点密度值，色散积分用中点求积。`basis.py` 的 `CardinalSourceRows`（basis.py:52-208）是另一条分支：把节点密度经正弦-cardinal 变换 `analytic.cardinal_transform` (analytic.py:16-21) 变成 z^n 多项式系数，得到一个真正的解析函数族 `CardinalAmplitude` (analytic.py:51-130)。这条分支是仓库为了诊断 PV 离散在节点间是否解析而自建的，论文没有。

### 1.3 色散/交叉的实现："PV 核"是什么、为什么需要主值

论文只说"用 (2.9) 计算分波"（TeX 929），没写离散密度下 s-道 Cauchy 积分 (1/π)∫σ(x)dx/(x−s) 在切割上 s=x_k+i0 怎么算。仓库采用（`kernels.physical_node_row`, kernels.py:136-143）：

- 直道（s 在节点 x_k 上）：Re 部分 = 离散 Hilbert 核 K_{kj}（`pv_matrix`, kernels.py:105-113，即论文 (3.67) 的 K̃ cot 公式；test_kernels.py:31-38 验证它把 sin(nφ) 映到 cos(nφ)），Im 部分 = i·δ_{kj}（节点跳跃 = 密度本身）。这就是"主值＋半留数"，是 Sokhotski–Plemelj 的离散版；**主值是需要的**，因为在切割上 Cauchy 核有极点，论文自己对形状因子 ReF 就用同一个 K（(3.66)-(3.67)，TeX 970-985）。仓库把它推广到散射密度，并用 `COLLOCATION_KERNEL_MAP_ZH.md` 核对了与作者 2403 公开代码一致。
- 阈下（0<s<4，用于手征点 s=1/2,1,3/2,2 和目标点 s=3）：`offcut_row` (kernels.py:128-134) 直接用有理核 w_j/(x_j−s)，无需主值。
- 交叉道 (t,u)：角积分 ∫P_ℓ(μ)/(x−t)dμ 解析求出，`operators.angular_kernels` (operators.py:216-242) 用外 Legendre Q 函数 `exterior_legendre_q` (kernels.py:115-126)；双密度交叉块 U_ij 用部分分式恒等式。test_kernels.py:14-30, 56-78 对 mpmath 直接积分做了独立校验。

**判定：合理离散化**（与 (3.58)-(3.67) 一致的最直接实现）。仓库文档指出节点 PV 值与阈下有理函数不是同一解析函数（留数 −w_i Im f(s_i) ≠ 0），这是有限离散的固有性质，论文实现同样存在；仓库为此另建 analytic-cardinal 族，不应算入主线。

### 1.4 幺正性的离散形式："1500 盘"、"hard-midpoint"

`canonical_waves` (kernels.py:162-163)：I=0,2 取 ℓ=0,2,…,18，I=1 取 ℓ=1,3,…,19，L=10 → 30 波；50 节点 × 30 = **1500** 个物理样本，每个是 2×2 PSD ⇔ 圆盘 2Im h − |h|² ≥ 0，h=κf（`model.unitarity_margin`, model.py:45-50）。这正是论文 (2.12)/(3.71)，M=50、"10 partial waves per isospin"（TeX 1068-1071）。其中 S0/P1 的 100 个样本进入 3×3 Gram (3.70)（`current_operators`, model.py:137-145，实对称化后 6 分量 svec；test_mainline.py:55-89 用复 U†BU 独立核对）。**判定：论文规定。**

"hard-midpoint" 与幺正无关，是 **FESR 求积截止**（`quotient._fesr_grid`, quotient.py:71-83）：权重 (π/M)·s'_i·1[s_i≤s0]，M=50 时 43 个低能节点计入、7 个高能节点进 FF 上界。"clipped-phi" 是仓库的端点修正变体。论文 (3.72) 的求和上限印为 M 但积分到 s0，未写掩码规则，所以 hard 是最字面读法。注意仓库自己的测试（test_linear.py:27）断言：常数谱下 hard 的 n=1 矩误差在 28–30%——这是一个已知且被保留的离散系统误差。

### 1.5 手征约束："separate-L2 χ=.002" / "combined-l2"

残差行（kernels.py:69-77）：r01(s_j)=f00(s_j)−[3(2s_j−1)/(s_j−4)]f11(s_j)，r21(s_j)=f02(s_j)−[3(2−s_j)/(s_j−4)]f11(s_j)，s_j=1/2,1,3/2,2 → 8 个实数。`model.chiral_slices` (model.py:16-19)：
- separate-l2：‖(r01)_{j=1..4}‖₂ ≤ ε 且 ‖(r21)_{j=1..4}‖₂ ≤ ε（两个 4 维球）；
- combined-l2：‖(r01,r21)‖₂ ≤ ε（一个 8 维球）。

论文 (3.64)（TeX 938-953）逐 s_j 写两条不等式，"with some norm"。字面读法是 8 条 |r|≤ε（L∞ 盒，LP 约束）；任何 L2 球都严格包含在该盒内。仓库两种选择都比字面读法**更紧**，且 combined 的依据是从 Fig.5 曲线反推的"指纹"（`SCIENCE.md` 自称 output-assisted）。当前合同用 separate。**判定：仓库选择，论文未定；影响 Fig.4/5 区域与代表点。** CLI 强制显式声明（model.py:21-27）。

### 1.6 FESR："printed 四 raw 矩盒 .002"

`fesr_targets` (quotient.py:54-69)：printed 分支直接取 (2.56) 的系数 J_n：S0 n=0: 3.09e−8(27.38/2+.61)，n=1: 3.09e−8·27.38/3；P1 n=−1: 4.34e−6·13.26，n=0: 4.34e−6(13.26/2−.41)；raw 目标 = J_n·s0^{n+2}，s0=(1.2/.14)²。约束（`current_operators` model.py:146-147；`JointProblem` operators.py:58-59；`joint_audit` certificates.py:154）：|Σ_i w_i s_i^n R_{ℓ,i} − J_n s0^{n+2}| ≤ 0.002，四条独立。论文 (3.74) 左边是 raw 积分、右边是"QCD value"、误差 ε_SR=2×10⁻³（TeX 1141），但 (2.56) 给的是 s0^{−n−2} 归一化矩，量纲不匹配，仓库取 raw 绝对误差为主读法，`--sr-error normalized-absolute` 为备选（`UVConfig`, __init__.py:37）。换算到归一化单位，四个容差约 (3.7e−7, 5.0e−9, 2.7e−5, 3.7e−7)，相对目标 (4.4e−7, 2.8e−7, 5.8e−5, 2.7e−5) 分别是 **84%、1.8%、47%、1.4%**——两个 n 较大的矩几乎被钉死，两个 n 较小的矩几乎自由。这直接决定 P1 谱的均方能标（仓库 DERIVATION §11 自己推出 ≈700–800 MeV）。**判定：论文口径未定；仓库读法字面合理但对 ρ 质量极敏感。**

### 1.7 形状因子 Gram 与高能上界

Gram：100 个实 3×3 PSD 块（2 通道 × 50 节点），条目 [[1+ReS, ImS, √2 k ReF],[ImS, 1−ReS, √2 k ImF],[…, R]]，F = 1 + K·ImF + i·ImF（(2.37)+(3.66)），k² 由 (2.33) 得 k0² = 3β/(256π⁵)，k1² = (s−4)β/(384π⁵)（`current_kinematic_squares`, model.py:56-58；与 (2.42) 一致）。FF 上界（model.py:148-154）：s_i>s0 的 7 个节点上 |kF|² ≤ 2m_q²ε_FF（S0）、ε_FF/2（P1），ε_FF=6e−5，m_q=(4+7.3)/2 MeV（`--mq-rule`）。**判定：论文规定 (3.70)、(3.75)；m_q 取平均是仓库补充的约定。**

### 1.8 密度正则化 B=377500，B(M)=100[M²+M(M+1)/2]

`prepare_amplitude` (operators.py:280-281)、`barrier_terms` (quotient.py:14-33)、`support_outer` (imaginary.py:202-208)：‖(ρ1_ij 全部, ρ2_ij i≤j)‖₄ ≤ B，M=50 时 d_M=3775，B=100·d_M=377500。**2309 v3 全文没有任何密度范数约束**（tex 中"regulari"只出现在被注释的书目）。它来自作者 2403 公开 MATLAB `norm(rho/Mrho,4)<=100`（仓库 `REGULATOR_METHOD_ZH.md` 逐行追踪）。为什么需要：离散色散参数化中密度可以在节点上大幅振荡而分波变化很小，无正则化时优化会把系数推到数值极端；2021 年 Kruczenski 等人的论文引入它并要求"随 B 增大寻找平台"。仓库自己的扫描（`COMMON_REGULATOR_DIAGNOSIS_ZH.md`）表明：五个手征方向在两个数量级窗口内**都排除 1% 平台**，pure +x 只有 B∈[1e3,1e5] 有平台；B 越大区域越大（+x 从 .0864 到 .1059，论文 Fig.4 右端约 .0826）。**判定：仓库移植（从后续工作）且可能改变物理**——它是可行集的一部分，而论文没有它。

### 1.9 目标函数：支撑函数、中心、+x tip、xref

- 支撑函数：`boundary` 命令最大化 d·(f00(3), f11(3))，H 的最后两行是 s=3 的目标行（kernels.py:80-82；`scattering.barrier_support` scattering.py:39）。论文 §3.1（TeX 931）从内点 (a,b) 沿角度 α 最大化 t。对凸集两者等价地刻画边界；不同处只在平坦边界段返回哪个点。
- "+x tip"：d=(1,0)，即 max f00(3)，对应论文 §4.3"tip of the shape (red dot)"（TeX ~1160）。
- "xref"：`analysis.XREF = 5/(16π²(92/140)²)`（analysis.py:10）= Weinberg 树级 f00(3) 在 fπ=92 MeV 的值（`weinberg_waves`, model.py:29-31，即 (2.18)），也就是论文 Fig.4/8 黑点的 x 坐标；y_ref=−xref/15。ref 点 = 在 x=xref 截面上最大化 y（`--fixed-x`，scattering.py:41-43, 87）；mid = 在 x=(x_tip+xref)/2 截面最大化 y（`PV_REPRESENTATIVE_RULE.json`；`select_gauge` 的 nearest 规则 analysis.py:231-239）。论文只说"two other points (pink, light pink) near the chiral point (black dot)"，无规则。**判定：仓库预注册规则，合理但非论文；无害的前提是承认"同一投影点对应的完整振幅不唯一"。**
- "中心"（`--require-center`, run.py:266；linear.py:202-206；scattering.py:221-224）：返回的不是极值点，而是加权障碍函数在最终 μ（1e−7…1e−10）处的收敛驻点（"数值中心"），χ 球权重默认 n_disks/2=750（ir.py:145、linear.py:50）。投影到 2D 后最优面是高维的，返回哪个完整振幅取决于障碍权重。论文用 MOSEK 得到的是（无权重）内点法解。**判定：仓库自创，改变所选代表振幅，从而改变相移。**

---

## 2. 求解器路线剖析

### 2.1 障碍函数与 Newton 步（linear.py / scattering.py / quotient.py）

障碍 Φ(z)（`quotient.barrier_terms` quotient.py:8-52 + `JointProblem.state` operators.py:71-91）：
−Σ_disks log(2y−x²−t y²)（t 是每行正缩放，operators.py:37-39，等价重排）
−w_χ Σ_groups log(ε²−‖r_g‖²)（w_χ=750）
−[密度 L4：−Σ log(v_i−ρ_i²)−log(1−Σv_i²) 对辅助 v 内层极小化，quotient.py:23-33 的 12 步牛顿]
−Σ_Gram log det G −Σ_FF log(1−‖kF‖²/cap) −Σ_moments log(1−r²)，r=(矩−目标)/误差。
目标：min −c·z/μ + Φ。Hessian 写成特征矩阵 FᵀF、梯度 −Fᵀb+c/μ（每个约束贡献若干行）。方向：正规方程用 PCG（linear.py:102-111），预条件子为列缩放 F 的 QR 上三角 R（每 3 步刷新，linear.py:98-101），fresh 64 步/复用 12 步，相对线性残差须 ≤1e−6；精确可行线搜索（倍增＋60 步二分，linear.py:128-143）；Newton 减量 <1e−14 判"centered"后 μ 减半（linear.py:11, 181, 205）。IR-only 版本 `barrier_support` 另加"径向居中"（整体缩放 z，scattering.py:54-78）。对偶乘子直接由障碍梯度读出（`JointProblem.duals` operators.py:127-134；scattering.py:187-195）。

### 2.2 Phase I、centering、support gap 证书

- Phase I（`PhaseOneProblem` operators.py:143-205；`gauge.initialize` gauge.py:85-141）：min τ，散射/χ/L4 保持硬约束，Gram+τI⪰0、FF 松弛 1+τ−‖kF‖²、矩 1+τ±r。种子是正弦型解析散射种子 c∝(0, v/8, v/8, v⊗v, …) 缩放到严格可行（gauge.py:105-124）。标准 epigraph Phase I。
- "support gap"证书（`imaginary.support_outer` imaginary.py:159-225；`certificates.joint_audit` certificates.py:84-219）：给定乘子 (kR,kI,y,Gram 对偶 Z,矩 m,FF h)，把 Lagrangian 残差 r = c − Σ乘子·约束行 用节点恒等式（Im f00 行、Im f20 行的自由列，`verify_scattering_layout` certificates.py:65-73）精确消去自由坐标，剩余密度残差用 Hölder 对偶 B·‖r‖_{4/3} 上界，加 disk 支撑项 Σ[a²/(√(a²+b²)−b) 或 √+b]、χ 项 ε‖y‖、Gram/FF/矩常数项，全部在 Arb 里对**保存的 float64 H** 精确计算，得 upper；lower = 可行点目标值。upper−lower ≤ `--gap`（默认 1e−4）即 `support_optimality_certified`（__init__.py:236）。数学含义：**该有限问题（float64 H 视为精确）上线性目标的弱对偶间隙**。它不是内点法的对偶残差意义上的最优性证书，也不含角积分/求积误差。
- Arb 局部下降证书（`merit.certified_local_decrease` merit.py:447-495）：在长双精度两个障碍值相减被舍入吞没时，用精确二进制有理数重算方向导数与自协调曲率，给出 Ψ(0)−Ψ(t) ≥ tD − t²h/(2(1−t√h)) 的严格正下界。这是为自研 Newton 的停机准则（减量 1e−14）**制造出来的问题**的修补。

### 2.3 Clarabel / SCS / highspy 用在哪里

- Clarabel：(a) `quotient.joint_hull_candidate` (quotient.py:85-165)：在 ≤~200 个已存父振幅的**凸包**内联合求电流变量（PSD(3)、SOC、非负锥），D 阶段初始化，0.26 s 解出（ATTEMPTS.json D3_hull_009）；(b) `conic.initialize_conic`/`solve_joint_cones` (conic.py:89-168)：全空间锥模型，Clarabel 或 SCS。
- SCS：仅 conic.py:92-107 与 toy 测试 test_mainline.py:147。
- highspy/`linprog`：IR-only 的"normal dual" LP（quotient.py:240-335；__init__.py:77-98），用于改进对偶上界，仍在用。

### 2.4 这是不是标准锥规划？规模能否用现成求解器？

锥类型（`conic.compile_joint_cones` conic.py:12-58 已经把它写成了标准形式）：1500 个 SOC(4)（幺正盘）、100 个 PSD(3)（Gram）、2 个 SOC(5)（或 1 个 SOC(9)，χ 球）、3775 个 SOC(3)+1 个 SOC(3776)（L4 球的精确提升）、14 个 SOC(3)（FF 上界）、8 条线性（矩盒）。变量 4076（+3775 辅助）。A 矩阵 28268×7846、4.7e7 非零（conic_model.json），主要来自稠密的 3011×3876 散射行。这是一个中等规模、**块很小**的 SOCP/SDP，MOSEK/Clarabel 的标准适用范围；作者后续代码就是用 CVX+MOSEK 解同规模的 M=50 问题。

仓库放弃现成求解器的实证依据（`SOLVER_DECISION_ZH.md`、`ATTEMPTS.json`、两次 interlaced_full_conic 报告）：
- 2026-09-08 D3_joint_001/002：Clarabel 全空间，`solver_seconds=300`，MaxTime 时已跑 26/25 次迭代（约 12 s/迭代）——即它**在正常迭代**，只是 300 s 不够一次内点法收敛（通常 50–100 次）；003–006 NumericalError（迭代 1 次，疑为未缩放的 1e−7 量级矩行）。
- 2026-09-11 interlaced_full_conic_01（Clarabel, 600 s）与 _scs_01（SCS, 600 s）：report 显示 `exit_code=-15`（SIGTERM）、`timeout=False`、无 conic_native.json、progress 停在"full_conic_compiled"——**两次都是在 690 s / 838 s 时被人为终止，没有任何求解结果**。
- 没有任何一次"同一冻结模型、同一精度、充足时限"的 MOSEK/Clarabel 对照；`SOLVER_DECISION_ZH.md` 自己也承认"不能宣称自研 Newton 更优"。

### 2.5 自研 barrier-Newton 带来的已知问题（结合仓库自述）

1. 中心化慢：每个 μ 层最多 100 步 Newton（linear.py:88），`NEWTON_REAUDIT_ZH.md` 记录 mid 点在 μ=6.25e−7 处 13 步减量仅从 12.28 降到 11.56；大量 run 以 `solver_seconds` 超时结束（`resolve_start_mu` gauge.py:263-295 专门为续算而设）。
2. 线性求解残差：`B_SOLVER_FIX_ZH.md` 承认旧逻辑在复用 QR 时 CG 实际相对残差达 0.04 仍采用该方向；修复后加了 1e−6 门槛（linear.py:111-116）。
3. 对数值抵消：近中心时 Ψ(0)−Ψ(t)≈1e−12 量级被舍入吞没，需要 merit.py 的 Arb 证书（`LOCAL_DESCENT_PROOF_ZH.md`）。
4. 无对偶证书：Newton 本身不产生最优性证明，靠 support_outer 的弱对偶界；gap 只到 1e−4…1e−6（内点法常规相对 1e−8…1e−10）。
5. 代表点不是极值点：`require_center` 返回的中心依赖 χ 权重 750（linear.py:50、ir.py:145），且 `SCIENCE.md` 承认"不同路径/目标面可能选择不同完整幅度"。
6. 复杂度：linear+scattering+quotient+conic+merit+ir+endpoints+gauge ≈ 1830 行求解/坐标变换代码，对比 kernels+operators(核部分)+model ≈ 900 行物理模型；额外的可逆换元（absorptive、subtracted、radial_coordinates __init__.py:129-157）各自引入失败模式（`stiff` 行、pivot 失败等）。
7. 尽管如此，Phase I 得到的 4076 变量见证（F_M50_L12）和 support 证书是真实产出；Arb 证书对**任何**来源的对偶都适用（hull 路径已把 Clarabel 对偶转成 `fiber_duals.npz`，quotient.py:156-163）。

### 2.6 明确回答

- Newton 路线是论文要求的吗？**否。** 2309 v3 只规定凸约束和"最大化线性泛函"（TeX 931），未指定任何算法；作者公开代码用 CVX/MOSEK。
- 是必要的吗？**否。** 问题是标准 SOCP+PSD(3)，`conic.py` 已能编译成 Clarabel/SCS 形式；被放弃的原因是 5–15 分钟的时限与人为终止，不是求解器不能解。
- 是自造复杂度吗？**是。** 自研路线引出 Phase I、QR/CG 修复、Arb 下降证书、径向居中、多套坐标变换、续算/缓存身份体系等一整层工程；其中真正有科学价值且可独立于 Newton 存在的只有 `support_outer`/`joint_audit` 的 Arb 对偶/原式审计。

---

## 3. 模型偏离清单

| 仓库要素 | 论文对应 | 判定 |
|---|---|---|
| 变量 {T0,σ,ρ1,ρ2} 3876 + {ImF,R} 200 (kernels.py:99-103, model.py:112) | (3.62),(3.69) TeX 926-929, 1010 | 论文规定 |
| 中点 φ 网格 x_j=4/cos²(φ_j/2)，w_j=x'_j/M (operators.py:207-214) | (3.58)-(3.61) TeX 913-923 | 论文规定 |
| 直道 PV+跳跃：Re=K_{kj}(cot 核)，Im=δ_{kj} (kernels.py:105-113, 136-143) | K 只在 FF (3.66)-(3.67) TeX 970-985 明示；散射未写 | 合理离散化（与作者后续代码一致） |
| 交叉道解析 Legendre-Q 角投影 (operators.py:216-242) | (2.9)-(2.10) TeX 561-568，无离散细节 | 合理离散化 |
| 1500 圆盘 = 50 节点 × 30 波，2×2 PSD (model.py:45-54) | (2.12)/(3.71)，M=50, L=10 TeX 1068-1071 | 论文规定 |
| S0/P1 100 个实 3×3 Gram，k² 由 (2.33) (model.py:56-65, 137-145) | (3.70) TeX 998-1006 | 论文规定 |
| 手征 8 残差 (kernels.py:69-77) | (3.63)-(3.64) TeX 938-953 | 论文规定 |
| separate-L2 / combined-L2 球 (model.py:16-19) | "some norm"；字面为 8 条标量 \|·\|≤ε（L∞ 盒） | 仓库自创且可能改变物理：两者都严格小于字面读法；combined 系 output-assisted |
| 四矩 raw 绝对 .002 盒 (quotient.py:54-69, model.py:146-147) | (3.72)-(3.74) TeX 1020-1033；ε_SR=2e−3 TeX 1141；量纲不闭合 | 论文口径未定；仓库读法可能改变物理（相对宽度 84/1.8/47/1.4%） |
| hard-midpoint 掩码 s_i≤s0 (quotient.py:78-79) | (3.72) 求和到 M、积分到 s0 | 合理离散化（含 28–30% 已知端点误差，test_linear.py:27） |
| clipped-phi (quotient.py:80-82) | 无 | 仓库自创但无害（变体、非主线） |
| FF 平方上界 7 节点×2 (model.py:148-154) | (3.75) TeX 1036-1041, ε_FF=6e−5 | 论文规定；m_q 取平均为仓库补充 |
| printed 矩 vs eq250 重算，mq-rule (quotient.py:57-68, __init__.py:73-76) | (2.50)/(2.56) | 论文规定；差 8.5% 为已声明输入不一致 |
| **B=377500，L4 双密度球，B(M)=100·d_M** (operators.py:280-281, quotient.py:14-33) | **无**（来自 2403 代码） | 仓库自创（移植）且可能改变物理 |
| infinity=free（T0 自由，无高能条件）(run.py:240, kernels.py:93-94) | (2.8) 含 T0 作为变量；论文无高能条件 | 论文规定 |
| T0=0 端点切片、FF 二阶零点、5 条渐近不等式 (endpoints.py:243-355) | 无 | 仓库自创（诊断分支，会改变可行集） |
| analytic-cardinal 族 (basis.py, analytic.py, quadrature.py) | 无 | 仓库自创（比较族） |
| 3123/3582 盘 = analytic-cardinal + 额外采样能量 (sampling.py:23-66) | 无（论文只在 M 个节点检查幺正） | 仓库自创（更强模型，非主线） |
| 目标 = 支撑函数 d·(x,y) / 固定 x 截面最大化 y (scattering.py:39-43) | 从内点沿射线最大化 t，TeX 931 | 合理等价改写 |
| +x tip (analysis.py:219-224) | "tip of the shape (red dot)" TeX ~1160 | 论文规定 |
| xref=Weinberg f00(3)@fπ=92 MeV；ref/mid 截面规则 (analysis.py:10, 231-239; PV_REPRESENTATIVE_RULE.json) | "near the chiral point (black dot)"，无规则 | 仓库自创但无害（预注册规则） |
| require-center：返回加权障碍中心，χ 权重 750 (linear.py:50, 202-206; ir.py:145) | 论文点为求解器极值解 | 仓库自创且可能改变物理（所选完整振幅/相移） |
| 相移 δ=unwrap(arg S)/2、η=\|S\|、E=.14√s (analysis.py:304-309) | (2.11) | 论文规定 |
| 90° 线性插值读 ρ 质量 (analysis.py:126-131) | "phase shift crosses π/2" TeX ~1163 | 合理离散化 |
| 五组 (M,L) 分辨率 (io.py:10) | 附录 A | 论文规定 |
| Watson 目标 `--objective watson` (imaginary.py:40-59) | 无（2403 方法） | 仓库自创（显式标记、非主线） |

---

## 4. 代码质量与可信度

### 4.1 测试覆盖（107 项 = 62 函数含参数化）

约 75 项是对小 M（2–4）玩具问题的**独立数学检查**：角核/有理行对 mpmath 直接积分（test_kernels.py:14-30, 56-78，18 项）、Hilbert 核 sin→cos（:31-45）、减除坐标恒等式（:87-149）、障碍梯度/Hessian 有限差分（test_kernels.py:150-167；test_mainline.py:121-161，12 项）、幺正性从 η e^{2iδ} 反推（:168-194）、FESR 单位与 (2.50) 的独立复围道积分（test_linear.py:9-52）、Gram 复矩阵 U†BU 对照（test_mainline.py:55-89）、Farkas/对偶界（:205-215, 269-278, 303-312）、cardinal 重构与 Arb 圆参数角积分（test_kernels.py:209-264）。约 30 项是数据身份/哈希/CLI 门禁/绘图拒绝等实现镜像（test_kernels.py:195-208, 311-350；test_linear.py:104-178, 248-280, 331-339；test_mainline.py:224-244, 313-350）。**没有任何测试在 M=50 生产问题上与外部求解器或论文数值对照**；唯一跨求解器检查是 M=4 玩具上 SCS 原始残差 <1e−7（test_mainline.py:147）。

### 4.2 明显的风险点（非 bug 的敏感开关与已知系统误差）

1. FESR 误差口径（§1.6）：raw .002 在归一化单位下是极不均匀的盒；这是决定 ρ 能标的最大隐藏旋钮。
2. hard-midpoint 端点：n=1 矩权重 28–30% 系统误差被保留为"打印规则"（quotient.py:78-79；test_linear.py:27）。
3. B 正则化无平台（§1.8）：区域大小随 B 单调变化。
4. χ 权重 750 与"中心"代表点：改变返回的完整振幅（linear.py:50；ir.py:145）。
5. 相位分支：`np.unwrap` 在 50 节点粗网格上（analysis.py:309），P1 附近 ΔE≈60 MeV，η 低至 0.24 时 S 接近 0，分支可能误判；文档已标注"winding scope"。
6. Arb 证书把 float64 舍入后的 H 当精确输入（operators.py:294 `float(v.mid())`）；PV 行的求积/角积分误差未进证书（kernels.py:24 `integral_enclosures=False`）。
7. `joint_audit` 会修改被审计点（提升自由高能 R，certificates.py:149-152），证书属于修改后的点。
8. 归一化/符号/单位逐项核对：S=1+iπ√(1−4/s)f ✓（model.py:39-43）；k0²,k1² 与 (2.42) ✓；FESR 雅可比 (π/M)s' ✓（quotient.py:77）；Weinberg (2.18) ✓；R01/R21 ✓；Q_ℓ(−z)=(−1)^{ℓ+1}Q_ℓ(z) ✓；ρ₂ 非对角 /2 打包 ✓。未发现符号或单位错误。

### 4.3 Arb 高精度证书证明的命题

- `support_outer`/`joint_audit`（imaginary.py:212-224；certificates.py:154-155, 210-218）：**仓库有限模型内**——保存的 float64 C/ImF/R 满足 1500 个节点圆盘、χ 球、L4、100 个 Gram 全部主子式、4 个矩盒、14 个 FF 上界；且线性目标在该有限可行集上的上界（弱对偶）。不是连续振幅可行，不是论文问题的解。
- `analytic_scattering_audit`/`analytic_fiber_audit`（certificates.py:244-316）：analytic-cardinal 族在**采样能量**上的区间可行（含行半径），仍是采样、非全能量。
- `native_rho_certificate`/`native_peak_bounds`（certificates.py:9-47）：相邻节点 P1 强度严格离散峰，明确"不证明极点"。
- 固定振幅 Farkas（gauge.py:76；certificates.py:206-208）：负上界 ⇒ 该 C 在有限模型内无电流扩展。
- `merit.certified_local_decrease`：缓存锥模型的局部下降，纯数值命题。

---

## 5. 结论

仓库主线与论文是**同一个凸问题的不同数值表示，外加三处改变可行集或代表点的仓库选择**：(1) 论文没有的双密度 L4 正则化 B=377500（移植自作者后续代码，且仓库自证无平台）；(2) 手征约束取 L2 球而非字面的逐点 L∞ 盒；(3) 代表点取 χ 权重 750 的障碍中心而非极值解。此外 FESR 误差口径（raw .002）虽是字面读法，但把四个矩的相对容差拉到 84%/1.8%/47%/1.4%，是 ρ 能标最敏感的开关。除此之外——变量、网格、PV+跳跃色散、解析角投影、1500 盘、100 个 Gram、FF 上界、+x tip、五组分辨率——都是论文规定或合理离散化；仓库对论文的逐项对照和自我标注是诚实的。

Newton 路线是"没苦硬吃"：论文只要求"最大化线性泛函"，问题是标准 SOCP+PSD(3)，`conic.py` 已把它编译成 Clarabel/SCS 形式；放弃现成求解器的证据只是 300–600 s 时限内未收敛、两次全空间运行被 SIGTERM 终止，从未做过充足时限的对照。自研 barrier-Newton 随后引发了线性残差、对数抵消、慢中心化、无最优性证书、代表点依赖障碍权重等一串已被仓库自己记录的问题，代码量约为物理模型的两倍。

要回到论文原问题，最少的改动是：(a) 用 Clarabel/MOSEK 在充足时限（小时级）和行缩放下直接解 `compile_joint_cones` 的锥模型，返回极值解，保留 `joint_audit` 作 Arb 证书；(b) 增加"8 条标量 |r|≤ε"手征选项并以此为默认；(c) 把 B 作为显式扫描参数或删除，报告 B 依赖；(d) 明确 FESR 误差归一化并与作者代码约定核对；(e) 去掉 require-center/χ 权重对代表点的影响；(f) 把 analytic-cardinal、T0=0、渐近不等式、3123/3582 盘继续留在诊断分支之外。
