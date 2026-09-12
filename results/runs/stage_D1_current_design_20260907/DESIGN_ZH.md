# D1：从同一散射振幅连接两类形状因子和电流

2026-09-07。当前只完成独立算子构造与数学控制；未运行 M50 物理准备或联合优化，D3 尚无联合可行点。用户随后要求先重新审计原方法，因此生产接线与阶段验收暂停。实现为 `operators.current_operators(H,kappa,energies,waves,M,config=None,bits=256)`；独立控制在 [independent_checks.json](independent_checks.json)。

## 论文公式及有限接口

令 p=1+2M+M²+M(M+1)/2，保存变量为

\[
w=[C_{\rm flat}[p],\ v_0[M],\ v_1[M],\ r_0[M],\ r_1[M]],
\]

其中 v=ImF，r 是**电流谱密度**，不是散射的单／双密度，也不是 r/k²。M50 共4076变量；切片依次为[0,3876)、[3876,3926)、[3926,3976)、[3976,4026)、[4026,4076)。原振幅的全部系数保留；其 T0 限制由调用方的散射问题决定。

论文 Eq.(2.35–37) 给 F1(0)=1，F0(0)≈mπ²=1。这里把标量的低能近似作为声明的归一化输入。Eq.(3.65–67) 在原生节点给

\[
F_{a,i}=1+\sum_jK_{ij}v_{a,j}+i v_{a,i},\quad
K_{ij}=\widetilde K_{i+j-2M-1}-\widetilde K_{i-j},\quad
\widetilde K_m=\frac{1-(-1)^m}{2M}\cot\frac{m\pi}{2M},\quad\widetilde K_0=0.
\]

当前 `kernels.pv_matrix` 与这项**明确给出的 FF 核**逐项相同；其名称不表示恢复了作者的整个散射 PV 实现。令 T_jn=sin(nφ_j)、W=T⁻¹，则有限函数 F_a(s)=1+Σ_jn W_nj z(s)^n v_aj，q(0)=0，包含 Nyquist 模式。独立控制核验了 Eq.(3.67) 与 mode1/Nyquist 的 sin→cos 关系。两类 FF 共用 K，不能用未减除 H= q+b 直接替代，否则 F(0)会改变。

仅在原生 M 节点建立当前 Gram：从传入的同一 H 和 (energy,I,ell) 元数据匹配 S0=(0,0)、P1=(1,1)，以 `Re f=H[id]C`、`Im f=H[N+id]C`、`S=1+iκf` 连接。不能沿用2403 MATLAB的不同分波排序，也不能把增加的散射能量误认成新增电流变量。其它 N−2M 行仍用普通散射约束；保留全部 N 行的散射约束也只产生冗余的原生主子块，不改变集合。电流谱在节点之间的插值尚未定义。

Eq.(2.33) 的正因子为

\[
k_0^2=\frac{3\sqrt{1-4/s}}{256\pi^5},\quad
k_1^2=\frac{(s-4)\sqrt{1-4/s}}{384\pi^5},\quad\mathcal F=kF.
\]

原生节点严格大于4，k 由 Arb 的正平方根后导出。`model.current_gram` 使用 diag(1,1,1/k) 对论文矩阵的合同变换，底右为 r/k²；该检查器可复用，不能在阈值 k=0 时除以 k。

## 实三阶 PSD，不必扩成六阶

记 S=a+ib、mathcal F=u+iv。对论文 Eq.(2.40)/(3.70) 使用固定酉矩阵

\[
U=\begin{pmatrix}1/\sqrt2&i/\sqrt2&0\\1/\sqrt2&-i/\sqrt2&0\\0&0&1\end{pmatrix},\quad
U^\dagger G U=\begin{pmatrix}1+a&b&\sqrt2u\\b&1-a&\sqrt2v\\\sqrt2u&\sqrt2v&r\end{pmatrix}.
\]

因此每块实 PSD 的按列上三角 svec 为

\[
[2-\kappa\Im f,\ \sqrt2\kappa\Re f,\ \kappa\Im f,\ 2k(1+Kv_a)_i,\ 2k v_{a,i},\ r_{a,i}].
\]

排列是[G11,√2G12,G22,√2G13,√2G23,G33]，可以直接接现有实 PSD 三角锥。M50 是100个3阶块，并无额外 Gram 决策变量。接口返回 `slack=constant+linear@w`；若求解器采用 A w+s=b，则 A=−linear、b=constant。不要漏掉 FF 的常数1或 svec 的√2。

独立控制取 S=η exp(2iδ)、F=|F|exp(iδ)，其精确谱下界为 r=2k²|F|²/(1+η)。η=1和.6的矩阵与 `model.current_gram` 的独立酉/合同变换一致，减半 r 后出现负特征值。没有强加 η=1、Watson 等式或两 pion 谱饱和；这些仅作为闭式控制。

## 四矩与高能 FF 的连接

直接复用 `model.fesr_data(M,UVConfig)`。四行依次为(S0,0)、(S0,1)、(P1,−1)、(P1,0)，`moment_linear` 只作用于 r0/r1，目标和误差都返回 raw 单位。默认目标约(.0023851042066,.1118381387347,.0042280457143,.1457112069971)，误差逐行.002。Jacobian 是(π/M) ds/dφ，不能用散射 Cauchy 核的无π权重替代。

默认硬节点 cutoff 下 M50有43个低能、7个高能节点。每个高能节点、每个电流增加一个三维 SOC：

\[
[\sqrt{c_a},\ k_a(1+Kv_a)_i,\ k_a v_{a,i}]\in Q_3,
\quad(c_0,c_1)=(2m_q^2\epsilon^{FF},\epsilon^{FF}/2).
\]

默认 c0=1.9544387755e−7、c1=3e−5，共14个 SOC。它们约束 |mathcal F|²；不能对 |F| 使用相同数字。raw/normalized 误差、cutoff 和 mq 分支全部沿用 UVConfig，不默默切换。不额外施加 F(∞)=0 的有限系数等式：论文这里给的是离散高能不等式，连续极限需另有明确论证。

返回布局：`matrices={gram_linear:CSR(600,4076), moment_linear:CSR(4,4076), ff_cap_linear:CSR(42,4076)}`；`arrays` 保存常数、K、k²/k、节点、原散射索引、四矩目标/误差及FF上限；`metadata` 保存变量切片、块顺序、锥大小、UVConfig和单位。以上 M50形状已做合成数据控制，未作为物理准备结果。

## 与本地后续作者代码的边界

v3来源为 `references/2309.12402v3-source/prd_submission_2.tex`：690–729行是 k、归一化和色散关系，728–781行是电流谱与Gram，962–1019行是节点核和两类PSD，1022–1041行是四矩及高能FF。

本地 `theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m` 确认作者用 CVX/MOSEK、仿射FF和Hermitian PSD，而其具体物理模型已经改变：nu0=−20、mpi=139.57MeV、匹配约2GeV、S0/P1/D0三电流、正 B-spline 谱系数、九矩、PSD只施加到s0，及不同FF高能界。`GTB_numerics.nb` 4480–4538行另有从nu0到0的减除修正 KFF=I⊗(K+1 K0)。这些不能无证明移植到本次v3两电流、nu0=0、全M节点设定。其目录README也明确将它标作2403的后续快照。D1当前实现选择的是v3明确节点FF核与已经声明的条件UV合同，而非复制后续三电流程序。
