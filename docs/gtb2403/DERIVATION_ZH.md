# GTB 2403 物理推导：从物理定义到离散核（Task A2）

本文件为 `PLAN_SDPB_2403_ZH.md` §3 所要求的推导文件。每节从物理定义出发，推到本仓库代码（`mma/gtb2403/GTBKernels.wl` 及 `src/smatrix_bootstrap/gtb/`）实现的离散公式，并以两行收尾："离散公式（本仓库实现）"与"与 2309 的关系"。主线 = 作者代码（PLAN §2）；正文与代码不一致处在对应节写明证据，§13 汇总。

## 0. 约定与引用

- **单位与运动学**：$m_\pi=1$（$m_\pi=139.57$ MeV），$s+t+u=4$，$q^2=(s-4)/4$，$t=-\tfrac{s-4}{2}(1-\mu)$，$u=-\tfrac{s-4}{2}(1+\mu)$，$\mu=\cos\theta$。
- **论文**：arXiv:2403.10772 v4 TeX `arxiv_submission_3_v10.tex`。式号一律用 v4 PDF 编号（TeX 标签→式号见 `paper/labels_v4.txt`），必要时附标签名，如 `HHtdef`=(8.3)。
- **作者代码**：notebook `GTB_numerics.nb`（submodule `third_party/gauge-theory-bootstrap` @801684d，`theories/qcd/papers/arxiv-2403.10772/`），按文本化文件 `GTB_numerics_2403.code.wl` 的单元序号写作 `#k`；MATLAB `GTB_numerics.m` 写作 `M3:行号`（去 CR 后）。本文件只写数学，不含作者代码原文。notebook 注释中的式号有几处沿用旧版编号（#30、#38、#153、#158），以 v4 为准。
- **2309**：`references/2309.12402v3-source/prd_submission_2.tex`，v3 式号：`curF`=(2.33)，`Bmt`=(2.41)，`h37`=(3.67)，`positiveB`=(3.70)，`FFasym`=(3.75)；我们的实现 `src/smatrix_bootstrap/sdp/{grid,hilbert,formfactor,constraints,projector}.py`。
- **矩阵目录**：`/playpen1/shiqiu/gtb2403_prep/notebooks/MATRICES.md`（下称 MATRICES），逐个导出矩阵的定义与索引约定。
- **"核算"**：标此字样的数字由本文件写作时独立计算（mpmath 40–60 位或 numpy 双精度），不读作者导出。
- **网格记号**：$\Delta=\pi/M$，$\phi_j=(j-\tfrac12)\Delta$，$x(\phi)=\nu_0+\frac{8-2\nu_0}{1+\cos\phi}$，$x_j=x(\phi_j)$，$x'_j=x'(\phi_j)$，$\omega_j=\pi\sqrt{(x_j-4)/x_j}$；$M=50$，$\nu_0=-20$，$s_0=(2000/139.57)^2=205.341$，$n_0=\#\{x_j\le s_0\}=39$，$X\equiv x_{n_0}=172.0853$；$c_n=\frac{1-(-1)^n}{2M}\cot\frac{n\pi}{2M}$（$c_0=0$）；$R=3.775\times10^5$。
- **变量（MATLAB 顺序）**：$v=(T_0,\sigma_1(x_{1..M}),\sigma_2(x_{1..M}))\in\mathbb R^{101}$；$\rho\in\mathbb R^{3775}$：先 $\rho_1(x_m,y_i)$（下标 $(m-1)M+i$，**首指标 $m$ 为 s 道**），后 $\rho_2(x_a,x_b)$ 的上三角 $a\le b$ 行优先；$\mathrm{Im}F\in\mathbb R^{150}$（S0、P1、D0 各 $M$ 个节点）；B 样条系数 $c\in\mathbb R^{120}$。
- **行顺序**：1500 行 = [$I=0$：$\ell=0,2,\dots,18$ | $I=2$：$\ell=0,\dots,18$ | $I=1$：$\ell=1,\dots,19$] × 50 节点（#89–#91，`M3:68`）。S0 = 行 1–50，D0 = 51–100，S2 = 501–550，D2 = 551–600，P1 = 1001–1050，F1 = 1051–1100。

**问题规模（代码主线）。**

| 对象 | 个数 | 本文节 |
|---|---|---|
| 振幅变量 $v$ / $\rho$ | 101 / 3775 | §2–§3 |
| $\mathrm{Im}F$ / B 样条系数（消去 $c_1$ 与 5 个不可见系数后） | 150 / 120（102） | §5、§7 |
| 幺正锥（全部 ℓ≤19、全部节点） | 1500 | §1、§10 |
| 复 3×3 Gram（S0、P1、D0，$j\le39$） | 117 | §6 |
| FESR L2 球（每流 3 个矩） | 3 | §7–§8 |
| FF 上界（每流 $j=40..50$） | 33 | §8 |
| 手征 L2 球（8 维、12 维）；等式 $\mathfrak F_0$ | 2；1 | §9 |
| ℓ4 正则化 | 1 | §12 |

## 1. 交叉、同位旋、分波与幺正

**同位旋。** $O(N)$（$N=3$）对称下 $\pi^a\pi^b\to\pi^c\pi^d$ 振幅为
$$\mathcal M=\delta^{ab}\delta^{cd}A(s,t,u)+\delta^{ac}\delta^{bd}A(t,s,u)+\delta^{ad}\delta^{bc}A(u,t,s),\qquad A(s,t,u)=A(s,u,t).$$
投影算符 $P^0=\tfrac1N\delta^{ab}\delta^{cd}$，$P^1=\tfrac12(\delta^{ac}\delta^{bd}-\delta^{ad}\delta^{bc})$，$P^2=\tfrac12(\delta^{ac}\delta^{bd}+\delta^{ad}\delta^{bc})-\tfrac1N\delta^{ab}\delta^{cd}$ 满足 $\delta^{ab}\delta^{cd}=NP^0$，$\delta^{ac}\delta^{bd}=P^0+P^1+P^2$，$\delta^{ad}\delta^{bc}=P^0-P^1+P^2$。代入得 (2.2)（`a10`）：
$$T^0=NA(s,t,u)+A(t,s,u)+A(u,t,s),\quad T^1=A(t,s,u)-A(u,t,s),\quad T^2=A(t,s,u)+A(u,t,s).$$
$t\leftrightarrow u$ 即 $\mu\to-\mu$：$T^0,T^2$ 偶、$T^1$ 奇，故 $I=0,2$ 只有偶 $\ell$、$I=1$ 只有奇 $\ell$（Bose 对称）。

**分波与归一。** (2.3)：$f^I_\ell(s)=\tfrac14\int_{-1}^1d\mu\,P_\ell(\mu)T^I$。逆展开为 $T^I(s,\mu)=\sum_\ell 2(2\ell+1)f^I_\ell(s)P_\ell(\mu)$（由 $\int P_\ell P_{\ell'}=\tfrac{2}{2\ell+1}\delta_{\ell\ell'}$）。与标准约定 $t^I_\ell=\tfrac{1}{64\pi}\int P_\ell T_{\rm phys}$、$S=1+2i\sqrt{1-4/s}\,t$ 比较 (2.5) $S=1+i\pi\sqrt{(s-4)/s}\,f$，得 $f=\tfrac{2}{\pi}t$，即本文的 $T=T_{\rm phys}/(8\pi^2)$。检验：Weinberg $A_{\rm phys}=(s-1)/f_\pi^2$ 给出 (2.14) $A=\frac{s-1}{8\pi^2f_\pi^2}$ 与 (2.15) $f^0_{0,\rm tree}=\frac{2s-1}{16\pi^2f_\pi^2}=\frac{2}{\pi}\frac{2s-1}{32\pi f_\pi^2}$。

**幺正。** 记 $h=\pi\sqrt{(s-4)/s}\,f$（(8.36)），$S=1+ih$：
$$|S|^2=(1-\mathrm{Im}\,h)^2+(\mathrm{Re}\,h)^2=1-2\,\mathrm{Im}\,h+|h|^2\ \Rightarrow\ |S|\le1\iff|h|^2\le2\,\mathrm{Im}\,h\iff(\mathrm{Re}\,h)^2+(\mathrm{Im}\,h-1)^2\le1 .$$
即 $h$ 落在以 $i$ 为心、半径 1 的 Argand 圆内，从而 $0\le\mathrm{Im}\,h\le2$，$|\mathrm{Re}\,h|\le1$，$|h|\le2$（§12 用）。SDP 形式由 Schur 补：$\begin{pmatrix}2\mathrm{Im}h&\mathrm{Re}h&\mathrm{Im}h\\ \mathrm{Re}h&1&0\\ \mathrm{Im}h&0&1\end{pmatrix}\succeq0$。

> **离散公式（本仓库实现）**　$\boxed{(\mathrm{Re}\,h_n)^2+(\mathrm{Im}\,h_n)^2\le2\,\mathrm{Im}\,h_n,\ \ \mathrm{Re}\,h_n=a_n\!\cdot v+A_n\!\cdot\rho,\ \ \mathrm{Im}\,h_n=b_n\!\cdot v+B_n\!\cdot\rho,\ \ n=1..1500}$（实际施加的是 §10 的 Λ 重标形式；#91，`M3:63–65`）
>
> **与 2309 的关系**：$T^I$、$f$ 的归一与 $|S|\le1$ 与 2309 (2.5)–(2.11) 完全相同（`projector.py` 的 `rows`）；2309 对非流道用 2×2 块 (3.71)，2403 每同位旋 10 个波（$\ell\le19$）。

## 2. 以 ν₀ 减除的 Mandelstam 表示

从未减除的 (2.1)（`Adef`）出发。令 $k(x,\nu)=\frac{1}{x-\nu}-\frac{1}{x-\nu_0}$，则 $\frac1{x-\nu}=k(x,\nu)+\frac1{x-\nu_0}$。

- 单重项：$\frac1\pi\!\int\!\frac{\sigma(x)}{x-\nu}=H(\nu)+\frac1\pi\!\int\!\frac{\sigma}{x-\nu_0}$，后者为常数，并入 $T_0$。
- 双重项：$\frac{1}{(x-\nu_1)(y-\nu_2)}=k(x,\nu_1)k(y,\nu_2)+\frac{k(x,\nu_1)}{y-\nu_0}+\frac{k(y,\nu_2)}{x-\nu_0}+\frac{1}{(x-\nu_0)(y-\nu_0)}$。第一项给 $\tilde H(\nu_1,\nu_2)$；交叉项 $\frac1\pi\!\int\!dx\,k(x,\nu_1)\tilde\sigma(x)$，$\tilde\sigma(x)=\frac1\pi\!\int\!dy\,\frac{\rho(x,y)}{y-\nu_0}$，是 $H$ 型单重项；末项为常数。

由此 $\sigma_1\to\sigma_1+\frac{2}{\pi}\!\int\!dy\frac{\rho_1(x,y)}{y-\nu_0}$；$\sigma_2\to\sigma_2+\frac1\pi\!\int\!dx\frac{\rho_1(x,\cdot)}{x-\nu_0}+\frac1\pi\!\int\!dy\frac{\rho_2(\cdot,y)}{y-\nu_0}$；$T_0$ 吸收全部常数；**$\rho_{1,2}$ 不变**。$t$ 与 $u$ 两处给 $\sigma_2$ 的平移相同（$(s,t)$、$(s,u)$ 用同一 $\rho_1$；$\rho_2$ 对称），所以 $H_2(t)+H_2(u)$ 的交叉结构保持，得到 (8.2)–(8.3)（`Arewrite`、`HHtdef`）。这是 $(T_0,\sigma,\rho)$ 的单位三角线性换元：$\rho$ 不动，$\sigma,T_0$ 平移 $\rho$ 的线性泛函。$T_0,\sigma$ 是无界自由变量，正则化只作用于 $\rho$，所以离散可行集（以分波计）不变，前提是减除常数与核用同一求积（作者代码的 $\mathcal W_0$ 即中点规则，#39）。连续层面减除只放宽收敛条件（$\sigma/x^2$ 可积即可）。

**为何取 ν₀ 为圆盘中心**：$H(\nu_0)=\tilde H(\nu_0,\cdot)=0$，每块在圆盘上于 $z=0$ 为零，于是割线上的实部恰为离散共轭算子作用于虚部（§4），无需另算常数。

> **离散公式（本仓库实现）**　$\boxed{A=T_0+H_1(s)+H_2(t)+H_2(u)+\tilde H_1(s,t)+\tilde H_1(s,u)+\tilde H_2(t,u),\ \ H,\tilde H\ \text{以}\ k(x,\nu)\ \text{为核 (8.3)}}$，变量 $T_0,\sigma_{1,2}(x_j),\rho_1(x_m,y_i),\rho_2(x_a,x_b)_{a\le b}$。
>
> **与 2309 的关系**：2309 (2.7) 不减除、圆盘中心在 $\nu=0$，我们的 2309 代码给每个 Cauchy 变换加常数行 $G(0)$（`hilbert.cauchy_oncut_real`，中点规则）。$\nu_0=0$ 时两者相差上述单位三角换元。

## 3. 分波投影与组装

**3.1 Legendre Q 闭式。** $x-t=x+\frac{s-4}{2}(1-\mu)=\frac{s-4}{2}(\zeta-\mu)$，$\zeta=1+\frac{2x}{s-4}$。故
$$\tfrac14\!\int_{-1}^1\!d\mu\,\frac{P_\ell(\mu)}{x-t}=\frac{1}{2(s-4)}\int_{-1}^1\!\frac{P_\ell(\mu)}{\zeta-\mu}d\mu=\frac{Q_\ell(\zeta)}{s-4},\qquad Q_\ell(\zeta)=\tfrac12\!\int_{-1}^1\!\frac{P_\ell}{\zeta-\mu}.$$
乘上 $H$ 的 $1/\pi$ 得 (8.6a) $\mathcal A(s,x)=\frac{Q_\ell(\zeta)}{\pi(s-4)}$。$s>4$ 时 $\zeta>1$；$0<s<4$ 时 $\zeta=1-\frac{2x}{4-s}<-1$（因 $x>4>4-s$），积分定义仍成立且为实数，$Q_\ell(\zeta)=(-1)^{\ell+1}Q_\ell(|\zeta|)$。实现用 $Q_\ell(\zeta)=\frac{\sqrt\pi\,\ell!}{\Gamma(\ell+\frac32)(2\zeta)^{\ell+1}}\,{}_2F_1\!\big(\tfrac{\ell+1}2,\tfrac{\ell+2}2;\ell+\tfrac32;\zeta^{-2}\big)$，对 $|\zeta|>1$ 两侧都给出正确宇称。核算：$(\ell,s,x)=(0,10,7),(3,10,7),(2,\tfrac32,5),(3,\tfrac12,9)$ 处直接数值积分、积分定义的 $Q$ 与 mpmath 第三类 `legenq(·,type=3)`（与作者 #45 所用 `LegendreQ[l,0,3,·]` 同一类型）三者一致到 10 位；$Q_\ell(-3.2)=(-1)^{\ell+1}Q_\ell(3.2)$，$\ell=0..3$ 精确成立。

**3.2 u 道。** $u(\mu)=t(-\mu)$，$P_\ell(-\mu)=(-1)^\ell P_\ell(\mu)$，故 $H(u)$ 的投影为 $(-1)^\ell\mathcal A$。

**3.3 $\delta_{\ell0}$ 常数项。** $\tfrac14\int P_\ell\,d\mu=\tfrac12\delta_{\ell0}$。于是 $H(s)$（与 $\mu$ 无关）只投到 $\ell=0$：$f^{(s)}=\tfrac12\delta_{\ell0}\!\int\!\mathcal K\sigma$；$H(t),H(u)$ 中 $-\frac1{x-\nu_0}$ 部分给 $-\tfrac12\delta_{\ell0}\!\int\!\mathcal W_0\sigma$；$T_0$ 每份给 $\tfrac12\delta_{\ell0}$。

**3.4 (t,u) 部分分式与 W 核。** $t+u=4-s\Rightarrow(x-t)+(y-u)=x+y+s-4>0$（$s>4$ 显然；$0<s<4$ 时 $x+y\ge8>4-s$），故
$$\frac{1}{(x-t)(y-u)}=\frac{1}{x+y+s-4}\Big[\frac1{x-t}+\frac1{y-u}\Big].$$
展开 $\tilde H(t,u)$ 的被积式为四项并投影：
$$\tfrac14\!\int\!P_\ell\,k(x,t)k(y,u)\cdot\tfrac{1}{\pi^2}=\frac{\mathcal A(s,x)+(-1)^\ell\mathcal A(s,y)}{\pi(x+y+s-4)}-\frac{\mathcal A(s,x)}{\pi(y-\nu_0)}-\frac{(-1)^\ell\mathcal A(s,y)}{\pi(x-\nu_0)}+\tfrac12\delta_{\ell0}\mathcal W_{00}(x,y)$$
$$=\mathcal A(s,x)\,\mathcal W(s,x,y)+(-1)^\ell\mathcal A(s,y)\,\mathcal W(s,y,x)+\tfrac12\delta_{\ell0}\mathcal W_{00}(x,y),\quad \mathcal W(s,y,x)=\tfrac1\pi\Big[\tfrac{1}{x+y+s-4}-\tfrac1{x-\nu_0}\Big].$$
第一项换名 $x\leftrightarrow y$ 即 (8.5c)。(8.5c) 与 (8.9) 把 $\mathcal W_{00}$ 写成 "$\mathcal W_0(x,y)$"，只是记号问题。$(s,t)$、$(s,u)$ 为乘积：$\mathcal K(s,x)[\pm\mathcal A(s,y)-\tfrac12\delta_{\ell0}\mathcal W_0(y)]$（(8.5a,b)）。

**3.5 物理区离散（行 $(\ell,j)$，$s=x_j$，整体乘 $\omega_j$）。** 记 $\alpha_\ell(s;i)=\Delta x'_i\frac{Q_\ell(1+2x_i/(s-4))}{\pi(s-4)}$，$w^0_i=\frac{\Delta x'_i}{\pi(x_i-\nu_0)}$，$\hat W_s(p\to q)=\frac{\Delta}{\pi}x'_q\big[\frac1{x_p+x_q+s-4}-\frac1{x_q-\nu_0}\big]$。割线上 $\mathcal K$ 的实部取 $K_{jm}$、虚部取 $\delta_{jm}$（作用于节点值，§4）；其余积分用 φ 中点规则（雅可比 $x'$ 并入核）。

| 片 | 作用于 | Re 系数 | Im 系数 | 单元 |
|---|---|---|---|---|
| $T_0$ | $T_0$ | $\tfrac12\delta_{\ell0}$ | 0 | #84 |
| (s) | $\sigma(x_m)$ | $\tfrac12\delta_{\ell0}K_{jm}$ | $\tfrac12\delta_{\ell0}\delta_{jm}$ | #65 |
| (t) / (u) | $\sigma(x_m)$ | $(\pm1)^\ell\alpha_\ell(x_j;m)-\tfrac12\delta_{\ell0}w^0_m$ | 0 | #45,#65 |
| (s,t) / (s,u) | $\rho(x_m,x_i)$ | $K_{jm}[(\pm1)^\ell\alpha_\ell(x_j;i)-\tfrac12\delta_{\ell0}w^0_i]$ | $\delta_{jm}[(\pm1)^\ell\alpha_\ell(x_j;i)-\tfrac12\delta_{\ell0}w^0_i]$ | #46,#67,#68 |
| (t,u) | $\rho(x_a,x_b)$ | $\alpha_\ell(x_j;a)\hat W_{x_j}(a\to b)+(-1)^\ell\alpha_\ell(x_j;b)\hat W_{x_j}(b\to a)+\tfrac12\delta_{\ell0}w^0_aw^0_b$ | 0 | #41,#43,#69 |

（"±" 中 + 对应 t/(s,t)，− 即 $(-1)^\ell$ 对应 u/(s,u)。）非物理点 $s_*$：去掉 $\omega$，$K_{jm}\to k^*_m=\frac{\Delta}{\pi}x'_m\big[\frac1{x_m-s_*}-\frac1{x_m-\nu_0}\big]$，无虚部，$\alpha_\ell(s_*;i)$ 用 $\zeta<-1$（#50–#61、#95–#98）。

**3.6 交叉组装。** 记各片为 $P^{(s)},P^{(t)},P^{(u)}$（作用于 σ）、$P^{(st)},P^{(su)},P^{(tu)}$（作用于 ρ，首指标对应括号内第一个变量），$\mathcal T$ 为 $\rho(x,y)\to\rho(y,x)$。三份振幅 $[A(s,t,u)],[A(t,s,u)],[A(u,t,s)]$ 中各变量的系数为：

| 变量 | $[A(s,t,u)]$ | $[A(t,s,u)]$ | $[A(u,t,s)]$ |
|---|---|---|---|
| $\sigma_1$ | $P^{(s)}$ | $P^{(t)}$ | $P^{(u)}$ |
| $\sigma_2$ | $P^{(t)}+P^{(u)}$ | $P^{(s)}+P^{(u)}$ | $P^{(s)}+P^{(t)}$ |
| $\rho_1$ | $P^{(st)}+P^{(su)}$ | $P^{(st)}\mathcal T+P^{(tu)}$ | $(P^{(tu)}+P^{(su)})\mathcal T$ |
| $\rho_2$ | $P^{(tu)}$ | $P^{(su)}$ | $P^{(st)}\mathcal T$ |

例：$[A(t,s,u)]$ 含 $\tilde H_1(t,s)=\frac1{\pi^2}\!\iint\!\rho_1(x,y)k(x,t)k(y,s)$，换名 $x\leftrightarrow y$ 后为 $P^{(st)}$ 作用于 $\mathcal T\rho_1$。这正是 (8.10)–(8.11)（#71、#73）。

**3.7 同位旋行。** $I=0$：$N[A(s,t,u)]+[A(t,s,u)]+[A(u,t,s)]$；$I=2$：$[A(t,s,u)]+[A(u,t,s)]$；$I=1$：$[A(t,s,u)]-[A(u,t,s)]$。$T_0$ 行因此为 $\tfrac12\delta_{\ell0}\omega_j\times(N+2,\ 2,\ 0)$，按代码行序 $[I{=}0|I{=}2|I{=}1]$（#75、#76、#82、#85）。
**论文笔误 (8.12)–(8.13)**：两式左侧标签为 $[\tilde\Phi^0;\tilde\Phi^1;\tilde\Phi^2]$，但第二行 $\tilde\Phi^{(t,s)}+\tilde\Phi^{(u,t)}$ 是 $T^2=A(t)+A(u)$，第三行 $\tilde\Phi^{(t,s)}-\tilde\Phi^{(u,t)}$ 是 $T^1$，即实际行序为 $[I{=}0;I{=}2;I{=}1]$。(8.14) 写 $\tfrac{\delta_{\ell0}}{2}[N+2;0;2]$，按 $[0;1;2]$ 标签读是**对的**（$T^1$ 不含 $T_0$），所以 (8.12)–(8.14) 三式行序互不一致。代码全程用 $[0;2;1]$（#85 给 $[5,2,0]$；#89 的 `listindex2` 取第二块偶 ℓ）；#128 的注释 "I=0,1,2" 同样误标。

**3.8 ρ₂ 对称化。** $f=\sum_{a,b}C_{ab}\rho_2(x_a,x_b)$，$\rho_2$ 对称，独立变量取 $a\le b$：系数 $C_{ab}+C_{ba}$（$a<b$）、$C_{aa}$（对角）。代码以 $2\cdot\tfrac12(I+\mathcal T)$ 作用后把对角列减半，再取上三角行优先（#78）。变量数 $1+2M+M^2+\tfrac{M(M+1)}2=101+3775$。

**3.9 单项锚点（核算，mpmath 40 位，供 A3/A5 的单项测试）。**

| 量 | 值 |
|---|---|
| $x_1$，$x'_1$，$\omega_1$ | 4.0059227368677782684，0.37711516955362880399，0.12079801770548657922 |
| $x_{39}$，$x_{50}$ | 172.08533939065126099，97256.336691443898924 |
| $K_{11}=K_{12}=-K_{21}=-K_{50,50}=c_1=\frac1{50}\cot\frac{\pi}{100}$ | 0.63641031907547916079 |
| $w^0_1$；$K_{0,1}$；$K_{0,50}$ | 3.1418510647329832647e−4；−1.5686029267701699871e−3；−2.6181015377878926974e−4 |
| $\alpha_0(x_1;1)$，$\alpha_1(x_1;1)$ | 9.40698778619304754e−4，2.3163187901303562331e−7 |
| $\alpha_0(\tfrac12;1)$，$\alpha_1(\tfrac12;1)$（$\zeta<-1$，符号检验宇称） | 2.2294459075099235473e−3，−7.1903201018980934001e−4 |
| $\Lambda_1(x_1)$，$\Lambda_2(x_{39})$，$\kappa_0(x_1)$ | 0.019232713130788359063，0.7354170156847934354，1.4724544226690316952e−6 |

例：行 $(I{=}2,\ell{=}0,j{=}1)$ 中 $T_0$ 的实部系数为 $\omega_1$（§3.7）；行 $(I{=}0,\ell{=}0,j{=}1)$ 中 $\sigma_{1,1}$ 的虚部系数为 $\tfrac32\omega_1$。

> **离散公式（本仓库实现）**　$\boxed{\mathrm{Re}/\mathrm{Im}\,h^I_{\ell}(x_j)=\omega_j\textstyle\sum(\text{3.5 表各片})\ \text{按 3.6 表交叉、按 3.7 组同位旋},\ \rho_2\ \text{系数}\ C_{ab}+C_{ba}\,(a<b),\,C_{aa}}$；行选 $I=0,2$ 偶 ℓ、$I=1$ 奇 ℓ，共 1500 行；$s_*$ 行同式去 $\omega$、无虚部。
>
> **与 2309 的关系**：闭式与 2309 代码相同（`projector.py`：$Q_c=\frac2aQ_\ell$，$a=\frac{s-4}2$；部分分式；`fold_sym` 同为 $C+C^{\mathsf T}$ 保留对角、`triu` 打包）。2403 多出 $-\tfrac12\delta_{\ell0}\mathcal W_0$、$+\tfrac12\delta_{\ell0}\mathcal W_{00}$ 减除项，以及 $\omega$、Λ 因子。

## 4. 圆盘映射、中点节点与离散共轭函数

**4.1 映射。** (8.15) $z(\nu)=\frac{a-\sqrt{4-\nu}}{a+\sqrt{4-\nu}}$，$a=\sqrt{4-\nu_0}$，$z(\nu_0)=0$。上岸 $\nu=x+i0$（$x>4$）：$\sqrt{4-\nu}=-i\sqrt{x-4}$，$z=\frac{a+ib}{a-ib}$，$b=\sqrt{x-4}$，$|z|=1$，$\phi=2\arctan(b/a)\in(0,\pi)$。反解 $x=4+a^2\tan^2\frac\phi2=\nu_0+\frac{8-2\nu_0}{1+\cos\phi}$（(8.16)），$x'=\frac{2(4-\nu_0)\sin\phi}{(1+\cos\phi)^2}$（(8.18)）。

**4.2 核 (8.20)。** 由 $\sqrt{4-\nu}=a\frac{1-z}{1+z}$ 得 $\nu-\nu_0=\frac{4a^2z}{(1+z)^2}$，又 $x-\nu_0=\frac{a^2}{\cos^2(\phi/2)}$，$x-\nu=a^2[\tau^2+w^2]$（$\tau=\tan\frac\phi2$，$w=\frac{1-z}{1+z}$），$x'=\frac{a^2\tau}{\cos^2(\phi/2)}$。于是
$$\tfrac1\pi\Big[\tfrac1{x-\nu}-\tfrac1{x-\nu_0}\Big]x'=\frac{4z\tau}{\pi[(1+z)^2\tau^2+(1-z)^2]}=\frac{2z\sin\phi}{\pi(1+z^2-2z\cos\phi)}.$$
核算：$\nu=3,-7.3,10+3i$ 处两边差 $<10^{-40}$。$z\to e^{i\psi}$ 时 $1+z^2-2z\cos\phi=2z(\cos\psi-\cos\phi)$，核变为 $\frac1\pi\frac{\sin\phi}{\cos\psi-\cos\phi}$（主值，(8.22)），另加 $\phi=\psi$ 处的 δ 函数（虚部）。

**4.3 函数类。** 每块振幅（§2）与 §5 的 $G$ 都是 $G(z)=\sum_{k\ge1}b_kz^k$，$b_k$ 实（实解析且 $G(0)=0$）。边界上 $\mathrm{Im}\,G=\sum b_k\sin k\phi$，$\mathrm{Re}\,G=\sum b_k\cos k\phi$。

**4.4 离散共轭算子（推导）。** 已知 $u_j=\mathrm{Im}\,G(\phi_j)$，用 $u(\theta)=\sum_{k=1}^{M}b_k\sin k\theta$ 插值（奇延拓到 $2M$ 个交错点，镜像 $\phi_{1-j}=-\phi_j$）。DST-II 正交性 $\sum_{j=1}^M\sin k\phi_j\sin k'\phi_j=\frac M2\delta_{kk'}$（$1\le k,k'<M$）、$=M$（$k=k'=M$）给出 $b_k=\frac2M\sum_ju_j\sin k\phi_j$（$k<M$），$b_M=\frac1M\sum_ju_j\sin M\phi_j$，插值唯一。在节点上 $\cos M\phi_i=\cos(i-\frac12)\pi=0$，第 $M$ 项消失：
$$\mathrm{Re}\,G(\phi_i)=\sum_ju_j\,\tfrac1M\!\sum_{k=1}^{M-1}\!\big[\sin k(\phi_j+\phi_i)+\sin k(\phi_j-\phi_i)\big].$$
用 $\sum_{k=1}^{M-1}\sin\frac{kn\pi}{M}=\frac{\sin\frac{n\pi}2\sin(\frac{n\pi}2-\frac{n\pi}{2M})}{\sin\frac{n\pi}{2M}}=\tfrac12(1-(-1)^n)\cot\frac{n\pi}{2M}$（$\sin\frac{n\pi}2\cos\frac{n\pi}2=0$），即 $\frac1M\sum_k\sin\frac{kn\pi}{M}=c_n$，且 $\phi_j\pm\phi_i$ 分别为 $(i+j-1)\Delta$、$(j-i)\Delta$，$c$ 为奇函数：
$$\mathrm{Re}\,G(\phi_i)=\sum_{j=1}^M K_{ij}\,\mathrm{Im}\,G(\phi_j),\qquad K_{ij}=c_{i+j-1}-c_{i-j}.$$
代码 #37 以 $c$ 的 $2M$ 周期性实现同一式。

**4.5 精确性与混叠。** 对 $\mathrm{Im}\,G=\sin k\phi$，$1\le k\le M$，插值精确，故 $K$ 精确给出 $\cos k\phi_i$。对 $M<k<2M$，节点上 $\sin k\phi_j=\sin(2M-k)\phi_j$（因 $2M\phi_j=(2j-1)\pi$），$K$ 返回 $\cos(2M-k)\phi_i=-\cos k\phi_i$，误差 $2|\cos k\phi_i|=O(1)$。核算（60 位，$M=50$）：$k=1..50$ 最大误差 $\le1.4\times10^{-59}$；$k=51$ 为 $1.999$，$k=52$ 为 $2.0$。$\max|K_{ij}|=0.6364$。

**4.6 论文记号问题**（不影响代码）：(8.21) 积分上限 "∞" 应为 π；(8.34) 求和 $\sum_{j_1=-M}^{M}$ 若按 $\phi_{-j}=-\phi_j$ 编号，(8.35) 的 cot 宗量会差 1，正确写法是 $j_1=-M+1..M$、$\phi_{1-j}=-\phi_j$；(8.33) 的 "$+\frac{i}{2\pi}\delta$" 与求和归一不符，操作含义是"实部 = $K$·虚部，虚部 = 节点值本身"。

**4.7 两种离散的混用。** 割线上用有限正弦共轭算子 $K$；割线外（$s_*$、§5 的 $K_0$、t/u 道积分）用 φ 中点求积。两者在内点不相同。核算：同一数据向量下中点行与"有限正弦插值在 $z(s_*)$ 的精确值"的最大差为 $1.8\times10^{-19}$（$s_*=\tfrac12$）、$5.1\times10^{-18}$（1）、$1.8\times10^{-16}$（$\tfrac32$）、$8.6\times10^{-15}$（2）、$9.9\times10^{-11}$（$s_*=3$，相对 $2.1\times10^{-9}$，约为 $|z|^{M}=0.661^{50}$）、$7.3\times10^{-21}$（$s=0$）。物理上可以忽略，但 G1 要求 (a)=(c) 到 $10^{-60}$，所以我们的核必须用作者的中点式，不能"改用更一致的"插值式。

> **离散公式（本仓库实现）**　$\boxed{K_{ij}=c_{i+j-1}-c_{i-j},\ \ c_n=\tfrac{1-(-1)^n}{2M}\cot\tfrac{n\pi}{2M},\ \ x_j=\nu_0+\tfrac{8-2\nu_0}{1+\cos\phi_j},\ \ x'_j=\tfrac{(8-2\nu_0)\sin\phi_j}{(1+\cos\phi_j)^2}}$
>
> **与 2309 的关系**：与 2309 (3.67) $K_{ij}=\tilde K_{i+j-2M-1}-\tilde K_{i-j}$ 恒等（$c$ 以 $2M$ 为周期；核算差为 0）。$K$ 与 $\nu_0$ 无关，只有节点 $x_j$ 依赖 $\nu_0$；$\nu_0=0$ 时 $x_j=8/(1+\cos\phi_j)$ 即 2309 (3.60)。2309 未减除，割线上的实部为 $K+\mathbf 1c(0)^{\mathsf T}$（`hilbert.py`）。

## 5. F(0)=1 的形状因子与 K₀

(2.7)：$F(s)=1+\frac1\pi\int_4^\infty dx\big[\frac1{x-s}-\frac1x\big]\mathrm{Im}F(x)$。拆 $\frac1{x-s}-\frac1x=k(x,s)+\big[\frac1{x-\nu_0}-\frac1x\big]$：
$$F(s)=F(\nu_0)+G(s),\quad G(s)=\tfrac1\pi\!\int\!k(x,s)\,\mathrm{Im}F,\quad F(\nu_0)=1+\tfrac1\pi\!\int\!\frac{\nu_0}{x(x-\nu_0)}\mathrm{Im}F\,dx .$$
$G(\nu_0)=0$，$\mathrm{Im}\,G(x+i0)=\mathrm{Im}F(x)$，故 $\mathrm{Re}\,G(x_j)=\sum_mK_{jm}\mathrm{Im}F_m$（§4）。常数项换到 φ：$\frac{x'}{x(x-\nu_0)}=\frac{(8-2\nu_0)\sin\phi}{(1+\cos\phi)^2}\cdot\frac{(1+\cos\phi)^2}{(8-\nu_0+\nu_0\cos\phi)(8-2\nu_0)}=\frac{\sin\phi}{8-\nu_0+\nu_0\cos\phi}$，即 (8.25) $K_0(\phi)=\frac{\nu_0\sin\phi}{\pi(8-\nu_0+\nu_0\cos\phi)}$。中点离散 $K_{0,m}=\Delta K_0(\phi_m)=\frac{\Delta}{\pi}\frac{x'_m\nu_0}{x_m(x_m-\nu_0)}$（#136）。**论文笔误 (8.41)**：写 $K_{0,j_1}=K_0(\phi_{j_1})$，漏了 $\Delta_\phi$（(8.24) 同样漏 $d\phi_1$）；代码有 $\Delta$。核算：$K_{0,m}\in[-2.04\times10^{-2},-2.6\times10^{-4}]$（$\nu_0<0$，全为负）；与插值精确值 $-(z_0^k)^{\mathsf T}S^{-1}$（$z_0=z(0)=0.42020$）差 $7.3\times10^{-21}$。三个流 S0、P1、D0 共用同一核（#138），作用于全部 50 个节点：$j\le39$ 进 Gram（§6），$j\ge40$ 进渐近上界（§8）。D0 的 $F^0_2(0)=1$ 来自 $\langle\pi|\!\int\!T^{00}|\pi\rangle=E\langle\pi|\pi\rangle$（(6.7) 后）；(2.23c) 写成 "$F^0_2(s)=1$" 是笔误，应为 $F^0_2(0)=1$。

> **离散公式（本仓库实现）**　$\boxed{\mathrm{Re}F_{c,j}=1+\textstyle\sum_{m=1}^{M}\big(K_{jm}+K_{0,m}\big)\mathrm{Im}F_{c,m},\ \ K_{0,m}=\tfrac{\Delta}{\pi}\tfrac{x'_m\nu_0}{x_m(x_m-\nu_0)},\ \ c\in\{S0,P1,D0\},\ j=1..M}$（`M3:78`）
>
> **与 2309 的关系**：2309 (3.66) $\mathrm{Re}F=1+K\,\mathrm{Im}F$（`formfactor.real_part_operator`）。$K_0\propto\nu_0$，2309 圆盘中心即减除点 0，故 $K_0=0$。D0 是 2403 新增的流。

## 6. 重标 Gram 矩阵：合同变换与饱和

**6.1 出发点。** 2309 (2.41) / 2403 (2.8)、(8.58)：$G_{ij}=\langle\psi_i|\psi_j\rangle$，$\psi_1=|\text{out}\rangle$，$\psi_2=|\text{in}\rangle$，$\psi_3=\mathcal O|0\rangle$：
$G=\begin{pmatrix}1&S&\mathcal F\\ S^*&1&\mathcal F^*\\ \mathcal F^*&\mathcal F&\rho\end{pmatrix}\succeq0$，其中 $\mathcal F=\kappa_\ell^{1/2}F$，$\kappa_\ell^{1/2}>0$（§7）。

**6.2 合同变换 (8.60)。** $\tilde\psi=T\psi$，$T=\begin{pmatrix}1&0&0\\ i\Lambda^{-1}&-i\Lambda^{-1}&0\\0&0&\kappa^{-1/2}\end{pmatrix}$，$\det T=-i\Lambda^{-1}\kappa^{-1/2}\ne0$（$s>4$ 时 $\Lambda,\kappa>0$），所以 $\tilde G=\bar TGT^{\mathsf T}\succeq0\iff G\succeq0$。逐项（bra 取复共轭系数）：
- $\tilde G_{12}=\Lambda^{-1}(-iS+i)=\Lambda^{-1}(-i)(ih)=h/\Lambda=\tilde h$；
- $\tilde G_{13}=\kappa^{-1/2}\mathcal F=F$；
- $\tilde G_{22}=\Lambda^{-2}\big[\langle\psi_2|\psi_2\rangle-\langle\psi_2|\psi_1\rangle-\langle\psi_1|\psi_2\rangle+\langle\psi_1|\psi_1\rangle\big]=\Lambda^{-2}(2-2\mathrm{Re}S)=2\,\mathrm{Im}\,h/\Lambda^2$；
- $\tilde G_{23}=\Lambda^{-1}\kappa^{-1/2}\,i(\mathcal F^*-\mathcal F)=2\,\mathrm{Im}\mathcal F/(\Lambda\kappa^{1/2})=2\,\mathrm{Im}F/\Lambda$；
- $\tilde G_{33}=\rho/\kappa=\hat\rho$。

得 (8.61)：$\tilde G=\begin{pmatrix}1&\tilde h&F\\ \tilde h^*&2\,\mathrm{Im}\hat h&2\,\mathrm{Im}\hat F\\ F^*&2\,\mathrm{Im}\hat F&\hat\rho\end{pmatrix}$，$\mathrm{Im}\hat h=\mathrm{Im}\,h/\Lambda^2$，$\mathrm{Im}\hat F=\mathrm{Im}F/\Lambda$。**论文笔误 (8.63)**：写 "$\mathrm{Im}\hat h=\Lambda^{-2}h$"，漏了 Im；代码为 $\hat B=B/\Lambda^2$（#129）。

求解器只接受实对称块，用 $\tilde G\succeq0\iff\begin{pmatrix}\mathrm{Re}\tilde G&-\mathrm{Im}\tilde G\\ \mathrm{Im}\tilde G&\mathrm{Re}\tilde G\end{pmatrix}\succeq0$（实 6×6 嵌入，特征值成对重复）；$\tilde h$ 与 $F$ 的相位都是变量，不能用对角酉变换把块化为实 3×3。

**6.3 推论（2×2 主子式）。** $\{1,2\}$：$|\tilde h|^2\le2\mathrm{Im}\hat h$，即幺正（Gram 块蕴含 S0/P1/D0 在 $j\le39$ 的幺正，代码两者都留）；$\{1,3\}$：$|F|^2\le\hat\rho$（(8.44)）；$\{2,3\}$：$2(\mathrm{Im}F)^2\le\mathrm{Im}\,h\,\hat\rho$。

**6.4 饱和与 Watson。** 秩一：$\tilde G=vv^\dagger$，$v_1=1$。$\tilde G_{12}=v_2^*=\tilde h$，$\tilde G_{13}=v_3^*=F$，$\tilde G_{23}=v_2v_3^*=\tilde h^*F=2\mathrm{Im}\hat F$（实），$\tilde G_{22}=|\tilde h|^2=2\mathrm{Im}\hat h$，$\tilde G_{33}=|F|^2=\hat\rho$。故
$$\tilde h=\frac{2\,\mathrm{Im}\hat F}{F^*}=\frac{2\,\mathrm{Im}\hat F\,F}{|F|^2}\ \Longleftrightarrow\ h=\frac{2\,\mathrm{Im}F\,F}{|F|^2}.$$
令 $F=|F|e^{i\alpha}$：$h=2\sin\alpha\,e^{i\alpha}$，$S=1+2i\sin\alpha\,e^{i\alpha}=e^{2i\alpha}$，正是 (2.25) 的 $S=\mathcal F/\mathcal F^*$，$\delta=\alpha$。**论文笔误 (8.76)/(8.78)**：(8.76) 写 $\tilde h=\mathrm{Im}\hat F/F^*$，(8.78) 也无因子 2。按该式 $h=\sin\alpha\,e^{i\alpha}$，$S=\cos\alpha\,e^{i\alpha}$，即 $\eta=|\cos\alpha|$、$\delta=\alpha/2$，两项都不是 Watson（核算 $\alpha=0.7$：$|S|=0.7648$，$\delta=0.35$）。代码 `M3:119–120` 带因子 2，主线按代码。

**6.5 Watson 泛函 (8.77)。** 每行 $\max\{\mathrm{Re}(t^*\tilde h)-\mathrm{Im}\hat h\}$ 在锥 $|\tilde h|^2\le2\mathrm{Im}\hat h$ 上于 $\tilde h=t$、$\mathrm{Im}\hat h=|t|^2/2$ 取到，值 $|t|^2/2$：极大点在饱和边界上，相位等于目标相位。代码（`M3:177–180`）：S0、P1、D0 用 $t=2\mathrm{Im}\hat F\,F/|F|^2$（旧解），S2、D2、F1 用旧解的 $\tilde h$，只取 $j\le n_0$ 的节点，只做一步；该步不再固定 $\mathfrak F_0$，其余约束不变（`M3:122–176`）。

> **离散公式（本仓库实现）**　$\boxed{\tilde G_{c,j}=\begin{pmatrix}1&\tilde h_n&F_{c,j}\\ \tilde h_n^*&2\mathrm{Im}\hat h_n&2\mathrm{Im}F_{c,j}/\Lambda_{\ell_c}(x_j)\\ F_{c,j}^*&2\mathrm{Im}F_{c,j}/\Lambda_{\ell_c}(x_j)&\hat\rho_{c,j}\end{pmatrix}\succeq0,\ j\le39}$，$n(S0,j)=j$，$n(D0,j)=50+j$，$n(P1,j)=1000+j$，$\ell_c=0,1,2$（#131，`M3:83–86`）；以实 6×6 嵌入施加，共 117 块。Watson 目标 $t_n=2\mathrm{Im}\hat F\,F/|F|^2$。
>
> **与 2309 的关系**：2309 (3.70) 直接施加未重标的 $G$（S0、P1）。我们的 2309 代码只做了 $\psi_3$ 的合同 $\mathrm{diag}(1,1,1/k)$（`formfactor.gram_scale`），得 $[[1,S,F],[S^*,1,F^*],[F^*,F,\rho/k^2]]$。2403 再做 $\psi_2\to\Lambda^{-1}(-i\psi_2+i\psi_1)$，加入 D0，且只取 $s_j<s_0$ 的节点。

## 7. κ_ℓ、B 样条与矩的归一

**7.1 κ_ℓ = |𝓕/F|²。** 由 (2.21e)、(2.22e)、(2.23e) 平方：
- $\kappa_0=\frac{6\pi}{256\pi^6}s^{-1/2}\big(\frac{s-4}4\big)^{1/2}=\frac{3}{256\pi^5}\sqrt{\frac{s-4}s}=\frac1{(2\pi)^4}\frac3{16\pi}\big(1-\frac4s\big)^{1/2}$；
- $\kappa_1=\frac{4\pi/3}{64\pi^6}s^{-1/2}\big(\frac{s-4}4\big)^{3/2}=\frac{(s-4)^{3/2}}{384\pi^5\sqrt s}=\frac1{(2\pi)^4}\frac{s}{24\pi}\big(1-\frac4s\big)^{3/2}$；
- $\kappa_2=\frac{\pi/5}{16\pi^6}s^{-1/2}\big(\frac{s-4}4\big)^{5/2}=\frac{(s-4)^{5/2}}{2560\pi^5\sqrt s}=\frac1{(2\pi)^4}\frac{s^2}{160\pi}\big(1-\frac4s\big)^{5/2}$。

与 (8.43)、#150 一致。一般地 $\mathcal F_\ell\propto q^{\ell+1/2}s^{-1/4}F$（$|p|^\ell$ 来自 $(-ip\cdot\Delta)^\ell$，$(q/\sqrt s)^{1/2}$ 来自两体相空间），故 $\kappa_\ell\propto q^{2\ell+1}/\sqrt s$。$F\approx1$ 时 $\rho\approx\kappa_\ell$，即自由 π 的低能谱 (2.20)。

**7.2 ρ̂ 与 ρ̂(0)=1。** $\rho=\kappa\hat\rho$（(8.42)）。$4\le s\le16$ 只有两 π 态，$\rho=|\mathcal F|^2\Rightarrow\hat\rho=|F|^2$。在 $0<s<4$ 上 $F$ 实，作者把 $\hat\rho$ 延拓为 $F^2$，由 $F(0)=1$ 得 $\hat\rho(0)=1$。这是延拓假设，不是约束推出的结论。

**7.3 B 样条。** 夹紧节点 $\{0^4,\ kX/37\ (k=1..36),\ X^4\}$，共 44 个节点、40 个三次基（#141–#145）；非负、单位分解，$B_1(0)=1$，其余在 0 处为零，故 $\hat\rho(0)=1\iff c_1=1$（(8.47)）。$c\ge0$ 是 $\hat\rho\ge0$ 的充分条件，不是必要条件，所以会缩小可行集。核算发现两点：
- $B_1(x_1)=2.67\times10^{-3}$：$c_1=1$ 只固定 $s=0$ 处的值，不固定阈值处；阈值幂律 $\rho\propto q^{2\ell+1}$ 由 $\kappa_\ell$ 给出。低端节点很密（$x_1..x_8$ 都在第一个节点区间 $[0,4.65]$ 内），在那里 $\hat\rho$ 被限制为一个三次多项式。
- **$B_{35},\dots,B_{39}$ 的支撑全落在 $(x_{38},x_{39})=(143.88,172.09)$ 内，不含任何节点**（$B_{39}$ 支撑 $[162.78,172.09]$，夹紧样条在 $X$ 处只有 $B_{40}=1$）。这五个系数不出现在任何约束和目标中，只受 $c\ge0$；它们应在程序中消去（置 0），这不改变其余变量上的可行集投影（§12 需要）。另 $B_{30},B_{31}$ 在节点上的最大值只有 $8.4\times10^{-4}$、$4.4\times10^{-5}$。节点上可见的系数只有 35 个，减去 $c_1$ 后为 34 个。

**7.4 矩。** 若 $\Pi\propto s^k\ln(-s)$，归一化矩为 $\frac{1}{s_0^{n+k+1}}\int\rho x^n\,dx$（§8），故 $n_\ell=k+1$：S0、P1 为 $k=1$，$n_\ell=2$；D0 为 $k=2$，$n_\ell=3$（(8.72)）。代码以 $X=x_{39}$ 同时作截断和归一（#152）：
$$\mu_{c,n}=\frac{1}{X^{n+n_c}}\sum_{j\le39}\Delta\,x'_j\,\kappa_c(x_j)\,x_j^n\,\hat\rho_{c,j}.$$
pQCD 右端 $2\pi C/(n+n_\ell)$ 与截断尺度无关，所以用 $X$ 代替 $s_0$ 自洽。**注意中点格的边界**（核算）：$j\le39$ 的中点格覆盖 $\phi\in[0,39\Delta]$，即 $x\in(4,189.16]$（1.920 GeV），而归一用 $X=172.085$（1.831 GeV）；$s_0=205.34$ 落在第 40 格 $[189.16,231.33]$ 内。对 pQCD 形状 $\rho\propto x^k$，离散归一矩为目标的 $(189.16/172.09)^{n+n_\ell}=1.208$（$n+n_\ell=2$）至 $1.605$（$=5$）倍。这是作者的离散选择，照搬，不修正。

> **离散公式（本仓库实现）**　$\boxed{\hat\rho_{c,j}=\textstyle\sum_{i=1}^{40}B_i(x_j)\,c_{c,i},\ c\ge0,\ c_{c,1}=1;\ \ \mu_{c,n}=X^{-(n+n_c)}\sum_{j\le39}\Delta x'_j\kappa_c(x_j)x_j^n\hat\rho_{c,j},\ n_c=(2,2,3)}$；$c_{c,35..39}$ 消去。
>
> **与 2309 的关系**：2309 用节点值 $\rho_{\ell,i}$（(3.68)），直接 $\rho\ge|\mathcal F|^2$，矩为 (3.72) $\frac\pi M\sum(ds/d\phi)s^n\rho$，按 $s_0$ 归一（`constraints.moment_row`）。2403 v1 也是节点 ρ；$\kappa$ 因子化与 B 样条是 v3/v4 代码新增。

## 8. SVZ 目标与 Brodsky–Lepage 上界

**8.1 FESR。** $\Pi$ 在沿 $[4,\infty)$ 切开的平面上解析，$\rho=2\mathrm{Im}\Pi(x+i0)$（7.1 节定义）。取围道：$|s|=s_0$ 的圆（逆时针）加割线两岸。若 $s^n\Pi$ 在 0 处正则，则
$$0=\oint_{|s|=s_0}\!s^n\Pi\,ds+\int_4^{s_0}\!x^n[\Pi(x+i0)-\Pi(x-i0)]dx\ \Rightarrow\ \int_4^{s_0}\!\rho\,x^ndx=i\oint_{|s|=s_0}\!s^n\Pi\,ds .$$
正则性决定 $n$ 的下限：S0 为 $n\ge0$；P1 由流守恒 $\Pi^1_1=s\,\Pi_V(s)$，$n\ge-1$；D0 由 $\Pi_2/s^2$ 解析（7.2.2 节末），$n\ge-2$。代入单位算符 OPE $\Pi=-C\,s^k\ln(-s/\mu^2)$，$s=s_0e^{i\theta}$，$\ln(-s/\mu^2)=\ln\frac{s_0}{\mu^2}+i(\theta-\pi)$，$N=n+k+1\ge1$：
$$\oint=-iCs_0^N\!\int_0^{2\pi}\!e^{iN\theta}\big[\ln\tfrac{s_0}{\mu^2}+i(\theta-\pi)\big]d\theta=-iCs_0^N\frac{2\pi}{N}\ \Rightarrow\ \frac{1}{s_0^N}\int\rho\,x^n=\frac{2\pi C}{N},$$
与 $\mu$ 无关，等价于 $\rho_{\rm pert}=2\pi Cs^k$ 从 0 积到 $s_0$。

**8.2 三个流的系数**（$a=\alpha_s/\pi$，$N_c=3$，$N_f=2$）：
- S0 $j=m_q(\bar uu+\bar dd)$：(2.21g) $C=\frac{\sum_fm_f^2}{(2\pi)^4}\frac{3}{8\pi^2}(1+\frac{13}3a)$（$N_fm_q^2\to m_u^2+m_d^2$，#157），$k=1$。
- P1 $j=\tfrac12(\bar u\slashed\Delta u-\bar d\slashed\Delta d)$：用 7.2.1 节的一般 ℓ 公式，$\ell=1$ 时 $\frac{N_c}{2\pi^2}\frac{2}{2}B(2,2)=\frac{N_c}{12\pi^2}$，乘味因子 $\frac14\mathrm{tr}(\tau^3\tau^3)=\frac12$，得 $\frac{1}{8\pi^2}$，与 (2.22g) 一致；再乘 $(1+a)$，$k=1$。
- D0 $j=T^{++}=T^{\mu\nu}\Delta_\mu\Delta_\nu$（能动张量流）：Zoller 结构 $P^{\mu\rho}P^{\nu\sigma}+P^{\mu\sigma}P^{\nu\rho}-\frac23P^{\mu\nu}P^{\rho\sigma}$ 在静止系用 $\bar\Delta\cdot\Delta=-1$、$\Delta^2=0$、$q\cdot\Delta=0$ 收缩，得 $2s^2$，即 $\Pi_2=\frac{2s^2C_0^T}{(2\pi)^4}$；$C_0^T$ 的系数在 $n_g=8,d_R=3,n_f=2,C_A=3,T_F=\frac12$ 下为 $-\frac{8}{10}-\frac{6}{20}=-\frac{11}{10}$ 与 $\frac{24}{18}-\frac{56}{144}=\frac{17}{18}$（核算），$k=2$。

结果即 (2.21h)、(2.22h)、(2.23h)。代码输入 $\alpha_s=0.314588744$，$m_u=3.6/139.57$，$m_d=6.5/139.57$（#155）。目标（核算，30 位，与 MATRICES.md 记录一致）：

| 流 | n | 目标 $T_{c,n}$ | 流 | n | 目标 | 流 | n | 目标 |
|---|---|---|---|---|---|---|---|---|
| S0 | 0 | 3.11259180672e−7 | P1 | −1 | 5.61716677265e−5 | D0 | −2 | 5.13358746336e−5 |
| S0 | 1 | 2.07506120448e−7 | P1 | 0 | 2.80858338632e−5 | D0 | −1 | 2.56679373168e−5 |
| S0 | 2 | 1.55629590336e−7 | P1 | 1 | 1.87238892422e−5 | D0 | 0 | 1.71119582112e−5 |

与 (8.71) 印出的两位数一致。容差（每流 L2，绝对）：代码 $\epsilon=(10^{-7},6\times10^{-6},5\times10^{-6})$（`M3:39`）；**正文 (8.73) 写 P1 $5\times10^{-6}$、D0 $6\times10^{-6}$，两者对调**。主线按代码。7.2.2 节末称 D0 "impose for $n=-2,-1$"，(8.70) 与代码为 $n=-2,-1,0$。输入注记（不算笔误，不调参）：S0 的 $O(\alpha_s)$ 系数 13/3 在正文和代码中一致；按跑动质量 $m(s)$ 的常见 $\overline{\rm MS}$ 系数是 17/3（2505 代码所用），若改用 17/3，S0 目标高 9.3%（3.40e−7）。

**8.3 Brodsky–Lepage。** 6.2 节一般式 $F^1_\ell\simeq-\frac{72\pi\alpha_sf_\pi^2}{t}\frac{C_2}{N_c^2}\{\frac1{\ell+1}+\frac{2(H_\ell-1)}{\ell+2}\}$，取 $\ell=1$、$C_2/N_c^2=4/9$，得 $72\cdot\frac49\cdot\frac12=16$；D0 为 $(24+84)\cdot\frac49=48$（胶子部分加夸克部分，核算）。所以 $|sF^1_1|\to16\pi\alpha_sf_\pi^2$，$|sF^0_2|\to48\pi\alpha_sf_\pi^2$。在 $\sqrt s\sim2$ GeV 尚未进入渐近区，作者乘放宽因子 $r=(0.05,2,6)$（(8.75)），用在 $j=40..50$ 共 11 个节点上（(8.74)，`M3:94`）。核算：$16\pi\alpha_sf_\pi^2=6.870756$，$48\pi\alpha_sf_\pi^2=20.612268$（$f_\pi=92/139.57$），代码常数 6.87076、20.6123 是其舍入（`M3:40`，主线照用）。S0 取常数上界而不取 $1/s$，是作者的选择（(8.74) 前的说明）。

> **离散公式（本仓库实现）**　$\boxed{\big\|(\mu_{c,n}-T_{c,n})_n\big\|_2\le\epsilon_c;\ \ |F_{S0,j}|\le0.05,\ \ x_j|F_{P1,j}|\le2(6.87076),\ \ x_j|F_{D0,j}|\le6(20.6123),\ \ j=40..50}$（`M3:89–94`）
>
> **与 2309 的关系**：2309 (2.56) 含凝聚，$\alpha_s(1.2\,{\rm GeV})=0.4$，$m_u=4$、$m_d=7.3$ MeV，S0 取 $n=0,1$、P1 取 $n=-1,0$，按 $s_0$ 归一；渐近条件 (3.75) 加在 $|\mathcal F|^2\le\epsilon^{FF}$ 上（`constraints.ff_asymptotic_bounds`）。2403 只用单位算符，每流 3 个矩，加 D0，BL 上界直接加在 $F$ 上。

## 9. 手征约束与泛函 𝔉₀、𝔉₁

**树图分波**（由 (2.14)，$\mathcal N=8\pi^2f_\pi^2$）：$T^0=\frac{2s-1}{\mathcal N}$，$T^1=\frac{t-u}{\mathcal N}=\frac{(s-4)\mu}{\mathcal N}$，$T^2=\frac{2-s}{\mathcal N}$，于是 $f^0_0=\frac{2s-1}{16\pi^2f_\pi^2}$，$f^1_1=\frac14\cdot\frac23\frac{s-4}{\mathcal N}=\frac{s-4}{48\pi^2f_\pi^2}$，$f^2_0=\frac{2-s}{16\pi^2f_\pi^2}$。$T^I$ 至多是 μ 的一次式，所以 $\ell\ge2$ 的树图分波全为零。比值 $R_{01}=\frac{3(2s-1)}{s-4}$，$R_{21}=\frac{3(2-s)}{s-4}$（(8.67)），与 $f_\pi$ 无关。

**约束 (8.69)。** 在 $s_*=\tfrac12,1,\tfrac32,2$ 上，第一范数取 8 维向量 $\|(f^0_0-R_{01}f^1_1,\ f^2_0-R_{21}f^1_1)\|_2\le\epsilon^\chi=2\times10^{-3}$；第二范数，**代码**取 12 维向量 $\|(f^0_2,f^2_2,f^1_3)\|_2\le\epsilon^\chi$（`M3:103`：D0、D2、F1），**正文**写 $f^1_1$。证据：(i) 第二范数收集的是树图为零的波，而 $f^1_1$ 树图不为零；(ii) 树图 $f^1_1(s_*)=(-0.01700,-0.01457,-0.01215,-0.00972)$，$\|\cdot\|_2=0.02727\gg2\times10^{-3}$（核算），按正文 Weinberg 点会被排除，这与 §8.4 "ε^χ 取为刚好不排除树图点的最小值" 矛盾；(iii) $f^1_3$ 列在 (8.65) 要计算的六个波中，正文别处却从未用到。主线按代码。

**泛函。** 由 §1 的逆展开，前向振幅为 $T^I(s,t=0)=\sum_\ell2(2\ell+1)f^I_\ell(s)$（$P_\ell(1)=1$）。故 $\mathfrak F_0=2(f^0_0+5f^0_2)(3)$ 是 $T^0(3,0)$ 截到 $\ell\le2$，$\mathfrak F_1=2(3f^1_1+7f^1_3)(3)$ 是 $T^1(3,0)$ 截到 $\ell\le3$（(2.24)）。Weinberg 值：$\mathfrak F_0=2f^0_0(3)=\frac{10}{16\pi^2f_\pi^2}=0.1457434$，$\mathfrak F_1=6f^1_1(3)=-\frac{6}{48\pi^2f_\pi^2}=-0.02914867$（$f_\pi=92/139.57$，核算）。选点 $\mathfrak F_0=\frac{n_i}{50}\times0.163108=0.143535,\ 0.150059,\ 0.156584$（$n_i=44,46,48$，蓝/绿/红，`M3:44–45`）；Weinberg 点对应 $n_i=44.68$。0.163108 不由任何发布的代码算出。

> **离散公式（本仓库实现）**　$\boxed{\|(f^0_0-R_{01}f^1_1,\,f^2_0-R_{21}f^1_1)_{s_*}\|_2\le\epsilon^\chi,\ \ \|(f^0_2,f^2_2,f^1_3)_{s_*}\|_2\le\epsilon^\chi,\ \ \mathfrak F_0=\mathfrak F_0^{\rm fix}\ (\text{消元}),\ \max\mathfrak F_1}$；$f(s_*)$ 行见 §3.5（无 ω、无 Λ、$\zeta<-1$）。
>
> **与 2309 的关系**：2309 (3.63)–(3.64) 的第一范数、取点与 $\epsilon^\chi$ 相同（`constraints.chiral_ratios`）；第二范数是 2403 代码新增的；2309 的泛函是 $(f^0_0(3),f^1_1(3))$（(2.24) 脚注）。

## 10. Λ_ℓ 重标

$\Lambda_\ell=\big(\frac{\sqrt s-2}{\sqrt s+2}\big)^{\ell/2}$（(8.62)，#125）。由 $\sqrt s-2=\frac{s-4}{\sqrt s+2}$ 得 $\Lambda_\ell=\big(\frac{2q}{\sqrt s+2}\big)^\ell=\big(\frac{q}{E+m}\big)^\ell$，阈值附近 $\to(q/2)^\ell$，$s\to\infty$ 时 $\to1$。**v1、v2 笔误**：分母写 $\sqrt s+1$（v1 `QCDsumrule2_30.tex:1965`，v2 `arxiv_submission_3.tex:2008`）；v3、v4 与代码为 $+2$，且只有 $+2$ 才满足 $(q/(E+m))^\ell$ 这一恒等式。

**阈值行为。** $T$ 在 $t=0$ 附近解析，$t=-2q^2(1-\mu)$，所以 $f_\ell=O(q^{2\ell})$，$h=O(q^{2\ell+1})\propto(s-4)^{\ell+1/2}$。从核看：$Q_\ell(1+\frac{x}{2q^2})\simeq\frac{(\ell!)^22^\ell}{(2\ell+1)!}\big(\frac{2q^2}{x}\big)^{\ell+1}$，故 $\mathcal A=O(q^{2\ell})$。重标后 $\mathrm{Im}\hat h=\mathrm{Im}h/\Lambda^2=O(q)$ 与 ℓ 无关，$\tilde h=O(q^{\ell+1})$；对远处的 $x$ 仍有 $(4/x)^\ell$ 压低，这是物理（Lehmann 椭圆）。重标是幺正块的合同变换 $\mathrm{diag}(1/\Lambda,1,1)$，可行集严格不变，只改善条件数。核算：$\rho_1$ 虚部核行 $\max_i|\omega_j\alpha_\ell(x_j;i)|$（原始 / ÷Λ / ÷Λ²）：

| 节点 | ℓ=0 | ℓ=1 | ℓ=2 | ℓ=5 | ℓ=10 | ℓ=19 |
|---|---|---|---|---|---|---|
| $j=1$（$s=4.006$） | 0.077 | 2.5e−7 / 1.3e−5 / 6.8e−4 | 5.7e−11 / 1.6e−7 / 4.2e−4 | 1.4e−21 / 5.2e−13 / 2.0e−4 | 5.0e−39 / 7.2e−22 / 1.0e−4 | 3.4e−70 / 1.4e−37 / 5.5e−5 |
| $j=10$（6.27） | 1.2 | 1.1e−3 / 3.2e−3 / 9.6e−3 | 7.4e−5 / 6.6e−4 / 5.9e−3 | 4.9e−8 / 1.2e−5 / 2.8e−3 | 4.5e−13 / 2.6e−8 / 1.5e−3 | 6.8e−22 / 7.3e−13 / 7.9e−4 |
| $j=39$（172.1） | 2.0 | 7.2e−3 / 8.4e−3 / 9.8e−3 | 1.9e−3 / 2.6e−3 / 3.6e−3 | 3.0e−4 / 6.5e−4 / 1.4e−3 | 3.3e−5 / 1.5e−4 / 7.2e−4 | 1.1e−6 / 2.1e−5 / 3.9e−4 |

结论：÷Λ² 使所有 ℓ 的 $\mathrm{Im}\hat h$ 行落在 $10^{-4}$–$10^{-2}$（原始最小到 $3\times10^{-70}$），锥的线性一侧在求解器容差下可见；$\tilde h$ 行仍小，因为它在锥里以平方出现。$\Lambda_{19}(x_1)^2=6.2\times10^{-66}$，$\Lambda_2(x_1)=3.7\times10^{-4}$，$\Lambda_1(x_{39})=0.858$。Gram 块中 $\mathrm{Im}\hat F=\mathrm{Im}F/\Lambda_{\ell_c}$，$\ell_c=0,1,2$（#131；v1 notebook 对 D0 误用 $\Lambda^0$，见 MATRICES §8.1）。

> **离散公式（本仓库实现）**　$\boxed{\tilde h_n=h_n/\Lambda_n,\ \ \mathrm{Im}\hat h_n=\mathrm{Im}h_n/\Lambda_n^2,\ \ |\tilde h_n|^2\le2\,\mathrm{Im}\hat h_n,\ \ \Lambda_n=\Lambda_{\ell(n)}(x_{j(n)})}$（#127–#129，`M3:63–65`）
>
> **与 2309 的关系**：2309 论文没有 Λ。我们的 2309 SDPB 代码有可选合同 `grid.lambda_rescale`（与行归一并列）；2403 把 Λ 写进问题本身，并用于 Gram 块。

## 11. 作者的离散选择与物理陈述

| 项 | 性质 | 内容（出处） |
|---|---|---|
| 交叉、同位旋、Bose 选择定则 | 物理 | (2.2)、(2.3)；§1 |
| Mandelstam（极大解析）表示 | 物理假设 | (2.1)；§2 |
| $\|S\|\le1$；KKP Gram 半正定 | 物理 | (2.6)、(2.8)；§1、§6 |
| $F(0)=1$（S0 $m_q\partial m_\pi^2/\partial m_q$；P1 电荷；D0 能量） | 物理 | (2.21c)–(2.23c)；§5 |
| $\kappa_\ell=\|\mathcal F/F\|^2$，$\rho\propto q^{2\ell+1}$ 阈值行为 | 运动学（精确） | (8.43)；§7 |
| $4\le s\le16$ 时 $\rho=\|\mathcal F\|^2$ | 物理（未作为等式施加，Watson 步逼近它） | 2.2 节末；§6.4 |
| 单位算符 OPE 到 $O(\alpha_s)$；忽略凝聚 | 物理输入 + 假设 | (2.21h)–(2.23h)、2.6 节；§8 |
| BL $16\pi$、$48\pi\,\alpha_sf_\pi^2/s$ | 物理（$s\to\infty$） | 6.2 节；§8.3 |
| Weinberg 树图比值 | 物理（LO χPT） | (2.15)；§9 |
| $\nu_0=-20$ 作圆盘中心与减除点 | 选择（物理中性，§2） | #19 |
| $M=50$ 个 φ 中点节点；ℓ≤19 截断 | 选择 | #17、#23 |
| 割线上用有限正弦共轭 $K$，割线外用中点求积（两者混用，差 $10^{-21}$–$10^{-9}$） | 选择 | #37、#39–#61；§4.7 |
| $\sigma,\rho,\mathrm{Im}F$ 取节点值；$\rho_2$ 上三角 | 选择 | #46、#78 |
| ρ̂ 用 40 个夹紧三次 B 样条，节点区间 $[0,x_{39}]$ 等距；$c\ge0$（充分条件）；$c_1=1$ 来自 $\hat\rho=F^2$ 的延拓；5 个系数不可见 | 选择 + 延拓假设 | #141–#145；§7.3 |
| Gram 只加在 $j\le39$、只用 S0/P1/D0 | 选择（物理要求所有 $s>4$、所有道） | `M3:83–86` |
| FESR 截在第 39 个节点：归一用 $x_{39}$（1.83 GeV），中点格延伸到 1.92 GeV；$s_0$ 只通过 $n_0$ 进入 | 选择 | #152；§7.4 |
| 每流 3 个矩；L2 绝对容差 $(10^{-7},6\times10^{-6},5\times10^{-6})$ | 选择（"按 $O(\alpha_s)$ 项估计"） | `M3:39` |
| FF 上界只加在 11 个节点上，$r=(0.05,2,6)$ 按可行性选取；S0 用常数界 | 选择 | (8.74)–(8.75) |
| 4 个手征点、L2 范数、$\epsilon^\chi=2\times10^{-3}$；第二范数 | 选择 | (8.66)、(8.69) |
| ℓ4 正则化 $\|\rho/3775\|_4\le100$ | 选择（非物理） | `M3:60` |
| Λ、κ 重标 | 数值中性（严格合同） | §6、§10 |
| 泛函 $\mathfrak F_{0,1}$、$\mathfrak F_0^{\rm fix}=n_i/50\times0.163108$ | 选点选择 | `M3:44–45` |
| 一步 Watson；S0/P1/D0 用 FF 相位目标，其余用旧 $\tilde h$ | 程序选择 | `M3:118–120` |
| 100 位计算核、双精度导出 | 数值 | #21、#162 |

## 12. 对偶证书的先验界（供 V2）

**用途。** 设程序为 $\max b\cdot y$，约束 $F_k(y)=F_{k0}+\sum_ny_nF_{kn}\succeq0$。对任意近似对偶 $Z_k$，记残差 $r_n=b_n+\sum_k\langle F_{kn},Z_k\rangle$，则对每个可行 $y$：
$$b\cdot y\le\sum_k\langle F_{k0},Z_k\rangle+\sum_n|r_n|\,\bar y_n+\sum_k\max\big(0,-\lambda_{\min}(Z_k)\big)\,\tau_k,$$
其中 $|y_n|\le\bar y_n$，$\mathrm{tr}F_k(y)\le\tau_k$（Jansson 型，PLAN §4.2）。下面给出每个 $\bar y_n$ 与 $\tau_k$ 的公式，全部只用核数据（$a_n,b_n,A_n,B_n$ 为 #91 的未重标行；带波浪的 $\tilde a_n,\tilde b_n,\tilde A_n,\tilde B_n$ 为 ÷Λ 后的行，即导出的 A1、B1、A2、B2；带帽的 $\hat b_n,\hat B_n$ 为 ÷Λ² 后的行，即 Bhat1、Bhat2；$\bar v$ 为 (B2)–(B3) 的 $|v|$ 上界向量）。顺序为 ρ → σ → $T_0$ → 各行 → $c,\hat\rho$ → $\mathrm{Im}F$ → 各块。

**(B0) 正则化。** $\|\rho\|_4\le R=100\times3775=3.775\times10^5\Rightarrow|\rho_i|\le R$，$\|\rho\|_2\le3775^{1/4}R=2.959\times10^6$，$\|\rho\|_1\le3775^{3/4}R=1.818\times10^8$。对任一行 $|A_n\!\cdot\!\rho|\le\|A_n\|_{4/3}R$（Hölder，$4/3$ 为对偶指数）。下记 $\alpha_n=\|A_n\|_{4/3}R$，$\beta_n=\|B_n\|_{4/3}R$。

**(B1) 幺正（每行）。** §1：$0\le\mathrm{Im}h_n\le2$，$|\mathrm{Re}h_n|\le1$，$|h_n|\le2$；重标后 $0\le\mathrm{Im}\hat h_n\le2/\Lambda_n^2$，$|\mathrm{Re}\tilde h_n|\le1/\Lambda_n$，$|\mathrm{Im}\tilde h_n|\le2/\Lambda_n$。高 ℓ 阈值处 $1/\Lambda^2$ 可达 $1.6\times10^{65}$，此时用核界（B4）。

**(B2) σ。** $\ell=0$ 虚部行中 σ 的系数只有（§3.5–3.7）：$\mathrm{Im}h^{I=2}_{0,j}=\omega_j\sigma_{2,j}+B_n\!\cdot\rho$，$\mathrm{Im}h^{I=0}_{0,j}=\omega_j(\tfrac32\sigma_{1,j}+\sigma_{2,j})+B_{n'}\!\cdot\rho$。由 (B1) 得
$$\sigma_{2,j}\in\Big[-\tfrac{\beta_n}{\omega_j},\ \tfrac{2+\beta_n}{\omega_j}\Big],\qquad \tfrac32\sigma_{1,j}\in\Big[-\tfrac{\beta_{n'}}{\omega_j}-\bar\sigma_{2,j}^{+},\ \tfrac{2+\beta_{n'}}{\omega_j}-\underline\sigma_{2,j}\Big],$$
其中 $n=(I{=}2,\ell{=}0,j)$，$n'=(I{=}0,\ell{=}0,j)$，$\bar\sigma^{+}_{2,j}$、$\underline\sigma_{2,j}$ 为前式的上下界。

**(B3) T₀。** $I=2,\ell=0$ 的实部行中 $T_0$ 的系数为 $\omega_j$：$|T_0|\le\min_j\big[1+\sum_m|a^{\sigma}_{n,m}|\bar\sigma_m+\alpha_n\big]/\omega_j$。

**(B4) 任一行的核界。** $|\mathrm{Re}\tilde h_n|\le\min\{1/\Lambda_n,\ |\tilde a_n|\!\cdot\!\bar v+\|\tilde A_n\|_{4/3}R\}$，$|\mathrm{Im}\tilde h_n|\le\min\{2/\Lambda_n,\ |\tilde b_n|\!\cdot\!\bar v+\|\tilde B_n\|_{4/3}R\}$，$\overline{\mathrm{Im}\hat h}_n=\min\{2/\Lambda_n^2,\ |\hat b_n|\!\cdot\!\bar v+\|\hat B_n\|_{4/3}R\}$。又 $|\tilde h_n|\le\sqrt{2\overline{\mathrm{Im}\hat h}_n}$。

**(B5) FESR → c 与 ρ̂。** 记 $I_{c,n,j}=\Delta x'_j\kappa_c(x_j)x_j^n/X^{n+n_c}>0$，$m_{c,n,i}=\sum_{j\le39}I_{c,n,j}B_i(x_j)\ge0$。$\|w\|_2\le\epsilon\Rightarrow\mu_{c,n}\le T_{c,n}+\epsilon_c$，且各项非负，故
$$c_{c,i}\le\bar c_{c,i}=\min_n\frac{T_{c,n}+\epsilon_c}{m_{c,n,i}},\qquad \hat\rho_{c,j}\le\bar{\hat\rho}_{c,j}=\min\Big\{\min_n\frac{T_{c,n}+\epsilon_c}{I_{c,n,j}},\ \sum_iB_i(x_j)\bar c_{c,i}\Big\}\ (c_{c,1}=1).$$
$m_{c,\cdot,i}\equiv0$ 的 $i=35..39$ 无界，按 §7.3 消去。核算（双精度）：$\bar{\hat\rho}$ 在 $j=1,20,39$ 处，S0 为 707、80.4、6.26，P1 为 $1.1\times10^4$、437、16.2，D0 为 $2.2\times10^5$、$3.6\times10^3$、99；$i\ge2$ 的 $\bar c$ 中位数为 74、319、2109（S0、P1、D0），最大者都在节点值很小的 $i=31$（$2.6\times10^5$、$6.9\times10^5$、$4.2\times10^6$）。

**(B6) 形状因子。** $j\le39$（Gram）：$|F_{c,j}|\le\sqrt{\bar{\hat\rho}_{c,j}}$，$|\mathrm{Im}F_{c,j}|\le\sqrt{\bar{\hat\rho}_{c,j}}$，$|2\mathrm{Im}\hat F_{c,j}|\le\sqrt{2\overline{\mathrm{Im}\hat h}_n\,\bar{\hat\rho}_{c,j}}$。$j\ge40$（上界）：$|\mathrm{Im}F_{c,j}|\le|F_{c,j}|\le\bar F_{c,j}=(0.05,\ 13.7415/x_j,\ 123.674/x_j)$。另有 $|\mathrm{Re}F_{c,j}|\le1+\sum_m|K_{jm}+K_{0,m}|\,\overline{\mathrm{Im}F}_{c,m}$。

**(B7) 各块的迹。**
- 幺正块 $U_n=\begin{pmatrix}2\mathrm{Im}\hat h_n&\mathrm{Re}\tilde h_n&\mathrm{Im}\tilde h_n\\ \cdot&1&0\\ \cdot&0&1\end{pmatrix}$：$\tau=2+2\overline{\mathrm{Im}\hat h}_n$。
- Gram 块：复 3×3 迹 $\le1+2\overline{\mathrm{Im}\hat h}_n+\bar{\hat\rho}_{c,j}$，实 6×6 嵌入加倍；非对角元界见 (B4)、(B6)。
- FESR 箭形块 $[[\epsilon,w^{\mathsf T}],[w,\epsilon I_3]]$：$\tau=4\epsilon_c$；手征 9×9 与 13×13：$9\epsilon^\chi$、$13\epsilon^\chi$；FF 上界块：$3\bar F_{c,j}$。
- $c_i\ge0$ 的 1×1 块：$\bar c_i$。ℓ4 的 E1 辅助量 $u_i,w_i\in[0,R]$，2×2 块迹 $\le2R$，$\sum w_i\le R$ 的松弛 $\le R$。

**(B8) 计算顺序（V2 实现清单）。**
1. 由导出的 A2、B2、Bhat2 行算 $\|\cdot\|_{4/3}$（只取 ρ 列），乘 $R$ 得 $\alpha_n,\beta_n$。
2. 由 $I{=}2$、$I{=}0$ 的 $\ell=0$ 行按 (B2) 得 $\bar\sigma$，再按 (B3) 得 $|T_0|$ 的界。
3. 按 (B4) 对 1500 行求 $\overline{\mathrm{Im}\hat h}_n$ 与 $|\tilde h_n|$ 的界。
4. 按 (B5) 求 $\bar c$、$\bar{\hat\rho}$；消去 $c_{\cdot,1}$ 与 $c_{\cdot,35..39}$。
5. 按 (B6) 求 $\overline{\mathrm{Im}F}$ 与 $|F|$ 的界。
6. 按 (B7) 汇总每块的 $\tau_k$，连同各变量的 $\bar y_n$ 代入本节开头的不等式。
全部步骤都是区间运算，可在 Arb 中带误差半径完成。

这些界很松（受 $R$ 主导），但全部有限，并可由导出的核逐项计算。若对偶界因此无用，按 PLAN §4.2 的退路：先做全空间求解，再修复对偶。

## 13. 与论文、notebook、计划的分歧汇总

| # | 位置 | 内容 | 证据 | 主线 |
|---|---|---|---|---|
| 1 | (8.12)–(8.13) | 标签 $[0;1;2]$，实际行序 $[0;2;1]$；(8.14) 按 $[0;1;2]$ 是对的 | §3.7；#75–#89 | 代码 $[0;2;1]$ |
| 2 | (8.62) v1/v2 | 分母 $\sqrt s+1$ | v3/v4、#125、$(q/(E+m))^\ell$ | $+2$ |
| 3 | (8.76)、(8.78) | 缺因子 2 | §6.4 秩一推导；`M3:119–120` | 有 2 |
| 4 | (8.73) | P1/D0 容差对调 | `M3:39` | $6\times10^{-6}$/$5\times10^{-6}$ |
| 5 | (8.56) | $\|\rho\|_4\le100$ | `M3:60` 为 $\|\rho/3775\|_4\le100$ | $R=3.775\times10^5$ |
| 6 | (8.68)–(8.69) | 第二范数写 $f^1_1$ | §9 三条证据；`M3:103` | $f^1_3$ |
| 7 | (8.63) | $\mathrm{Im}\hat h=\Lambda^{-2}h$ 缺 Im | §6.2；#129 | $\mathrm{Im}h/\Lambda^2$ |
| 8 | (8.41)、(8.24) | $K_0$ 缺 $\Delta_\phi$ / $d\phi_1$ | §5；#136 | 有 Δ |
| 9 | (8.21)、(8.33)、(8.34) | 上限 ∞→π；$i/2\pi$ 归一；镜像下标差 1 | §4.6 | 记号，无影响 |
| 10 | (8.55) 前后 | "$3ML/2$ 个幺正约束，$L=19$" 为 1425 | 代码 10 波 × 3 × 50 = 1500 | 1500 |
| 11 | (2.23c) | $F^0_2(s)=1$ | 应为 $F^0_2(0)=1$（6.1 节） | — |
| 12 | 7.2.2 节末 | D0 FESR 只取 $n=-2,-1$ | (8.70) 与代码为 $n=-2,-1,0$ | 3 个矩 |
| 13 | notebook #128 及若干式号 | "I=0,1,2" 误标；#30、#38、#153、#158 用旧编号 | §0、§3.7 | — |
| 14 | 计划/PLAN "求和规则截在 $x_{39}$" | 只对归一成立；中点格实际延伸到 $x(39\pi/50)=189.16$（1.920 GeV） | §7.4 核算 | 照搬代码 |
| 15 | 新发现 | $B_{35..39}$ 支撑内无节点，系数不可见；$c_1=1$ 不固定阈值值 | §7.3 核算 | 消去 5 个系数 |
| 16 | 新发现 | $s_*=3$ 的中点行与有限正弦插值行相对差 $2\times10^{-9}$ | §4.7 | 必须用中点式 |
| 17 | 输入注记 | S0 的 13/3 与常见 17/3 | §8.2 | 13/3，不调 |

计划 Task A2 给出的公式（$K_{ij}$、$K_0$、Q 闭式、部分分式、合同变换、因子 2、$n_\ell$、Weinberg 值）经本文推导全部成立；唯一的细化是第 14 条。
