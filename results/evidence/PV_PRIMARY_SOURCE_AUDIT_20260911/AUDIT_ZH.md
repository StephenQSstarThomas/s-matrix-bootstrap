# PV主线的一手资料审计：2309代码可用性与2403数值定义

任务日期：2026-09-11。核查范围限于2309.12402v3官方源包、对应APS页面、作者公开仓库的当前分支/标签/release，以及最早2403公开快照的FESR、散射核、χ约束和密度正则化。没有读取输出曲线来反推参数，没有执行Mathematica/MATLAB/CVX或优化，没有修改本仓库源码或当前物理输入。小数表仅为明确公式的独立算术还原。

核心结论：

1. 作者明确说明2309工作当时未提供公开代码；本次范围内未找到原2309数值producer。这是有范围的检索结论，不是“互联网中确实不存在原代码”的证明。
2. 最早2403公开代码确实采用PV余切矩阵、交叉有理核的节点求和，以及实际独立双密度变量上的四范数界。它为本仓库mixed-PV结构和B=100Mρ的具名选择提供后续一手依据，但不证明2309当时采用同一完整实现。
3. 2403代码的.002是χ容差，FESR另用各通道三矩的归一化残差L2球。能量网格、截止、通道/矩数、误差域及额外IR条件都有变化，不能回填2309参数。不存在统一线性单位换元把当前四条raw盒变成该后续误差域。

来源字节与检索范围见[source_manifest.json](source_manifest.json)；结构化结论、公式及数值见[manifest.json](manifest.json)。两个数值源码仅以只读、不可执行的[GTB_numerics.m.txt](sources/GTB_numerics.m.txt)和[GTB_numerics.nb.txt](sources/GTB_numerics.nb.txt)保存，原字节与换行均未改。

## 1. 原2309数值producer的可用性

[arXiv官方记录](https://arxiv.org/abs/2309.12402v3)确认v3为2024-12-04，原工作首次提交于2023-09-21。本次直接下载[官方TeX源包](https://arxiv.org/src/2309.12402v3)，在内存中只列档：共20个普通文件，为TeX、bbl和图件，没有数值程序。原源包也已按原字节只读保存，未解包执行。SHA-256为：

f8fe776bc387766604c42c60557452abeb08696b6dc9b1ef4eb17bc553b9adfd

对应[PRD出版页](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.096001)和[PRL出版页](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.191601)可访问的链接中，本次未定位数值代码或相关附件入口；这不构成对未链接资源的排除。

作者仓库[README的Paper snapshots](https://github.com/hyfysics/gauge-theory-bootstrap/blob/801684d9ece3de2918a20145178f569459081098/README.md#paper-snapshots)对2309明确说：“No public code was provided at the time.” 同页把2403列为首次公开程序。此次GitHub API返回唯一分支main；完整且未截断的目录树、两个标签与两个release，仅对应2403和2505，没有2309快照。当前main及2505标签commit为801684d9ece3de2918a20145178f569459081098；2403标签commit为c8409331d5df911f592e53e5a43628418a39b213。API响应与README原字节均已存档。

[最早2403 release](https://github.com/hyfysics/gauge-theory-bootstrap/releases/tag/arxiv-2403.10772)发布于2025-09-30，包含GTB_numerics.m/.nb；它是2403论文的程序身份，不能写成2309的原producer。

原2309正文可确认节点与离散变量（Eqs.3.61–62）、FF核（Eqs.3.66–67）及一般支持求解思路；χ约束Eq.3.64只写some norm。§4.3给出tip及手征点附近另两点的选择思路，没有完整精确的代表点算法。此次没有找到能够唯一核定原2309散射交叉离散、χ打包、密度正则化和全部代表选择的数值producer。[原文§3–4](https://arxiv.org/html/2309.12402v3#S3)

## 2. 固定2403快照：变量和归一化

本节源码全部固定于commit c8409331d5df911f592e53e5a43628418a39b213：

- [MATLAB脚本](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m)：9398字节，SHA-256 b12bf411407fc4cb8ae56aea41c88aa336e0f73aa7ae44ebef2b1bca3fc08ad2。
- [Mathematica矩阵producer](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb)：224480字节，SHA-256 9a3a81d7e99ef35355aa5a2bbf10ce9a21e6af929e0cb8371c4afaff133922bd。

下文用R表示电流谱，避免与脚本rho所表示的散射双谱密度混淆。

| 代码变量 | 数学身份 |
|---|---|
| sig | T0及两组单密度，共2M+1个实数 |
| rho | 第一双密度的M²项、第二对称双密度的M(M+1)/2个独立项 |
| specC | S0、P1、D0三通道各NB个非负B-spline系数 |
| spec | 节点值 \(\widehat R_c(s_j)\)，由B-spline矩阵乘specC获得 |
| ImFF | \(F_c(0)=1\)归一化下的FF虚部；ReFF由相同色散核计算 |
| fesr | 已乘回谱运动学因子、再按截止幂归一化的三个矩/通道 |
| w | fesr减去QCD目标；没有再除以QCD目标 |

令\(\beta(s)=\sqrt{1-4/s}\)。NB L4729–4771明确给出

\[
R_c(s)=k_c^2(s)\widehat R_c(s),\qquad
k_{S0}^2(s)=\frac{3\beta(s)}{256\pi^5},
\]
\[
k_{P1}^2(s)=\frac{(s-4)\beta(s)}{384\pi^5},\qquad
k_{D0}^2(s)=\frac{(s-4)^2\beta(s)}{2560\pi^5}.
\]

这里的\(k_c^2\)是谱/FF运动学因子，不是散射幺正性中的\(\pi\beta\)。[NB L4729起](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4729)

在所有量以mπ为单位时，上述变量无量纲。若\(t=m_\pi^2s\)为物理能量平方，则电流谱R的物理量纲为S0/D0的质量四次方、P1的质量二次方；因此物理谱分别乘回\(m_\pi^4,m_\pi^2,m_\pi^4\)。标量物理FF相应为\(m_\pi^2F_{S0}\)，向量及D0的归一化FF为无量纲。

## 3. FESR：误差球如何还原为raw矩

脚本L38的.002是echi。L39的实际FESR容差为

\[
e_{S0}=10^{-7},\qquad e_{P1}=6\times10^{-6},\qquad e_{D0}=5\times10^{-6}.
\]

L89–91将每通道的三个残差作一个Euclidean L2球。矩阶为

\[
N_{S0}=(0,1,2),\quad N_{P1}=(-1,0,1),\quad N_{D0}=(-2,-1,0).
\]

[脚本L38–39](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L38-L39)、[L88–91](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L88-L91)。

设\(\Delta\phi=\pi/M\)，\(s_*=s_{n_0}\)为最后一个不超过名义截止的节点。NB的ints0[n,nl]及随后逐通道乘子，给出精确的代码换元：

\[
I^{\rm raw}_{c,n}=\Delta\phi\sum_{j\le n_0}s'_j\,s_j^nR_c(s_j),
\quad
w_{c,n}=\frac{I^{\rm raw}_{c,n}}{s_*^{n+d_c}}-q_{c,n},
\quad d_{S0}=d_{P1}=2,\ d_{D0}=3.
\]

NB L4798–4857定义这些行，L5149–5157直接导出到intS0/intP1/intD0和SVZSR文件；不存在导出后另一层隐藏误差缩放。[积分定义](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4798-L4857)、[直接导出](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L5149-L5157)。

这里的通道运动学prefactor已经包含在积分行内，不能再乘/除一次。QCD目标同样含原谱定义中的π因子。具体地，令\(\bar m_{u,d}=m_{u,d}/m_\pi\)、\(a=\alpha_s\)，NB L4965–5030给出

\[
q_{P1,n}=\frac{1+a/\pi}{64\pi^5(n+2)},
\quad
q_{S0,n}=\frac{3(\bar m_u^2+\bar m_d^2)(1+13a/(3\pi))}
{64\pi^5(n+2)},
\]
\[
q_{D0,n}=\frac{11/10-17a/(18\pi)}{64\pi^5(n+3)}.
\]

残差是直接相减，不是相对误差\((I-q)/q\)。[QCD目标定义](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4965-L5030)

因此raw域为

\[
\sum_{n\in N_c}
\left[\frac{I^{\rm raw}_{c,n}-s_*^{n+d_c}q_{c,n}}
{e_c\,s_*^{n+d_c}}\right]^2\le1.
\]

各矩的单轴半径为\(e_cs_*^{n+d_c}\)，但其他残差占据预算后，该矩不能任意取满轴界。还原物理单位则有

\[
I^{\rm phys}_{S0,n}=m_\pi^{2n+6}I^{\rm raw}_{S0,n},\quad
I^{\rm phys}_{P1,n}=m_\pi^{2n+4}I^{\rm raw}_{P1,n},\quad
I^{\rm phys}_{D0,n}=m_\pi^{2n+6}I^{\rm raw}_{D0,n}.
\]

按该快照默认值M=50、ν0=−20、mπ=.13957 GeV和名义截止2 GeV，独立公式算术给\(n_0=39\)、\(s_*=172.08533939065126099\)、\(m_\pi\sqrt{s_*}=1.8308974290695215335\) GeV。共有的四个矩对应：

| 矩 | raw单轴半径 | 物理单轴半径 | 半径/QCD目标 |
|---|---:|---:|---:|
| S0,n=0 | .002961336403319563 | 2.188972050400521e−8 GeV⁶ | .321275664171 |
| S0,n=1 | .5096025800151375 | 7.337840139108745e−8 GeV⁸ | .481913496257 |
| P1,n=−1 | .001032512036343908 | 2.011311237464030e−5 GeV² | .106815415010 |
| P1,n=0 | .1776801841991738 | 6.742288156581814e−5 GeV⁴ | .213630830020 |

表格属于后续默认模型，不是对当前1.2 GeV输入的重新定标。该后续口径并非统一“收紧raw .002”：在P1的n=−1方向半径较小，在n=0方向却大得多，而且有三矩L2耦合。

还有一项必须保留的来源差异：2403v4印刷Eq.8.73（PDF第66页）写P1为5e−6、D0为6e−6，与固定MATLAB脚本对调。本证据还原的是固定脚本，未擅自选择或“纠正”二者。[论文Eq.8.72–73](https://arxiv.org/html/2403.10772v4#S8.SS4)

### 为什么不能把当前四盒当作后续误差域

当前声明的四条raw约束是\(|I_n-Q_n|\le.002\)。做同样的单位换元后，它们仍是四个独立区间，只是半宽变为\(.002/s_*^{n+d_c}\)。后续约束则是每通道三维L2球；即使只看每通道前两矩，其误差球投影也是椭圆盘。

仅讨论误差集合时，非退化线性单位换元将盒变为多面体，不能变成球/椭球的光滑边界；单个统一标量更不可能同时改变各矩尺度和分组范数。额外矩、D0及其他物理约束又是进一步的差别。因此没有一个共同数值归一化能将两份完整误差合同等同。

两矩比可确定正谱测度的平均能量平方，不是峰位公式。本核查既不证明原2309的.002必须按后续归一化解释，也不证明当前raw四盒一定就是作者原实现；更未把约822 MeV的矩均值尺度当作修改参数或认定ρ根因的依据。

## 4. 散射kernel、χ范数和密度L4

### 4.1 确认了PV余切＋交叉有理节点求和结构

NB L943–1018构造周期余切序列及M×M核。在一基索引下可写为

\[
K_{ij}=\widetilde K_{i+j-2M-1}-\widetilde K_{i-j},
\quad
\widetilde K_m=\frac{1-(-1)^m}{2M}\cot\frac{m\pi}{2M},
\]

零/可去奇点按离散定义置零，序列按2M周期使用。交叉有理求积块之一是

\[
W_{(j,n),(m,n)}
=\frac{\Delta\phi}{\pi}s'_m
\left[\frac1{s_m+s_n+s_j-4}-\frac1{s_m-\nu_0}\right].
\]

NB还逐项组合其余核生成分波矩阵；不是只导入了一个未说明的外部矩阵。[PV核](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L943)、[交叉块](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L1130)

FF另有减除点转换列
\[
(K_0)_j=\frac{\Delta\phi}{\pi}
\frac{s'_j\nu_0}{s_j(s_j-\nu_0)},\qquad
K_{\rm FF}=K+\mathbf1K_0.
\]
ν0=0时该附加列为零；ν0=−20时不能漏掉。[NB L4464–4525](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4464-L4525)

MATLAB的A1/A2/B1/B2/Bhat1/Bhat2还有\(\Lambda_\ell\)行缩放：
\[
h=\pi\beta f,\quad
\widetilde h=h/\Lambda_\ell,\quad
\operatorname{Im}\widehat h=\operatorname{Im}h/\Lambda_\ell^2,\quad
\Lambda_\ell=\left(\frac{\sqrt{s}-2}{\sqrt{s}+2}\right)^{\ell/2}.
\]
每个有限物理节点上\(\Lambda_\ell>0\)，对应Gram基变换及幺正锥的可逆数值缩放。[NB L4245起、L4332起](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4245)

这些一手定义支持“后续公开程序具有相同类型的mixed-PV结构”。ν0改变实际采样点和有限矩阵，减除坐标、核组合及输入也必须共同对照；没有证明当前保存H等于作者2309 H，更没有转移可行性、支持或连续解析性证书。

### 4.2 χ不是两个独立四维球

MATLAB L97–103把四个S0/P1残差与四个S2/P1残差串成一个八维向量，施加L2≤.002；另对D0、D2、F1在四个阈下点的十二个值施加第二个L2≤.002。后一个条件相对2309原两组比值约束是额外IR输入，前一个分组也不能当作当前separate-L2的数值重标度。[脚本L97–103](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L97-L103)

### 4.3 B=100Mρ确实来自实际独立密度变量

脚本L10、L53和L60给
\[
M_\rho=M^2+\frac{M(M+1)}2,\qquad
\|\rho/M_\rho\|_4\le100
\quad\Longleftrightarrow\quad
\|\rho\|_4\le100M_\rho.
\]

这里除以Mρ，不是除以\(M_\rho^{1/4}\)。M50时Mρ=3775，实际界为377500。NB L2837–2922将ρ2的矩阵列先做两倍对称化，再将对角列除2、取上三角；于是非对角贡献的2在核矩阵内，优化变量保留实际独立密度。若用本仓库C_flat记录，必须先还原非对角实际ρ2再算该四范数。[脚本密度约束](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L53-L60)、[NB对称化](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L2837-L2922)

该代码事实与本仓库B(M)配方相符，是后续来源支持；不是2309作者同样使用B(M)的追认证明。正则化本身限定有限可行域，也不能被视为不影响模型的纯单位变化。

## 5. 可以转写的归一化与不能移用的后续定义

| 类别 | 可确认范围 |
|---|---|
| 纯数值等价 | 固定节点上的正Gram基变换、h/Λ缩放、R与hatR的可逆乘子；矩除以截止幂亦可还原，但必须同时搬运完整误差集合。 |
| 有限实现/约束改变 | ν0=−20导致不同节点；B-spline谱及首系数=1；合并χ球和额外高波IR匹配；三通道各三矩；不同误差球和归一化截止；新的投影泛函。 |
| 物理/近似输入改变 | mπ=139.57 MeV、名义匹配2 GeV、αs=.3145887440、mu/md=3.6/6.5 MeV；代码QCD目标使用上述质量平方和及所列扰动项，未包含原2309冻结的凝聚中心值贡献；FF界亦另有定义。 |

后续脚本的初始投影还使用包含D0/F1的泛函，并非2309的单纯\((f_{00}(3),f_{11}(3))\)平面。本次只辨认其模型身份，没有展开后续选点迭代。[设置与目标](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L5-L45)、[泛函L105–109](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L105-L109)、[QCD输入NB L4924起](https://github.com/hyfysics/gauge-theory-bootstrap/blob/c8409331d5df911f592e53e5a43628418a39b213/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4924)

当前主线物理输入没有改变。本资料适合在报告中说明哪些选择已有后续一手依据、哪些2309来源缺口仍存在；不把一组更接近ρ的输出当作改范数、误差域或质量的证据。
