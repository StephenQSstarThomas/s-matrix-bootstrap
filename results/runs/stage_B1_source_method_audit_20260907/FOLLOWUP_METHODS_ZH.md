# 后续作者程序与方法审计：2403、2505 对 B1 原型机意味着什么

2026-09-07。这里只审计公开来源，没有运行作者优化程序，也没有修改当前核心。**后续程序提供了可复核的成熟锥求解基线，但没有证据表明它们排除了全部连续能量、全部自旋上的幺正违例；也不能把它们的参数当作2309.12402v3的原始设置。** 逐条事实、出处及适用范围保存在 [followup_methods.json](followup_methods.json)。

## 来源与版本

| 来源 | 本次核对版本 | 可确认的代码关系 |
|---|---|---|
| [2403.10772](https://arxiv.org/abs/2403.10772v4) | v4，2025-08-08；v3为2025-02-27 | [官方tag](https://github.com/hyfysics/gauge-theory-bootstrap/releases/tag/arxiv-2403.10772)发布于2025-09-30；本地`.m/.nb`与该tag相同。tag未逐项标明对应哪次arXiv修订，不能擅称v1代码 |
| [2505.19332](https://arxiv.org/abs/2505.19332v3) | v3，2026-02-04；v2为2025-10-02 | [官方tag](https://github.com/hyfysics/gauge-theory-bootstrap/releases/tag/arxiv-2505.19332)发布于2025-10-08；本地两份MATLAB及生成notebook与tag相同。该tag早于v3，不能宣称完整复现v3全部修订 |
| [作者README](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/README.md#L81-L100) | paper snapshots说明 | 明确2309当时没有公开代码；2403是首次公开程序，2505是后续版本。“当时无代码”不是“永远没有代码”，但当前没有2309v3原producer的证明 |

因此，本报告的“源码事实”来自可定位的官方tag；论文说明使用明确的arXiv版本。二者吻合之处可交叉确认，差异则保留，不用后续事实填补2309空白。

## 求解器和数值精度：公开基线并不是自写Newton

2403的执行入口明确是`cvx_begin sdp`和`cvx_solver mosek`；2505保持相同链路。CVX负责把二阶锥、L4范数及Hermitian PSD约束编译给MOSEK。2403v4参考文献列CVX2.1和MOSEK MATLAB手册10.1；这些是文献版本，脚本没有锁定具体安装build，也没有`cvx_precision`或`cvx_solver_settings`覆盖。不能凭空给作者程序补上某个求解容差。[2403源码51–65行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L51-L65)、[2505源码27–45行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m#L27-L45)。

两个生成notebook都设置`p=100`，即Mathematica的100位十进制预计算；2403§8.4明确说明这些核随后乘机器精度优化变量。**核生成精度与优化精度必须分开。** 2505导出JSON前还做`Chop[...,tol]`，其中`tol=10^-50`，并非默认的`10^-10`。这不是严格区间数据产物。[2403本地源码](/home/shiqiu/s-matrix-bootstrap/references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb:520)、[2505精度设置](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L224-L250)、[2505导出清理](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L4373-L4426)。

## 散射核、网格和变量：哪些细节必须先对齐

后续程序将色散核减除于`nu0=-20`，物理采样点为

\[
\phi_j=(j-\tfrac12)\pi/M,\qquad
s_j=\nu_0+\frac{8-2\nu_0}{1+\cos\phi_j}.
\]

角向投影先用解析Legendre第二类函数，例如
\(\mathcal A_\ell(s,x)=Q_\ell(1+2x/(s-4))/[\pi(s-4)]\)；再对谱积分作离散求和。物理节点的奇异部分用离散PV cot核处理，阈下/一般复点则用对应核的有限求和。它不是当前repo“所有直接及交叉槽均采用同一个有限sine-cardinal Cauchy变换，再做角积分”的既定函数族。**同一节点Hilbert矩阵不足以证明这两种有限M处方等价。** [2403§8.1–8.3](https://arxiv.org/html/2403.10772v4#S8.SS1)、[2505 PV构造](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L505-L586)、[Legendre-Q源码](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L751-L790)。

| 源码事实 | 对原型机的直接影响 |
|---|---|
| MATLAB `M=50,l=10`；notebook `L=20`表示`ell<L` | 每个同位旋10个允许奇偶波，最大spin19，共1500个物理disk；不能把三个不同L记号混用 |
| 作者矩阵按`I=0,2,1`分块，再按波、能量排列 | 本repo按能量及`I=0,1,2`组织；导入矩阵必须显式置换。MATLAB的P1索引从`2*l*M+1`开始可直接验证 |
| `sig`有`2M+1`个自由分量 | 减除常数保留自由；公开程序没有另设未减除常数为零 |
| `rho`有`Mrho=M²+M(M+1)/2`个分量 | 全部双密度方向保留，未施加逐分量振幅密度正性 |

出处：[2403设置5–14行、变量53–57行、通道索引68–69行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L5-L69)，[2505同类设置](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m#L25-L38)。

作者还作了与阈值/自旋相关的解析重缩放：

\[
\Lambda_\ell(s)=\left(\frac{\sqrt s-2}{\sqrt s+2}\right)^{\ell/2},\quad
\widetilde h=h/\Lambda,\quad \operatorname{Im}\widehat h=\operatorname{Im}h/\Lambda^2,
\]
\[
\| (\operatorname{Re}\widetilde h,\operatorname{Im}\widetilde h)\|_2
\le\sqrt{2\operatorname{Im}\widehat h}.
\]

这是同一离散幺正锥的重缩放，适合优先作为数值对照。电流3×3 Gram也先换状态基，再消除已知运动学小因子；不是把接近1的原始`S`直接塞给机器精度PSD程序。[2403§8.4式8.60–8.64](https://arxiv.org/html/2403.10772v4#S8.SS4)、[notebook4221–4382行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L4221-L4382)。

## 密度正则化：执行含义已明确，但不能倒灌到2309

两套源码都写

```matlab
Mrho = M^2 + M*(M+1)/2;
norm(rho/Mrho,4) <= 100;
```

因此M50执行的是ordinary L4上界`100*3775=377500`。2403论文式8.56文字称`Mreg=100`，而代码显式对向量先除`Mrho`；不能把论文那个数直接解释成未归一化向量上界100。要复现公开程序，就保留代码表达；要追溯2309，仍须单独证明其归一化。[2403源码10、60行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L10-L60)、[2505源码37–38行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m#L37-L38)。

ρ₂的度量可以从生成器核对：`2*symtr`先把核变成`Kij+Kji`，对角列再除2，最后取上三角。故优化向量存的是**actual ρ₂上三角**，不是本repo的off-diagonal `C_flat=2ρ₂`。也不要混淆振幅的`rho`与电流的`specC/a`；2403的`specC>=0`约束的是电流谱的非负BSpline系数。[2403对称化源码2840–2933行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L2840-L2933)、[电流系数74–75行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L74-L75)。

## 后续改进的范围：不能当作2309参数表

| 项目 | 2403公开快照 | 2505公开快照 |
|---|---|---|
| UV匹配能量 | 2GeV，代码`s0=(2000/mpi)^2` | 同为2GeV |
| 电流 | S0/P1/D0三通道 | 同三通道 |
| 电流谱表示 | 40个三次BSpline，系数非负、首系数1 | 每通道49系数的阈值适配解析参数化；leading系数[1,10]、其余系数绝对值界 |
| FESR | 每通道3矩，组内l2误差 | 6矩、逐分量误差；更新的pinched/高阶pQCD信息 |
| χ约束 | 合并8维S/P误差，另加12维D0/D2/F1误差 | 合并8维S/P误差；本快照未再加上一列的12维约束 |
| 选点/目标 | 固定`F0=2(f00+5f02)`，最大化`F1=2(3f11+7f13)`；一次Watsonian改进 | 先`minimize(0)`，默认5次Watsonian迭代 |

代码位置：[2403 73–111及117–182行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m#L73-L182)、[2505参数44–50行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m#L44-L50)、[2505约束58–121行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/optimize_core.m#L58-L121)。此外，2403式8.73的P1/D0误差标注为5e−6/6e−6，快照变量却赋`eP1=6e−6,eD0=5e−6`；应保留这一待对齐差异，不能自行挑一个版本拟图。

2505的Watsonian步骤用上一解构造新的线性目标，将幅度模长与形状因子相位联系起来；每一步仍解凸问题。它旨在选出趋于饱和的解，并不把`η=1`当作所有能量上的原始强制条件，更不是对任意固定支持方向的最优性证明。论文对“固定点代表物理解”也保留系统不确定性。[2505§2.2](https://arxiv.org/html/2505.19332v3#S2.SS2)、[迭代源码91–121行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/matlab/gtb_qcd_02_optimize.m#L91-L121)、[结论](https://arxiv.org/html/2505.19332v3#S7)。

## 幺正性覆盖与已知局限

源码`all partial waves at all energies`的注释，实际对应有限矩阵的所有行：M50、每同位旋10波；电流PSD仅在`s_j<=s0`的节点上施加。所审查的生成器、优化器和分析notebook中，没有找到额外离节点违例分离、无限高spin包络、global-FG PSD或连续区间认证流程。这是**所审查实现的证据范围**，不是对作者未公开工作作断言。

2403附录B.1做M60、最大spin23、BSpline46、Mreg500等变参稳定性图；这与“所有未采样点都满足幺正性”是不同结论。公开程序没有给出我们现在要求的原式primal/dual间隙证书。可以忠实复现作者有限问题并同时记录外部检查失败，不能把失败隐藏，也不能把加回全部失败约束后的新问题仍称同一个有限原型。[有限离散实现](https://arxiv.org/html/2403.10772v4#S8.SS3)、[附录B.1](https://arxiv.org/html/2403.10772v4#A2.SS1)。

## 方法来源与原型机路线建议

作者的[2103.11484v3，§3](https://arxiv.org/html/2103.11484v3#S3)讨论了双谱反问题的病态、高频近零方向、范数正则化及对偶外界。它解释了为什么“方程全写对”仍可能数值困难，也支持独立保留可验证外界；它不是只检查少量primal波即可获得全域可行性的论证。

另一个离散基实现及其代码约定由[专项审计](DISCRETE_BASIS_ZH.md)单独处理；本报告不将它的截断、速度或参数结论混入官方2403/2505快照。

对当前原型机最直接的建议是：先把目标论文的有限算子、减除约定、度量、行顺序和锥重缩放核对到可运行的成熟conic基线，再讨论求解器强化。把“作者数值原型复现”“加入发现的违例后的加强问题”“连续认证”分层记录。后续Watsonian与电流谱改进应排在其各自科学阶段；不能用它们替代B1固定方向支持验收，亦不能把当前自研Newton的耗时当成作者原型机固有复杂度。
