# 2309v3 定义重审：χ、regulator、FF/FESR与减除

本次重新读本地v3 TeX的定义及当前model/kernels/operators/gauge，不以既有E投影一致性或测试通过替代定义审查。未修改src、未读取archive、未运行优化；仅做一个M3减除代数控制和密度范数控制。

**结论：所列定义未发现足以直接判定为主要数值差异成因的运算错误。当前实现与“已声明的有限模型”一致，不等于2309所有未公开数值选择已唯一恢复。可直接确认的一项错误是准备报告的χ范数元数据硬编码；两种χ分组及regulator来源则是实质性的定义身份问题。**

## 1. χ约束：原文写的是线性残差，不是除以P1后的相对误差

v3 TeX938–951先给比值

\[
R_{01}(s)=3(2s-1)/(s-4),\quad R_{21}(s)=3(2-s)/(s-4),
\]

随后实际施加Eq.(3.64)：

\[
\|f_{00}(s_j)-R_{01}(s_j)f_{11}(s_j)\|\le\epsilon_\chi,
\quad\|f_{20}(s_j)-R_{21}(s_j)f_{11}(s_j)\|\le\epsilon_\chi,
\quad s_j=\tfrac12,1,\tfrac32,2.
\]

`model.py:39–45`、`operators.py:29–37`、`gauge.py:40,74–76` 对应此式，未错除P1。若改成 `|f00/f11−R|≤ε`，就会改变问题，尤其P1接近0时；不是本次审计支持的修正。

把 `f→λf` 后线性残差也乘λ，说明这不是归一化比值误差。原文TeX951明确“some norm”，953强调只匹配四点，不强制全区间线性或精确散射长度。代码也不把fπ固定进χ残差；fπ用于Weinberg参照及选点，这是正确区分。

当前 `model.chiral_slices:11–15` 提供：
- separate-l2：两个四维L2球；
- combined-l2：一个拼接八维L2球。

二者有 `K_combined(ε)⊂K_separate(ε)⊂K_combined(√2 ε)`，不能同名混用。后续官方2403代码 `references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m:25,97–103` 明确采用combined8，另加高波χ条件；该证据支持具名后续分支，不独自证明2309原图也是相同全部设置。

**可直接修正的元数据错误：** `operators.prepare_amplitude:302` 将 `restrictions.chiral_norm` 固定写为 `two separate l2 four-vectors`，未读取实际请求的chiral_norm。准备的H只保存8个线性残差行，不施加球范数；实际gauge通过 `self.groups=chiral_slices(chiral_norm)` 选择组别（gauge:16）。因此这会误标combined准备报告，但不证明现有separate生产计算施加了错误约束，也不能解释其数值差异。

## 2. regulator：当前只管双密度，且是普通向量L4

当前定义准确为

\[
\left(\sum_{ij}|\rho_{1,ij}|^4+\sum_{i\le j}|\rho_{2,ij}|^4\right)^{1/4}\le B,
\qquad B=100d,\quad d=M^2+M(M+1)/2.
\]

没有T0、σ1、σ2，也没有电流谱ρ_current；不是Schatten范数，不是加求积权重的积分范数，不是逐元素上限。`model.density_fourth_power:163–174` 从 `1+2M` 后开始，rho2非对角C_flat除2后取四次幂；`gauge.py:18–24,76` 只把原振幅部分交给此规范；`operators.py:265–266,300–301` 明确B配方。

后续官方MATLAB `GTB_numerics.m:10,53–60` 分开声明 `rho(Mrho)` 与 `sig(Msigma)`，只约束 `norm(rho/Mrho,4)<=100`。由于这是向量L4，等价上限就是100Mrho，而非100Mrho^(1/4)。当前“不含单密度”与这项后续代码一致。

**来源未定：** 2309v3正文给出了T0/σ/ρ变量（TeX553–562、925–929），但没有给出这个regulator范数或数值上限。TeX416只有一个未在正文约束中使用的Mreg宏；1337的regularization是注释掉的书目条目。不能把后续配方称为已唯一恢复的2309原定义，也没有来源根据把σ或T0临时并入范数。

相同配方在M=45/50/60时分别给B=306000/377500/543000；它是同一离散配方，不是固定数值B。不同M结果因此也带有这一数值正则化依赖，不能把变化全部归因于节点分辨率。

独立小控制：M3所有actual双密度为1、C_flat非对角为2，任意大T0/σ时，代码的第四次幂恰为d=15。这确认了实际执行的分量范围和rho2计数，并不证明原作者用了此配方。

## 3. FF与电流谱：归一化对象及质量幂正确

原文TeX703–725给

\[
F_1(0)=1,\quad F_0(0)\simeq m_\pi^2=1,
\quad F_\ell(s)=1+\pi^{-1}\int_4^\infty[1/(x-s)-1/x]\Im F_\ell(x)dx.
\]

F0的1是论文低能近似，F1的1是电荷归一化；都不是F(4)=1或F(s0)=1。当前 `gauge.py:181–185,192–194` 的FF常数1和减除Hilbert核对应此式，不能额外乘2mq或把常数1移到阈值。

按TeX703–704平方，`model.py:81–93` 的

\[
\mathcal F=kF,\quad k_0^2=3\beta/(256\pi^5),\quad k_1^2=(s-4)\beta/(384\pi^5)
\]

正确。原文Gram的共轭排列见TeX762–766；gauge的实PSD只是固定酉变换及正对角合同。其右下变量是电流谱，不能与散射双密度混同。

原文TeX1037–1039明确约束高节点的 **|mathcal F|²**。代码使用 `caps=(2mq² εFF, εFF/2)`、`ffscale=sqrt(caps)`，而 `gauge.py:62–65` 的内部ff已经是mathcal F/sqrt(cap)。因此 `1−||ff||²≥0` 正好是原数值式；给裸F直接套同一个cap或再多除一次cap才是错误，但当前链没有这样做。

ρ_current的物理维数是S0为4、P1为2；四raw矩为6、8、2、4次质量幂。P1使用的是本文带s的真空极化对应谱，不是把常见无量纲矢量谱直接搬入。`model.py:149–160` 明确这些单位，未发现丢掉s或mπ幂的定义错误。

## 4. 四FESR：当前是raw绝对分量误差，不是0.2%精度

原文TeX829–835选(S0,n=0,1)、(P1,n=−1,0)；TeX878–879打印

\[
J_n=s_0^{-n-2}I_n,\qquad I_n=\int_4^{s_0}\rho(x)x^n dx.
\]

当前 `model.py:95–114` 正确将打印目标还原为raw：约(.002385104207,.111838138735,.004228045714,.145711206997)。`gauge.py:59–60,70,77` 定义r=(I−target)/error并逐项要求1−r²≥0，故是四个独立绝对误差区间。

对固定channel、n，TeX1031写的是raw积分减QCD值的标量残差，按字面执行成四项raw绝对误差是相符读法；不能只因容差相对较大就断言程序少除某个量。文中引用的是normalized打印表，转换I=s0^(n+2)J不可省略；代码已转换。未公开作者执行细节可能存在额外归一化／分组，但这属于未恢复的信息，不是已证明的漏因子。

raw εSR=.002对应四目标的相对半宽83.85%、1.79%、47.30%、1.37%。因此不能将其解释成“FESR相对精度0.2%”，也不能由小support gap推出相移／ρ质量相同比例精度。若改成normalized-absolute=.002，会把raw误差乘s0^(n+2)，产生另一个更宽问题；相对误差或联合四维范数也都不是同一问题。

`s0=(1.2/.140)²`、FESR Jacobian=(π/M)ds/dφ与散射Cauchy权重ds/dφ/M均正确（model:116–137；原文TeX1022–1026）。hard-midpoint与clipped-phi是不同截止求积；当前生产选择须以保存的current preparation为准，不能从UVConfig的默认值猜测。clipped仅调整截止单元宽度，不是对未知谱的精确积分。

当前主线使用打印目标，因此算术平均mq与RMS重算Eq.(2.50)的约8.48%／.68%标量差，不直接改变主线已冻结的打印矩；mq_rule仍影响标量FF cap。这不能被混成“当前P1 FESR少了8.48%”。原文先假设等质量再列不同mu/md，如何合成有效mq不是唯一明确给出的规则。

## 5. 减除：σ的移项、T0二次项及返回符号一致

从原式Eq.(2.7)令D=q+b，独立展开得到

\[
\sigma_1^s=\sigma_1^u+2R_1b,\quad
\sigma_2^s=\sigma_2^u+R_1^Tb+R_2b,
\]
\[
T_0^s=T_0^u+b^T\sigma_1^u+2b^T\sigma_2^u+2b^TR_1b+b^TR_2b.
\]

逆式的T0二次项仍为 **+2bᵀR1b+bᵀR2b**：T0u=T0s−bᵀσ1s−2bᵀσ2s+这些二次项。`model.convert_subtraction:297–314` 与之相同；不能把逆变换的全部项统统改负。

`quotient.subtraction_rows:136–158` 接收actual密度行，`subtract_amplitude_coordinates:160–179` 明确先把actual rho2 off乘2交给C_flat函数，再除2返回；`gauge.py:22–24` 同时对行／系数作互逆ds转换。因此没有发现混用这两个坐标约定的实际接线。

M3非对称R1、对称R2、全部非零dyadic系数的独立手算展开与forward helper逐项差包含0，inverse回原系数也逐项包含0。这个控制专门检查σ移项方向、转置和T0符号，不以旧E投影测试替代。密度本身在减除中不变；T0u=0仍是一条额外约束，不是减除恒等式。当前PV生产guard要求free，未发现把T0s=0偷偷当成同义坐标操作的情况。

## 判定边界

本轮支持修正的是χ准备元数据误标；支持继续明确的是χ范数、regulator身份、FESR误差及截止选择。没有找到应直接改动χ残差、π、rho2因子、FF常数或减除符号以“修好曲线”的依据。原作者未公开的定义选择仍可能影响主要差异，但本审计不把可能性当成已证实原因，也不宣称排除了所有算子／优化／选点错误。
