# 2309→2403散射归一化与密度regulator的换算

**结论：两篇论文对散射A、T、f、S及振幅双谱rho使用同一整体归一化。将2403代码的`norm(rho/Mrho,4)<=100`表达为本方actual-rho L4上限B=100Mrho时，不需要额外乘除pi、2或8pi²。** 唯一须保留的2是rho2对称非对角的存储换元，不是整体幅度单位转换。

本次只核对相关公式和官方notebook的打包，没有查看拟合相移、研究新的regulator值、修改src或运行优化。

## 1. 两篇公式逐项一致

| 内容 | 2309v3本地TeX | 2403v4 |
|---|---|---|
| Mandelstam A | 555–557：单谱前1/pi、双谱前1/pi²，rho1跨s/t及s/u，rho2跨t/u | (2.1)逐项相同 |
| isospin T | 548–550：T0=3A(s,t,u)+A(t,s,u)+A(u,t,s)，T1=A(t,s,u)−A(u,t,s)，T2=A(t,s,u)+A(u,t,s) | (2.2a–c)相同 |
| 分波 | 567：f=1/4∫P_l T | (2.3)相同 |
| S/幺正性 | 575–579：S=1+i*pi*sqrt((s−4)/s)*f=eta exp(2i delta)，abs(S)<=1 | (2.5–6)相同 |
| 双谱对称性 | 560：rho2(x,y)=rho2(y,x) | (2.1)之后同样规定 |

来源：[2309本地TeX](../../../references/2309.12402v3-source/prd_submission_2.tex)、[2403v4 §2.1](https://arxiv.org/html/2403.10772v4#S2.SS1)。例如两套定义都给常数A=T0的f00=5T0/2、f20=T0；不存在全局二倍或pi倍的分波差异。

2403用h=pi*sqrt((s−4)/s)*f作数值变量，随后再令htilde=h/Lambda。对应[notebook786–819](../../../references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb)和4330–4346。pi及Lambda作用在**输出分波/核矩阵**上，不是把振幅双谱rho整体重新定义成pi*rho。

## 2. 减除不改变双谱rho的整体尺度

2403 (8.2–3)把Cauchy核换成`1/(x−nu)−1/(x−nu0)`。展开两个减除核的乘积，原双Cauchy项仍有相同的rho/pi²；新增项至多是单变量积分或常数，可吸收到sigma及T0中。因而减除可以改变T0/sigma坐标，但不产生rho的整体倍率。[2403 §8.1](https://arxiv.org/html/2403.10772v4#S8.SS1)

nu0=−20与本方nu0=0还会改变取样位置/离散问题，这与整体幅度归一化是不同事项；不能由此给B补一个pi或2。

## 3. 官方rho2打包消除了非对角二倍疑问

本地2403 `GTB_numerics.nb`的执行链为：

- 673–701构造转置置换；703–709定义symtr=(I+transpose)/2。
- 2851–2886对rho2核作2*symtr，即非对角列合并为Kij+Kji。
- 2862–2874及2888–2900再把对角列除2，使对角只计一次。
- 2903–2926取含对角的上三角索引。
- 2951–2962把完整rho1和上三角rho2的核拼接，乘MATLAB的rho向量；阈下f核在3693–3735也执行同一规则。

所以官方优化变量rho包含的是全部rho1_ij及一次计数的**实际rho2_ij（i<=j）**，其非对角二倍已经放进核列，不是把变量本身都乘2。

本方输出C_flat对rho2非对角储存2*rho2_ij，因此施加范数前要除2回actual-rho。这是已声明的逐分量坐标变换：

\[
\|\rho\|_4^4=\sum_{ij}|\rho_{1,ij}|^4+
\sum_i|C_{2,ii}|^4+\sum_{i<j}|C_{2,ij}/2|^4.
\]

它不支持把整个B乘2或除2。

## 4. regulator语句的直接含义

[GTB_numerics.m第10、53、60行](../../../references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m)分别定义Mrho=M²+M(M+1)/2、声明rho(Mrho)、施加`norm(rho/Mrho,4)<=100`。等价为

\[
\|\rho_{\mathrm{actual}}\|_4\le100\,M_\rho.
\]

M50因此为377500，无额外整体换算。这里的rho是**振幅双谱**；电流谱rho_l^I、其运动学因子和FESR单位不进入这项正则。

本结论解决的是两篇论文/代码之间的归一化转移问题，不额外证明2309未公开的原图regulator恰好采用2403这一数值配方。没有来源支持通过补一个整体pi或2来改变当前B。
