# 从闭式振幅证明全能区IR非空，以及UV约束的非平凡性

这是一个辅助定理，**不替代Fig.7的绿色边界代表，也不计作B–F的数值复现**。它回答底层函数族是否连IR约束本身都不自洽，并说明归一化FF、高能界和谱矩为什么会排除足够弱的无ρ振幅。结论关于全部分波的弹性块收缩，不构造其它散射通道的完整解析幺正S矩阵。

## 1. 一个属于当前解析族的明确函数

取

\[
h(s)=1+z(s)=\frac4{2+\sqrt{4-s}},
\]
\[
A_\lambda(s,t,u)=\lambda\left[
\frac{h(s)}8+\frac{h(t)+h(u)}8+
h(s)(h(t)+h(u))+\frac12h(t)h(u)\right],\qquad\lambda>0.
\]

对任意M，当前完整cardinal族都包含此函数：令v_j=sinφ_j，T0=0，σ₁=σ₂=λv/8，R=λvvᵀ，Q=λvvᵀ/2。完整DST的第一谐波恒等式给∑v_jh_j=h；Q非对角的原存储系数是2Q。这里使用明确的闭式函数，不从有限H证书转移连续结论。

h在规定割平面解析、满足实解析性；A在t↔u下对称，按原同位旋交叉组合得到完整pion振幅。对t<4有正Stieltjes表示

\[
h(t)=\int_4^\infty\frac{d\mu(x)}{x-t},\qquad
d\mu(x)=\frac{4\sqrt{x-4}}{\pi x}\,dx.
\]

令x=4+u²，积分化为(8/π)∫u²/[(u²+4)(u²+4−t)]du；部分分式积分即得上述h。它给出所需正性及绝对可积控制，随后由解析延拓固定物理边界值。

## 2. 对任意分波统一控制交叉投影

设d=s−4>0，t=−d(1−x)/2，定义

\[
J_l=\frac14\int_{-1}^1P_l(x)h(t)dx,
\qquad U_l=\frac14\int_{-1}^1P_l(x)h(t)h(u)dx.
\]

正Stieltjes表示给J_l=∫dμ(a) Q_l(1+2a/d)/d>0。这里Q_l(y)>0可直接由Rodrigues分部积分证明：

\[
Q_l(y)=2^{-l-1}\int_{-1}^1\frac{(1-x^2)^l}{(y-x)^{l+1}}dx>0\quad(y>1).
\]

由|P_l|≤1，又有0<J_l≤J₀。部分分式

\[
\frac1{(a-t)(b-u)}=\frac{(a-t)^{-1}+(b-u)^{-1}}{a+b+d}
\]

表明奇l时U_l=0，偶l时0≤U_l≤2J_l，因为内部剩余积分是h(−a−d)≤1。交换积分由正测度和|P_l|≤1的绝对控制保证。

物理割线上

\[
|h(s)|=\frac4{\sqrt s},\quad
\operatorname{Im}h(s)=\frac{4\sqrt{s-4}}s,\quad
J_0=\frac4{s-4}\left[\sqrt s-2-2\log\frac{\sqrt s+2}4\right].
\]

所以J₀≤4/(√s+2)，且J₀≥h(−d)/2=2/(√s+2)。因此

\[
\frac{\operatorname{Im}h}{\kappa J_l}\ge\frac1\pi,
\qquad |h|\le2,\qquad |h|\le4J_0.
\]

## 3. 全部允许分波的幺正不等式

直接使用三种交叉组合可得（I=0,2取偶l，I=1取奇l）

\[
f_l^0=\lambda\left[\frac5{16}\delta_{l0}h+
(\frac54+9h)J_l+\frac72U_l\right],
\]
\[
f_l^2=\lambda\left[\frac18\delta_{l0}h+
(\frac12+3h)J_l+2U_l\right],
\qquad f_l^1=\lambda hJ_l.
\]

从上节界限，统一有

\[
|f_l^I|\le\lambda C_IJ_l,
\qquad \operatorname{Im}f_l^I\ge\lambda c_I\operatorname{Im}h\,J_l,
\quad (C_0,C_1,C_2)=(55/2,2,11),\quad(c_0,c_1,c_2)=(9,1,3).
\]

S波中额外的直接h项，其虚部为正；实模界用|h|≤4J₀吸收。于是

\[
1-|S_l^I|^2\ge
\kappa\lambda J_l\left(2c_I\operatorname{Im}h-\kappa\lambda C_I^2J_l\right)>0
\]

只需

\[
0<\lambda<\min_I\frac{2c_I}{\pi C_I^2}
=\frac{72}{3025\pi}.
\]

同一个λ因此对所有s>4、所有允许l有效；不是逐个l取一个不同的小参数。阈值、s→∞或l→∞时裕量可以趋零，未声称存在统一正裕量。

同时|f_l^1|≤λ，故|S_l^1−1|≤πλ<1。P1的S始终在右半平面，连续相移没有90°上穿。这是明确的无ρ式相移族，而不是用图像挑出的无峰点。

## 4. 同时满足给定IR范数

在四个手征输入点及其角积分中，s,t,u∈[0,4)，所以1≤h<2，从而|A_λ|≤43λ/4。利用同位旋系数及∫|P₁|=1，有

\[
|f_0^0|\le215\lambda/8,\quad
|f_1^1|\le43\lambda/8,\quad |f_0^2|\le43\lambda/4.
\]

四点上|R₀₁^χ|≤9/2、|R₂₁^χ|≤9/7，故两条独立L2残差分别不超过817λ/8和989λ/28。取λ<8εχ/817即可使两球严格满足。实际双密度又满足‖r‖₄≤λ[M²+M(M+1)/2]^(1/4)，因此再取λ小于B除以此因子即可。

这证明当前解析函数族在原IR不等式层面有明确的全能区、全允许分波例子；不是完整QCD或完整多通道S矩阵的存在定理，也没有得到物理fπ位置的边界代表。

## 5. 原UV输入为何排除该族的弱耦合部分

只用向量通道。记F_j=1+∑K_ji v_i+i v_j，v_i=ImF_i，\(\mathcal F_j=k_jF_j\)。低能集合L取s_j≤s0，高能集合H取s_j>s0；向量高能上界给|v_j|、|ReF_j|≤c_j=√ε_F/k_j（ε_F=3×10⁻⁵）。

对任意h∈H，归一化与三角不等式给

\[
\left|\sum_{j\in L}K_{hj}v_j\right|\ge
b_h:=1-c_h-\sum_{j\in H}|K_{hj}|c_j.
\]

若b_h>0，这要求低能ImF不能全小。另一方面，本族有1−ReS_j≤πλ。原Gram相应2×2主子式给

\[
R_j\ge\frac{2k_j^2v_j^2}{1-\operatorname{Re}S_j}
\ge\frac{2k_j^2v_j^2}{\pi\lambda}.
\]

向量逆矩上界U及正权重w_j由此限制∑_L w_j k_j²v_j²≤πλU/2。加权Cauchy–Schwarz于是得到必要条件

\[
\lambda\ge\frac{2b_h^2}{\pi U D_h},\qquad
D_h=\sum_{j\in L}\frac{K_{hj}^2}{w_jk_j^2}.
\]

对原M50、打印矩/raw绝对.002/hard-midpoint、s0=(1.2/.14)²及向量FF上界，使用精确原生节点、K和
k_j²=(s_j−4)√(1−4/s_j)/(384π⁵)、w_j=(π/M)tan(φ_j/2)，逐项区间计算给第49号（零基）高节点的

\[
b_h>0.85464142575,\qquad
\lambda\ge2.35634303968\times10^{-5}.
\]

而λ=2⁻³⁰满足上述全部严格IR条件，却违反此UV必要下界。因此该明确完整解析IR振幅没有原有限UV扩展；实际上整个上述保证手征可行的小λ区间都被这个界排除。全部高节点、公式输入与包络见[计算记录](../runs/physical_consistency_20260911/GLOBAL_IR_UV_EXCLUSION.json)。没有输入任何ρ质量或论文输出曲线。

## 6. 此定理闭合了什么

它严格连接了“共同解析、交叉及分波收缩的IR族存在”“IR不强制ρ式上穿”和“归一化FF＋高能FF界＋有限谱预算会排除其中的弱相互作用部分”。它没有证明所有无ρ振幅都被UV排除，更没有从一个逆矩推出唯一ρ质量。

实际Fig.7绿色边界振幅的全部C另有原式及解析Farkas排除证据；完整多变量UV代表、图形比较和分辨率检查也另行交付。这个辅助解析例子不替代那些主线结果。
