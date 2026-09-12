# 不借助输出曲线的UV→P波必要下界

本节把[明确IR族的UV排除](GLOBAL_IR_SEED_AND_UV_OBSTRUCTION_20260911.md)中的关键步骤推广到任意原有限UV可行点。只使用原向量电流Gram、F(0)=1所固定的离散色散关系、向量高能FF界和逆矩上界；不使用代表坐标、实验ρ质量、相移或正则化数值拟合。

记低能节点集L={s_j≤s0}，高能集H={s_j>s0}，F_j=1+Σ_i K_ji v_i+i v_j，v_i=ImF_i，mathcal F_j=k_jF_j，R_j为向量电流谱。记w_j>0为低能逆矩权重，Σ_L w_jR_j≤U。高能界给|F_j|≤c_j=√(3×10⁻⁵)/k_j。

对每个h∈H，

\[
\left|\sum_{j\in L}K_{hj}v_j\right|
\ge b_h:=1-c_h-\sum_{j\in H}|K_{hj}|c_j.
\]

取b_h>0。完整Gram的一个2×2主子式给

\[
2k_j^2v_j^2\le R_j y_j,\qquad
 y_j=1-\Re S_j=\kappa_j\Im f_j,\quad0\le y_j\le2.
\]

这条不等式在y_j=0或R_j=0也成立；此时相关v_j=0，无需除以零。加权Cauchy–Schwarz于是严格推出

\[
 b_h^2\le
 \left(\sum_{L}w_jR_j\right)
 \left(\sum_L\frac{K_{hj}^2 y_j}{2w_jk_j^2}\right)
 \le\frac U2\sum_L\frac{K_{hj}^2y_j}{w_jk_j^2}.
\]

因此任何可行P1振幅都满足一个显式线性必要条件

\[
\boxed{\sum_L\frac{K_{hj}^2\kappa_j\Im f_j}{w_jk_j^2}
\ge\frac{2b_h^2}{U}}.
\]

令D_h=Σ_L K_hj²/(w_jk_j²)，则至少有一个低能节点满足

\[
\max_{j\in L}(1-\Re S_j)
\ge\frac{2b_h^2}{U D_h}>0.
\]

先前[精确K／节点的区间记录](../runs/physical_consistency_20260911/GLOBAL_IR_UV_EXCLUSION.json)已经给出原M50、printed/raw/hard条件下全部高节点的b_h、D_h和U。对h=49，右端为该记录lambda_lower乘π，严格大于7.40×10⁻⁵。因此至少一个原生低能P1节点有

\[
\frac{|S_j-1|^2}{4}\ge\frac{(1-\Re S_j)^2}{4}
>1.369\times10^{-9}.
\]

这里的“吸收”指振幅的虚部κImf，不是非弹性损失1−|S|²。上述结果是很弱但严格、且不依赖任何边界选择的吸收／散射强度下界。它说明UV输入为何不能由无限减弱整个P波来满足：FF归一化和高能衰减要求低能ImF承担有限色散贡献，而Gram及有限正谱预算把这转成P波吸收的需求。

这个下界仍远小于90°上穿所需的1−ReS>1（S≠0时）。所以**该不等式自身不证明ρ共振、峰位、峰数或极点**；也不能据此声称完整UV可行域包含或不包含无ρ振幅。论文的几何预选边界、实际相移曲线和五组分辨率检验仍是另一层必须完成的数值结果。

它适用于原有限UV模型及其加强子集。原生K、k、权重和Γ的身份必须相同；不同插值、矩单位或密度定义之间不能转移该数值下界。
