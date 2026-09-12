# 有限解析族的领先高能幺正条件

范围：当前M模式未减除cardinal振幅，T0=0；所有系数实，R任意、Q对称，Q的非对角原存储C等于2Q。结论对每个固定分波成立，严格条件可给任意有限保留分波集合的最终高能幺正性；不提供无限自旋统一界，也未给出有限 crossover 能量。

## 极限与积分换序

写h_j=q_j+b_j，q_j(−1)=−b_j，令a_j=4q'_j(−1)。物理割线上

\[
1+z(s)=(8+4i\sqrt{s-4})/s,
\quad\sqrt{s}\,\operatorname{Im}h_j(s)\to a_j,
\quad\sqrt{s}\,\operatorname{Re}h_j(s)\to0.
\]

对t=−(s−4)(1−x)/2，有

\[
\sqrt{s}\,h_j(t)\to {\sqrt2 a_j\over\sqrt{1-x}}.
\]

因为h_j(z)=(1+z)\widetilde q_j(z)，有限多项式在闭单位盘上一致有界。s≥8时有√s|1+z(t)|≤8/√(1−x)，故右端存在与s无关的可积控制函数。两个交叉核之积乘s也被常数/√(1−x²)控制，可积分。因此相应积分换序合法，而非逐角极限的形式代入。

Rodrigues公式或Beta积分给

\[
\int_{-1}^{1}{P_l(x)\over\sqrt{1-x}}dx={2\sqrt2\over2l+1}.
\]

具体地，对Rodrigues公式分部积分l次，每一边界项在x→1时均为O(√(1−x))、在x→−1时至少有一个(1+x)因子，因此为零。剩余积分为(1/2)_l/(2^l l!)乘∫(1−x)^(-1/2)(1+x)^l dx。换元y=(1+x)/2后，用B(l+1,1/2)，得到√2 Γ(l+1/2)/Γ(l+3/2)=2√2/(2l+1)。

用J_lj=(1/4)∫P_l h_j(t)dx，得到√s J_lj→a_j/(2l+1)。两个交叉核乘积的积分为O(1/s)，故不贡献分波实部的1/√s领先项；它们在物理区为实，也不贡献虚部。

## 由交叉组合得到五项条件

定义

\[
\alpha_1=a\cdot\sigma_1,\quad\alpha_2=a\cdot\sigma_2,
\qquad r=a^TRa,\quad q=a^TQa.
\]

两条S波虚部有

\[
\sqrt{s}\,\operatorname{Im}f_0^0\to\tfrac32\alpha_1+\alpha_2,
\qquad\sqrt{s}\,\operatorname{Im}f_0^2\to\alpha_2.
\]

对非零保留分波，令d_l=2l+1，则

\[
\sqrt{s}\,\operatorname{Re}f_l^I\to{2A_I\over d_l},
\qquad s\,\operatorname{Im}f_l^I\to{2K_I\over d_l},
\]

\[
(A_0,A_1,A_2)=(\alpha_1+4\alpha_2,\alpha_1-\alpha_2,\alpha_1+\alpha_2),
\qquad(K_0,K_1,K_2)=(4r+q,r-q,r+q).
\]

例如I=0的直接×交叉R项给虚部8r/d_l，Q项给2q/d_l；I=1给2(r−q)/d_l。因κ→π，Δ=2κImf−κ²|f|²满足

\[
s\Delta_l^I\to {4\pi\over d_l^2}(d_lK_I-\pi A_I^2),\qquad l>0.
\]

因此全高能幺正性的必要条件是

\[
3\alpha_1+2\alpha_2\ge0,\quad\alpha_2\ge0,
\]
\[
5(4r+q)\ge\pi(\alpha_1+4\alpha_2)^2,
\quad3(r-q)\ge\pi(\alpha_1-\alpha_2)^2,
\quad5(r+q)\ge\pi(\alpha_1+\alpha_2)^2.
\]

最低非零偶波l=2和奇波l=1最严格；满足它们即满足其它固定更高波的领先必要条件。若五式全部严格，由极限定义存在足够大能量，使任意给定有限保留波集合全部Δ>0。等号时次领先项仍可能改变结论，不能用≥0冒称严格的最终高能证书。

## 凸性、种子与对偶

r和q对双密度系数线性，α对单密度线性。后三式都是线性项不小于线性项平方，构成凸集合；不是按ρ峰加入的经验约束。可用五个形式2y−x²≥0表示：前两个取x=0，后面取x=√(π/d_l)A_I、y=K_I/2。正的尺度变换x→x/√w、y→y/w不改集合。

对闭抛物域2y≥x²，支持函数为

\[
\sup (u x+v y)=-{u^2\over2v}\quad(v<0),
\]

v=0时仅u=0给有限值0，其余无有限支持。由−log(2y−x²)的梯度产生v=−2μ/(2y−x²)<0。验回应把实际u、v乘原解析行计入站立性残差，并加上述支持常数；不能沿用圆盘的支持公式。

正的IR种子已存在：取σ1=σ2=λ sinφ/8、R=λ vvᵀ、Q=λ vvᵀ/2、T0=0。由于a·v=4，α1=α2=λ/2、r=16λ、q=8λ，对足够小λ>0五式均严格。再结合有限所有散射盘及χ/L4的小缩放构造，可作为完整Phase I的硬IR内点。它不证明加入UV后必定可行；UV仍须完整求解、原式及解析核验。

端点切片的第一个完整联合见证确实违反后三式，数值包络见`results/runs/physical_consistency_20260911/endpoint_witness_analytic/source_audit.json`。这是当前点的严格高能障碍，不是整个端点切片无解。新增五式后必须另立模型记录并重新求见证。

## 独立有限能量检查

`results/runs/physical_consistency_20260911/TAIL_LIMIT_CONTROL.json`使用手工定义的M3首谐波v=(1/2,1,1/2)、λ=1/1024，完整解析投影分别在s=10^6、10^8求值。α1=α2=λ/2、r=16λ、q=8λ由h=1+z直接给定，未从同一导数矩阵计算期望值。六个保留波的相应缩放虚部接近所推导极限，误差随能量提高明显下降。这只是独立数值交叉检查；上述控制收敛与积分计算才是极限证明，不由两个能量样本代替。
