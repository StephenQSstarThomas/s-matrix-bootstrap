# 有限解析族的阈值幺正必要条件

本结论针对本仓库的未减除有限 cardinal 解析族；不是原文另行给定的输入。目前只推导并独立检查，**没有加入生产优化约束**。高能必要条件成立不能替代这里的阈值条件，也不能替代中间能区的幺正性。

设 d=s−4↓0，h_j(z)=b_j+∑_{n=1}^M D_{nj}z^n，e_j=h'_j(1)。物理割线上 z=(8−s+4i√(s−4))/s，所以

\[
\operatorname{Im}h_j(s)=e_j\sqrt d+O(d^{3/2}).
\]

实系数、有限多项式保证余项在任意固定有限核集合上一致。对交叉变量 t=−d(1−x)/2，h_j(t)在t=0附近解析，其Taylor级数在全部x∈[−1,1]上一致收敛。因此逐项投影合法。

由 t=16z/(1+z)^2 和 Lagrange 反演，对 l≥1有

\[
c_{lj}:=[t^l]h_j(t)=16^{-l}\sum_{n=1}^{\min(M,l)}D_{nj}\frac nl {2l\choose l-n}.
\]

定义 J_{lj}=(1/4)∫_{−1}^1P_l(x)h_j(t)dx。Legendre正交性消掉所有次数小于l的项，而

\[
\int_{-1}^1P_l(x)(1-x)^l dx
=(-1)^l\frac{2^{l+1}(l!)^2}{(2l+1)!}.
\]

此恒等式直接由 Rodrigues 公式分部积分l次得到；(x²−1)^l使所有边界项为零。因此

\[
J_{lj}=g_{lj}d^l+O(d^{l+1}),\qquad
g_{0j}=b_j/2,\qquad
g_{lj}=\frac{(l!)^2}{2(2l+1)!}c_{lj}\quad(l\ge1).
\]

以下仅讨论允许的分波：I=0、2取偶数l，I=1取奇数l；禁戒奇偶的振幅恒等为零，不能套用下面的非零波系数。所有交叉核在物理区为实。按原交叉组合，令R为非对称双密度矩阵、Q为对称双密度矩阵（存储的非对角系数是2Q），得到

\[
\begin{aligned}
B_l^0&=6e^TRg_l+2g_l^TRe+2e^TQg_l,\\
B_l^1&=2g_l^TRe-2e^TQg_l,\\
B_l^2&=2g_l^TRe+2e^TQg_l.
\end{aligned}
\]

对l=0，第一式还应加(3/2)e·σ₁+e·σ₂，第三式还应加e·σ₂。于是 Im f_l^I=B_l^I d^{l+1/2}+O(d^{l+3/2})。对l>0，正交性同时给 Re f_l^I=O(d^l)；对S波令A_I=f_0^I(4)，它是完整C的实线性泛函。

由于κ=π√d/2+O(d^{3/2})，原幺正裕量Δ=2κIm f−κ²|f|²满足

\[
\lim_{d\downarrow0}\frac{\Delta_0^I}{d}
=\pi B_0^I-\frac{\pi^2}{4}A_I^2\quad(I=0,2),
\qquad
\lim_{d\downarrow0}\frac{\Delta_l^I}{d^{l+1}}
=\pi B_l^I\quad(l>0).
\]

因此阈值邻域幺正必需 B₀ᴵ≥πAᵢ²/4，以及每个非零保留波Bₗᴵ≥0。S波条件是凸抛物域，高波条件线性。对任意固定有限波集合，若这些条件全部严格，则由上述极限存在共同的δ>0，使4<s<4+δ上的全部保留波幺正。等号情形必须考察更高阶；本结论不给出δ数值，也不处理无限自旋一致性。

## 独立首谐波检查及适用范围

取M=3、v=(1/2,1,1/2)、λ=1/1024、σ₁=σ₂=λv/8、R=λvvᵀ、Q=λvvᵀ/2、T₀=0。此时∑v_jh_j=1+z，不需要上述D矩阵就可直接展开：

\[
g_l=\frac1{2(2l+1)(l+1)16^l},\quad
A_0=12\lambda,\quad A_2=\tfrac92\lambda,
\quad B_0^0=77\lambda/16,\quad B_0^2=13\lambda/8,
\]

非零波有B_l^0=9λg_l、B_l^1=λg_l、B_l^2=3λg_l。
单一CLI计算 `results/runs/physical_consistency_20260911/threshold_limit_control_evaluate` 在 s=4.000001、4.00000001 直接投影六个波；实际binary64能量作为精确输入，并保存解析包络。它用于检查极限系数和幂次，不是QCD联合可行见证。对应数值比较见 `results/runs/physical_consistency_20260911/THRESHOLD_LIMIT_CONTROL.json`：两能量上的最大相对虚部系数误差约6.0×10⁻⁷、6.0×10⁻⁹；裕量系数误差约7.25×10⁻⁷、7.25×10⁻⁹，与O(d)余项一致。

这里的证明解释为何极小的高波负裕量也可能是严格物理违例：其自然量级本来就是d^(l+1)。不能将统一的绝对浮点容差当作物理可行性标准。
