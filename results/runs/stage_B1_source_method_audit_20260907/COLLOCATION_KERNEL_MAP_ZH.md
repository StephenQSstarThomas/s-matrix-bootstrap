# ν₀=0 的 author-collocation 核映射与有限sine核的精确差

2026-09-07。结论：**现有`angular_kernels`、`physical_node_row`、`offcut_row`及`assemble_density_row`，在ν₀=0时，逐项实现了后续作者Legendre-Q／PV-midpoint有限离散式。它们可以直接生成M50/L10的3011×3876算子；数学核本身没有发现需要修正的归一化或ρ₂因子。** 这不证明它们是2309v3未公开producer的唯一恢复，也不转移任何历史rank或可行性结论。当前finite-sine族与这个collocation方案的差可以显式写出，不能再称为有限M下等价。

## 约定与单核

令
\[
\phi_j=(j+\tfrac12)\pi/M,\quad x_j=4/\cos^2(\phi_j/2),\quad
w_j=\frac{\Delta\phi}{\pi}x'(\phi_j)=\frac{x'(\phi_j)}M,\quad b_j=w_j/x_j.
\]

因此`weights`已含谱积分的`1/pi`，无需再乘或除π。现函数返回未乘κ的
\(f_\ell^I=\frac14\int_{-1}^{1}P_\ell T^I\,d\mu\)；幺正性仍在外层使用\(h=\kappa f\)、\(\kappa=\pi\sqrt{1-4/s}\)。

后续作者减除核为\(\pi^{-1}[(x-s)^{-1}-(x-\nu_0)^{-1}]\)。在ν₀=0，crossed／offcut的离散单核就是
\[
g_j(v)=\frac{w_j}{x_j-v},\qquad q_j(v)=g_j(v)-b_j.
\]

`offcut_row(...,subtracted=True)`直接使用q；`False`使用g。物理**原生节点**的direct核由规定的PV／jump代替极点：
\[
q_j(s_k+i0)\rightsquigarrow K_{kj}+i\delta_{kj},\qquad
g_j(s_k+i0)\rightsquigarrow K_{kj}+b_j+i\delta_{kj}.
\]

这里`pv_matrix`与作者周期cot矩阵及奇反射完全一致。作者源码的direct虚部是`(1/2) DiagonalMatrix[vOmega]`，对应S-wave投影的\(c=1/2\)，而非额外的`1/(2pi)`节点因子。纸面δ分布/离散矩阵的写法不能忽略其积分权重。[作者§8.1–8.3](https://arxiv.org/html/2403.10772v4#S8.SS1)、[官方PV源码505–586行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L505-L586)、[direct核及jump源码1537–1624行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2505.19332/theories/qcd/papers/arxiv-2505.19332/src/mathematica/gtb_qcd_01_generate_matrices.nb#L1537-L1624)。

## 角投影与双密度

记\(d=s-4\)、\(\epsilon=(-1)^\ell\)、\(c=\delta_{\ell0}/2\)，则
\[
A_i=\frac{w_i}{d}Q_\ell\!\left(1+\frac{2x_i}{d}\right)
=\frac14\int_{-1}^{1}P_\ell(\mu)g_i(t)\,d\mu,
\qquad J_i=A_i-cb_i.
\]

阈下\(0<s<4\)的Q实分支满足\(Q_\ell(-z)=(-1)^{\ell+1}Q_\ell(z)\)、z>1。`exterior_legendre_q`采用该积分定义的延拓；不是把负轴的支割虚部误加到阈下实振幅。

令\(D_{ij}=x_i+x_j+s-4\)。部分分式恒等式给出
\[
U^{u}_{ij}=\frac{w_jA_i+\epsilon w_iA_j}{D_{ij}},
\]
\[
U^{s}_{ij}=A_iw_j\!\left(\frac1{D_{ij}}-\frac1{x_j}\right)
+\epsilon A_jw_i\!\left(\frac1{D_{ij}}-\frac1{x_i}\right)+cb_ib_j
=\frac14\int P_\ell q_i(t)q_j(u)\,d\mu.
\]

这正是`angular_kernels`的`U`；其中`W[j][i]`的索引是必要的，不可误换成`W[i][j]`。`U_ji=epsilon*U_ij`与奇自旋对角零均是精确交换对称性。mixed direct/crossed的双密度块是\(D_iJ_j\)，其中D表示当前direct核，并不是再插入一次角积分。[作者式8.4–8.12](https://arxiv.org/html/2403.10772v4#S8.SS1)、[现有解析核](/home/shiqiu/s-matrix-bootstrap/src/smatrix_bootstrap/kernels.py:51)。

`assemble_density_row`按\(T^0=3A(s,t,u)+A(t,s,u)+A(u,t,s)\)、\(T^1=A(t,s,u)-A(u,t,s)\)、\(T^2=A(t,s,u)+A(u,t,s)\)组装。结果包含正确的常数系数\(f_{00}=5T_0/2\)、\(f_{20}=T_0\)，以及各单／双密度的转置和奇偶符号。

函数输出的ρ₂列乘**actual对称上三角**：off列为两项之和、diag只计一次。与作者`2*symtr`、diag再/2、取upper的packing一致。写成本repo的C_flat时，系数向量off=2ρ₂，所以**算子off列应除2**；可直接调用`density_row_to_cflat`。rho1保持完整row-major，不作对称化。[作者packing源码2840–2933行](https://github.com/hyfysics/gauge-theory-bootstrap/blob/arxiv-2403.10772/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb#L2840-L2933)、[当前组装](/home/shiqiu/s-matrix-bootstrap/src/smatrix_bootstrap/kernels.py:79)。

## 减除常数与原型行数

g=q+b是每个槽相同的代数平移，故两种坐标由完整可逆三角变换联系：
\[
\sigma_1^s=\sigma_1^u+2R_1b,\qquad
\sigma_2^s=\sigma_2^u+R_1^Tb+R_2b,
\]
\[
T_0^s=T_0^u+b^T\sigma_1^u+2b^T\sigma_2^u+2b^TR_1b+b^TR_2b.
\]

现有`convert_subtraction`正是这项变换；它也适用于本套离散行。作者减除常数是自由变量。**设置未减除T₀=0是额外条件，不是坐标变换。**

M50/L10的原型可用现有函数直接组装：

1. `midpoint_grid(50)`、`pv_matrix(50)`产生原生50点及PV矩阵。
2. 对每节点及30个允许波调用`physical_node_row`，输出1500个复f行，转C_flat后分别写Re/Im，共3000行。角核可按(s,ell)缓存，供I0/I2复用。
3. 用`offcut_row`在s=1/2,1,3/2,2构成8个χ残差行。
4. 加1行未减除常数坐标诊断、2行s=3的f00/f11目标，共3011×3876。

第4步的常数行只是一项坐标诊断；对collocation方案不能未加说明地将它宣称为某个已证明连续振幅的高能极限。若存未减除C_flat，该行是(1,0,…,0)；若存减除坐标，则用上述逆变换的T₀ᵘ线性泛函。

所需最小实现变动是**具名算子生成器、元数据和调度**，不是修改现有数学核：例如`author-collocation-nu0-0`，明确spectral midpoint求和、解析角投影、native PV/jump及actualρ₂度量。保存作者行排列与本repo排列的置换。此处未生成完整M50算子，也未做优化。

## 有限sine与midpoint核的精确桥

取当前同一conformal变量z及\(D_j(z)=1-2z\cos\phi_j+z^2\)。含Nyquist的正弦插值给出
\[
q^{sine}_j(z)=\frac1M\left[2\sum_{n=1}^{M-1}\sin(n\phi_j)z^n+\sin(M\phi_j)z^M\right].
\]

利用\(\cos(M\phi_j)=0\)求有限几何和，得到
\[
q^{sine}_j(z)=\frac{2z\sin\phi_j+\sin(M\phi_j)z^M(z^2-1)}{M D_j(z)},
\qquad q^{mid}_j(z)=\frac{2z\sin\phi_j}{M D_j(z)}.
\]

加回同一个b_j：
\[
H^{mid}_j=\frac{w_j}{x_j-s}=\frac{b_j(1+z)^2}{D_j(z)},\qquad
\boxed{H^{sine}_j-H^{mid}_j=
\frac{\sin(M\phi_j)z^M(z^2-1)}{M D_j(z)}}.
\]

因此单个nodal basis的差一般是\(O(z^M)\)，不能笼统说成\(O(z^{2M})\)。M1、s=3给出可立即检查的反例：\(H_{sine}=4/3\)，\(H_{mid}=8/5\)。M50的几个offcut参数下，本次数学核对得到：

| 参数v | max_j abs(H_sine−H_mid) |
|---|---:|
| 3 | 5.568e−26 |
| −1 | 4.500e−65 |
| −100 | 2.381e−10 |
| −10^6 | .1306211 |

低能目标行的差很小，不能据此推断全能区约束给出的可行域相同；高能散射的crossed参数可落在后两种区域。以上也不证明这一个差异独自解释Fig.4端点的全部偏移。

在物理节点上，H_mid具有极点，而有限sine核有限；桥式是offcut恒等式，不能把其未消去的极点值当作作者的PV／jump。作者collocation的native direct行按上述离散规定与sine一致，crossed/offcut核仍不同。也不能用现有SineSourceRows的连续evaluate来不加声明地验收另一处方的系数；要进行离节点物理检查，必须明确选定并审计该处方的延拓规则。

数学核对：现有19项rational/angular/subtraction独立积分测试通过；M=1,3,50、v=3,−1,−100,−10^6的全部216个桥式差球在384-bit Arb下包含0。只支持以上恒等式与离散行映射，不支持旧rank、旧可行点或2309未公开实现的迁移。
