# Original-center midpoint/PV amplitude rows

This chapter specifies the amplitude operator suggested by the archived
2024 Mathematica kernels, with the **2023 map center restored to zero**.
It retains all 3876 M50 amplitude variables and the original normalization
`f=(1/4) integral P_l T`, `S=1+i pi beta f`. It imports no later density
regularizer, three-current setup, B-spline parametrization or matching scale.

The result is a declared **discrete midpoint/PV operator**. Its assigned
physical-node values are not literal limits of its off-cut rational
functions at their poles. It must remain separate from the exact
finite-sine analytic family of [chapter 69](69_full_analytic_rank_and_finite_target_plane.md).
No M50 operator-generation, rank or optimization worker was run here.

## 1. Nodes, weights and primitive density coordinates

Use zero-based i,j,k=0,...,M-1 and

\[
\phi_i=\frac{\pi(i+1/2)}M,\quad x_i=4\sec^2(\phi_i/2),\quad
w_i=\frac{\Delta\phi}{\pi}x'_i=\frac{x'_i}M,
\quad b_i=\frac{w_i}{x_i}=\frac1M\tan(\phi_i/2). \tag{PK1}
\]

The weights already include the dispersion factor 1/pi. Density variables
are unweighted sampled values, in the order

```
T0; sigma1_i; sigma2_i; rho1_ij (row major); rho2_ij (i <= j).
```

R1=(rho1_ij) is unrestricted and R2=(rho2_ij) is symmetric. The upper
entries are the actual symmetric-density values: an off-diagonal R2_ij
multiplies the **sum** of its two ordered products, while a diagonal entry
multiplies one product. Thus there are
`1+2M+M^2+M(M+1)/2=3876` variables at M50.

For the subtracted convention used below, the off-cut scalar kernel is

\[
h_i(v)=w_i\left(\frac1{x_i-v}-\frac1{x_i}\right). \tag{PK2}
\]

The unsubtracted kernel is `g_i(v)=w_i/(x_i-v)=h_i(v)+b_i`.
The formal amplitude construction uses

\[
\begin{aligned}
A(s,t,u)={}&T_0+\sigma_1^T h(s)+\sigma_2^T[h(t)+h(u)]\\
&+\sum_{ij}R_{1,ij}h_i(s)[h_j(t)+h_j(u)]
+\sum_{ij}R_{2,ij}h_i(t)h_j(u). \tag{PK3}
\end{aligned}
\]

The last sum is over the full symmetric matrix, even when stored in upper
coordinates. This is why importing a half-symmetrized polynomial row
requires multiplying its off-diagonal rho2 columns by two.

## 2. Physical PV and jump versus off-cut evaluation

Define the periodic sequence, modulo 2M,

\[
a_m=\frac{1-(-1)^m}{2M}\cot\frac{\pi m}{2M},\qquad a_0=0,
\qquad K_{ki}=-a_{k-i}+a_{k+i+1}. \tag{PK4}
\]

Even m entries vanish exactly. The folded half-circle matrix obeys

\[
K\sin(n\phi_i)=\cos(n\phi_k),\quad n=1,\ldots,M,
\]

including the highest mode, whose cosine samples vanish. The highest
imaginary mode is retained. In general the diagonal K_kk is nonzero.

The integrated direct-channel rows are

| Convention | Off cut, including 0<s<4 | Assigned physical node k |
|---|---|---|
| Subtracted at zero | `d_i(s)=w_i[(x_i-s)^-1-x_i^-1]` | `d_i=K_ki+i delta_ki` |
| Unsubtracted | `d_i(s)=w_i/(x_i-s)` | `d_i=K_ki+b_i+i delta_ki` |

The unit imaginary jump follows directly from the 2023 dispersion
normalization: `(1/pi) integral sigma(x)/(x-s-i0)` has imaginary part
sigma(s). The notebook independently confirms it: `his` is one half
times `DiagonalMatrix[vOmega]`, with no extra M or pi divisor. The
unintegrated delta notation in a manuscript kernel must not be inserted
again after the quadrature measure has already been absorbed.

The folded PV matrix is the notebook's K at its physical-node step. Its
off-cut `Kunphys` is the midpoint expression (PK2). At a physical pole
the module uses the assigned node row; it never evaluates (PK2) there or
silently substitutes a small imaginary regulator.

## 3. Closed angular kernels and real branches

Write `d=s-4`, `epsilon=(-1)^ell`, `c_ell=delta_(ell,0)/2`, and

\[
z_i(s)=1+\frac{2x_i}{d},\qquad
\mathcal A_i^\ell(s)=\frac{w_i}{d}Q_\ell(z_i). \tag{PK5}
\]

Q is defined by the real exterior integral

\[
Q_\ell(z)=\frac12\int_{-1}^1\frac{P_\ell(\mu)}{z-\mu}\,d\mu,
\quad |z|>1.
\]

For physical s>4, z_i>1. For every 0<s<4, including all five source
energies 1/2,1,3/2,2,3, z_i<-1. The continuation is
`Q_l(-z)=(-1)^(l+1) Q_l(z)`, with no added imaginary Ferrers term.

The single crossed-channel moments are

\[
J_i^\ell=\frac14\int P_\ell h_i(t)\,d\mu
=\mathcal A_i^\ell-c_\ell b_i,\qquad
\frac14\int P_\ell h_i(u)\,d\mu=\epsilon J_i^\ell. \tag{PK6}
\]

Define the notebook subtraction kernels, including all midpoint weights,

\[
W_{ij}(s)=w_i\left[\frac1{x_i+x_j+s-4}-\frac1{x_i}\right],
\quad W_{0,i}=b_i,\quad W_{00,ij}=b_i b_j. \tag{PK7}
\]

Then the double crossed moment is

\[
U_{ij}^\ell=\frac14\int P_\ell h_i(t)h_j(u)\,d\mu
=\mathcal A_i^\ell W_{ji}
+\epsilon\mathcal A_j^\ell W_{ij}+c_\ell b_i b_j. \tag{PK8}
\]

This follows by partial fractions from

\[
\frac1{(x_i-t)(x_j-u)}
=\frac{(x_i-t)^{-1}+(x_j-u)^{-1}}{x_i+x_j+s-4}.
\]

There is no singular i=j case. `U_ji=epsilon U_ij`; odd diagonal entries
are identically zero. For the unsubtracted convention simply use

\[
J_i^\ell=\mathcal A_i^\ell,\qquad
U_{ij}^\ell=\frac{w_j\mathcal A_i^\ell+
\epsilon w_i\mathcal A_j^\ell}{x_i+x_j+s-4}. \tag{PK9}
\]

At s=4 the unsubtracted single limit is `c_l w_i/x_i`. All subtracted
crossed single and double moments vanish exactly. The implementations
treat this limit separately instead of dividing by s-4.

The Arb implementation evaluates Q for positive exterior z using

\[
Q_\ell(z)=\frac{2^\ell(\ell!)^2}{(2\ell+1)!}\,z^{-\ell-1}
{}_2F_1\left(\frac{\ell+1}2,\frac{\ell+2}2;
\ell+\frac32;z^{-2}\right), \tag{PK10}
\]

and uses the exact negative-argument parity above. At ell=0 this reduces
to `atanh(1/z)`, fixing the normalization. No source angular quadrature
appears in these closed-form rows.

## 4. Complete isospin rows

Let Abar be the partial-wave projection of A(s,t,u) for one density
coordinate and Bbar that of A(t,s,u). The third permutation contributes
epsilon Bbar. The entire column table is

| Density coordinate | Abar | Bbar |
|---|---|---|
| T0 | c_l | c_l |
| sigma1_i | c_l d_i | J_i |
| sigma2_i | `(1+epsilon) J_i` | `c_l d_i+epsilon J_i` |
| rho1_ij | `(1+epsilon) d_i J_j` | `d_j J_i+U_ij` |
| rho2_ii | U_ii | `epsilon d_i J_i` |
| rho2_ij, i<j | `(1+epsilon) U_ij` | `epsilon(d_i J_j+d_j J_i)` |

Apply the original 2023 isospin definitions:

\[
f^0=3\overline A+(1+\epsilon)\overline B,\quad
f^1=(1-\epsilon)\overline B,\quad
f^2=(1+\epsilon)\overline B. \tag{PK11}
\]

Only even ell enter I0/I2 and odd ell enter I1; every forbidden-parity
row is zero. In particular the constant column is 5/2 in f00, 1 in f20,
and zero in the other retained waves. Multiplication of a physical row
by `Omega_j=pi sqrt((s_j-4)/s_j)` gives G, not f. The notebook's physical
Ap/Am arrays include Omega; their unphysical counterparts do not.
The later Lambda cone rescalings are not applied here.

Equations (PK4)–(PK11) provide every physical-node and subthreshold row
without changing variable count or angular normalization. They also avoid
depending on potentially confusing I1/I2 labels of a plus/minus stack:
the original definitions and the explicit row permutation control it.

## 5. Subtraction and polynomial-coordinate maps

Writing superscripts u,s for unsubtracted/subtracted parameters and
keeping R1,R2 unchanged, the finite triangular map is

\[
\begin{aligned}
\sigma_1^s&=\sigma_1^u+2R_1b,\\
\sigma_2^s&=\sigma_2^u+R_1^Tb+R_2b,\\
T_0^s&=T_0^u+(\sigma_1^u)^Tb+2(\sigma_2^u)^Tb
+2b^TR_1b+b^TR_2b. \tag{PK12}
\end{aligned}
\]

This works for both off-cut and assigned node rows when the b shift in
the direct channel is carried consistently. It fixes no parameter and
is invertible. For example the unsubtracted infinity functional in
subtracted coordinates is

\[
A_\infty=T_0^s-(\sigma_1^s)^Tb-2(\sigma_2^s)^Tb
+2b^TR_1b+b^TR_2b.
\]

An exact identification with the old polynomial **coefficient labels** is
especially simple in the subtracted nodal coordinates. Let
`T_in=sin(n phi_i)`, n=1,...,M. Then

\[
T_0^s=c_z,\quad \sigma_1^s=T a_z,\quad\sigma_2^s=T b_z,
\quad R_1=T B_zT^T,\quad R_2=T C_zT^T, \tag{PK13}
\]

where the old C_z matrix has off-diagonal entries equal to half its flat
coefficient. This is a block-diagonal invertible parameter transformation.
It does not identify the two finite operators: under the midpoint map,
the off-cut z^n function becomes

\[
r_n(z)=\frac{z^n+z^{2M-n}}{1+z^{2M}},\qquad 1\le n\le M,
\]

while the separately assigned physical-node value remains z_j^n. Between
nodes the rational expression is real on the physical circle and has
poles at the nodes. Thus the new operator is a product-integration/PV
prescription, not a newly established globally analytic finite amplitude.
The theorem about the old exact finite-sine operator does not transfer
without a new rank/operator analysis.

## 6. Source-code trace, implementation and tests

[The exact input-cell excerpt archive](../pv_midpoint_kernel_source_excerpts.json)
contains byte offsets, literal source fragments and hashes for 18 notebook
symbols. The pinned notebook SHA256 is
`9a3a81d7e99ef35355aa5a2bbf10ce9a21e6af929e0cb8371c4afaff133922bd`, at author
commit `801684d9ece3de2918a20145178f569459081098`. Relevant cells are:

| Source cells | What was checked |
|---|---|
| K; W0/W/W00 | Folded PV signs and already absorbed midpoint/dispersion weights |
| Ap1/Am1; Ap1unphys | Q branches, angular normalization and physical Omega factor |
| Kunphys/Wunphys/W00unphys | Direct off-cut and subtraction kernels |
| hrs/his; hrst/hist | PV and unit imaginary jump after integration |
| hrI1; hrI2b | Isospin combinations; twice symmetrization followed by halving diagonals |
| hrIsigma; listindex; fsigma | Constant/single column order and physical/unphysical row ordering |

The notebook orders physical rows by I0-even, I2-even, I1-odd, then spin
and energy. The repository's canonical order is energy, then I0-even,
I1-odd, I2-even. The permutation must be explicit. The notebook uses
L=20 before parity selection; the 2023 ten-waves-per-isospin convention
corresponds to even 0–18 and odd 1–19.

The new [kernel module](../runs/repo_reorganization_20260906/recovery.json) (historical member: `src/smatrix_bootstrap/pv_midpoint_kernels.py`)
exposes `midpoint_grid`, `pv_matrix`, `angular_kernels`, `density_labels`,
`physical_node_row`, and `offcut_row`. The last API deliberately accepts
only 0<s<4; off-grid physical boundary continuation is not invented.
All physical row values have Acb type; subthreshold rows have Arb type.

Twenty tests passed against independent high-precision direct angular
integration, including all five subthreshold energies, full raw isospin
assembly, both kinds of rho2 entries, independent sine-sum PV values,
the unit jump, the Nyquist mode and threshold limits. A further test
passes the complete unsubtracted/subtracted coordinate transformation
for physical and subthreshold rows. These are small numerical regression
checks of analytically derived formulas, not a full-source rank proof.

An independent angular implementation uses a positive series/log
recurrence and separately validated Gauss integration. Its bounded M2
comparison supplies 868 interval overlap/identity checks over the five
subthreshold energies and s=9, ell=0,...,5; all pass in 0.615 seconds.
[Independent review artifact](../midpoint_angular_independent_review.json).

The next full-size calculation needs its own declared precision, finite
operator label, variable/row permutation, exact jump and rho2 checks, and
an early target-rank test. Root owns that execution. No full M50 worker,
regularizer, gauge-parameter substitution, or optimizer was introduced by
this construction.
