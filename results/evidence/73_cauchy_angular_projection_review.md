# Angular projection of the original-center midpoint Cauchy kernels

The raw Cauchy-kernel projections, their subtraction terms, and the
isospin/density assembly have been independently checked. The resulting
formulas agree with the angular part of the authors' later notebook after
restoring `nu0=0`, and with the new
[midpoint/PV implementation](../runs/repo_reorganization_20260906/recovery.json) (historical member: `src/smatrix_bootstrap/pv_midpoint_kernels.py`).
Seventeen small validated tests pass, and a separate M2 comparison gives
868 passing interval checks. No M50 matrix, rank calculation, or optimizer
was run for this review.

## 1. Keep dispersion and angular normalization separate

The original [2023 TeX](../../references/2309.12402v3-source/prd_submission_2.tex)
specifies `1/pi` for each dispersion integration, `1/pi^2` for the double
integration, and `f_l^I=(1/4) integral P_l T^I`. With the original map center,

\[
 x(\phi)=4\sec^2(\phi/2),\qquad
 x'(\phi)=\frac{8\sin\phi}{(1+\cos\phi)^2}
          =x(\phi)\tan(\phi/2).
\]

Midpoint angles are `phi_j=(j+1/2)pi/M`. Thus one integrated Cauchy kernel is

\[
 h_j(\nu)=\frac{w_j}{x_j-\nu},\qquad
 w_j=\frac{\Delta\phi}{\pi}x'_j=\frac{x'_j}{M}.
\]

The `1/pi` has already canceled the midpoint factor pi. The independent
angular factor remains `1/4`. There is no additional `1/pi` to attach to
h or its products, and no physical `pi beta` factor in an f-normalized
row.

## 2. Exact single and crossed-double projections

Let `d=s-4`, `epsilon=(-1)^ell`, and

\[
 t=-\frac d2(1-\mu),\quad u=-\frac d2(1+\mu),\qquad
 z_i=1+\frac{2x_i}{d},\qquad c_\ell=\delta_{\ell0}/2.
\]

Use the integral-defined exterior Legendre function

\[
 Q_\ell(z)=\frac12\int_{-1}^1\frac{P_\ell(\mu)}{z-\mu}\,d\mu.
\]

For d nonzero,

\[
 \boxed{a_i^\ell(s):=\frac14\int_{-1}^1P_\ell(\mu)h_i(t)\,d\mu
       =\frac{w_i}{d}Q_\ell(z_i),}\qquad
 \frac14\int P_\ell h_i(u)=\epsilon a_i^\ell.
\]

The factor follows from `x_i-t=(d/2)(z_i-mu)`; the factor two in the
definition of Q then cancels half of the angular quarter.

For the double projection, the two denominators add to the positive
constant `D_ij=x_i+x_j+d`. Partial fractions give

\[
 \frac1{(x_i-t)(x_j-u)}=
 \frac1{D_{ij}}\left(\frac1{x_i-t}+\frac1{x_j-u}\right).
\]

Consequently

\[
 \boxed{u_{ij}^{\ell,\mathrm{raw}}
 =\frac14\int P_\ell h_i(t)h_j(u)
 =\frac{w_j a_i^\ell+\epsilon w_i a_j^\ell}{x_i+x_j+s-4}.}
\]

This implies `u_ji=epsilon u_ij`. The diagonal is regular:

\[
 u_{ii}^{\ell,\mathrm{raw}}=
 \frac{(1+\epsilon)w_i a_i^\ell}{2x_i+s-4}.
\]

In particular every odd-spin diagonal vanishes exactly. No division by
`x_i-x_j`, derivative limit, or special diagonal factor one half is needed
for the crossed product. For every `s>0` and `x_i,x_j>=4`, its denominator
is strictly greater than four.

Any direct-channel factor independent of mu can be pulled outside the
projection: `integral(P_l h_i(s))/4=c_l h_i(s)` and
`integral(P_l h_i(s)h_j(t))/4=h_i(s)a_j`. This statement does not assign a
finite value to the pole at `s=x_i`; the physical PV prescription is a
separate operation.

## 3. Subthreshold branch and threshold limits

For physical s>4, `z_i>1`. For `0<s<4`, one has `z_i<-1`, because
`x_i>4-s`. All crossed denominators stay positive on the real angular
interval. The real continuation required by the integral is

\[
 Q_\ell(-Z)=(-1)^{\ell+1}Q_\ell(Z),\qquad Z>1.
\]

It follows directly by reflecting the integration variable. For example,

\[
 Q_0(z)=\tfrac12\log\frac{z+1}{z-1},\qquad Q_1(z)=zQ_0(z)-1,
\]

where the displayed logarithm has a positive real argument on both
exterior intervals. No Ferrers continuation adding an imaginary multiple
of `pi P_l` is appropriate below threshold. The notebook explicitly uses
`LegendreQ[ell,0,3,z]`, consistent with this exterior branch and the
integral definition in [the later paper, Eq.(8.7)](https://arxiv.org/html/2403.10772v2#S8.SS1).

The apparent d=0 singularity is removable. In particular,

\[
 a_i^0(s)=\frac{w_i}{2d}\log(1+d/x_i),\qquad
 a_i^\ell(4)=c_\ell w_i/x_i.
\]

The double limit is `c_l w_i w_j/(x_i x_j)`. After subtraction at zero,
the t and u kernels vanish at threshold, and so do their projected
products. These limits are implemented as exact identities when s=4 is
specified exactly.

## 4. Stable analytic evaluation with an explicit tail

Forward recurrence for Q at `z>>1` can cancel hundreds of decimal digits.
A direct moment expansion avoids that issue. Set

\[
 A_i=x_i+d/2>0,\qquad r_i=\frac{d}{2x_i+d},\qquad |r_i|<1.
\]

Expanding `1/(A_i-d mu/2)` and using Rodrigues' formula gives

\[
 \frac14\int_{-1}^1P_\ell(\mu)\mu^{\ell+2m}\,d\mu
 =\frac{(\ell+2m)!}{2(2m)!!(2\ell+2m+1)!!}.
\]

All lower-degree or opposite-parity powers integrate to zero. Hence

\[
 a_i^\ell=\frac{w_i}{2A_i}r_i^\ell
 \sum_{m=0}^{\infty}c_{\ell m}r_i^{2m},\qquad
 c_{\ell m}=\frac{(\ell+2m)!}{(2m)!!(2\ell+2m+1)!!}>0.
\]

Here `c_0m=1/(2m+1)` and
`c_(ell+1),m/c_ell,m=(ell+2m+1)/(2ell+2m+3)<1`, so `c_ell,m<=1`.
Truncating before m=K therefore has the rigorous absolute tail

\[
 \boxed{|\mathrm{tail}|\le
 \frac{|w_i|}{2A_i}\frac{|r_i|^{\ell+2K}}{1-r_i^2}.}
\]

The standalone [reference evaluator](../runs/repo_reorganization_20260906/recovery.json) (historical member: `cauchy_angular_review.py`) uses this
positive-coefficient series for `|r|<=1/2`, with the tail added outward.
For larger |r| it uses the real logarithm and recurrence only on
`1<|z|<2`, with interval arithmetic. The new production prototype instead
uses the exterior hypergeometric identity

\[
 Q_\ell(z)=\frac{2^\ell(\ell!)^2}{(2\ell+1)!}\,z^{-\ell-1}
 {}_2F_1\left(\frac{\ell+1}{2},\frac{\ell+2}{2};
 \ell+\frac32;z^{-2}\right)
\]

on positive z and the proved parity continuation for negative z. The
different evaluation paths provide useful independent controls. Their
ell=0 case reduces to `atanh(1/z)`, fixing the normalization.

## 5. What W, W0, and W00 mean in the notebook

The [archived notebook](../../references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb)
uses center `nu0=-20`; the present derivation restores zero. Its subtracted
kernel is `H_i(nu)=h_i(nu)-b_i`, with `b_i=w_i/x_i`. At zero center its
arrays implement

\[
 W0_i=b_i,\quad W_{i|j}(s)=\frac{w_i}{x_i+x_j+s-4}-b_i,
 \quad W00_{ij}=b_i b_j.
\]

These are subtraction components, not independent angular factors. The
complete projected kernels are

\[
 j_i=a_i-c_\ell b_i,\qquad
 u_{ij}=a_iW_{j|i}+\epsilon a_jW_{i|j}+c_\ell b_i b_j.
\]

Expanding `(h_i(t)-b_i)(h_j(u)-b_j)` proves this formula immediately.
The `W00` term carries `c_l=delta_l0/2`; omitting that factor is incorrect.
The index reversal in the first W factor is essential.

Notebook `Ap/Am` and `Ap1/Am1` are `Omega(s)a_i` and
`Omega(s)epsilon a_i` in the physical region, where the actual code sets
`Omega=pi sqrt((s-4)/s)`. The unphysical versions omit Omega and evaluate
real f rows. This agrees with the source relation `S=1+i pi beta f`.
Equations (8.3)–(8.9) of the later paper describe these subtractions; the
2023 representation is unsubtracted. A variable redefinition between them
must be explicit when comparing complete finite operators.

## 6. Independent full-row and symmetric-density review

Write d_i for the chosen direct-channel kernel, j_i for the crossed single,
and u_ij for the matching crossed double. Use either all raw or all
consistently subtracted kernels. Let a and b denote projections of one
basis contribution to `A(s,t,u)` and `A(t,s,u)`:

| Density variable | a | b |
|---|---|---|
| T0 | c_l | c_l |
| sigma1_i | c_l d_i | j_i |
| sigma2_i | (1+epsilon)j_i | c_l d_i+epsilon j_i |
| rho1_ij | (1+epsilon)d_i j_j | d_j j_i+u_ij |
| rho2_ii | u_ii | epsilon d_i j_i |
| rho2_ij, i<j | (1+epsilon)u_ij | epsilon(d_i j_j+d_j j_i) |

The I0/I1/I2 combinations are respectively `3a+(1+epsilon)b`,
`(1-epsilon)b`, and `(1+epsilon)b`. Forbidden parities vanish. This
matches the original definitions `T1=B-C` and `T2=B+C`. The plus/minus
rows in the later paper's displayed Eq.(8.12)/(8.13) should not be read
using their apparently swapped I1/I2 labels; the notebook stacks the
corresponding blocks in I0,I2,I1 order.

For a genuine symmetric sampled density, rho2_ij is the common value of
the two full-array entries. Its off-diagonal basis is therefore
`h_i(t)h_j(u)+h_j(t)h_i(u)`, with **no half**. The notebook implements
`2 symtr`, then halves only the diagonal before extracting the upper
triangle. That gives precisely the table above. A half-sum basis would
require a compensating coefficient redefinition, including any density
norm; it is not the same sampled-density convention.

The actual code in `pv_midpoint_kernels.py` was reviewed against all these
formulas. No angular, branch, transposition, threshold, or rho2-packing
blocker was found. Assigning physical `PV+iI` values remains a discrete
prescription: they are not literal finite limits of the off-cut rational
kernels at their source-node poles. This angular review does not replace
that separate consistency analysis.

## 7. Bounded tests and provenance

The [17 tests](../runs/repo_reorganization_20260906/recovery.json) (historical member: `tests/test_cauchy_angular_review.py`) compare the formulas
with direct Arb Gauss-Legendre sums and an explicit Bernstein-ellipse
remainder. If the chosen ellipse excludes the rational poles, each
denominator has a positive real-part gap. Bounding the rational kernel
by B and using `|P_l|<=rho^ell` gives normalized quadrature error
`2 B rho^(ell-2q)/(1-1/rho)`. The tests include both energy regions,
crossed diagonals, subtractions, very large Q arguments, and threshold
limits; their agreement is not a Float64 tolerance check.

The separate [review artifact](../midpoint_angular_independent_review.json)
records 868 passing interval comparisons of A/J/U/W/W00 at x=(5,13),
w=(2,3), all five source subthreshold energies and s=9, with ell=0,...,5.
It also checks the exact M2 grid and Jacobian using square-root-two formulas.
That comparison took 0.615 seconds at 256 bits and snapshots both
implementations and the primary-source hashes. The checked kernel-module
SHA256 is `926493cab1b7eaf80c1fa50f3dc379e361bc33b89b02d293b328ea9885c18512`.

These results establish and test the angular ingredients. They do not
establish a full source-sized matrix rank, target kernel, regularizer
choice, or reproduction of the published bounded region.
