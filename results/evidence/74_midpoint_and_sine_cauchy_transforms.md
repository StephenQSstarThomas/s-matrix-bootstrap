# Exact midpoint and finite-sine Cauchy transforms

The literal midpoint Cauchy sum and exact finite-sine transform use the same
M density samples but define different off-cut functions. Their exact
correction is generally full rank across evaluation points. The midpoint
sum has real source-node poles; assigning a finite nodal PV matrix and an
imaginary jump is a separate prescription.

This addresses part of the [source decision](71_source_boundedness_decisions.md).
It does not transfer the finite-sine rank-3011 theorem to a PV/midpoint map.
No large rank calculation, optimizer, or frozen-source edit was performed.

## 1. Normalization and definitions

The original conformal-map center is nu0=0:

\[
 z(\nu)=\frac{2-\sqrt{4-\nu}}{2+\sqrt{4-\nu}},\quad
 \nu=\frac{16z}{(1+z)^2},\quad
 \phi_j=\frac{\pi(j+1/2)}M,\quad x_j=4\sec^2(\phi_j/2).
\]

Put \(x'_j=x_j\tan(\phi_j/2)\) and \(w_j=x'_j/M\). The weight includes
the Jacobian: the dispersion factor 1/pi times midpoint spacing pi/M is
1/M. These are samples of a density with respect to dx, not dphi.
For real samples rho_j, define the **unsubtracted** transforms

\[
 H_{\rm mid}(\nu)=\sum_j\frac{w_j\rho_j}{x_j-\nu},\qquad
 H_{\rm sine}(\nu)=\sum_{n=1}^M a_n\{z(\nu)^n-(-1)^n\},
 \quad a=T^{-1}\rho,\quad T_{jn}=\sin(n\phi_j). \tag{MC1}
\]

The second is the exact continuous dispersion integral of the sine
interpolant, by [chapter 3](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/03_discretization.md`). All M modes are kept:
\(T^TT=\operatorname{diag}(M/2,\ldots,M/2,M)\). A charge subtraction is
H(nu)-H(0); it is distinct from choosing the map center nu0=0.

The square root is positive below threshold. On the upper physical lip,
\(z=e^{i\phi}\) and \(\sqrt{4-s}=-i\sqrt{s-4}\); lower values are conjugate.

## 2. Exact midpoint response of every sine mode

For rho_j=sin(n phi_j), 1<=n<=M,

\[
 \boxed{H_{{\rm mid},n}(\nu)
 =\frac{z^n+z^{2M-n}}{1+z^{2M}}-(-1)^n
 =\frac{T_{M-n}(8/\nu-1)}{T_M(8/\nu-1)}-(-1)^n.} \tag{MC2}
\]

The last T symbols are Chebyshev polynomials. The singularity at nu=0
is removable, and the expression is rational in nu.

To prove (MC2), compare poles and residues. The source cosines c_j are
the roots of T_M. With sigma_j=(-1)^j,

\[
 T'_M(c_j)=\frac{M\sigma_j}{\sin\phi_j},\qquad
 T_{M-n}(c_j)=\sigma_j\sin(n\phi_j).
\]

Since \((8/\nu-1)'=-8/\nu^2\), the coefficient of \(1/(x_j-\nu)\)
in (MC2) is \(x_j^2\sin\phi_j\sin(n\phi_j)/(8M)=w_j\sin(n\phi_j)\).
Both rational functions vanish at infinity, so matching all residues
proves equality. Zero residues simply cancel the corresponding pole.

The Nyquist mode is

\[
 H_{{\rm mid},M}=\frac{2z^M}{1+z^{2M}}-(-1)^M,\qquad
 H_{{\rm sine},M}=z^M-(-1)^M. \tag{MC3}
\]

The factor two is essential. At the map center both transforms obey
\(H_n(0)=-(-1)^n\), hence

\[
 H(0)=\ell^T\rho,\qquad \ell_j=\tan(\phi_j/2)/M.
\]

For \(P(z)=\sum a_nz^n\), the charge-subtracted midpoint transform is
\([P(z)+z^{2M}P(1/z)]/(1+z^{2M})\), whereas the subtracted sine transform
is P(z).

## 3. The exact correction is generally full rank

Mode by mode,

\[
 H_{{\rm mid},n}-H_{{\rm sine},n}
 =\frac{z^{2M-n}(1-z^{2n})}{1+z^{2M}}. \tag{MC4}
\]

In nodal coordinates the cardinal identity of
[chapter 51](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/51_analytic_imaginary_structure.md`) gives

\[
 \boxed{H_{\rm mid}(\nu)-H_{\rm sine}(\nu)
 =\frac{z(\nu)^M\sqrt{4-\nu}}{2M}
       \sum_j\frac{\sigma_jx_j\rho_j}{x_j-\nu}.} \tag{MC5}
\]

This difference is unchanged by subtracting the common value at zero.
Across K distinct nonzero off-cut points, its matrix is a Cauchy matrix
\(1/(x_j-\nu_i)\) multiplied by nonzero row and column factors. Its rank
is min(K,M), so it is not generally rank one or low rank. Evaluation at
nu=0 gives an exact zero row and is an explicit exception.

For M points, the Cauchy determinant has numerator
\(\prod_{i<k}(\nu_i-\nu_k)\prod_{j<l}(x_l-x_j)\) and denominator
\(\prod_{i,j}(x_j-\nu_i)\), which proves nonsingularity under these conditions.

There is a low-rank **leading Taylor term**: the correction to mode n
starts at \(z^{2M-n}\), so only the Nyquist coefficient contributes at
order \(z^M\). This is not the rank of the complete correction. On |z|<=r<1,

\[
 |H_{{\rm mid},n}-H_{{\rm sine},n}|
 \le\frac{r^{2M-n}(1+r^{2n})}{1-r^{2M}}.
\]

Since z has a simple zero at nu=0, every nodal density has identical
midpoint and sine derivatives through order M-1 at that center. The
first possibly different derivative, order M, depends only on the
Nyquist coefficient. Their continuous values also agree at nu=4;
this does not equate their threshold derivatives or physical jumps.

A general-density bound must include its mode coefficients. A small
off-cut operator error does not uniformly control unbounded coefficient
families or preserve a sampled target kernel; see
[chapter 50](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/50_rank_transfer_limits.md`).

Both transforms retain all M density directions as global off-cut
functions. Midpoint residues recover each rho_j; sine polynomial
coefficients recover each a_n. Their unsubtracted spans intersect only
at zero: a polynomial in z that is also rational in nu must be invariant
under z->1/z, hence constant, and the unsubtracted transform vanishes at
infinity. With a free constant, the intersection is precisely constants.
Equal dimensions therefore do not imply an invertible reparameterization
of the same function family.

## 4. Physical boundaries, PV, and the pole's finite part

Away from source nodes, the midpoint boundary is

\[
 H_{{\rm mid},n}(s\pm i0)
 =\frac{\cos((M-n)\phi)}{\cos(M\phi)}-(-1)^n\in\mathbb R. \tag{MC6}
\]

The two lips agree there. The sine boundaries instead have imaginary
parts plus/minus sin(n phi), and

\[
 H_{{\rm mid},n}-H_{{\rm sine},n}^{\pm}
 =\sin(n\phi)\{\tan(M\phi)\mp i\}. \tag{MC7}
\]

At x_i, H_mid has the pole \(w_i\rho_i/(x_i-\nu)\) unless rho_i=0.
Its distributional jump is \(2\pi i\sum_jw_j\rho_j\delta(s-x_j)\);
it is not the continuous sine-density jump. An ordinary boundary value
at a nonzero source density is therefore not a finite PV term plus/minus
i rho_i.

The exact finite Hilbert matrix

\[
 K_{ij}=\frac2M\sum_{n=1}^{M-1}\cos(n\phi_i)\sin(n\phi_j)
\]

sends sin(n phi) to cos(n phi) for n<M and annihilates sin(M phi).
It is the real nodal map for the **subtracted** sine transform. The
unsubtracted sine boundary map is

\[
 (K+\mathbf1\ell^T)\rho\ \pm i\rho. \tag{MC8}
\]

The highest density mode cannot be deleted merely because its real K
image is zero: its imaginary values and off-cut transform are nonzero.

The midpoint pole's finite part is also distinct. Set
\(A_{{\rm off},ij}=w_j/(x_j-x_i)\) for j!=i and set its diagonal to zero.
Let

\[
 D_\phi=\dot T T^{-1},\quad \dot T_{jn}=n\cos(n\phi_j),\qquad
 d_i=\frac{x_i''}{2Mx_i'}=\ell_i+\frac{\cot\phi_i}{2M}.
\]

The last column of dot-T is zero. The exact finite-part identity is

\[
 A_{\rm off}=K+\mathbf1\ell^T-\frac1M D_\phi-\operatorname{diag}(d).
 \tag{MC9}
\]

For a sine mode, expansion of (MC6) at phi_i gives the phi-finite part
\((1-n/M)\cos(n\phi_i)-(-1)^n\). The self-pole itself has phi-finite
part \(\rho_i x_i''/(2Mx_i')\), which must also be removed to obtain
the finite part in nu. This proves (MC9), including its diagonal term.
For charge-subtracted transforms, remove the common rank-one constant
matrix from both sides. Neither finite part is simply K.

A PV matrix and midpoint off-cut sums can be numerical approximations to
a continuum dispersion problem. They are not the ordinary boundary and
off-cut values of the same finite rational function (MC1). Another shared
analytic interpretation requires its own construction and proof.

## 5. Analytically projected midpoint kernels

For \(t=-(s-4)(1-\mu)/2\), \(u=-(s-4)(1+\mu)/2\), put

\[
 C_\ell(s,x)=\frac14\int_{-1}^1\frac{P_\ell(\mu)}{x-t}\,d\mu.
\]

The real Legendre-Q integral gives

\[
 C_\ell(s,x)=
 \begin{cases}
 (s-4)^{-1}Q_\ell(1+2x/(s-4)),&s>4,\\
 (-1)^\ell(4-s)^{-1}Q_\ell(2x/(4-s)-1),&0<s<4.
 \end{cases} \tag{MC10}
\]

All Q arguments displayed here exceed one. The single-cardinal projection
is w_j C_l(s,x_j). At s=4 its continuous extension is
C_0(4,x)=1/(2x) and C_l(4,x)=0 for l>0. Exact partial fractions give

\[
 \frac14\int_{-1}^1P_\ell h_i^{\rm mid}(t)h_j^{\rm mid}(u)\,d\mu
 =\frac{w_iw_j}{x_i+x_j+s-4}
   \{C_\ell(s,x_i)+(-1)^\ell C_\ell(s,x_j)\}. \tag{MC11}
\]

These formulas retain the projection factor one quarter. They multiply
the full density matrix; symmetric upper-triangle packing weights must
be applied separately and consistently.

The projected single-cardinal correction is

\[
 \frac14\int_{-1}^1P_\ell
       (h_j^{\rm mid}(t)-h_j^{\rm sine}(t))\,d\mu
 =\frac{\sigma_jx_j}{8M}\int_{-1}^1
       \frac{P_\ell(\mu)\sqrt{4-t}\,z(t)^M}{x_j-t}\,d\mu. \tag{MC12}
\]

This is not a universally low-rank angular operator. For example, at
fixed s>4, consecutive spins 0,...,L-1 and any L distinct density nodes
(L<=M) give a nonzero determinant: after fixed signs are removed the
weight sqrt(4-t)|z(t)|^M is positive, and the moment-determinant integral
contains a Legendre Vandermonde and a Cauchy determinant of fixed sign.
Parity-restricted and fully crossed source families need their own rank
analysis.

For one fixed argument pair the double-density coefficient correction
\(h_{\rm mid}(\nu)h_{\rm mid}(\omega)^T-
h_{\rm sine}(\nu)h_{\rm sine}(\omega)^T\) has rank at most two. Collecting
many pairs or integrating angularly need not preserve that small rank.

## 6. Rank scope and small exact checks

The literal off-cut midpoint amplitude also retains its full density
dimension as a rational function family. At s=x_i,t=x_j the unrestricted
s-t double pole isolates that coefficient; u=4-x_i-x_j<0 is not another
source pole. The t-u double poles isolate the symmetric family. Once
those residues vanish, single poles isolate the two single families and
then the free constant. This is not a rank theorem for a hybrid nodal
PV/midpoint operator or a claim about its physical boundary interpretation.

The [small exact check](../midpoint_sine_transform_small_exact.json)
and its [script](../runs/repo_reorganization_20260906/recovery.json) (historical member: `scripts/check_midpoint_sine_transform.py`) verify 21
mode identities for M=1,...,6, including every Nyquist mode. This uses
Chebyshev-root traces independently of (MC2). Fourteen small-M finite-part
checks verify (MC9), and a two-node off-cut correction has rank two.
The run took about 1.27 seconds and used no quadrature or optimization.

No full M50 PV/midpoint rank, target image, bounded contour, or coefficient
restriction is inferred. The previous analytic rank and gauge-plane
corollary remain statements about the declared finite-sine completion.
A newly specified physical-PV/off-cut-midpoint row map needs a fresh
kernel test while retaining every original variable.
