# Full analytic rank and the finite source-sized target plane

The declared M50/L10, degree-50 analytic amplitude completion now has
**augmented row rank 3011**. This is a theorem about its true analytic
operators, proved with enclosing intervals; the earlier theorem about a
rounded operator is no longer the only available rank result.

The [rank certificate](../analytic_full_source_rank3011_certificate.json)
gives

\[
\|I-R_{Q,0}B_Q\|_\infty
\le 2.206960154905338\times10^{-68}<1. \tag{AR1}
\]

Here \(B_Q\) is the selected 1510 by 1510 minor of the actual analytic
quotient after an exact constant pivot, with the recorded dyadic scales.
The rounded number displayed in (AR1) is an upper bound; the artifact
retains the exact rational bound and all 1510 rational row bounds.

Combined with the previously certified exact source-current sector, the
rank theorem proves that the **declared finite analytic model allows every
target pair** \((f^0_0(3),f^1_1(3))\). This applies even when all eight raw
chiral residuals and the additional infinity functional are zero, and when
the stated shared-node current, form-factor and moment constraints hold.
It does not assert this for the authors' unarchived implementation or for
continuum scattering unitarity.

## 1. The complete operator and numerical certificate

Let I denote the 1500 **raw imaginary f** rows of the analytic source
registry. Let C contain its 1500 real f rows, eight chiral rows, infinity
row, and two targets, in that order. Exact positive multiplication by
\(\kappa_j=\pi\sqrt{1-4/s_j}\) converts f to G; it preserves rank and
the imaginary kernel. The earlier imaginary rank theorem therefore gives
\(\operatorname{rank}I=1500\).

The [exact minimal quotient](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/65_minimal_imaginary_quotient.md`) gives
\(Q\in\mathbb R^{1511\times2376}\), with

\[
\operatorname{rank}\begin{bmatrix}I\\C\end{bmatrix}
=1500+\operatorname{rank}Q. \tag{AR2}
\]

Every entry of Q was enclosed in the
[complete calculation](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/68_complete_analytic_quotient.md`). Its first row has
exact constant component a=5/2. The exact row operations described below
remove this one pivot and leave \(H\in\mathbb R^{1510\times2375}\).
One column-pivoted QR chooses 1510 columns P; all 1510 rows remain. The
actual selected entries are assembled from the archived interval data,
not from the Float64 selection proxy:

\[
B_Q=D_rH_{:P}D_c.
\]

Both D matrices are exact positive dyadic diagonal matrices. The single
inverse candidate is computed at 1024 bits and frozen entrywise as an
exact dyadic matrix \(R_{Q,0}\). A separate 2048-bit full interval product
then evaluates (AR1). Every outward input conversion contains its original
ball. For any true analytic matrix represented by these enclosures,

\[
B_Q^{-1}=\sum_{k=0}^{\infty}(I-R_{Q,0}B_Q)^kR_{Q,0}
\]

converges, proving nonsingularity. This is independent of any numerical
rank cutoff or the success flag of the approximate solve.

| Stage | Verified result |
|---|---:|
| Imaginary analytic rank | 1500 |
| Complete quotient | 1511 by 2376, all entries enclosed |
| Exact pivot and selected minor | 1510 by 1510, all remaining rows retained |
| Candidate / verification precision | 1024 / 2048 bits |
| Strict residual bound | approximately 2.2069601549e-68 |
| Candidate infinity norm | approximately 5.31363221275e61 |
| Certified inverse norm of scaled minor | at most \(\|R_{Q,0}\|_\infty/(1-\delta_Q)\) |
| Successful supervisor time | 489.67 seconds |
| Cap / timeouts / retries / reselections | 900 seconds / none / none / none |
| Source or input changes during the run | none |

All 1510 rational row bounds and their strict maximum were rechecked after
termination, together with every direct input and candidate hash. The
independent [full-quotient data audit](../source_full_quotient_independent_audit.json)
checks all 3,590,136 exact encodings, every row metric, all 98 artifact
hashes and all constant anchors. Its comparisons with the independent
control run find all 4752 tested intervals exactly equal. The rational
bound checks and data audit do not claim to be a second computation of the
1510-order residual product.

The resulting rank chain is

\[
\boxed{\operatorname{rank}H=1510,\qquad
\operatorname{rank}Q=1511,\qquad
\operatorname{rank}M_{\rm aug}=3011.} \tag{AR3}
\]

The complete 3876-coefficient basis and every source row are retained.
The resulting dimensions are:

| Fixed linear data | Rank | Amplitude kernel dimension |
|---|---:|---:|
| All physical scattering samples, real and imaginary | 3000 | 876 |
| Those samples, eight chiral residuals, and infinity | 3009 | 867 |
| Those values and both targets | 3011 | 865 |

The target map has rank two on the 867-dimensional kernel in the middle
row. Thus it can vary both targets while all previously fixed data stay
exactly unchanged.

## 2. An explicit linear amplitude lift

The proof can specify a right inverse without finding a new 3011-column
minor in the original amplitude coordinates. Let A be the exact linear
map taking a coefficient **dual** r to its scaled quotient q(r), and put

\[
Z=A^T.
\]

Then \(IZ=0\), \(CZ=Q\), and Z is injective. It is an implicit basis map
onto the imaginary amplitude kernel. Z must not be confused with K, the
map that reconstructs the projected dual from its quotient in chapter 65.

Write \(Q=[c\mid V]\), with first row \([a\mid v_0]\), a=5/2. The
constant-pivot matrix is

\[
H_i=aV_i-c_i v_0\quad(i=1,\ldots,1510).
\]

For arbitrary \(t=(t_0,t_{\rm rest})\), set

\[
\begin{aligned}
y_{P}&=D_cB_Q^{-1}D_r(a t_{\rm rest}-c_{\rm rest}t_0),
&y_{\overline P}&=0,\\
z_0&=(t_0-v_0y)/a,
&R_Qt&=(z_0,y).
\end{aligned} \tag{AR4}
\]

Direct substitution gives \(QR_Q=I_{1511}\). All inverses in this
definition are the true inverses certified above, not their approximate
candidates.

Let R be the right inverse of the raw imaginary sample matrix I obtained
from [chapter 61](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/61_imaginary_right_inverse.md`) by the exact G/f row
scaling: if \(I_G=D_\kappa I\), then \(R=R_GD_\kappa\).
Given arbitrary desired imaginary samples g and remaining data h,
define the real amplitude coefficient vector

\[
\boxed{x=Rg+ZR_Q(h-CRg).} \tag{AR5}
\]

Indeed \(Ix=g\), while \(Cx=CRg+QR_Q(h-CRg)=h\). Equation (AR5) is an
explicit exact construction from the already defined small-factor maps
and the certified selected inverse. No dense right inverse, numerical
amplitude coefficient vector, or rational-coefficient claim is required.
The analytic coefficients in this model need not be rational.

## 3. Pure and chiral finite target projections

For every physical retained wave choose S=0. With
\(S=1+i\kappa_j f\), the desired f samples are

\[
\operatorname{Re}f=0,\qquad \operatorname{Im}f=1/\kappa_j.
\]

Choose all eight chiral residuals and the optional infinity functional
equal to zero, and let the two target values be any point of
\(\mathbb R^2\). Equation (AR5) realizes all these linear data
simultaneously. Every scattering disk is strict because \(|S|=0<1\).
Every chiral norm ball of nonnegative radius contains the zero residuals,
independently of whether the displayed families are grouped or separate.

Consequently the pure and chiral finite models in this completion have
the entire target plane. Changing only the norm or tolerance of these raw
chiral residuals cannot produce the bounded published regions while these
zero-residual lifts remain admissible.

## 4. The shared-node exact gauge model

The independent exact-source certificate
[source_M50_exact_current_fiber_256.json](../source_M50_exact_current_fiber_256.json)
establishes the necessary current-sector premise using the exact nodes,
phase-space factors, Jacobian and printed Hilbert-transform relation. It
keeps the saved rational spectral samples, targets and bounds fixed.
The earlier certificate explicitly did not assert an analytic amplitude
lift; (AR3)–(AR5) now supply that missing premise.

Use

\[
\theta_j=\frac{(2j+1)\pi}{200},\quad \phi_j=2\theta_j,\quad
s_j=4/\cos^2\theta_j,\quad F_j=(1+e^{i\phi_j})^3.
\]

The form factor is analytic with F(0)=1. In the S0 and P1 current waves,
choose \(S_j=\eta F_j/F_j^*=\eta e^{3i\phi_j}\), with \(\eta=9/10\).
All other 28 retained waves at each node have S=0. The saved positive
rational spectra satisfy the 100 strict current Gram inequalities, 14
high-energy FF inequalities and four raw scalar moment inequalities at
the exact source constants. Phase alignment reduces the current lower
density requirement to

\[
\rho_j>\frac{20}{19}k_j^2|F_j|^2,
\]

as derived in [chapter 48](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/48_shared_G_gauge_fiber.md`). The source-current
artifact certifies these inequalities and the exact affine Hilbert relation.

Convert each chosen S value to desired raw f samples:

\[
\operatorname{Re}f_j=\frac{\operatorname{Im}S_j}{\kappa_j},\qquad
\operatorname{Im}f_j=\frac{1-\operatorname{Re}S_j}{\kappa_j}. \tag{AR6}
\]

Set chiral and infinity values to zero and choose arbitrary targets d.
The exact lift (AR5) realizes (AR6). Its current S values are the **same
sampled amplitude values** used for scattering, so the fixed current
Grams and FF/moment inequalities continue to hold as d varies. The disk
margin is at least \(1-\eta^2=19/100\).

This proves the entire target plane for this explicitly defined analytic
gauge variant, with the source-current artifact's 0.140 GeV pion mass,
1.2 GeV matching energy, 43-node hard mask, printed normalized moment
targets converted to raw units, raw tolerance .002, and declared scalar
quark-mass completion. The proof does not silently change those inputs.

## 5. Implications and remaining reproduction work

The missing analytic-rank step is now resolved for this completion. No
rank-transfer inference from a nearby rounded matrix is used. The very
large inverse norm also explains why ordinary floating-point numerical
behavior alone was not a reliable guide to the finite feasible set.

A further support scan of this same unrestricted finite model cannot
produce a bounded target region: its projection is provably all of
\(\mathbb R^2\). The published bounded plots therefore require a different
operator completion, additional restrictions, or an effective numerical
restriction not established by these declared inputs alone. The result
does not identify which of those possibilities applies to the authors'
unarchived calculation.

The next reproduction work must resolve or bracket that distinction and
then verify the relevant amplitudes and phase curves. Finite samples give
no scattering guarantee between energies or at omitted spins. No continuum
unitarity, published pure/chiral/gauge region, rho phase crossing, or
support optimality is proved here. The original end-to-end goal remains
active.
