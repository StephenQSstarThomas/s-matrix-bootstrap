# Source-sized imaginary rank for the discrete PV/midpoint operator

The declared M50/L10 PV/midpoint operator now has **imaginary row rank
1500**. A new interval certificate proves rank 450 of its own residual
map. The finite-sine rank certificate was not transferred to this operator.
The real-row rank, target-row rank and interpretation as one global analytic
amplitude remain separate questions.

The [certificate index](../pv_midpoint_imaginary_rank1500_certificate.json)
records the exact rational residual bound, all artifact links and scope.
One 300-second attempt finished successfully in 30.42 seconds, with

\[
\|I-R_0B\|_\infty
\le 3.150777397246630\times10^{-262}<1.
\]

The displayed number rounds upward; the full report retains the exact
rational bound and all 450 rational row-sum bounds. No optimizer was run.

## 1. The operator and coefficient pairing

This is the original-center discrete prescription of
[chapter 72](72_midpoint_PV_kernel_construction.md): exact M50 nodes,
w=x'/M, charge-subtracted off-cut midpoint kernels, and assigned physical
values K+iI. Partial waves use the factor one quarter and are f-normalized.
Physical G rows differ by the exact positive factor pi beta.

The unweighted amplitude density variables are T0, two independent
50-vectors, one unrestricted 50-square rho1, and the upper triangle of
symmetric rho2: 3,876 real coordinates. In the stored physical density
convention, rho2 off-diagonal entries multiply a sum of two products.
To use the generic flat-C algebra, define

\[
C_{\rm flat,ii}=\rho_{2,ii},\qquad
C_{\rm flat,ij}=2\rho_{2,ij}\quad(i<j).
\]

The corresponding physical-row conversion divides only those off-diagonal
columns by two. This is an invertible column change; it preserves every
left-null relation and row rank. It is not a density restriction.

The physical imaginary jump is the nodal unit vector u=e_j. Thus the
phase-coordinate matrix in the generic algebra is **T=I50**, not the
trigonometric sine matrix used by the previous analytic completion.
The off-cut crossed moments J are real. In the full symmetric matrix
C=rho2 convention, the imaginary factorization is

\[
\begin{aligned}
\Im f^0_\ell&=c_\ell(3\sigma_1+2\sigma_2)^Tu
 +2J^T(3\rho_1^T+\rho_1+C)u,\\
\Im f^2_\ell&=2c_\ell\sigma_2^Tu+2J^T(\rho_1+C)u,\\
\Im f^1_\ell&=2J^T(\rho_1-C)u,\qquad
c_\ell=\tfrac12\delta_{\ell0}.
\end{aligned}
\]

The first two use even spins, the last odd spins. The real PV matrix and
real double-crossed moments drop out of this imaginary factorization.
Their absence here does not remove them from the full discrete operator.

## 2. Why the same residual algebra applies

Use D=Im(f0-f2), E=Im f2 and V=Im f1. In a left-null functional, the
independent single-density coefficient equations force the 50 D S-wave
weights and 50 E S-wave weights to zero.

At each node, the 50 by 19 positive-spin J matrix has columns in the order
even2,...,18 followed by odd1,...,19. The new local inverse certificates
prove these columns independent. Let E_j and O_j denote the two column
blocks. They have dimensions 9 and 10 with zero intersection.

The remaining coefficient equations are exactly

\[
3P^T+Q+R=0,\qquad \operatorname{Sym}(Q-R)=0,
\]

with P[:,j],Q[:,j] in the image of E_j and R[:,j] in the image of O_j.
Writing P[:,j]=E_j alpha_j gives 450 unknowns. For each node, a recorded
19-row selector S_j defines the true inverse

\[
A_j=(S_jJ_j)^{-1}S_j,\qquad
M_j=I-J_jA_j,\qquad \Delta_j=[E_j\ {-O_j}]A_j.
\]

Membership and symmetric residual rows are

\[
\{M_jP^T_{:j}\}_j,\qquad
-3\operatorname{pack}\operatorname{Sym}
  (\{\Delta_jP^T_{:j}\}_j).
\]

They form F of shape 3775 by 450. Packing is the ordinary unweighted
upper triangle, with the usual one-half in off-diagonal Sym entries.
There are 950 exact membership identities, determined by the **actual
arbitrary selectors**, leaving 2,825 active rows. No assumption that the
first 19 density rows were selected is made.

Local independence makes this residual kernel equivalent to the original
imaginary left-null space. Consequently

\[
\operatorname{rank}(\Im f_{\rm PV/midpoint})=1050+\operatorname{rank}F.
\]

This uses the generic algebra of
[chapters 51](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/51_analytic_imaginary_structure.md`) and
[57](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/57_imaginary_residual_rank_review.md`), applied to the newly certified
nodal-jump data. It does not reuse any full rank of the sine-based F.

## 3. Authenticated inputs and preflight

The source-factor archive is
[indexed here](../pv_midpoint_source_factors_2048.json), with
[factor report](../pv_midpoint_storage/source_factors_2048_attempt01/report.json)
SHA256:

    9cc07c41216785d7ac85b435ab3fa7dc264b661b9f343a2bfd2f62c05402ac8b

It contains the new 2048-bit moments and validated local inverse factors.
The producer keeps the true J spin columns unscaled; its preliminary
column scaling is only a Float64 selector device.

The separate
[input preflight](../pv_midpoint_storage/imaginary_rank_preflight_2048_attempt01/report.json)
completed in 2.74 seconds. It read 103 direct inputs: the report and
source snapshot, grid, 50 master files and 50 local-factor files. It checked:

- Matching PV/midpoint, nu0=0, f, subtraction and flat-C conventions.
- All 2,500 entries of the exactly identity nodal-jump matrix.
- All 50 distinct source nodes, complete column and coefficient labels.
- All 50 recorded local rational residual bounds below one, and the
  infinity norms of their saved 19 by 19 residual matrices.
- 47,500 E/O-to-master interval overlaps and all 950 selected exact-zero
  membership rows, using each node's own selector.

The preflight generated no F rows, performed no global QR or inverse,
and recomputed no local residual matrix product. The validated producer
definitions and authenticated archives supply the local inverse enclosures.

Seven focused tests passed, including 396 exact synthetic comparisons of
the PV density rows with the generic flat-C imaginary formulas, arbitrary
selector coverage, wrong-convention/hash rejection and timeout flags.
The tests do not use an M50 physical rank computation.

## 4. One complete interval rank gate

The new [driver](../runs/repo_reorganization_20260906/recovery.json) (historical member: `scripts/certify_pv_midpoint_imaginary_rank.py`) evaluates
all active residual rows from the authenticated interval factors. A
Float64 midpoint proxy receives exactly two dyadic equilibration passes.
One QR of its transpose selects 450 distinct rows, retaining all 450
unknowns. There is no numerical rank threshold or discarded coordinate.

The true selected minor B is reassembled at 2048 bits from the original
interval factors with the recorded row and column powers. It is saved
before the single 1024-bit approximate inverse solve. The returned finite
midpoint matrix R0 is fixed as exact dyadics and saved before verification.
A separate full 2048-bit interval product computes I-R0 B.

An approximate inverse alone has no certification role. The outward
infinity norm is the maximum row sum of absolute interval bounds, and its
strict rational maximum must be below one. All 450 rational row bounds
are retained. No large residual-matrix blob is needed because B, R0 and
the row bounds are archived.

| Quantity | Result |
|---|---:|
| Full residual / active residual shape | 3775x450 / 2825x450 |
| Selected minor | 450x450 |
| Proxy QRs / inverse candidates / verification products | 1 / 1 / 1 |
| Candidate / verification precision | 1024 / 2048 bits |
| Residual infinity upper bound | approximately 3.1507773972466293e-262 |
| Candidate infinity norm | approximately 3.5837364741175884e38 |
| Supervisor time / cap | 30.41946s / 300s |
| Retries, reselections, source/input changes | none |

The
[complete rank report](../pv_midpoint_storage/source_imaginary_rank_1024_2048_attempt01/report.json)
has SHA256:

    68f2ec9e36e6bb849e949a50bfa735a3cee38ccdb24a958a4b6d3ea0650f2dff

The [selection](../pv_midpoint_storage/source_imaginary_rank_1024_2048_attempt01/selection.json),
[interval B](../pv_midpoint_storage/source_imaginary_rank_1024_2048_attempt01/selected_minor_B.json.gz)
and [exact R0](../pv_midpoint_storage/source_imaginary_rank_1024_2048_attempt01/R0_exact_dyadic.json.gz)
are immutable outputs.

After termination, all 103 direct-input and five output-artifact hashes
were rechecked. The exact maximum of all 450 rational row bounds equals
the reported delta and is strictly below one. The stated inverse-norm
bound also satisfies ||B^-1||<=||R0||/(1-delta). These scalar and hash
checks are not a second residual-matrix multiplication.

## 5. What is established and what remains separate

The selected minor is nonsingular, so rank(F)=450 and the new discrete
imaginary sample matrix has rank 1500. In particular all retained imaginary
sample values can be prescribed by some coefficient vector. This alone
does not fix their real partners or show that an arbitrary unitary set of
complex samples, chiral values, or target pair can be realized.

The M4, M5 and M8 direct-angular controls were useful checks but supplied
no source-M50 rank transfer. Conversely, the new source-M50 imaginary result
is not an augmented-rank or target-plane theorem. The fully symmetric
capacity obstruction in the smaller M4 case and the independently proved
finite-sine results retain their own scopes.

The physical-node PV+jump and off-cut midpoint assignments are the
declared discrete operator of chapter 72. This certificate does not
establish them as one global analytic amplitude, infer continuum
unitarity, or identify the authors' unarchived 2023 numerical operator.
No density regularizer, coefficient bound, phase fit or optimizer was
introduced by this rank gate.
