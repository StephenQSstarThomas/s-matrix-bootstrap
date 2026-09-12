# Isospin capacity and complete small PV/midpoint rank controls

The alternative PV/midpoint prescription now passes independent complete
small-matrix checks. A first M4/L1 attempt also exposed an exact isospin
capacity obstruction: its 35 rows cannot all be independent, despite
having 35 amplitude coefficients. The obstruction is algebraic and is
not specific to midpoint integration or an arithmetic precision.

The operator conventions are those of
[chapter72](72_midpoint_PV_kernel_construction.md): original-center midpoint
nodes, assigned subtracted PV+iI physical values, analytically integrated
crossed Cauchy kernels, and the full unrestricted density families.
The mathematical distinction from the finite-sine analytic completion
is established in [chapter74](74_midpoint_and_sine_cauchy_transforms.md).

## 1. A representation-sensitive dimension bound

For M density coordinates the amplitude has

\[
d=1+2M+M^2+M(M+1)/2
\]

real coefficients. The following subspace gives a fully permutation-symmetric
scalar amplitude, for any consistently used one-variable kernels:

\[
\sigma_1=\sigma_2,\qquad \rho_1=\rho_1^T=\rho_2.
\]

Its coefficient dimension is

\[
d_{\rm sym}=1+M+M(M+1)/2.
\]

This count uses actual symmetric density entries: an off-diagonal rho2
parameter supplies both matrix orientations. In the repository's alternate
half-symmetrized coefficient convention its flat value would be twice the
off-diagonal matrix entry, with no change of dimension.

For a fully symmetric scalar amplitude, the isospin amplitudes obey
\(T^0=5A\), \(T^2=2A\), and \(T^1=0\). Consequently the functionals

\[
f^0_\ell-\tfrac52 f^2_\ell,\qquad f^1_\ell
\]

annihilate all \(d_{\rm sym}\) coefficient directions. Their combined
row space therefore has dimension at most

\[
d-d_{\rm sym}=M^2+M. \tag{IC1}
\]

At L allowed waves per isospin and M physical nodes, the complete augmented
family has 6ML+11 real rows. Apply the invertible row replacement
\(f^0\mapsto f^0-(5/2)f^2\) at every physical sample, separately to
real and imaginary parts. Do the same with the four chiral pairs:

\[
\chi_{01}(s)-\tfrac52\chi_{21}(s)
=f^0_0(s)-\tfrac52 f^2_0(s)
 -\bigl(r_{01}(s)-\tfrac52r_{21}(s)\bigr)f^1_1(s).
\]

There are **4ML+5 mixed-sector rows**: 4ML physical rows, four chiral
combinations, and the f11 target. They all obey (IC1). The other
2ML+6 rows require no symmetry assumption for the bound. Thus

\[
\boxed{\operatorname{rank}M_{\rm aug}
\le\min\{d,\ \min(4ML+5,M^2+M)+2ML+6\}.} \tag{IC2}
\]

A necessary condition for full augmented row rank is
\(M(M+1)\ge4ML+5\), in addition to the total variable count. At M4/L1,
21 mixed rows fit into only 20 available directions; (IC2) gives rank at
most 34. At the actual source M50/L10, the corresponding numbers are
2005 and 2550, so this argument supplies no rank obstruction.

## 2. Independent angular construction

The small controls integrate the rational crossed functions directly with
Arb's validated complex quadrature at 512 bits and a declared 1e-120
integration goal. They use

\[
q_i(v)=\frac{w_i v}{x_i(x_i-v)},\qquad
\frac14\int_{-1}^1P_\ell(\mu)q_i(t)q_j(u)\,d\mu,
\]

including the single-factor moments. Reflection parity supplies the paired
matrix entry and proves odd diagonal zeros. No angular quadrature node
approximation or numerical rank truncation is imposed.

The physical direct matrix is constructed independently through the exact
sine Gram identity, \(K=C(T^TT)^{-1}T^T\), rather than the new module's
cotangent formula. Its imaginary jump is the unit node selector. The
controls use the existing generic crossing assembler with the required
factor two on off-diagonal actual-rho2 columns.

All physical real and imaginary rows, eight raw chiral rows, the explicit
infinity functional, and both targets are retained. The infinity row in
subtracted density coordinates is

\[
\bigl(1,\ -b_i,\ -2b_i,\ 2b_i b_j,\
b_i^2\text{ on rho2 diagonals},\ 2b_i b_j\text{ off diagonal}\bigr),
\qquad b_i=w_i/x_i.
\]

## 3. Verified ranks and the retained first attempt

| Control | Complete matrix | Strict certificate | Consequence |
|---|---:|---:|---|
| M4/L1 initial full-square attempt | 35 by35 | residual bound about1.328e4; inconclusive | No full-rank claim; rank35 is subsequently excluded by (IC2) |
| Same M4 data, one34-square minor | First34 rows, selected34 of35 columns | residual below3.296e-70 | Rank exactly34 by the independent upper bound |
| M5/L1 direct angular control | 41 by51 | residual below3.792e-70 | Full row rank41 |
| M8/L2 direct angular control | 107 by117 | residual below8.942e-67 | Full row rank107, including non-S0 even waves |

Each successful lower-bound certificate uses one 256-bit inverse candidate
and a separate 512-bit full interval residual. Numerical QR only chooses
columns, retaining every declared row of that certificate. The complete
original matrices, selected matrices, candidates, row bounds and sources
are archived under `results/pv_midpoint_storage/`.

The M4 rank34 certificate excludes only the final f11 target row from its
minor; it retains all other 34 rows. Its first 33 base rows are therefore
independent. Adding both targets increases rank by exactly one. The
21 mixed rows also enclose zero on all 15 exact symmetric coefficient
directions in 315 independent consistency checks. Those interval zeros
confirm the implementation; the algebra above proves the rank upper bound.

No precision increase or repeated angular integration was used to address
the initial M4 result. The structural argument changed the appropriate
rank question. M5 and M8 were separate predeclared complete controls whose
sector capacities permit full row rank; neither result is transferred to M50.

The M4 and M5 angular constructions were additionally compared coefficient
by coefficient against the new closed-form kernel module. All **3316 real
entries overlap**, including physical, chiral, infinity and target rows.
[The comparison](../pv_midpoint_small_complete_comparison.json)
retains its source snapshots and input hashes. The broader independent
[angular review](73_cauchy_angular_projection_review.md) checks the exterior
branches, subtractions and nonzero-spin formulas separately.

## 4. Scope and next gate

For the M5 and M8 finite prescriptions, full augmented row rank implies
rank two of the target map after fixing all other retained rows. Pure
scattering samples can be assigned S=0 with zero chiral/infinity residuals,
so their two-target projections are the entire plane. This is a small-model
statement. It neither establishes the current-sector premise at these
different node counts nor reproduces the original source-sized model.

The source M50/L10 rank remains the next gate. Its local Cauchy moment
spaces and the full 3876-coordinate operator must be evaluated and certified
at that size. The PV/midpoint rows remain a discrete prescription; their
finite physical assignments are not asserted to be the boundary values
of their off-cut rational functions.
