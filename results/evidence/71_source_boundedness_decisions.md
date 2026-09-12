# Next source decision after the finite target-plane theorem

**Next derive the original-node version of the authors' principal-value /
midpoint Cauchy-kernel operator, then check its target kernel before any
boundary optimization.** Keep amplitude regularization as a separate,
explicit hypothesis. The original 2023 text does not identify its use or
size, but prior and later primary sources show that it is a substantive
part of the authors' numerical methodology.

The [rank-3011 theorem](69_full_analytic_rank_and_finite_target_plane.md)
settles the unregularized, exact finite-sine completion currently audited.
Changing a solver tolerance or chiral norm cannot turn its entire target
plane into the published bounded region. This does not settle a different
finite amplitude operator or continuum unitarity.

## A. What 2309.12402v3 explicitly fixes

The [archived TeX](../../references/2309.12402v3-source/prd_submission_2.tex)
and [prescription ledger](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/44_source_prescription.md`) give the following.

| Item | Printed restriction or remaining omission |
|---|---|
| Amplitude family | Adef and Svars retain T0, two independent M-vectors, one unrestricted M-square density and one symmetric M-square density: 3876 real variables at M50. A total-degree cutoff is not this count. |
| Double support | The `intxy` macro, TeX418, is the full rectangle `[4,infinity)^2`. No Landau-domain mask, coefficient deletion or positivity condition on these amplitude densities is supplied. Current spectral positivity is a different constraint. |
| Analytic and crossing construction | Adef specifies the dispersion representation; fdef specifies angular projection by `integral P_l/4`. The source maps zero to the disk center and uses midpoint angles. It does not print the amplitude-density interpolation or complete physical/off-cut integration matrices. |
| Unitarity | `uni` states the continuum condition for every physical energy and spin. The numerical results specify M50 and ten waves per isospin; Appendix A compares L8/10/12 and M45/50/60. No additional between-node grid or analytic large-spin/tail enforcement is described. |
| Subtractions/asymptotics | T0 is a free amplitude parameter. The explicit charge subtraction and high-energy suppression concern form factors. No extra amplitude sum rule, fixed subtraction constant or coefficient norm is stated. The added `A_infinity=0` condition is already included in the rank theorem. |
| Low energy | Four raw linear chiral residuals per family, with an unspecified norm. The text explicitly avoids fixed scattering lengths and does not fix the full amplitude through f_pi. |
| Numerical restrictions | No active specification of regularizer, rank truncation, coefficient bounds, solver precision or optimizer tolerances. The Tikhonov and CVX bibliography entries are commented out. The active 2021 regularization-paper citation occurs in the scalar-bound footnote, not as a stated prescription for the gauge figures. |

Thus support masks, amplitude-density positivity, suppressed highest modes,
fixed absolute chiral amplitudes and additional high-energy amplitude
conditions remain new restrictions unless further original-code evidence
is obtained. No such evidence was found in this bounded review.

## B. What cannot restore boundedness by itself

For the audited operator, the target-changing kernel preserves every
retained scattering value, all eight chiral residuals, and the infinity
functional. Therefore none of the following alone removes it:

- Changing the raw chiral norm, its grouping or its nonnegative tolerance:
  the theorem realizes the zero residual vector.
- Changing current/FF/FESR constraints that depend on the amplitude only
  through those same sampled S0/P1 values, provided a compatible current
  fiber remains nonempty. Such changes may empty the model; they cannot
  bound a nonempty target fiber while leaving those amplitude directions.
- Choosing another radial center, vertical slice or nonzero linear target
  objective; changing physical units; or applying invertible coefficient
  coordinates and positive cone rescalings.
- Rewriting subtractions while retaining all free constant/single/double
  families. This is an invertible redefinition when the underlying
  function family is unchanged. Fixing an additional quantity would be a
  different operation.
- Improving arithmetic or angular integration toward this same exact
  finite-sine operator. The obstruction is now analytic, not a Float64
  rank artifact.

Changing the node set, amplitude completion, extra amplitude observables
or continuum constraints requires a fresh argument. A variable/constraint
count alone is not a substitute for that argument.

## C. Two concrete primary-source leads

**A different finite operator.** The later primary paper explicitly
constructs angular-projected Legendre-Q kernels, midpoint off-cut sums,
and a discrete principal-value Hilbert matrix in sections 8.1–8.3.
The local notebook implements these as K, W0, W, W00, Ap/Am and their
unphysical-region counterparts. This is a concrete construction to derive,
with the original map center restored to zero; it is not evidence that
the 2023 finite operator was identical.
[Later paper, section 8](https://arxiv.org/html/2403.10772v2#S8),
[archived notebook](../../references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb).

The distinction is already visible without a numerical scan. With
zero-based source nodes, `sigma_j=(-1)^j`, the exact identity from
[chapter 51](../runs/repo_reorganization_20260906/recovery.json) (historical member: `docs/51_analytic_imaginary_structure.md`) is

\[
h_j^{\rm sine}(\nu)=\frac{x_j}{M(x_j-\nu)}
\left[\tan\frac{\phi_j}{2}
-\frac{\sigma_j}{2}z(\nu)^M\sqrt{4-\nu}\right].
\]

Literal midpoint integration of the unweighted sampled density instead
gives `h_j^mid(nu)=x'_j/[M(x_j-nu)]`, exactly the first term above.
The second term is not a constant subtraction or an invertible change of
coefficients. It also cancels the apparent source-node poles in the exact
sine expression. Consequently a nodal PV rule plus off-cut Cauchy sums
must be checked for their shared analytic interpretation; one cannot
declare them the same finite analytic amplitude merely because their
sample counts agree.

The next derivation should retain all 3876 variables and original nodes,
apply the angular projection analytically before midpoint sums, and
derive the physical diagonal jump/PV factors from the original dispersion
normalization. It should track rho2 symmetry weights and use the same
sampled scattering operator in ordinary and current constraints. First
compare simple density-mode identities and determine the target kernel;
do not optimize or transplant the later regularizer simultaneously.

**An actual amplitude restriction.** The pinned MATLAB snapshot has
`norm(rho/Mrho,4) <= 1e2` in both optimization passes, with
`Mrho=3775`. In ordinary vector-norm units this is
`||rho||_4 <= 377500`, not 100. The later paper prints an unnormalized
norm bound of 100 in Eq.(8.51); the norm convention must therefore be
recorded rather than equating those numerical statements. Its distinct
htilde/hhat arrays in Eqs.(8.48)–(8.50) are only positive cone rescalings,
not additional physical restrictions.
[Executable snapshot](../../references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m),
[later paper, section 8.4](https://arxiv.org/html/2403.10772v2#S8.SS4).

The prior authors' paper discusses why double-density directions require
regularization and presents both interpolation-point and sine-basis
methods. Its Eqs.(5.26)–(5.27) give precisely the sine densities and
`Phi_n=z^n-(-1)^n` used in the current completion. That supplies historical
support for both methods and for regularization, but does not establish
the missing 2023 implementation parameters.
[Prior primary paper, sections 3 and 5](https://arxiv.org/html/2103.11484v3).

There is a direct reason such a bound can remove the present obstruction.
Within the audited sine completion at fixed M, bounding both double-density families bounds their B/C
coefficients under the finite tensor transform. The sampled imaginary
D/E S waves then bound `T a` and `T b` after subtracting bounded double
terms. Invertibility of T bounds all 100 singles; one real S0 row, whose
constant coefficient is 5/2, bounds the final constant. Thus any finite
full double-density norm bound makes the pure/chiral amplitude feasible
set bounded, and its closed finite constraints make that set compact.
This argument is ours; it does not determine a source regulator value.
It also does not make the entire gauge feasible set compact: unconstrained
high-energy current-density directions can remain unbounded, so closedness
of its target projection and attainment need separate treatment.

## D. Decision and provenance limits

The source-supported next step is the **PV/midpoint amplitude-kernel
derivation and its early target-kernel test**. In parallel at the level of
definitions, specify the later-code density norm in the original sampled
density coordinates, including normalization and symmetric-array weights.
Do not choose its value by fitting the bounded plot. If that restriction
is eventually used, label it a regularized source variant until its 2023
provenance is established.

Continuum scattering constraints are physically motivated by the paper's
stated problem and remain a separate route to boundedness. The existing
edge, fixed-spin and joint-limit analyses can guide that route, but their
addition would not establish which finite scheme generated the published
M50 figures. Solver rank truncation or implicit suppression of large
directions is an unprinted numerical possibility, not an acceptable
unannounced replacement for those restrictions.

The [author repository README](../../references/upstream-gauge-theory-bootstrap/README.md)
identifies the 2023 paper as the original proposal without public code at
the time and supplies later snapshots. The local commit was revalidated
as `801684d9ece3de2918a20145178f569459081098`. The original TeX still has
SHA256 `39ec9e10da96f9c09c8db9e7d3ac63f342b1214df592896281adda727fb6be84`;
the inspected MATLAB and notebook hashes are respectively
`b12bf411407fc4cb8ae56aea41c88aa336e0f73aa7ae44ebef2b1bca3fc08ad2` and
`9a3a81d7e99ef35355aa5a2bbf10ce9a21e6af929e0cb8371c4afaff133922bd`.
The primary 2021 and 2024 HTML versions were consulted only to resolve
their additional method details. No outreach, optimizer, parameter scan,
phase fit, or source-code modification was performed for this memo.
