# Scientific contract — He–Kruczenski 2309.12402v3

Updated: 2026-09-09. The target is [Bootstrapping gauge theories, v3](../../../references/2309.12402v3.pdf), especially Fig.3–11. [The version record](https://arxiv.org/abs/2309.12402v3) dates v3 to 2024-12-04; “2023 paper” alone does not identify this numerical target. Snowmass is background. Progress and the remaining deliverables are tracked only in [STATUS.md](../../../STATUS.md).

## Current source-identified chiral grouping (2026-09-09)

C selection is now decoupled from unrelated global polygon gaps. For its audited complete point (xc,yc), any valid upper plane nx*x+ny*y<=U with ny>0 bounds ymax(xc)-yc by (U-nx*xc)/ny-yc. `select` computes the minimum bound with Arb and scales by epsilon; the same .01 metric budget can certify this local representative while B3's global geometry remains unfinished. Both scopes are recorded; no phase data enter the choice. See [the derivation](../../../results/runs/mainline_alignment_20260909/C_LOCAL_SELECTION_ZH.md).

The [original Fig.5 residual reconstruction](../../../results/runs/mainline_alignment_20260909/CHIRAL_SOURCE_RESIDUALS_ZH.md) identifies a much stronger source fingerprint than the words “some norm”: the three colored curves have combined eight-component L2/epsilon approximately1.000007,1.000020,1.000146. Interpolation uncertainty is retained; this is not a strict violation test or exact producer recovery. Together with the author's explicit later stack8 implementation, it motivates `combined-l2` as the current CLI default and mainline choice. Historical separate-L2 replays must request it explicitly.

[PAPER_MAINLINE](../../../results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json) now fixes combined-L2 and the uniform hard-midpoint quadrature below. Pure and H/current operators are unchanged; B/C and joint points must be re-audited under the new ball. Reusing a dual requires evaluating its support function for the new norm, and scaling a scattering-only amplitude requires rechecking the resulting full coefficients. A complete gauge point is never rescaled as if F(0)=1 were homogeneous. The [normalization audit](../../../results/runs/mainline_alignment_20260909/NORMALIZATION_TRANSFER_ZH.md) confirms that the actual-amplitude rho bound requires no additional pi or overall factor2 conversion.

## Printed-quadrature source decision (2026-09-09)

The new [cutoff decision](../../../results/runs/mainline_alignment_20260909/CUTOFF_DECISION_ZH.md) adopts `hard-midpoint` for D–F: the uniform pi/M weights printed in v3 Eq.(3.72), masked to native s<=s0. Previous `clipped-phi` results remain a named endpoint-corrected variant. Masses, s0, chi, actual-density bound, moment targets/errors and FF caps are unchanged. This is a choice of the direct published discretization, not recovery of the unpublished producer. At M50 the exact-weight hard (S,F) projection is contained in clipped by monotonically increasing the latter's endpoint current spectrum; floating matrices still require fresh audits. M45/M60 have different positive-weight supports, so that proof does not transfer.

The [hard D witness](../../../results/runs/mainline_alignment_20260909/D_RESULT_ZH.md) passes all original joint constraints. New E near-black x positions are declared in [HARD_MAINLINE](../../../results/runs/mainline_alignment_20260909/HARD_MAINLINE.json) before new phase evaluation. `--selection-reference` reads the original Fig.8 marker geometry and scales its two near-black x positions by XREF/x_black; it does not use phase data or claim original amplitude identity. The default selection remains available for historical replays.

Numerical trials of unit chi-barrier weight, wider QR reuse and a unit step cap did not close the old mid problem; the original QR/line-search strategy was restored. `--chiral-barrier-weight` is explicitly numerical; omitted values preserve checkpoint weights or the established n_disks/2 default. A convex mixture of boundary points is not an analytic center. Production now uses complete center paths from appropriate joint interiors. Converged full vectors are retained in center_path.npz and their native S/F changes in center_observables.json: projection gap does not bound those changes. The author explicitly permits repeating2403 Watsonian updates until convergence; old-Q objective optimality remains distinct from that iteration's fixed-point behavior.

## Earlier A–F claim alignment (2026-09-09)

The [18-step review](../../../results/runs/claims_alignment_20260909/CLAIMS_ZH.md) supersedes any reading of “finite chain delivered” as completion of all physical claims. Original E three-point robustness, mid/high-energy S0 and parts of Fig.11 stability remain unresolved. The [definition audit](../../../results/runs/claims_alignment_20260909/DEFINITION_REAUDIT_ZH.md) found no supported correction to chi residuals, conjugations, rho2, FESR units or squared FF caps. Prepare produces eight chi residual rows and imposes no norm; its formerly hardcoded two-ball metadata has been corrected. The production solver grouping and historical amplitudes are unchanged.

The [source review](../../../results/runs/claims_alignment_20260909/SOURCE_REAUDIT_ZH.md) distinguishes explicit v3 equations, unpublished numerical choices and later methods. In particular s0=(1.2GeV)^2 is explicit; f_pi is used to locate the target region; v3 does not publish the regulator or Newton algorithm. The spectral-moment-ratio mass estimate occurs only in a commented TeX line. A moment ratio remains a derivable spectral diagnostic, not a published v3 rho-mass definition.

`compare --comparison-manifest ...` in `selection.py` compares saved curves with declared reference filters. The [48 comparisons](../../../results/runs/claims_alignment_20260909/CURVE_COMPARISON_ZH.md) separate paper bootstrap from phenomenology, original E projection points from later Watsonian points, and raw-GeV from identical native-s comparisons. A native-prefix mapping is inferred from reference x coordinates alone, never by minimizing phase errors. The inferred display mass near139.57MeV and black point near139.57/93 do not change the production140/92 inputs. RMS values are descriptive marker-sample differences; no statistical or continuum accuracy is inferred. Modulo180 differences diagnose branch ambiguity without changing any reported phase.

The [selection audit](../../../results/runs/claims_alignment_20260909/SELECTION_REAUDIT_ZH.md) shows that the original mid borrowed a valid ref support plane but did not converge its own fixed-x objective. The later Watson-ref precision check cannot exclude this upstream issue. Its separate [fixed-model precision result](../../../results/runs/claims_alignment_20260909/MID_PRECISION_ZH.md) preserves the original position and all constraints. F's existing five supports test the original UV+x selection only; they do not certify the later Watsonian trio's resolution stability.

## E closure audit and the named follow-up method (2026-09-09)

The [independent operator audit](../../../results/runs/stage_E_closure_audit_20260909/OPERATOR_AUDIT_ZH.md) reconstructs Eq.(2.7) angular projections without calling the production projection kernels. Representative M50 rows agree across all3876 coefficients within5.33e-15. The original complex Gram determinant, conjugations, current dimensions, FESR Jacobian and squared FF cap have independent checks; no major normalization error was found. This is an audit of the declared finite model, not recovery of unpublished v3 settings.

The [primary-source audit](../../../results/runs/stage_E_closure_audit_20260909/SOURCE_AUDIT_ZH.md) identifies the targeted continuation in2403.10772v4 Eq.(2.29) and the authors' public MATLAB snapshot. For an old feasible full joint solution define Q=F_old/F_old* in S0/P1 and Q=S_old in S2. Optimize the linear functional

\[
W=\frac1{3N_{\rm low}}\sum_{I=0,1,2}\sum_{4<s_j\le s_0}
\Lambda_{\ell}(s_j)^{-2}\Re[Q_{Ij}^*(S_{Ij}^{\rm new}-1)],
\qquad \Lambda_0=1,\quad \Lambda_1=\frac{\sqrt{s-4}}{\sqrt s+2}.
\]

This is the three-primary-wave specialization for2309 (S0/P1/S2); the later six-wave objective and D0 current are not claimed to be reproduced. Every original scattering disk, chiral ball, actual-density bound, current Gram, FESR and FF cap remains imposed. The positive overall normalization changes no maximizer. As in the author's second CVX block, the original projection equalities are released: later curves must carry both their initial geometric role and their new projected coordinates. Q comes from the prior computed FF, with no experimental phase input, branch choice or fitted resonance mass. For any feasible S the elementary upper bound is sum(w*(abs(Q)-Re(Q)))/(3*N_low); a smaller jointly audited dual bound may certify the actual optimum even when this ideal disk ceiling cannot be reached.

This is an explicitly named follow-up algorithm, not a step documented by2309v3. It improves the search for the Watson/saturation behavior underlying the physical argument; PSD alone permits eta<1 and does not force S=F/F*. Changing the optimization objective cannot change the Fig.8 feasible region. The2505 update has a different iteration and UV inputs and is not silently substituted. `--objective watson` starts a new frozen functional; `--resume-objective` continues the same one. A solver time limit is not a completed Watson iteration.

A useful exact consequence of the original Gram determinant is

\[
|S-Q|^2\le (1-|S|^2)\left(\frac{\rho_{\rm current}}{|\mathcal F|^2}-1\right),
\qquad Q=F/F^*,\quad F\ne0.
\]

Thus numerical eta close to1 does not by itself give a uniform Watson phase bound when the two-pion fraction r=|mathcal F|^2/rho_current is very small. In the completed mid continuation the first P1 node has eta=0.99999999853 but r=2.75e-9 and |S-Q|=.9184. This satisfies the Gram constraint and exposes the role of spectral dominance in the paper's physical argument. It is not evidence of a lost unitarity disk or justification for replacing the input spectrum. Exact eta=1 with nonzero F still forces Q=S.

For a uniform additional phase diagnostic, define delta_FF=unwrap(arg F_new)+arg(S_new/(F_new/F_new*))/2, with the real nonzero F(4) fixing the overall integer-pi threshold anchor. This exactly preserves eta*exp(2i*delta_FF)=S at the evaluated nonzero S/F samples. The known FF family gives F(4)=1+sum_j cot(phi_j/2)*ImF_j/M. Both lifts and their integer-pi difference are reported; the primary Fig.9–10 curves retain the nearest native S lift. Exact elastic Watson equality, continuity and treatment of zeros/branch crossings are additional requirements for a unique physical interpretation. The2403 discussion treats4<s<16 as the elastic region; above it saturation is approximate. Native inequalities alone do not impose exact saturation, and a completed old-Q optimization does not imply a new-S/new-F fixed point. FF analyticity cannot supply the missing PV scattering winding between nodes. This diagnostic is fixed uniformly before viewing the new phases, not selected by a desired rho mass.

The later author's combined eight-component chiral L2 ball is available as `--chiral-norm combined-l2`, with matching primal/dual math and result identity. For the same epsilon, K_comb(epsilon) is a subset of K_sep(epsilon), which is a subset of K_comb(sqrt(2)*epsilon). The default and current A–E baseline remain `separate-l2`; the new branch is a source diagnostic, not a retroactive correction to v3. A new complete M50 joint witness satisfies the combined constraint. No old two-ball feasibility certificate is transferred, and no new six-epsilon or C scan is required while that branch is not adopted. F is now following the fixed resolution protocol below.

## F resolution protocol and numerical Phase I (2026-09-09)

The [Fig.11 source audit](../../../results/runs/stage_F_mainline_20260909/SOURCE_AUDIT_ZH.md) fixes exactly (M,L)=(50,8),(50,10),(50,12),(45,10),(60,10). L counts waves per isospin. The original appendix does not specify the exact representative functional or a regulator rescaling rule. Before inspecting new phases, this reproduction chooses the full joint +f00(3) support for all five, with no Watsonian step, and the declared discrete recipe B(M)=100[M²+M(M+1)/2]. This is a conditional common recipe, not a recovered continuum norm. All chiral/current/UV inputs remain the E projection baseline. Each support is checked to gap<=1e-4.

Cross-M density projection uses the sine Galerkin map in subtracted coordinates, then returns to original C_flat; this avoids losing large subtraction cancellations when modes are dropped. It is only a numerical initializer. Neither a PV amplitude identity nor feasibility is transferred. Cross-L amplitudes retain their complete variables; every added disk is enforced. At fixed M all shared old/new rows agree exactly, but the old L10 tip violates some added L12 high-energy/high-spin disks, so simple inheritance is invalid.

The initializer now follows B→D explicitly. Construct the same strict analytic scattering seed used in B on the actual M/L grid, with the chiral and actual-density bounds. Every scattering disk, chiral ball and L4 bound remains hard throughout Phase I. Only the current subsystem receives one scalar tau: normalized Grams G+tau*I, FF slacks 1+tau−||F_scaled||², and moment halfspaces 1+tau±r where r=(moment−target)/error. A same-resolution B amplitude can initialize this step; a small analytic-seed mixture preserves strict scattering slack. Current Grams use the same positive diagonal equilibration as the full joint support solver. This is a fixed invertible congruence from the numerical reference, not an extra physical condition. Earlier unbalanced current-only runs and their checkpoints remain separate. Minimizing tau is a standard convex epigraph Phase I. Negative tau guarantees all original current constraints; acceptance uses the independent original joint audit. No Phase-I dual is reused for support. Failed projected pools and earlier all-constraint slack attempts remain numerical failures, not original-model infeasibility proofs.

For a fresh Phase I, start with mu=max(.001,tau_initial/nu), nu=3*N_Gram+N_FF+8, balancing the initial objective and current barrier; decrease it in a Newton neighborhood (decrement squared below1); it does not require support-level central accuracy. Its checkpoint is distinct from the unrelaxed support checkpoint. Both stages use the existing feasible line-minimization expansion; a trial cap at a full Newton step did not improve the calculation and was removed. Support resumes follow half-mu continuation from a center, avoiding a premature tenfold jump. These are numerical choices, without altered physical constraints or phase data. Independent directional first/second derivative tests cover both barriers.

A Phase-I checkpoint also supplies candidate covectors for the original zero-objective problem: the hard amplitude/chiral covectors, positive mu*(G+tau*I)^(-1) with its saved congruence, the existing FF covector using shifted FF slack, and moment weights mu*(1/ms_minus−1/ms_plus)/error. The original joint audit includes all residual support terms and no tau constants. Only a strictly negative original Arb upper bound certifies finite infeasibility; a positive bound remains inconclusive. These covectors are saved separately and never presented as +x support weights.

F compares fresh native S0/S2/P1 values and eta on each actual M grid. Across-M interpolation is a plot/comparison diagnostic only; it does not evaluate or certify an off-node PV amplitude. The original Fig.11 vector marker CSVs are reference readouts, never optimization inputs. Finite support error, resolution sensitivity and the remaining E physical discrepancies are reported separately.

## F delivery (2026-09-09)

[The F result](../../../results/runs/stage_F_mainline_20260909/F_RESULT_ZH.md) completes all five full joint +x supports to original gap<=1e-4 and all669 fresh native primary-wave checks. Across the five finite models,7650 scattering disks,510 current Grams,10 chiral balls,5 actual-density bounds,20 FESRs and74 FF caps pass. All five P1 curves cross90 degrees, with descriptive readings840.59/812.36/798.34/752.12/773.79MeV in the declared configuration order. Fixed-M L variation spans42.25MeV versus12.54 in the original plot; fixed-L M variation spans60.25 versus54.35. S2 is stable and low-energy S0 remains close across resolutions, while the absolute mid/high-energy S0 discrepancy persists. No all-paper-claims or continuum-certification assertion follows.

The unused cross-resolution projected hull/SOC4 initializer is removed; the original same-model D hull remains. Failed numerical inputs and producer source snapshots are retained. The independent zero-objective Farkas control exposed and fixed a NaN when a derived current-spectrum upper box was negative: replace that box by the weaker nonnegative bound before taking a square root, retaining the original FESR in both primal checks and dual constants. Positive physical-input runs are unchanged. All120 tests pass; the active source has15 focused modules within the existing size limits.

## Acceptance update and C mainline (2026-09-08)

The user's acceptance criterion is alignment with the paper's principles and core scientific comparisons, not exact recovery of every displayed point. A declared implementation of the unspecified norms is admissible; measured numerical differences remain visible. No new parameter/solver scans or regulator-platform prerequisite is introduced. See [the Fig.4 audit](../../../results/runs/stage_C_mainline_20260908/FIG4_AUDIT_ZH.md).

C1 uses three preselected upper-boundary approximations at the same physical x_ref=5/[16*pi^2*(92/140)^2], epsilon=.006/.004/.002. Each is a convex mixture of two stored feasible amplitudes; all3876 coefficients remain and the rounded mixture is checked against the original1500 disks, density ball and chiral balls. The existing regional error controls the projected boundary approximation, not uniqueness or phase stability of an exact extremizer. Point selection never uses phase data.

For S0,S2,P1, the Weinberg reference is ((2s-1),(2-s),(s-4)/3)/(16*pi^2*(92/140)^2). Direct polynomial-amplitude angular projection independently tests these factors. On the81-point subthreshold grid the deviations decrease with epsilon; the blue sample has no observed S0 sign change, orange/green zeros approach s=.5. These are results for the selected amplitudes, not a theorem about every allowed amplitude.

The actual Fig.6 magenta coincides with the Fig.4 green highlight within readout precision. C2 therefore locks the same green coefficient vector used in C1. The caption's word closest does not specify a distance metric; it does not justify switching to a different branch by a newly invented minimization rule.

C3 uses that same vector at43 original M50 nodes through1.2GeV plus threshold. The phase is unwrap(arg S)/2, anchored by S(4)=1, and eta=abs(S); eta=1 is not imposed. No off-node physical PV continuation or sine substitution is made. The three primary waves pass the132 listed checks; the original1500-disk check remains available from selection. S0/S2 broadly reproduce the low-energy pattern while P1 lacks the observed rho behavior. S0 deviations, including the region near1GeV, remain in the comparison. One wrong P1 representative shows that these IR inputs do not force the correct P1; it does not exclude resonant amplitudes throughout the chiral feasible set.

[The C delivery](../../../results/runs/stage_C_mainline_20260908/C_RESULT_ZH.md) contains coefficients, phases/eta, plots and replay parameters. PDF fit comparisons are descriptive native-sample errors, not chi-squared or tuning objectives; reference curves are not extrapolated. No new bootstrap optimization was used. D subsequently completes the joint witness described below; E/F retain the same declared scattering setting for the IR/UV comparison.

## E delivery update (2026-09-09)

[The closure delivery](../../../results/runs/stage_E_closure_20260909/E_RESULT_ZH.md) completes the bounded E calculation and audit; all-paper-claims remains false. A new lower support and Arb feasible-segment/support-plane section prove that the upper shrink is larger at xref, with ratio[1.026389,1.167740]. The earlier unresolved sign below is historical. The original display readout near5.35 is still unmatched, and no dense UV geometry or continuum certificate is claimed.

All three frozen initial roles have completed one named2403 Watsonian optimization in the same4076-variable finite model. Tip/mid/ref support gaps are .000621638/.000663697/.000920794. Their new coordinates and complete vectors are saved; the original projection-only vectors are retained. The three nearest-S phase curves now all cross90 degrees, with descriptive readings .807565/.728386/.697596GeV. Their separation, and S0/S2 differences, remain explicit. These moved points do not retroactively reproduce the original three boundary amplitudes' claim.

A tighter ref check completes the current mu center with gap .000646807. Maximum P1 phase change is .08271 degrees, its crossing changes .029MeV, and minimum eta changes .84699 to .84949. The frozen primary trio is not replaced. This excludes that particular unfinished center as an explanation of the large discrepancy, not every numerical/discretization effect. Watson fixed-point diagnostics are reported without making exact saturation at every energy a new paper acceptance gate. F is now following the fixed resolution protocol below.

New Watsonian problems default to mu=.001; same-objective resumes read their saved mu unless explicitly overridden. Existing valid duals are carried forward, and a candidate close enough to an earlier valid upper bound is sent to the original audit before demanding another full center. QR reuse is bounded and checked against the actual linear residual; all physical constraints stay unchanged. The full native conic trial did not improve M50 and was removed from active code; its source and failed candidates remain in the run records.

## Earlier E projection-only delivery (2026-09-08)

[E delivery](../../../results/runs/stage_E_mainline_20260908/E_RESULT_ZH.md) completes the full finite calculation chain while explicitly retaining failed physical comparisons. D's PV/M50/L10, all3876 amplitude coefficients and200 currents, two radius-.002 chiral balls, ordinary actual-density L4 bound377500, printed/raw FESR and clipped-phi/FF settings stay fixed. The +x joint interval [.09924805847073,.09925920163675] is strictly below the same-setting IR endpoint near.1059. At xref, however, the upper/lower shrink ratio is [.984,1.202], compared with about5.35 from original Fig.8 marker hulls. Sparse inner/outer regions are delivered; dense UV contour convergence is not certified.

The fixed numerical coordinates are z0=x=f00(3) and zj=y+x/15, with actual amplitude rho unchanged. For a fixed-x support of d0*x+d1*y, eliminate only the z0 increment. If phi is the joint barrier, lambda=c0-mu*partial0(phi) gives the global normal n=(d0-lambda,d1). An original-operator audit of U implies d0*xref+d1*y<=U+(d0-n0)*xref. On a projection ray y=a*x, eliminate dzj=(a+1/15)*dz0; its multiplier gives n=d+lambda*(a,-1), whose objective agrees with d on that ray. No amplitude density direction is discarded beyond the prescribed section/ray equality.

A saved valid joint dual can audit another complete candidate under the same operators. For n_y>0, (U-n dot p)/n_y bounds the candidate's distance below the upper boundary at its own x. E uses this for the midpoint, whose y-distance bound is3.3244e-5; its internal center did not fully converge. The reference point has a6.77e-6 bound. These are geometric bounds, not phase-error estimates. The best verified primal and the best valid dual are retained independently, then audited together. Lower-side comparison uses feasible points and valid outer bounds; it does not impose an identical arbitrary support gap on every point.

Initialization combines complete audited joint points, including all currents, or first constructs the D-type hull witness and then centers the full joint problem. For a complete convex mixture, q(S_mix)=alpha*q(S1)+(1-alpha)*q(S2)+alpha*(1-alpha)*|S1-S2|²; Grams are affine. A substantial contribution from an existing joint center preserves its weak-constraint slack scale. A zero objective at fixed x minimizes the same full barrier and stops only after centering; it is never reported as support optimality. This numerical initializer is distinct from a final phase representative.

The joint Newton solver uses QR-preconditioned longdouble CG, an invertible increment change isolating the two stiff chiral radial rows, and exact quadratic slack updates. At M50/L10 each chiral log barrier has weight750, equivalent to repeated identical constraints: radius .002 is unchanged, and value, gradient, Hessian and exported dual all carry the weight. The feasible line minimization expands a derivative bracket then bisects; a center is solved before decreasing mu. Actual M50 directional gradient/curvature contractions and independent finite differences agree. Step caches retain coordinates, slacks and the Gram reference; same-section objectives may change, while operators and physical tolerances must agree. These are numerical choices, not a claim that v3 used this Newton implementation.

E2's tip, xref and midpoint rules preceded the new phase calculations. The three complete vectors remain frozen, with no phase ranking or final amplitude mixing. All original joint constraints pass, and three fresh primary-wave evaluations pass396 listed native checks. The tip's first P1 90-degree bracket is [.792136,.864418]GeV, with descriptive linear reading .812362GeV (+5.50% vs770MeV). Ref crosses in [.680414,.731675]GeV, reading .699955GeV (-9.10%); mid has no upward90 crossing in the fixed threshold-anchored nearest-branch unwrap. Thus the tip reproduces the paper's characteristic example, while three-point robustness is not reproduced.

Native S fixes delta only modulo pi. Mid's absent crossing does not prove absence of a continuous rho pole, and the PV provider supplies no physical off-node winding to settle that ambiguity. The unwrapped P1 high-end spread near190.7 degrees partly reflects a pi branch difference. Branch-invariant pairwise angles are about89.5 degrees around.792GeV and11.8 degrees at the high end; complex-S separations are also retained. Eta is not forced to one: its P1 minima are approximately .818/.332/.174 for tip/ref/mid. S0/S2 spread more at high energy, but also have material quantitative discrepancies, including S0 already around.51GeV. See the [phase audit](../../../results/runs/stage_E_mainline_20260908/E3_AUDIT_ZH.md).

A concrete allowed freedom appears in the P1 moments. Raw absolute error .002 is about47.3% of the printed M_minus1 and1.373% of M0. The three solutions put M0 near its common upper limit but choose different M_minus1; .14*sqrt(M0/M_minus1) is about .810/.694/.741GeV for tip/ref/mid. This is a positive rho(s)/s-weighted spectral scale, not a pole or 90-degree mass, and not a proof of a unique cause of the phase differences. The raw reading was fixed before D; no tolerance was retuned. Original setting identity, extremizer/phase stability and resolution remain open scientific issues rather than grounds for silently fitting the output.

## D current/UV delivery (2026-09-08)

D1/D2 now have the actual PV/M50/L10 production operators and independent real/complex Gram wiring checks. Before any joint solve, D2 fixes clipped-phi cutoff-cell overlap, the printed moments, four raw absolute errors .002, epsilonFF=6e-5 and the arithmetic-mean quark-mass rule. This is an explicit finite quadrature choice, not recovery of an unpublished unique v3 prescription; the constant-spectrum diagnostic is not an error bound for an unknown spectrum. See [the D derivation and result](../../../results/runs/stage_D_mainline_20260908/D_RESULT_ZH.md).

D3 has a complete4076-variable joint feasible witness. The candidate generator uses the convex hull of all204 saved same-model B amplitudes, including pure parents, and explicitly imposes the requested two chiral balls while solving the current SDP. This is a constructive existence calculation in the original M50 model: convexity preserves ordinary scattering disks and the original actual-density L4 ball, and the final independent audit nevertheless checks all original constraints on the rounded full coefficients. Nonnegative normalized weights are selected without phase or experimental data. It is not an extremizer or a replacement for the full-amplitude E1 support problem.

Only14 unconstrained high-energy current spectra are eliminated, then recovered by the Schur complement for strict |S|<1; all100 original Grams are checked afterward. Candidate generation uses a declared1e-6 inner cone margin; acceptance uses the original unshrunk constraints. Unitarity, chi, density, current PSD, moments and FF must all pass jointly. An amplitude-only feasible flag is insufficient. The zero direction used for witness audits is never reported as support optimality.

The saved-operator768-bit audit passes1500 disks, both chi balls, the L4 bound, all700 principal minors of100 current Grams, four moments and14 FF caps. Fresh768-bit native PV evaluation also passes all1500 waves and chi; its largest complex-f difference from the saved H contraction is about5.10e-11 for this coefficient vector. These arithmetic checks do not certify continuum unitarity, physical phase precision or discretization error. E1–E3 now deliver full-variable supports, three representative boundary approximations and the rho comparison, with the physical differences and limitations stated above.

The earlier C audit also quantifies the convex-mixing contribution to inelasticity: its maximum extra1-|S|² in S0 is about.12159. This is consistent with unitarity but prevents interpreting the C representative as a high-precision unique extremizer. Parent-phase variation and regional geometry errors are not fixed-point phase error bars. See [the readiness audit](../../../results/runs/stage_DF_readiness_audit_20260908/AUDIT_ZH.md).

## Source audit and prototype scope (2026-09-07)

The primary deliverable is the complete finite prototype supporting v3's IR-to-UV comparisons, regions, phase shifts and resolution study. The [original-paper audit](../../../results/runs/stage_B1_source_method_audit_20260907/ORIGINAL_PAPER_ZH.md) separates the all-energy/all-spin axiom from the finite M/L calculation and empirical Fig.11 evidence. The [method review](../../../results/runs/stage_B1_source_method_audit_20260907/OUR_METHOD_REVIEW_ZH.md) identifies the present zero-infinity, five-tail, global-FG and expanded-grid model as our strengthened branch. Those conditions have valid derivations in our finite family; that does not identify them as the author's finite implementation.

A closed finite support problem, numerical stability of the paper's physical observables, and a continuum certificate are distinct deliverables. We will retain all existing violation evidence without treating a continuum certificate as a prerequisite for every finite current/UV implementation step. Both source and scope branches are implemented. `--prescription auto` selects sine for preparation and otherwise resolves the saved source. PV boundary/dual require `--unitarity-scope sampled --infinity free`; they do not apply sine FG or tail formulas. This executable prototype does not establish the unspecified v3 norm or regulator. The conditional B prototype and regional geometry are complete; v3 quantitative agreement remains unresolved. Stage completion is tracked in STATUS.

v3 does not specify a Newton solver. The later public MATLAB snapshot invokes CVX/MOSEK; MOSEK uses a [homogeneous primal-dual interior-point method](https://docs.mosek.com/latest/toolbox/solving-conic.html). Our barrier already analytically minimizes standard SOC density auxiliaries, uses a diagonal-plus-rank-one Schur complement, QR preconditioning and longdouble iterative refinement. These are numerical methods, not additional QCD input. Solver choice cannot establish equivalence of different interpolation families, regularization norms or constraint sets.

The [independent 2511.11513v1 implementation audit](../../../results/runs/stage_B1_source_method_audit_20260907/DISCRETE_BASIS_ZH.md) gives relevant evidence for separating basis truncation from energy sampling and measuring convergence. Its reported speed is not a reproduced benchmark here, and its altered UV inputs and density norm cannot be substituted for v3.

## B plan audit: finite delivery and source identity

The current B plan keeps one declared PV/M50/L10 contract and all original finite constraints. The later-added common 1%/two-decade regulator certification is retained as a failed sensitivity study, but removed as a prerequisite for this conditional finite prototype. Six epsilon regions still use the same B and norm. Existing strict-gap failures remain failures; valid wider bounds remain usable for the unchanged .01 regional geometry requirement. Original figure settings and quantitative agreement are reported separately, never inferred from calculator readiness. See [the bounded plan](../../../results/runs/stage_B_plan_audit_20260907/PLAN_ZH.md).

The feasible coefficient sets are convex: all sampled unitarity disks, both chiral L2 balls and the actual-density L4 ball are convex inverse images under fixed linear maps. Thus the convex hull of saved feasible amplitudes supplies an inner projected polygon, while each validated support inequality supplies an outer half-plane. With the same H and B, epsilon1 <= epsilon2 implies K(epsilon1) subset K(epsilon2) subset K(pure) directly from the constraints. The aggregate also checks all available smaller-set witnesses against larger-set upper bounds; no observed contradiction is a numerical consistency check in addition to this nesting argument.

For pure the declared metric is q=(x,y). For chiral it is q=(x/x_ref,(y+x/15)/epsilon), where x_ref=5/[16*pi^2*(92/140)^2]. The second coordinate resolves the thickness around the Weinberg line. An Arb upper bound on the distance of every outer vertex to the inner polygon bounds the entire outer-to-inner Hausdorff distance by convexity. The retained .01 threshold therefore measures the accuracy of the finite projected region, not its distance to the paper's marker hull. Support and geometry enclosures apply to the saved float64 finite operator; they do not enclose changes of discretization, unresolved source choices or all physical energies and spins.

Positive-span detection orders directions by half-plane and Arb cross products. Floating atan2 values remain diagnostics: two distinct, nearly parallel directions can round to the same angle and previously produced a false unboundedness report. Both ordering variants of a square with two redundant near-parallel inequalities now preserve its known enclosure. The [full M50 reaggregation](../../../results/runs/stage_B_delivery_20260907/geometry_order_check/regions.json) retains every support half-plane and restores closure; no optimizer, amplitude constraint or geometric tolerance was changed.

## Stage A: preserved finite-sine contract

A1–A3 were completed with the conditional finite-sine formulation below. It remains available alongside the active PV source prototype; neither establishes recovery of the unpublished v3 implementation. The local [author README](../../../references/upstream-gauge-theory-bootstrap/README.md) states that 2309 had no public code; the preserved 2403 snapshot is a later calculation.

| Choice | A finite-sine definition | Comparison / scope |
|---|---|---|
| Scattering family | all M sine-cardinal densities, including Nyquist; one analytic function for direct and crossed arguments | active PV/midpoint is a distinct collocation prescription, defined below |
| Coordinates / constant | unsubtracted Eq.(2.7); full operator retains T0; strengthened default uses `--infinity zero` | free-T0 is a sampled relaxation of this sine family; PV explicitly uses sampled/free |
| Double-density bound | ordinary l4 on actual rho1 and rho2 upper triangle, default B(M)=100[M²+M(M+1)/2] | `--density-limit` fixes B; a common stable B platform for the regions is not yet selected |
| Chiral norm | two separate l2 four-vectors, epsilon=.002 | combined eight-vector norm in later code is not silently substituted |
| FESR targets | the four printed (2.56) moments | `--moment-source eq250` recomputes from UV inputs |
| FESR errors | four independent raw-integral absolute errors <=.002 | `--sr-error normalized-absolute` changes the feasible set |
| Cutoff quadrature | `hard-midpoint`, using nodes <=s0 | `--fesr-cutoff clipped-phi` retains partial endpoint-cell width |
| Effective mq | arithmetic mean, with printed targets held fixed | `--mq-rule rms`, also affecting the scalar squared FF cap |

The sine transform has a direct derivation. For T_jn=sin(n phi_j), T^T T=diag(M/2,...,M/2,M). Define

\[
W_{nj}=(T^{-1})_{nj},\quad q_j(v)=\sum_{n=1}^M W_{nj}z(v)^n,\qquad
H_j(v)=\sum_{n=1}^M W_{nj}[z(v)^n-(-1)^n]=q_j(v)+b_j,
\]

where b_j=tan(phi_j/2)/M. H is the exact continuous Cauchy transform of the sine density interpolant. Hence q(0)=0 and q_j(x_k+i0)=K_kj+i delta_kj, while the same H is evaluated in all direct and crossed slots of Eq.(2.7). Every density variable is retained; the finite family assumes O(sqrt(x−4)) threshold and O(x^(−1/2)) high-energy density behavior in each argument. No midpoint rational crossed kernel is mixed into this family.

Writing R1 and R2 as actual density matrices, the change H=q+b is

\[
\sigma_1^s=\sigma_1^u+2R_1b,\quad
\sigma_2^s=\sigma_2^u+R_1^Tb+R_2b,
\]
\[
T_0^s=T_0^u+b^T\sigma_1^u+2b^T\sigma_2^u+2b^TR_1b+b^TR_2b.
\]

It leaves both double densities and their norm unchanged. Inverting gives

\[
A_\infty=T_0^u=T_0^s-b^T\sigma_1^s-2b^T\sigma_2^s+2b^TR_1b+b^TR_2b.
\]

Thus infinity=0 removes a genuine direction and is not a subtraction identity. It is absent from v3's explicit finite setup and the later MATLAB CVX constraints, but is a necessary consequence of full-energy unitarity **within this finite sine family**. Here each H_j tends to zero at large physical and crossed arguments and is uniformly bounded. Dominated convergence in the angular projection gives

\[
\lim_{s\to\infty}f^0_0(s)=\tfrac52T_0^u,\qquad
\lim_{s\to\infty}f^2_0(s)=T_0^u.
\]

The coefficients are real, so these limits are real. Since kappa tends to pi, either nonzero limit would give |S|²→1+pi²f_infinity²>1. Therefore T0^u=0 is necessary; in subtracted coordinates impose the full A_infinity functional above, not T0^s=0. The original A operator deliberately retains the constant so both branches remain representable. Current `boundary` defaults to zero; free-T0 support is a sampled relaxation, not a candidate for all-energy unitarity in this family. Zero is not sufficient: the subsequent 150-energy candidate still fails independent physical evaluation. This derivation does not establish the authors' finite prescription or transfer to a different high-energy density family.

A finite density bound already makes the sampled pure/chiral amplitude problem compact with free T0: the two imaginary S-wave disks bound sigma2 and then sigma1, because the nodal jump map is the identity and the bounded double densities contribute bounded rows; one real S-wave disk then bounds T0. All remaining constraints are closed, and the zero amplitude is feasible. This establishes existence of finite support optima without discarding the constant; it does not establish coupled gauge feasibility.

`SineSourceRows(M,L,bits,order)` implements this family, `convert_subtraction` implements the triangular map, and `prepare` emits the complete operator plus the above physical contract. The infinity row is a diagnostic row; producing an operator imposes no constraints by itself. For this sine prescription, `evaluate` reads full C_flat coefficients with their declared identity and evaluates partial waves at arbitrary positive s, including the threshold. Physical upper/lower lips follow the same conformal map.

The [M50/L10 run](../../../results/runs/stage_A_M50_L10_20260906/report.json) generated all 3011×3876 entries. Angular integrals use Arb arithmetic and Gauss panels placed near the crossed-cut singularities. Four representative rows compared order24 with order32: max absolute change 2.01e−40, max relative change 6.76e−27, including the first-node spin19 wave. This is an empirical quadrature comparison. The NPZ matrix is rounded to float64 for optimization; it has neither that many stored significant digits nor certified integral radii. Independent tests integrate the full crossed amplitude directly and verify the Cauchy transform and coordinate change. No continuum unitarity statement follows.

## Physical question and source hierarchy

For B3, v3 §4.2 and Fig.4's caption explicitly display only x=f00(3)>0. Region aggregation therefore declares `published_positive_x` and intersects the outer support half-planes with x>=0. Its inner polygon uses retained feasible points, their convex mixtures at x>=0, and the exact zero amplitude. Mixture weights are rounded toward the positive endpoint and their coordinate balls are checked. Negative-x coefficient files remain available; the underlying support problems and unitarity constraints are unchanged. The generic full-plane geometry remains available, and B2 uses it. The .01 metric distance budget is unchanged; a passed positive-window calculation is never called a result for the entire chiral plane.

The paper combines IR chiral symmetry breaking, pion S-matrix/form-factor positivity, and UV pQCD finite-energy SVZ sum rules for Nc=3, Nf=2. Its central comparison is chiral-only S0/S2 agreement with inadequate P1 (Fig.7), followed by improved P1 and a rho resonance after UV information is coupled (Fig.8–10). It neither determines a unique QCD amplitude nor performs a precision fit to experiment.

Source order: v3 equations and figures → explicitly identified author implementation → declared local finite prescription. Later author code is evidence about possible numerical choices, not proof of the prescription that generated v3's figures. Source anchors below refer to [v3 TeX](../../../references/2309.12402v3-source/prd_submission_2.tex); no plot is used to tune missing density bounds or UV inputs.

## Fixed scattering conventions (A1/A3)

Use pion-mass units, s+t+u=4, t=−(s−4)(1−μ)/2, and

\[
T^0=3A(s,t,u)+A(t,s,u)+A(u,t,s),\quad
T^1=A(t,s,u)-A(u,t,s),\quad T^2=A(t,s,u)+A(u,t,s),
\]
\[
f^I_\ell(s)=\tfrac14\int_{-1}^{1}P_\ell(\mu)T^I(s,t,u)d\mu,
\qquad S^I_\ell=1+i\pi\sqrt{1-4/s}\,f^I_\ell=\eta e^{2i\delta}.
\]

These are (2.7–12), TeX lines 554–579. Only even spins appear for I=0,2 and odd spins for I=1. L=10 means 0,2,…,18 in each even channel and 1,3,…,19 in I=1, not all spins through 10. The original M=50 energy grid gives 1500 retained scattering samples; the sine provider allows a separate scattering grid, so 150 energies with the same M50 basis give 4500 samples. The current PV provider keeps its original physical nodes. There is one constant, two M-vectors, an unrestricted M² double density, and a symmetric M² double density: 3876 real amplitude coefficients before the optional zero equality is imposed. A gauge discretization adds two ImF M-vectors and two current-density M-vectors, for 4076 raw real variables before eliminations.

The full scattering operator obeys S_full^dagger S_full=1. Let P project onto
the normalized two-pion elastic states in a fixed I,ell sector. Then the
elastic block satisfies

\[
(PS_{\rm full}P)^\dagger(PS_{\rm full}P)
=P-PS_{\rm full}^\dagger(1-P)S_{\rm full}P\preceq P.
\]

In the one-channel pion sector this is |S_ell^I|²<=1. The missing probability
is the sum/integral of transitions from that incoming channel to other
states. Eta=1 requires no such transition at that energy and partial wave;
it is not a condition imposed at every energy in v3. For every physical
s>=4 and every allowed I,ell, write kappa=pi sqrt(1-4/s). The concrete
constraint is Eq.(2.12):

\[
|S|\le1\quad\Longleftrightarrow\quad
\Delta=1-|S|^2=2\kappa\operatorname{Im}f-\kappa^2|f|^2\ge0
\quad\Longleftrightarrow\quad
\begin{pmatrix}1&S\\S^*&1\end{pmatrix}\succeq0.
\]

The Gram matrix is Eq.(3.71), PDF p.24, TeX 1009–1017; its eigenvalues are 1±|S| and determinant is Delta. The normalization and unitarity equations (2.9–12) are on PDF p.9. For s>4, this is equivalent to Im f>=kappa|f|²/2; Im f>=0 alone is insufficient. Eta<1 is allowed, eta=1 is an elastic limiting case. At exact threshold s=4 and finite f, kappa=0, S=1 and Delta=0; the Gram is positive semidefinite, not positive definite. Subthreshold amplitudes do not obey this physical disk condition.

`model.scattering_matrix`, `unitarity_margin(f,s)` and `scattering_gram(f,s)` share these conventions. The margin uses kappa(2 Im f−kappa|f|²) to avoid subtracting two numbers near one. Independent analytic tests start from S=eta exp(2i delta), with eta=1,.6,1.1, and check the margin and Gram determinant; the eta>1 example has Im f>0 and is correctly rejected. These physical helpers reject nonfinite f and energies below, or uncertain across, threshold; the general amplitude evaluator still accepts subthreshold energies. A returned sign is about the listed energy/wave pairs only. Arb rounding enclosures for the evaluated quadrature do not include angular-integration error, and neither a sampled pass nor additional finite spins proves continuum or full S-matrix unitarity.

Amplitude density and current spectral density are different variables. In the current flat coefficient representation `C_flat_ii=rho2_ii`, `C_flat_ij=2*rho2_ij` off diagonal. The density norm must be evaluated on the actual amplitude densities, never on quotient coordinates.

The target plane is x=f00(3), y=f11(3), below the physical threshold s=4. Its coordinates are not cross sections. A boundary projection need not identify a unique full amplitude. Save the amplitude selected by each optimization, including any selection tie-break.

## Parameters directly stated in v3 (A2)

| Input | v3 value and source |
|---|---|
| Mass and reference scale | mπ=140 MeV; fπ≈92 MeV, (2.2), TeX 507–513 |
| Matching scale | sqrt(s0)=1.2 GeV, hence dimensionless s0=(1.2/.140)²=3600/49; some prose incorrectly writes s0 with energy units |
| UV theory | Nc=3, Nf=2, αs=.4; mu=4 MeV, md=7.3 MeV, (2.54), TeX 865–869 |
| Condensates | gluon .023 GeV⁴; scalar −(.1 GeV)⁴, (2.55), TeX 870–874 |
| Chiral samples | s=1/2,1,3/2,2, (3.63–64), TeX 938–953 |
| Chiral tolerances | .006,.004,.002,.001,.0006,.0002 in Fig.4; .002 for subsequent plots, TeX 1079–1083 |
| Current channels | S0 and P1; S2 retains ordinary scattering unitarity, (3.70–71) |
| FESR moments | S0: n=0,1; P1: n=−1,0; four moments, (2.50)/(2.56) |
| UV tolerances | εSR=.002, εFF=.00006, §4.3, TeX 1141; SR scaling remains an open source question |
| Resolution study | (M,L)=(50,8),(50,10),(50,12),(45,10),(60,10), Appendix A; five configurations |

The chiral linear residuals are

\[
r_{01}(s)=f^0_0(s)-\frac{3(2s-1)}{s-4}f^1_1(s),\qquad
r_{21}(s)=f^2_0(s)-\frac{3(2-s)}{s-4}f^1_1(s).
\]

v3 says “some norm”; two separate Euclidean four-vector balls are our explicit hypothesis. The ratios remove fπ. The black reference point uses fπ separately: x≈.07332, y≈−.004888 and y=−x/15. Four residual constraints do not impose exact linear functions or chiral scattering lengths. `boundary --mode chiral` applies the selected ε and `--mode pure` removes these balls for either declared source. In the strengthened sine branch, barrier uses FG separation/cuts and native Clarabel can lift the global polynomial cone to PSD; zero-infinity also imposes the five tails. Sampled/free does not impose these conditions. Two L2 balls are not equivalent to the later combined eight-vector norm. The B finite regions and C method reproduction are complete. The regulator sensitivity and original-setting identity questions remain documented, not promoted into new prerequisites.

## Current / UV contract (A2, then D1–D3)

The scalar normalization is F0(0)≈mπ²=1, whereas the vector charge gives F1(0)=1. Use the paper's subtracted dispersion relation (2.37), and ReF_i=1+K_ij ImF_j from (3.66–67). The latter explicitly defines a form-factor Hilbert kernel; it does not specify the entire scattering double-density implementation.

Write mathcal F=k F. Equation (2.33) gives

\[
k_0^2=\frac{3\sqrt{1-4/s}}{256\pi^5},\qquad
k_1^2=\frac{(s-4)\sqrt{1-4/s}}{384\pi^5}.
\]

For each retained node in S0 and P1, impose (3.70)

\[
\begin{pmatrix}1&S&\mathcal F\\S^*&1&\mathcal F^*\\
\mathcal F^*&\mathcal F&\rho_{\rm current}\end{pmatrix}\succeq0.
\]

`model.current_gram` uses the congruent matrix with F and rho_current/k²; k>0 makes the two PSD conditions equivalent. This congruence applies strictly above threshold; k=0 at threshold cannot be divided out. Its S must come from the same amplitude variables as all the scattering constraints. Expected two-pion saturation and Watson's theorem explain the physics; do not add eta=1 or rho_current=|mathcal F|² as strict extra constraints.

For normalized moments J_n=s0^(−n−2)∫₄^s0 rho_current(x)x^n dx, the printed (2.56) gives approximately

| Channel, n | J_n from printed coefficients |
|---|---|
| S0, 0 | 4.41870e−7 |
| S0, 1 | 2.82014e−7 |
| P1, −1 | 5.75484e−5 |
| P1, 0 | 2.69948e−5 |

These numbers are rounded paper inputs, not exact UV identities. A comparison branch should begin with the printed moments; independently recompute them from (2.50) to audit the UV conventions. The paper starts with equal mq but reports unequal mu,md; selecting mq=(mu+md)/2 is an additional choice, not an explicitly recovered rule. With the printed masses it produces scalar moments 8.484% below the printed values. Replacing Nf mq² by mu²+md² (the RMS choice) reduces the difference to about .680%. The latter suggests a possible flavor convention, but does not prove the original prescription. Both are explicit audit branches; the scalar FF cap follows the chosen mq rule.

The FESR Jacobian is (π/M)s'_i, while a scattering Cauchy kernel already includes 1/π and therefore uses s'_i/M. Confusing these weights changes the moments by π. Equation (3.72) integrates only to s0 but writes a full-M sum: a finite implementation must specify the low-energy mask and the partial endpoint cell. Our historical current witness used a hard mask (43 low, 7 high nodes at M50). This is a declared numerical choice. Equation (3.74) writes raw integral errors while (2.56) reports normalized moments, leaving the scale and grouping of εSR=.002 unresolved. Raw per-moment absolute errors, normalized errors and a combined vector norm are different feasible sets. A2 fixes raw per-moment errors as the main conditional reading and keeps normalized absolute errors as a named comparison. In normalized units the four raw .002 tolerances become approximately (3.70525e−7,5.04325e−9,2.72222e−5,3.70525e−7).

For s_i>s0, (3.75) bounds the **squared rescaled** form factors:

\[
|\mathcal F_0|^2\lesssim2m_q^2\epsilon^{FF},\qquad
|\mathcal F_1|^2\lesssim\epsilon^{FF}/2.
\]

It is not a bound of the same numerical size on |F| or |F|². With the conditional arithmetic-mean mq=5.65/140 and epsilonFF=6e−5, the caps on |mathcal F0|² and |mathcal F1|² are 1.9544387755e−7 and 3e−5; the RMS branch changes the scalar cap explicitly. Four FESR integrals and these high-energy inequalities do not insert a rho resonance by hand. Active code contains these UV helpers plus `current_operators`, re-exported by `operators`, and the `prepare-current` entry. The affine layout has 4076 variables, 100 real 3-by-3 PSD blocks, four moments and 14 high-FF cones at M50; independent algebraic controls have passed. D has completed this M50 current preparation and a fully coupled feasible solution, as recorded above. The [independent current witness](../../../results/source_M50_exact_current_fiber_256.json) remains historical evidence; it is not used as a substitute for the joint solution.

The units are explicit in `fesr_data`: physical scalar and vector spectral densities are divided by mπ⁴ and mπ², respectively; raw moments in the order (S0,0),(S0,1),(P1,−1),(P1,0) are divided by mπ^(6,8,2,4). F0_phys is divided by mπ², while F1 is dimensionless. Direct integration of the complex OPE around the FESR circle reproduces the Eq.(2.50) implementation independently.

For rho_current=1, exact raw integrals are (s0−4,(s0²−16)/2,log(s0/4),s0−4). At M50, hard-midpoint relative errors are approximately (14.64%,28.86%,4.57%,14.64%); clipped-phi reduces them to (.917%,2.855%,.106%,.917%). The latter uses the overlap of each phi cell with [0,phi0], keeping its original node and Jacobian; in some M it uses a partial cell whose center is above s0. The high-FF mask always remains node>s0. This example demonstrates a real discretization sensitivity, not an estimate for an optimized unknown spectrum. Later D/F calculations must compare the declared rules.

## Distinct finite scattering prescriptions and the completed rank diagnosis

The active `PVSourceRows` follows the audited later-author PV/midpoint and analytic Legendre-Q collocation kernels with the target paper's center nu0=0. It retains all coefficients and has defined native physical values, subthreshold values and the threshold limit. Other physical points, extra-row caches and analytic `embed` are rejected. Its coordinate label is `unsubtracted PV-midpoint C_flat`; sine retains `unsubtracted sine-cardinal C_flat`. The [kernel mapping](../../../results/runs/stage_B1_source_method_audit_20260907/COLLOCATION_KERNEL_MAP_ZH.md) establishes the local formula correspondence, not the identity of all settings that produced v3's figures.

Source and constraint scope are independent declarations. PV solving is restricted to sampled/free, with sine-tail and FG diagnostics marked `not_applicable`. Sine sampled/free reports these diagnostics without enforcing them; explicit sine sampled/zero retains the five tails. Sine strengthened enforces global FG and, when zero is selected, the five tails. The old working-set/constraint-generation dispatcher is retired; active solvers use all supplied rows. Saved coefficients/duals carry their source identity; cross-source coefficients may initialize a new full-H check, but do not transfer the analytic amplitude, feasibility, old dual bounds or resampling proof. `regions` rejects mixed source/scope/B/input identities.

Original-center midpoint nodes and subtracted off-cut kernels are

\[
\phi_j=\pi(j+\tfrac12)/M,\quad x_j=4/\cos^2(\phi_j/2),\quad
w_j=x'_j/M,\qquad q_j(v)=\frac{w_jv}{x_j(x_j-v)}.
\]

PV/midpoint combines these rational functions with assigned physical values K+iI, which are not limits of those rational functions at their poles. Finite sine densities instead have analytic modes z(v)^n with z=(2−sqrt(4−v))/(2+sqrt(4−v)). For sampled sin(nφ), the midpoint transform is exactly

\[
\frac{z^n+z^{2M-n}}{1+z^{2M}},\qquad1\le n\le M,
\]

including the doubled Nyquist numerator. It is not z^n up to a constant or an invertible finite-basis change. [Transform derivation](../../../results/evidence/74_midpoint_and_sine_cauchy_transforms.md).

The angular kernel is J_i^ell=w_i Q_ell(1+2x_i/(s−4))/(s−4) before subtraction; below threshold use its real exterior continuation. [Construction](../../../results/evidence/72_midpoint_PV_kernel_construction.md), [independent review](../../../results/evidence/73_cauchy_angular_projection_review.md). The paper's scattering representation (2.7) is unsubtracted and §3 gives nodes and variables, but does not recover our full physical/off-cut prescription. Exact algebraic changes of coordinates must be distinguished from changing the function family or adding infinity=0.

Both declared unregularized finite models have augmented rank3011 at M50/L10. The imaginary block has 1500 rows, reduced through 450 columns; the remaining block contains 1500 real scattering rows, eight chiral residuals, infinity, and two targets. An exact constant pivot reduces its 1511×2376 matrix to a 1510×2375 selection problem, retaining every row. Certification requires ||I−R0 B||∞<1, not a numerical inverse alone. PV's residual is <9.984e−159, independently checked with integer arithmetic. Fixing the 3009 nontarget rows leaves 867 coefficient directions and both targets independently variable. This proves the entire target plane only for those finite prescriptions, with their compatible fixed current witness; it proves neither the authors' implementation nor continuum unboundedness.

These are historical proofs: the quotient reconstruction is an isometry only on surviving flat coordinates, and the M4/L1 rank34 ceiling has an isospin capacity proof. The old rank and PV-support executables have been retired; the present PV provider uses the unified calculation entry, while the [pre-global-SDP core and tests](../../../results/evidence/pre_global_sdp_core_20260907.tar.gz) and original proof inputs remain. No rank calculation is a prerequisite to current region work. [Finite-sine proof](../../../results/evidence/69_full_analytic_rank_and_finite_target_plane.md), [integer audit](../../../results/evidence/70_full_rank_integer_audit.md), [small-control capacity](../../../results/evidence/75_isospin_capacity_and_small_PV_rank.md), [PV imaginary proof](../../../results/evidence/76_PV_midpoint_imaginary_rank.md).

## B1: bounded hypothesis and independent support audit

Later author code contains an ordinary vector l4 bound on the two amplitude double-density families. Its use in the target calculation is unresolved; see [source decision](../../../results/evidence/71_source_boundedness_decisions.md). In our flat coordinates,

\[
\|\rho\|_4^4=\sum_{ij}(\rho_1)_{ij}^4+
\sum_i(C_{\rm flat})_{ii}^4+\tfrac1{16}\sum_{i<j}(C_{\rm flat})_{ij}^4.
\]

At M50, later code `norm(rho/3775,4)<=100` corresponds to B=377500. This is an unweighted ordinary vector norm, not a continuum L4 integral. A finite bound on both double-density families bounds the singles through imaginary S waves and the constant through a real S wave, making this finite amplitude problem bounded. It does not prove compactness of the entire current sector.

Both active sources retain all 3775 double-density variables at M50. The code's default B377500 is the later-author recipe, not a selected stability platform or a recovered v3 regulator. The [regulator protocol](../../../results/runs/stage_B_regulator_protocol_20260907/protocol.json) fixed a logarithmic B sequence; its failed chiral stability windows are retained. The bounded B delivery now fixes B377500 for all seven regions and reports that dependence separately. Closeness to a PDF curve is not a selection rule. In particular, the current PV chiral right endpoint remains materially different from Fig.4 despite a small support gap. See STATUS for the numerical records.

The predeclared numerical requirements are project choices: <=1% support variation across two regulator decades, B1 support gap <=1e-4 in stated objective units, and regional Hausdorff distance <=.01 in the declared metric. v3 gives its finite M/L setup and empirical comparisons, not these thresholds. The regulator study remains failed and is no longer a prerequisite for the conditional B prototype; historical pointwise failures also remain failures. Regional queries use a .0025 metric budget and aggregate all valid bounds under the unchanged .01 geometry criterion. Neither a failed added requirement nor a completed conditional region establishes agreement or disagreement with an exact unpublished author calculation.

The sine 10062-row strengthened runs have stopped with their saved state and stop reasons; they are not current ongoing jobs. Historical [PV seeds](../../../results/runs/pv_regularized_seed_20260906/report.json), failures and [retired source/tests](../../../results/evidence/pre_global_sdp_core_20260907.tar.gz) remain separate evidence. `quotient` supplies exact coordinate maps, Gram/moment cones and normal duals; `imaginary` handles high-spin checks, native scattering cones and the original joint audit. No active optimization uses the retired CVXPY path.

The quartic ball has the exact lift v_i>=(rho_i/B)², ||v||2<=1. Positive
column rescaling changes units only. Any positive disk margin is an inner
approximation; radial contraction may improve feasibility but lower the
objective. Only the final independently checked coefficients supply a lower
support bound. A pure/chiral support audit does not establish a coupled
gauge bound: D/E must add the current PSD and UV contributions explicitly.

## Strengthened finite-sine implementation and retained controls

The native solver can use unsubtracted or subtracted numerical coordinates.
`operators.subtraction_rows` and `subtract_amplitude_coordinates` implement
the exact H=q+b change derived above on complete actual-density coordinates.
With `--native-coordinates subtracted`, T0_sub remains a variable; only if
`--infinity zero` is selected does a linear equality impose T0_unsub=0. Setting T0_sub=0 would change
the problem and is not done. Both density families and their norm are
unchanged. Results are converted back and audited in original C_flat
coordinates; positive numerical row/variable scales do not change kappa,
chiral tolerances or the physical constraints.

A nonzero ideal analytic witness exists inside the full finite sine family,
including M50; it does not restrict the model's variable family to rank one.
Let w_j=sin(phi_j), h=H·w=1+z, and choose T0=0,
sigma1=sigma2=lambda a w, R1=lambda b ww^T, R2=lambda b ww^T/2,
with a=1/8, b=1. All coefficients here are exact analytic definitions.
For s>4 write h_s=8/s+i v, v=4sqrt(s−4)/s, and define
J_ell=(1/4)integral P_ell h(t), U_ell=(1/4)integral P_ell h(t)h(u).
Positive Stieltjes density and the exterior Q kernel give
0<J_ell<=J0<=min(1/2,4/sqrt(s)). Explicitly,
J0=4[sqrt(s)−2−2log((sqrt(s)+2)/4)]/(s−4).
Partial fractions give, for even ell,
U_ell=2/[pi(s−4)] integral rho(x)Q_ell(1+2x/(s−4))h(−x−s+4) dx,
so 0<=U_ell<=2J_ell; odd U_ell=0 by reflection symmetry.

The unscaled seed waves are exactly

\[
\widehat f^0_\ell=(\tfrac{5a}{2}\delta_{\ell0}+9bJ_\ell)h_s
 +10aJ_\ell+\tfrac{7b}{2}U_\ell,\qquad
\widehat f^2_\ell=(a\delta_{\ell0}+3bJ_\ell)h_s+4aJ_\ell+2bU_\ell,
\qquad \widehat f^1_\ell=b h_sJ_\ell.
\]

Their imaginary parts are strictly positive on every s>4 and allowed spin.
Using kappa J_ell/v<=pi, kappa|h_s|²/v<=2pi and
|D h_s+R|²<=2D²|h_s|²+2R² gives the uniform ratio bounds

\[
\frac{\kappa|\widehat f^0_\ell|^2}{2\operatorname{Im}\widehat f^0_\ell}
\le\pi\left[5a+9b+\frac{(10a+7b)^2}{9b}\right]=\frac{275\pi}{16},
\quad
\frac{\kappa|\widehat f^2_\ell|^2}{2\operatorname{Im}\widehat f^2_\ell}\le10\pi,
\quad
\frac{\kappa|\widehat f^1_\ell|^2}{2\operatorname{Im}\widehat f^1_\ell}\le2\pi.
\]

Thus f=lambda f_hat satisfies every physical scattering disk strictly for
0<lambda<=.01: its margin is at least
2 kappa lambda Im(f_hat)[1−lambda(275pi/16)]>0; at lambda=.01 the
bracket is 1−11pi/64. The bound 275pi/16≈53.9961 is independent of both
s and ell. Threshold has equality by the finite-f limit. There is no
uniformly positive absolute margin as threshold or asymptotic limits are
approached. If prescribed chiral or density bounds require it, lambda can
be reduced further because those residuals/norms are homogeneous; no
current-sector feasibility is asserted.

This proves continuum **scattering** feasibility for this ideal mode-one
witness, not for rounded nodal coefficient files or a numerical optimum.
Rounded coefficients retain only their separately checked scope. The full
problem therefore has no intrinsic absence of nontrivial scattering
solutions, but this witness does not replace B1 support optimization or
Fig.3–11. The previous numerical seed failure and precision correction
remain separate implementation evidence below.

For numerical scaling, each disk can be written exactly as
\[
\left\|(\sqrt{2/t_i}\Re G_i,\sqrt{2/t_i}\Im G_i,
\Im G_i/t_i-1)\right\|_2\le\Im G_i/t_i+1.
\]
An explicitly positive margin shifts this to an inner approximation;
zero margin is the original disk. A full QR change uses all constraint
columns and the finite density scale, without dropping singular directions.
The density norm has the equivalent power-cone form
\(t_j\ge |\rho_j/B|^4\), \(\sum t_j\le1\).

With arbitrary support weights kR,kI, chiral vectors y and tail weights a,b, form the residual
of the requested objective in the original coefficient coordinates.
In the free branch, one real S0 row removes its T0 component; the two imaginary S-wave
rows at each node remove the remaining single-density components. Add
these correction weights back into kR,kI. Only the actual-density residual
remains, giving
\[
h\le\sum_i(\sqrt{kR_i^2+kI_i^2}+kI_i)
+\epsilon\sum_g\|y_g\|_2+B\|r_{\rho,\mathrm{actual}}\|_{4/3}
+\sum_j h_{\mathrm{tail}}(a_j,b_j).
\]
For x=P_j z, y=Q_j z and x²<=2y, the finite parabolic support is
h_tail(a,b)=−a²/(2b) when b<0, and zero when a=b=0. Other weights have
unbounded support and are rejected. These weighted tail rows are removed
from the objective residual before the remaining free-direction correction.
For the linear rows P=0 the implementation uses a=0, b<=0 and hence zero
tail support. With the native global cone, the original-coordinate residual
also includes +nu D/B, where D maps complete C_flat coefficients to
Chebyshev-T coefficients of p. Arb verifies the moment Toeplitz matrix
T(nu)_ij=nu_|i-j| as PSD; a corrected nu is used in the residual again,
never silently substituted into an old bound. Internal p/unit scaling is
converted back to p/B for this audit. For pure mode omit the chiral term. Negative scattering kI
uses kR²/(sqrt(kR²+kI²)−kI) to avoid cancellation.
C_flat off-diagonal residuals receive a factor2 in the actual-density dual
norm. In the zero branch T0 has no free residual to correct. The correction is algebraic; a small floating residual in an
unbounded free direction cannot just be ignored. This check can certify
bounds for the saved float64 operator, with angular quadrature uncertainty
reported separately.

The saved M50 sampled-support failures and their independent physical violations
remain linked in [STATUS.md](../../../STATUS.md). They motivated the two distinct
asymptotic conditions below; no old feasibility conclusion is transferred.

### Necessary fixed-spin tail constraints in the finite sine family

These are consequences of the declared finite family and Eq.(2.12), not
additional parameters read from a plot or recovered from the authors' code.
Assume fixed finite M, real amplitude coefficients, and T0^u=0. Define

\[
d_j=4\sum_{n=1}^M(-1)^{n+1}nW_{nj}
=\frac{(-1)^{M+j}x_j}{2M},\qquad j=1,\ldots,M,
\qquad e_j=8\sum_{n=1}^M(-1)^{n+1}n^2W_{nj}.
\]

The closed expression includes the Nyquist mode with its correct inverse
sine-transform weight. Expanding the same analytic H on its physical lip
and crossed ray gives

\[
H_j(s+i0)=\frac{i d_j}{\sqrt{s}}+\frac{e_j}{s}+O(s^{-3/2}),\qquad
H_j(-X)=\frac{d_j}{\sqrt X}-\frac{e_j}{X}+O(X^{-3/2}).
\]

For fixed ell, set a=(1−mu)/2 and t=−(s−4)a. The angular kernel obeys

\[
J_{\ell,j}=\frac14\int_{-1}^1P_\ell(\mu)H_j(t)\,d\mu
=\frac{d_j}{(2\ell+1)\sqrt s}+O\!\left(\frac{\log s}{s}\right),
\]

because integral_0^1 P_ell(1−2a)/sqrt(a) da=2/(2ell+1). This expansion
does not apply pointwise at a=0: the endpoint layer a=O(1/s), where H is
bounded, contributes O(1/s). Away from it, the O(1/(sa)) remainder integrates
to O(log(s)/s). Thus endpoints do not change the displayed 1/sqrt(s) term.
The crossed product H_i(t)H_j(u) is real and its full angular integral is
O(1/s); its endpoint layers contribute O(s^(−3/2)). It cannot supply the
leading imaginary term below.

Use the unsubtracted Eq.(2.7) single-density coefficients; another coordinate
system must first apply the exact subtraction conversion. Write R1,R2 as
actual density matrices, retaining the unrestricted R1, and set

\[
p_a=\sigma_a^Td,\quad r_a=d^TR_ad,\qquad
(q_0,q_1,q_2)=(4r_1+r_2,r_1-r_2,r_1+r_2),
\]
\[
(b_0,b_1,b_2)=(p_1+4p_2,p_1-p_2,p_1+p_2).
\]

Only terms containing physical H(s) can be imaginary. Isospin crossing
groups their double-density matrices as 3R1+R1^T+R2 in I=0, R1^T−R2 in
I=1, and R1^T+R2 in I=2. Projection of H(t)±H(u) supplies 2J_ell for
the allowed parity. Consequently

\[
\operatorname{Im}f^0_0(s)=\frac{3p_1+2p_2}{2\sqrt s}
 +\frac{2q_0}{s}+O\!\left(\frac{\log s}{s^{3/2}}\right),\qquad
\operatorname{Im}f^2_0(s)=\frac{p_2}{\sqrt s}
 +\frac{2q_2}{s}+O\!\left(\frac{\log s}{s^{3/2}}\right).
\]

Thus 3p1+2p2>=0 and p2>=0 are two necessary linear S-wave tail conditions.
For every fixed allowed ell>0, the direct single-density term projects to zero,
but its crossed part is real and cannot be discarded:

\[
\operatorname{Re}f^I_\ell(s)=\frac{2b_I}{(2\ell+1)\sqrt s}
 +O\!\left(\frac{\log s}{s}\right),\qquad
\operatorname{Im}f^I_\ell(s)=\frac{2q_I}{(2\ell+1)s}
 +O\!\left(\frac{\log s}{s^{3/2}}\right).
\]
\[
\Delta^I_\ell(s)=\frac{4\pi}{(2\ell+1)^2s}
 \left[(2\ell+1)q_I-\pi b_I^2\right]
 +O\!\left(\frac{\log s}{s^{3/2}}\right).
\]

The strongest leading conditions among fixed nonzero allowed spins are

\[
5q_0\ge\pi b_0^2\quad(\ell=2),\qquad
3q_1\ge\pi b_1^2\quad(\ell=1),\qquad
5q_2\ge\pi b_2^2\quad(\ell=2).
\]

They are three convex second-order-cone inequalities, since each q and b is
linear in the original amplitude coefficients. Their linear consequences
q1>=0 and q2>=0 imply q0=(3q1+5q2)/2>=0. Merely adding these linear
consequences misses the Re(f)^2 term: positive leading Im f is not sufficient.
For any larger fixed allowed ell, the displayed leading inequality follows
from the respective lowest-spin one; this does not justify a uniform limit
over infinitely many spins or control finite-energy behavior.

If a scalar leading coefficient vanishes, its next order must also be
examined: 3p1+2p2=0 requires q0>=pi b0², and p2=0 requires q2>=pi b2².
These conditional S-wave requirements must not be imposed unconditionally
when the corresponding 1/sqrt(s) coefficient is positive. Vanishing leading
gaps can require further orders. Therefore the two linear conditions and
three cones are necessary cuts, not a sufficient continuum certificate.

For the saved M50 candidates, direct 384-bit Arb contractions of the JSON
coefficients and exact d isolate these tails without any angular quadrature.
The d formula agrees with its independent finite-sum form to below 6.4e−111.
R2 uses C_flat/2 off diagonal, hence r2=sum_i d_i² C_ii+sum_(i<j) d_i d_j C_ij.

| Asymptotic quantity | [Chiral coefficients](../../../results/runs/stage_B_line_M50_20260906/coefficients.json) | [Pure coefficients](../../../results/runs/stage_B_pure_line_20260906/coefficients.json) |
|---|---:|---:|
| sqrt(s) Im f00 limit | 774.468020 | 767.309023 |
| sqrt(s) Im f20 limit | 876.457951 | 750.418408 |
| lim s Im f^0_2, I=0 and ell=2 | −162705.049524 | −171464.259470 |
| lim s Im f^1_1, I=1 and ell=1 | 404559.055165 | 324525.879939 |
| lim s Im f^2_2, I=2 and ell=2 | −210723.279669 | −185415.020566 |
| 5q0−pi b0² | −39163459.971300 | −30661963.603075 |
| 3q1−pi b1² | −981747.644532 | −256056.964095 |
| 5q2−pi b2² | −4687433.430296 | −4140297.254801 |

Both candidates pass the two displayed S-wave leading signs but violate all three cone
conditions. Even P1, with positive leading Im f, has a negative leading
unitarity margin. For these fixed saved amplitudes this proves eventual
high-energy violation; it does not locate its first onset or certify another
candidate. These five necessary conditions are now implemented in both the
barrier and the independent primal/outer audit. The common solver form is
2Q_a z−(P_a z)^2>=0: three rows encode (2ell+1)q_I−pi b_I², and P_a=0
for the two linear S-wave rows. `operators.asymptotic_margins` separately
checks the analytic coefficient formulas, rather than treating a pass for
rounded tail maps as the full analytic test.

An independent control has the exact optimum 1/sqrt(2): a disk together with
v²<=2rho and ||rho||4<=1/4. The direct outer bound encloses this value,
detects a disk-feasible candidate that violates the parabola, and rejects
tail dual weights with unbounded support. The barrier and its saved duals
also recover a bracket of width below 1e−5. The complete
[M3/L2, ten-energy CLI control](../../../results/runs/stage_B_tail_control_20260906/report.json)
gives lower .05739694297648645 and upper .05739759949863588, a gap about
6.56522e−7, with all five analytic tail conditions passed. This validates
the finite method; it neither reproduces the M50 problem nor establishes
continuum unitarity; independent mathematical tests retain this control.

The [new M50/L12 preparation](../../../results/runs/stage_B_unitarity_extended_20260906/report.json)
has 167 scattering energies, 6012 disk constraints and a 12035×3876 matrix.
It retains the M50 density family and is an additional coverage tool, not a
substitute for the derived necessary tail conditions. The separate old
[150-energy precision run](../../../results/runs/stage_B_extended_precision_20260906/report.json)
has ended at lower x≈.0735049050 with outer≈2.51360e12, before the new
tail conditions were applied. It has no useful support gap or accepted
physical result.

The earlier full-size [chiral](../../../results/runs/stage_B_tail_chiral_M50_20260906/report.json)
and [pure](../../../results/runs/stage_B_tail_pure_M50_20260906/report.json) attempts
stopped after a nonpositive initializer scale led to invalid log2 and NaN;
their failures remain preserved. The fail-fast
[diagnostic](../../../results/runs/stage_B_seed_representation_20260906/progress.jsonl)
located s≈4.000109664275457, I=2, ell=22, with computed
Im G≈−5.93257e−120. The cause is now identified as angular-projection
cancellation at 384-bit precision. For h(t)=1+z(t), tau=s−4→0+,

\[
J_\ell(s)=\frac14\int_{-1}^1P_\ell(\mu)h(t)\,d\mu
\sim\frac{(\tau/16)^\ell}{2(\ell+1)(2\ell+1)}.
\]

The coefficient follows by expanding h(−tau a) in Catalan coefficients and
using integral_0^1 P_ell(1−2a)a^ell da=(−1)^ell(ell!)²/(2ell+1)!.
At that first near-threshold point J22 and J23 are about 2^(−388.4) and
2^(−405.7). Thus 384-bit summation cannot preserve these signs reliably.
768/1024-bit projections agree with an independent analytic evaluation and
are positive. Increasing Gauss order alone does not repair arithmetic
cancellation.

`model.projection_precision` estimates this loss from tau and the maximum
spin, reserves 128 target bits plus 64 guard bits, rounds upward in 64-bit
blocks and respects any larger requested precision. Exact threshold is
handled by its known limit. Basis matrices and caches are rebuilt when
working precision changes; an old 384-bit cache cannot become a higher
precision operator merely by changing the arithmetic context. This is a
numerical precision budget based on the analytic threshold scale, not a
rigorous quadrature-error enclosure for every possible coefficient vector.

The [resolved M50/L12 preparation](../../../results/runs/stage_B_unitarity_resolved_20260906/report.json)
uses 640-bit arithmetic and order32 for the complete 12035×3876 operator
at 167 energies. Earlier sampled candidates remain historical controls;
the accepted next initializer and its limited scope are stated below.

`embed` changes the sine basis without changing the exact amplitude, unless
an optional interior is supplied for a separately reported convex FG repair.
Neither operation transfers norm bounds, feasibility or optimality; saved
coefficients require new checks in the target problem.

Earlier M50 sampled candidates and normal-LP failures remain in
[STATUS.md](../../../STATUS.md); in particular, the old radial candidates fail the
fixed-energy large-spin test below. NumericalError returns, including zero
vectors, are not new results. Recovery retains a better certified incumbent,
uses genuine coefficient scaling or convex mixing, records objective loss,
and checks the final saved coefficients with the original operator and FG
polynomials. It never clips S or polynomial values.

## Fixed-energy large-spin necessary condition (B1)

This derivation and its global-polynomial tests concern the finite-sine family. They are not applied to PV collocation; that provider reports `not_applicable`.

This is the other unitarity limit: fix a finite s>4 and take allowed
ell→infinity, instead of fixing ell and taking s→infinity. It follows from
the same Mandelstam representation (2.7), projection (2.9), and unitarity
(2.12), specialized to the declared finite sine family. It is not an
assumption of componentwise double-density positivity and is not claimed
to have been implemented in the unpublished algorithm behind v3.

For the density basis r_j(x)=sum_n W_nj sin(n phi(x)), define

\[
a_j=\sum_{n=1}^M nW_{nj}
=\frac{(-1)^{j+1}}{2M\sin^2(\phi_j/2)},\qquad
r_j(x)=a_j\sqrt{x-4}+O((x-4)^{3/2}),\quad
v_j(s)=\operatorname{Im}H_j(s+i0).
\]

This threshold vector a differs from the high-energy vector d used above.
With actual R2 entries and the unrestricted R1, set

\[
C_0=3R_1+R_1^T+R_2,\qquad C_1=R_1^T-R_2,\qquad
C_2=R_1^T+R_2,\qquad \alpha_I(s)=v(s)^TC_Ia.
\]

At nonzero allowed ell, direct single-density imaginary terms project out;
the entirely crossed products are real. Hence Im f_ell^I=2v^TC_I J_ell
exactly, with

\[
J_{\ell,j}=\frac1{\pi(s-4)}\int_4^\infty r_j(x)
Q_\ell\!\left(1+\frac{2x}{s-4}\right)dx.
\]

The factor 1/[pi(s−4)] includes the paper's 1/4 angular normalization.
Let N=ell+1/2, beta=sqrt(1−4/s) and eta0=arcosh(1+8/(s−4)). The positive
Q kernel is Laplace-dominated by its endpoint x=4. Its N^(−1/2) factor
together with the sqrt(x−4) endpoint integral gives

\[
J_{\ell,j}=N^{-2}e^{-N\eta_0}
\left[\frac{a_j}{2\beta}+O(N^{-1})\right],\qquad
\operatorname{Im}f_\ell^I=N^{-2}e^{-N\eta_0}
\left[\frac{\alpha_I(s)}{\beta}+O(N^{-1})\right].
\]

The real partial wave has the same or faster exponential decay, so its
square decays faster than a nonzero displayed imaginary term. Consequently

\[
1-|S_\ell^I|^2=N^{-2}e^{-N\eta_0}
\left[2\pi\alpha_I(s)+O(N^{-1})\right].
\]

Thus alpha_I(s)>=0 is necessary for every fixed s>4. A negative alpha
proves eventual violation along the allowed spin parity; zero leading
coefficient requires the next endpoint order. This does not fix the first
violating spin or justify an exchange of the large-s and large-ell limits.

Writing c=cos(phi), sin(phi)>0 on the physical open cut gives

\[
\alpha_I(s)=\sin\phi\,p_I(c),\qquad
p_I(c)=\sum_{n=1}^M(WC_Ia)_n U_{n-1}(c),\qquad
s=\frac8{1+c}.
\]

Each p_I is a degree-at-most M−1 Chebyshev-U polynomial, linear in the
amplitude double-density coefficients. Therefore p_I(c)>=0 throughout
[−1,1] is an equivalent necessary sign condition, including the endpoint
limits by polynomial continuity. A chosen c gives a linear separating
constraint on the original coefficients; it does not require every rho1 or
rho2 component to be positive. `imaginary.high_spin_diagnostic` now isolates
all derivative roots with Arb/Acb intervals and evaluates the endpoint and
stationary intervals to bound the global minimum. A conservative coefficient
bound remains available if isolation is incomplete; inconclusive bounds
are not reported as positivity. Negative witnesses and a certified global
lower bound are distinct outputs, both independent of angular quadrature.

For degree at most 49 (M50), the native full-interval representation is

\[
p_I(c)=(1+c)v(c)^T G_{I,+}v(c)+(1-c)v(c)^T G_{I,-}v(c),\qquad
v=(T_0(c),\ldots,T_{24}(c))^T,\qquad G_{I,\pm}\succeq0.
\]

Thus the three channels use six 25-by-25 Gram matrices. This is the standard
weighted-SOS equivalence, not a density truncation; see
[Roh–Vandenberghe (2006), §5.1, Eq.(45)](../../../references/roh_vandenberghe_2006.pdf)
(and its trigonometric/cosine forms in §§3–4). The Gram coefficient map is
exact in the chosen polynomial basis. Solver residuals are still numerical:
primal polynomials and the moment-dual PSD/residual are independently checked.

The following witnesses use the saved
[chiral radial coefficients](../../../results/runs/stage_B_radial_chiral_M50_20260906/coefficients.json)
and [pure radial coefficients](../../../results/runs/stage_B_radial_pure_M50_20260906/coefficients.json).
Floating stationary points and endpoints were searched, then 512-bit Arb
confirmed the negative polynomial values at those finite physical points.
All six displayed witnesses lie within the paper's E<=1.2 GeV window.

| Candidate | I | cos(phi), approximately | E (GeV) | p_I(c), approximately |
|---|---:|---:|---:|---:|
| chiral | 0 | −.879938878892578 | 1.14280421 | −1595764.806897 |
| chiral | 1 | −.709800696708287 | .73506348 | −226997.492870 |
| chiral | 2 | −.891111111111111 | 1.2 | −725520.511779 |
| pure | 0 | −.485140948611078 | .55185992 | −860991.988070 |
| pure | 1 | −.891111111111111 | 1.2 | −612984.653323 |
| pure | 2 | −.491825045216312 | .55547741 | −442216.570207 |

Thus both candidates have a necessary infinite-spin obstruction even though
they passed the five fixed-spin high-energy conditions. They have been
excluded from physical acceptance. The observed dense-audit violations at
ell=24–27 motivate this test, but a particular finite-spin violation should
only be attributed to its leading asymptotic term after checking that limit
at the same energy. The p_I diagnostic and FG linear cuts are now active:
`analysis.high_spin_rows` builds the coefficient rows, and retained cuts are
passed through `--additional-constraints`. Imposing finitely many cuts does
not by itself certify the full interval; that is checked separately on the
returned polynomial. This necessary-condition gate is not full unitarity.

Active `evaluate` checks independent energy/spin grids and the global
high-spin polynomials; `violations.json` and cumulative `constraints.json`
record finite-wave violations and FG cuts. `boundary` and separate `dual`
consume the retained rows, with original-coordinate support audits.
`__init__` contains the public Model/UVConfig types, still importable from
`model`; the ten-file structure and physical conventions are unchanged.

The strengthened-sine [M50 initializer v2](../../../results/runs/stage_B_FG_repaired_M50_v2_20260907/report.json)
has target (.080550644094,−.005276824472). Its
[independent 832-bit/order40 audit](../../../results/runs/stage_B_FG_repaired_M50_audit_20260907/report.json)
passes 707×42 checks, with global p_I lower bounds approximately
.4740803, .1563960 and .1616019. Its [full LP upper bound](../../../results/runs/stage_B_FG_primal_dual_M50_20260907/report.json)
was 21478.2379 at that stage; this record is an initializer, not a boundary.
The [767-energy/L14 preparation](../../../results/runs/stage_B_dense_M50_prepare_20260907/report.json)
generated 32214 scattering rows and a 64439×3876 matrix in about385 seconds
at832-bit/order40; this is an operator artifact, not a support result.

Equivalent subtracted native coordinates give Solved for both
[M3 qdldl](../../../results/runs/stage_B_native_subtracted_qdldl_M3_20260907/report.json)
and [M3 faer](../../../results/runs/stage_B_native_subtracted_faer_M3_20260907/report.json),
with original-coordinate gaps about1.12e−7 and1.60e−7. These controls do not
transfer to M50. The full subtracted-native
[chiral](../../../results/runs/stage_B_native_subtracted_M50_round1_20260907/report.json)
and [pure](../../../results/runs/stage_B_native_subtracted_pure_M50_round1_20260907/report.json)
runs ended with MaxTime after 13 and 12 iterations. Original-coordinate
recovery gives chiral lower .04174070348279 versus upper 289961049.6036,
and pure lower .04779318187414 versus upper 2525930.456843. These are
saved finite feasible points with unusable gaps; neither run completes B1.
The native objectives 2.21484 and 8.56481 preceded feasibility recovery
and are not valid lower bounds.

The retired strengthened-sine `--constraint-generation` dispatcher solved a subset of the declared finite physical
rows, checks the candidate against every original H row, and adds violations.
The full M50 coefficient family, density/chiral constraints and FG cone remain.
A subset dual gives an upper bound after omitted multipliers are padded by
zero; its candidate becomes a lower bound only after full-H feasibility.
This changes the solve order, not Eq.(2.12) or the reproduction criterion.
The first [native M50 exchange run](../../../results/runs/stage_B_CG_soc_M50_20260907/report.json)
ended without a useful bound. The direct-margin barrier with global FG
separation closes the [M3/L2 control](../../../results/runs/stage_B_barrier_CG_M3_control_20260907/report.json)
to gap 8.14e-5. The first full M50 chiral exchange gives an independently checked interval
[.080550643289,.092990467638], still wider than the required gap.
Subsequent finite support closure and newly observed violations are tracked in STATUS; neither a solver success flag nor a strengthened finite model identifies the published implementation.

The [reference extraction](../../../results/runs/stage_B_reference_20260906/README.md)
contains the original Fig.3 and six Fig.4 marker sets recovered from PDF
vector paths. These are reference geometry, not bootstrap results. The
Fig.4 black marker is approximately (.07131334,−.00475418), whereas the
stated 140/92 parameters give (.07332139,−.00488809). Keep both labels.
The paper's ray angle and a support normal are different parameters of
the same convex boundary; compare shapes and support values, not direction
indices.

## Equivalent coordinates and historical constraint generation (B1)

The optional `absorptive` numerical coordinates retain the complete amplitude.
After the subtraction change, let D0 and D2 be the actual-density columns of
the original native-node Im f00 and Im f20 rows. Their single-density blocks
are exactly (3I/2,I) and (0,I). Define u0=Im f00 and u2=Im f20 at those nodes.
Then the inverse triangular map is

\[
\sigma_1^s=\tfrac23(u_0-u_2)-\tfrac23(D_0-D_2)\rho,
\qquad \sigma_2^s=u_2-D_2\rho.
\]

T0_sub and every actual double-density coordinate remain independent variables;
the single-density determinant is (2/3)^M. All amplitude and tail rows,
including the T0_unsub=0 equality when declared, undergo the same substitution. Output
first reverses the numerical units, then this map, then the subtraction map.
Original-coordinate residuals include numerical errors in these operations.
No density direction or physical constraint is removed.

The two native L4 representations are also equivalent: the power-cone version
uses u_j>=|rho_j/B|^4 and sum u_j<=1. It permits both signs of rho_j, just as
the SOC representation. Neither representation changes the norm or B.

The retired strengthened-sine `--constraint-generation` path retained the original M and complete physical row
universe. A working set starts with the native nodes; all supplied physical
rows are checked and violated inactive rows are added. The barrier variant
also separates the global polynomials before recovering a feasible point:
otherwise a convex recovery could hide the very violations needing new cuts.
It adds p_I(x_j)/2>=0 at negative witnesses and independently checks the final
polynomials on the whole interval. The native variant uses the exact PSD lift.

For nonpositive multipliers b_j of the polynomial halfspaces, define

\[
\nu_{I,k}=-\frac B2\sum_{j:I_j=I} b_j T_k(x_j).
\]

This is a positive moment functional. After outward moment repair, the
independent outer calculation uses +nu D/B in the original residual, and
keeps only the original analytic tail multipliers separately. Thus neither
the FG contribution nor a rounded cut row is silently counted twice.
Each working feasible set contains the full constrained set, so its valid
dual support remains an upper bound. A lower bound still requires the full
physical-row, tail and global-polynomial checks; an iteration count or a
solver success label is insufficient.

`cg_state.npz` retains physical indices and polynomial points for continuation.
Reusing indices requires the same M and an exact ordered energy/wave prefix;
only appended physical rows are allowed. Saved duals are padded with zero on
new rows and their support is recomputed on the new H. No old feasibility or
optimality claim is transferred. The better original feasible point is
kept separately from the slightly interior reference used for numerical units.
Squares of real Arb balls use multiplication, avoiding the installed library's
NaN result for generic powers of some zero-centered nonzero-radius balls.
This is an arithmetic correction, not a change to the unitarity inequality.

## Result and verification contract (C/E/F)

Recover phases from the saved S: δ=unwrap(arg S)/2, eta=|S|, E_GeV=.140 sqrt(s); record branch choice and threshold behavior. Do not force eta=1 or interpolate assigned PV node values into an allegedly certified global amplitude. A3 supplies the sine off-node definition. Active and historical PV values are defined only on their stated domain; connecting plotted native points does not define or certify a physical off-node amplitude.

Fig.3 is a baseline, explicitly not used further in §4.1. Fig.4–7 establish the chiral comparison; Fig.8 adds UV data without changing the common scattering assumptions. Fig.9–10 require the tip plus nearby representative amplitudes, a rho crossing near 90° with the reported roughly 6% displacement, and S0/S2 checks. Fig.11 contains five resolution runs and is empirical sensitivity evidence, not a continuum theorem.

Save full coefficients/current variables, source and configuration records, individual constraint diagnostics, primal/outer support bounds and plot data under `results/runs`. Separate solver error, truncation, UV-input uncertainty and unresolved prescriptions. The paper's M/L comparisons and stability of claim-relevant observables remain required. Additional checks must respect the source's defined domain; the PV off-node limitation is not evidence of a violation and cannot be filled by a sine continuation. Rigorous interval coverage and tail estimates are additional requirements only for a continuum-certification claim. If source choices remain unresolved, final plots describe a conditional reproduction even when they resemble the paper.

For B1, the five predeclared regulator supports are summarized through `run regulators`; the published source and conclusions are in [the regulator results](../../../results/runs/stage_B_regulator_protocol_20260907/RESULTS_ZH.md). B1 retained the requested 1e-4 support gaps. B2/B3 use valid brackets and the unchanged .01 whole-region metric budget; new query objectives are scaled to unit normals in that metric with a .0025 stopping budget. Positive objective scaling changes neither the feasible amplitude set nor its exact extremizer.

## B solver correction: support-aware central stopping

A matched M50 epsilon=.0002 checkpoint exposed a concrete stopping error. The computed Newton decrement squared was 7.19e-12, but the gradient error after eliminating the free coordinates had density-support cost 3.58725. The cached-state, factorization-gradient and original-H back-transformation discrepancies were below 1e-9 in that same support scale. Evaluating the already-computed full Newton step reduced the density defect to 5.36e-5 and gave an original-H feasible support interval [.7811192483,.7852060190] at mu=1e-6. No new Newton solve or changed physical constraint was needed for this diagnostic.

Let e=c-mu*grad(phi). The nodal free-coordinate elimination writes e as disk covectors (delta kR, delta kI) plus an actual-density remainder r. Each G=kappa*f disk has radius one and the density ball has radius B. Therefore the error functional's width over the declared feasible set obeys:

W(e) <= 2 sum_j hypot(delta kR_j, delta kI_j) + 2 B ||r||_(4/3).

Ignoring the additional chiral inequalities in this bound is conservative. The support_gradient_width routine estimates this quantity in the working arithmetic; it is a numerical stopping check, not the final certificate. The near-center branch now requires both the existing small-decrement condition and W(e)<=gap/4. Otherwise it takes the already-computed full Newton step only if the existing quadratic slack update stays strictly feasible, retains that update without a fresh low-precision matrix product or radial rescaling, and recomputes the gradient on the next iteration.

The same-input full M50 check then reached [.7842651767,.7847169067], gap4.52e-4, in about122 solver seconds, stopping at mu=1e-7. Original-H Arb feasibility and support verification remained mandatory. [The matched diagnostic](../../../results/runs/stage_B1_source_method_audit_20260907/CACHED_STATE_DIAGNOSTIC.json) and [the successful run](../../../results/runs/stage_B3_pv_M50_eps0002_central_fixed_20260907/report.json) preserve the evidence. Matched first/best center records are now saved; their coefficient JSON is explicitly only an initializer. This removes an identified numerical error and does not settle the remaining differences from the paper's displayed regions.
