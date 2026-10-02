# Pre-registration: arXiv:2403.10772 v4 reproduction (frozen 2026-10-01, before any 2403 solve)

This file is committed before any SDPB solve of the 2403 model exists. Its git commit and its SHA-256 in
`/playpen1/shiqiu/gtb2403_runs/GATE_LOG.md` are the freeze record. Changes after the first solve are appended as dated
amendments with the reason; the original text is never edited.

## 0. Rules that apply to every claim

- **Inputs** are the authors' code values (spec §2, plan Global Constraints). Nothing is tuned to any published curve.
- **Points**: blue / green / red = ni 44 / 46 / 48, 𝔉₀ = ni/50 × 0.163108, maximise 𝔉₁, then one Watson step
  (factor-2 target, S2/D2/F1 targets = their own initial waves). No other point selection is used for C4–C21.
- **Comparison data** come only from the paper: Fig. 3 dots (`fig3_shape_dots_extracted.txt`, vector, ±3e-5), Figs. 4,
  5, 6, 8, 10, 11, 18, 19, 20 vector paths extracted by a script frozen before it is run on our results, and the tables.
  The paper's Fig. 4 dots sit on the same 39 nodes as ours (M=50, ν₀=−20), so phase comparisons are node by node.
- **Phase convention**: S = 1 + iΛ_ℓ h̃ at each node; δ = ½ unwrap(arg S) from threshold, η = |S|; 90° crossings by
  linear interpolation in √s between adjacent nodes. Branch choices are not tuned per solution.
- **Verification level** (docs/gtb2403/VERIFICATION_AUDIT_ZH.md): any value is reported with V1 (Arb primal); any
  statement of an upper bound, exclusion or infeasibility needs V2 on the full-space problem.
- **Tolerance for tables**: each blue and red value within max(5 %·|q_paper|, ½|q_blue − q_red|_paper).
- **Regulariser encoding**: a number obtained with a relaxation (ℓ∞, ℓ2 inner/outer) is labelled as such and never
  used alone for a verdict; verdicts use the exact ℓ4 encoding or a bracket that decides the inequality.

## 1. Claims

| ID | Quantity (paper value: blue, red unless noted) | Rule | Level |
|---|---|---|---|
| C1 | 𝔉₁^max pre-Watson at ni 44/46/48: −0.027229 / −0.028770 / −0.030397 | \|Δ𝔉₁\| ≤ 1e-4 for each; a bracket wholly outside the band is FAIL, a bracket straddling it with width > 2e-4 is INCONCLUSIVE | V1 + V2 |
| C2 | Allowed region (𝔉₀, 𝔉₁): 50 upper + 50 lower boundary dots; tip ≈ (0.1632, −0.0324) | RMS of 𝔉₁ differences over the 100 dots ≤ 3e-4; tip 𝔉₀ within 1 %; tree-level dot (0.145743, −0.0291487) inside | V1 + V2 |
| C3 | ε^χ = 2e-3 is the minimal tolerance keeping the tree-level point | tree point feasible at 2e-3; bisected exclusion threshold in [1e-3, 2e-3] | V1 (feasible), V2 (infeasible) |
| C4 | Six phase shifts and η after Watson (Fig. 4) | per-wave RMS ≤ 5° against the same-colour paper dots on nodes below 1.4 GeV for S0, P1, D0, S2, D2, F1; min η ≥ 0.98 below 1.83 GeV for every wave | V1 |
| C5 | δ¹₁ = 90° at ≈ 0.83–0.85 GeV | crossing within ±20 MeV of the same-colour paper crossing | V1 |
| C6 | δ⁰₂ = 90° at ≈ 1.35–1.40 GeV | within ±30 MeV of the same-colour paper crossing | V1 |
| C7 | Second P1 rise (ρ(1450)) | δ¹₁ at the last node below s₀ within ±20° of the same-colour paper dot | V1 |
| C8 | Table 1 (a⁰₀ 0.174/0.191, a²₀ −0.0357/−0.0403, b⁰₀ 0.285/0.282, b²₀ −0.061/−0.071, a¹₁ 28.7/30.7e-3, b¹₁ 12.4/7.28e-3, a⁰₂ 14.6/14.9e-4, a²₂ 3.41/2.53e-4) | table tolerance; a and b from the representation at threshold (analytic derivative), threshold fit as control | V1 |
| C9 | Table 2 radii S 0.64/0.55, V 0.394/0.391, T 0.153/0.140 fm² | table tolerance; F′(0) from the discrete dispersion relation | V1 |
| C10 | Table 3 l̄₁ 0.96/−0.10, l̄₂ 4.5/4.8, l̄₄ 4.8/4.2, l̄₆ 14.6/14.5 | table tolerance (algebra on C8, C9 with (4.5)–(4.8)) | V1 |
| C11 | Table 9 GTBII | table tolerance (algebra (A.1)–(A.8), l̄₃ = 0) | V1 |
| C12 | \|F⁰₀\| broad peak ≈ 0.48–0.50 GeV, heights 2.0 / 1.68 | peak node within ±40 MeV; height within 10 % | V1 |
| C13 | ρ pole Re√s 829/844, Im√s 71/82 MeV (Table 4) | Re within ±20 MeV, Im within ±15 MeV; continuation S(s)S*(s̄)=1 with complex-s kernels | V1 |
| C14 | m_ρ 838/859, Γ_ρ 114/128 MeV (Table 5) | Gounaris–Sakurai fit of \|F¹₁\| on nodes 0.5–1.1 GeV (BW as control); mass ±20 MeV, width ±20 % | V1 |
| C15 | g_ρππ 4.9/5.2 | ±5 % (from C14 via (5.5)) | V1 |
| C16 | 4π/g_ρ² 0.61/0.57 | ±10 %; integration over nodes 0.6–1.0 GeV | V1 |
| C17 | Two-BW fit of \|𝓕¹₁\|² (Table 6): m₁ 824/841, Γ₁ 137/164, m₂ 1635/1645, Γ₂ 248/272, A₁ 1.29/1.39, A₂ 0.27/0.28 | fit on nodes 0.5–1.83 GeV; m₁ ±20, m₂ ±40 MeV; Γ ±25 %; A ±15 % | V1 |
| C18 | f₂ BW of \|F⁰₂\| (Table 7): m 1366/1388, Γ 164/188 MeV | fit on nodes 1.0–1.83 GeV; m ±25 MeV; Γ ±25 % | V1 |
| C19 | g_f₂ππ 8.5/8.9 | ±5 % (from C18 via (5.13)) | V1 |
| C20 | Forward amplitudes (Fig. 5, red) | Im T₀₊, T₀₀, T^{I_t=1} peak positions within ±30 MeV, heights within 10 % of the red paper curves | V1 |
| C21 | Watson step fixes D0 (Fig. 20) | before: min η⁰₂ < 0.7 for at least two colours; after: η⁰₂ ≥ 0.98 and δ⁰₂ 90° crossing within ±30 MeV of the paper's | V1 |
| C22 | Stability at green (Fig. 17): M=60, L=23, N_B=46, M_reg=5×10² | per-wave RMS difference from the base run ≤ 10° for S0, P1, S2 below 1.4 GeV | V1 |
| C23 | Minimal feasible caps r⁰₀ 2e-2, r¹₁ 1.5, r⁰₂ 4; r¹₁ = r⁰₂ = 1 infeasible | feasible values found feasible (V1); infeasibility certified by a dual ray (V2) | V1 + V2 |
| C24 | ρ mass varies ≈ 80 MeV for r¹₁ = 1.5–2.5 | spread of the \|F¹₁\| peak node in 80 ± 30 MeV | V1 |
| C25 | n_SR = 3, 4, 5 give the same \|F\| | peak positions within ±10 MeV, heights within 5 % | V1 |
| C26 | Below 1.2 GeV the phases agree with the s₀ = 1.2 GeV results of 2309 | S0, P1, S2 at green vs the 2309 paper's Fig. 9/10 red curves: RMS ≤ 10° below 1.2 GeV | V1 |
| C27 | Removing P1 FESR removes the ρ; removing D0 FESR does not completely remove the f₂ | no δ¹₁ 90° crossing below 1.2 GeV without P1 FESR; max δ⁰₂ ≥ 45° below 1.83 GeV without D0 FESR | V1 |
| C28 | Uncertainty model = blue–red spread | definitional; our P8 extension reports V2-backed [min, max] ranges alongside | — |
| C29 | Runtime | informational | — |

## 2. Gates (plan Phase A)

- **G1**: our kernels = authors' reference export to ≤ 1e-30 relative (each exception explained in the derivation), and
  = Arb rebuild within Arb radii.
- **G2**: a production encoding with a measured cost; the choice is recorded before G3 runs.
- **G3** = C1. On FAIL: stop production; diagnose (regulariser scale, full space vs reduced, ‖ρ‖₄ at the optimum,
  omitted-wave unitarity) without changing any input; report.

## Amendment 1 (2026-10-01, before any 2403 analysis quantity is evaluated): analysis prescriptions


- Off-node and complex-s evaluation of the amplitude uses the disk interpolant w_j(z) = (i/2M)(1 − (z/z_j)^M)(z_j + z)/(z_j − z) − (same with z_j → 1/z_j), z_j = e^{iφ_j}, for every s-channel discontinuity (σ and the s-index of ρ); it satisfies Im w_j(z_k) = δ_jk and Re w_j(z_k) = K_kj, i.e. it is the unique band-limited analytic function that the authors' K represents at the nodes. Crossed channels use Q_ℓ at the same (real or complex) s. Control: the authors' midpoint sum (Δ/π) x'_m [1/(x_m − s) − 1/(x_m − ν₀)] for the s-channel; the difference is reported for every quantity that uses off-node evaluation.
- Threshold parameters (C8): a and b from a least-squares fit of Re f^I_ℓ(s)/k^{2ℓ} = a + b k² + c k⁴ on 12 evenly spaced k² ∈ (0, 0.25] using the interpolant evaluation; the same fit on the first 4 nodes is the control.
- Radii (C9): F′(0) = Σ_m (Δ/π) x'_m Im F_m / x_m² (the authors' midpoint construction of K₀, differentiated); ⟨r²⟩ = 6 F′(0)/m_π² in fm² with ħc = 197.327 MeV fm; control: derivative of the interpolant at s = 0.
- ρ pole (C13): zero of S¹₁(s) = 1 + iπ√((s−4)/s) f¹₁(s) (S = 1 + ih, h = π√((s−4)/s) f) continued to complex s in the lower half plane of the first sheet (equivalently the pole of S on the second sheet, since S_II = 1/S_I), found by Newton from √s = 0.83 − 0.04i GeV; reported as √s_p in MeV.
- Fits (C14, C17, C18): windows as in PREREGISTRATION.md; Gounaris–Sakurai and Breit–Wigner forms of (5.2), (5.3), (5.9) of v4; unweighted least squares on node values.
- Forward amplitudes (C20): T^I(s, 0) = 2 Σ_{ℓ≤19} (2ℓ+1) f^I_ℓ(s) at the nodes; T₀₊ = (T¹ + T²)/2, T₀₀ = (T⁰ + 2T²)/3, T^{I_t=1} = T⁰/3 + T¹/2 − 5T²/6.


## Amendment 2 (2026-10-01, from the blind figure extraction only; no 2403 phase result of ours exists yet)
- C4 and C21: the η threshold 0.98 rejects the paper's own post-Watson curves (Fig. 4 P1: η = 0.978 at blue node 29,
  0.986 at red node 35). Changed to min η ≥ 0.97 so that the rule accepts the published figure itself.
- C4–C7 comparisons are node by node against `/playpen1/shiqiu/gtb2403_runs/paper_data/` (extraction script copied to
  `scripts/gtb/extract_2403_figures.py`; dots on the nodes to ≤ 3e-5 GeV). Missing paper dots (Fig. 4 D0 δ at nodes
  2–12 even, red S0 node 39) are skipped, not interpolated.
- C23–C25: the Fig. 18–19 runs match none of the blue/green/red points (differences 0.17–1.3 in |F|); those claims are
  compared in shape (peak positions/heights relative to the same run's baseline), never point to point.

## Amendment 3 (2026-10-01, before any 2403 phase result of ours exists)
- Phase comparisons (C4, C7, C21) use d = ((δ_ours − δ_paper + 90°) mod 180°) − 90°: S = η e^{2iδ} fixes δ only modulo
  180°, and the authors' plotting notebook applies hand-made branch shifts (P5 #8). Crossings (C5, C6) use our δ
  unwrapped from threshold, as frozen in §0.
- Every claim rule is first run with the paper's own extracted curves in place of ours; a rule that rejects the paper
  itself is a defect of the rule and is amended (with this record) before our leaves are evaluated.

## Amendment 4 (2026-10-01, before any C2 solve exists): compute budget for the Fig. 3 scan
- C2 is evaluated on k ∈ {2, 6, 10, …, 50} (13 values of 𝔉₀ = k/50 × 0.163108) × {max, min} 𝔉₁ = 26 solves, compared
  with the paper's dot at the same k (the paper's Fig. 3 dots exist for every k). Rule unchanged otherwise: RMS of the
  𝔉₁ differences over these 26 dots ≤ 3e-4; tip 𝔉₀ within 1 % (k = 50 feasible and k = 51 infeasible is not required);
  tree-level dot inside. Reason: ~2 h per reduced-basis solve on the shared machine (load ~100 on 80 cores).
- These scan solves use the reduced basis (an inner approximation: max values are lower bounds, min values upper bounds);
  the encoding is the one selected at G2. They are reported as inner bounds, never as V2-certified boundaries.

## Amendment 5 (2026-10-01 22:55Z, before any M = 50 certificate, omitted-wave η or phase result exists)
State when written: the unconverged SDPB iterates of the three ℓ∞-inner leaves (ni 44/46/48, reduced basis) were
visible; their feasible-side objective was ≈ 1e-3 above the paper's C1 values. Nothing below changes a tolerance.
- C1 by one-sided certificates: a V1-certified, ℓ4-feasible leaf with 𝔉₁ > paper + 1e-4 is FAIL at that ni whatever
  the upper bound (the bracket [lower, ∞) is disjoint from the band); a V2 bound below paper − 1e-4 is FAIL likewise.
  This only makes explicit that no full-space V2 solve is needed for such a verdict.
- On a C1 FAIL, Phase B is computed for two regulariser variants and both are reported side by side, never chosen by
  agreement with the paper:
  (A) the authors' scale M_reg = 1e2: representative = the ℓ∞-inner leaf (|ρ_i| ≤ R/3775^{1/4}, inside the ℓ4 ball),
      control = the ℓ∞-relax leaf (|ρ_i| ≤ R); their spread is reported as the encoding uncertainty of (A).
  (B) the physics-selected scale, the 2309 rule transferred unchanged: ladder M_reg ∈ {100, 30, 10, 3, 1}
      (R = M_reg · Mrho), ℓ∞-relax box |ρ_i| ≤ R, ni = 46 pre-Watson 𝔉₁ max, reduced basis, 128 bit. Omitted-wave
      η = max over all M nodes of |1 + i h| for ℓ = 20…23 (from the L = 24 export of the same kernels).
      M_reg* = the largest rung with η ≤ 1.02 whose L = 24 re-solve gives 𝔉₁ within 2 % of the L = 20 value.
  For each variant: initial + Watson leaves at ni 44/46/48 and claims C4–C21. The C2 scan (amendment 4) uses (A).

## Amendment 6 (2026-10-01 23:25Z, after C1 = FAIL and the rung-100 omitted-wave η, before any other rung is solved)
- Seen when written: C1 FAIL (V1-certified inner leaves 1.5–1.8e-3 above the paper at ni 44/46/48); omitted-wave η at
  M_reg = 100: 2.14 (inner, ni 46) and 2.83 (relax, ni 46), maximum at node 47.
- Variant (B) ladder: if no rung of {100, 30, 10, 3, 1} has η ≤ 1.02, it continues with the same factor
  (0.3, 0.1, 0.03, …) until a rung passes or the problem becomes infeasible; infeasible before passing = variant (B)
  does not exist and is reported so. η is evaluated at every node (no energy cut), exactly as in 2309.

## Amendment 7 (2026-10-02 04:00Z, before any text-reading solve): the text's regulariser normalisation as a control
- Seen when written: C1 FAIL; ladder F1(ni 46) for M_reg 100/30/10/3/1 = -0.02546/-0.02651/-0.02732/-0.02793/-0.02823;
  Fig. 3 dots lie on F0 = k/50 x 0.1632 (the shipped script uses 0.163108; calibration by the tree dot to 6e-6).
- Control (not a mainline, not a tuning): the text's (8.56) reading ||rho||_4 <= 100 (no division by Mrho), bracketed
  by l-inf inner (|rho_i| <= 100/3775^(1/4)) and relax (|rho_i| <= 100) at ni 46, reduced basis, 128 bit, V1.
  Verdict wording fixed now: if the paper's -0.028770 lies inside the certified bracket, "the published point is
  consistent with the text's normalisation (the figure run differs from the shipped code)"; if outside, "neither the
  code's nor the text's normalisation reproduces the published point".

## Amendment 8 (2026-10-02, before any P_W kernel, solve or analysis exists): the well-posed version P_W
- Definition, principles and gates: `docs/gtb2403/WELLPOSED_ZH.md` (sha256 4c12ec437ae94228), frozen with this amendment.
- P_W = the authors' code problem with exactly four changes: W1 FESR moments integrated and normalised to the cell
  boundary x_b = x(n0 pi/M) = 189.16 (both sides); W2 the same unitarity cone for all l < L_u = 40 at all nodes
  (control L_u = 60 at ni 46; if |dF1| > 1e-4 P_W becomes L_u = 60); W3 the same cone at all M-1 phi-midpoints for
  l < 20 (interpolant), by constraint generation until the certified max |S| - 1 <= 1e-6; W4 the code's l4 regulariser
  kept as a technical bound (l-inf relax), required inactive (||rho||_4 <= R at the optimum and |dF1| < 1e-6 for +10 % R),
  otherwise reported as active. Full space, 128 bit, V1 + V2. No physical input of the paper changes.
- New prescription (P_W and as a control for variant A): threshold parameters from the interpolant below threshold,
  12 points s in [2, 3.8], fit of f_l/k^(2l) = a + b k^2 + c k^4 with k^2 = (s-4)/4, extrapolated to k^2 = 0.
- C1-C27 rules are unchanged; P_W verdicts are reported next to variant A and variant B, never used to pick a setting.
