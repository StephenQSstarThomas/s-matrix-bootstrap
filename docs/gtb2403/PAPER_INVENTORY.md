<!-- Derived inventory (our own text) of arXiv:2403.10772 v1-v4, prepared 2026-10-01 from the arXiv sources; author code is referenced by file/line, not reproduced. Source: /playpen1/shiqiu/gtb2403_prep/paper/PAPER_2403_INVENTORY.md -->
# PAPER_2403_INVENTORY — He & Kruczenski, arXiv:2403.10772 ("GTB: Pion amplitudes and low energy parameters")

Prepared 2026-10-01 (read-only scoping). All paths below are under `/playpen1/shiqiu/gtb2403_prep/paper/` unless absolute.
Equation/figure/table numbers refer to **v4** (numbers verified against the arXiv v4 PDF and a local tectonic build, `labels_v4.txt`).
Code references: `CODE` = released ancillary `GTB_numerics.m` of v3/v4 (= GitHub snapshot
`/home/shiqiu/s-matrix-bootstrap/references/upstream-gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.m`,
which differs from `src_2403v4/anc/GTB_numerics.m` only by a trailing newline). `NB` = `src_2403v4/anc/GTB_numerics.nb` (matrix generator).
`CODE-v1` = v1/v2 ancillary (identical md5 in v1 and v2).

Saved materials: `src_2403v{1..4}.tar.gz` + extracted `src_2403v{1..4}/`, `2403.10772v{1..4}.pdf`, `v{1..4}.txt` (pdftotext),
`v2_v3_worddiff.txt`, `abs_2403v*.html`, `inspire_2403.json`, `PoS_CD2024_040.pdf/.txt` (proceedings),
`src_2505v3.tar.gz` + `src_2505v3/`, `2505.19332v3.pdf`, `fig_png/` (renders), `fig3_shape_dots_extracted.txt` (vector read-off of Fig. 3).

---------------------------------------------------------------------------------------------------

## 1. Version history

| ver | date (UTC) | arXiv comment | pages/figs | TeX file | ancillary code |
|---|---|---|---|---|---|
| v1 | 2024-03-16 | "74 pages, 16 figures, ancillary files contain program..." | 74 / 16 | `QCDsumrule2_30.tex` | CODE-v1 (md5 b0825c0f) |
| v2 | 2024-09-15 | "section 5.2.3 updated with analysis of two Breit-Wigner resonances, appendix B added with test on numerical parameters" | 76 / 17 | `arxiv_submission_3.tex` | identical to v1 |
| v3 | 2025-02-27 | "Numerical method improved with low energy kinematic constraints on spectral densities (sec 2.2) and Watsonian unitarization (sec 2.5); updated figures and calculation results; included pion spin-2 charge radius (sec 4.2); ancillary files updated. Overall idea and results remain essentially the same" | 84 / 20 | `arxiv_submission_3_v5.tex` | CODE (md5 c966f9bd) |
| v4 | 2025-08-08 | "(labelled v3:) 85 pages, 20 figures. Definitions in section 6.1.1 updated in agreement with the rest of the paper" | 85 / 20 | `arxiv_submission_3_v10.tex` | identical to v3 |

**Journal status:** no journal version. INSPIRE record 2769436 has no `publication_info` (record updated 2025-08-12; still none as of 2026-10-01).
The only refereed-venue item is the proceedings *PoS(CD2024)040* (DOI 10.22323/1.479.0040, 2026; saved as `PoS_CD2024_040.pdf`). Its tables are
identical to v3/v4 (a, b, r², l̄, ρ/ρ′/f2 masses). By comparison, the companion 2309.12402 appeared in PRL 133 191601 and PRD 110 096001, and 2505.19332 in PTEP 2026 063B05.

### What changed, by step

**v1 → v2 (results unchanged; numerics unchanged)**
- §5.2.3: new two-Breit–Wigner fit to |𝓕¹₁|² (new eq. for F¹₁, new Table "ρ(1450)", new fig `rho1450twopole.pdf` replacing `rho1450.pdf`; v2 values m1 825/828, Γ1 128/126, m2 1632/1636, Γ2 180/210, A1 1.23/1.21, A2 0.19/0.20, A1−A2 1.04/1.01).
- New App. B "Test on numerical parameters M, L" (fig `psandetaconverge.png`: M=50,L=19 red; L=23 green; M=60 purple). Text notes that the f2 rise shifts for M=60 "because changing M shifts n0, so the FF asymptotic bounds differ".
- Cosmetic (author order in citation, "elasticity" → "inelasticity").

**v2 → v3 (major numerical revision; all results changed)**
- Spectral densities: point values → cubic B-spline basis with N_B = 40, ρ = κ_ℓ ρ̂, ρ̂(0) = 1 (§2.2, §8.3.2, eqs. (8.42)–(8.49)).
- Positivity matrix: v1/v2 [[1,S,𝓕],[S*,1,𝓕*],[𝓕*,𝓕,ρ]] at points → rescaled matrix (8.61) with entries h̃, 2Im ĥ, F, 2Im F̂, ρ̂.
- FF asymptotics: v1/v2 caps on |𝓕|² (rescaled FF): |𝓕⁰₀|² ≤ 3e-8, |𝓕¹₁|² ≤ 2e-6, |𝓕⁰₂|² ≤ 4e-2 (v1 eq. (8.65)) → v3 caps on |F⁰₀| ≤ r⁰₀ and |s F| ≤ r·(Brodsky–Lepage), r = (0.05, 2, 6) (eqs. (8.74)–(8.75)).
- Unitarisation: v1/v2 plain linearised unitarity functional on 6 waves → "Watsonian" unitarisation (§2.5, eqs. (2.26)–(2.29), (8.76)–(8.78)).
- Fig. 3 now has three points (blue, green, red) instead of two (blue, red). The v1/v2 remark that lower-edge points "display numerical instabilities where the sum rules are satisfied by discontinuous jumps near threshold" was removed.
- ε^χ in the text corrected from 3×10⁻² (v1/v2 text) to 2×10⁻³ (the value the v1 code already used).
- New spin-2 (mechanical) radius ⟨r²⟩_T; new Table 8 (summary of numerical parameters); App. B expanded to B.1–B.4 (Figs. 17–20); the discussion of particle production at f2 and ρ(1450) was removed; SR-moment display (8.71) corrected (P1, D0 moment weights written with the right powers 1/x, 1/x²).
- Λ_ℓ typo fixed: v1/v2 text ((√s−2)/(√s+1))^{ℓ/2} → ((√s−2)/(√s+2))^{ℓ/2} (eq. (8.62); NB uses +2).
- All tables updated (see §3). **Fig. 7 (`rhopole.pdf`) was NOT updated: same md5 (b9f620) in v1–v4, while Table 4 changed (v1: 829/832, 63/64 → v3: 829/844, 71/82).**

**v3 → v4 (no numerical change)**
- §6.1.1 general-ℓ form-factor definitions rewritten (isospin-index notation; normalisation of 𝓕^{I}_ℓ for general ℓ). √(3π/15) → √(π/5) in (2.23e) (identical value).
- (2.20a) gains the m_π⁴ factor; T_{μν} split into gluon and fermion parts; t = −(s−4)(1−cosθ)/2 sign fixed.
- Footnotes citing 2505.19332 (B-spline replaced by better parametrisation; iterated Watson converges, "largely independent of the functional chosen"); GitHub link added.
- `psandeta.pdf` re-rendered: only the δ⁰₀ value of the last point (1.83 GeV) is drawn on another branch mod π (blue/green no longer wrap to ~10–50°). Same data.
- Tables, ancillary code and all other figures are identical to v3.

### Known text-vs-code mismatches (status in each version)

| item | v1/v2 text | CODE-v1 | v3/v4 text | CODE (v3/v4, GitHub) | fixed? |
|---|---|---|---|---|---|
| FESR tolerance P1 | 5×10⁻⁶ (v1/v2 eq. (8.62)) | esrP1 = 6e-6 | ε¹₁ = 5×10⁻⁶ (8.73) | eP1 = 6e-6 | **No, present in all 4 versions** |
| FESR tolerance D0 | 6×10⁻⁶ | esrD0 = 5e-6 | ε⁰₂ = 6×10⁻⁶ (8.73) | eD0 = 5e-6 | **No** |
| FESR tolerance S0 | 1×10⁻⁷ | 1e-7 | 1×10⁻⁷ | 1e-7 | consistent |
| ε^χ | 3×10⁻² | 2e-3 | 2×10⁻³ | 2e-3 | fixed in v3 |

Evidence that the **code values were used to make the figures**: the Fig. 19 legends (`SRerr11.pdf`, `SRerr02.pdf`) show the central (orange) curves
at ε¹₁ = **6×10⁻⁶** and ε⁰₂ = **5×10⁻⁶**, with variations ε¹₁ ∈ {1e-6, 6e-6, 1e-5} and ε⁰₂ ∈ {3e-6, 5e-6, 7e-6}. The text's own prescription
("ε estimated from the O(α_s) term", below (8.71)) gives: ‖O(α_s) part of the 3 moments‖₂ = 1.23e-7 (S0), **5.97e-6 (P1)**, 5.63e-6 (D0),
computed with NB inputs. So P1 = 6e-6 also follows from the stated rule. Conclusion: eq. (8.73) has P1 and D0 swapped. Use the code values.

Other text/code discrepancies found in this pass (v4 text vs CODE):
1. **Regulator normalisation.** (8.56) states ‖[ρ₁,ρ₂]‖₄ ≤ M_reg = 10². CODE has `norm(rho/Mrho,4) <= 10^2` with Mrho = M² + M(M+1)/2 = 3775, i.e. ‖ρ‖₄ ≤ 3.775×10⁵. The 2505 main code is the same.
2. **Second chiral norm.** (8.69) lists ‖[f⁰₂, f²₂, **f¹₁**]‖ ≤ ε^χ. CODE constrains `[fchiD0; fchiD2; fchiF1]`, i.e. f¹₃ (F1), at s* = ½, 1, 3/2, 2. The text is a typo: f¹₁ cannot be made small, because it is the tree-level wave.
3. **Watson target normalisation.** (8.76)/(8.78) give h̃ = Im F̂ / F* (Re h̃_fnl = Re F·Im F̂/|F|²). CODE uses `2*real(FF3).*ImFFh3./|FF3|²` (factor 2). The factor 2 is the one implied by the rank-1 saturation of matrix (8.61), whose (2,3) entry is 2 Im F̂. So the code is self-consistent and the text drops the 2.
4. **Unitarity-cone count.** (8.55) text says 3ML/2 cones (= 1425 for L = 19). CODE has 3·l·M = 1500: l = 10 waves per isospin, ℓ = 0..19. Also, "L" means maximal spin 19 in the text but the Fig. 17 legend means waves per isospin (L = 10, 12 ↔ ℓ_max = 19, 23).
5. Variable count "**M** ∼ 5,000" (§2.1, §2.4) is really **M** = 1 + 2M + M² + M(M+1)/2 + 3M + 3N_B = 4146.
6. P1 FESR n = 1 value printed 1.8×10⁻⁵ in (8.71); NB inputs give 1.872×10⁻⁵, which should round to 1.9. The other eight values round correctly.
7. "g_ρ = 4.8 ≃ g_ρππ" (after (5.8)) is stale from v1 (4π/g_ρ² = 0.55 → 4.78). The v3/v4 values 0.61/0.57 give g_ρ = 4.54/4.70.
8. CODE plotting bug: `deltaP1=angle(...)` lacks `*deg`, so P1 is plotted in radians. This affects plotting only.
9. α_s, m_u, m_d are not given in the text; they appear only in NB (see §2). The S0 FESR uses m_u² + m_d² where (2.21h) writes N_f m_q².

---------------------------------------------------------------------------------------------------

## 2. Numerical setup of the latest version (v4), with code agreement

Units m_π = 1 (m_π = 139.57 MeV in NB/CODE). Two flavours, N_c = 3.

| # | Item | Paper (v4) value / definition | Eq./loc. | CODE / NB | Agree? |
|---|---|---|---|---|---|
| 1 | Low-energy data | m_π = 139.57 MeV, f_π = 92 MeV | §2.2, (2.11) | NB `Mπ=139.57`; CODE `mpi=139.57; fpi=92` | yes |
| 2 | α_s | not stated numerically | (2.21)–(2.23) | NB `αs2GeV = 0.3145887440`, equal to one-loop α_s(s0 = (2 GeV)²) of the 2505 NB formula 1.3963/ln(0.41218 s); CODE `FP1asym = 6.87076 = 16π α_s f_π²` gives α_s = 0.31459 | implied |
| 3 | Quark masses | "m_q" only; N_f m_q² in (2.21h) | (2.21h) | NB `mu2GeV = 3.6/Mπ`, `md2GeV = 6.5/Mπ`; S0 FESR ∝ (m_u² + m_d²), equivalent to m_q = 5.25 MeV (rms) | implied |
| 4 | Condensates | omitted (suppressed by s0); no condensate inputs | §2.6 | none | yes |
| 5 | Currents | S0: j⁰₀ = m_q(ūu + d̄d); P1: ½(ūγ₊u − d̄γ₊d); D0: T^{++} | (2.21a), (2.22a), (2.23a) | 3 FFs + 3 spectral densities (S0, P1, D0) | yes |
| 6 | FF ↔ rescaled FF | 𝓕⁰₀ = √(6π)/(16π³) s^{-1/4}((s−4)/4)^{1/4} F; 𝓕¹₁ = √(4π/3)/(8π³) s^{-1/4}((s−4)/4)^{3/4} F; 𝓕⁰₂ = (1/4π³)√(π/5) s^{-1/4}((s−4)/4)^{5/4} F | (2.21e), (2.22e), (2.23e) | in NB (κ_ℓ, Λ3) | yes |
| 7 | Low-energy spectral densities | κ₀ = (2π)⁻⁴ (3/16π)(1−4/s)^{1/2}; κ₁ = (2π)⁻⁴ (s/24π)(1−4/s)^{3/2}; κ₂ = (2π)⁻⁴ (s²/160π)(1−4/s)^{5/2} | (2.20), (8.43) | ρ = κ_ℓ ρ̂ | yes |
| 8 | Amplitude ansatz | Mandelstam double-dispersive with T₀, σ₁,₂, ρ₁,₂; kernels subtracted at ν₀ | (2.1), (8.2)–(8.14) | NB → B1, B2, A1, A2, Bhat1, Bhat2 | yes |
| 9 | Map / grid | z(ν) (8.15), x(φ) = ν₀ + (8 − 2ν₀)/(1 + cos φ) (8.16); φ_j = (j − ½)π/M (8.28) | (8.15)–(8.29) | `vnu = nu0+(8-2*nu0)./(1+cos(((1:M)-1/2)*pi/M))` | yes |
| 10 | ν₀ | −20 ("tested a few values") | footnote to (8.28) | `nu0=-20` | yes |
| 11 | M | 50 | Table 8, after (8.55) | `M=50` | yes |
| 12 | L (max spin) | 19 | after (8.50), Table 8 | NB `L=20` (ℓ < 20); CODE `l=10` per isospin | yes |
| 13 | s0 | √s0 = 2 GeV, s0 = (2000/139.57)² = 205.34; PSD imposed for s_j < s0 | (8.54), Table 8 | `s0=(2000/mpi)^2` | yes |
| 14 | Grid facts (derived) | footnote: last point s₃₉ ≈ 1.83 GeV | footnote after (8.71) | n0 = 39 (√s₁ = 0.2793 GeV … √s₃₉ = 1.8309 GeV); 11 points above s0: √s = 2.016 … 43.53 GeV | yes |
| 15 | Discrete principal value | K̂_{j1j2} = −(1/2M)(1 − (−)^{j1−j2}) cot(π(j1 − j2)/2M) | (8.33)–(8.35) | NB | yes |
| 16 | FF dispersion | F(0) = 1, one subtraction; Re F_j = 1 + (K̂ + K₀) Im F | (2.7), (8.24), (8.40) | `ReFF=KFF*ImFF+1` | yes |
| 17 | Unitarity | Re h̃² + Im h̃² ≤ 2 Im ĥ, all waves ℓ ≤ 19, all M points | (8.52), (8.64) | `norms([Reht,Imht],2,2)<=sqrt(2*Imhh)` | yes (count: see mismatch 4) |
| 18 | Rescaling | Λ_ℓ(s) = ((√s − 2)/(√s + 2))^{ℓ/2}; h̃ = Λ⁻¹h, Im ĥ = Λ⁻² Im h, Im F̂ = Λ⁻¹ Im F; kernels built in high precision, then multiplied by machine-precision variables | (8.60)–(8.64) | NB `Λ[s,l]`; CODE `ImFFh3=.../Lambda3` | yes |
| 19 | PSD matrix | [[1, h̃, F],[h̃*, 2Im ĥ, 2Im F̂],[F*, 2Im F̂, ρ̂]] ⪰ 0, for S0, P1, D0 at s_j < s0 (3×39 = 117 matrices) | (8.61), (8.54) | `Bmat(:,:,i) == [...]; Bmat>=0`, i = 1..3n0 | yes |
| 20 | B-spline basis | ρ̂(s) = Σ_{i=1}^{N_B} ρ_{ℓ,i} B_i(s), cubic (p = 3), knots {0,0,0,0, 36 equally spaced interior, s_{n0}×4} on [0, s_{n0}]; ρ_{ℓ,i} ≥ 0; ρ̂ extended to 0 < s < 4 as F²(s) | (8.45)–(8.48) | NB `Nb=40; pd=3; innerKnots=i*vx[[n0]]/(Nb-pd), i=1..Nb-pd-1`; CODE `specC>=0` | yes |
| 21 | N_B | 40 | Table 8 | `nB=40` | yes |
| 22 | Low-energy ρ behaviour | ρ_{ℓ,1} = 1 (ρ̂(0) = 1) | (8.47) | `specC(1)==1; specC(nB+1)==1; specC(2*nB+1)==1` | yes |
| 23 | Chiral matching points | s* = {½, 1, 3/2, 2} | (8.66) | `vsunphys=1/2+((1:4)-1)/2` | yes |
| 24 | Chiral ratios | R₀₁ = 3(2s − 1)/(s − 4), R₂₁ = 3(2 − s)/(s − 4) | (8.67) | same | yes |
| 25 | Chiral constraints | ‖[f⁰₀ − R₀₁ f¹₁, f²₀ − R₂₁ f¹₁]‖₂ ≤ ε^χ; ‖[f⁰₂, f²₂, f¹₁(sic)]‖₂ ≤ ε^χ | (8.69) | second norm is over [f⁰₂, f²₂, **f¹₃**] | text typo (mismatch 2) |
| 26 | ε^χ | 2×10⁻³ ("minimal tolerance before the tree-level point is excluded"; from 2309 §4.2) | after (8.69), Table 8 | `echi=2*10^(-3)` | yes |
| 27 | FESR moments | n_SR = 3: S0 n = 0, 1, 2; P1 n = −1, 0, 1; D0 n = −2, −1, 0; normalisation 1/s0^{n+2} (S0, P1), 1/s0^{n+3} (D0) | (2.21h), (2.22h), (2.23h), (8.70) | NB `SVZSR` with SVZS0[0..2], SVZP1[−1..1], SVZD0[−2..0] | yes |
| 28 | FESR RHS | S0: 3.1e-7, 2.1e-7, 1.6e-7; P1: 5.6e-5, 2.8e-5, 1.8e-5; D0: 5.1e-5, 2.6e-5, 1.7e-5 (O(α_s), identity operator only) | (8.71) | from NB inputs: S0 3.113e-7, 2.075e-7, 1.556e-7; P1 5.617e-5, 2.809e-5, **1.872e-5**; D0 5.134e-5, 2.567e-5, 1.711e-5 | yes (rounding, mismatch 6) |
| 29 | FESR discretisation | Σ_{j ≤ n0} Δφ x′(φ_j) κ_ℓ(s_j) ρ̂_j s_j^n / s_{n0}^{n+n_ℓ}; integration ends at s_{n0} (1.83 GeV) | (8.49), (8.72) | NB `intS0/intP1/intD0` | yes |
| 30 | FESR tolerances (2-norm of 3 residuals) | ε⁰₀ = 1e-7, ε¹₁ = 5e-6, ε⁰₂ = 6e-6 ("≈ 10–30 % of moments") | (8.73) | eS0 = 1e-7, **eP1 = 6e-6, eD0 = 5e-6** | **No** (swap; Fig. 19 agrees with code) |
| 31 | FF asymptotic form | |F⁰₀| ∼ 4π α_s f_π² m_q² ln(s/m_π²)/s ("crude estimate"); |F¹₁| ≃ 16π α_s f_π²/s; |F⁰₂| ≃ 48π α_s f_π²/s | (2.21i), (2.22i), (2.23i) | FP1asym = 6.87076, FD0asym = 20.6123 (= 3×FP1asym) | yes |
| 32 | FF caps above s0 (j = n0+1..M, 11 points) | |F⁰₀(s_j)| ≤ r⁰₀ (a constant: "surprisingly ... constant instead of 1/s"); |s_j F¹₁| ≤ r¹₁·16πα_s f_π²; |s_j F⁰₂| ≤ r⁰₂·48πα_s f_π² | (8.74) | `abs(FF(n0+1:M))<=rS0; abs(vnuhigh'.*FF(P1))<=rP1*FP1asym; ...` | yes |
| 33 | Relaxation factors | r⁰₀ = 5×10⁻², r¹₁ = 2, r⁰₂ = 6 ("r = 1 infeasible"; App. B.2: minimal feasible 2e-2, 1.5, 4) | (8.75), Fig. 18 | rS0 = 5e-2, rP1 = 2, rD0 = 6 | yes |
| 34 | Regulariser | ‖[ρ₁, ρ₂]‖₄ ≤ M_reg = 10² ("avoids highly oscillating double spectral densities, renders dual feasible"; from He:2021eqn) | (8.56), Table 8 | `norm(rho/Mrho,4) <= 10^(2)` | **normalisation differs** (mismatch 1) |
| 35 | Functionals | 𝔉₀ = 2(f⁰₀(3) + 5f⁰₂(3)), 𝔉₁ = 2(3f¹₁(3) + 7f¹₃(3)) | (2.24) | `F0=2*(f3S0+5*f3D0); F1=2*(3*f3P1+7*f3F1)` | yes |
| 36 | Point selection | Fig. 3 made by scanning 𝔉₀ and max/min 𝔉₁; blue/green/red = 3 upper-edge points nearest the tree-level black dot; spread = error estimate | §2.4, §8.4 "Maximization..." | `ni=46; F0fix=ni/50*0.163108` (ni = 44, 46, 48 = blue, green, red → 𝔉₀ = 0.14354, 0.15006, 0.15658), `F0==F0fix; maximize(F1)` | **numbers only in code**; Fig. 3 vector read-off confirms the grid k/50·0.1632 and the 3 points |
| 37 | Tree-level black dot (derived) | 𝔉₀ = 0.145746, 𝔉₁ = −0.029152 (f_π = 92 MeV) | (2.15), (2.24) | read-off (0.14574, −0.02915) | yes |
| 38 | Watson step (one iteration) | maximise Σ[Re h̃_fnl Re h̃ + Im h̃_fnl Im h̃] − Σ Im ĥ over S0, P1, D0 and S2, D2, F1 at s_j < s0; for S0/P1/D0 the target is h̃_fnl from the old FF (Watson, (8.76)); for S2/D2/F1 the target is the old partial waves; "can be repeated until convergence" (not done in the paper) | (2.29), (8.76)–(8.78) | second CVX solve, `maximize(v1+v2)`, with all constraints repeated (F0 not fixed in step 2) | yes (except factor 2, mismatch 3) |
| 39 | Solver | CVX + MOSEK, double precision; no precision settings given | §8 footnote, §8.4 | `cvx_begin sdp; cvx_solver mosek` (default cvx_precision) | yes |
| 40 | Runtime | "a few minutes on an average PC"; Matlab "∼20 min on an average laptop" | §1, footnote §1 | – | – |
| 41 | Kernel files | built by NB | §8.3 | 17 `.mat` files (A1, A2, B1, B2, Bhat1, Bhat2, fsigmaF, frhoF, fsigmachi, frhochi, Lambda3, KFF3, MBSpline3, intS0, intP1, intD0, SVZSR) are **not shipped**; Mathematica is required | gap |

Not in the paper or the code: no η/inelasticity constraint beyond unitarity; ρ̂ = |F|² is **not** imposed for 4 ≤ s ≤ 16 (it is stated physically but enforced only through the PSD matrix); no condensates; no D0 3-loop.

---------------------------------------------------------------------------------------------------

## 3. Figures and tables (v4 numbering)

"Blue/green/red" = ni = 44/46/48 (𝔉₀ = 0.14354/0.15006/0.15658), each after one Watson step unless stated. Figures marked **vector** have dot positions recoverable from the PDF (colour-coded `<path>` fills; counts given). Figure numbers are v4.

| Fig | file(s) in tarball (`figures/`) | content / axes | runs | exp. overlay | digitise | key read-offs (eye ±1 px unless "vector") |
|---|---|---|---|---|---|---|
| 1 | UVIR.pdf | schematic UV–IR bridge | – | – | n/a | – |
| 2 | GTB_schematic.pdf | schematic of GTB | – | – | n/a | – |
| 3 | shape.pdf (Inkscape, vector) | 𝔉₀ (0–0.165) vs 𝔉₁ (−0.033 to 0.001); yellow upper/lower boundary (50 + 50 dots), black tree-level line and dot; blue/green/red | 50 𝔉₀ values × max & min 𝔉₁, initial maximisation only | no | **vector, done**: `fig3_shape_dots_extracted.txt` | 𝔉₁^max at 𝔉₀ = 0.14362/0.15015/0.15668: **−0.027229 / −0.028770 / −0.030397** (±3e-5); black (0.14574, −0.02915); tip 𝔉₀ ≈ 0.1632, 𝔉₁ ≈ −0.0324 |
| 4 | psandeta.pdf (Inkscape, vector; 462–463 dots per colour) | 6 panels δ & η vs √s (0.28–2.0 GeV): δ⁰₀ (0–220°), δ⁰₂ (0–220°), δ²₀ (−45–0°), δ²₂ (−8 to +1.3°), δ¹₁ (0–350°), δ¹₃ (≈ −8 to +8°); η panels 0–1.1 | blue/green/red, 39 points each | yes: gray dots [Protopopescu 73, Losty 74, Hyams 75]; gray lines PY [Pelaez:2004vs] (S0 with K K̄ removed) | vector | after Watson η = 1 for all GTB points; δ¹₁ = 90° near √s ≈ 0.83–0.85 GeV (data ≈ 0.77); second P1 rise from ≈ 1.5 GeV (δ¹₁ ≈ 310° at 1.83); δ⁰₂ passes 90° near 1.35–1.40 GeV; δ²₂ positive (blue) or negative (red) above 1 GeV, i.e. not fixed; δ²₀ ≈ −30 to −45° at 1.5–1.8 GeV |
| 5 | forwardRe.pdf, forwardIm.pdf (vector) | Re/Im T₀₊, T₀₀, T^{I_t=1}(s, t=0) vs √s 0.28–1.5 GeV (Re −1.9 to 2; Im 0–3) | red point; needs all ℓ ≤ 19 | PY dashed (to 1.4 GeV, K K̄ removed) | vector | Im T peaks: T₀₊ ≈ 2.1 at 0.85 GeV (PY ≈ 2.1 at 0.77); T₀₀ ≈ 2.75 at ≈ 1.40 (PY 2.7 at 1.27); peaks shifted up by ≈ 80–130 MeV |
| 6 | S0ff.pdf (vector) | |F⁰₀| vs √s (0.28–2.0 GeV; 0–2.1) | blue/green/red | no | vector | broad peak at √s ≈ 0.48–0.50 GeV, heights ≈ 2.0 (blue), 1.81 (green), 1.68 (red); |F⁰₀(4m_π²)| ≈ 1.4; red has a second bump ≈ 0.96 at 1.53 GeV |
| 7 | rhopole.pdf (raster 3-D) | |S¹₁(s)| on the second sheet, axes √Re s 650–980 MeV, √Im s 0–600 MeV | – | no | raster, poor | **stale: identical file since v1** (not regenerated for v3 numbers) |
| 8 | P1FF.pdf (vector, 33 pts per colour) | |F¹₁| vs √s 0.28–1.2 GeV (0–8.5); dashed BW (5.2), solid GS (5.3) fits | blue, red | no | vector | data peaks ≈ 7.4 at 0.85 GeV (both); fitted peaks: blue ≈ 7.9 at 0.83, red ≈ 7.3 at 0.85 |
| 9 | eehadron.pdf (vector) | |F¹₁| GTB (blue, one fit) vs GS with m_ρ = 770, Γ_ρ = 149 MeV (black; caption says "PDG data") vs rescaled R(e⁺e⁻) (green, PDG) | blue | yes (PDG R-ratio) | vector | GTB peak ≈ 7.9 at 0.83 GeV vs GS ≈ 5.8 at 0.75 GeV; ω, φ spikes only in data |
| 10 | rho1450twopole.pdf (vector, 39 pts per colour) | |𝓕¹₁|² vs √s 0.28–2.0 (0–0.015) with two-BW fit (5.9) | blue, red | no | vector | peak 1: blue 0.0149 at 0.84, red 0.0131 at 0.85 GeV; peak 2 ≈ 0.0037 at ≈ 1.61–1.62 GeV |
| 11 | D0FF.pdf (vector) | |F⁰₂| vs √s 0.28–2.0 (0–8.5); BW fits blue/red; PDG f2 BW (black) | blue, red | PDG BW curve | vector | fit peaks: blue ≈ 8.35 at 1.365, red ≈ 7.35 at 1.385, PDG ≈ 6.8 at 1.275 GeV |
| 12 | Vertex0/1/2.png | Feynman vertices | – | – | n/a | – |
| 13 | FFdiag1/2.png | Brodsky–Lepage diagrams | – | – | n/a | – |
| 14 | FFdiag3.png | gluon contribution to the EMT form factor | – | – | n/a | – |
| 15 | jjdiag1.png | 2-pt function diagram | – | – | n/a | – |
| 16 | jjdiag2/3.png | four-quark current 2-pt diagrams | – | – | n/a | – |
| 17 | psandetaconverge.png (**raster**, 4997×6631) | as Fig. 4, green point with variations: base M=50, L=10(=ℓ≤19), N_B=40 (green); L=12 (orange); N_B=46 (purple); M=60 (olive); M_reg=5×10² (black) | 5 runs | same as Fig. 4 | raster only (high res) | all variants overlap within ≈ 5–10° except δ²₂, δ¹₃, and δ⁰₂ around 1.3–1.4 GeV (M=60 shifts later) |
| 18 | Fasymp00A/B, Fasymp11A/B, Fasymp02A/B (vector) | |F| and phase α for S0, P1, D0 vs √s (0.28–2.0/2.02) for r⁰₀ ∈ {2, 5, 7}×10⁻²; r¹₁ ∈ {1.5, 2, 2.5}; r⁰₂ ∈ {4, 6, 8} (blue = minimal feasible, orange = used) | 9 runs; point not stated | no | vector | |F¹₁| peak: 7.4 at 0.80 (r = 1.5), 6.9 at 0.85 (r = 2), 6.25 at 0.85–0.90 (r = 2.5), giving the quoted "∼80 MeV ρ-mass variation"; |F⁰₂| peak: 6.85 at 1.32 (r = 4), 5.85 at 1.32 (r = 6), ≈ 5.3 at 1.42 (r = 8). Orange r¹₁ = 2 peak (6.9) is lower than the Fig. 8 peaks, so these runs are probably pre-Watson or another point (unstated) |
| 19 | SRerr00/11/02 and nSR00/11/02 (vector) | |F| for S0, P1, D0 with ε⁰₀ ∈ {5e-8, 1e-7, 5e-7}, **ε¹₁ ∈ {1e-6, 6e-6, 1e-5}, ε⁰₂ ∈ {3e-6, 5e-6, 7e-6}**, and n_SR ∈ {3, 4, 5} | 15 runs | no | vector | n_SR curves coincide; ε¹₁ = 1e-6 moves the P1 peak to ≈ 0.80 GeV with height 6.1; legends confirm code tolerances |
| 20 | D0beforeW.pdf, D0afterW.pdf (vector) | δ⁰₂, η⁰₂ before/after Watson for blue/green/red (0.28–2.0 GeV) | 3 points × 2 stages, **in the released code** | gray data and PY | vector | before: δ⁰₂ rises early (blue ≈ 90° at 1.0 GeV, red at ≈ 1.3) with η⁰₂ down to 0.3–0.5; after: η = 1 and δ⁰₂ = 90° at ≈ 1.35–1.40 GeV |

v1 had Figs. 1–16 with the same roles, minus Figs. 17–20 (v2 added the convergence figure, v3 added 18–20). v1 also had `rho1450.pdf` (replaced in v2).

### Tables (verbatim from v4; blue, red)

**Table 1 — scattering lengths / effective ranges, units of m_π (label table:scatlength)**

| | W | GTB (blue, red) | CGL | PY |
|---|---|---|---|---|
| a⁽⁰⁾₀ | 0.16 | 0.174, 0.191 | 0.220 ± 0.005 | 0.230 ± 0.010 |
| a⁽²⁾₀ | −0.046 | −0.0357, −0.0403 | −0.0444 ± 0.0010 | −0.0422 ± 0.0022 |
| b⁽⁰⁾₀ | 0.18 | 0.285, 0.282 | 0.280 ± 0.001 | 0.268 ± 0.010 |
| b⁽²⁾₀ | −0.092 | −0.061, −0.071 | −0.080 ± 0.001 | −0.071 ± 0.004 |
| a⁽¹⁾₁ (×10⁻³) | 31 | 28.7, 30.7 | 37.0 ± 0.13 | 38.1 ± 1.4 |
| b⁽¹⁾₁ (×10⁻³) | 0 | 12.4, 7.28 | 5.67 ± 0.13 | 4.75 ± 0.16 |
| a⁽⁰⁾₂ (×10⁻⁴) | 0 | 14.6, 14.9 | 17.5 ± 0.3 | 18.0 ± 0.2 |
| a⁽²⁾₂ (×10⁻⁴) | 0 | 3.41, 2.53 | 1.70 ± 0.13 | 2.2 ± 0.2 |

Definition (4.1): Re f^I_ℓ(s) ≃ k^{2ℓ}(a + b k² + …), k = √(s − 4)/2. v1/v2 values: a⁰₀ 0.178/0.182, a²₀ −0.0369/−0.0378, b⁰₀ 0.287/0.290, b²₀ −0.064/−0.066, a¹₁ 28.0/28.4, b¹₁ 2.86/3.37, a⁰₂ 12.6/12.3, a²₂ 2.87/2.81.

**Table 2 — radii (fm²) (table:rSV)**: ⟨r²⟩_S = 0.64, 0.55 (exp. 0.61 ± 0.04); ⟨r²⟩_V = 0.394, 0.391 (0.439 ± 0.008); ⟨r²⟩_T = 0.153, 0.140 (no direct exp.; text quotes PDF-based 0.67–0.77 and 2.1–2.4 fm²). Definition (4.2)–(4.4). v1/v2: S 0.64/0.61, V 0.388/0.381, no T.

**Table 3 — χPT constants (table:lj)** (from Tables 1–2 via (4.5)–(4.8))

| | GTB (blue, red) | GL | Bij | CGL |
|---|---|---|---|---|
| l̄₁ | 0.96, −0.10 | −2.3 ± 3.7 | −1.7 ± 1.0 | −0.4 ± 0.6 |
| l̄₂ | 4.5, 4.8 | 6.0 ± 1.3 | 6.1 ± 0.5 | 4.3 ± 0.1 |
| l̄₄ | 4.8, 4.2 | 4.3 ± 0.9 | 4.4 ± 0.3 | 4.4 ± 0.2 |
| l̄₆ | 14.6, 14.5 | 16.5 ± 1.1 | 16.0 ± 0.5 ± 0.7 | |

I re-derived these from Tables 1–2 (f_π = 92 MeV, m_π = 139.57 MeV, ħc = 197.327 MeV·fm) and got l̄₁ = 0.96/−0.12, l̄₂ = 4.49/4.83, l̄₄ = 4.75/4.23, l̄₆ = 14.5/14.4, consistent. v1/v2: 0.92/0.93, 4.1/4.0, 4.7/4.6, 14.3/14.1.

**Table 4 — ρ T-matrix pole (MeV) (table:rhopole)**: Re√s_ρ = 829, 844 (PDG 761–765 ± 0.23); Im√s_ρ = 71, 82 (PDG 71–74 ± 0.8). Method: |S¹₁| continued with S(s)S*(s̄) = 1 (5.1). v1/v2: 829/832, 63/64.

**Table 5 — ρ mass/width from FF fit (MeV) (table:mrho)**: m_ρ = 838, 859 (PDG 775 ± 0.23); Γ_ρ = 114, 128 (149.1 ± 0.8). The text does not say whether these come from the BW (5.2) or the GS (5.3) fit; GS is the one it calls better. v1/v2: 836/839, 111/111.

**Table 6 — two-BW fit to |𝓕¹₁|² up to 2 GeV (table:1450)**: m₁ = 824, 841; Γ₁ = 137, 164; m₂ = 1635, 1645 (PDG 1465 ± 25); Γ₂ = 248, 272 (400 ± 60); A₁ = 1.29, 1.39; A₂ = 0.27, 0.28. Also (5.10) A₁ − A₂ = 1.03, 1.11. v2: 825/828, 128/126, 1632/1636, 180/210, 1.23/1.21, 0.19/0.20; A₁ − A₂ = 1.04/1.01.

**Table 7 — f2 from BW fit of |F⁰₂| (MeV) (table:mf2)**: m_f2 = 1366, 1388 (PDG 1275.4 ± 0.6); Γ_f2 = 164, 188 (186.6 ± 2.3). v1/v2: 1165/1283, 151/201.

**Table 8 — numerical parameters (numpara)**: Discretization M = 50, L = 19, N_B = 40; Regularization M_reg = 10²; χSB (2.15) ε^χ = 2×10⁻³; FF asymptotics (2.21i), (2.22i), (2.23i) r^I_ℓ (8.75); SVZ (2.21h), (2.22h), (2.23h) n_SR = 3, ε^I_ℓ (8.73); UV scale s0 = 2 GeV.

**Table 9 — GTB vs GTBII (scattering lengths recomputed from Table 3 via (A.1)–(A.8) with l̄₃ = 0) (table:GTBII)**

| | GTB | GTBII | PY |
|---|---|---|---|
| a⁽⁰⁾₀ | 0.174, 0.191 | 0.21, 0.21 | 0.230 ± 0.010 |
| a⁽²⁾₀ | −0.0357, −0.0403 | −0.043, −0.043 | −0.0422 ± 0.0022 |
| b⁽⁰⁾₀ | 0.285, 0.282 | 0.26, 0.26 | 0.268 ± 0.010 |
| b⁽²⁾₀ | −0.061, −0.071 | −0.080, −0.079 | −0.071 ± 0.004 |
| a⁽¹⁾₁ | 28.7, 30.7 | 36.3, 36.7 | 38.1 ± 1.4 ×10⁻³ |
| b⁽¹⁾₁ | 12.4, 7.28 | 2.62 3.43 | 4.75 ± 0.16 ×10⁻³ |
| a⁽⁰⁾₂ | 14.6, 14.9 | 14.8, 15.0 | 18.0 ± 0.2 ×10⁻⁴ |
| a⁽²⁾₂ | 3.41, 2.53 | 3.4, 2.5 | 2.2 ± 0.2 ×10⁻⁴ |

Inline numbers: g_ρππ = 4.9, 5.2 (5.6) from (5.5) (PDG-based 6). I verified 4.95/5.15 from Table 5. 4π/g_ρ² ≃ 0.61, 0.57 vs 0.507 ± 0.011 (after (5.8); the ρ-peak integration window is unspecified). g_f2ππ = 8.5, 8.9 (5.13) vs 9.26 (EWERZ); I verified 8.48 from Table 7, which assumes 100 % ππ.

---------------------------------------------------------------------------------------------------

## 4. Quantitative claims (pre-registration list)

Path classes:
- **R** = output of the released path (one fixed-𝔉₀ maximisation + one Watson step) for ni = 44/46/48 (three runs).
- **R+A** = released path plus analysis not shipped (threshold fits, FF derivative, BW/GS fits, pole continuation, Legendre sums).
- **N** = needs runs not in the released code (scans, parameter variations, ablations, other s0).

Values are (blue, red) unless noted. "PDG/exp" are the paper's comparators, not targets.

| ID | Claim | Quoted value | Location | Class |
|---|---|---|---|---|
| C1 | 𝔉₁^max at fixed 𝔉₀ = (44, 46, 48)/50 × 0.163108, before Watson | −0.02723, −0.02877, −0.03040 (vector read-off ±3e-5) | Fig. 3 | R (first SDP only) |
| C2 | Allowed region in (𝔉₀, 𝔉₁): upper and lower edges, tip at 𝔉₀ ≈ 0.163, 𝔉₁ ≈ −0.0324; tree-level dot (0.14575, −0.02915) is inside | 100 boundary points | Fig. 3 | N (50 × max/min) |
| C3 | ε^χ = 2e-3 is the minimal tolerance keeping the tree-level point allowed | 2×10⁻³ | after (8.69), §8.4 summary | N (feasibility scan) |
| C4 | Six phase shifts and inelasticities in "reasonable agreement" with data up to 2 GeV except S0 K K̄ region; η = 1 after Watson | qualitative, Fig. 4 | §3, Fig. 4 | R |
| C5 | ρ(770) in P1 at "near its physical value"; δ¹₁ crosses 90° at ≈ 0.83–0.85 GeV | eye-read | §3, Fig. 4 | R |
| C6 | f2(1270) appears in D0 ("strong resonance near 1.2 GeV"); δ⁰₂ = 90° at ≈ 1.35–1.40 GeV | eye-read | §3, §5.3, Fig. 4 | R |
| C7 | Second P1 phase rise at ≃ 1.6 GeV, identified with ρ(1450) | ≃ 1.6 GeV | §3, §5.2.3 | R |
| C8 | Scattering lengths/effective ranges | Table 1 (8 rows × 2) | Table 1 | R+A |
| C9 | Radii ⟨r²⟩_S, ⟨r²⟩_V, ⟨r²⟩_T | 0.64/0.55; 0.394/0.391; 0.153/0.140 fm² | Table 2 | R+A |
| C10 | l̄₁, l̄₂, l̄₄, l̄₆ | 0.96/−0.10; 4.5/4.8; 4.8/4.2; 14.6/14.5 | Table 3, (4.5)–(4.8) | R+A (algebra on C8, C9) |
| C11 | GTBII recomputation agrees "reasonably" | Table 9 | App. A | R+A (algebra) |
| C12 | Scalar FF has a broad peak at ∼ 500 MeV (σ-like) | ≈ 0.48–0.50 GeV; heights 2.0/1.68 | §5.1, Fig. 6 | R |
| C13 | ρ pole on sheet II | Re√s = 829/844, Im√s = 71/82 MeV | Table 4, (5.1) | R+A (complex continuation not in code) |
| C14 | ρ mass/width from FF fit | 838/859, 114/128 MeV | Table 5, Fig. 8 | R+A |
| C15 | g_ρππ | 4.9/5.2 | (5.6) | R+A |
| C16 | 4π/g_ρ² | 0.61/0.57 | after (5.8) | R+A (window unspecified) |
| C17 | Two-BW fit: ρ and ρ′ | m₁ 824/841, Γ₁ 137/164, m₂ 1635/1645, Γ₂ 248/272, A₁ 1.29/1.39, A₂ 0.27/0.28, A₁ − A₂ 1.03/1.11 | Table 6, (5.10) | R+A |
| C18 | f2 BW fit | 1366/1388, 164/188 MeV | Table 7, Fig. 11 | R+A |
| C19 | g_f2ππ | 8.5/8.9 | (5.13) | R+A |
| C20 | Forward amplitudes (red point) agree with PY "reasonably", with shifted peaks | Fig. 5 | §3.1 | R+A (needs ℓ ≤ 19 sums) |
| C21 | Watson step greatly improves D0: before, η⁰₂ ≈ 0.3–0.6 and early δ rise; after, η = 1 and f2 shape | Fig. 20 | App. B.4 | R (both stages in code) |
| C22 | Discretisation stability at green point: L = 23 (l = 12), N_B = 46, M = 60, M_reg = 5×10² change results little | Fig. 17 | App. B.1 | N |
| C23 | Minimal feasible caps r⁰₀ = 2e-2, r¹₁ = 1.5, r⁰₂ = 4; r = 1 (P1, D0) infeasible | Fig. 18, after (8.75) | App. B.2, §8.4 | N (feasibility bisection) |
| C24 | ρ mass varies by ∼ 80 MeV for r¹₁ = 1.5–2.5 | ∼ 80 MeV | App. B.2 | N |
| C25 | n_SR = 3, 4, 5 give the same |F|; small ε changes have little effect | Fig. 19 | App. B.3 | N (n_SR > 3 needs new NB moments) |
| C26 | s0 independence: below 1.2 GeV the phase shifts agree with the s0 = 1.2 GeV results of 2309.12402; ρ "largely unchanged" | qualitative | §2.3, §5.2, §9 | N (cross-paper; 2309 runs) |
| C27 | Ablation: removing P1 FESR removes the ρ; removing D0 FESR does not completely remove the f2 | qualitative | §1 | N |
| C28 | Uncertainty model = blue–red spread ("main source of error is identification of the point") | – | §2.4, §4 | definitional |
| C29 | Runtime ∼ 20 min laptop (Matlab) | – | §1 footnote | R (informational) |

The released code computes one point per run (default ni = 46, green). The blue/red numbers in Tables 1–7 need ni = 44 and 48, which is a one-line edit. Nothing in Tables 1–7 is computed in the released code beyond the raw arrays (ht, FF, spec).

---------------------------------------------------------------------------------------------------

## 5. Context: what 2505.19332 (v3, 2026-02-04; PTEP 2026 063B05) changes, and what the main-branch code covers

**Method changes vs 2403**
- **Spectral densities:** B-splines are replaced by a (1−z)-power ("atorho") basis that builds in κ_ℓ threshold behaviour, with leading coefficient A, B, C ∈ [1, 10] and bounded higher coefficients |a_n| ≤ 10, |b_n| ≤ 100, |c_n| ≤ 1000 (App. B of 2505; `aS0b=10, aP1b=100, aD0b=1000`, Na = 49).
- **Sum rules:** "pinched" kernel K(s) = sⁿ(1 − s/s0)e^{−s/s0} for S0 and P1 with 3-loop (α_s²) coefficients, n_SR = 6; D0 keeps sⁿ weights (App. C, Tables of c_{i,j}). Code tolerances eS0 = 0.05, eP1 = 0.1, eD0 = 0.1, applied componentwise (abs). These look like relative tolerances, but the paper does not state them.
- **FF asymptotics:** Brodsky–Lepage bounds |F| ≤ χ·BL/s with χ⁰₀ = 8 (on m_π⁴ ln(s)/s), χ¹₁ = 2.2, χ⁰₂ = 9 (eq. chichoice; "close to minimal feasible"). Main code: rS0 = 8, rP1 = 2.2, **rD0 = 8** (≠ 9). Code α_s = 0.327631 = one-loop α_s at s_{n0} (1.83 GeV), versus 0.31459 at 2 GeV in 2403. The bounds start at index n0 (the last point below s0) rather than n0 + 1. s0 itself is still (2 GeV)², with M = 50, L = 19, ν₀ = −20 and the same M_reg normalisation.
- **Low energy:** only the two ratio constraints (no D/F-wave vanishing constraint), ε^χ = 2e-3 (scanned in Fig. 13).
- **Selection:** no boundary point and no 𝔉₀/𝔉₁ functional. Pure feasibility (`minimize(0)`) is followed by iterated Watson; the Watson target is |h_old|·F_old/|F_old| (modulus of the old wave, phase of the FF). The paper's figures use 10 iterations; the main code has `nWatson = 5`.
- **Diagnostic:** λ = (3π/4)A(4/3, 4/3, 4/3) ("pion coupling") is printed every iteration; the paper also maximises and minimises λ, a's and r's.

**2505 headline numbers (converged solution):** λ = 0.02 (initial spread 0.01–0.022); f_π = 101 MeV inferred; l̄₁ … l̄₆ = 1.6, 5.5, 7.8 (l̄₃), 4.7, 14.3; a⁰₀ = 0.172, a²₀ = −0.0366, b⁰₀ = 0.251, b²₀ = −0.063, a¹₁ = 32.1e-3, b¹₁ = 2.69e-3, a⁰₂ = 12.9e-4, a²₂ = −1.1e-4; ⟨r²⟩_S, V, T = 0.55, 0.441, 0.146 fm²; m_ρ = 758, Γ_ρ = 137, m_ρ′ = 1514, Γ_ρ′ = 162, m_f2 = 1180, Γ_f2 = 249 MeV.

**Main-branch code** (`theories/qcd/src/matlab/gtb_qcd_02_optimize.m` + `optimize_core.m`; identical to the 2505 snapshot; Mathematica `gtb_qcd_01_generate_matrices.nb`; `gtb_qcd_03_analyze_plot.nb` only plots and contains no fits):
- **Reproduces (modulo nWatson = 5 vs 10 and χ⁰₂ = 8 vs 9):** the "no functional" chain, giving Fig. 6 (pwplot/ffplot), the green-triangle series of Fig. 2 (λ per iteration), and the raw arrays behind Tables (masses, LECs, a, r).
- **Not reproduced without modification:** the max/min λ, a, r chains (Figs. 2–5 red/blue; the objective must change); the χ scans (Figs. 9–12); the ε^χ scan (Fig. 13); the χPT one-loop overlay (Fig. 7); the R-ratio comparison (Fig. 8); thermodynamics (Figs. 14–15); all BW fits and LEC extraction (no analysis code).
- 2505 Fig. 1 is a schematic.

---------------------------------------------------------------------------------------------------

## 6. Statements on precision, solver, instabilities, regulariser, convergence

**2403 (v4 unless noted)**
- Solver: CVX with MOSEK "for fast numerical optimizations" (§8 footnote); no solver tolerances or precision settings are given. CODE uses defaults (`cvx_begin sdp`, no `cvx_precision`).
- Conditioning: kernels Ã, B̃, B̂ are computed "with high precision before importing them into the optimization program and multiplying them with the machine precision variables" (§8.4, after (8.64)). The rescaled matrix (8.61) has "order 1" entries and is "better scaled ... for machine precision semi-definite programming". Λ_ℓ is taken from the pure S-matrix work He:2021eqn. v1 wording: the rescaling is meant to avoid "potential numerical instability".
- Regulariser: ‖[ρ₁, ρ₂]‖₄ ≤ M_reg = 10² "avoids unphysical highly oscillating double spectral densities in the primal problem and renders the dual problem feasible"; "value tested with the pure S-matrix bootstrap to ensure convergence"; M_reg = 5×10² checked in Fig. 17 (§8.4 after (8.56); App. B.1). v1: "For different numerical problems, one would need to test this". The code applies it to ρ/Mrho (mismatch 1).
- FF caps: order-one relaxation factors are "necessary since otherwise the problem becomes unfeasible" (§8.4). r⁰₀ is "close to the smallest one that makes the problem feasible", although Fig. 18 shows 2e-2 is still feasible. v1/v2: the D0 Brodsky–Lepage estimate made the program infeasible, so the D0 cap was "gradually raised until the program just becomes feasible".
- Instabilities (v1/v2 only, removed in v3): lower-edge points "display numerical instabilities where the sum rules are satisfied by discontinuous jumps near threshold". v2 App. B: with M = 60 the f2 rise shifts because n0 moves and the FF bounds change.
- Convergence in M, L, N_B, M_reg: "good overall stability" for M = 60, L = 23, N_B = 46, M_reg = 5×10² at the green point (App. B.1, Fig. 17). ν₀ convergence "tested for a few values" (footnote to (8.28)). No convergence in M beyond 60 and no error bars.
- s0: "going to very large energies is challenging numerically due to the larger number of interpolation points required" (footnote §2.3). s0 = 3 GeV "will require new improved techniques that we have not found yet" (§8.4 summary). The B-spline and Watson techniques were "required" for s0 = 2 GeV.
- Watson: "can be repeated until convergence" but only one step is used. v4 footnote: iterating converges to a solution "largely independent of the functional" (from 2505).
- FESR: integrals truncated at s_{n0} ≈ 1.83 GeV rather than s0 (footnote after (8.71)); ε set from the O(α_s) term.

**2505**
- No solver or precision statements. χ = 1 is infeasible ("not in the fully asymptotic region"). Coefficient regulators a_max, b_max, c_max are "chosen in a range that does not affect the results" (values only in code). The iteration "converges to a unique solution" from max, min or feasible starts (Figs. 2–5). "In all cases we always do one initial Watsonian iteration to smoothen out the functions" (footnote §3.1). The Watson iteration removes the ε^χ dependence except for I = 2 and F1 (Fig. 13).

---------------------------------------------------------------------------------------------------

## 7. Gaps between the released code and what the 2403 figures need (summary)

1. **Inputs not shipped.** Running needs Mathematica (NB → 17 `.mat`) plus MATLAB + CVX + MOSEK. Under this project's SDPB-only constraint, NB and the SOCP/SDP must be re-implemented. Tolerances are absolute and tiny (ε⁰₀ = 1e-7 on moments ∼ 3e-7, |F⁰₀| caps 0.05), so scaling matters for SDPB.
2. **One point per run.** Default ni = 46 (green); blue/red need ni = 44/48. CODE-v1 instead used 𝔉₀ = (49, 50)/60 × 0.1776.
3. **Fig. 3** needs a 50-point 𝔉₀ scan with a max **and min** branch. There is no minimise branch in the code, and these values are pre-Watson.
4. **No analysis code at all.** Missing pieces:
   - threshold expansion (Table 1)
   - FF slope at s = 0 (Table 2)
   - l̄ algebra (Table 3)
   - second-sheet continuation of S¹₁ at complex s, which needs kernels at complex s that NB does not export (Table 4, Fig. 7)
   - BW/GS and two-BW fits, including the unspecified fit windows (Tables 5–7, Figs. 8, 10, 11)
   - g_ρ peak integral
   - forward amplitudes from all ℓ ≤ 19 (Fig. 5)
   - η plotting
   - all experimental overlays: Protopopescu/Losty/Hyams data, PY curves with K K̄ removed, PDG R-ratio, PDG f2 BW
5. **App. B variations need regenerated kernels** (M = 60, L = 23, N_B = 46, n_SR = 4, 5). Only M_reg, r and ε are one-line edits. The point used in Figs. 18–19 is unstated.
6. **Fix list before trusting text values:** SR ε (P1 6e-6, D0 5e-6), M_reg normalisation (/Mrho), second chiral norm on F1 (not P1), Watson target factor 2, α_s = 0.3145887, m_u = 3.6 MeV, m_d = 6.5 MeV.
7. **Fig. 7 is stale since v1** (it does not match Table 4). There is no pole figure to compare against for v3/v4 numbers.
8. **Target version:** v4 text with the v3/v4 ancillary code (= GitHub snapshot). There is no journal version. PoS(CD2024)040 repeats the v3/v4 numbers.
