# Gate log (trial-and-error record, kept verbatim)

- 11:53 launched Mreg = 1e3, 1e4, 1e5, 1e6 with the first encoding of the regulariser
  blocks: [[Mreg*y0 - rho_i]], [[Mreg*y0 + rho_i]] (constant entry = Mreg).
- 12:10 Mreg=1e6 terminated after 15 iterations with `maxComplementarity exceeded`:
  primal objective grew 1e14 -> 1e100 while the dual objective stayed ~3.6e3 and the
  dual error at 9.9e5.  Diagnosis: SDPB initialises every block matrix at
  initialMatrixScale = 1e4, so the 7550 regulariser blocks contribute
  7550 * 1e4 * Mreg = 7.5e13 to the initial primal objective at Mreg=1e6; the step
  control never recovered.  This is an encoding/scaling artefact, not a property of
  the regularised problem.
- Fix: divide each regulariser block by Mreg -> [[y0 -/+ rho_i/Mreg]] (constant entry 1,
  same inequality).  Tests updated (7 pass).  Mreg=1e5 and 1e6 relaunched with the
  unit-scaled encoding as `gate_tip_linf_1e5_unit`, `gate_tip_linf_1e6_unit`.
  Mreg=1e3 and 1e4 (constant entries 1e3/1e4) were left running; their encoding is
  recorded in each report's pmp block sizes and in git history (commit 51ce1e1 vs the
  unit-scaled commit that follows).  A plateau statement will only compare runs whose
  acceptance passed the full Arb replay, regardless of encoding.
- 12:36 Mreg=1e3 (raw encoding) accepted: tip 0.0839635, rho_linf 999.997 (bound active), ||rho||_4 5.80e3, SDPB gap < 1e-6, full Arb replay passed.
       Mreg=1e4 (raw encoding) accepted: tip 0.0885125, rho_linf 9999.995 (active), ||rho||_4 5.79e4.
  Paper 0.08257.  1e3 is inside the pre-registered band (+1.7 %), 1e4 is outside (+7.2 %); the two differ by 4.6e-3, so no plateau
  between these decades.  Launched Mreg = 1e2 and 3e2 (unit-scaled) to map the low side; 1e5/1e6 unit-scaled still running.
- 12:44 Diagnostics of the two accepted solutions (independent float re-projection of the saved amplitude):
    Mreg=1e3: T0=-10.1, |sigma|max 3.8e3, 18.3 % of the 3775 double-density values sit on the bound,
              omitted waves (I0/I2 l=20,22; I1 l=21,23) worst eta = 1.03-1.12, all at the last node (5.94 GeV);
              tip phases: S0 56 deg at 0.9 GeV, P1 no 90-deg crossing below 1.2 GeV (chiral-only tip, as expected).
    Mreg=1e4: T0=15.3, |sigma|max 4.6e4, 17.9 % on the bound, omitted-wave worst eta = 1.14-1.46.
    Without any regulariser the same omitted waves had eta = 53 (earlier SDPB runs).
  Observation to be confirmed: the Mreg=1e4 solution has ||rho||_4 = 5.8e4 < 3.775e5, i.e. it is feasible for the
  authors' released l4 bound, yet its tip 0.0885 exceeds the paper's 0.0826 by 7 %.  If confirmed by an explicit
  l4-bounded run, the exact optimum of the authors' stated problem lies above their plotted endpoint, and the plotted
  value reflects the double-precision solver they used (the repo's own Clarabel ladder reproduced 0.0826 with
  ||rho||_4 ~ 6e3-1.5e4).

- 12:46 PRE-REGISTERED fallback rule, written before the Mreg=1e2/3e2/1e5/1e6 results are known:
  (1) If the tip agrees within 1e-3 over two consecutive decades of Mreg, that is the plateau (2103 sec. 3.3) and the
      smallest Mreg on it is Mreg*.
  (2) If no such plateau exists, Mreg* is the LARGEST scanned Mreg for which the regularised amplitude satisfies the
      unitarity that the truncation dropped: max eta <= 1.02 over all nodes for the six first-omitted waves
      (I0/I2 l=20,22; I1 l=21,23), AND the L=8/10/12 tips at that Mreg agree within 2 % (the paper's Appendix-A
      stability).  Both are physics conditions independent of Fig. 4.  The full tip-vs-Mreg curve is reported either way.
  (3) The paper's endpoint band (0.0826 +- 5 %) is then a *test* of the chosen Mreg*, never an input to it.
- 13:16 unit-scaled runs accepted: Mreg=1e5 tip 0.1036422, ||rho||_4 5.70e5; Mreg=1e6 tip 0.1464423, ||rho||_4 5.61e6.
  Series (Mreg -> tip): 1e3 0.08396, 1e4 0.08851, 1e5 0.10364, 1e6 0.14644: monotone, accelerating (+5 %, +17 %, +41 % per decade).
  No plateau on the high side; rule (1) cannot apply there.  The authors' l4 bound 3.775e5 is crossed between the 1e4 and 1e5
  solutions, so the exact optimum of their stated problem lies between 0.0885 and ~0.10 (log-interpolation ~0.10), i.e. 7-20 %
  above their plotted 0.0826.  Decision now depends on the low side (1e2, 3e2, still solving) via rule (2).
- 13:19 Norm-independence control started at M=30/L=8 (the ladder resolution that reproduced 0.0826 in double precision):
  linf at Mreg=1e2 and 1e3 first; l2/l4 follow once the linf scale that satisfies the omitted-wave rule is known.
  Mathematica formula audit skipped for these control runs (--skip-mma-audit): the audit is M-independent in content and
  passed at M=45/50; the control compares norms against each other at fixed M.
- 13:32 M=30/L=8 linf controls accepted:
    Mreg=1e2: tip 0.0826714 (paper 0.0825728, +0.12 %), ||rho||_2 2619, ||rho||_4 488, 30 % of values on the bound,
              omitted waves (l=16..21) max eta 1.0072  -> satisfies the pre-registered eta <= 1.02.
    Mreg=1e3: tip 0.0843413 (+2.1 %), ||rho||_4 4801, omitted-wave max eta 1.107 -> fails eta <= 1.02.
  Rule (2) therefore selects Mreg = 1e2 at M=30 (largest scanned value satisfying the omitted-wave condition), and the
  paper's endpoint is reproduced to 0.1 % without using it as input.  Independent corroboration: the double-precision
  Clarabel ladder at M=30 gave 0.082598 with ||rho||_4 of the same order.
- 13:33 Norm-independence control launched at M=30/L=8 with anchors taken from the linf Mreg=1e2 solution:
  l2 with B2 = 2600 (~||rho||_2 of that solution), l4 with B4 = 500 (~||rho||_4); plus linf Mreg=3e1 to map the low side.
  Each norm's scale is then chosen by the same omitted-wave rule; a x3 variant per norm follows to confirm it is the largest.
- 13:41 M=30/L=8 linf Mreg=3e1 accepted: tip 0.0820470 (paper -0.64 %), ||rho||_4 148.  Low-side series at M=30:
  3e1 0.08205, 1e2 0.08267, 1e3 0.08434 -> the change per half-decade on the low side is 6e-4, i.e. the curve flattens
  below Mreg ~ 1e2 while the omitted-wave condition is met; rule (2) selection Mreg=1e2 stands.
- 13:50 M=50 low side accepted: Mreg=1e2 tip 0.0827570 (paper +0.22 %), ||rho||_4 598, six-wave omitted eta 1.0106 (<= 1.02);
  Mreg=3e2 tip 0.0832653 (+0.84 %), ||rho||_4 1745, eta 1.046 (> 1.02).
  Full M=50 series: 1e2 0.08276 | 3e2 0.08327 | 1e3 0.08396 | 1e4 0.08851 | 1e5 0.10364 | 1e6 0.14644.
  Rule (2), first condition: Mreg* = 1e2 at M=50 (recorded in MREG_STAR_M50.json).  Second condition (L=8/10/12 tips
  within 2 %) now being tested: launched M50/L8 and M50/L12 at Mreg=1e2.
- 13:55 Consolidated table (tip = +x end of the eps=0.002 chiral region; phases at the tip, 0.9 GeV, degrees):
    run           tip      vs paper  ||rho||4  eta_omit  S0    S2    P1
    M50 L10 1e2   0.08276  +0.22 %   598       1.011     59.8  -26.4  33.1
    M50 L10 3e2   0.08327  +0.84 %   1.7e3     1.050     57.5  -27.1  32.4
    M50 L10 1e3   0.08396  +1.68 %   5.8e3     1.145     55.7  -27.7  32.0
    M50 L10 1e4   0.08851  +7.19 %   5.8e4     1.637     49.8  -31.1  30.4
    M50 L10 1e5   0.10364  +25.5 %   5.7e5     1.703     44.1  -35.2  31.8
    M50 L10 1e6   0.14644  +77.4 %   5.6e6     3.049     34.7  -39.5  47.3
    M30 L8  3e1   0.08205  -0.64 %   148       1.006     62.1  -25.6  35.3
    M30 L8  1e2   0.08267  +0.12 %   488       1.007     59.8  -26.1  34.1
    M30 L8  1e3   0.08434  +2.14 %   4.8e3     1.107     57.5  -27.7  33.3
  At the physics-selected scale Mreg=1e2, M=30 and M=50 agree to 0.1 % in the tip and to ~1 degree in all three
  tip phase shifts (M-stability); at Mreg >= 1e4 the tip and the phases drift strongly (S0 at 0.9 GeV: 60 -> 35 deg).
- 14:05 M=30/L=8 l2 control, B2=2600 (anchored on ||rho||_2 of the linf Mreg=1e2 solution): accepted, tip 0.0830076
  (paper +0.53 %, linf-selected value +0.41 %), ||rho||_2 2599 (active), ||rho||_4 875.  Omitted-wave eta recorded in
  omitted_waves.json; l2 at B2=8000 launched to confirm whether 2600 is the largest scale satisfying eta <= 1.02.
- 14:10 l2 B2=2600: six-wave omitted eta = 1.0211, marginally above the 1.02 threshold -> not the l2 physics scale.
  Stopped the B2=8000 run (a larger scale cannot pass) and launched B2=1800 as the next candidate below.
- 14:40 M=50/L=8, Mreg=1e2 accepted: tip 0.0828558 (L=10: 0.0827570, +0.12 %).  L=12 at iteration 58 with dual objective
  0.08250 (gap 1.9e-3), i.e. heading to about -0.3 % of L=10.  Rule (2) second condition (<= 2 %) will be met by a wide
  margin once L=12 is accepted.
- 14:55 M=50/L=12 Mreg=1e2 accepted: tip 0.0827039 (L=10 -0.06 %).  l2 B2=1800 accepted: tip 0.0828355, six-wave eta 1.0138
  (passes) -> l2 physics scale = 1800; linf vs l2 tips differ by 0.20 %, tip phases by < 0.5 degree.
  GATE PASSED under rule (2); verdict in GATE_RESULT.md.  Phase 2 started: eps=0.002 x_ref sections (support reuse of the
  accepted tip PMP) and the eps=0.004 / 0.006 tips at Mreg=1e2.
- 15:00 l4 control (M=30/L=8, B4=500, anchored on ||rho||_4 of the linf Mreg=1e2 solution) accepted: tip 0.0828383,
  ||rho||_4 = 499.8 (active).  See the norm table appended below; the omitted-wave eta of this run is in its
  omitted_waves.json.  Norm-independence: l-inf 0.082671 / l2 0.082836 / l4 0.082838, spread 0.20 %, tip phases at
  0.9 GeV within 0.5 degree -> the pre-registered <= 2 % / <= 2 degree criterion is met; the l-inf form is validated
  as the mainline regulariser.  Phase 0b closed.
- 15:58 Phase 2, eps=0.002 x_ref sections accepted (support reuse of the Mreg=1e2 tip PMP):
    upper y = -0.0043725 (digitised Fig.4/8 green upper -0.0043678, +0.11 %), lower y = -0.0051343 (paper -0.0051276, +0.13 %),
    width 0.0007618 (paper 0.0007598, +0.26 %).  Omitted-wave eta: upper 1.028, lower 1.003.
  Preview of the Fig.7 comparison at the upper x_ref point (the frozen rule of 15:10 restricts the representative to the
  upper branch and the +-0.001/0.002 fan cannot beat the x_ref point, whose 2-D distance to P is 0.00052 < 0.001):
    S0 RMS 0.4 deg (endpoint +1.0), S2 RMS 0.1 deg (+0.7), P1 RMS 0.5 deg (+2.0) over 0.28-1.196 GeV; no P1 90-degree
    crossing below 1.2 GeV.  The lower-branch point is qualitatively different (S0 38 deg at 0.9 GeV), as the paper's
    choice of the upper boundary anticipates.  Formal C4 verdict waits for the fan receipt; C3 evaluation next.
- 16:05 eps ladder tips accepted at Mreg=1e2: eps=0.004 -> 0.1255315 (Clarabel M30 ladder 0.125637, -0.08 %);
  eps=0.006 -> 0.1607729 (ladder 0.161054, -0.17 %).  C3 at eps=0.002 (upper x_ref point, 131-point subthreshold grid):
  RMS vs the digitised Fig.5 green curves 0.96 % (S0), 0.27 % (S2), 0.19 % (P1) of f00(3) [budget 8 %]; S0 chiral zero at
  s=0.4262 (paper ~0.425).  Queued the x_ref upper sections for eps=0.004 and 0.006 (Fig.5 orange/blue) ahead of the rest.
- 16:20 eps=0.001 tip accepted at Mreg=1e2: 0.0550677 (Clarabel M30 ladder 0.055038, +0.05 %).  Ladder so far, monotone:
  0.006 0.16077 | 0.004 0.12553 | 0.002 0.08276 | 0.001 0.05507; 0.0006 and 0.0002 queued.
- 16:35 fan: x_ref-0.001 upper y = -0.0042911 (distance to P 0.001165 vs 0.000516 for the x_ref point).
- 16:48 fan: x_ref+0.001 upper y = -0.0044553 (distance to P 0.001090).
- 17:02 eps=0.004 / 0.006 upper x_ref sections accepted: y = -0.0035180 / -0.0027953.  C3 (Fig.5) evaluated on the three
  upper x_ref points: RMS/f00(3) S0 0.96 / 2.61 / 6.43 %, S2 0.27 / 0.75 / 1.93 %, P1 0.19 / 0.51 / 1.28 % (budget 8 %);
  S0 chiral zero 0.4262 / 0.2929 / none (paper 0.425 / 0.305 / none).  C3 PASS -> C3_RESULT.json.
- 17:20 fan: x_ref-0.002 upper y = -0.0042108 (distance to P 0.002112).
- 18:20 fan complete (x_ref+0.002 upper y = -0.0045398, distance 0.00203).  Frozen rule selects the x_ref upper point
  (distance 0.000516); near points = x_ref +- 0.001.  IR_REPRESENTATIVE_SELECTION.json written before any further phase use.
  C4 formal: S0 RMS 0.4 / S2 0.1 / P1 0.5 deg vs digitised Fig.7, endpoints +1.0 / +0.7 / +2.0 deg, no P1 crossing
  below 1.2 GeV -> PASS (C4_RESULT.json).
- 18:35 eps=0.0006 tip 0.0410551, eps=0.0002 tip 0.0222912 accepted.  Six-epsilon ladder complete and monotone:
  0.006 0.16077 | 0.004 0.12553 | 0.002 0.08276 | 0.001 0.05507 | 0.0006 0.04106 | 0.0002 0.02229.
  C2 formal verdict PASS (C2_RESULT.json): +x end +0.22 % (tol 5 %), x_ref width +0.26 % (tol 20 %), monotone.
- 18:55 chi-c control (two separate 4-dim L2 balls, eps=0.002, Mreg=1e2): tip 0.0860468, i.e. +4.0 % over chi-b
  (0.0827570) and +4.2 % over the paper.  The Newton-era +28 % (0.1059) therefore was not caused by the norm grouping
  alone; chi-b stays the mainline as pre-registered.
- 19:10 Fig.3 pure-unitarity tip at Mreg=1e3 accepted: 2.2244538 (digitised 2.2328885, -0.38 %).  Mreg=1e4 running, 1e2 queued;
  the omitted-wave rule decides the pure-problem scale separately (amplitude scale ~2 vs ~0.08 for the chiral region).
- 19:35 Fig.3 pure tip at Mreg=1e2 accepted: 2.2115961 (-0.95 % vs 2.2328885), ||rho||_4 658.  Mreg=1e3 gave 2.2244538 (-0.38 %),
  omitted eta 1.017.  Waiting for Mreg=1e4 before applying rule (2) to the pure problem.
- 19:50 Fig.3 pure Mreg scan complete: 1e2 2.2116 (eta 1.005) | 1e3 2.2245 (eta 1.016) | 1e4 2.2305 (eta 1.138).
  Rule (2) -> pure Mreg* = 1e3 (MREG_STAR_PURE.json); tip -0.38 % vs the digitised 2.2328885 (C1 tolerance 2 %).
  24-direction support sweep (15-degree steps) queued from the Mreg=1e3 tip PMP.
- 19:55 PRE-REGISTERED before the UV tips are accepted: the UV +x ends are heading to ~0.075 (SR-a) and ~0.072 (SR-d),
  below the digitised cyan tip 0.0811 (chiral tip 0.0828 -> paper shrinks by 2 %, we by 9-13 %).  The paper's text on the
  form-factor bound (3.75) is ambiguous: "(3.75) ||F(s_i)||^2 <= ... for s_i > s0" with per-node kinematic factors, versus
  the sentence "due to the factor between F and calF as in (2.33) (which we evaluate at s = s0)".  The per-node factor
  tightens the P1 bound by up to 15x at the highest node.  Both readings are therefore run as declared variants:
  uv_SRa_ffs0_tip and uv_SRd_ffs0_tip use the factor frozen at s0.  Decision rule (extends PLAN 2.2, FESR row): the
  packaging/reading pair is chosen by the C5 criteria (UV +x end 0.0811 +- 5 %, x_ref shrink ratio) -- source
  identification against Fig.8 only, never against phase shifts; all four variants are reported.
- 20:10 UV SR-a tip accepted: x=0.0760129 (digitised cyan tip 0.0811, -6.3 %), f11/f00=-0.0650.  Phases at the tip:
  P1 90-degree crossing 773 MeV but min eta(P1)=0.585 (inelastic dip at 0.8 GeV); S0 134 deg at 0.9 GeV, 181 deg at 1.196
  (paper 84-103 / 87-106).  FESR moments: S0 n=0 sits on the upper edge of the raw box (+84 %), S0 n=1 +1.8 %, P1 n=-1 +13 %,
  P1 n=0 +1.4 %.  Spectral mean energies sqrt(M1/M0)*140: S0 713 MeV (targets give 959), P1 777 MeV (targets give 822):
  the rho position and the S0 rise track the moment ratios exactly as REVISED_CLAIMS 4b predicted for the raw box -- now
  confirmed on the regularised model.  The raw-absolute packaging leaves the S0 n=0 and P1 n=-1 moments effectively free.
  Consequence: the FESR packaging, a genuine 2309 under-specification, decides C5-C7.  Variants queued per the pre-registered
  rule (decided by Fig.8 numbers only): SR-b (relative 10 % per moment, the authors' released-code magnitude) with node and
  frozen FF factors, in addition to SR-a/SR-d node and frozen.
- 20:25 UV SR-d tip accepted: x=0.0759477 (-6.4 %), rho 773 MeV, min eta(P1) 0.585 -- indistinguishable from SR-a: the
  per-wave L2 ball of radius 2e-3 constrains the small moments (S0 n=0, P1 n=-1) as little as the box does.
- 23:55 UV SR-a upper section at x_ref accepted: y=-0.0046498 (paper cyan upper -0.0046002, +1.1 %); upper shrink 2.77e-4
  (paper 2.32e-4, +19 %, inside the +-25 % band).  Phases at this point: no P1 90-degree crossing below 1.2 GeV,
  min eta(P1)=0.244 -- under the raw-box packaging the near-black-point representative has no rho, like every earlier
  raw-box result.  Lower section running; SR-b and frozen-factor variants queued.
- 00:50 UV SR-a lower section at x_ref accepted: y=-0.0048571 -> lower RISE 2.77e-4 (paper 3.7e-5): the lower boundary moves
  up by as much as the upper moves down (ratio 1.0; paper 6.3).  C5 under SR-a: FAIL on the +x end (-6.3 %) and on the
  asymmetry (C5_RESULT_SRa.json).  The lower-branch point has a rho crossing at 809 MeV with min eta 0.63.
- 01:20 Frozen-factor UV tips accepted: SR-a/frozen 0.0796228, SR-d/frozen 0.0795407 (paper 0.0811, -1.8 / -1.9 %, inside
  the 5 % band) -- but at the tip P1 has min eta 0.026 and no 90-degree crossing below 1.2 GeV.  Diagnostics of the
  SR-a/node solutions: all 14 form-factor caps are ACTIVE at all three boundary points (used/cap = 1.000), S0 Gram is
  two-pion saturated (share 0.95-1.00, Watson phase within 1 deg), P1 |F| peaks at 5.2-6.7 near 0.7-0.8 GeV (physical
  |F_pi(rho)| ~ 6) with Watson holding at the tip and lower point but violated by 73-180 deg at the upper x_ref point,
  where P1 turns absorptive (share 0.67, eta 0.46) instead of resonant.  Interpretation: the finite problem as written
  (PSD Gram + moment boxes + FF caps) does not enforce two-pion saturation or Watson; which branch a boundary point takes
  is not fixed by the paper's constraints.  Sections for the frozen variants queued.
- 01:30 PRE-REGISTERED prediction for the SR-b (relative 10 % per moment) variants, written while uv_SRb_tip is still
  solving: with all four moments pinned to +-10 %, the spectral mean energies are forced to 959 +- ~5 % (S0) and
  822 +- ~5 % (P1); therefore the P1 90-degree crossing at the tip should land in 780-870 MeV and S0 at 0.9 GeV should
  drop from ~134 deg toward 85-110 deg.  If instead the tip still shows S0 > 120 deg at 0.9 GeV or no rho in
  760-880 MeV, the moment packaging is NOT the decisive missing convention and the degenerate-face (saturation) issue is.
- 02:05 Full-space check (no SVD basis reduction, 3876 coordinates, Mreg=1e2, M50/L10, eps=0.002): tip 0.0827768 vs the
  reduced-basis 0.0827570 -> +0.02 % (plan tolerance 1 %; the sign is the expected one, full >= reduced).  The basis
  reduction is validated for the mainline.
- SR-a/frozen upper section at x_ref accepted: y=-0.0044257 -> upper shrink 5.3e-5 (paper 2.3e-4): with the weaker
  (frozen) form-factor bound the upper boundary barely moves; at this point P1 crosses 90 deg at 1.03 GeV with min eta 0.65.
- SR-a/frozen lower section accepted: y=-0.0048853 -> lower rise 2.49e-4, upper shrink 5.3e-5, ratio 0.21 (paper 6.3):
  C5 FAIL for the frozen reading too (C5_RESULT_SRa_ffs0.json); the +x end passes (-1.8 %) but the asymmetry is inverted.
  Both FF readings with the raw box move the LOWER boundary by ~2.5-2.8e-4.
- SR-b (relative 10 % per moment) tips, node and frozen: after ~95 iterations SDPB is stalled -- step sizes 1e-4,
  primal error fixed at 3.7, mu ~1e7, primal objective diverging to -1e38 while the dual objective sits at 0.0797 / 0.0786.
  This is the signature of a primal-infeasible (or weakly infeasible) problem: pinning the S0 n=0 moment to +-10 % is
  incompatible with the (3.75) form-factor caps, the Gram blocks and unitarity at M=50 (all SR-a solutions sit at +84 %
  on that moment with all FF caps active).  Consistent with the 2026-09-12 Clarabel study (eps_FF = 6e-5 infeasible
  under SR-b at M=20/25).  The runs are left to hit their 7200 s solver budget so the failure is recorded.
  Consequence: the loose raw box is the only feasible reading of eps_SR = 2e-3 with eps_FF = 6e-5 in this finite model,
  i.e. the paper's stated tolerances do not pin the moments; the 822 MeV rho and the 85-110 deg S0 of Fig.9/10 are then
  properties of the particular optimal amplitude the authors' solver returned on a degenerate optimal face.

## 2026-09-14 02:25Z  Fig.3 sweep stalled on a source-hash check; face diagnostic pre-registered

- fig3_dir00..07 all failed at restore(): `Source operator changed: precision.py`.  The SR-d commit (bad9f7a) changed one
  line of PrecisionRows.uv_data (the SR-d tolerance); the pure Fig.3 model has no UV sector, so its rows are unaffected,
  but the whole-file hash check on the six operator modules is deliberate and is respected rather than bypassed.  The
  stubs were moved to stale_fig3_dirs_precision_hash/ (each holds only a status=assembling report).  The runner had been
  launching one stale direction every 90 s since 00:52Z whenever the live count dropped below 5; runner stopped 02:25Z.
  Fix: the pure tip is re-solved with the current source as fig3_pure_linf_1e3_tip_v2 (launched 02:27Z, same spec
  M50 L10 linf Mreg=1e3); the 24-direction sweep is re-queued from it once accepted.  All four UV tip sources and the
  x_ref sections were built after bad9f7a and pass the hash check (verified 02:25Z).

- Degenerate-face diagnostic, pre-registered before any result (user's choice A, 02:20Z).
  Construction: from an accepted support leaf with direction d, section fix_f00 and objective value v*, keep the section,
  add the 1x1 slab block d.(f00,f11) >= v* - 2e-6 and replace the objective by +-1 node functional.  Margin 2e-6 is
  twice the SDPB duality-gap threshold (1e-6 absolute for objectives below 1), 1 % of the Fig.8 upper shrink 2.3e-4 and
  2.6e-5 relative on the tip x; an exact equality at the optimum would leave no interior for the barrier method.
  Functionals (all exactly linear in the PMP variables): ImKH = Im(kappa h) = 1 - Re S at one node, ImS = Im S,
  ImF_ell(s_k), rho_hat_ell(s_k).  Nodes 37/38/39 sit at 0.732/0.792/0.864 GeV and bracket the paper's 822 MeV rho.
  Every leaf must pass the unchanged acceptance rules (solver optimal, Arb primal feasibility including the section and
  the slab); the functional value is re-evaluated in Arb from S at the node, not read from the solver objective.
  Reference (paper's Fig.9 P1: Breit-Wigner M=822 MeV, Gamma=150 MeV, eta=1): node 38 delta=68 deg, 1-ReS=1.72, ImS=0.69;
  node 39 delta=119 deg, 1-ReS=1.52, ImS=-0.85; node 37 delta=40 deg, 1-ReS=0.82.
  Rules:
  F1 (exclusion, exact up to the margin): max over the face of ImKH_P1(38) < 1.5  <=>  Re S_P1(0.792 GeV) > -0.5 on the
     whole near-optimal face; then no elastic resonance with delta >= 60 deg at 0.79 GeV lies on it and the paper's P1
     shape is NOT a property of this finite problem's optimal face at that representative point.
  F2 (degeneracy width): W = max - min of ImKH_P1(38) on the face.  W >= 0.5: the observable is undetermined by the
     extremal problem (representative curves are solver-dependent); W <= 0.05: determined (any Fig.9 discrepancy comes
     from the constraint set, not from the choice of representative); between: partially determined.
  F3 (post-resonance node): max ImKH_P1(39) >= 1.4 is necessary for a Fig.9-like curve; below it the paper's shape is
     excluded at that node as well.
  F4 (S0): the range of ImKH_S0(38) is reported as a delta range at eta=1 (delta = arccos(1-v)/2); C7 itself keeps its
     own rule.  F5 (form factor): the range of ImF_1(38) on the x_ref upper face is reported without a threshold.
  Jobs (12 solves, SR-a raw box, eps_FF 6e-5): node-factor tip (uv_SRa_tip/tip): ImKH P1 38 max/min, ImKH P1 39 max,
  ImKH S0 38 max/min; node-factor x_ref upper (uv_SRa_section_hi): ImKH P1 38 max/min, ImKH P1 39 max, ImF P1 38
  max/min; frozen-factor tip (uv_SRa_ffs0_tip/tip): ImKH P1 38 max, ImKH P1 39 max.

## 2026-09-14 02:55Z  SR-b: what "stalled" means, and the pre-registered moment-range test that replaces the inference

- The two SR-b tips (uv_SRb_tip, uv_SRb_ffs0_tip; launched 00:43/00:45Z, SDPB started 00:50Z) are the ORIGINAL runs, never
  restarted.  SDPB log at iteration 128 (6717 s): mu 1.8e50, primal objective -6.5e54, dual objective +0.0802 (node) /
  +0.0790 (frozen), gap 1.00, P-err 3.8e-4, p-err 1.4e3, D-err 3.0, steps ~1e-4.  In SDPB's pairing our variables are its
  dual; a primal objective running to -infinity with a bounded dual objective and vanishing steps is the interior-point
  signature of an infeasible (or weakly infeasible) dual, i.e. of an infeasible SR-b problem.  This is an inference from
  the iterate trajectory, NOT a certificate: SDPB 3.1.0 emits no infeasibility certificate.
- Note for readers of these directories: out/out.txt written 81 s after start ("maxIterations exceeded", runtime 81,
  matching the iteration-3 values) is the residue of SDPB's internal timing run ("Start timing run" in sdpb.log; the
  2-iteration iterations.0.json and c_minus_By.0 are the same residue, present in every leaf).  It is overwritten when the
  real solve terminates and must not be read as the result of the run.
- Accepted SR-a leaves (Arb FESR rows): S0 n=0 moment sits at target + 2e-3 exactly (slack 1e-10 at both tips, +83.9 %
  of the target 2.3851e-3); S0 n=1 at +1.8 % and P1 n=0 at +1.4 %, both at the raw box's upper edge (slack < 1e-8);
  P1 n=-1 free inside the box (-33 % .. +32 %).  So three of four raw boxes are active at their UPPER edge: every optimal
  amplitude wants larger moments than QCD's.  Under SR-b the S0 n=0 box is 8.4x tighter (2.39e-4) and P1 n=-1 4.7x
  tighter (4.23e-4) while S0 n=1 and P1 n=0 become 5.6x/7.3x looser; the two sets are not nested.
- Pre-registered moment-range test (rule SR-1), to replace the inference by a solver-precision proof: solve the SR-b
  model with the S0 n=0 box removed (spec.sr_free = ((S0,0),)) and the objective replaced by the S0 n=0 moment itself
  (functional SRmom S0 0 min, then max if needed).  The feasible set is convex and the moment is linear, so the
  achievable values form an interval [m_lo, m_hi].  SR-b is feasible iff [m_lo, m_hi] meets [0.9 t, 1.1 t],
  t = 2.3851e-3.  If m_lo > 1.1 t the min solve alone proves SR-b infeasible and (m_lo - t)/t is the smallest relative
  tolerance the S0 n=0 moment admits with the other three SR-b boxes; if m_lo <= 1.1 t the max solve is run.  Node-factor
  reading first (uv_SRb_tip spec + sr_free), frozen reading only if the node result is not decisive.  Same acceptance
  rules; the moment is re-evaluated in Arb from rho_hat (verification.functional), not read from the solver objective.

## 2026-09-14 03:35Z  SR-b drivers crashed at verification (operator error), post-processing re-run

- Both SR-b tips ended by SDPB's maxRuntime (7148 s node / 7134 s frozen; out.txt and y.txt complete).  The Python
  drivers then crashed with `NameError: name 'sr_free' is not defined`: arbaudit.py was edited in two steps (call site
  02:49Z, signature 02:51Z) while the drivers, which import arbaudit lazily at verification time, hit that window.
  Cause is the operator editing src/ under running drivers; nothing about the solve.  Rule from now on: no edits to
  src/smatrix_bootstrap/sdp while any driver is alive, except in one atomic write with tests passing first.
- Repair: scripts/sdp/finish_leaf.py re-runs exactly the readback/Arb-verification/convergence tail of run_once on the
  untouched SDPB outputs (Pmp.from_saved restores the saved basis; nothing is re-solved) and records the incident in the
  report under `postprocess`.  Expected outcome: status not_accepted, terminate reason maxRuntime exceeded, Arb
  verification listing which constraints the terminal iterate violates.
- Running drivers checked for the same hazard: the SR-d lower section (started 02:14Z) and the Fig.3 v2 tip (02:27Z)
  hold the pre-face precision_pmp in memory and call the audit without sr_free (default ()); the two face jobs (02:52/
  02:53Z) hold the face-era precision_pmp; all four resolve consistently against the current arbaudit.py.

## 2026-09-14 03:48Z  C5 with the SR-d packaging (per-wave L2 ball 2e-3, frozen FF factor): FAIL, identical to SR-a frozen

- uv_SRd_ffs0_section_hi y = -0.0044260268 (SR-a frozen: -0.0044257449), uv_SRd_ffs0_section_lo y = -0.0048844421
  (SR-a frozen: -0.0048852868), tip x = 0.0795407 (SR-a frozen 0.0796228).  C5_RESULT_SRd_ffs0.json: +x end -1.95 %
  (ok), upper shrink 5.35e-5 (paper 2.32e-4), lower rise 2.50e-4 (paper 3.7e-5), ratio 0.21 (paper 6.3): FAIL.
- Reading: replacing the per-moment raw box by the authors' per-wave L2 ball with the same 2e-3 changes the x_ref
  sections by < 1e-6 and the tip by 8e-5.  The FESR packaging is NOT the missing ingredient of Fig.8; the asymmetry
  (upper shrinks 6x more than the lower rises) is absent under every feasible reading of eps_SR = 2e-3 tried (SR-a
  node, SR-a frozen, SR-d frozen), and the relative reading (SR-b) is not feasible (03:35Z entry; proof pending SR-1).
  Node-factor SR-d tip (uv_SRd_tip) x = 0.0759477 vs SR-a node 0.0760129: same -6.3 % gap to the paper's 0.0811.

## 2026-09-14 04:40Z  First face result (tip, max ImKH_P1 at 0.792 GeV); eta-free reading and rule F6 pre-registered

- face_SRa_tip_ImKH_P1_38_max accepted (6079 s, gap 7.3e-7): max over the near-optimal face = 1.49349; the slab is
  active (slack 3e-14, f00 = x* - 2e-6), so the value is the ceiling permitted by the 2e-6 margin, not an interior
  optimum.  The tip's own amplitude has 1.470 at that node (delta 108 deg, eta 0.585): the maximiser moves it by +0.023
  only.  Pre-registered F1 (threshold 1.5, elastic Breit-Wigner reference): 1.4935 < 1.5 -> the elastic Fig.9 shape is
  excluded from the tip face, but by 0.4 % of the threshold, so F1 alone is not a comfortable verdict.
- Refinement, fixed before any further leaf is read: the paper never shows eta, so the reference must not assume eta=1.
  From the digitised Fig.9 curves (references/figure9_p1_phases.csv, linear interpolation to the node energies):
     node 37 (0.732 GeV): red 32.5, pink 31.3, light pink 32.8 deg
     node 38 (0.792 GeV): red 68.7, pink 69.9, light pink 76.2 deg      (90-deg crossings 0.827 / 0.824 / 0.813 GeV)
     node 39 (0.864 GeV): red 112.8, pink 114.5, light pink 124.1 deg
  With S = eta e^{2i delta}: 1 - Re S = 1 - eta cos(2 delta) and Im S = eta sin(2 delta).  At node 38 the paper's
  2 delta lies in [137, 152] deg, at node 39 in [226, 248] deg: cos(2 delta) < 0 at both nodes for all three curves, so
  ANY amplitude with the paper's phase there has 1 - Re S > 1 whatever eta > 0 is, and Im S > 0 at node 38.
  Rules (eta-free, necessary conditions; a violation excludes the paper's phase at that node for every eta):
  F1' : max_face (1 - Re S_P1)(node 38) < 1  =>  paper's delta_P1(0.792) excluded on that face.
  F3' : max_face (1 - Re S_P1)(node 39) < 1  =>  paper's delta_P1(0.864) excluded on that face.
  F6  : max_face Im S_P1(node 38) < 0  =>  the paper's pre-resonance sign at 0.792 GeV (crossing at 0.81-0.83 GeV)
        is excluded on that face for every eta > 0; runs face_*_ImS_P1_38_max queued (lines 3-4) for tip and x_ref upper.
  The original F1/F3 (elastic reference, thresholds 1.5 / 1.4) stay recorded and are reported alongside.
- x_ref upper face, max ImKH_P1(38): converging at 0.936 (gap 6.8e-6 at 04:35Z) -> F1' will exclude the paper's
  delta_P1(0.792) at the x_ref upper point for every eta (final value recorded when the leaf lands).

## 2026-09-14 04:45Z  x_ref upper face result; UV-stage Mreg control queued

- face_SRa_xrefhi_ImKH_P1_38_max accepted: max over the near-optimal face of 1 - Re S_P1(0.792 GeV) = 0.93559 (slab
  active, slack 7e-16).  F1': Re S_P1 >= 0.064 on the whole face, so |delta_P1(0.792)| < 43 deg for every eta there; the
  paper's 68.7-76.2 deg at that energy is excluded at the x_ref upper representative for every eta.  (Elastic F1 also
  fails.)  The source amplitude itself has 0.563 at that node (eta 0.71, delta -26 deg).
- Control pre-registered: the UV stage inherited Mreg = 1e2 from the chiral gate (rule (2)); its sensitivity was checked
  for Fig.3/Fig.4 but not with the UV blocks.  uv_SRa_mreg1e3_tip (SR-a node reading, Mreg = 1e3, otherwise the
  mainline) is queued.  Reading rule: if its +x end, P1 90-degree crossing and min eta differ from uv_SRa_tip by more
  than the Fig.4 Mreg spread (2 %) / 20 MeV / 0.1, the UV conclusions above are re-examined with an Mreg scan before
  any C5-C7 verdict is finalised; otherwise Mreg is recorded as not the missing ingredient.

## 2026-09-14 04:55Z  What the authors' released code (2403) and the Córdoba code (2511) actually do (read locally, not executed)

- GTB_numerics.m (arXiv:2403.10772 ancillary; local copy results/evidence/PV_PRIMARY_SOURCE_AUDIT_20260911/sources):
  nu0 = -20, s0 = (2 GeV)^2, three currents S0/P1/D0 with three moments each, sum rules as per-current L2 balls
  eS0 = 1e-7, eP1 = 6e-6, eD0 = 5e-6 on the normalised moments (i.e. 10-25 % of the QCD values), chiral L2 2e-3 (our
  chi-b), form-factor caps with "relaxation factors" (|F_S0| <= 0.05 above s0; |s F_P1| <= 2 x 6.87), and the
  regulariser  norm(rho/Mrho, 4) <= 1e2  -- an l4 M-regularisation, so the regulariser we had to reintroduce is the
  authors' own device.  Solver: MOSEK (CVX).  Selection: block 1 fixes F0 = ni/50 x 0.163108 (ni = 44/46/48) and
  maximises F1; block 2 ("Watsonian unitarization") re-imposes every constraint WITHOUT fixing F0 or F1 and maximises a
  linearised saturation functional v1+v2 built from the block-1 form factors, "repeated until convergence".  The
  reported amplitude is therefore the most Watson-saturated point of the whole feasible set, not the amplitude of the
  boundary point; the boundary point only seeds the linearisation.
- Córdoba et al. (2511.11513, MATLAB): s0 = (2 GeV)^2, Mreg = 10 with norm(x, 4) <= Mreg x 10^power, per-current L2
  sum-rule tolerances MSVZ0 = 1e-7 / MSVZ = 4e-6 / MSVZ2 = 6e-6, squared FF caps 3e-8 / 2e-6 / 4e-2, MOSEK; they
  report the rho(770) and f2(1270) of 2403 and a sigma pole.  Same family of choices as 2403, none of them in 2309.
- Consequence for the 2309 reproduction: every later implementation (i) regularises the double density (l4), (ii) uses
  s0 = 2 GeV, a third current and 3 moments per current with ~10 % relative tolerances, and (iii) selects the reported
  amplitude by a saturation functional over the full set.  None of (ii)-(iii) is stated in 2309; our face diagnostic
  shows that with the 2309-stated set the boundary-point amplitudes cannot carry the paper's P1 phase at 0.79 GeV
  (x_ref upper: excluded for every eta; tip: face rigid at 1.47-1.49).  The "10 % precision" the collaborator mentions
  matches the 2403/2511 tolerances on normalised moments, and in the 2309 model that reading (SR-b) is not feasible
  (SR-1 solve running).
- Our tip amplitude: P1 crosses 90 deg at 773 MeV with eta = 0.585 at 0.79 GeV (chiral-only stage: no crossing below
  1.2 GeV, C4).  Our x_ref upper amplitude: |S_P1| dips to 0.24 at 0.73 GeV; the node values admit two phase branches,
  the nearest-node lift (delta swinging to -47 deg) and the +180 deg branch (delta rising through 90 deg at 0.708 GeV);
  the finite problem does not decide between them.  Both readings are recorded; neither is the paper's 0.82 GeV.

## 2026-09-14 05:05Z  Tip face: min ImKH_P1(38) = 1.44112  ->  range [1.441, 1.493], width 0.052

- face_SRa_tip_ImKH_P1_38_min accepted.  Together with the max (1.49349): on the tip's near-optimal face (margin 2e-6)
  1 - Re S_P1(0.792 GeV) lies in [1.441, 1.493]; the tip's own amplitude sits at 1.470, mid-range.  F2: width 0.052,
  just above the 0.05 "determined" threshold -> "partially determined" by the letter of the rule, in substance rigid:
  the extremal problem fixes the P1 S-matrix at 0.79 GeV to +-0.026 at the tip.  The paper's 68.7-76.2 deg would need
  > 1 for any eta and 1.70-1.85 at eta = 1; neither is inside the range, so no representative choice on the tip face
  reaches the Fig.9 tip curve.

## 2026-09-14 05:21Z  x_ref upper face: min ImKH_P1(38) = 0.02279  ->  range [0.023, 0.936], width 0.91

- F2 at the x_ref upper representative: undetermined (width 0.91 >= 0.5).  The extremal problem leaves Re S_P1(0.792 GeV)
  anywhere in [0.064, 0.977] on the near-optimal face: the paper's premise that a boundary point carries one set of
  partial waves fails there, while the whole face still has Re S > 0, i.e. the paper's 68.7-76.2 deg remains excluded
  for every eta (F1').  Contrast with the tip, where the same observable is pinned to [1.441, 1.493].

## 2026-09-14 05:55Z  SR-1 verdict: the relative (10 %) reading of the sum rules is infeasible in the 2309 model -- proven; tip node 39

- srmom_SRb_free_S0n0_min accepted (gap 9.3e-7; Arb: unitarity, chiral, regulariser, Gram, FF and the three imposed
  SR-b boxes all pass).  Minimum of the S0 n=0 moment over the SR-b feasible set with its own box removed:
  m_lo = 3.74537e-3 = 1.570 x target (t = 2.38510e-3).  Since the feasible set is convex and the moment linear, every
  point of that set has moment >= m_lo > 1.1 t = 2.6236e-3, so the SR-b box |m - t| <= 0.1 t is empty: SR-b is
  infeasible, to solver precision, independently of the 03:35Z stall inference.  At the minimiser S0 n=1 sits at its
  lower edge (-10.0 %) and P1 n=-1 at its upper edge (+10.0 %): lowering the S0 n=0 moment is paid for by the other two
  low moments, and the smallest relative tolerance the S0 n=0 sum rule admits in this model is +57 %.
  The minimiser is an interior amplitude, (f00, f11) = (0.06304, -0.00396); P1 crosses 90 deg at 756 MeV, S0 at 619 MeV.
- Reading: the "10 %" tolerances of the later codes (2403: eS0 = 1e-7 on the normalised S0 moments, i.e. ~23 %; Córdoba
  MSVZ0 = 1e-7) live in a model with s0 = 2 GeV, three moments per current and a third current; transplanted into the
  2309 model at s0 = 1.2 GeV they cannot be met.  The literal 2e-3 raw box is the only feasible reading found, and under
  it the S0 n=0 moment is pinned to its UPPER edge (+83.9 %) at every UV representative.  Both facts say the same thing:
  the 2309-stated finite problem wants a low-energy S0 spectral integral 1.6-1.8 times the QCD value.
- face_SRa_tip_ImKH_P1_39_max accepted: max over the tip face of 1 - Re S_P1(0.864 GeV) = 0.63457 < 1.  F3': the
  paper's 112.8-124.1 deg at that energy (cos 2 delta < 0 for all three curves) is excluded on the tip face for every
  eta; elastic F3 (threshold 1.4) fails as well.  With F1' (x_ref upper, node 38) and the tip range at node 38, the
  paper's Fig.9 P1 curve is now excluded on the near-optimal faces of both representatives examined, at one node each
  without any assumption on eta.

## 2026-09-14 06:03Z  x_ref upper face, node 39: max ImKH_P1(0.864 GeV) = 0.70098 < 1

- F3' at the x_ref upper representative: the paper's 112.8-124.1 deg at 0.864 GeV is excluded for every eta (Re S_P1 >=
  0.299 on the whole face).  Both representatives examined are now excluded at both nodes 38 and 39 without any
  assumption on eta; the remaining P1 runs (Im S sign at node 38, frozen-factor tip) can only add to this, not reverse it.

## 2026-09-14 06:30Z  Code-versus-paper conformance audit of the UV sector (re-derived from the 2309 source, not from our own reading)

Checked against prd_submission_2.tex equations (h22)-(h25), (srexpression), (srnumbers), (h39)-(h43), (FFasym):
- Kinematic factors: k0^2 = 3/(16 pi (2pi)^4) sqrt((s-4)/s) and k1^2 = (s-4)^{3/2}/(24 pi (2pi)^4 sqrt s) from (h24)/(h25);
  formfactor.kinematic_factor squares agree to 6e-16 at s = 10, 30.25, 73.47, 200.  Gram block (positiveB) is imposed
  with the exact congruence rho_hat = rho/k^2; S in the block is the same S as in the unitarity disks.
- FESR targets: printed (srnumbers) x s0^{n+2} with s0 = (1200/140)^2 = 3600/49.  Recomputing (srexpression) with the
  inputs (qcddata1-2) reproduces the printed S0 numbers to -0.7 % with m_q = rms(m_u, m_d) and the P1 numbers to +0.1 %;
  the printed numbers are what the paper imposed, so they are used (recomputed_targets kept as a control).
- Discretisation (srnum): (pi/M) sum (ds/dphi)_i s_i^n rho_i over the nodes with s_i <= s0 (43 nodes; last node s =
  73.40 vs s0 = 73.47, first excluded s = 97.3).  quad_weights = (pi/M) ds/dphi and PrecisionRows.w = (1/M) ds/dphi with
  the explicit pi in uv_data and in the Arb audit: the float and Arb rows are the same weight (ratio pi checked).
  The authors' 2403 code uses the same hard cutoff (is0 = indexes(1)-1, with the alternative commented out).
  Caveat recorded: for a smooth integrand this midpoint sum overshoots int_4^{s0} x^n dx by +14.6 % (n=0) and +28.9 %
  (n=1), with the last node carrying 24 % / 41 % of the weight; in the accepted solutions the last node contributes
  only 1-6 % (tip) and up to 16 % (x_ref, S0 n=1) / 31 % (frozen tip, P1 n=0) of the moments, so the artefact is
  present but not dominant.  A cell-clipped quadrature would be a deviation from both the paper and the 2403 code; it is
  noted as an optional sensitivity control, not part of the mainline.
- Form-factor caps (FFasym): |cF_0(s_i)|^2 <= 2 m_q^2 eps_FF and |cF_1(s_i)|^2 <= eps_FF/2 on the 7 nodes above s0,
  cF = k F with F normalised by Re F = 1 + K Im F (F(0) = 1; for j_S this is the LO sigma-term normalisation m_pi^2 = 1).
  With eps_FF = 6e-5 this means |F_1(s0)| <= 0.228 and |F_0(s0)| <= 0.072 (node reading).  The frozen reading keeps
  k(s0) for all seven nodes.  Both were run.
- Spectral density above s0: no constraint in the paper besides the Gram block, and none in our model; the solver
  returns rho_hat ~ 1e13-1e15 there (a free direction).  Harmless for every imposed constraint, recorded for the reader.
- Chiral (chiralconstraints): 4 points s = 1/2, 1, 3/2, 2, ratios R01, R21, one combined L2 norm 2e-3 (the 2403 code
  uses exactly this packaging, echi = 2e-3).
Conclusion: every quantity the paper defines is implemented as written; the open items are exactly the ones the text
leaves open (eps_SR norm, FF factor evaluation point, regulariser, representative selection, solver) and each has been
run in at least two readings with the controls recorded above.

## 2026-09-14 06:31Z  Frozen-factor tip face, node 38: max ImKH_P1(0.792 GeV) = 0.32497 < 1

- Under the frozen (3.75) reading the UV tip (uv_SRa_ffs0_tip, x = 0.07962) has no P1 crossing below 1.2 GeV (first at
  1.49 GeV), and its near-optimal face cannot supply one at 0.79 GeV either: Re S_P1 >= 0.675 on the whole face, the
  paper's 68.7-76.2 deg excluded for every eta (F1').  So the rho at the tip is a property of the per-node (3.75)
  reading (mainline), not of the frozen one; C6 keeps the mainline reading, the frozen result is reported as the
  declared alternative.

## 2026-09-14 06:37Z  Tip face, F6: max Im S_P1(0.792 GeV) = -0.27605 < 0

- On the whole tip near-optimal face Im S_P1 at 0.792 GeV is negative (<= -0.276), i.e. the P1 wave is already past its
  90-degree crossing at that energy for every amplitude on the face; the paper's curves cross at 0.813-0.827 GeV and have
  Im S > 0 there.  F6 excludes the paper's pre-resonance sign at the tip for every eta > 0, independently of F1/F3'.
  The tip's own value is -0.349; the range [.., -0.276] is again narrow (rigid face).

## 2026-09-14 06:45Z  F6 at the tip, frozen-factor tip node 38, and what the authors' CURRENT (2505-era) code does

- face_SRa_tip_ImS_P1_38_max accepted: max over the tip face of Im S_P1(0.792 GeV) = -0.27605 < 0.  F6: every amplitude
  on the tip's near-optimal face is already past the resonance at 0.79 GeV (Im S < 0), for every eta > 0; the paper's
  curves (90-deg crossing at 0.813-0.827 GeV, Im S > 0 at 0.79) are excluded there independently of eta.  Combined with
  the node-38 range [1.441, 1.493] and node-39 max 0.635, the tip face is characterised: a P1 resonance that has
  crossed 90 deg below 0.79 GeV (crossing 0.772-0.774 GeV on every extremal leaf), eta 0.56-0.59 at the crossing.
- face_SRa_ffs0_tip_ImKH_P1_38_max accepted: max 1 - Re S_P1(0.792) = 0.32497 < 1 on the frozen-factor tip face
  (source tip x = 0.07962, its own P1 crossing at 1.49 GeV): the paper's phase at 0.79 GeV is excluded for every eta
  under the frozen reading as well; that reading has no rho below 1.2 GeV at the tip at all.
- x_ref upper Im S max converging at +0.124 (gap 4.7e-5): F6 will NOT exclude at x_ref (Im S may be positive there), but
  F1' already does (Re S >= 0.064): with Im S <= 0.124 the phase on that face is bounded by ~31 deg where Im S > 0.
- Authors' current code (upstream main, gtb_qcd_02_optimize.m / optimize_core.m, read only): nu0 = -20, s0 = 2 GeV,
  three currents with nSR = 6 moments each, spectral densities as Na = 49 coefficients with explicit bounds, sum rules
  as per-moment boxes abs(w) <= eS0 = 0.05 / eP1 = 0.1 / eD0 = 0.1 on residuals that are evidently normalised (an
  absolute 0.05 on the raw normalised moments ~1e-7..1e-5 would be vacuous), FF asymptotics with multipliers 8 / 2.2 / 8,
  chi L2 2e-3, the same l4 regulariser norm(rho/Mrho,4) <= 1e2, MOSEK.  Objective: ni = 0 is a pure feasibility solve
  (minimize 0); ni = 1..5 maximise the linearised Watson functional.  There is NO (f00, f11) support functional and no
  boundary point anywhere in the current method: the reported amplitude is the Watson fixed point of the feasible set.
- Reading: the 5-10 % relative sum-rule boxes that the collaborator calls "10 % precision" are the 2505 setting at
  s0 = 2 GeV; at s0 = 1.2 GeV in the 2309 model the same reading is provably infeasible (SR-1, 05:55Z).  The evolution
  2309 -> 2403 -> 2505 (larger s0, more moments, Watson selection, no boundary point) is consistent with the 2309 fine
  features not being properties of the 2309-stated finite problem.

## 2026-09-14 06:45Z  x_ref upper face, F6: max Im S_P1(0.792 GeV) = +0.12405 > 0  (sign test inconclusive there)

- F6 does not by itself exclude the paper's sign at the x_ref upper point: Im S can reach +0.124 on the face.  Combined
  with F1' (Re S >= 0.064 on the same face) the admissible S-matrix points at 0.792 GeV have 2 delta <= 63 deg, i.e.
  delta <= 31 deg for every eta -- still far from the paper's 68.7-76.2 deg, so the exclusion at x_ref rests on F1'
  (and F3' at node 39), not on the sign.  Recorded as such; no rule is changed after the fact.

## 2026-09-14 07:07Z  x_ref upper face: max Im S_P1(0.792 GeV) = +0.12405

- F6 does not exclude at the x_ref upper representative (Im S may be positive there, up to 0.124); F1' already excludes
  the paper's phase there (Re S >= 0.064 on the face).  Combined: on that face S_P1(0.792) lies in the sector
  Re S in [0.064, 0.977], Im S <= 0.124, i.e. |delta| <= 31 deg wherever Im S > 0 -- far from the paper's 69-76 deg.
- Frozen-factor tip node 39 converging at 0.645 (< 1): F3' excludes there as well (recorded when the leaf lands).
- Mreg = 1e3 control at iteration 60: dual objective 0.0661 rising (tip x at Mreg = 1e2: 0.07601); verdict when accepted.

## 2026-09-14 07:16Z  Frozen-factor tip face, node 39: max ImKH_P1(0.864 GeV) = 0.64519 < 1  (F3' excluded for every eta)

## 2026-09-14 11:50Z  (operator-18) C8 ladder: M=60 hit the 7200 s budget one step short; resumed from its checkpoint

- uv_M60_L10_tip ended by maxRuntime at iteration ~67 (97 s per iteration at M=60), duality gap 4.9e-3, terminal
  f00 = 0.075033 (not accepted, not used).  Resumed with `sdp refine --resume` from the saved sdp.ck under the same
  PMP, precision and rank layout (uv_M60_L10_tip_resume, budget 7200 s); acceptance rules unchanged.
- Accepted so far (SR-a node reading, Mreg = 1e2): (50,8) x = 0.076176, P1 crossing 771 MeV, min eta 0.584;
  (50,10) 0.076013, 773 MeV, 0.585; (50,12) 0.075944, 773 MeV, 0.585; (45,10) 0.074275, 736 MeV, 0.226.
  L spread at M=50: 2 MeV (rule <= 20).  S0 at 1 GeV is ~150 deg at every configuration (rule 85-105): C8 will fail on
  that criterion whatever M=60 gives; the L-stability part of the paper's statement is reproduced.

## 2026-09-14 12:19Z  (operator-18) C8 verdict (Fig.11, M/L stability of the UV tip): FAIL by the letter; L-stability reproduced

- uv_M60_L10_tip_resume accepted from the checkpoint (checkpoint_loaded true, 13 further iterations): x = 0.075097,
  P1 crossing 739 MeV, min eta 0.597, S0 crossing 667 MeV.  Full ladder (SR-a node reading, Mreg = 1e2):
     (50,8)  x 0.076176  P1 771 MeV  eta 0.584  S0(1 GeV) 151 deg
     (50,10) x 0.076013  P1 773 MeV  eta 0.585  S0(1 GeV) 150 deg
     (50,12) x 0.075944  P1 773 MeV  eta 0.585  S0(1 GeV) 150 deg
     (45,10) x 0.074275  P1 736 MeV  eta 0.226  S0(1 GeV) 165 deg
     (60,10) x 0.075097  P1 739 MeV  eta 0.597  S0(1 GeV) 156 deg
- Pre-registered rule (claims.c8): L spread 1.5 MeV (<= 20: pass); M spread 36.6 MeV (40-70: fail by 3 MeV); ordering
  M60 < M45 < M50 required, observed M45 (736) < M60 (739) < M50 (773): fail; S0 at 1 GeV 150-165 deg (85-105: fail).
  Verdict FAIL.  Reading: the paper's qualitative statement -- "for M=50 ... L=8,10,12 ... reasonably good convergence"
  and "the P1 phase shifts show a rho resonance with the peak slightly shifted depending on M" -- is reproduced (2 MeV
  in L, 37 MeV in M); what is not reproduced is the absolute S0 level at 1 GeV (ours 150 deg vs the paper's ~95 deg),
  the same UV-stage S0 discrepancy as in C7, and the paper's M-ordering.

## 2026-09-14 12:30Z  (operator-18) Consolidated record of results not yet logged by the peer session (07:30Z-09:42Z)

- Mreg control (rule 04:45Z): uv_SRa_mreg1e3_tip accepted: x = 0.076713 (+0.92 % vs Mreg = 1e2; rule <= 2 %), P1
  crossing 769 MeV (-4 MeV; rule <= 20), min eta 0.560 (-0.025; rule <= 0.1), S0 crossing 696 MeV, f11 -0.004982.
  All three within threshold: the regulariser scale is NOT the missing ingredient; the Fig.8 asymmetry and the 822 MeV
  rho do not appear at Mreg = 1e3 either.
- Tip face, S0 node 38: 1 - Re S_S0(0.792 GeV) in [1.58975, 1.72638] (width 0.137); the phase is pinned at 111.0-111.4
  deg on all three leaves and only eta varies (0.79-0.99).  Paper Fig.10 at 0.792 GeV: red 76.1, pink 84.8,
  light pink 98.6 deg; at 1.0 GeV 91-106 deg vs our 150 deg.  F4 reported; the S0 discrepancy is UV-induced (C4 passed).
- x_ref upper face, Im F_1 node 38: [0.07843, 3.70991] (F5): the vector form factor at 0.79 GeV is essentially
  undetermined by the extremal problem at that point.  FACE_RESULT.json is the final face summary (14 leaves).
- UV upper sections (SR-a node, Mreg 1e2), x_ref + {-0.002, -0.001, 0, +0.001, +0.002}: y = -0.004485, -0.004566,
  -0.004650, -0.004738, -0.004836.  Frozen rule (PLAN 2.2) applied: representative = x_ref upper (distance to P
  0.00024), mid = x_ref + 0.001 upper, tip = +x end; receipt UV_REPRESENTATIVE_SELECTION.json.
- C6 (Fig.9), C67_RESULT.json: FAIL -- tip crossing 773 MeV (min eta 0.585), mid 739 MeV (min eta 0.077), ref no
  crossing on the nearest-node lift (|S| dips to 0.24 at 0.73 GeV; +180 deg branch crosses at 0.708 GeV); paper
  813-827 MeV with three nearly coincident curves.  A rho-like structure is present at all three points, 708-773 MeV,
  strongly inelastic; the pre-registered band 795-845 MeV and min eta >= 0.9 are not met.
- C7 (Fig.10), same file: FAIL on S0 only -- S2 RMS 0.35-1.44 deg and endpoints -28 to -32 deg (paper -18 to -36) pass;
  S0 RMS 5.7-16.6 deg with endpoints at 1.196 GeV of 181-192 deg (paper 87-106 deg).

## 2026-09-14 15:40Z  (operator-18) C1 verdict (Fig.3, 24 supports at Mreg = 1e3): FAIL on three of four extrema; +x end passes

- All 24 directions accepted (fig3_dir00-23, reusing the fig3_pure_linf_1e3_tip_v2 PMP; dir00 reproduces the tip
  bit for bit).  Extrema vs the digitised Fig.3 (rule: each within 2 %):
     f00_max  2.2245 vs 2.2329  (-0.38 %)  pass
     f00_min -2.2470 vs -2.9020 (+22.6 %)  fail (ours inside)
     f11_min -0.6243 vs -0.7336 (+14.9 %)  fail (ours inside)
     f11_max  0.0663 vs 0.0793  (-16.4 %)  fail (ours inside)
- Reading: our region is strictly inside the paper's, and the deficit sits on the -x and +-y extrema, where the
  amplitudes are large and the |rho_ij| <= Mreg cap is the binding constraint; the +x end, the only extremum that
  feeds the chiral stage, is Mreg-insensitive (2.2116 / 2.2245 / 2.2305 at Mreg = 1e2 / 1e3 / 1e4, MREG_STAR_PURE).
  The paper states no regulariser (a double-precision solver on the unregularised problem returns whatever its implicit
  regularisation allows), so the three failing extrema measure exactly the under-specified item.  Control queued:
  fig3_pure_linf_1e4_tip_v2 (fresh source, current hashes), then -x/+y/-y supports at Mreg = 1e4, to report how far
  the extrema move per decade of Mreg.  C1 stays FAIL under the pre-registered rule; the control is a sensitivity
  annotation, not a re-tuning.

## 2026-09-17 18:05Z  (operator-18) C1 Mreg = 1e4 control landed on 2026-09-14 (recorded now)

- fig3_pure_linf_1e4_tip_v2 accepted, +x end 2.23046 (paper 2.23289, -0.11 %).  Supports from it: -x end -2.6927
  (paper -2.9020; at Mreg = 1e3 it was -2.2470), +y 0.07562 (paper 0.07932; 1e3: 0.06629), -y -0.71764 (paper
  -0.73359; 1e3: -0.62430).  Per decade of Mreg the three failing extrema move from 15-23 % inside the paper's region
  to 2-7 % inside; the +x end moves by 0.3 %.  The paper's Fig.3 therefore corresponds to a larger effective
  regulariser (or none, i.e. whatever the double-precision solver tolerated).  C1 verdict unchanged (rule fixed at
  Mreg* = 1e3); this is the sensitivity annotation promised in the comparison document.

## 2026-09-17 18:40Z  (operator-18) Watson / unitarity-saturation iteration: pre-registration

- Trigger: the authors' reply of 2026-09-17 -- "|S| does not saturate unitarity. We have an iterative procedure that
  improves saturation of unitarity ... Our results are always shown after that."  The 2309 text does not describe this
  step; the released 2403 code does (block 2).  It was deliberately not applied in our main line (PLAN 2.1).
- Implementation (src/smatrix_bootstrap/sdp/watson.py, `support --functional watson`): every constraint of the finite
  model kept; the section f00(3)=x released (as in the authors' block 2; `--watson-keep-section` as a control); the
  objective replaced by  max sum_{waves S0,P1,S2} sum_{s_k <= s0} [Re(conj(t_k) h_k) - Im h_k]  with h = kappa f
  (S = 1 + i h) and targets from the previous accepted leaf: S0/P1 t_k = -i(F_k/conj(F_k) - 1) (Watson's theorem at
  |S| = 1, F the previous form factor), S2 t_k = radial projection of the previous h_k onto |S| = 1.  On the unitarity
  disc |h|^2 <= 2 Im h this functional is maximised exactly at h = t, so the fixed point is the saturated,
  Watson-aligned amplitude; the functional value at the fixed point is sum Im t_k.  Unit node weights; no weights,
  no D0 current, no s0 = 2 GeV: only the 2309 model with a different selection of the amplitude.
- Verification unchanged: SDPB optimal + Arb feasibility; the functional is re-evaluated in Arb from the node S.
- Runs: two chains, from the accepted SR-a node-reading leaves uv_SRa_tip/tip (tip) and uv_SRa_section_hi (x_ref
  upper), 5 rounds each, stop earlier if max_k |h_k - t_k| < 0.02 or a round is not accepted.  Start distances:
  see WATSON_ITER.json round 0.
- Readings, fixed now: per round, (f00, f11) drift, min |S_P1| below 1.2 GeV, P1 90-degree crossing, S0 phase at 0.792
  and 1.0 GeV, Watson residual and two-pion fraction.  Afterwards C6/C7 are re-evaluated on the final round of each
  chain as a separately labelled "post-Watson" reading (the original verdicts stand for the text-only model).  The
  face-diagnostic conclusions are unaffected (they concern the boundary points themselves).
- Equivalence checked against the authors' text (arXiv:2403.10772 Section 2.5, eq. 2.29, F_W = sum int ds
  Re[e^{-2 i alpha}|old (S-1)], alpha = old form-factor phase, or the old phase shift for waves without a current):
  with S - 1 = i h and t = 2 sin(alpha) e^{i alpha}, Re(conj(t) h) - Im h = Re[e^{-2 i alpha}(S - 1)] identically, so
  our per-node term IS (2.29) at the node; the released 2403 code likewise sums nodes up to s0 without ds weights,
  keeps every constraint and drops the F0 fix.  Their 2505 code differs only in the target magnitude (|h_old| instead
  of 2 sin alpha), which does not change the maximiser on the unitarity disc.  Declared differences: waves S0, P1, S2
  (the ones 2309 plots) instead of their six; no D0 current; s0 = 1.2 GeV; unit weights in the physical h.

## 2026-09-17 19:50Z  (operator-18) Watson chains: the keep-section control was vacuous at the tip; replaced by pinned-point controls

- watson_tip_keepsec stopped after 40 minutes: the tip leaf has no section (it is a +x support), so "keep section"
  changed nothing and the chain duplicated watson_tip iteration for iteration (identical SDPB trajectories).  Directory
  moved to stale_watson_tip_keepsec_vacuous/.
- The question it was meant to answer -- does the plotted amplitude stay at the boundary point? -- is answered by a
  pinned-point control instead: constraints and section from the start leaf plus the near-optimal slab
  d.(f00,f11) >= v* - 2e-6 (the face-diagnostic slab), Watson targets from the previous round.  Two such chains were
  launched at 19:50Z: watson_ref_pinned (x_ref upper, whose face is wide) and watson_tip_pinned (tip, whose face is
  rigid).  If saturation improves little under pinning while the released chains saturate, the authors' plotted
  amplitudes cannot be the boundary-point amplitudes.  Released chains watson_tip / watson_ref / watson_mid continue
  (round 1 at SDPB iteration ~45 of an expected ~85).

## 2026-09-17 20:35Z  (operator-18) Watson chain, tip, round 1 (released section)

- Accepted (gap < 1e-6, Arb feasible).  The amplitude left the tip: (f00, f11) = (0.076013, -0.004944) ->
  (0.073874, -0.004857), i.e. inward by 2.8 % in f00, to the x_ref neighbourhood.  S0 became elastic at all 43 nodes
  below s0 (Watson residual 0, two-pion fraction >= 0.996); S2 unchanged (elastic); P1 elastic everywhere except the
  0.792 GeV node, where |S| rose from 0.585 to 0.932 (|F|^2/rho_hat 0.80 there).  P1 crossing 773 -> 770 MeV,
  delta_P1(0.792) 108 -> 111 deg; S0 crossing 694 -> 688 MeV, S0(1 GeV) 150 -> 148 deg: the fast S0 rise is not a
  saturation effect.  FESR: S0 n=0 still at the upper edge (+84 %), P1 n=-1 moved from +13 % to +34 %.
- Artefact noted: after the round the P1 form factor at the threshold node (0.280 GeV) is |F| = 18 with phase 75 deg
  while |F|^2/rho_hat = 0 there -- a free direction the objective does not see.  The target built from that phase for
  the next round is noise, but the node has kappa ~ 0.15 and cannot move S; its only effect is on the raw convergence
  metric.  The driver now also records max |h - t| over the two-pion-saturated nodes (|F|^2/rho_hat >= 0.5) and uses
  that for the stop rule in future launches (tip: 0.418 -> 0.183); the running chains keep their 5-round budget.

## 2026-09-17 20:50Z  (operator-18) Watson chain, x_ref upper (ref), round 1 (released section)

- Accepted.  (f00, f11): (0.073321, -0.004650) -> (0.069295, -0.004497), inward by 5.5 % in f00.  S0 elastic at all
  nodes below s0 (was |S| = 0.52 at 1.2 GeV); S0 crossing 684 -> 686 MeV, S0(1 GeV) 131 deg.  P1: the |S| = 0.24 dip
  at 0.73 GeV became a clean, nearly elastic resonance (|S| = 0.883 at 0.73 GeV, 1.000 elsewhere) with the phase rising
  14 -> 49 -> 117 -> 157 deg over 0.60-0.79 GeV and the 90-degree crossing at 711 MeV (the +180-degree branch of the
  round-0 amplitude had given 708 MeV).  Distance to the fixed point over two-pion-saturated nodes 0.856 -> 0.377.
- Reading so far: one round of the authors' selection removes the inelasticity almost completely at both points, but
  the rho position it leaves (770 MeV at the tip chain, 711 MeV at the ref chain) is below the paper's 813-827 MeV, and
  the two chains do not (yet) agree with each other.  The authors' statement that the iteration converges to a solution
  independent of the starting functional is what the remaining rounds test.

## 2026-09-17 21:30Z  (operator-18) Pinned-point controls, round 1

- watson_ref_pinned (x_ref upper, section kept + slab 2e-6, Watson objective): accepted, point unchanged
  ((0.073321, -0.004652), slab active).  Saturation improves only partly: P1 min |S| 0.244 -> 0.367 at 0.73 GeV
  (0.625 at 0.68 GeV, 0.889 at 0.79 GeV), mean(1-|S|^2) over the nodes below s0 0.050 -> 0.039; S0 min |S| 0.52 -> 0.88.
  The released chain from the same start reached P1 min |S| 0.883 and S0 fully elastic in one round.
- watson_tip_pinned (tip, slab 2e-6): accepted, point unchanged; P1 min |S| 0.585 -> 0.593, essentially nothing --
  as the face diagnostic predicted (the tip face is rigid in the P1 observables).
- Reading: at the boundary points themselves the saturation the authors describe is not available; the iteration
  saturates only because it leaves the boundary (2.8 % inward at the tip, 5.5 % at x_ref after one round).  Whatever
  coordinates the paper's dots mark, the amplitudes it plots after the iteration cannot be the boundary-point
  amplitudes of this finite problem.

## 2026-09-17 23:05Z  (operator-18) Correspondence with the released code re-checked; first-round inputs verified; weighted control launched

- Variable mapping (audit of the 2403 notebook, PV_PRIMARY_SOURCE_AUDIT_20260911): h-tilde = h/Lambda_l, Im h-hat =
  Im h/Lambda_l^2, F-hat = F/Lambda, Lambda_l(s) = ((sqrt s - 2)/(sqrt s + 2))^(l/2); phase shifts read as arg(h-tilde)
  for the current waves, i.e. S = 1 + i h-tilde with h-tilde = h -- our convention.  Their block-2 term
  Re(conj(fnl) h-tilde) - Im h-hat with fnl = 2 Im F-hat F/|F|^2 equals (1/Lambda_l^2)[Re(conj(t) h) - Im h] with our
  t = 2 sin(alpha) e^{i alpha}: our functional with node weights 1/Lambda_l^2 (1 for S0, S2; (sqrt s+2)/(sqrt s-2) for
  P1, ~1800 at the threshold node, 2.2 at the rho).  Same per-node maximiser, hence the same fixed point when the
  saturated amplitude is feasible; a different compromise otherwise.  Not a transcription: our objective is built from
  our own physical-h rows and Arb-verified; theirs from their scaled kernels.
- Inputs of the first rounds verified on watson_tip/round_01, watson_ref/round_01, watson_ref_pinned/round_01,
  watson_tip/round_02: recorded targets identical to a recomputation from the linearisation source; released leaves
  carry the base blocks only (ref 9175 -> 9173, section removed; tip 9173), pinned leaf 9176 (section + slab, slab
  active); the objective row in pmp.json equals an independent float construction from the unscaled operator rows
  projected on the saved basis to 1e-11 (row scale 1e2); SDPB dual objective = Arb recomputation = float check.
- Control launched: watson_ref_authorsweights (released, weights 1/Lambda_l^2), to measure how much the node weighting
  moves the compromise at the rho node where full saturation is not feasible.

## 2026-09-18 00:50Z  (operator-18) watson_tip_pinned stopped at its fixed point

- Rounds 2 and 3 of the pinned tip chain are identical to 1e-6 in every recorded quantity (f00 0.0760109, f11
  -0.0049434, P1 min |S| 0.5926, W 25.4888): the pinned iteration has converged, without saturating, to the tip
  face's own amplitude.  Further rounds would repeat it; the chain was stopped and its round-4 directory (just
  started) removed.  Result on record: at the tip the authors' step cannot saturate unless the point is released.

## 2026-09-18 01:55Z  (operator-18) Sensitivity runs pre-registered (the authors' "depends on the parameters")

- Five UV tips (SR-a node reading, Mreg 1e2, M50 L10), everything else as the main line:
  eps_SR = 5e-3 and 1e-2 (main line 2e-3); the S0 n=0 box dropped (all other boxes kept); eps_FF = 2e-4 and 1e-3
  (main line 6e-5).  s0 variation is deferred until the Watson chains finish, because s0 lives in the hash-checked
  operator module and editing it would break the chains' restore step.
- Readings (no verdict changes; annotation of the "parameter dependence" the authors accept): UV +x end, P1 90-degree
  crossing, min |S_P1| below 1.2 GeV, S0 phase at 1 GeV, the FESR moments relative to their targets.  If loosening the
  S0 n=0 box (or dropping it) slows the S0 rise and/or moves the rho towards 820 MeV, the sum-rule tolerance is the
  parameter behind both discrepancies; if not, the FF caps are tested next by the eps_FF pair.
- Watson chains: after their 5 rounds each released chain is continued for 4 more rounds from its last leaf, to see
  whether the inward drift in (f00, f11) stops and whether the three chains approach one solution.

## 2026-09-18 04:00Z  (operator-18) Sensitivity: eps_SR at the tip (SR-a node reading)

- eps_SR 2e-3 -> 5e-3 -> 1e-2: UV +x end 0.07601 -> 0.07631 -> 0.07672 (+0.9 %); P1 crossing 773 -> 779 -> 785 MeV
  (+12 MeV); min |S_P1| 0.585 -> 0.571 -> 0.567; S0 crossing 694 -> 696 -> 706 MeV, S0(1 GeV) 150 -> 155 -> 154 deg.
  FESR moments: S0 n=0 +84 % -> +89 % -> +92 % of its target -- at 5e-3 and 1e-2 the box is no longer active (its edges
  would be +210 % and +419 %); the amplitude settles near +90 % on its own.  S0 n=1 +2 -> +4 -> +9 %, P1 n=-1 +13 ->
  +14 -> +17 %, P1 n=0 +1 -> +3 -> +7 %.
- Reading: the sum-rule tolerance is not the parameter behind either discrepancy: a five-fold loosening moves the rho
  by 12 MeV (towards the paper, but 30-40 MeV short) and leaves the S0 rise unchanged.  The S0 n=0 moment wants to be
  1.9 times the QCD value in this model whatever the box; with the 2e-3 box it is merely clipped at +84 %.
  Pending: the box dropped altogether, and the two eps_FF values.

## 2026-09-18 05:05Z  (operator-18) Pinned-point control at x_ref completed (5 rounds)

- P1 min |S| below 1.2 GeV by round: 0.244, 0.367, 0.409, 0.455, 0.505, 0.493 -- a plateau near 0.5 at the point
  (0.073321, -0.004652); the released chain from the same start reached 0.975 by round 5 while drifting to
  (0.06534, -0.00431).  Round 5's driver wall clock expired 250 s after SDPB's optimal termination (machine load ~120);
  the leaf was post-processed from its complete outputs (finish_leaf, note in report.json) and is accepted.
- Together with the tip control (fixed point at |S| = 0.593): at the boundary points of this finite problem the
  saturation the authors describe cannot be reached; it needs the point released.

## 2026-09-18 05:20Z  (operator-18) Sensitivity: S0 n=0 box dropped (tip)

- Accepted.  x_tip 0.07603 (main line 0.07601); P1 crossing 774 MeV (773), min |S_P1| 0.583 (0.585); S0 crossing
  unchanged, S0(1 GeV) 155 deg (150).  The freed S0 n=0 moment settles at about +90 % of the QCD value, the same value
  the loosened boxes gave.  The S0 n=0 sum rule is therefore inactive in effect: neither the rho position nor the S0
  rise depends on it.  What remains of the "parameter dependence" is eps_FF (two values queued) and s0.

## 2026-09-18 06:00Z  (operator-18) Sensitivity: eps_FF = 2e-4 at the tip -- the sensitive parameter

- eps_FF 6e-5 (paper's value) -> 2e-4, everything else the main line: UV +x end 0.07601 -> 0.08078 (paper 0.0811,
  now -0.4 %); S0 phases at 0.68/0.79/0.86/0.95/1.06/1.20 GeV 86/111/126/145/157/181 deg -> 63/73/82/88/111/152 deg
  against the paper's red curve 67/76/83/86/98/100 deg: the S0 wave now follows the paper to within 3-5 deg up to
  1 GeV (S0 crossing 694 -> 960 MeV).  P1: the rho moves from 773 to 972 MeV and becomes nearly elastic (min |S| 0.973);
  the paper has 813-827 MeV.  FESR: S0 n=0 +84 % -> +42 %, P1 n=-1 +13 % -> -22 %.
- Reading: the form-factor cap is the parameter behind both discrepancies.  The paper's stated 6e-5 with the factor at
  each node gives the fast S0 and a low rho; a 3.3 times looser cap gives the paper's S0 and +x end but a rho that is
  too high.  The paper's figures sit between the two, so its effective constraint on the form factors above s0 is
  looser than our node reading of (3.75) with 6e-5 and tighter than 2e-4 -- possibly with a different balance between
  the S0 and P1 caps (the ratio 2 m_q^2 : 1/2 in (3.75)).  Next: eps_FF 1e-4 and 1.4e-4 (both currents), and 2e-4 on
  one current at a time; eps_FF 1e-3 is already queued and will show the overshoot direction.  This is a sensitivity
  study, recorded as such: the main-line inputs and verdicts are not changed.

## 2026-09-18 08:50Z  (operator-18) Watson chain from the tip converged (8 rounds)

- Stop rule met at round 8 (watson_tip_cont round 3): max |h - t| over the two-pion-saturated nodes 0.0194 < 0.02.
  Trajectory (f00, f11): (0.07601, -0.00494) -> ... -> (0.06782, -0.00449), inward steps shrinking 2.8, 1.9, 1.6, 1.4,
  1.2, 1.0, 0.8, 0.6 %; P1 min |S| 0.585 -> 0.998; P1 crossing 773 -> 775 MeV; S0 crossing 694 -> 676 MeV.
  The converged amplitude is saturated and Watson-aligned, has the rho at 775 MeV (paper 813-827) and the S0 wave as
  fast as before.  The authors' weighting (ref start) gave the same rho (714 MeV) as the unit weighting (713).

## 2026-09-18 08:55Z  (operator-18) eps_FF scan: 1e-4 lands the rho just above the paper's; S0 needs a looser cap

- eps_FF 1e-4 (both currents): x_tip 0.07835; rho 849 MeV (paper 813-827), min |S_P1| 0.785; S0 at 0.79/0.86/0.95/
  1.06 GeV 93/107/123/146 deg (paper light pink 99/103/104/109, red 76/83/86/98), S0 crossing 776 MeV; FESR S0 n=0
  +66 %, P1 n=-1 -1 %.  Sequence so far (6e-5 / 1e-4 / 2e-4 / 1e-3): rho 773 / 849 / 972 / 1604 MeV, x_tip 0.0760 /
  0.0784 / 0.0808 / 0.0828, S0(0.95 GeV) 145 / 123 / 88 / (no crossing) deg.  By interpolation the paper's rho needs
  eps_FF about 8e-5 on the P1 current, while the paper's S0 needs about 2e-4 on the S0 current: no common value fits
  both, so the effective constraint of the paper differs between the two currents from the 2 m_q^2 : 1/2 split of (3.75).
- Pre-registered follow-up (hypothesis test, not a tuning of the main line): after the one-current runs (S0-only and
  P1-only 2e-4, queued) confirm that the two caps act separately, run the tip with eps_FF_S0 = 2e-4 and eps_FF_P1 = 8e-5.
  Reading: the hypothesis "the paper's effective caps are ~2e-4 (S0) and ~8e-5 (P1)" is supported if that run gives
  x_tip within 1 % of 0.0811, the rho crossing within 813-827 MeV and S0 within 10 deg of the paper's curves up to
  1 GeV -- all three at once; then the x_ref sections with the same pair test whether Fig.8's asymmetry appears.
  If not all three, the sensitivity study ends with the table above.

## 2026-09-18 10:15Z  (operator-18) eps_FF 1.4e-4 (recovered after a driver wall-clock expiry; SDPB optimal at 6523 s)

- x_tip 0.07964; rho 905 MeV, min |S_P1| 0.845; S0 at 0.79/0.86/0.95/1.06 GeV 82/92/105/130 deg; FESR S0 n=0 +53 %,
  P1 n=-1 -14 %.  The scan 6e-5 / 1e-4 / 1.4e-4 / 2e-4 / 1e-3 is monotone in every quantity: rho 773 / 849 / 905 /
  972 / 1604 MeV; x_tip 0.0760 / 0.0784 / 0.0796 / 0.0808 / 0.0828; S0(0.95 GeV) 145 / 123 / 105 / 88 / -- deg.
  (Correction: the two one-current tips had already started with the 7200 s budget when this was written; a
  wall-clock expiry after an optimal termination is recovered with finish_leaf, as above.)

## 2026-09-18 10:50Z  (operator-18) One-current cap: loosening the S0 cap alone fixes the S0 wave and leaves the rho untouched

- eps_FF_S0 = 2e-4 with eps_FF_P1 = 6e-5 (main line): x_tip 0.07758; rho 775 MeV, min |S_P1| 0.585 -- identical to the
  main line's P1; S0 at 0.79/0.86/0.95/1.06 GeV 75/84/90/111 deg against the paper's red curve 76/83/86/98 deg (within
  3 deg up to 0.95 GeV), S0 crossing 951 MeV; FESR S0 n=0 +44 %, P1 n=-1 +13 % (unchanged).
- The two caps act independently: the S0 form-factor cap sets the S0 rise (and the S0 n=0 moment), the P1 cap sets the
  rho position and its inelasticity.  The matched-pair test (S0 2e-4, P1 8e-5) is queued as pre-registered; the
  P1-only 2e-4 run is still solving.

## 2026-09-18 10:55Z  (operator-18) One-current cap, the other direction: loosening the P1 cap alone moves the rho and leaves S0 untouched

- eps_FF_P1 = 2e-4 with eps_FF_S0 = 6e-5: x_tip 0.07891; rho 967 MeV (both-currents 2e-4 gave 972), min |S_P1| 0.956;
  S0 at 0.79/0.86/0.95/1.06 GeV 110/125/148/164 deg -- the main line's fast S0, essentially unchanged.  Decoupling confirmed in
  both directions: S0 cap <-> S0 wave, P1 cap <-> rho position and inelasticity; x_tip responds to both
  (+0.0016 from the S0 cap, +0.0029 from the P1 cap at 2e-4).

## 2026-09-18 11:05Z  (operator-18) Comparison document extended (interim; chains still running)

- REPRODUCTION_COMPARISON_EN.{tex,pdf} regenerated with two new sections, both data-driven from the results root:
  "After the unitarity-saturation iteration" (Table tab:watson: last accepted round per chain, min|S_P1|, rho, S0 at 1 GeV,
  S2 at 1.2 GeV; pinned-control sentence; fig9_watson) and "Sensitivity to the form-factor cap" (Table tab:epsff: eps_FF scan,
  one-current rows, paper row; fig_epsff).  The request paragraph now asks how the caps on F0, F1 above s0 were normalised
  and for the iteration details.  15 pages, no overfull boxes.  Committed (branch sdpb-2309-regularised) and pushed.
- Tables will be regenerated once watson_ref_cont round 4, watson_mid_cont3 round 2 and the matched-pair tip
  (eps_FF S0 2e-4 / P1 8e-5) land; verdicts are not touched by these annotations.
- Status at 11:00Z: watson_ref_cont/round_04 SDPB at iteration 52 (gap 0.55, 3985 s of the 7200 s budget, load ~400);
  watson_mid_cont3/round_02 in the early iterations (9000 s budget); sens_epsFF_S0_2e-4_P1_8e-5_tip solving since 10:47Z.

## 2026-09-18 11:50Z  (operator-18) watson_ref_cont round 4: optimal at 6212 s, driver wall clock expired during the solution write; finish_leaf recovery started

- SDPB: "found primal-dual optimal solution", gap 6.2e-7, primal error 3e-22, dual error 1e-46, 81 iterations, 6212 s
  (load ~450; round 3 needed 6078 s at lower load).  Driver marked the leaf solver_failed at 11:49Z (wall clock 7260 s
  passed while the 9179 solution files were being written).  Same pattern as watson_ref_pinned round 5 and the
  watson_mid_cont rounds; recovered with scripts/sdp/finish_leaf.py (no re-solve; hash-checked restore, identical
  readback / Arb verification / acceptance logic).  Result is logged when the verification finishes.

## 2026-09-18 11:53Z  (operator-18) watson_ref_cont round 4 recovered and accepted; the x_ref chain ends after 9 rounds (5 + 4)

- round 4 (9th overall): f00 0.06266, f11 -0.00416 (round 0 of the chain 0.07332 / -0.00465: -14.5 % / -10.5 % in total,
  -0.9 % in this round); Watson objective 25.1117 (SDPB dual = Arb = float check); min |S_P1| 0.991; rho crossing 714 MeV
  (unchanged over the whole chain, paper 813-827); S0 at 1 GeV 135.8 deg; S2 at 1.2 GeV -18.1 deg.
- Convergence metric max |h - t| over two-pion-saturated nodes: 0.194 -> 0.172 -> 0.153 -> 0.137 -> 0.122 over the four
  rounds (about -11 % per round; the pre-registered tolerance 0.02 is not reached, extrapolation needs ~15 more rounds),
  so the chain is stopped by the round budget, not by convergence.  All-node metric stays 1.66 (the threshold node of P1,
  as noted on 09-17).  WATSON_ITER.json round-4 entry updated from the recovered leaf (accepted, metrics, note).
- No change to the picture: saturation is reached by moving inwards, the rho does not move, S0 stays fast.

## 2026-09-18 12:12Z  (operator-18) watson_mid_cont3 round 2 accepted; all three chains finished; post-Watson C6/C7 recorded

- x_ref+0.001 chain, 9th round overall: f00 0.06292, f11 -0.00419 (chain start 0.07432 / -0.00474); objective 24.9487;
  min |S_P1| 0.913; rho 741 MeV; S0 at 1 GeV 138 deg; S2 at 1.2 GeV -18.0 deg.  Saturated-node metric 0.119 -> 0.112 -> 0.112
  over the last three rounds: this chain has stopped improving its saturation while still drifting inwards (~ -0.9 %/round).
- Final leaves: tip watson_tip_cont/round_03 (8 rounds), x_ref watson_ref_cont/round_04 (9), mid watson_mid_cont3/round_02 (9).
- C67_POSTWATSON.json (scripts/sdp/c67_eval.py on the three final leaves, label post-Watson):
     C6 FAIL: crossings 775 / 714 / 741 MeV, spread 61 MeV, band [795, 845] not met; min eta >= 0.9 now TRUE at all three
              (0.998 / 0.991 / 0.913) -- the unitarity half of the rule is what the iteration fixes.
     C7 FAIL: S0 r.m.s. vs the paper's red curve 16.5 / 8.3 / 13.1 deg (before: 16.6 / 5.7 / 10.7); S2 r.m.s. 3.3 / 1.6 / 3.0
              (before 1.0 / 1.4 / 0.4); delta00(1.196) = 180.3 / 180.3 / 180.2 deg at all three points.
  Annotation only: C6/C7 main-line verdicts (C67_RESULT.json) unchanged, as pre-registered on 09-17 18:40Z.
- Remaining: the matched-pair tip (eps_FF S0 2e-4, P1 8e-5) is at iteration ~50; documents are regenerated once it lands.

## 2026-09-18 12:25Z  (operator-18) Matched-pair tip (eps_FF S0 2e-4, P1 8e-5) accepted: pre-registered criteria NOT all met (2 of 3)

- sens_epsFF_S0_2e-4_P1_8e-5_tip/tip: SDPB optimal, Arb verified.  x_tip 0.07855 (paper 0.0811: -3.15 %), f11 -0.00515;
  rho crossing 816.9 MeV (band 813-827: MET); min |S_P1| below 1.2 GeV 0.751 (unsaturated at the rho, before any iteration);
  S0 vs the paper's red curve up to 1 GeV: r.m.s. 1.9 deg, max 4.1 deg (criterion 10 deg: MET); S0 at 0.79/0.86/0.95/1.06 GeV
  74/83/90/111 (paper 76/83/86/98); S2 vs red: r.m.s. 0.3 deg, max 1.1 deg; S2 at 1.2 GeV -35.4 deg.
- Rule pre-registered on 09-18 (x_tip within 1 % of 0.0811 AND rho in 813-827 AND S0 within 10 deg up to 1 GeV, all three):
  NOT MET -- the +x end misses by 3 %.  Consequence as pre-registered: the x_ref sections with the same pair (Fig. 8 asymmetry
  test) are NOT run.  Reading: with the two caps decoupled, the pair reproduces the paper's Fig. 9 (red) and Fig. 10 (red)
  curves at the tip within digitisation accuracy, while the region's +x end stays 3 % short; since x_tip rises with either
  cap and the rho rises with the P1 cap, no pair of caps gives all three at once under our reading (per-node factor,
  everything else main line).  The residual 3 % therefore sits elsewhere (the eps_SR reading, the (3.75) normalisation
  above s0, or the authors' iteration acting on the region itself).  Annotation only: main-line inputs and verdicts unchanged.

## 2026-09-18 12:40Z  (operator-18) Post-reply programme complete; documents regenerated

- Done since the authors' reply of 09-17: saturation iteration implemented (eq. 2.29 form, unit and authors' weights, released
  and pinned), three chains of 8-9 rounds each finished, post-Watson C6/C7 recorded (both FAIL, unitarity half now passes),
  eps_SR scan (inactive), eps_FF scan and one-current tests (decisive, decoupled), matched pair (2 of 3 criteria; follow-up
  not triggered).  Verdicts C1-C8 unchanged; all main-line inputs unchanged.
- Regenerated from the results root: final_figures (fig9_watson, fig_epsff), REPRODUCTION_COMPARISON_EN.{tex,pdf,md,html}
  (sections 5-6 new, request renumbered 7; 15 pages), index.html (new section "作者回复后的补充"), REPORT_SDPB_2309_ZH.md
  section 8, docs README; receipts refreshed (C67_POSTWATSON.json, GATE_LOG.md, queue.log).
- Still open, deliberately: s0 sensitivity (requires editing the hash-checked operator module; nothing is running now, so it
  can be scheduled), and the authors' answers on the (3.75) normalisation and the eps_SR norm.

## 2026-09-18 13:20Z  (operator-18) Comparison document, second version (replaces the first, not incremental)

- REPRODUCTION_COMPARISON_EN.{tex,pdf,md,html} regenerated as a second version dated 18 September: "What was done" now
  states plainly that, after reading the 2403 / 2505 code and after Martin's replies of 17 September, three things were
  added (saturation iteration, eps_SR / eps_FF scan with separate caps, matched-pair test); Table 2 rows on (3.75) and on
  the boundary-point amplitude updated; a data-driven paragraph after the summary table; one sentence each in the Fig. 9
  and Fig. 10 subsections on the iterated amplitudes; Section 5 says what was taken from the two code versions and how it
  was adapted; Section 6 closes with the joint reading; the request names the three items that would settle the rest.
  15 pages, no overfull boxes.  Verdicts and main-line inputs unchanged.
- TODO_ZH.md created at the repository root (s0 sensitivity, the three answers awaited from the authors, share pin of the
  online ledger, o/operator-cache); linked from REPORT_SDPB_2309_ZH.md section 6 and the docs README.
