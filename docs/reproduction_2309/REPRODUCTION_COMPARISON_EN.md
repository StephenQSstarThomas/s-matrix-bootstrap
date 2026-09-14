# Reproduction of arXiv:2309.12402 with SDPB: comparison against the paper's figures

Prepared for Yifei He and Martin Kruczenski. Date: 2026-09-14. Code and receipts: branch `sdpb-2309-regularised` of https://github.com/StephenQSstarThomas/s-matrix-bootstrap (folder `docs/reproduction_2309`).

## 1. What was done

We reproduced the finite S-matrix / form-factor bootstrap of *Bootstrapping gauge theories* (arXiv:2309.12402 v3) from the equations of Sections 2 and 3 of the paper, and compared the result with each figure of the paper. The optimiser is SDPB 3.1.0 (192-bit arithmetic, duality gap 1e-6, primal and dual error 1e-10); no other solver was used. The discretised operators (conformal map and grid (3.58)-(3.62), the cot kernel (3.67), the angular projections of (2.7), the current kernels (3.68)-(3.75)) were derived from the paper and checked independently: a Mathematica script recomputes the formulas at M=50, an independent reconstruction of the (2.7) projections agrees with the production operators to 5e-15 over all 3876 coefficients, and every accepted solution is re-verified against all original constraints in Arb interval arithmetic. The comparison thresholds for each figure were fixed before the corresponding runs and are recorded in the log (`receipts/GATE_PREREGISTRATION.md`, `receipts/GATE_LOG.md`).

The paper leaves five numerical items unstated. We list them here together with the readings we ran, because they decide the outcome for Figures 8-10.

| Item | Paper text | Readings run |
|---|---|---|
| Norm of the sum-rule tolerance (3.73), eps_SR = 2e-3 | "with some norm" | per-moment box on the raw moments (main line); per-wave L2 ball (structure of your 2403 code); relative 10 % per moment |
| Evaluation point of the factor between F and script-F in (3.75) | "which we evaluate at s = s0" appears in the estimate of eps_FF | factor at each node s_i > s0 (main line); factor frozen at s0 |
| Regularisation of the double spectral density | not mentioned in 2309; your 2103.11484 Section 3 and your 2403 code use an M-bound | |rho_ij| <= Mreg with Mreg fixed by an in-model rule (omitted-wave unitarity and L-stability), Mreg = 1e2 for the chiral and UV stages, 1e3 for the pure stage; l2 and l4 controls; Mreg x 10 control |
| Selection of the amplitude at a boundary point | "only points at the boundary have partial waves associated with them" | the solver's optimal point, plus a diagnostic of the whole near-optimal face (Section 4) |
| Chiral norm in (3.64) | "with some norm" | one combined 8-dimensional L2 norm (as in your 2403 code); two separate 4-dimensional norms as a control |

Everything else is taken as printed: M = 50, phi_i = (i - 1/2) pi / M, nu_0 = 0, L = 10 waves per isospin, s0 = (1.2 GeV)^2, alpha_s = 0.4, m_u = 4 MeV, m_d = 7.3 MeV, the condensates (2.54), m_pi = 140 MeV, f_pi = 92 MeV, the printed sum-rule numbers (2.56) times s0^(n+2), eps_chi = 2e-3, eps_FF = 6e-5, chiral points s = 1/2, 1, 3/2, 2, moments n = 0, 1 (S0) and -1, 0 (P1).

## 2. Summary table

| Figure | Paper statement | Our result | Agreement |
|---|---|---|---|
| Fig. 3 | pure-unitarity region, M=50, L=10 | +x end 2.22445 vs 2.23289 (-0.4 %); -x end -2.247 vs -2.902; f11 range [-0.6243, 0.06629] vs [-0.7336, 0.07932] | +x end agrees; our region is inside yours on the -x and +-y sides (regulariser-sensitive directions) |
| Fig. 4 | chiral constraints collapse the region onto f11 = -f00/15 | eps_chi = 2e-3: +x end 0.082757 vs 0.0825728 (+0.2 %); x_ref section width 0.0007618 vs 0.0007598; six-tolerance ladder monotone | agrees |
| Fig. 5 | subthreshold waves nearly linear; S0 chiral zero moves with eps_chi | RMS <= 6.4 % of f00(3) against the digitised curves at eps_chi = 2e-3, 4e-3, 6e-3; S0 zero at 0.426 / 0.293 / absent | agrees |
| Fig. 7 | chiral-only phases: S0, S2 agree with experiment, P1 has no rho | RMS vs digitised curves S0 0.388 deg, S2 0.144 deg, P1 0.456 deg; no P1 crossing below 1.2 GeV | agrees |
| Fig. 8 | with the sum rules the upper boundary shrinks much more than the lower | upper shrink / lower rise at x_ref: 0.000277 / 0.000277 (node factor), 5.32e-05 / 0.000249 (frozen factor); paper 2.32e-4 / 3.7e-5. +x end 0.076013 / 0.079623 vs 0.08112 | the asymmetry does not appear in any feasible reading |
| Fig. 9 | P1 phase shift with the rho, 90-degree crossing about 6 % above 770 MeV | crossings: tip 772.7 MeV (|S| min 0.585), mid 738.9 MeV (|S| min 0.077), ref: |S| dips to 0.24 at 0.73 GeV (crossing at 708 MeV on the +180 deg branch); paper 813-827 MeV | a rho-like structure appears at every point once the UV constraints are added; it sits 40-110 MeV lower than yours and is strongly inelastic |
| Fig. 10 | S0 and S2 coincide at low energy and spread above | S2: RMS 0.354-1.44 deg, endpoints -28 to -32 deg (paper -18 to -36); S0: RMS 5.7-16.6 deg, endpoint at 1.196 GeV 181-192 deg (paper 87-106) | S2 agrees; S0 rises faster than yours above 0.6 GeV |
| Fig. 11 | good convergence in L; rho peak shifts with M | P1 crossing at M=50: L=8/10/12 771.4 / 772.7 / 772.9 MeV; L=10: M=45/50/60 736.1 / 772.7 / 738.8 MeV; S0 at 1 GeV 150-165 deg (paper about 95) | L-stability agrees; the M shift is 37 MeV; S0 level differs |

## 3. Figure-by-figure comparison

In each figure below the left (or upper) panel is the paper's figure rendered from the arXiv source, the right (or lower) panel shows our accepted solutions with the digitised paper curve or boundary overlaid.

### Fig. 3, pure unitarity region

Paper: "we have used M=50 in (3.61) for discretization and imposed unitarity for 10 partial waves per isospin."

Ours: 24 support directions at 15-degree steps, all accepted, M=50, L=10, Mreg = 1e3 (chosen by the omitted-wave rule: the six first omitted waves stay within |S| <= 1.02 at every node and the L=8/10/12 endpoints agree to 2 %). Extrema against the digitised figure: f00 max 2.22445 vs 2.23289 (-0.4 %), f00 min -2.247 vs -2.902, f11 min -0.6243 vs -0.7336, f11 max 0.06629 vs 0.07932. The +x end is insensitive to Mreg (2.2116 / 2.2245 / 2.2305 at Mreg = 1e2 / 1e3 / 1e4). The three other extrema lie inside your region by 15-23 %; these are the directions where the amplitudes are large and the |rho_ij| bound is active. We read this as a measurement of the unstated regulariser rather than of the model; a Mreg = 1e4 control for these three directions is running.

![Fig. 3: paper (left) and our 24 support points with the digitised boundary (right).](figures/fig3.png)
*Fig. 3: paper (left) and our 24 support points with the digitised boundary (right).*

### Fig. 4, chiral constraints

Paper: "restricted by the chiral constraints (3.64) with tolerances 6e-3, 4e-3, 2e-3, 1e-3, 6e-4, 2e-4 (from the outer shape inward) ... with some norm".

Ours: one combined 8-dimensional L2 norm on the four-point residuals (the packaging of your 2403 code), Mreg = 1e2. At eps_chi = 2e-3 the +x end is 0.082757 (paper 0.0825728, +0.2 %) and the vertical section at x_ref has width 0.0007618 (paper 0.0007598); the +x ends of the six tolerances decrease monotonically (0.1608, 0.1255, 0.0828, 0.0551, 0.0411, 0.0223 for eps_chi = 6e-3 ... 2e-4). The two-norm control (chi-c) moves the +x end by less than the ladder spacing.

![Fig. 4: paper (left) and our six +x ends plus the eps_chi = 2e-3 sections with the digitised eps_chi = 2e-3 boundary (right).](figures/fig4.png)
*Fig. 4: paper (left) and our six +x ends plus the eps_chi = 2e-3 sections with the digitised eps_chi = 2e-3 boundary (right).*

### Fig. 5, subthreshold partial waves

Paper: "For the larger tolerances the partial waves are not approximately linear ... the chiral zero of the amplitude disappears for the blue points. ... The value eps_chi = 0.002 ... allows the physical value of f_pi".

Ours: the curves are evaluated from the Arb coefficients of the accepted leaves on 0.05 <= s <= 3.95, at the upper-branch x_ref representative of each tolerance. RMS against your digitised curves is at most 6.4 % of f00(3) (threshold 8 %); the S0 zero sits at s = 0.426 (2e-3), 0.293 (4e-3) and is absent at 6e-3, as in the figure.

![Fig. 5: paper's three panels (top) and ours (bottom).](figures/fig5.png)
*Fig. 5: paper's three panels (top) and ours (bottom).*

### Fig. 7, chiral-only phase shifts

Paper: "We choose a point closest to the black dot ... The partial waves at the magenta point and other nearby agree very well with experimental values for the S0 and S2 waves but not for the P1".

Ours: the representative is chosen by a rule fixed before any phase was read (nearest upper-branch endpoint to the black dot among the vertical sections at x_ref + {-0.002, -0.001, 0, +0.001, +0.002}). RMS against your digitised curves: S0 0.388 deg, S2 0.144 deg, P1 0.456 deg; P1 has no 90-degree crossing below 1.2 GeV.

![Fig. 7: paper (top) and ours (bottom).](figures/fig7.png)
*Fig. 7: paper (top) and ours (bottom).*

### Fig. 8, region after the sum rules and the form-factor bounds

Paper: "the plots are produced by taking eps_SR = 2e-3 in (3.73) and eps_FF = 6e-5"; caption: "The upper boundary shrinks notable but the lower not so much."

Ours, three readings of eps_SR and two of the (3.75) factor (all values at the x_ref section, paper values from the digitised figure):

| Reading | UV +x end | upper shrink | lower rise | ratio |
|---|---|---|---|---|
| raw box 2e-3 per moment, factor at each node | 0.076013 | 0.000277 | 0.000277 | 1 |
| raw box, factor frozen at s0 | 0.079623 | 5.32e-05 | 0.000249 | 0.214 |
| per-wave L2 ball 2e-3, factor frozen at s0 | 0.079541 | 5.35e-05 | 0.00025 | 0.214 |
| paper | 0.08112 | 2.32e-4 | 3.7e-5 | 6.3 |

The relative reading (10 % of each moment) is infeasible in this model. We proved this by removing the S0 n=0 box, keeping the other three 10 % boxes and every other constraint, and minimising that moment: the minimum is 0.0037454 = 1.57 x the QCD value, so no point of the feasible set lies within 10 % of it. In the raw-box reading the accepted solutions sit at the upper edge of three of the four boxes (S0 n=0 at +84 % of its target, S0 n=1 at +1.8 %, P1 n=0 at +1.4 %): the sum rules bind, and they bind in the direction of larger spectral integrals.

![Fig. 8: paper (left) and our chiral and chiral+UV boundary points with the digitised boundaries (right, with a zoom on the x_ref region).](figures/fig8.png)
*Fig. 8: paper (left) and our chiral and chiral+UV boundary points with the digitised boundaries (right, with a zoom on the x_ref region).*

### Fig. 9, P1 phase shift

Paper: "The resonance energy where the phase shift crosses pi/2 is slightly shifted from the real world data on the mass of the rho at 770 MeV by roughly 6 %." Digitised crossings: 0.827, 0.824, 0.813 GeV.

Ours, at the three representatives of the same frozen rule (tip = +x end; ref = nearest upper-branch point to the black dot; mid = its neighbour towards the tip): tip crossing 772.7 MeV with |S_P1| = 0.585 at 0.79 GeV; mid 738.9 MeV with |S_P1| down to 0.077; ref: |S_P1| falls to 0.24 at 0.73 GeV, and the node values admit two phase branches (the nearest-node lift swings to -47 deg, the +180 deg branch crosses 90 deg at 708 MeV). With the chiral constraints alone (Fig. 7) there is no crossing below 1.2 GeV, so the resonance is produced by the sum rules and the form-factor bounds, as in the paper; its position is 40-110 MeV below yours and the wave is far from elastic there. The paper does not show |S|; we add it in the second panel.

![Fig. 9: paper (left) and our P1 phases at the representatives with the three digitised curves (right).](figures/fig9.png)
*Fig. 9: paper (left) and our P1 phases at the representatives with the three digitised curves (right).*

![|S_P1| at the same points (not shown in the paper).](figures/fig9_eta.png)
*|S_P1| at the same points (not shown in the paper).*

### Fig. 10, S0 and S2 phase shifts

Paper: "the three bootstrap (red, pink, light pink) curves ... coincide at low energy but they spread at high energy."

Ours: S2 agrees with your curves (RMS 0.354-1.44 deg, endpoints -28 to -32 deg at 1.196 GeV against your -18 to -36 deg). S0 agrees below about 0.6 GeV and then rises faster: it crosses 90 deg at 0.69 GeV and reaches 181-192 deg at 1.196 GeV, where your curves are at 87-106 deg. With the chiral constraints alone our S0 matches your Fig. 7 to 0.4 deg, so the difference is introduced by the UV constraints.

![Fig. 10: paper (top) and ours (bottom).](figures/fig10.png)
*Fig. 10: paper (top) and ours (bottom).*

### Fig. 11, dependence on M and L

Paper: "For M=50, we have plotted the results with L=8, 10, 12. In all three partial waves there is a reasonably good convergence. ... the P1 phase shifts show a rho resonance with the peak slightly shifted depending on M."

| (M, L) | UV +x end | P1 crossing (MeV) | min |S_P1| below 1.2 GeV | S0 at 1 GeV (deg) |
|---|---|---|---|---|
| (50, 8) | 0.076176 | 771.4 | 0.584 | 151 |
| (50, 10) | 0.076013 | 772.7 | 0.585 | 150 |
| (50, 12) | 0.075944 | 772.9 | 0.585 | 150 |
| (45, 10) | 0.074275 | 736.1 | 0.226 | 165 |
| (60, 10) | 0.075097 | 738.8 | 0.597 | 156 |

The L dependence is 1.5 MeV in the crossing and 0.3 % in the endpoint. The M dependence is 37 MeV, with M=45 and M=60 both below M=50; in your figure the ordering is M=60 < M=45 < M=50. The S0 level at 1 GeV is 150-165 deg at every configuration against about 95 deg in your figure.

![Fig. 11: paper (top) and ours (bottom).](figures/fig11.png)
*Fig. 11: paper (top) and ours (bottom).*

## 4. Why the fine features of Figs. 8-10 do not follow from the stated problem

Two diagnostics were run on the accepted UV solutions. Both are exact up to the solver tolerance and were registered before the runs.

(a) Near-optimal face. At a boundary point with support direction d and optimal value v*, we add the constraint d.(f00, f11) >= v* - 2e-6 (twice the duality-gap threshold), keep the section, and maximise or minimise one node functional that is linear in the primal variables: 1 - Re S(s_i), Im S(s_i), Im F(s_i). The range of the functional over this slab shows how far the amplitude at that boundary point is determined by the extremal problem. Your Fig. 9 phases at the nodes 0.792 GeV (68.7-76.2 deg) and 0.864 GeV (112.8-124.1 deg) have cos 2 delta < 0, so any amplitude with those phases has 1 - Re S > 1 whatever |S| is.

| Representative | Functional | Range over the face | Consequence |
|---|---|---|---|
| tip (+x end) | 1 - Re S_P1 (0.792 GeV) | [1.441, 1.493] | pinned to +-0.03 around the tip's own value 1.470; below the 1.72-1.85 of your curves at |S| = 1 |
| tip | 1 - Re S_P1 (0.864 GeV) | max 0.6346 | below 1: your 112.8-124.1 deg excluded for every |S| |
| tip | Im S_P1 (0.792 GeV) | max -0.276 | negative on the whole face: every amplitude there is past the resonance at 0.79 GeV |
| x_ref upper (ref) | 1 - Re S_P1 (0.792 GeV) | [0.02279, 0.9356] | wide (undetermined), yet below 1 everywhere: your 68.7-76.2 deg excluded for every |S| |
| x_ref upper | 1 - Re S_P1 (0.864 GeV) | max 0.701 | below 1: excluded for every |S| |
| x_ref upper | Im F_1 (0.792 GeV) | [0.07843, 3.71] | the vector form factor is essentially free at this point |
| tip | 1 - Re S_S0 (0.792 GeV) | [1.59, 1.726] | the S0 phase is fixed at 111 deg on the face; only |S| varies (0.79-0.99) |

The tip face is rigid in the P1 observables and does not contain your P1 shape; the x_ref face is degenerate (the paper's premise that a boundary point carries one set of partial waves does not hold there) and still excludes your P1 phase at both nodes. Choosing a different point on these faces, or a different solver, cannot reproduce Fig. 9 from this constraint set.

(b) Regulariser scale. Raising Mreg from 1e2 to 1e3 moves the UV +x end by +0.9 %, the P1 crossing by -3.6 MeV and min |S_P1| by -0.025; the Fig. 8 asymmetry does not appear at either scale.

![Ranges of the node functionals over the near-optimal faces (bars) against the floor 1 - Re S = 1 required by the paper's phases (dashed).](figures/face_ranges.png)
*Ranges of the node functionals over the near-optimal faces (bars) against the floor 1 - Re S = 1 required by the paper's phases (dashed).*

## 5. Questions

The following items would let us close the gap or confirm that it is real. We have your 2403 and 2505 codes and can see how they treat these points, but we are asking specifically about the 2309 runs.

1. Which norm was used in (3.73), and is the QCD value compared with the raw moment or with the normalised one (divided by s0^(n+2))? The raw per-moment box is the only reading with eps_SR = 2e-3 that is feasible in our implementation.
2. In (3.75), is the factor between F and script-F evaluated at each node above s0 or frozen at s0?
3. Was a bound on the double spectral density rho_ij (as in 2103.11484 Section 3, or the l4 bound of your 2403 code) applied in the 2309 runs, and at what scale? Without any bound the finite problem is ill-posed at 192-bit precision; with our in-model rule the +x ends agree with your Figs. 3 and 4 to 0.4 %.
4. Were the Fig. 9-10 amplitudes the solver's optimal point of the support problem, or was a saturation / Watson step (as in your later code) already used to select the amplitude at each of the three points?
5. Could you share the (f00(3), f11(3)) coordinates and |S_P1| of the red, pink and light pink points, and the value of s0 used for Figs. 8-11?

## 6. Receipts

`receipts/C1_RESULT.json` ... `C8_RESULT.json` (per-figure verdicts on the pre-registered rules), `receipts/FACE_RESULT.json` (face diagnostic), `receipts/UV_REPRESENTATIVE_SELECTION.json` and `IR_REPRESENTATIVE_SELECTION.json` (frozen representative choices), `receipts/GATE_LOG.md` (time-stamped log of every run and decision). Each accepted leaf keeps its PMP, the SDPB output and the Arb verification (`report.json`); these are available on request (several GB).
