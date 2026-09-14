# Gate experiment: regularised M=50 chiral tip (pre-registered before any result)

Written 2026-09-13 before the first solve finished.

Model (identical in the four runs except Mreg): M=50, L=10, mixed-pv (the authors'
nodal discretisation: cot kernel, midpoint rule, Legendre-Q), nu0=0, chiral only,
chi-b (one combined 8-dim L2 norm), eps_chi = 0.002, no UV, T0 free, reduce-basis
(SVD tol 1e-12), regulariser |rho_{a,ij}| <= Mreg (2103.11484 sec. 3 M-regularisation),
SDPB 192 bit, operator rows 40 digits, printed 30 digits, gap 1e-6, 8 ranks.
Objective: maximise f00(3) (the +x tip of the Fig. 4 eps=0.002 region).

Runs: Mreg in {1e3, 1e4, 1e5, 1e6}.

Prediction (falsifiable): a plateau exists over at least two consecutive decades of
Mreg on which the tip agrees within 1e-3 absolute, and the plateau value is
0.0826 +- 0.004 (5 % of the digitised paper value 0.08257; the repo's own MOSEK
ladder gave 0.08246-0.08260 for M=25..40).  On the plateau ||rho||_4 is expected
far below 1e6 (the un-regularised SDPB solutions had 1e16-1e17).

Pass: prediction holds.  Fail: no plateau, or plateau value outside the band.
On fail the root-cause judgement is refuted and the plan must be revised before
any further physics runs.  Nothing here is tuned to phase shifts or to Fig. 4
beyond the pre-stated acceptance band.
