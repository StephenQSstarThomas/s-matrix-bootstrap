# Gate verdict (2026-09-13 14:55) — PASS under the pre-registered rules

Model: M=50, L=10, mixed-pv, chi-b, eps_chi=0.002, no UV, T0 free, SVD-reduced basis, SDPB 192 bit, gap 1e-6,
regulariser |rho_{a,ij}| <= Mreg (l-inf, unit-scaled 1x1 blocks).

Rule (1) (two-decade plateau within 1e-3): NOT met — the tip grows monotonically with Mreg
(1e2 0.08276, 3e2 0.08327, 1e3 0.08396, 1e4 0.08851, 1e5 0.10364, 1e6 0.14644).

Rule (2) (physics selection, fixed at 12:46 before the deciding runs):
  (i)  largest Mreg with the six first-omitted waves at eta <= 1.02 on every node: Mreg = 1e2
       (eta 1.011; 3e2 gives 1.046, 1e3 gives 1.145).
  (ii) L-stability at Mreg=1e2, M=50: L=8 0.082856, L=10 0.082757, L=12 0.082704 — spread 0.18 % (<= 2 %).
  => Mreg* = 1e2.

Post-hoc test against the paper (never an input): tip 0.082757 vs digitised 0.082573 = +0.22 % (band +-5 %);
f11/f00 = -0.0649 vs the chiral line -1/15 = -0.0667.
M-stability at Mreg*=1e2: M=30/L=8 tip 0.082671 (0.1 % from M=50), tip phase shifts at 0.9 GeV within 1 degree.
Norm control at M=30/L=8 (each norm at its own largest scale passing eta <= 1.02):
  l-inf Mreg=1e2: 0.082671, eta 1.007;  l2 B2=1800: 0.082836, eta 1.014 (B2=2600 fails at 1.021).
  Difference 0.20 % (<= 2 %); tip phases at 0.9 GeV within 0.5 degree (<= 2 degrees).  l4 control still running.

Finding to carry into the report: the Mreg=1e4 solution satisfies the authors' released l4 bound
(||rho||_4 = 5.8e4 < 3.775e5) yet exceeds the plotted endpoint by 7 %; the exact optimum of the authors' stated
problem is >= 0.0885, so the plotted 0.0826 reflects their double-precision solver, which implicitly stopped at the
scale where the omitted partial waves are unitary.  All numbers, times and the rule text are in GATE_LOG.md.
