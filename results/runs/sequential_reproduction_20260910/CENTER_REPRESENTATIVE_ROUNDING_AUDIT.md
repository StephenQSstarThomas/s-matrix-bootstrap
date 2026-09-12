# Terminal representative rounding audit

Created: 2026-09-11T05:39:03.846917+00:00

All six requested terminal representatives preserve every C, ImF and moment-bounded R coordinate exactly from the reconstructed `JointProblem.raw(saved z)`. Each direct original audit passes after its permitted free-spectrum Schur lift and reproduces the delivered `joint.npz` point exactly. No delivered center flag in this six-run set is affected by segment restoration or seed fallback.

| Run | M / L | Changed C / ImF / bounded R | Free R entries lifted | Maximum free R lift | Last Newton decrement |
|---|---:|---:|---:|---:|---:|
| E1_tip_center | 50 / 10 | 0 / 0 / 0 | 14 | 0.0553704712396 | 1.81893e-14 |
| E2_ref_01 | 50 / 10 | 0 / 0 / 0 | 14 | 0.790725421909 | 1.61532e-11 |
| E3_mid_02 | 50 / 10 | 0 / 0 / 0 | 14 | 0.90754658928 | 1.23835e-10 |
| F_M50_L8_support_05 | 50 / 8 | 0 / 0 / 0 | 14 | 0.00613258886566 | 3.4083e-16 |
| F_M50_L12_support_02 | 50 / 12 | 0 / 0 / 0 | 14 | 0.104779879054 | 7.48651e-13 |
| F_M45_L10_support_01 | 45 / 10 | 0 / 0 / 0 | 14 | 0.152420781628 | 2.84263e-10 |

Each reconstruction uses its saved preparation H and energies, current matrices, dimensions, bound, epsilon, separate-L2 norm, objective/section, and `reference_point` from `barrier_state.npz`, with `absorptive=True`. The constructor, raw transform, absorptive change, original audit and restoration function were compared to each archived producer by AST hash. The JSON records those results and SHA-256 hashes for all input files. All calculations were sequential, single-thread BLAS, nice 10, and read-only except these two reports. The three E producers predate the optional analytic-preparation extension in `joint_audit`; its full-function hashes differ, while the normalized None/float native path matches. Replayed points and both support bounds equal the saved results in all six cases. No archived code was executed.

The saved checkpoint mu matches the final Newton entry, final center history, and final center-path mu in every case. The final center-path point equals the delivered point exactly. The 14 differing coordinates are outside all moment supports and are the free high-energy R values removed from the numerical barrier. The original auditor restores their Schur floors without changing the active center. The active segment fit is the identity endpoint (fraction 1 with zero residual); direct audit success establishes that restoration/fallback was unnecessary.

This establishes provenance of the numerical center claims. It does **not** independently recompute the delivered point’s barrier decrement or provide a mathematical center certificate. The Newton decrement was logged before the terminal step, in the transformed representation with incrementally propagated values. Exact agreement of the saved active float64 coordinates is stronger lineage evidence than a small coordinate difference, but feasibility and support optimality alone still do not prove centeredness.

The latent control-flow defect remains in the audited producer: `linear.center_joint_support` restores an infeasible candidate or falls back to the run input seed, then can retain `converged=True` and record that changed point in `center_path` and `representative`. Fix the delivered-center eligibility flag when active coordinates change, keep the numerical convergence and finite-support records separate, and require actual re-centering before promoting a changed point. Permit free high-R lifts because those coordinates are absent from the barrier. Regression checks should cover direct Schur lift, segment recovery, seed fallback, and checkpoint/representative mismatch.

Detailed evidence and input hashes: [CENTER_REPRESENTATIVE_ROUNDING_AUDIT.json](CENTER_REPRESENTATIVE_ROUNDING_AUDIT.json).
