# Scientific mainline — He–Kruczenski 2309.12402v3

## Scientific progression and lessons from B

- Follow the paper's principles, constraints and scientific comparisons.
  Exact recovery of every displayed point is not an acceptance requirement.
  Declare underspecified numerical choices and retain measured differences;
  do not tune them to the published or experimental curves.
- Reuse existing amplitudes and operators before launching optimization.
  State the next physical decision and the smallest calculation that answers it.
  Do not turn invented regulator platforms or continuum certificates into
  prerequisites for the paper's finite prototype.
- Default to one scientific calculation at a time with a bounded runtime.
  Add concurrency only for concrete independent work; avoid broad parameter,
  solver or boundary scans without a mainline reason.
- Fix representative-point selection before inspecting phase shifts.
  Preserve the same complete amplitude through the linked figures and report
  approximation choices explicitly. Never relax fundamental unitarity to match a plot.
- Use published output curves for comparison, not for choosing new physical
  inputs, uncertainty norms or representative coordinates. Label historical
  output-assisted choices; do not retroactively call them independent.
- Distinguish an exact analytic function family from a mixed collocation
  approximation. Native matrix certificates do not establish global
  amplitude existence or transfer to another interpolation prescription.

## Repository and evidence

- Reproduce Fig.3–11 of `references/2309.12402v3.pdf`.
  Snowmass is background; scalar calibration and completed rank proofs
  do not replace the pion/gauge results. Map work to STATUS steps A1–F3.

- Start with `STATUS.md` and `SCIENCE.md`. Read only needed core files.
  Read `REPRODUCTION_GUIDE_ZH.md` for interfaces and dependencies.
  Do not scan historical results or compressed source snapshots by default.
- Keep a small set of modules organized by scientific responsibility.
  A few additional modules are allowed; ten is not a hard limit.
  Keep joint gauge optimization and figure delivery behind the same CLI;
  avoid redundant branches, experiment-specific scripts and architecture.
- Each core file remains limited to350 lines and24KiB. Split only at a
  clear scientific responsibility; keep calculations behind the same CLI.
- Keep at most three active test files, each at most350 lines. Preserve
  independent mathematical checks; avoid redundant implementation-mirroring tests.
  Delete redundant legacy files only
  after checking their historical proof/source recovery records.
- Use `python -m smatrix_bootstrap.run` as the single calculation entry.
  Runs write data and provenance to `results/runs`, not new source files.
- Focus on derivations, source fidelity and the next physical decision.
  A smaller passing model never replaces the requested published result.
- Preserve historical proof inputs and failures. The finite-sine and
  PV/midpoint prescriptions have distinct scopes; do not transfer rank or
  feasibility results without proof.
- Keep the removed legacy tree removed. Essential derivations live in
  `results/evidence`; original producer bytes remain in authenticated
  non-executable snapshots. Do not restore old modules to run experiments.
- Record density/infinity assumptions, chiral norm, FESR units and cutoff,
  and form-factor normalization explicitly. Do not tune unresolved
  prescriptions to the published or experimental curves.
- Sampled feasibility is not continuum certification. Published regions,
  phase curves, convergence and support verification remain required.
