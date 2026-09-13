# legacy/ — retired Newton mainline (evidence only)

`smatrix_bootstrap_newton/` is the former `src/smatrix_bootstrap/*.py` package (23 modules, self-written
barrier–Newton solver, barrier-centre representative points, analytic-cardinal / finite-sine / T0=0 / five-tail
prescriptions) that produced the results under `results/runs/{stage_*,mainline_alignment_20260909,
physical_consistency_20260911,sequential_reproduction_20260910,major_claims_20260912,...}`.  It is kept so that
those records can be re-verified; it is **not** part of the SDPB reproduction and is not installed.  Why it was
retired, module by module, is in `../TRIAGE_ZH.md` §2.2.

Still reusable from it as independent comparison targets: `kernels.py` (PV rows, exterior Legendre-Q),
`operators.py` lines 207–312 (grid, angular kernels), `model.py` lines 29–172 (unitarity, kinematics, Gram, FESR,
FF operators).  `src/smatrix_bootstrap/sdp/crosscheck.py` loads them from here when they are present.

Run its own tests (107, ~35 s, need python-flint) with `python -m pytest legacy/tests`.
To run the retired CLI: `PYTHONPATH=legacy python -m smatrix_bootstrap_newton.run --help` (requires the
`legacy` extra: highspy, scs, clarabel, cvxpy).
