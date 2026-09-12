# SDP-route reproduction of arXiv:2309.12402v3 — run of 2026-09-12

Produced by the new subpackage `src/smatrix_bootstrap/sdp/`, which re-derives
every operator from the paper's equations and shares no code with the
repository's historical Newton mainline.

## What is here

| file | content |
|---|---|
| `preregistration.json` | the three under-determined calibers and the C1–C8 acceptance thresholds, frozen before P3 |
| `reports/**/report.json` | one per solver invocation: argv, git rev, code hashes, model spec, solver status/iterations/seconds, objective, a posteriori verification, solution hash |
| `manifest.json` | index from every reported number to the `report.json` that produced it |
| `verdicts.json` | the C1–C8 verdicts computed by `smatrix_bootstrap.sdp.claims` |
| `tables.json` | comparison tables against `references/figure*.csv` |
| `solver_study_M20.json` | the solver-behaviour study behind the stop-rule (c) finding |
| `figures/` | figure PDFs |
| `REPORT_SDP_ZH.md` | the report |

Heavy artefacts (solution vectors, per-worker logs) stay on the scratch disk at
`/playpen1/shiqiu/sdp-work/runs/sdp_reproduction_20260912/` and are not copied
into the repository.

## How to recompute

```bash
export PYTHONPATH=/tmp/collocation_arb:src          # python-flint for the Arb audit
python -m smatrix_bootstrap.sdp selfcheck           # the section 5a self-checks
python -m smatrix_bootstrap.sdp prereg  --out <run>
python scripts/sdp/drive.py fig3        --root <run> --jobs 8
python scripts/sdp/solver_study.py      --M 20 --L 6 --out <run>/solver_study_M20.json
python scripts/sdp/finalize.py          --root <run> --repo-out results/runs/sdp_reproduction_20260912
```

Nothing in this directory overwrites any pre-existing `results/runs` entry.
