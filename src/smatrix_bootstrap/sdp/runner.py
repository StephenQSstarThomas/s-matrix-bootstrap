"""Job runner: one solve -> one ``report.json`` from which every number recomputes.

Each record carries the argv, the input hashes, the full model spec, the solver
status/iterations/residuals, the objective, the a posteriori verification of the
*unmodified* constraints, and the solution hash.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import asdict

import numpy as np

from . import constraints as C
from .observables import resonance_report, subthreshold_curves
from .problem import Model, ModelSpec
from .verify import full_report

CODE_FILES = ("grid.py", "legendreq.py", "hilbert.py", "projector.py",
              "formfactor.py", "constraints.py", "problem.py", "verify.py",
              "runner.py", "observables.py")


def code_hash() -> dict:
    here = os.path.dirname(__file__)
    out = {}
    for f in CODE_FILES:
        p = os.path.join(here, f)
        if os.path.exists(p):
            out[f] = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
    return out


def provenance() -> dict:
    try:
        rev = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                      cwd=os.path.dirname(__file__),
                                      stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        rev = None
    import clarabel, cvxpy, scipy
    return {"argv": sys.argv, "git_rev": rev, "python": sys.version.split()[0],
            "numpy": np.__version__, "scipy": scipy.__version__,
            "cvxpy": cvxpy.__version__, "clarabel": clarabel.__version__,
            "host": platform.node(), "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "code_sha256": code_hash()}


def run_job(spec: ModelSpec, jobs, outdir: str, solver: str = "CLARABEL", **kw) -> dict:
    """``jobs``: list of ``(name, direction, extra_constraint_builder)``.

    ``extra`` is called as ``extra(model, ctx)`` where ``ctx`` accumulates the
    results already obtained in this process, so the representative points can
    chain: ``tip`` is solved first and ``ref``/``mid`` read ``ctx["tip"]`` for
    the section they are constrained to.
    """
    os.makedirs(outdir, exist_ok=True)
    model = Model(spec)
    records, ctx = [], {}
    built_with_extra = None      # None = never built
    for name, direction, extra in jobs:
        # The objective direction is a cvxpy Parameter, so a pure direction
        # sweep canonicalises once and only re-solves.  A job that adds an extra
        # constraint has to rebuild, and so does the first plain job after one.
        if extra is not None or built_with_extra is not False:
            model.finalize(extra(model, ctx) if extra else ())
            built_with_extra = extra is not None
        t0 = time.time()
        res = model.solve(direction, solver=solver, **kw)
        rec = {"job": name, "spec": asdict(spec), "result": res,
               "wall_seconds": time.time() - t0}
        sol = model.solution()
        if sol.get("c") is not None:
            rec["verification"] = full_report(model, sol)
            rec["observables"] = resonance_report(model.ops, sol["c"])
            rec["subthreshold"] = subthreshold_curves(
                model.ops, sol["c"], np.linspace(0.05, 3.95, 40))
            np.savez_compressed(os.path.join(outdir, f"sol_{name}.npz"), **{
                k: v for k, v in sol.items() if v is not None})
        ctx[name] = res
        records.append(rec)
        _write(outdir, records)
    return {"records": records}


def _write(outdir, records):
    doc = {"provenance": provenance(), "records": records,
           "fesr_target_audit": C.fesr_target_audit()}
    with open(os.path.join(outdir, "report.json"), "w") as fh:
        json.dump(doc, fh, indent=1, default=float)


def sweep_directions(n: int, half: bool = False):
    """``n`` unit directions; ``half`` keeps only the upper half plane."""
    ang = np.linspace(0.0, np.pi if half else 2 * np.pi, n, endpoint=half)
    return [(float(np.cos(a)), float(np.sin(a))) for a in ang]
