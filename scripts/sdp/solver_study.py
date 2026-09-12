#!/usr/bin/env python
"""Reproduce the solver-behaviour study quoted in the report (task stop rule (c)).

The conic programme of section 3 is, at the paper's resolution, at the edge of
what Clarabel handles.  This script records the evidence, so the claim
"the solver, not the model, is the limiting factor" is checkable rather than
asserted.  It runs at M = 20 (seconds per solve) where every configuration can
be compared against a common reference.

Reported per configuration: the objective, whether the returned point passes the
a posteriori feasibility check on the unmodified operators, and ||rho||_4.

The monotonicity test is the decisive one: enlarging the density ball can only
enlarge the feasible set, so the optimum must be non-decreasing in B.  Where the
reported optimum decreases instead, the solver is returning a suboptimal point.
"""
from __future__ import annotations

import argparse
import json
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from smatrix_bootstrap.sdp.problem import Model, ModelSpec
from smatrix_bootstrap.sdp.verify import full_report


def one(M, L, tag, **kw):
    opts = {k: kw.pop(k) for k in list(kw) if k.startswith("equilibrate")}
    meta = {"B": kw.get("B"), "B_norm": kw.get("B_norm", "l2"),
            "reduce_basis": kw.get("reduce_basis", False)}
    m = Model(ModelSpec(M=M, L=L, cone_scaling=kw.pop("cone_scaling", "rownorm"), **kw))
    m.finalize()
    r = m.solve((1.0, 0.0), max_iter=1000, time_limit=900.0, **opts)
    out = {"tag": tag, "status": r["status"], "seconds": r["seconds"],
           "iterations": r.get("iterations"), "n_vars": r.get("n_reduced"), **meta}
    if r.get("f00_3") is None:
        return out
    v = full_report(m, m.solution())
    out.update({"f00_3": r["f00_3"], "feasible": v["unitarity"]["feasible"],
                "max_rel_violation": v["unitarity"]["max_relative_violation_active"],
                "rho_l4": v["rho_l4"], "rho_l2": v["rho_l2"]})
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--M", type=int, default=20)
    p.add_argument("--L", type=int, default=6)
    p.add_argument("--out", default=None)
    a = p.parse_args()
    rows = []
    for scale in ("none", "centrifugal", "rownorm"):
        rows.append(one(a.M, a.L, f"cone_scaling={scale}, B=377500 l4",
                        cone_scaling=scale, B=3.775e5, B_norm="l4"))
    for B in (1e4, 1e5, 3.775e5, 1e6, 1e7, 1e9):
        rows.append(one(a.M, a.L, f"B={B:.3e} l4", B=B, B_norm="l4"))
    for B in (3.775e5, 1e7, 1e9):
        rows.append(one(a.M, a.L, f"B={B:.3e} l2", B=B, B_norm="l2"))
    rows.append(one(a.M, a.L, "no B, exact basis reduction",
                    B=None, reduce_basis=True, basis_tol=1e-12))
    rows.append(one(a.M, a.L, "no B, full space", B=None))

    fine = [r for r in rows if r.get("feasible")]
    best = max(fine, key=lambda r: r["f00_3"]) if fine else None
    # monotonicity is only meaningful within one norm and one parametrisation
    l4 = [r for r in rows if r.get("B_norm") == "l4" and r.get("B") is not None
          and not r["reduce_basis"] and r.get("f00_3") is not None
          and r["tag"].startswith("B=")]
    l4.sort(key=lambda r: r["B"])
    mono = all(l4[i]["f00_3"] <= l4[i + 1]["f00_3"] + 1e-9 for i in range(len(l4) - 1))
    doc = {"M": a.M, "L": a.L, "rows": rows,
           "monotonicity_series": [{"B": r["B"], "f00_3": r["f00_3"]} for r in l4],
           "best_verified_feasible": best,
           "objective_is_monotone_in_B": mono,
           "note": "the optimum must be non-decreasing in B; a decrease proves the "
                   "solver returned a suboptimal point"}
    print(json.dumps(doc, indent=1, default=float))
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(doc, fh, indent=1, default=float)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
