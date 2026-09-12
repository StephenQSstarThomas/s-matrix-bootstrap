#!/usr/bin/env python
"""Resolution ladder: the certified extremum against M, for one constraint set.

The paper's figures are at M = 50, L = 10.  Where that resolution does not
converge, quoting a result at a smaller M is only meaningful next to the trend,
so the ladder is an artefact in its own right rather than a footnote.

Every entry is certified by constraint generation: the relaxation optimum whose
returned point satisfies all 3 L M unitarity disks of the full problem.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings

warnings.filterwarnings("ignore")

from smatrix_bootstrap.sdp.problem import ModelSpec
from smatrix_bootstrap.sdp.runner import solve_generated
from smatrix_bootstrap.sdp.verify import full_report

GRID = [(20, 6), (25, 8), (30, 8), (35, 10), (40, 10), (45, 10), (50, 10)]
TARGETS = {"pure": ("Fig. 3 +x tip", 2.2328885260009974),
           "chiral": ("Fig. 4 +x end at eps=2e-3", 0.0826)}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["pure", "chiral", "uv"], default="pure")
    p.add_argument("--direction", nargs=2, type=float, default=[1.0, 0.0])
    p.add_argument("--eps-chi", type=float, default=2e-3)
    p.add_argument("--out", default=None)
    p.add_argument("--grid", default=None,
                   help="comma list of M:L, e.g. 20:6,30:8,40:10")
    p.add_argument("--time-limit", type=float, default=3600.0)
    a = p.parse_args()

    label, target = TARGETS.get(a.mode, (None, None))
    grid = ([tuple(int(x) for x in g.split(":")) for g in a.grid.split(",")]
            if a.grid else GRID)
    rows = []
    for M, L in grid:
        t0 = time.time()
        sp = ModelSpec(M=M, L=L, cone_scaling="rownorm", B=3.775e5, B_norm="l4",
                       chiral=a.mode in ("chiral", "uv"), eps_chi=a.eps_chi,
                       uv=a.mode == "uv")
        res, model, rounds = solve_generated(sp, tuple(a.direction), max_iter=800,
                                             time_limit=a.time_limit)
        row = {"M": M, "L": L, "status": res["status"], "seconds": time.time() - t0,
               "rounds": len(rounds),
               "relaxation_upper_bounds": [r.get("objective") for r in rounds]}
        if res.get("f00_3") is not None:
            v = full_report(model, model.solution())
            row.update({"objective": res["objective"], "f00_3": res["f00_3"],
                        "f11_3": res["f11_3"], "certified": res.get("certified"),
                        "n_disks_imposed": res.get("n_disks_imposed"),
                        "feasible": v["unitarity"]["feasible"],
                        "c_norm_inf": v["c_norm_inf"], "rho_l4": v["rho_l4"],
                        "B_active": v.get("B_active")})
            if target:
                row["rel_to_paper"] = res["objective"] / target - 1.0
        rows.append(row)
        print(json.dumps(row, default=float), flush=True)
        if a.out:                      # write after every M so a long ladder is
            _dump(a, grid, label, target, rows)   # usable while it is still running
    if a.out:
        _dump(a, grid, label, target, rows)
    return 0


def _dump(a, grid, label, target, rows):
    doc = {"mode": a.mode, "direction": a.direction, "eps_chi": a.eps_chi,
           "grid": grid, "paper_reference": label, "paper_value": target,
           "rows": rows}
    with open(a.out, "w") as fh:
        json.dump(doc, fh, indent=1, default=float)


if __name__ == "__main__":
    raise SystemExit(main())
