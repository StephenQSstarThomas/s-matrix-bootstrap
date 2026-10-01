"""G2 (plan Task A10): cost and outcome of each trial encoding.  Usage: python scripts/gtb/g2_report.py LEAF_DIR..."""
import json
import re
import sys
from pathlib import Path


def iters(log):
    rows = [l for l in Path(log).read_text().splitlines() if re.match(r"^\d+\s+\d+\s", l)]
    if not rows:
        return 0, None
    last = rows[-1].split()
    return int(last[0]), float(last[1])


out = {}
for d in map(Path, sys.argv[1:]):
    r = json.loads((d / "report.json").read_text())
    n, t = iters(d / "sdpb.log") if (d / "sdpb.log").exists() else (0, None)
    s, o = r.get("sdpb_result", {}), r.get("observables", {})
    R = r["spec"]["reg_bound"] * 3775
    out[d.name] = {"reg": r["regulariser"], "reduced": r["reduce_basis"], "precision": r["settings"]["precision"],
                   "nproc": r["settings"]["nproc"], "n_vars": r["n_vars"], "blocks": r["pmp"]["n_blocks"],
                   "pmp_GB": round(r["pmp"]["bytes"] / 1e9, 2), "pmp2sdp_s": round((r.get("pmp2sdp") or {}).get("seconds", 0)),
                   "iterations": n, "sdpb_s": t, "s_per_iter": round(t / n, 1) if n and t else None,
                   "terminate": s.get("terminateReason"), "F1_primal": s.get("dualObjective"),
                   "F1_dual": s.get("primalObjective"), "gap": s.get("dualityGap"),
                   "rho_linf": o.get("rho_linf"), "rho_l4": o.get("rho_l4"),
                   "l4_feasible_float": (o.get("rho_l4") or 1e99) <= R}
print(json.dumps(out, indent=1))
