"""Re-run the post-solve tail of ``sdpb.run_once`` on a leaf whose SDPB finished but whose driver died afterwards.

Usage: python scripts/sdp/finish_leaf.py LEAF_DIR --note "why the driver died"

Nothing is re-solved: the saved pmp.json/basis.npy are restored through ``Pmp.from_saved`` (hash-checked), out/out.txt and
out/y.txt are read as they are, and the same readback, Arb verification, convergence and acceptance logic as run_once is
applied.  The report keeps the original assembly record and gains ``postprocess`` with the reason and timestamp.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from smatrix_bootstrap.sdp.pmp import Pmp
from smatrix_bootstrap.sdp.sdpb import Settings, convergence_report, read_out, write_json


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("leaf"); p.add_argument("--note", required=True)
    a = p.parse_args(argv)
    dest = Path(a.leaf).resolve(); report = dest / "report.json"
    rec = json.loads(report.read_text())
    out_txt, yp = dest / "out/out.txt", dest / "out/y.txt"
    wall_clock = (rec.get("status") == "solver_failed" and rec.get("sdpb", {}).get("returncode") == 124
                  and out_txt.exists() and "found primal-dual optimal solution" in out_txt.read_text())
    if (rec.get("status") not in ("solving",) and not wall_clock) or "verification" in rec:
        raise SystemExit(f"{dest}: status {rec.get('status')!r}; only an unfinished 'solving' leaf, or one whose driver "
                         "wall clock expired after SDPB had terminated optimally, is post-processed")
    if not (out_txt.exists() and yp.exists() and yp.stat().st_size):
        raise SystemExit("SDPB outputs incomplete; nothing to post-process")
    if wall_clock:
        rec.pop("solver_failure", None); rec.pop("failure_trace", None); rec.pop("untrusted_partial_output", None)
    cfg = Settings(**rec["settings"])
    proc = json.loads((dest / "sdpb_process.json").read_text())
    rec["sdpb"] = {"returncode": None, "command": proc.get("command"),
                   "seconds": out_txt.stat().st_mtime - proc.get("started_unix", out_txt.stat().st_mtime),
                   "scope": ("driver wall clock expired after SDPB had terminated optimally and written out.txt and y.txt" if wall_clock
                             else "driver died after SDPB wrote out.txt and y.txt; return code not recorded, outputs complete")}
    w = Pmp.from_saved(str(report), tuple(rec["direction"]), rec.get("fix_f00"), rec.get("face"), rec.get("functional"))
    rec["sdpb_result"] = read_out(out_txt)
    sol = w.read_solution(yp)
    rec["y_sha256"] = hashlib.sha256(yp.read_bytes()).hexdigest()
    np.savez_compressed(dest / "solution.npz", **{k: v for k, v in sol.items() if k != "y_text"})
    rec["solution_sha256"] = hashlib.sha256((dest / "solution.npz").read_bytes()).hexdigest()
    rec["verification"] = w.verify(sol)
    rec["convergence"] = convergence_report(rec["sdpb_result"], 0, rec["verification"].get("objective_recomputed"), cfg)
    rec["accepted"] = bool(rec["convergence"]["solver_optimal"] and rec["verification"].get("primal_feasible", False))
    rec["status"] = "numerically_accepted" if rec["accepted"] else "not_accepted"
    from smatrix_bootstrap.sdp.accuracy import trace_report
    rec["iterate_trace"] = trace_report(dest / "out/iterations.json", cfg.precision)
    rec["observables"] = rec["verification"].pop("observables")
    rec["subthreshold"] = rec["verification"].pop("subthreshold")
    rec["postprocess"] = {"note": a.note, "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                          "tool": "scripts/sdp/finish_leaf.py", "scope": "readback/verification only; no re-solve"}
    write_json(report, rec)
    print(dest.name, rec["status"], "terminate", rec["convergence"]["terminate_reason"],
          "primal_feasible", rec["verification"]["primal_feasible"], "objective", rec["verification"].get("objective_recomputed"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
