"""Run the linearised unitarity-saturation (Watson) iteration from an accepted UV leaf.

Usage: python scripts/sdp/watson_iterate.py --start LEAF/report.json --out ROOT/watson_NAME [--rounds 5] [--waves S0+P1+S2]
       [--keep-section] [--timeout 7200] [--tol 0.02]

Round r solves the same finite problem (every constraint kept, section released unless --keep-section) with the
objective of watson.py built from round r-1 (round 0 = the start leaf).  Each round is a full SDPB solve through the
usual `support` path, so it inherits the acceptance rules (solver optimal, Arb feasibility).  The chain stops when a
round is not accepted, when max_k |h_k - t_k| over the listed waves falls below --tol (the fixed point is reached
within the tolerance), or after --rounds rounds.  Every round's metrics are appended to WATSON_ITER.json in --out.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from smatrix_bootstrap.sdp import watson  # noqa: E402

PY = sys.executable
REPO = Path(__file__).resolve().parents[2]


def distance_to_targets(report_path, waves, saturated_only=False):
    """max_k |h_k - t_k| where t is the target that THIS leaf would define (how far from its own fixed point).

    With ``saturated_only`` the maximum runs over the nodes where the current is two-pion saturated
    (|F|^2 / rho_hat >= 0.5) or, for waves without a current, over all nodes: at nodes where rho_hat >> |F|^2 the
    form-factor phase is not constrained by anything physical and the target built from it is noise (seen at the
    threshold node of P1, where |F| = 18 with |F|^2/rho_hat = 0 after the first round).
    """
    t = watson.targets_from_leaf(report_path, waves)
    r = json.loads(Path(report_path).read_text())
    sol = np.load(Path(report_path).parent / "solution.npz")
    M = r["spec"]["M"]; nodes = np.array(t["nodes"])
    worst = 0.0
    for w in waves:
        S = watson.node_S(r, w); h = -1j * (S - 1.0)
        tt = np.array([complex(a, b) for a, b in t["targets"][w]])
        keep = np.ones(len(nodes), dtype=bool)
        if saturated_only and w in ("S0", "P1") and r["spec"].get("uv"):
            ell = watson.WAVE_INDEX[w][1]
            F = watson.form_factor(sol, ell, M); rh = np.asarray(sol["rho_hat"][ell], dtype=float)
            frac = np.abs(F[nodes]) ** 2 / np.where(rh[nodes] > 0, rh[nodes], np.inf)
            keep = frac >= 0.5
        if keep.any():
            worst = max(worst, float(np.max(np.abs(h[nodes][keep] - tt[keep]))))
    return worst


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--start", required=True); p.add_argument("--out", required=True)
    p.add_argument("--rounds", type=int, default=5); p.add_argument("--waves", default="S0+P1+S2")
    p.add_argument("--keep-section", action="store_true"); p.add_argument("--timeout", type=int, default=7200)
    p.add_argument("--tol", type=float, default=0.02)
    p.add_argument("--pin", type=float, default=None, help="pinned-point control: keep the start leaf's section and add the slab "
                   "d.(f00,f11) >= v* - PIN; constraints always from the start leaf, targets from the previous round")
    a = p.parse_args(argv)
    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    waves = tuple(a.waves.split("+"))
    log_path = out / "WATSON_ITER.json"
    log = {"start": str(Path(a.start).resolve()), "waves": list(waves), "keep_section": a.keep_section, "tol": a.tol, "pin": a.pin,
           "rounds": [], "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    src = Path(a.start).resolve()
    m0 = watson.iteration_metrics(src, waves); m0["distance_to_own_targets"] = distance_to_targets(src, waves)
    m0["distance_to_own_targets_saturated_nodes"] = distance_to_targets(src, waves, saturated_only=True)
    log["rounds"].append({"round": 0, "leaf": str(src), "accepted": True, "metrics": m0})
    log_path.write_text(json.dumps(log, indent=1, default=float))
    for r in range(1, a.rounds + 1):
        dest = out / f"round_{r:02d}"
        base = Path(a.start).resolve() if a.pin is not None else src      # pinned: constraints always from the start leaf
        cmd = [PY, "-m", "smatrix_bootstrap.run", "sdp", "support", "--source-report", str(base), "--out", str(dest),
               "--functional", "watson", a.waves, "0", "max", "--timeout", str(a.timeout)]
        if a.pin is not None:
            cmd += ["--face-margin", str(a.pin), "--watson-targets", str(src)]
        if a.keep_section:
            cmd.append("--watson-keep-section")
        t0 = time.time()
        proc = subprocess.run(cmd, cwd=str(REPO), env={**__import__("os").environ, "PYTHONPATH": "src",
                              "OPENBLAS_NUM_THREADS": "4", "OMP_NUM_THREADS": "4"}, capture_output=True, text=True)
        rep_path = dest / "report.json"
        rep = json.loads(rep_path.read_text()) if rep_path.exists() else {}
        entry = {"round": r, "leaf": str(rep_path), "accepted": bool(rep.get("accepted")), "status": rep.get("status"),
                 "seconds": time.time() - t0, "returncode": proc.returncode, "stderr_tail": proc.stderr[-600:]}
        if rep.get("accepted"):
            m = watson.iteration_metrics(rep_path, waves); m["distance_to_own_targets"] = distance_to_targets(rep_path, waves)
            m["distance_to_own_targets_saturated_nodes"] = distance_to_targets(rep_path, waves, saturated_only=True)
            m["objective_value"] = rep["verification"].get("functional", {}).get("value")
            entry["metrics"] = m
        log["rounds"].append(entry); log_path.write_text(json.dumps(log, indent=1, default=float))
        print(json.dumps({k: v for k, v in entry.items() if k != "stderr_tail"}, default=float)[:600], flush=True)
        if not rep.get("accepted"):
            log["stopped"] = f"round {r} not accepted"; break
        src = rep_path
        if entry["metrics"]["distance_to_own_targets_saturated_nodes"] < a.tol:
            log["stopped"] = (f"converged at round {r} (max |h - t| over two-pion-saturated nodes = "
                              f"{entry['metrics']['distance_to_own_targets_saturated_nodes']:.4f} < {a.tol})"); break
    log["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    log_path.write_text(json.dumps(log, indent=1, default=float))
    print("done:", log.get("stopped", "rounds exhausted"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
