"""Degenerate-face diagnostic verdicts (pre-registered rules F1-F5, GATE_LOG 2026-09-14 02:25Z).

Usage: python scripts/sdp/face_eval.py --root RESULTS_ROOT --out FILE

Reads every ``face_*/report.json`` under ROOT (leaves written by ``sdp support --face-margin --functional``),
groups them by source leaf, and reports for each registered node functional the accepted extreme values on the
near-optimal face, the width of the range, and the F1/F2/F3 verdicts.  Values come from ``verification.functional``
(Arb re-evaluation of the node S-matrix), never from the solver objective.  A leaf that is not accepted is listed
as such and does not enter a verdict.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import time
from pathlib import Path

BW = {"M_MeV": 822.0, "Gamma_MeV": 150.0}          # Fig.9 P1 reference used in the pre-registration
NODE_E = {37: 0.7317, 38: 0.7921, 39: 0.8644}      # GeV, M=50 ladder


def bw_reference(node):
    """1 - Re S and Im S of an elastic Breit-Wigner at the node energy."""
    e = 1000 * NODE_E[node]
    d = math.atan2(BW["Gamma_MeV"] / 2, BW["M_MeV"] - e)
    return {"delta_deg": math.degrees(d), "ImKH": 1 - math.cos(2 * d), "ImS": math.sin(2 * d)}


def delta_at_eta1(v):
    """delta from 1 - Re S = v for an elastic wave (eta = 1), in degrees; None outside [0, 2]."""
    return None if not 0 <= v <= 2 else math.degrees(math.acos(1 - v) / 2)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    leaves = {}
    for path in sorted(glob.glob(str(Path(a.root) / "face_*" / "report.json"))):
        r = json.loads(Path(path).read_text())
        f = r.get("functional")
        if f is None:
            continue
        src = r.get("reused_source_report", "?")
        v = r.get("verification", {})
        rec = {"leaf": path, "status": r.get("status"), "accepted": bool(r.get("accepted")),
               "value": (v.get("functional") or {}).get("value"), "f00_3": v.get("f00_3"), "f11_3": v.get("f11_3"),
               "face": v.get("face"), "seconds": r.get("seconds"),
               "terminate": r.get("convergence", {}).get("terminate_reason")}
        key = (f["kind"], f["wave"], f["node"])
        leaves.setdefault(src, {}).setdefault("%s_%s_%d" % key, {})[f["sense"]] = rec
    out = {"rules": "GATE_LOG 2026-09-14 02:25Z F1-F5", "bw_reference": {n: bw_reference(n) for n in NODE_E},
           "sources": {}, "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    for src, funcs in leaves.items():
        s = json.loads(Path(src).read_text()) if Path(src).exists() else {}
        sv = s.get("verification", {})
        entry = {"source_point": {"f00_3": sv.get("f00_3"), "f11_3": sv.get("f11_3"),
                                  "direction": s.get("direction"), "fix_f00": s.get("fix_f00")},
                 "source_observables": {w: {"delta_deg_at_38": (s.get("observables", {}).get(w, {}).get("delta_deg") or [None] * 39)[38],
                                            "eta_at_38": (s.get("observables", {}).get(w, {}).get("eta") or [None] * 39)[38],
                                            "crossing_90_GeV": s.get("observables", {}).get(w, {}).get("crossing_90_GeV")}
                                        for w in ("S0", "P1")},
                 "functionals": {}, "verdicts": {}}
        for name, senses in funcs.items():
            kind, wave, node = name.split("_"); node = int(node)
            mx = senses.get("max", {}); mn = senses.get("min", {})
            hi = mx.get("value") if mx.get("accepted") else None
            lo = mn.get("value") if mn.get("accepted") else None
            row = {"max": mx or None, "min": mn or None, "range_hi": hi, "range_lo": lo,
                   "width": (hi - lo) if (hi is not None and lo is not None) else None}
            if kind == "ImKH":
                row["delta_deg_at_eta1"] = {"hi": delta_at_eta1(hi) if hi is not None else None,
                                            "lo": delta_at_eta1(lo) if lo is not None else None}
                row["bw_reference"] = bw_reference(node) if node in NODE_E else None
            entry["functionals"][name] = row
            if kind == "ImKH" and wave == "P1" and node == 38:
                if hi is not None:
                    entry["verdicts"]["F1_exclusion"] = {"max_ImKH_P1_38": hi, "threshold": 1.5,
                        "paper_shape_excluded_from_face": hi < 1.5,
                        "reading": "Re S_P1(0.792 GeV) > %.3f on the whole near-optimal face" % (1 - hi)}
                if row["width"] is not None:
                    w = row["width"]
                    entry["verdicts"]["F2_width"] = {"width": w, "verdict":
                        "undetermined (>=0.5)" if w >= 0.5 else "determined (<=0.05)" if w <= 0.05 else "partially determined"}
            if kind == "ImKH" and wave == "P1" and node == 39 and hi is not None:
                entry["verdicts"]["F3_post_resonance"] = {"max_ImKH_P1_39": hi, "threshold": 1.4, "necessary_condition_met": hi >= 1.4}
            if kind == "ImKH" and wave == "S0" and node == 38 and (hi is not None or lo is not None):
                entry["verdicts"]["F4_S0_range_deg_at_eta1"] = row["delta_deg_at_eta1"]
            if kind == "ImF" and wave == "P1" and node == 38:
                entry["verdicts"]["F5_ImF1_38_range"] = [lo, hi]
        out["sources"][src] = entry
    Path(a.out).write_text(json.dumps(out, indent=1, default=float))
    for src, e in out["sources"].items():
        print(Path(src).parent.parent.name if Path(src).parent.name == "tip" else Path(src).parent.name, e["source_point"])
        for name, row in e["functionals"].items():
            print("   %-14s lo %-12s hi %-12s width %s" % (name, row["range_lo"], row["range_hi"], row["width"]))
        for k, v in e["verdicts"].items():
            print("   ", k, v)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
