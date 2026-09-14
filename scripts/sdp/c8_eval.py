"""C8 (Fig.11) verdict: rho position stability across (M,L) at the UV tip.

Usage: python scripts/sdp/c8_eval.py --root RESULTS_ROOT --out FILE
Pre-registered claims.c8 rule: L spread <= 20 MeV at fixed M=50; M spread 40-70 MeV at fixed L=10 with the
paper's ordering M60 < M45 < M50; S0 at 1 GeV within 85-105 deg.  Only the 'tip' representative is available for
every configuration; the rule is applied to it.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from smatrix_bootstrap.sdp import claims

RUNS = {(50, 10): "uv_SRa_tip/tip", (50, 8): "uv_M50_L8_tip/tip", (50, 12): "uv_M50_L12_tip/tip", (45, 10): "uv_M45_L10_tip/tip", (60, 10): "uv_M60_L10_tip/tip"}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", required=True); p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    by_ml, rows = {}, []
    for ml, rel in RUNS.items():
        path = Path(a.root) / rel / "report.json"
        alt = Path(a.root) / (rel.split("/")[0] + "_resume") / "report.json"      # checkpoint resume of a timed-out leaf
        if alt.exists() and json.loads(alt.read_text()).get("accepted"):
            path, rel = alt, alt.parent.name
        if not path.exists(): continue
        r = json.loads(path.read_text())
        if not r.get("accepted"): continue
        by_ml[ml] = {"tip": r["observables"]}
        o = r["observables"]
        rows.append({"M": ml[0], "L": ml[1], "x_tip": r["verification"]["f00_3"], "P1_crossing_MeV": None if o["P1"]["crossing_90_GeV"] is None else 1000 * o["P1"]["crossing_90_GeV"],
                     "min_eta_P1": o["P1"]["min_eta_below_1p2GeV"], "S0_at_1GeV_deg": claims._at(o["S0"], 1.0), "leaf": rel})
    res = claims.c8(by_ml)
    if isinstance(res.get("rho_MeV"), dict):
        res["rho_MeV"] = {f"M{k[0]}_L{k[1]}": v for k, v in res["rho_MeV"].items()}
    res.update(rows=rows, recorded_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), claim="C8 Fig.11 M/L stability of the UV tip")
    Path(a.out).write_text(json.dumps(res, indent=1, default=float))
    print("C8:", res["verdict"], "|", res.get("evidence"))
    for row in rows: print("   ", row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
