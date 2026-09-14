"""C6 (Fig.9, rho at the three UV representatives) and C7 (Fig.10, S0/S2) verdicts from accepted leaves.

Usage: python scripts/sdp/c67_eval.py --tip LEAF --ref LEAF --mid LEAF --label NAME --out FILE

Roles follow the frozen representative rule (PLAN 2.2): tip = +x end (paper red), ref = nearest upper-branch
point to the Weinberg point (paper light pink), mid = the neighbouring upper-branch point on the tip side
(paper pink).  Thresholds are the pre-registered ones in smatrix_bootstrap.sdp.claims.c6 / c7.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from smatrix_bootstrap.sdp import claims


def obs(path):
    r = json.loads(Path(path).read_text())
    if not r.get("accepted"):
        raise SystemExit(f"{path}: not an accepted leaf")
    return r["observables"], r["verification"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--tip", required=True); p.add_argument("--ref", required=True); p.add_argument("--mid")
    p.add_argument("--label", required=True); p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    points, coords = {}, {}
    for role, path in (("tip", a.tip), ("ref", a.ref), ("mid", a.mid)):
        if path:
            o, v = obs(path)
            points[role] = o; coords[role] = {"x": v["f00_3"], "y": v["f11_3"], "leaf": path}
    c6 = claims.c6(points)
    c7 = claims.c7(points)
    out = {"packaging": a.label, "points": coords, "C6": c6, "C7": c7,
           "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    Path(a.out).write_text(json.dumps(out, indent=1, default=float))
    print("C6:", c6["verdict"], "|", c6.get("evidence"))
    for r in c6.get("rows", []):
        print("   ", r)
    print("C7:", c7["verdict"], "|", c7.get("evidence"))
    for r in c7.get("rows", []):
        print("   ", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
