"""C1 (Fig.3) verdict: 24 accepted support directions of the pure-unitarity region versus the digitised figure.

Usage: python scripts/sdp/c1_eval.py --root RESULTS_ROOT --out FILE
Uses the pre-registered claims.c1 rule (extrema within 2 % of the digitised Fig.3, 24 verified directions).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from smatrix_bootstrap.sdp import claims, figures


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--root", required=True); p.add_argument("--out", required=True)
    p.add_argument("--source", default="fig3_pure_linf_1e3_tip_v2", help="tip run whose basis the directions reuse")
    a = p.parse_args(argv)
    recs = []
    for r in figures.load_reports(a.root):
        f = r["_file"]
        if f.startswith("fig3_dir") and f.count("/") == 1:
            r["job"] = "dir" + f[len("fig3_dir"):f.index("/")]; recs.append(r)
        elif f == f"{a.source}/tip/report.json":
            r["job"] = "dir_tip"; recs.append(r)
    res = claims.c1(recs)
    res.update(n_records=len(recs), recorded_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), claim="C1 Fig.3 pure-unitarity region")
    Path(a.out).write_text(json.dumps(res, indent=1, default=float))
    print("C1:", res["verdict"], "|", res.get("evidence"))
    for row in res.get("rows", []):
        print("   ", row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
