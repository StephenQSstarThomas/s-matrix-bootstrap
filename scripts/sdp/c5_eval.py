"""C5 (Fig.8) verdict for one UV packaging: +x end and the asymmetric shrink of the x_ref section.

Usage: python scripts/sdp/c5_eval.py --chi-hi LEAF --chi-lo LEAF --uv-tip LEAF --uv-hi LEAF --uv-lo LEAF --label NAME --out FILE

Pre-registered thresholds (PLAN 1 / claims.c5): UV +x end within 5 % of the digitised cyan tip 0.0811249;
upper shrink = y_hi(chiral) - y_hi(UV) within 25 % of 2.324e-4; lower rise = y_lo(UV) - y_lo(chiral) in [0, 1e-4];
ratio upper/lower >= 4.  Paper reference numbers are read from references/figure8_boundary.csv through
smatrix_bootstrap.sdp.claims.fig8_reference and the digitised x_ref values recorded in GATE_LOG.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

PAPER = {"uv_x_end": 0.0811249, "upper_shrink": 2.324e-4, "lower_rise": 3.687e-5, "chi_hi": -0.0043678,
         "chi_lo": -0.0051276, "uv_hi": -0.0046002, "uv_lo": -0.0050907}


def leaf(path):
    r = json.loads(Path(path).read_text())
    if not r.get("accepted"):
        raise SystemExit(f"{path}: not an accepted leaf")
    v = r["verification"]
    return r, v["f00_3"], v["f11_3"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("chi-hi", "chi-lo", "uv-tip", "uv-hi", "uv-lo"):
        p.add_argument("--" + k, required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    _, _, chi_hi = leaf(a.chi_hi); _, _, chi_lo = leaf(a.chi_lo)
    tip, x_end, _ = leaf(a.uv_tip); _, _, uv_hi = leaf(a.uv_hi); _, _, uv_lo = leaf(a.uv_lo)
    upper, lower = chi_hi - uv_hi, uv_lo - chi_lo
    ratio = upper / lower if lower > 0 else float("inf")
    rows = [
        {"quantity": "UV +x end", "ours": x_end, "paper": PAPER["uv_x_end"], "rel": x_end / PAPER["uv_x_end"] - 1,
         "pass": abs(x_end / PAPER["uv_x_end"] - 1) <= 0.05},
        {"quantity": "upper shrink at x_ref", "ours": upper, "paper": PAPER["upper_shrink"],
         "rel": upper / PAPER["upper_shrink"] - 1, "pass": abs(upper / PAPER["upper_shrink"] - 1) <= 0.25},
        {"quantity": "lower rise at x_ref", "ours": lower, "paper": PAPER["lower_rise"], "pass": 0 <= lower <= 1e-4},
        {"quantity": "upper/lower ratio", "ours": ratio, "paper": PAPER["upper_shrink"] / PAPER["lower_rise"], "pass": ratio >= 4},
    ]
    asym = all(r["pass"] for r in rows[1:])
    out = {"claim": "C5 Fig.8", "packaging": a.label, "verdict": "PASS" if all(r["pass"] for r in rows) else
           ("PARTIAL (asymmetry ok, +x end off)" if asym else "FAIL"),
           "rows": rows, "sections": {"chi_hi": chi_hi, "chi_lo": chi_lo, "uv_hi": uv_hi, "uv_lo": uv_lo},
           "paper_sections": {k: PAPER[k] for k in ("chi_hi", "chi_lo", "uv_hi", "uv_lo")},
           "uv_tip_phase_summary": {w: {"crossing_90_GeV": tip["observables"][w].get("crossing_90_GeV"),
                                        "min_eta_below_1p2": tip["observables"][w].get("min_eta_below_1p2GeV")}
                                    for w in ("S0", "S2", "P1")},
           "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    Path(a.out).write_text(json.dumps(out, indent=1))
    for r in rows:
        print(f"{r['quantity']:22s} ours {r['ours']:.6g}  paper {r['paper']:.6g}  {'ok' if r['pass'] else 'FAIL'}")
    print("verdict:", out["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
