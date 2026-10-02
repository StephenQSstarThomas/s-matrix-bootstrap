"""Variant B of PREREGISTRATION amendments 5-6: the 2309 regulariser rule on the M_reg ladder.

Usage: python scripts/gtb/ladder_rule.py OUT_JSON RUNG_DIR... [--l24 MREG:LEAF ...]
Each RUNG_DIR is an ni = 46 pre-Watson F1max leaf (l-inf relax, R = M_reg * Mrho) with report.json/solution.npz.
eta = max over every node of |1 + i h| for l = 20..23 from the L = 24 export of the same kernels.  M_reg* = the
largest rung with eta <= 1.02 whose L = 24 re-solve (--l24 MREG:LEAF) gives F1 within 2 % of the L = 20 value.
"""
import json
import sys
from pathlib import Path

import numpy as np

from smatrix_bootstrap.gtb import kernels
from smatrix_bootstrap.gtb.leafanalysis import omitted_waves

ETA_MAX, L_STAB = 1.02, 0.02


def rung(leaf, big):
    rep = json.loads((Path(leaf) / "report.json").read_text())
    M = rep["spec"]["M"]
    v = np.load(Path(leaf) / "solution.npz")["v"]
    eta = omitted_waves(v[:1 + 2 * M + M * M + M * (M + 1) // 2], big, rep["spec"]["L"])["eta_max"]
    cert = Path(leaf) / "certificate.json"
    f1 = float(json.loads(cert.read_text())["lower"]["value"]) if cert.exists() else rep["observables"]["F1"]
    return {"leaf": str(leaf), "M_reg": rep["spec"]["reg_bound"], "F1": f1, "V1": cert.exists(),
            "eta_max": max(eta.values()), "eta": eta, "eta_pass": max(eta.values()) <= ETA_MAX}


if __name__ == "__main__":
    out, args = Path(sys.argv[1]), sys.argv[2:]
    l24 = {}
    if "--l24" in args:
        i = args.index("--l24")
        for item in args[i + 1:]:
            m, leaf = item.split(":", 1)
            l24[float(m)] = leaf
        args = args[:i]
    big = kernels.load("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M50_L24")
    rungs = sorted((rung(d, big) for d in args), key=lambda r: -r["M_reg"])
    star = None
    for r in rungs:
        if not r["eta_pass"]:
            continue
        leaf = l24.get(r["M_reg"])
        if leaf is None:
            r["L_stability"] = "missing"
            break                                  # the largest eta-passing rung needs its L = 24 check first
        f24 = json.loads((Path(leaf) / "report.json").read_text())["observables"]["F1"]
        r["F1_L24"], r["L_rel"] = f24, abs(f24 - r["F1"]) / abs(r["F1"])
        if r["L_rel"] <= L_STAB:
            star = r["M_reg"]
            break
    res = {"rule": f"largest M_reg with eta <= {ETA_MAX} and L=24 within {L_STAB:.0%}", "M_reg_star": star,
           "rungs": rungs}
    out.write_text(json.dumps(res, indent=1) + "\n")
    for r in rungs:
        print(f"M_reg {r['M_reg']:>7}: F1 {r['F1']:.7f}  eta_max {r['eta_max']:.4f}  pass {r['eta_pass']}"
              + (f"  L24 {r.get('L_rel')}" if "L_rel" in r or "L_stability" in r else ""))
    print("M_reg* =", star)
