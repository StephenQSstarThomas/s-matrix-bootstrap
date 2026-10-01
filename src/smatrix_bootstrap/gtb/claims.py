"""Pre-registered claim verdicts (docs/gtb2403/PREREGISTRATION.md, amendment 1)."""
from __future__ import annotations

C1_PAPER = {44: -0.027229, 46: -0.028770, 48: -0.030397}       # Fig. 3 vector read-off (+-3e-5)
C1_TOL = 1e-4


def c1_verdict(lower: float, upper: float, paper: float, tol: float = C1_TOL) -> str:
    """Certified bracket [lower, upper] of F1max against the band paper +- tol: PASS only if inside, FAIL only if
    disjoint, otherwise INCONCLUSIVE (a verdict needs a bracket that decides the inequality)."""
    if lower > upper:
        raise ValueError("inverted bracket")
    lo, hi = paper - tol, paper + tol
    if lo <= lower and upper <= hi:
        return "PASS"
    if upper < lo or lower > hi:
        return "FAIL"
    return "INCONCLUSIVE"


# ------------------------------------------------------------------ phase and form-factor claims (C4-C7, C12)
import csv
import json
from pathlib import Path

import numpy as np

from .analysis import crossing

PAPER_DATA = Path("/playpen1/shiqiu/gtb2403_runs/paper_data")
COLOUR = {44: "GTB_blue_ni44", 46: "GTB_green_ni46", 48: "GTB_red_ni48"}
WAVES6 = ("S0", "D0", "S2", "D2", "P1", "F1")
N0 = 39


def _panel(stem, series):
    """{node_j (1-based): y} for the dots of one series in one extracted panel."""
    out = {}
    with open(PAPER_DATA / f"{stem}.csv") as fh:
        for r in csv.DictReader(fh):
            if r["series"] == series and r["kind"] == "dot" and r.get("node_j"):
                out[int(r["node_j"])] = float(r["y"])
    return out


def _nodes():
    return [n["sqrt_s_gev"] for n in json.loads((PAPER_DATA / "index.json").read_text())["nodes"]][:N0]


def paper_as_analysis(ni):
    """The paper's own curves in the leafanalysis format (missing dots -> NaN); used to test the rules."""
    E, col = _nodes(), COLOUR[ni]
    an = {"n0": N0, "E_GeV": E, "phases": {}, "FF_abs": {}}
    for w in WAVES6:
        d, e = _panel(f"fig04_{w}_delta", col), _panel(f"fig04_{w}_eta", col)
        an["phases"][w] = {"delta_deg": [d.get(j, np.nan) for j in range(1, N0 + 1)],
                           "eta": [e.get(j, np.nan) for j in range(1, N0 + 1)]}
    f = _panel("fig06_S0ff", col)
    an["FF_abs"]["S0"] = [f.get(j, np.nan) for j in range(1, N0 + 1)]
    return an


def _wrap(d):
    return (np.asarray(d) + 90.0) % 180.0 - 90.0


def _cross(an, w):
    d, E = np.array(an["phases"][w]["delta_deg"], float), np.array(an["E_GeV"], float)
    ok = ~np.isnan(d)
    return crossing(E[ok], d[ok], 90.0)


def phase_claims(analyses: dict) -> dict:
    """Pre-registered C4, C5, C6, C7, C12 for {ni: leafanalysis result} (post-Watson leaves)."""
    E = np.array(_nodes())
    rows, out = {}, {}
    c4 = c5 = c6 = c7 = c12 = True
    for ni, an in analyses.items():
        ref = paper_as_analysis(ni)
        r = rows[ni] = {}
        for w in WAVES6:
            ours, pap = np.array(an["phases"][w]["delta_deg"][:N0], float), np.array(ref["phases"][w]["delta_deg"])
            m = (E < 1.4) & ~np.isnan(pap) & ~np.isnan(ours)
            rms = float(np.sqrt(np.mean(_wrap(ours[m] - pap[m]) ** 2))) if m.any() else None
            eta = float(np.nanmin(an["phases"][w]["eta"][:N0]))
            r[w] = {"rms_deg": rms, "min_eta": eta}
            c4 &= rms is not None and rms <= 5.0 and eta >= 0.97
        for w, tol, key in (("P1", 0.020, "c5"), ("D0", 0.030, "c6")):
            a, b = _cross(an, w), _cross(ref, w)
            r[f"{w}_cross_GeV"] = (a, b)
            ok = a is not None and b is not None and abs(a - b) <= tol
            if key == "c5":
                c5 &= ok
            else:
                c6 &= ok
        last = float(_wrap(an["phases"]["P1"]["delta_deg"][N0 - 1] - ref["phases"]["P1"]["delta_deg"][N0 - 1]))
        r["P1_last_node_diff_deg"] = last
        c7 &= abs(last) <= 20.0
        fo, fp = np.array(an["FF_abs"]["S0"][:N0], float), np.array(ref["FF_abs"]["S0"], float)
        m = ~np.isnan(fp)
        jo, jp = int(np.nanargmax(np.where(m, fo, np.nan))), int(np.nanargmax(fp))
        r["S0FF_peak"] = {"ours_GeV": float(E[jo]), "paper_GeV": float(E[jp]), "ours": float(fo[jo]), "paper": float(fp[jp])}
        c12 &= abs(E[jo] - E[jp]) <= 0.040 and abs(fo[jo] / fp[jp] - 1) <= 0.10
    for cid, ok in (("C4", c4), ("C5", c5), ("C6", c6), ("C7", c7), ("C12", c12)):
        out[cid] = {"verdict": "PASS" if ok else "FAIL"}
    out["rows"] = rows
    return out
