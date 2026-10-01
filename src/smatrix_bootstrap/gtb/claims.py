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


# ------------------------------------------------------------------ table claims (C8-C10, C13-C19), blue = 44, red = 48
TABLES = {   # v4 Tables 1-7 and inline numbers: (blue, red)
    "a00": (0.174, 0.191), "a20": (-0.0357, -0.0403), "b00": (0.285, 0.282), "b20": (-0.061, -0.071),
    "a11": (28.7e-3, 30.7e-3), "b11": (12.4e-3, 7.28e-3), "a02": (14.6e-4, 14.9e-4), "a22": (3.41e-4, 2.53e-4),
    "r2S": (0.64, 0.55), "r2V": (0.394, 0.391), "r2T": (0.153, 0.140),
    "l1": (0.96, -0.10), "l2": (4.5, 4.8), "l4": (4.8, 4.2), "l6": (14.6, 14.5),
    "pole_re": (829, 844), "pole_im": (71, 82), "m_rho": (838, 859), "G_rho": (114, 128), "g_rho": (4.9, 5.2),
    "4pi_grho2": (0.61, 0.57), "m1": (824, 841), "G1": (137, 164), "m2": (1635, 1645), "G2": (248, 272),
    "A1": (1.29, 1.39), "A2": (0.27, 0.28), "m_f2": (1366, 1388), "G_f2": (164, 188), "g_f2": (8.5, 8.9)}


def _table_tol(key):
    b, r = TABLES[key]
    return lambda colour: max(0.05 * abs(TABLES[key][colour]), 0.5 * abs(b - r))


def _abs_tol(x):
    return lambda colour: x


def _rel_tol(f):
    return lambda key: (lambda colour: f * abs(TABLES[key][colour]))


def _get(an, key):
    th, ra, lb, fi = an.get("threshold", {}), an.get("radii_fm2", {}), an.get("lbar", {}), an.get("fits", {})
    pole = (an.get("pole") or {}).get("interp")
    m = {**th, "r2S": ra.get("S"), "r2V": ra.get("V"), "r2T": ra.get("T"), **lb,
         "pole_re": pole["sqrt_s_MeV"][0] if pole else None, "pole_im": -pole["sqrt_s_MeV"][1] if pole else None}
    gs, tb, f2 = fi.get("rho_gs", {}), fi.get("two_bw", {}), fi.get("f2_bw", {})
    m.update({"m_rho": gs.get("m_MeV"), "G_rho": gs.get("Gamma_MeV"), "g_rho": gs.get("g_rho_pi_pi"),
              "4pi_grho2": fi.get("four_pi_over_g_rho2"), "m1": tb.get("m1_MeV"), "G1": tb.get("G1_MeV"),
              "m2": tb.get("m2_MeV"), "G2": tb.get("G2_MeV"), "A1": tb.get("A1"), "A2": tb.get("A2"),
              "m_f2": f2.get("m_MeV"), "G_f2": f2.get("Gamma_MeV"), "g_f2": f2.get("g_f2_pi_pi")})
    return m.get(key)


CLAIM_KEYS = {
    "C8": {k: _table_tol(k) for k in ("a00", "a20", "b00", "b20", "a11", "b11", "a02", "a22")},
    "C9": {k: _table_tol(k) for k in ("r2S", "r2V", "r2T")},
    "C10": {k: _table_tol(k) for k in ("l1", "l2", "l4", "l6")},
    "C13": {"pole_re": _abs_tol(20), "pole_im": _abs_tol(15)},
    "C14": {"m_rho": _abs_tol(20), "G_rho": _rel_tol(0.20)("G_rho")},
    "C15": {"g_rho": _rel_tol(0.05)("g_rho")},
    "C16": {"4pi_grho2": _rel_tol(0.10)("4pi_grho2")},
    "C17": {"m1": _abs_tol(20), "m2": _abs_tol(40), "G1": _rel_tol(0.25)("G1"), "G2": _rel_tol(0.25)("G2"),
            "A1": _rel_tol(0.15)("A1"), "A2": _rel_tol(0.15)("A2")},
    "C18": {"m_f2": _abs_tol(25), "G_f2": _rel_tol(0.25)("G_f2")},
    "C19": {"g_f2": _rel_tol(0.05)("g_f2")},
}


def table_claims(blue: dict, red: dict) -> dict:
    out = {}
    for cid, keys in CLAIM_KEYS.items():
        rows, ok = {}, True
        for key, tol in keys.items():
            for colour, an in ((0, blue), (1, red)):
                v = _get(an, key)
                good = v is not None and abs(v - TABLES[key][colour]) <= tol(colour)
                rows[f"{key}:{'blue' if colour == 0 else 'red'}"] = {"ours": v, "paper": TABLES[key][colour],
                                                                    "tol": tol(colour), "pass": good}
                ok &= good
        out[cid] = {"verdict": "PASS" if ok else "FAIL", "rows": rows}
    return out
