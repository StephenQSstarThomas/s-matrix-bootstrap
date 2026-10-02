"""Pre-registered rules for the N-class claims C2, C3, C22-C25, C27 (PREREGISTRATION.md table and amendment 4).

Inputs are V1-certified objective values (C2, C3), infeasibility certificates (C3) and `gtb analyze` outputs of the
Watson leaves (C22-C27).  Every rule is checked against the paper's own numbers in the tests (amendment 3).
"""
from __future__ import annotations

import re

import numpy as np

from . import analysis as A

F1_TREE, F0_TREE = -0.0291487, 0.145743
C2_RMS, TIP_REL = 3e-4, 0.01


def fig3_dots(path) -> dict:
    """{'edges': {k: {'max', 'min'}} on the paper's grid F0 = k/50 * 0.1632, 'points': {colour: (F0, F1)}}."""
    edges, points = {}, {}
    for line in open(path):
        m = re.match(r"(\w+)\s+F0=([-\d.]+)\s+F1=([-\d.]+)", line)
        if not m:
            continue
        c, f0, f1 = m.group(1), float(m.group(2)), float(m.group(3))
        if c != "yellow":
            points[c] = (f0, f1)
            continue
        k = int(round(f0 / 0.1632 * 50))
        e = edges.setdefault(k, {"max": f1, "min": f1})
        e["max"], e["min"] = max(e["max"], f1), min(e["min"], f1)
    return {"edges": edges, "points": points}


def c2(ours: dict, tree: dict, fig3_path) -> dict:
    """ours: {k: {'max': F1, 'min': F1}} (None = infeasible); tree: {'max', 'min'} at F0 = 0.145743."""
    paper = fig3_dots(fig3_path)["edges"]
    d = [ours[k][s] - paper[k][s] for k in sorted(ours) for s in ("max", "min") if ours[k].get(s) is not None]
    rms = float(np.sqrt(np.mean(np.square(d)))) if d else float("inf")
    tip = all(ours.get(50, {}).get(s) is not None for s in ("max", "min"))      # 0.163108 is 0.06 % below 0.1632
    inside = tree.get("min") is not None and tree["min"] <= F1_TREE <= tree["max"]
    ok = rms <= C2_RMS and tip and inside
    return {"verdict": "PASS" if ok else "FAIL", "rms": rms, "n_dots": len(d), "tip_feasible": tip,
            "tree_inside": inside, "differences": {f"{k}:{s}": ours[k][s] - paper[k][s] for k in sorted(ours)
                                                   for s in ("max", "min") if ours[k].get(s) is not None}}


def c3(at2e3: dict, excl1e3: dict, excl2e3=None, at1e3=None) -> dict:
    """PASS: tree inside at eps 2e-3 (V1 max/min straddle it) and certified infeasible at 1e-3.
    FAIL: certified exclusion at 2e-3, or the tree inside at 1e-3 (threshold below 1e-3).  Else INCONCLUSIVE."""
    inside = lambda r: r is not None and r.get("min") is not None and r["min"] <= F1_TREE <= r["max"]
    if (excl2e3 or {}).get("certified_infeasible") or inside(at1e3):
        v = "FAIL"
    elif inside(at2e3) and excl1e3.get("certified_infeasible"):
        v = "PASS"
    else:
        v = "INCONCLUSIVE"
    return {"verdict": v, "inside_2e3": inside(at2e3), "excluded_1e3": bool(excl1e3.get("certified_infeasible"))}


def _unwrapped(delta_deg):
    return np.degrees(np.unwrap(np.radians(2 * np.asarray(delta_deg, dtype=float))) / 2)


def c22(base: dict, var: dict, waves=("S0", "P1", "S2"), emax=1.4, tol=10.0) -> dict:
    """Per-wave RMS of (delta_var - delta_base) modulo 180 deg at the base nodes below emax GeV."""
    Eb, Ev = np.asarray(base["E_GeV"]), np.asarray(var["E_GeV"])
    sel = Eb <= emax
    rms = {}
    for w in waves:
        dv = np.interp(Eb[sel], Ev, _unwrapped(var["phases"][w]["delta_deg"]))
        d = (dv - np.asarray(base["phases"][w]["delta_deg"])[sel] + 90) % 180 - 90
        rms[w] = float(np.sqrt(np.mean(d ** 2)))
    return {"verdict": "PASS" if all(r <= tol for r in rms.values()) else "FAIL", "rms_deg": rms}


def peak(an: dict, ch: str):
    """(E_GeV, |F|) of the largest |F| among the nodes below s0."""
    E = np.asarray(an["E_GeV"])
    f = np.asarray(an["FF_abs"][ch])[:E.size]
    i = int(np.argmax(f))
    return float(E[i]), float(f[i])


def c24(analyses: list, lo=50.0, hi=110.0) -> dict:
    """Spread of the |F^1_1| peak node over r11 = 1.5, 2, 2.5 in 80 +- 30 MeV."""
    e = [peak(a, "P1")[0] * 1000 for a in analyses]
    spread = max(e) - min(e)
    return {"verdict": "PASS" if lo <= spread <= hi else "FAIL", "peaks_MeV": e, "spread_MeV": spread}


def c25(base: dict, variants: list, de=10.0, dh=0.05) -> dict:
    """n_SR = 4, 5 against 3: |F| peak positions within 10 MeV and heights within 5 % for S0, P1, D0."""
    rows, ok = [], True
    for v in variants:
        for ch in ("S0", "P1", "D0"):
            (eb, hb), (ev, hv) = peak(base, ch), peak(v, ch)
            good = abs(ev - eb) * 1000 <= de and abs(hv - hb) <= dh * hb
            ok &= good
            rows.append({"ch": ch, "dE_MeV": (ev - eb) * 1000, "rel_dh": (hv - hb) / hb, "ok": good})
    return {"verdict": "PASS" if ok else "FAIL", "rows": rows}


def c27(no_p1: dict, no_d0: dict) -> dict:
    """Without the P1 FESR no delta^1_1 90-deg crossing below 1.2 GeV; without the D0 FESR max delta^0_2 >= 45 deg
    below 1.83 GeV."""
    E1 = np.asarray(no_p1["E_GeV"])
    x = A.crossing(E1, _unwrapped(no_p1["phases"]["P1"]["delta_deg"]), 90.0)
    rho_gone = x is None or x > 1.2
    E2 = np.asarray(no_d0["E_GeV"])
    d0 = _unwrapped(no_d0["phases"]["D0"]["delta_deg"])[E2 <= 1.83]
    f2_kept = bool(d0.size and d0.max() >= 45.0)
    return {"verdict": "PASS" if rho_gone and f2_kept else "FAIL", "p1_crossing_GeV": x,
            "max_d0_deg": float(d0.max()) if d0.size else None}


def red_2309(ref_dir) -> dict:
    """{wave: (E_GeV, delta_deg)} of the red series of 2309.12402v3 Figs. 9 (P1) and 10a/b (S0, S2)."""
    import csv
    from pathlib import Path
    files = {"P1": "figure9_p1_phases.csv", "S0": "figure10_s0_phases.csv", "S2": "figure10_s2_phases.csv"}
    out = {}
    for w, name in files.items():
        rows = [r for r in csv.DictReader(open(Path(ref_dir) / name)) if r["group"] == "red"]
        rows.sort(key=lambda r: float(r["energy_gev"]))
        out[w] = (np.array([float(r["energy_gev"]) for r in rows]), np.array([float(r["phase_deg"]) for r in rows]))
    return out


def c26(green: dict, ref: dict, emax=1.2, tol=10.0) -> dict:
    """Our green Watson phases against the 2309 red curves: RMS (mod 180 deg) at the 2309 energies below 1.2 GeV."""
    rms = {}
    for w, (E, d) in ref.items():
        base = {"E_GeV": list(E), "phases": {w: {"delta_deg": list(d)}}}
        rms[w] = c22(base, green, waves=(w,), emax=emax, tol=tol)["rms_deg"][w]
    return {"verdict": "PASS" if all(r <= tol for r in rms.values()) else "FAIL", "rms_deg": rms}
