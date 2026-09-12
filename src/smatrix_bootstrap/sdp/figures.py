"""Figures and the C1-C8 comparison tables, built from the report.json files.

The digitised originals in ``references/figure*.csv`` are the comparison
target.  Their metadata records that the abscissa was calibrated with
m_pi = 139.57 MeV whereas the paper's text sets m_pi = 140 MeV; the resulting
0.3% shift of the energy axis is carried explicitly (``ENERGY_AXIS_NOTE``) and
never absorbed into a fit.
"""
from __future__ import annotations

import csv
import glob
import json
import os

import numpy as np

REF = os.path.join(os.path.dirname(__file__), "..", "..", "..", "references")
ENERGY_AXIS_NOTE = ("digitised figures use m_pi = 139.57 MeV, the paper's text "
                    "m_pi = 140 MeV; 0.31% shift on the energy axis")


def load_csv(name: str) -> dict:
    path = os.path.join(REF, name)
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]
            if _isnum(rows[0][k])}
    for k in rows[0]:
        if not _isnum(rows[0][k]):
            cols[k] = np.array([r[k] for r in rows])
    return cols


def _isnum(x):
    try:
        float(x)
        return True
    except Exception:
        return False


def load_reports(root: str) -> list[dict]:
    out = []
    for f in sorted(glob.glob(os.path.join(root, "**", "report.json"), recursive=True)):
        with open(f) as fh:
            doc = json.load(fh)
        for r in doc.get("records", []):
            r["_file"] = os.path.relpath(f, root)
            out.append(r)
    return out


def boundary_points(records, want=None, verified_only: bool = True) -> np.ndarray:
    """Points in the projection plane.

    Only points that pass the a posteriori feasibility check on the unmodified
    operators are kept by default: the solver's own "optimal" status is not
    trusted on this problem (see the report), so a returned point counts as
    evidence only if it is verifiably inside the allowed region -- which makes
    every number an inner bound on the support function.
    """
    pts = []
    for r in records:
        if want and not r["job"].startswith(want):
            continue
        res = r.get("result", {})
        if res.get("f00_3") is None:
            continue
        if verified_only:
            u = r.get("verification", {}).get("unitarity", {})
            if not u.get("feasible", False):
                continue
        pts.append((res["f00_3"], res["f11_3"]))
    return np.array(pts) if pts else np.zeros((0, 2))


def best_support(records) -> dict:
    """Best verified-feasible objective per direction, over all configurations."""
    out = {}
    for r in records:
        res = r.get("result", {})
        u = r.get("verification", {}).get("unitarity", {})
        if res.get("objective") is None or not u.get("feasible", False):
            continue
        d = tuple(np.round(res.get("direction", [np.nan, np.nan]), 9))
        if d not in out or res["objective"] > out[d]["objective"]:
            out[d] = {"objective": res["objective"], "f00_3": res["f00_3"],
                      "f11_3": res["f11_3"], "file": r.get("_file"), "job": r["job"]}
    return out


def support_values(records) -> dict:
    """max d.(f00,f11) per direction, comparable with the Fig.3 metadata."""
    out = {}
    for r in records:
        res = r.get("result", {})
        if res.get("objective") is None:
            continue
        d = res.get("direction")
        out[tuple(np.round(d, 9))] = res["objective"]
    return out


def c1_table(records) -> dict:
    """C1: extrema of the pure-unitarity region against the digitised Fig. 3."""
    pts = boundary_points(records)
    ref = load_csv("figure3_boundary.csv")
    got = {"f00_min": float(pts[:, 0].min()), "f00_max": float(pts[:, 0].max()),
           "f11_min": float(pts[:, 1].min()), "f11_max": float(pts[:, 1].max())}
    tgt = {"f00_min": float(ref["f00_s3"].min()), "f00_max": float(ref["f00_s3"].max()),
           "f11_min": float(ref["f11_s3"].min()), "f11_max": float(ref["f11_s3"].max())}
    rows = [{"quantity": k, "ours": got[k], "paper": tgt[k],
             "rel_diff": (got[k] - tgt[k]) / abs(tgt[k]),
             "pass": abs((got[k] - tgt[k]) / tgt[k]) <= 0.02} for k in got]
    return {"claim": "C1", "rows": rows, "pass": all(r["pass"] for r in rows),
            "n_directions": int(len(pts))}


def phase_comparison(ours: dict, csv_name: str, group: str | None = None) -> dict:
    """Node-wise rms between our phase shift and one digitised curve, in degrees.

    ``group`` selects the colour series of the original figure; the pairing with
    our representative points is fixed in advance by the task:
    tip <-> red, mid <-> pink, ref <-> light_pink.
    """
    ref = load_csv(csv_name)
    x, y = ref["energy_gev"], ref["phase_deg"]
    if group is not None and "group" in ref:
        m = ref["group"] == group
        if m.sum() == 0:
            return {"error": "group not present", "groups": sorted(set(ref["group"]))}
        x, y = x[m], y[m]
    o = np.argsort(x)
    x, y = x[o], y[o]
    E = np.asarray(ours["E_GeV"])
    m = (E >= x.min()) & (E <= x.max())
    if m.sum() == 0:
        return {"n": 0}
    pred = np.interp(E[m], x, y)
    d = np.asarray(ours["delta_deg"])[m] - pred
    return {"n": int(m.sum()), "rms_deg": float(np.sqrt((d ** 2).mean())),
            "max_abs_deg": float(np.abs(d).max()), "group": group,
            "E_range": [float(E[m].min()), float(E[m].max())]}


def build_all(root: str) -> dict:
    """Assemble every table the report needs and write ``tables.json``."""
    recs = load_reports(root)
    out = {"n_records": len(recs), "energy_axis_note": ENERGY_AXIS_NOTE}
    pure = [r for r in recs if not r["spec"]["chiral"] and not r["spec"]["uv"]]
    if pure:
        out["C1"] = c1_table(pure)
    with open(os.path.join(root, "tables.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    _plots(root, recs)
    return out


def _plots(root: str, recs) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig_dir = os.path.join(root, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    pure = [r for r in recs if not r["spec"]["chiral"] and not r["spec"]["uv"]]
    if pure:
        pts = boundary_points(pure)
        ref = load_csv("figure3_boundary.csv")
        fig, ax = plt.subplots(figsize=(7, 4.4))
        ax.plot(ref["f00_s3"], ref["f11_s3"], ".", ms=2.5, color="0.6",
                label="digitised Fig. 3")
        if len(pts):
            o = np.argsort(np.arctan2(pts[:, 1] - pts[:, 1].mean(),
                                      pts[:, 0] - pts[:, 0].mean()))
            ax.plot(pts[o, 0], pts[o, 1], "o-", ms=4, lw=1, color="C0",
                    label="this work (SDP)")
        ax.set_xlabel(r"$f^0_0(s=3)$")
        ax.set_ylabel(r"$f^1_1(s=3)$")
        ax.legend(frameon=False, fontsize=8)
        ax.set_title("Pure S-matrix bootstrap region", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "fig3_region.pdf"))
        plt.close(fig)

    # ---- chiral region(s) against the digitised Fig. 8 boundaries
    chi = [r for r in recs if r["spec"]["chiral"] and not r["spec"]["uv"]
           and r["spec"]["chi_caliber"] == "chi-b"]
    if chi:
        by_eps = {}
        for r in chi:
            by_eps.setdefault(round(r["spec"]["eps_chi"], 10), []).append(r)
        ref = load_csv("figure8_boundary.csv")
        pts = load_csv("figure8_selected_points.csv")
        fig, ax = plt.subplots(figsize=(7, 4.4))
        m = ref["group"] == "chiral_only_green"
        ax.plot(ref["f00_s3"][m], ref["f11_s3"][m], ".", ms=2.5, color="0.55",
                label=r"digitised Fig. 8, chiral only ($\epsilon^\chi=2\times10^{-3}$)")
        for k, (eps, rs) in enumerate(sorted(by_eps.items(), reverse=True)):
            q = boundary_points(rs, want="dir")
            if not len(q):
                continue
            o = np.argsort(np.arctan2(q[:, 1] - q[:, 1].mean(),
                                      q[:, 0] - q[:, 0].mean()))
            ax.plot(q[o, 0], q[o, 1], "o-", ms=3, lw=1, color=f"C{k}",
                    label=r"this work, $\epsilon^\chi=%.0e$" % eps)
        for g, mk in (("chiral_reference_black", "k*"), ("red", "r^"),
                      ("pink", "v"), ("light_pink", "s")):
            mm = pts["group"] == g
            if mm.sum():
                ax.plot(pts["f00_s3"][mm], pts["f11_s3"][mm], mk, ms=6,
                        label=g.replace("_", " "))
        ax.set_xlabel(r"$f^0_0(s=3)$")
        ax.set_ylabel(r"$f^1_1(s=3)$")
        ax.legend(frameon=False, fontsize=7)
        ax.set_title("Chiral-constrained region", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "fig4_chiral_region.pdf"))
        plt.close(fig)
