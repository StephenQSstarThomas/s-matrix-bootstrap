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
import hashlib
import json
import os

import numpy as np
from .assembly import basis_identity, recorded_basis_key

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


def load_reports(root: str, subthreshold_reports=()) -> list[dict]:
    out = []
    for f in sorted(glob.glob(os.path.join(root, "**", "report.json"), recursive=True)):
        with open(f) as fh:
            doc = json.load(fh)
        # Only leaf solve reports, not aggregate copies of the same amplitudes.
        if doc.get("solver") == "SDPB" and "spec" in doc:
            v = doc.get("verification", {})
            job = os.path.basename(os.path.dirname(f))
            if job.startswith("round_"):
                job = os.path.basename(os.path.dirname(os.path.dirname(f)))
            if job == "support" and doc.get("fix_f00") is None:
                job = "dir_support"
            out.append({"job": job, "spec": doc["spec"], "solver": "SDPB",
                "fix_f00": doc.get("fix_f00"),
                "accepted": doc.get("accepted", False), "verification": v,
                "convergence": doc.get("convergence", {}),
                "basis":doc.get('basis'),"basis_identity":basis_identity(doc,f),
                "source_sha256":doc.get('source_sha256',{}),
                "result": {"objective": v.get("objective_recomputed"),
                           "f00_3": v.get("f00_3"), "f11_3": v.get("f11_3"),
                           "direction": doc.get("direction")},
                "observables": doc.get("observables", {}),
                "subthreshold": doc.get("subthreshold", {}),
                "_file": os.path.relpath(f, root), "_source_report":os.path.realpath(f)})
        for r in doc.get("records", []):
            r["_file"] = os.path.relpath(f, root)
            out.append(r)
    if subthreshold_reports:
        from .evaluation import attach_subthreshold_overlays
        out = attach_subthreshold_overlays(out,subthreshold_reports)
    return out


def accepted_support(record):
    """Complete SDPB numerical acceptance; historical inner points are not extrema."""
    return (record.get("solver") == "SDPB" and record.get("accepted", False) and
            recorded_basis_key(record) is not None and
            record.get("verification", {}).get("primal_feasible", False) and
            record.get("convergence", {}).get("solver_optimal", False))


def model_key(record):
    spec = {k: v for k, v in record["spec"].items() if k not in ("tag", "disk_mask")}
    basis=recorded_basis_key(record)
    if basis is None:raise ValueError('Reduced finite problem has no authenticated basis identity')
    return hashlib.sha256(json.dumps({'spec':spec,'basis':basis}, sort_keys=True).encode()).hexdigest()[:16]


def unique_supports(records):
    """Deduplicate repeated supports, never merge distinct finite model contracts."""
    out = {}
    for r in records:
        if not accepted_support(r) or r.get("fix_f00") is not None:
            continue
        d = np.asarray(r.get("result", {}).get("direction"), dtype=float)
        if d.shape != (2,) or not np.all(np.isfinite(d)) or np.linalg.norm(d) == 0:
            continue
        key = (model_key(r), tuple(np.round(d/np.linalg.norm(d), 12)))
        out[key] = r
    return list(out.values())


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
            if not accepted_support(r):
                continue
        pts.append((res["f00_3"], res["f11_3"]))
    return np.array(pts) if pts else np.zeros((0, 2))


def best_support(records) -> dict:
    """Best verified-feasible objective per direction, over all configurations."""
    out = {}
    for r in records:
        res = r.get("result", {})
        if res.get("objective") is None or not accepted_support(r):
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
    records = unique_supports(records)
    if len({model_key(r) for r in records}) > 1:
        return {"claim": "C1", "status": "mixed_models", "pass": False, "rows": [],
                "n_directions": 0, "reason": "Different model contracts require separate comparisons"}
    pts = boundary_points(records)
    if len(pts) < 24:
        return {"claim": "C1", "status": "incomplete", "pass": False,
                "n_directions": len(pts), "rows": [],
                "reason": "24 accepted directions required by the registered comparison"}
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


def build_all(root: str, subthreshold_reports=()) -> dict:
    """Assemble every table the report needs and write ``tables.json``."""
    recs = load_reports(root,subthreshold_reports)
    out = {"n_records": len(recs), "energy_axis_note": ENERGY_AXIS_NOTE}
    from .evaluation import subthreshold_claim_tables
    out['C3_by_source_objective'] = subthreshold_claim_tables(recs)
    pure = [r for r in recs if not r["spec"]["chiral"] and not r["spec"]["uv"]]
    if pure:
        out["C1"] = c1_table(pure)
        groups = {}
        for r in pure:
            if accepted_support(r):groups.setdefault(model_key(r), []).append(r)
        out["C1_by_model"] = {k: c1_table(rs) for k, rs in groups.items()}
    with open(os.path.join(root, "tables.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    groups = {}
    for r in recs:
        if accepted_support(r):
            groups.setdefault(model_key(r), []).append(r)
    for key, records in groups.items():
        path = os.path.join(root, "by_model", key)
        os.makedirs(path, exist_ok=True)
        _plots(path, unique_supports(records))
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
    chi = [r for r in recs if r["spec"]["chiral"] and not r["spec"]["uv"]]
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


def ir_profile_artifacts(outdir,record):
    """Saved subthreshold profiles and the data behind a selected IR comparison."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    dest=Path(outdir);data=record['subthreshold'];s=np.asarray(data['s'],float)
    x=record['verification']['f00_3']
    lines={'S0':(2*s-1)*x/5,'S2':(2-s)*x/5,'P1':(s-4)*x/15}
    fig,axes=plt.subplots(1,3,figsize=(12,3.6));subrows=[];phaserows=[]
    for ax,wave in zip(axes,('S0','S2','P1')):
        values=np.asarray(data[wave],float)
        ax.plot(s,values,label='selected IR');ax.plot(s,lines[wave],'k--',label='same-x Weinberg line')
        ax.axhline(0,color='.7',lw=.6);ax.set(title=wave,xlabel='s / m_pi^2',ylabel='partial wave f')
        ax.grid(alpha=.15)
        subrows.extend([wave,float(q),float(v),float(w)] for q,v,w in zip(s,values,lines[wave]))
        obs=record['observables'][wave]
        phaserows.extend(['selected_IR',wave,e,d,eta] for e,d,eta in zip(obs['E_GeV'],obs['delta_deg'],obs['eta']))
        ref=load_csv(f'figure7_{wave.lower()}_phases.csv');m=ref['group']=='ir_magenta'
        phaserows.extend(['paper_Fig7',wave,float(e),float(d),''] for e,d in zip(ref['energy_gev'][m],ref['phase_deg'][m]))
    axes[0].legend(frameon=False,fontsize=8)
    fig.suptitle(f"Selected IR subthreshold curves; epsilon_chi={record['spec']['eps_chi']:g}")
    fig.tight_layout();fig.savefig(dest/'subthreshold.pdf');fig.savefig(dest/'subthreshold.png',dpi=180);plt.close(fig)
    for filename,header,rows in [('subthreshold.csv',['wave','s','f','Weinberg_same_x'],subrows),
                                  ('phase_eta.csv',['source','wave','energy_GeV','phase_deg','eta'],phaserows)]:
        with (dest/filename).open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(header);writer.writerows(rows)
