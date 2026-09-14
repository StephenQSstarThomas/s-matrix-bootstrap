"""Side-by-side comparison figures: the paper's figure (rendered PDF) on the left, our reproduction with the
digitised paper curve overlaid on the right.  Every number plotted comes from an accepted leaf report.json.

Usage: python scripts/sdp/final_figures.py --root RESULTS_ROOT --paper-png RESULTS_ROOT/paper_figures --out RESULTS_ROOT/final_figures
Missing leaves are skipped and listed in figures.json.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

REF = Path(__file__).resolve().parents[2] / "references"
GREY = "0.45"


def load_csv(name):
    rows = list(csv.DictReader(open(REF / name)))
    out = {}
    for k in rows[0]:
        col = [r[k] for r in rows]
        try:
            out[k] = np.array([float(v) for v in col])
        except ValueError:
            out[k] = np.array(col)
    return out


def leaf(root, rel):
    p = Path(root) / rel / "report.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    return r if r.get("accepted") else None


def two_panel(paper_png, title):
    fig = plt.figure(figsize=(12, 4.6))
    axl = fig.add_axes([0.02, 0.05, 0.46, 0.85]); axr = fig.add_axes([0.56, 0.13, 0.42, 0.77])
    if paper_png and Path(paper_png).exists():
        axl.imshow(mpimg.imread(paper_png)); axl.set_title("paper (2309.12402 v3)", fontsize=10, color=GREY)
    axl.axis("off"); axr.set_title(title, fontsize=10)
    return fig, axr


def polygon_order(pts):
    c = pts.mean(axis=0)
    return np.argsort(np.arctan2(pts[:, 1] - c[1], pts[:, 0] - c[0]))


def fig3(root, png, out, info):
    pts = []
    for rel in ["fig3_pure_linf_1e3_tip_v2/tip"] + sorted(p.replace(root + "/", "") for p in glob.glob(f"{root}/fig3_dir??")):
        r = leaf(root, rel)
        if r:
            v = r["verification"]; pts.append((v["f00_3"], v["f11_3"]))
    ref = load_csv("figure3_boundary.csv")
    fig, ax = two_panel(png / "purSplot.png", "Fig.3 pure unitarity region, M=50 L=10, Mreg=1e3 (ours)")
    ax.plot(ref["f00_s3"], ref["f11_s3"], ".", ms=2.5, color=GREY, label="digitised Fig.3")
    if pts:
        P = np.array(pts)
        if len(P) > 2:
            o = polygon_order(P); ax.plot(P[o, 0], P[o, 1], "o-", ms=4, lw=1, color="C0", label=f"this work, {len(P)} support points")
        else:
            ax.plot(P[:, 0], P[:, 1], "o", ms=6, color="C0", label="this work (+x tip)")
    ax.set_xlabel("$f_0^0(3)$"); ax.set_ylabel("$f_1^1(3)$"); ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "fig3.png", dpi=110); plt.close(fig)
    info["fig3"] = {"n_points": len(pts), "tip_x": pts[0][0] if pts else None, "paper_tip_x": float(ref["f00_s3"].max())}


def fig4(root, png, out, info):
    ref = load_csv("figure8_boundary.csv"); sel = load_csv("figure8_selected_points.csv")
    fig, ax = two_panel(png / "chiplot.png", "Fig.4 chiral-constrained region (ours, Mreg=1e2)")
    m = ref["group"] == "chiral_only_green"
    ax.plot(ref["f00_s3"][m], ref["f11_s3"][m], ".", ms=2.5, color=GREY, label=r"digitised $\epsilon^\chi=2\times10^{-3}$ (Fig.8 green)")
    tips = {}
    for tag, eps in (("0002", 2e-4), ("0006", 6e-4), ("001", 1e-3), ("002", 2e-3), ("004", 4e-3), ("006", 6e-3)):
        r = leaf(root, f"fig4_eps{tag}_tip/tip")
        if r:
            tips[eps] = (r["verification"]["f00_3"], r["verification"]["f11_3"])
    for k, (eps, (x, y)) in enumerate(sorted(tips.items())):
        ax.plot([x], [y], "o", ms=5, color=f"C{k}", label=r"+x end, $\epsilon^\chi$=%g" % eps)
    sec = []
    for tag in ("xm002", "xm001", "", "xp001", "xp002"):
        for side in ("hi", "lo"):
            rel = f"fig4_eps002_section_{tag + '_' if tag else ''}{side}"
            r = leaf(root, rel)
            if r:
                sec.append((r["verification"]["f00_3"], r["verification"]["f11_3"]))
    if sec:
        S = np.array(sec); ax.plot(S[:, 0], S[:, 1], "s", ms=4, color="C3", label=r"$\epsilon^\chi=2\times10^{-3}$ sections (ours)")
    xs = np.linspace(0, 0.09, 10); ax.plot(xs, -xs / 15, "k-", lw=0.8, label="Weinberg line $f_1^1=-f_0^0/15$")
    mm = sel["group"] == "chiral_reference_black"; ax.plot(sel["f00_s3"][mm], sel["f11_s3"][mm], "k*", ms=9, label="black dot ($f_\\pi$=92 MeV)")
    ax.set_xlabel("$f_0^0(3)$"); ax.set_ylabel("$f_1^1(3)$"); ax.legend(frameon=False, fontsize=7)
    fig.savefig(out / "fig4.png", dpi=110); plt.close(fig)
    info["fig4"] = {"tips": {str(e): v for e, v in tips.items()}, "n_sections": len(sec)}


def fig5(root, png, out, info):
    ref = load_csv("figure5_subthreshold.csv")
    fig = plt.figure(figsize=(12, 4.6)); fig.suptitle("Fig.5 subthreshold partial waves: paper (left, three panels) vs ours (right)", fontsize=10)
    for j, name in enumerate(("lin1", "lin2", "lin3")):
        axp = fig.add_axes([0.01 + 0.155 * j, 0.08, 0.15, 0.8]); p = png / f"{name}.png"
        if p.exists(): axp.imshow(mpimg.imread(p))
        axp.axis("off")
    axes = [fig.add_axes([0.50 + 0.165 * j, 0.14, 0.15, 0.72]) for j in range(3)]
    rows = {}
    for j, wave in enumerate(("S0", "P1", "S2")):
        ax = axes[j]
        for k, (tag, eps) in enumerate((("002", 0.002), ("004", 0.004), ("006", 0.006))):
            m = (ref["wave"] == wave) & (np.isclose(ref["epsilon"], eps))
            ax.plot(ref["s"][m], ref["f"][m], "-", lw=2.2, color=f"C{k}", alpha=0.35, label=f"paper $\\epsilon$={eps}" if j == 0 else None)
            r = leaf(root, f"fig4_eps{tag}_section_hi")
            if r:
                sub = r["subthreshold"]; ax.plot(sub["s"], sub[wave], "-", lw=1, color=f"C{k}", label=f"ours $\\epsilon$={eps}" if j == 0 else None)
                rows[f"{wave}_{eps}"] = True
        ax.set_title(wave, fontsize=9); ax.set_xlabel("s"); ax.set_xlim(0, 4)
    axes[0].legend(frameon=False, fontsize=6)
    fig.savefig(out / "fig5.png", dpi=110); plt.close(fig); info["fig5"] = {"curves": sorted(rows)}


def phases_panel(ax, r, wave, csv_name, groups, label_ours, eta=False):
    o = r["observables"][wave]; E = np.array(o["E_GeV"]); d = np.array(o["delta_deg"])
    ref = load_csv(csv_name)
    for g, col in groups:
        m = (ref["group"] == g) if "group" in ref else np.ones(len(ref["energy_gev"]), bool)
        ax.plot(ref["energy_gev"][m], ref["phase_deg"][m], ".", ms=3, color=col, alpha=0.7, label=f"paper {g}")
    m = E <= 1.3; ax.plot(E[m], d[m], "o-", ms=3, lw=1, color="C0", label=label_ours)
    if eta:
        ax2 = ax.twinx(); ax2.plot(E[m], np.array(o["eta"])[m], "s--", ms=2, lw=0.8, color="C3", label=r"$\eta$ (ours)"); ax2.set_ylim(0, 1.05); ax2.set_ylabel(r"$\eta$", color="C3")
    ax.set_xlabel("E (GeV)"); ax.set_ylabel("phase shift (deg)"); ax.set_xlim(0.28, 1.3)


def fig7(root, png, out, info):
    r = leaf(root, "fig4_eps002_section_hi")
    fig = plt.figure(figsize=(12, 4.6)); fig.suptitle("Fig.7 chiral-only phase shifts at the representative near the black dot", fontsize=10)
    for j, (wave, name, csvn) in enumerate((("S0", "chiS0plot", "figure7_s0_phases.csv"), ("S2", "chiS2plot", "figure7_s2_phases.csv"), ("P1", "chiP1plot", "figure7_p1_phases.csv"))):
        axp = fig.add_axes([0.01 + 0.155 * j, 0.08, 0.15, 0.8]); p = png / f"{name}.png"
        if p.exists(): axp.imshow(mpimg.imread(p))
        axp.axis("off")
        ax = fig.add_axes([0.50 + 0.165 * j, 0.14, 0.15, 0.72]); ax.set_title(wave, fontsize=9)
        if r and (REF / csvn).exists():
            phases_panel(ax, r, wave, csvn, [("ir_magenta", "m")], "ours (x_ref upper)")
        if j == 0: ax.legend(frameon=False, fontsize=6)
    fig.savefig(out / "fig7.png", dpi=110); plt.close(fig); info["fig7"] = {"leaf": "fig4_eps002_section_hi" if r else None}


def fig8(root, png, out, info):
    ref = load_csv("figure8_boundary.csv"); sel = load_csv("figure8_selected_points.csv")
    fig, ax = two_panel(png / "chiSRplot.png", "Fig.8 chiral (green) vs chiral+UV (cyan): ours, SR-a raw box, node FF factor")
    for g, col in (("chiral_only_green", "g"), ("gauge_cyan", "c")):
        m = ref["group"] == g; ax.plot(ref["f00_s3"][m], ref["f11_s3"][m], ".", ms=2.5, color=col, alpha=0.5, label=f"digitised {g}")
    chi, uv = [], []
    for tag in ("xm002", "xm001", "", "xp001", "xp002"):
        for side in ("hi", "lo"):
            r = leaf(root, f"fig4_eps002_section_{tag + '_' if tag else ''}{side}")
            if r: chi.append((r["verification"]["f00_3"], r["verification"]["f11_3"]))
            r = leaf(root, f"uv_SRa_section_{tag + '_' if tag else ''}{side}")
            if r: uv.append((r["verification"]["f00_3"], r["verification"]["f11_3"]))
    r = leaf(root, "fig4_eps002_tip/tip"); chi.append((r["verification"]["f00_3"], r["verification"]["f11_3"])) if r else None
    r = leaf(root, "uv_SRa_tip/tip"); uv.append((r["verification"]["f00_3"], r["verification"]["f11_3"])) if r else None
    if chi: C = np.array(chi); ax.plot(C[:, 0], C[:, 1], "s", ms=5, color="g", label="ours chiral (sections + tip)")
    if uv: U = np.array(uv); ax.plot(U[:, 0], U[:, 1], "^", ms=6, color="c", markeredgecolor="k", label="ours chiral+UV (sections + tip)")
    for g, mk in (("chiral_reference_black", "k*"), ("red", "r^"), ("pink", "v"), ("light_pink", "s")):
        mm = sel["group"] == g
        if mm.sum(): ax.plot(sel["f00_s3"][mm], sel["f11_s3"][mm], mk, ms=7, label=f"paper {g}")
    ax.set_xlabel("$f_0^0(3)$"); ax.set_ylabel("$f_1^1(3)$"); ax.legend(frameon=False, fontsize=6, ncol=2, loc="upper right")
    ins = ax.inset_axes([0.08, 0.08, 0.42, 0.5])
    for g, col in (("chiral_only_green", "g"), ("gauge_cyan", "c")):
        m = ref["group"] == g; ins.plot(ref["f00_s3"][m], ref["f11_s3"][m], ".", ms=3, color=col, alpha=0.5)
    if chi: ins.plot(C[:, 0], C[:, 1], "s", ms=5, color="g")
    if uv: ins.plot(U[:, 0], U[:, 1], "^", ms=6, color="c", markeredgecolor="k")
    for g, mk in (("chiral_reference_black", "k*"), ("red", "r^"), ("pink", "v"), ("light_pink", "s")):
        mm = sel["group"] == g
        if mm.sum(): ins.plot(sel["f00_s3"][mm], sel["f11_s3"][mm], mk, ms=7)
    ins.set_xlim(0.069, 0.083); ins.set_ylim(-0.0056, -0.0040); ins.tick_params(labelsize=6); ins.set_title("zoom: x_ref region", fontsize=7)
    fig.savefig(out / "fig8.png", dpi=110); plt.close(fig); info["fig8"] = {"n_chi": len(chi), "n_uv": len(uv)}


def fig9_10(root, png, out, info):
    pts = [("uv_SRa_tip/tip", "tip (+x end)"), ("uv_SRa_section_hi", "x_ref upper")]
    for tag in ("xm001", "xp001", "xm002", "xp002"):
        pts.append((f"uv_SRa_section_{tag}_hi", f"x_ref{tag[1:].replace('m','-').replace('p','+')} upper"))
    fig = plt.figure(figsize=(12, 4.6)); fig.suptitle("Fig.9 P1 phase shift after the QCD sum rules (ours: phase and inelasticity)", fontsize=10)
    axp = fig.add_axes([0.01, 0.05, 0.47, 0.85]); p = png / "P1ps.png"
    if p.exists(): axp.imshow(mpimg.imread(p))
    axp.axis("off"); ax = fig.add_axes([0.55, 0.13, 0.36, 0.77])
    ref = load_csv("figure9_p1_phases.csv")
    for g, col in (("red", "r"), ("pink", "orchid"), ("light_pink", "pink")):
        m = ref["group"] == g; ax.plot(ref["energy_gev"][m], ref["phase_deg"][m], ".", ms=3, color=col, alpha=0.8, label=f"paper {g}")
    used = []
    for k, (rel, lab) in enumerate(pts):
        r = leaf(root, rel)
        if not r: continue
        o = r["observables"]["P1"]; E = np.array(o["E_GeV"]); m = E <= 1.3; d = np.array(o["delta_deg"]); eta = np.array(o["eta"])
        ax.plot(E[m], d[m], "o-", ms=3, lw=1, color=f"C{k}", label=f"ours {lab}")
        i = int(np.argmin(eta[30:45])) + 30
        if d[i] < d[i - 1] and eta[i] < 0.6:      # |S| dip with a downward jump: the +180 deg branch is equally consistent
            alt = d.copy(); alt[i:] += 180
            ax.plot(E[m], alt[m], "--", lw=0.9, color=f"C{k}", alpha=0.7, label=f"ours {lab}, +180° branch")
        used.append({"leaf": rel, "crossing_90_GeV": o["crossing_90_GeV"], "min_eta": o["min_eta_below_1p2GeV"]})
    ax.axhline(90, color=GREY, lw=0.5, ls=":"); ax.set_xlabel("E (GeV)"); ax.set_ylabel("$\\delta_{P1}$ (deg)"); ax.set_xlim(0.28, 1.3); ax.legend(frameon=False, fontsize=6)
    ax2 = fig.add_axes([0.93, 0.13, 0.05, 0.77]); ax2.axis("off")
    fig.savefig(out / "fig9.png", dpi=110); plt.close(fig); info["fig9"] = used
    # eta panel (our extra information)
    fig, ax = plt.subplots(figsize=(6, 3.4))
    for k, (rel, lab) in enumerate(pts):
        r = leaf(root, rel)
        if not r: continue
        o = r["observables"]["P1"]; E = np.array(o["E_GeV"]); m = E <= 1.3
        ax.plot(E[m], np.array(o["eta"])[m], "s-", ms=3, lw=1, color=f"C{k}", label=lab)
    ax.set_xlabel("E (GeV)"); ax.set_ylabel(r"$\eta_{P1}=|S_{P1}|$ (ours; not shown in the paper)"); ax.set_ylim(0, 1.05); ax.legend(frameon=False, fontsize=7); fig.tight_layout()
    fig.savefig(out / "fig9_eta.png", dpi=110); plt.close(fig)
    fig = plt.figure(figsize=(12, 4.6)); fig.suptitle("Fig.10 S0 and S2 phase shifts at the same points", fontsize=10)
    for j, (wave, name, csvn) in enumerate((("S0", "S0ps", "figure10_s0_phases.csv"), ("S2", "S2ps", "figure10_s2_phases.csv"))):
        axp = fig.add_axes([0.01 + 0.24 * j, 0.08, 0.23, 0.8]); p = png / f"{name}.png"
        if p.exists(): axp.imshow(mpimg.imread(p))
        axp.axis("off"); ax = fig.add_axes([0.53 + 0.24 * j, 0.14, 0.21, 0.72]); ax.set_title(wave, fontsize=9)
        ref = load_csv(csvn)
        for g, col in (("red", "r"), ("pink", "orchid"), ("light_pink", "pink")):
            m = ref["group"] == g; ax.plot(ref["energy_gev"][m], ref["phase_deg"][m], ".", ms=3, color=col, alpha=0.8, label=f"paper {g}" if j == 0 else None)
        for k, (rel, lab) in enumerate(pts):
            r = leaf(root, rel)
            if not r: continue
            o = r["observables"][wave]; E = np.array(o["E_GeV"]); m = E <= 1.3
            ax.plot(E[m], np.array(o["delta_deg"])[m], "o-", ms=3, lw=1, color=f"C{k}", label=f"ours {lab}" if j == 0 else None)
        ax.set_xlabel("E (GeV)"); ax.set_xlim(0.28, 1.3)
        if j == 0: ax.legend(frameon=False, fontsize=6)
    fig.savefig(out / "fig10.png", dpi=110); plt.close(fig)


def fig11(root, png, out, info):
    models = [(50, 10, "uv_SRa_tip/tip"), (50, 8, "uv_M50_L8_tip/tip"), (50, 12, "uv_M50_L12_tip/tip"), (45, 10, "uv_M45_L10_tip/tip"), (60, 10, "uv_M60_L10_tip/tip")]
    fig = plt.figure(figsize=(12, 4.6)); fig.suptitle("Fig.11 dependence on M and L (UV tip, ours vs paper)", fontsize=10); used = []
    for j, (wave, name, csvn) in enumerate((("S0", "S0ML", "figure11_s0_phases.csv"), ("P1", "P1ML", "figure11_p1_phases.csv"), ("S2", "S2ML", "figure11_s2_phases.csv"))):
        axp = fig.add_axes([0.01 + 0.155 * j, 0.08, 0.15, 0.8]); p = png / f"{name}.png"
        if p.exists(): axp.imshow(mpimg.imread(p))
        axp.axis("off"); ax = fig.add_axes([0.50 + 0.165 * j, 0.14, 0.15, 0.72]); ax.set_title(wave, fontsize=9)
        if (REF / csvn).exists():
            ref = load_csv(csvn)
            for (M, L) in sorted(set(zip(ref["model_M"].astype(int), ref["model_L"].astype(int)))):
                m = (ref["model_M"] == M) & (ref["model_L"] == L); ax.plot(ref["energy_gev"][m], ref["phase_deg"][m], ".", ms=2, alpha=0.5, label=f"paper M{M} L{L}" if j == 0 else None)
        for k, (M, L, rel) in enumerate(models):
            r = leaf(root, rel)
            if not r: continue
            o = r["observables"][wave]; E = np.array(o["E_GeV"]); m = E <= 1.3
            ax.plot(E[m], np.array(o["delta_deg"])[m], "-", lw=1.2, color=f"C{k}", label=f"ours M{M} L{L}" if j == 0 else None)
            if j == 1: used.append({"M": M, "L": L, "P1_crossing_GeV": o["crossing_90_GeV"], "x_tip": r["verification"]["f00_3"]})
        ax.set_xlabel("E (GeV)"); ax.set_xlim(0.28, 1.3)
        if j == 0: ax.legend(frameon=False, fontsize=5, ncol=2)
    fig.savefig(out / "fig11.png", dpi=110); plt.close(fig); info["fig11"] = used


def diagnostics(root, out, info):
    face = Path(root) / "FACE_RESULT_partial.json"
    if not face.exists(): return
    F = json.loads(face.read_text()); fig, ax = plt.subplots(figsize=(8, 3.6)); y = 0; labels = []
    for src, e in F["sources"].items():
        name = Path(src).parent.parent.name if Path(src).parent.name == "tip" else Path(src).parent.name
        for fn, row in e["functionals"].items():
            lo, hi = row["range_lo"], row["range_hi"]
            if hi is None: continue
            lo = hi if lo is None else lo
            ax.plot([lo, hi], [y, y], "-", lw=6, color="C0", alpha=0.6); ax.plot([hi], [y], "|", color="k")
            labels.append(f"{name}: {fn}"); y += 1
    ax.axvline(1.0, color="r", ls="--", lw=1, label="1 - Re S = 1: floor required by the paper's phase for every eta (nodes 38, 39)")
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7); ax.set_xlabel("range over the near-optimal face (margin 2e-6)"); ax.legend(frameon=False, fontsize=7); fig.tight_layout()
    fig.savefig(out / "face_ranges.png", dpi=110); plt.close(fig); info["face_ranges"] = labels


def main(argv=None):
    p = argparse.ArgumentParser(); p.add_argument("--root", required=True); p.add_argument("--paper-png", required=True); p.add_argument("--out", required=True)
    a = p.parse_args(argv); out = Path(a.out); out.mkdir(parents=True, exist_ok=True); png = Path(a.paper_png); info = {}
    for fn in (fig3, fig4, fig5, fig7, fig8, fig9_10, fig11):
        try:
            fn(a.root, png, out, info)
        except Exception as exc:  # keep going; record the failure
            info[fn.__name__ + "_error"] = repr(exc)
    diagnostics(a.root, out, info)
    (out / "figures.json").write_text(json.dumps(info, indent=1, default=float))
    print(json.dumps({k: (v if not isinstance(v, (list, dict)) else "...") for k, v in info.items()}, indent=0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
