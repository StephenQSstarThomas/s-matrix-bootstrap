"""Comparison figures for the 2403 report: ours (filled, thin line) vs the paper's blind-extracted dots (hollow).

Points ni = 44/46/48 keep the paper's blue/green/red identities (validated triple, CVD warn band -> secondary encoding
by marker shape: circle/square/triangle).  One y-axis per panel; legends always present.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STYLE = {44: ("#2a78d6", "o", "ni=44 (blue)"), 46: ("#008300", "s", "ni=46 (green)"), 48: ("#e34948", "^", "ni=48 (red)")}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
WAVE_TITLES = {"S0": r"$\delta^0_0$", "D0": r"$\delta^0_2$", "S2": r"$\delta^2_0$", "D2": r"$\delta^2_2$",
               "P1": r"$\delta^1_1$", "F1": r"$\delta^1_3$"}


def _axes_style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(GRID)


def phases_figure(ours: dict, paper: dict, path, quantity="delta_deg", ylabel=r"$\delta$ (deg)"):
    """ours/paper: {ni: analysis dict with 'E_GeV' and 'phases'[wave][quantity]}; six panels."""
    fig, axes = plt.subplots(3, 2, figsize=(9, 9.5), sharex=True)
    for ax, w in zip(axes.flat, ("S0", "D0", "S2", "D2", "P1", "F1")):
        for ni, (col, mk, lab) in STYLE.items():
            if ni in paper:
                p = paper[ni]
                ax.plot(p["E_GeV"], p["phases"][w][quantity], ls="none", marker=mk, ms=5, mfc="none", mec=col, mew=1.0,
                        label=f"paper {lab}")
            if ni in ours:
                o = ours[ni]
                n = len(o["phases"][w][quantity])
                ax.plot(o["E_GeV"][:n], o["phases"][w][quantity], lw=1.2, color=col, marker=mk, ms=3, label=f"ours {lab}")
        ax.set_title(WAVE_TITLES[w] if quantity == "delta_deg" else WAVE_TITLES[w].replace(r"\delta", r"\eta"),
                     color=INK, fontsize=10)
        _axes_style(ax)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\sqrt{s}$ (GeV)", color=MUTED, fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def region_figure(ours: dict, c1: dict, dots: dict, path, f0_scale=0.163108):
    """Fig. 3: allowed region in (F0, F1).  ours: {k: {'max', 'min'}} (V1-certified, reduced basis = inner bounds);
    c1: {ni: certified F1max}; dots: appclaims.fig3_dots of the paper.  Paper hollow, ours filled."""
    from pathlib import Path
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ks = sorted(dots["edges"])
    for s, lab in (("max", "paper edges (Fig. 3)"), ("min", None)):
        ax.plot([k / 50 * 0.1632 for k in ks], [dots["edges"][k][s] for k in ks], ls="none", marker="o", ms=4,
                mfc="none", mec=MUTED, mew=0.8, label=lab)
    kk = sorted(k for k in ours if ours[k].get("max") is not None)
    for s, lab in (("max", "ours: max / min F1 (V1, inner bounds)"), ("min", None)):
        ax.plot([k / 50 * f0_scale for k in kk], [ours[k][s] for k in kk], lw=1.2, color=INK, marker="o", ms=4,
                label=lab)
    for ni, (col, mk, lab) in STYLE.items():
        p = dots["points"].get({44: "blue", 46: "green", 48: "red"}[ni])
        if p:
            ax.plot(*p, ls="none", marker=mk, ms=8, mfc="none", mec=col, mew=1.4, label=f"paper {lab}")
        if ni in c1:
            ax.plot(ni / 50 * f0_scale, c1[ni], ls="none", marker=mk, ms=7, color=col, label=f"ours {lab}")
    if "black" in dots["points"]:
        ax.plot(*dots["points"]["black"], ls="none", marker="*", ms=11, color=INK, label="tree level")
    ax.set_xlabel(r"$\mathfrak{F}_0$", color=MUTED)
    ax.set_ylabel(r"$\mathfrak{F}_1$", color=MUTED)
    _axes_style(ax)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)
