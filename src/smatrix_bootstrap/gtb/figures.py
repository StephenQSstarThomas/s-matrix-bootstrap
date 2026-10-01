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
