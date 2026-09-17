"""Linearised unitarity-saturation (Watson) step, after the second CVX block of the authors' released code.

The authors select the amplitude they plot by an iterative procedure that improves the saturation of unitarity
(|S| -> 1 with arg F = delta at the nodes below s0).  Their released code (arXiv:2403.10772, block 2) keeps every
constraint, drops the projection equalities, and maximises a functional that is linear in the new variables and built
from the previous solution.  We implement the same idea for the 2309 model in the following declared form.

For a wave with node values h_k = kappa f(s_k) (so S_k = 1 + i h_k) the unitarity disk is |h_k|^2 <= 2 Im h_k, the
disc of radius 1 centred at i in the h-plane.  Given a target t_k on its boundary, maximising

    Re(conj(t_k) h_k) - Im h_k  =  Re(conj(t_k - i) h_k)

over the disc returns h_k = t_k exactly, so the fixed point of the iteration is the saturated amplitude with the
target phases.  The targets come from the previous accepted solution:

* S0, P1 (waves with a current): S_k^t = F_k / conj(F_k), i.e. Watson's theorem arg F = delta at |S| = 1, with
  F_k = 1 + (K Im F)_k + i Im F_k the previous form factor; t_k = -i (S_k^t - 1).
* S2 (no current): the previous h_k projected radially onto the boundary, t_k = i + (h_k - i)/|h_k - i|,
  i.e. |S| -> 1 at the previous phase.

The objective is the plain sum over the nodes with s_k <= s0 of the chosen waves; no weights.  Every other
constraint of the model is kept; the section f00(3) = x is released unless asked for (the authors' block 2 releases
the projection equalities).  This is a selection rule for the amplitude, not a change of the feasible region.

Relation to the authors' Watsonian functional, eq. (2.29) of arXiv:2403.10772,
F_W = sum_{l,I} int ds Re[ e^{-2 i alpha(s)}|_old (S(s) - 1) ], alpha = phase of the old form factor (old phase
shift where no form factor exists): with S - 1 = i h and t = 2 sin(alpha) e^{i alpha} one has identically
Re(conj(t) h) - Im h = Re[(conj(t) + i) h] = Re[(sin 2alpha + i cos 2alpha) h] = Re[e^{-2 i alpha} (S - 1)],
so the node functional here is (2.29) evaluated at the collocation nodes below s0, as in their released code
(node sums, no ds weights).  The S2 target (radial projection of the old S) gives alpha = old delta, their rule for
waves without a current.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import constraints as C
from . import formfactor as FFM
from .grid import s_nodes

WAVE_INDEX = {"S0": (0, 0), "P1": (1, 1), "S2": (2, 0)}
DEFAULT_WAVES = ("S0", "P1", "S2")


def _report_and_solution(report_path):
    report_path = Path(report_path).resolve()
    r = json.loads(report_path.read_text())
    sol = np.load(report_path.parent / "solution.npz")
    return r, sol


def node_S(report, wave):
    o = report["observables"][wave]
    eta, d = np.asarray(o["eta"], dtype=float), np.radians(np.asarray(o["delta_deg"], dtype=float))
    return eta * np.exp(2j * d)


def form_factor(sol, ell, M):
    ImF = np.asarray(sol["ImF"][ell], dtype=float)
    return 1.0 + FFM.hilbert_kernel(M) @ ImF + 1j * ImF


def targets_from_leaf(report_path, waves=DEFAULT_WAVES):
    """Targets ``t_k`` per wave for the nodes below s0, plus saturation diagnostics of the source."""
    r, sol = _report_and_solution(report_path)
    M = r["spec"]["M"]
    nodes = [int(k) for k in np.flatnonzero(C.below_s0(M))]
    out = {"source": str(Path(report_path).resolve()), "nodes": nodes, "targets": {}, "source_saturation": {}}
    for wave in waves:
        I, ell = WAVE_INDEX[wave]
        S = node_S(r, wave)
        h = -1j * (S - 1.0)
        rad = h - 1j
        norm = np.abs(rad)
        radial = 1j + np.where(norm > 1e-14, rad / np.where(norm > 1e-14, norm, 1.0), 0.0)   # old phase, |S| -> 1
        if wave in ("S0", "P1") and r["spec"].get("uv"):
            F = form_factor(sol, ell, M)
            phase_ok = np.abs(F) > 1e-8                       # a vanishing form factor has no phase: fall back to the old delta
            St = np.where(phase_ok, F / np.where(phase_ok, np.conj(F), 1.0), 1.0)
            t = np.where(phase_ok, -1j * (St - 1.0), radial)
        else:
            t = radial
        out["targets"][wave] = [[float(t[k].real), float(t[k].imag)] for k in nodes]
        elastic = 1.0 - np.abs(S[nodes]) ** 2
        out["source_saturation"][wave] = {
            "mean_1_minus_S2": float(np.mean(elastic)), "min_abs_S": float(np.min(np.abs(S[nodes]))),
            "n_nodes_abs_S_below_0p99": int(np.sum(np.abs(S[nodes]) < 0.99)),
            "max_abs_h_minus_target": float(np.max(np.abs(h[nodes] - t[nodes])))}
        if wave in ("S0", "P1") and r["spec"].get("uv"):
            F = form_factor(sol, ell, M)
            res = np.abs(F - S * np.conj(F)) / np.maximum(np.abs(F), 1e-12)
            out["source_saturation"][wave]["watson_residual_max"] = float(np.max(res[nodes]))
            out["source_saturation"][wave]["watson_residual_median"] = float(np.median(res[nodes]))
    return out


def functional_from_leaf(report_path, waves=DEFAULT_WAVES, keep_section=False):
    """The registered ``functional`` dict for one Watson step."""
    t = targets_from_leaf(report_path, waves)
    return {"kind": "watson", "waves": list(waves), "sense": "max", "node": 0, "wave": "all",
            "nodes": t["nodes"], "targets": t["targets"], "linearisation_source": t["source"],
            "source_saturation": t["source_saturation"], "keep_section": bool(keep_section),
            "convention": "maximise sum_k Re(conj(t_k) h_k) - Im h_k over the nodes with s_k <= s0 of the listed waves; "
                          "t from the previous leaf (S0/P1: F/F*, S2: radial projection of the previous S); "
                          "unit node weights; section released unless keep_section"}


def objective_value(functional, S_by_wave):
    """The functional evaluated on node S-matrix values ``S_by_wave[wave][k]`` (any float/complex-like)."""
    total = 0.0
    for wave in functional["waves"]:
        for (tr, ti), k in zip(functional["targets"][wave], functional["nodes"]):
            S = S_by_wave[wave][k]
            total += tr * float(np.imag(S)) + (ti - 1.0) * (1.0 - float(np.real(S)))
    return total


def iteration_metrics(report_path, waves=DEFAULT_WAVES):
    """Saturation and Watson metrics of one accepted leaf (for the iteration log)."""
    r, sol = _report_and_solution(report_path)
    M = r["spec"]["M"]; nodes = np.flatnonzero(C.below_s0(M))
    out = {"f00_3": r["verification"]["f00_3"], "f11_3": r["verification"]["f11_3"], "waves": {}}
    for wave in waves:
        I, ell = WAVE_INDEX[wave]
        S = node_S(r, wave); o = r["observables"][wave]
        row = {"min_abs_S_below_s0": float(np.min(np.abs(S[nodes]))), "mean_1_minus_S2": float(np.mean(1 - np.abs(S[nodes]) ** 2)),
               "crossing_90_GeV": o.get("crossing_90_GeV"), "delta_deg_at_0p792": float(o["delta_deg"][38]) if M == 50 else None,
               "abs_S_at_0p792": float(abs(S[38])) if M == 50 else None}
        if wave in ("S0", "P1") and r["spec"].get("uv"):
            F = form_factor(sol, ell, M); rh = np.asarray(sol["rho_hat"][ell], dtype=float)
            res = np.abs(F - S * np.conj(F)) / np.maximum(np.abs(F), 1e-12)
            frac = np.abs(F) ** 2 / np.where(rh > 0, rh, np.nan)
            row.update(watson_residual_max=float(np.max(res[nodes])), watson_residual_median=float(np.median(res[nodes])),
                       two_pion_fraction_min=float(np.nanmin(frac[nodes])), two_pion_fraction_median=float(np.nanmedian(frac[nodes])),
                       abs_F_at_0p792=float(abs(F[38])) if M == 50 else None)
        out["waves"][wave] = row
    return out
