"""Phase shifts, inelasticities and the rho position from a solution.

Definitions follow (2.11):  ``S^I_ell = 1 + i kappa f^I_ell = eta exp(2 i delta)``
with ``kappa = pi sqrt((s-4)/s)``, so

    eta   = |S| ,      delta = (1/2) unwrap(arg S)   along the node ladder.

Energies are reported in GeV as ``E = m_pi sqrt(s)`` with ``m_pi = 0.140 GeV``,
the paper's convention (2.2).  The digitised figures were calibrated with
``m_pi = 139.57 MeV``; that 0.3% difference in the abscissa is carried
explicitly rather than absorbed.
"""
from __future__ import annotations

import numpy as np

from .grid import M_PI_MEV, kappa, s_nodes

WAVES = {"S0": (0, 0), "S2": (2, 0), "P1": (1, 1)}


def energy_gev(s, m_pi_mev: float = M_PI_MEV) -> np.ndarray:
    return (m_pi_mev / 1000.0) * np.sqrt(np.asarray(s, dtype=float))


def wave_values(ops, c: np.ndarray, wave: str) -> np.ndarray:
    """Complex ``f^I_ell(s_i)`` on the node ladder."""
    return _wave_h(ops, c, wave) / kappa(ops.s)


def _wave_h(ops, c, wave):
    """Native kappa*f directly, preserving exact zeros of 1+i*h."""
    I, ell = WAVES[wave]
    a = ops.index.index((I, ell))
    lo, hi = a * ops.M, (a + 1) * ops.M
    return ops.h_re[lo:hi] @ c + 1j * (ops.h_im[lo:hi] @ c)


def phase_shift(ops, c: np.ndarray, wave: str) -> dict:
    h = _wave_h(ops, c, wave)
    f, S = h / kappa(ops.s), 1.0 + 1j * h
    # The threshold convention is S(4)=1, delta(4)=0 when f stays finite.
    # Only delta -> delta + k*pi preserves S.  A pi/2 shift changes S's sign.
    # This is a nearest-node lift, not a claim about winding between PV nodes.
    delta = 0.5 * np.unwrap(np.angle(np.r_[1.0 + 0j, S]))[1:]
    zeros = np.flatnonzero(S == 0)
    first_zero = int(zeros[0]) if zeros.size else None
    if first_zero is not None:
        delta[first_zero:] = np.nan          # phase cannot be anchored through a zero
    return {"s": ops.s.copy(), "E_GeV": energy_gev(ops.s), "f": f, "S": S,
            "eta": np.abs(S), "delta_deg": np.degrees(delta),
            "first_zero_node": first_zero,
            "phase_convention": "nearest-node lift from S(4)=1; no off-node winding certificate"}


def crossing_energy(ph: dict, target_deg: float = 90.0) -> float | None:
    """Linear interpolation of the first ``delta = target`` crossing, in GeV."""
    d, E = ph["delta_deg"], ph["E_GeV"]
    for i in range(len(d) - 1):
        if (d[i] - target_deg) * (d[i + 1] - target_deg) <= 0 and d[i] != d[i + 1]:
            t = (target_deg - d[i]) / (d[i + 1] - d[i])
            return float(E[i] + t * (E[i + 1] - E[i]))
    return None


def modulus_peak(ph: dict) -> float:
    """Node energy at which |f| peaks (reported next to the 90 degree crossing)."""
    return float(ph["E_GeV"][int(np.argmax(np.abs(ph["f"])))])


def resonance_report(ops, c: np.ndarray) -> dict:
    out = {}
    for w in WAVES:
        ph = phase_shift(ops, c, w)
        out[w] = {"delta_deg": [float(x) for x in ph["delta_deg"]],
                  "first_zero_node": ph["first_zero_node"],
                  "phase_convention": ph["phase_convention"],
                  "eta": [float(x) for x in ph["eta"]],
                  "E_GeV": [float(x) for x in ph["E_GeV"]],
                  "crossing_90_GeV": crossing_energy(ph, 90.0),
                  "modulus_peak_GeV": modulus_peak(ph),
                  "min_eta_below_1p2GeV": float(np.min(ph["eta"][ph["E_GeV"] <= 1.2]))}
    return out


def interpolate_at(ph: dict, e_gev: float, key: str = "delta_deg") -> float | None:
    E, y = ph["E_GeV"], ph[key]
    if e_gev < E[0] or e_gev > E[-1]:
        return None
    return float(np.interp(e_gev, E, y))


def subthreshold_curves(ops, c: np.ndarray, s_values) -> dict:
    """Real partial waves in 0 < s < 4 (Fig. 5), evaluated from the same rows."""
    out = {"s": list(map(float, s_values))}
    for w, (I, ell) in WAVES.items():
        out[w] = [float(ops.op.rows(I, ell, float(s))[0] @ c) for s in s_values]
    return out
