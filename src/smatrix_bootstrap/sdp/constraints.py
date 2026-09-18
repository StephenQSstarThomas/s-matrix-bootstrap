"""FESR moments, form-factor asymptotics and chiral residuals.

FESR (3.72)-(3.73): the moments are discretised with the same midpoint rule in
phi that (3.61) induces, restricted to the nodes below ``s_0``,

    int_4^{s0} rho(x) x^n dx  ->  (pi/M) sum_{i: s_i <= s0} (ds/dphi)_i s_i^n rho_i .

Targets.  The paper prints normalised values in (2.56); the raw targets are
those times ``s0^{n+2}``.  :func:`recomputed_targets` redoes them from (2.50)
with the QCD inputs of (2.53)-(2.54), which is how the quark-mass convention is
pinned down only approximately (see ``fesr_target_audit``): the printed S0 numbers are closer to
``m_q = sqrt((m_u^2+m_d^2)/2)`` (the combination that actually appears in
``N_f m_q^2 -> sum_f m_f^2``), not to the arithmetic mean.

Form-factor asymptotics (3.75), imposed on the 7 nodes above ``s_0``:

    |cF_0(s_i)|^2 <= 2 m_q^2 eps_FF ,     |cF_1(s_i)|^2 <= eps_FF / 2 .

Chiral constraints (3.63)-(3.64) at ``s_j = 1/2, 1, 3/2, 2``:

    r01(s) = f^0_0(s) - [3(2s-1)/(s-4)] f^1_1(s),
    r21(s) = f^2_0(s) - [3(2-s)/(s-4)] f^1_1(s).
"""
from __future__ import annotations

import numpy as np

from .grid import (ALPHA_S, M_D, M_PI_MEV, M_Q, M_U, N_F, S0, quad_weights,
                   s_nodes)

CHIRAL_POINTS = (0.5, 1.0, 1.5, 2.0)
EPS_CHI_MAIN = 2.0e-3
EPS_CHI_GRID = (6e-3, 4e-3, 2e-3, 1e-3, 6e-4, 2e-4)
EPS_SR = 2.0e-3
EPS_FF = 6.0e-5
MOMENTS = {0: (0, 1), 1: (-1, 0)}          # n per wave, as stated below (2.52)

_GEV4 = (1000.0 / M_PI_MEV) ** 4           # 1 GeV^4 in units m_pi = 1
COND_G2 = 0.023 * _GEV4                    # <alpha_s G^2 / pi>      (2.54)
COND_JS = -(0.1 ** 4) * _GEV4              # <j_S>                   (2.54)
M_Q_RMS = float(np.sqrt((M_U ** 2 + M_D ** 2) / 2.0))


def below_s0(M: int) -> np.ndarray:
    return s_nodes(M) <= S0


def moment_row(M: int, n: int) -> np.ndarray:
    """Row ``m`` with ``M_n = m . rho`` over the M spectral-density variables."""
    s, w = s_nodes(M), quad_weights(M)
    return np.where(below_s0(M), w * s ** float(n), 0.0)


def printed_targets() -> dict[tuple[str, int], float]:
    """Raw (un-normalised) FESR targets from (2.56), i.e. printed * s0^{n+2}."""
    out = {}
    for n in MOMENTS[0]:
        out[("S0", n)] = 3.09e-8 * (27.38 / (n + 2) + (0.61 if n == 0 else 0.0)) * S0 ** (n + 2)
    for n in MOMENTS[1]:
        out[("P1", n)] = -4.34e-6 * (-13.26 / (n + 2) + (0.41 if n == 0 else 0.0)) * S0 ** (n + 2)
    return out


def recomputed_targets(m_q: float) -> dict[tuple[str, int], float]:
    """Raw FESR targets recomputed from (2.50) with the inputs (2.53)-(2.54)."""
    out = {}
    for n in MOMENTS[0]:
        d = 1.0 if n == 0 else 0.0
        out[("S0", n)] = (S0 ** (n + 1) * N_F * m_q ** 2 / (2 * np.pi) ** 4) * (
            3 * S0 / (4 * np.pi * (n + 2)) * (1 + 13.0 / 3.0 * ALPHA_S / np.pi)
            + d * np.pi / (4 * S0) * COND_G2 + d * 3 * np.pi / S0 * COND_JS)
    for n in MOMENTS[1]:
        d = 1.0 if n == 0 else 0.0
        out[("P1", n)] = -(S0 ** (n + 1) / (2 * np.pi) ** 4) * 0.5 * (
            -S0 / (2 * np.pi * (n + 2)) * (1 + ALPHA_S / np.pi)
            + d * np.pi / (6 * S0) * COND_G2 + d * 2 * np.pi / S0 * COND_JS)
    return out


def fesr_target_audit() -> dict:
    """Self-check 5a.9: printed (2.56) vs an independent evaluation of (2.50)."""
    pr = printed_targets()
    rows = []
    for key, pv in pr.items():
        ra = recomputed_targets(M_Q)[key]
        rr = recomputed_targets(M_Q_RMS)[key]
        rows.append({"wave": key[0], "n": key[1], "printed_raw": pv,
                     "recomputed_mq_mean": ra, "ratio_mean": pv / ra,
                     "recomputed_mq_rms": rr, "ratio_rms": pv / rr})
    return {"m_q_mean": M_Q, "m_q_mean_MeV": M_Q * M_PI_MEV,
            "m_q_rms": M_Q_RMS, "m_q_rms_MeV": M_Q_RMS * M_PI_MEV,
            "condensate_alphaG2_over_pi": COND_G2, "condensate_jS": COND_JS,
            "rows": rows}


def sr_tolerances(caliber: str, eps_sr: float = EPS_SR) -> dict[tuple[str, int], float]:
    """Half-width of the FESR box for each pre-registered caliber (task section 4)."""
    t = printed_targets()
    if caliber in ("SR-a", "SR-d"):
        # SR-a: per-moment box |M_n - T_n| <= eps_SR.  SR-d: per-wave L2 ball
        # ||(M_n - T_n)_n||_2 <= eps_SR over the two moments of each wave -- the
        # packaging of the authors' released code (norm(wS0) <= eS0) with the
        # 2309 value eps_SR = 2e-3; the per-moment entry here is the ball radius.
        return {k: float(eps_sr) for k in t}
    if caliber == "SR-b":
        return {k: 0.10 * abs(v) for k, v in t.items()}
    if caliber == "SR-c":
        return {k: 0.20 * abs(v) for k, v in t.items()}
    raise ValueError(caliber)


def ff_asymptotic_bounds(M: int, m_q: float = M_Q, eps_ff: float = EPS_FF,
                         frozen_at_s0: bool = True, eps_ff_s0=None, eps_ff_p1=None):
    """(3.75): indices above ``s_0`` and the bound on ``|cF_ell(s_i)|``.

    The paper writes the bound on the rescaled current ``cF`` of (2.33), but the
    numerical implementation explicitly evaluates its factors at each s_i.
    The later at-s0 parenthetical estimates epsilon; it does not establish a
    frozen factor in the inequality. The CLI therefore passes False.
    True is retained only to replay the historical, weaker finite problem.
    That choice was motivated by an old numerical feasibility result which
    cannot establish infeasibility of the corrected source model.

    Returns ``(indices, bound_on_cF, kinematic_factor_to_use)``.
    """
    from .formfactor import kinematic_factor
    idx = np.flatnonzero(~below_s0(M))
    e0 = eps_ff if eps_ff_s0 is None else eps_ff_s0
    e1 = eps_ff if eps_ff_p1 is None else eps_ff_p1
    bound = {0: float(np.sqrt(2.0 * m_q ** 2 * e0)),
             1: float(np.sqrt(0.5 * e1))}
    s = s_nodes(M)
    kin = {ell: (np.full(idx.size, float(kinematic_factor(ell, np.array([S0]))[0]))
                 if frozen_at_s0 else kinematic_factor(ell, s[idx]))
           for ell in (0, 1)}
    return idx, bound, kin


def chiral_ratios(s: float) -> tuple[float, float]:
    """(3.63):  R01 = 3(2s-1)/(s-4),  R21 = 3(2-s)/(s-4)."""
    return 3.0 * (2.0 * s - 1.0) / (s - 4.0), 3.0 * (2.0 - s) / (s - 4.0)


def chiral_reference_point() -> tuple[float, float]:
    """Weinberg-model point (2.13) at the physical f_pi: (f00(3), f11(3))."""
    from .grid import F_PI
    x = (2.0 / np.pi) * (2.0 * 3.0 - 1.0) / (32.0 * np.pi * F_PI ** 2)
    return x, -x / 15.0
