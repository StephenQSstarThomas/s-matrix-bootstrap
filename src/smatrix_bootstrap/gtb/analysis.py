"""Analysis of 2403 leaves (Phase B, Task B3; prescriptions in PREREGISTRATION amendment 1).  Units m_pi = 1."""
from __future__ import annotations

import numpy as np

M_PI_MEV, F_PI_MEV, HBARC = 139.57, 92.0, 197.327        # MeV, MeV, MeV fm


def phases(S):
    """delta (degrees) = 1/2 unwrap(arg S) from threshold, eta = |S| (pre-registered convention)."""
    S = np.asarray(S, dtype=complex)
    return np.degrees(0.5 * np.unwrap(np.angle(np.r_[1.0 + 0j, S]))[1:]), np.abs(S)


def crossing(E, delta_deg, target=90.0):
    """First crossing of delta = target by linear interpolation in E, or None."""
    for i in range(len(E) - 1):
        a, b = delta_deg[i] - target, delta_deg[i + 1] - target
        if a * b <= 0 and delta_deg[i] != delta_deg[i + 1]:
            return float(E[i] + (target - delta_deg[i]) * (E[i + 1] - E[i]) / (delta_deg[i + 1] - delta_deg[i]))
    return None


def threshold_fit(k2, re_f, ell):
    """Re f = k^(2l) (a + b k^2 + c k^4): least squares on the given k^2 points -> (a, b)."""
    k2, y = np.asarray(k2, float), np.asarray(re_f, float) / np.asarray(k2, float) ** ell
    a, b, _ = np.linalg.lstsq(np.vstack([np.ones_like(k2), k2, k2 ** 2]).T, y, rcond=None)[0]
    return float(a), float(b)


def radius_fm2(fprime0):
    """<r^2> = 6 F'(0) (m_pi^-2) in fm^2."""
    return 6.0 * fprime0 * (HBARC / M_PI_MEV) ** 2


def lbar(a02, a22, r2S_fm2, r2V_fm2, fpi=F_PI_MEV / M_PI_MEV):
    """One-loop chiPT, v4 (4.5)-(4.8): D-wave scattering lengths and the S, V radii -> lbar_1, 2, 4, 6."""
    c = 1440 * np.pi ** 3 * fpi ** 4
    x, y = a02 * c + 53 / 8, a22 * c + 103 / 40           # l1 + 4 l2,  l1 + l2
    l2 = (x - y) / 3
    l1 = y - l2
    to_mpi = (M_PI_MEV / HBARC) ** 2                      # fm^2 -> m_pi^-2
    l4 = 16 * np.pi ** 2 * fpi ** 2 * r2S_fm2 * to_mpi / 6 + 13 / 12
    l6 = 96 * np.pi ** 2 * fpi ** 2 * r2V_fm2 * to_mpi / 6 + 1
    return {"l1": float(l1), "l2": float(l2), "l4": float(l4), "l6": float(l6)}


# ------------------------------------------------------------------------------------------ resonance fits
def _bw(s, m, g):
    """v4 (5.2): F = -m^2 / (s - m^2 + i m Gamma theta(s - 4))."""
    s = np.asarray(s, float)
    return -m ** 2 / (s - m ** 2 + 1j * m * g * (s > 4))


def _gs(s, m, g):
    """v4 (5.3) as printed: F = m^2 (1 + d Gamma/m) / ((m^2 - s) - i m Gamma (q/q_m)^3 (m/sqrt s))."""
    s = np.asarray(s, float)
    q, qm = np.sqrt(np.maximum(s - 4, 0) / 4), np.sqrt((m ** 2 - 4) / 4)
    d = 3 / np.pi / qm ** 2 * np.log((m + 2 * qm) / 2) + m / (2 * np.pi * qm) - m / (np.pi * qm ** 3)
    return m ** 2 * (1 + d * g / m) / ((m ** 2 - s) - 1j * m * g * (q / qm) ** 3 * (m / np.sqrt(s)))


FF_FORMS = {"bw": _bw, "gs": _gs}


def _lsq(fun, p0, bounds=(-np.inf, np.inf)):
    from scipy.optimize import least_squares
    return least_squares(fun, p0, bounds=bounds, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000).x


_MG = ([2.01, 1e-4], [40.0, 20.0])                       # physical box for (m, Gamma) in m_pi units


def fit_rho(s, absF, form, guess):
    """Unweighted least squares of |F(s)| with the BW or GS form -> (m, Gamma) in m_pi units."""
    f = FF_FORMS[form]
    return tuple(float(v) for v in _lsq(lambda p: np.abs(f(s, *p)) - absF, guess, _MG))


def two_bw(s, m1, g1, m2, g2, A1, A2):
    """v4 (5.9): F = -A1 m1^2/(s - m1^2 + i m1 G1) + A2 m2^2/(s - m2^2 + i m2 G2) (theta(s-4) on the widths)."""
    s = np.asarray(s, float)
    th = s > 4
    return -A1 * m1 ** 2 / (s - m1 ** 2 + 1j * m1 * g1 * th) + A2 * m2 ** 2 / (s - m2 ** 2 + 1j * m2 * g2 * th)


def fit_two_bw(s, absF2, guess):
    lo, hi = [2.01, 1e-4, 2.01, 1e-4, -10, -10], [40, 20, 40, 20, 10, 10]
    return tuple(float(v) for v in _lsq(lambda p: np.abs(two_bw(s, *p)) ** 2 - absF2, guess, (lo, hi)))


def g_rho_pi_pi(m, g):
    """Gamma = g^2 m/(48 pi) (1 - 4/m^2)^(3/2)  (m_pi = 1)."""
    return float(np.sqrt(48 * np.pi * g / (m * (1 - 4 / m ** 2) ** 1.5)))


def g_f2_pi_pi(m, g, M0=1000 / M_PI_MEV):
    """v4 (5.13) with Gamma_pipi = Gamma(pi0 pi0) + Gamma(pi+ pi-) = 3 Gamma(pi0 pi0) taken as the fitted width."""
    return float(np.sqrt((g / 3) * 960 * np.pi / (m * (m / M0) ** 2 * (1 - 4 / m ** 2) ** 2.5)))


def forward(f_by_wave_isospin):
    """T^I(s,0) = 2 sum_l (2l+1) f^I_l(s) and the s<->u (anti)symmetric combinations of v4 section 3.1."""
    T = {I: 2 * sum((2 * l + 1) * f for l, f in waves.items()) for I, waves in f_by_wave_isospin.items()}
    return {"T0+": (T[1] + T[2]) / 2, "T00": T[0] / 3 + 2 * T[2] / 3, "TIt1": T[0] / 3 + T[1] / 2 - 5 * T[2] / 6}
