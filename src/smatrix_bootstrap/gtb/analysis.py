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
