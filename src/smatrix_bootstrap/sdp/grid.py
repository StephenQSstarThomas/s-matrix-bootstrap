"""Conformal map, collocation grid and kinematic factors -- paper (3.58)-(3.62).

All quantities are in units m_pi = 1.  Derived here from scratch:

* ``z(nu) = (2 - sqrt(4-nu)) / (2 + sqrt(4-nu))``                    -- (3.58)
* ``z(nu + i0) = exp(i phi)`` for ``nu > 4``                         -- (3.59)
* ``nu(phi) = 8 / (1 + cos phi)``                                    -- (3.60)
* ``phi_i = (pi/M) (i - 1/2)``, ``i = 1..M``                         -- (3.61)
* ``ds/dphi = 8 sin phi / (1 + cos phi)^2``                          -- (3.71)

The quadrature rule induced by (3.61) is the *midpoint rule in phi*, which is
exactly what the paper uses for the FESR discretisation (3.72):

    int_4^inf dx g(x)  ->  (pi/M) sum_i (ds/dphi)_i g(s_i)
"""
from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------- QCD inputs
M_PI_MEV = 140.0          # (2.2): m_pi = 140 MeV sets the unit
F_PI = 92.0 / 140.0       # f_pi = 92 MeV in units m_pi = 1
S0_GEV2 = 1.2 ** 2        # (2.53): s_0 = (1.2 GeV)^2
S0 = S0_GEV2 * 1.0e6 / (M_PI_MEV ** 2)   # s_0 in units m_pi^2
ALPHA_S = 0.4             # (2.53)
M_U = 4.0 / M_PI_MEV      # (2.53): m_u = 4 MeV
M_D = 7.3 / M_PI_MEV      # (2.53): m_d = 7.3 MeV
M_Q = 0.5 * (M_U + M_D)   # average light quark mass
N_C, N_F = 3, 2


def z_map(nu: np.ndarray | complex) -> np.ndarray | complex:
    """Conformal map (3.58) of the cut nu-plane onto the unit disk."""
    r = np.sqrt(4.0 - np.asarray(nu, dtype=complex))
    return (2.0 - r) / (2.0 + r)


def phi_nodes(M: int) -> np.ndarray:
    """Collocation angles (3.61)."""
    return (np.pi / M) * (np.arange(1, M + 1) - 0.5)


def s_nodes(M: int) -> np.ndarray:
    """Collocation energies s_i = nu(phi_i), eq. (3.60)+(3.61)."""
    return 8.0 / (1.0 + np.cos(phi_nodes(M)))


def dsdphi(M: int) -> np.ndarray:
    """Jacobian (3.71) evaluated on the collocation grid."""
    p = phi_nodes(M)
    return 8.0 * np.sin(p) / (1.0 + np.cos(p)) ** 2


def quad_weights(M: int) -> np.ndarray:
    """Weights w_i with  int_4^inf dx g(x) ~ sum_i w_i g(s_i)  (midpoint in phi)."""
    return (np.pi / M) * dsdphi(M)


def kappa(s: np.ndarray) -> np.ndarray:
    """Phase-space factor of (2.11):  S = 1 + i kappa f,  kappa = pi sqrt((s-4)/s)."""
    s = np.asarray(s, dtype=float)
    return np.pi * np.sqrt((s - 4.0) / s)


def lambda_rescale(s: np.ndarray, ell: int, floor: float = 0.0) -> np.ndarray:
    """Centrifugal factor Lambda_ell(s) = ((sqrt(s)-2)/(sqrt(s)+2))^(ell/2).

    The exact value spans an enormous range -- Lambda_19(s_1)^2 ~ 1e-80 at M=50,
    because s_1 - 4 ~ 1e-3 -- but it stays far above the double-precision
    underflow threshold (min normal ~2.2e-308), so ``floor`` defaults to 0 and
    the value returned is the true one.  Pass a positive ``floor`` only where a
    *division* by Lambda would otherwise be unbounded; do not use a floored
    Lambda to report the analytic row scale, which is what it means.

    This removes the angular momentum threshold order, not all node/column
    normalization differences. Without clipping it leaves 6.24 decades in
    the M50 cone rows, versus 1.25 for rownorm. Both are valid congruences.
    """
    s = np.asarray(s, dtype=float)
    r = (np.sqrt(s) - 2.0) / (np.sqrt(s) + 2.0)
    with np.errstate(under="ignore"):
        lam = r ** (0.5 * ell)
    return lam if floor <= 0.0 else np.maximum(lam, floor)


def n_amplitude_vars(M: int) -> int:
    """1 (T0) + 2M (sigma_1,2) + M^2 (rho_1) + M(M+1)/2 (rho_2 symmetric)."""
    return 1 + 2 * M + M * M + M * (M + 1) // 2


def sym_pack_index(M: int) -> tuple[np.ndarray, np.ndarray]:
    """Row/col indices of the upper triangle (i<=j) used to pack rho_2."""
    return np.triu_indices(M)
