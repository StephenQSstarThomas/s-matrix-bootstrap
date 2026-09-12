"""Form factors, rescaled currents and the 3x3 Gram matrix -- (2.33),(2.41),(2.45)-(2.47),(3.65)-(3.70).

Variables (3.69): ``ImF_{ell,i}`` and the current spectral densities
``rho_{ell,i} >= 0``, ``ell = 0`` (scalar current j_S, I=0) and ``ell = 1``
(vector current j_V, I=1), ``i = 1..M``.

Dispersion relation (2.41) with ``F_ell(0) = 1``, mapped to the disk, gives the
paper's (3.66):  ``ReF_i = 1 + K_ij ImF_j`` with the *same* kernel (3.67)
derived in :mod:`smatrix_bootstrap.sdp.hilbert`.

Rescaled currents (2.33):

    cF_0(s) = sqrt(6 pi)/(16 pi^3)  s^{-1/4} ((s-4)/4)^{1/4} F_0(s)
    cF_1(s) = sqrt(4 pi/3)/(8 pi^3) s^{-1/4} ((s-4)/4)^{3/4} F_1(s)

whose squares reproduce (2.46)-(2.47) exactly -- checked in ``tests/sdp``:

    |cF_0|^2 = (2 pi)^{-4} (3/16 pi) sqrt((s-4)/s) |F_0|^2
    |cF_1|^2 = (2 pi)^{-4} (1/24 pi) (s-4)^{3/2}/sqrt(s) |F_1|^2

Positivity (2.37)/(3.68) is imposed per node as

    [[1, S, cF], [S*, 1, cF*], [cF*, cF, rho]] >= 0 .
"""
from __future__ import annotations

import numpy as np

from .hilbert import hilbert_kernel

def gram_scale(ell, s):
    """Per-node congruence scale for (3.68):  V_i = diag(1, 1, 1/g_i).

    Congruence by an invertible matrix preserves positive semidefiniteness
    exactly, so the feasible set is untouched; only the conditioning changes.
    Taking ``g_i = k_ell(s_i)``, the (2.33) kinematic factor itself, turns the
    block into

        [[1, S, F], [S*, 1, F*], [F*, F, rho / k^2]]

    i.e. the (1,3) entry is the *form factor* rather than the rescaled current.
    F is O(1) at low s and decays at high s, whereas cF_1 = k_1 F_1 grows like
    sqrt(s) and spans four decades across the grid -- which is what makes a
    single constant scale fail.  (The repository's historical
    ``model.current_gram(S, F, rho, k_squared)`` uses the same form.)
    """
    return kinematic_factor(ell, s)


def kinematic_factor(ell: int, s: np.ndarray) -> np.ndarray:
    """The (2.33) prefactor relating ``cF_ell`` to ``F_ell``."""
    s = np.asarray(s, dtype=float)
    if ell == 0:
        return np.sqrt(6.0 * np.pi) / (16.0 * np.pi ** 3) * s ** -0.25 * ((s - 4.0) / 4.0) ** 0.25
    if ell == 1:
        return np.sqrt(4.0 * np.pi / 3.0) / (8.0 * np.pi ** 3) * s ** -0.25 * ((s - 4.0) / 4.0) ** 0.75
    raise ValueError("only the S0 (ell=0) and P1 (ell=1) currents are used")


def kinematic_square_reference(ell: int, s: np.ndarray) -> np.ndarray:
    """Independent form of ``|cF_ell / F_ell|^2`` taken from (2.46)-(2.47)."""
    s = np.asarray(s, dtype=float)
    if ell == 0:
        return (2 * np.pi) ** -4 * (3.0 / (16.0 * np.pi)) * np.sqrt((s - 4.0) / s)
    if ell == 1:
        return (2 * np.pi) ** -4 * (1.0 / (24.0 * np.pi)) * (s - 4.0) ** 1.5 / np.sqrt(s)
    raise ValueError(ell)


def real_part_operator(M: int) -> tuple[np.ndarray, np.ndarray]:
    """``(const, K)`` with ``ReF = const + K @ ImF``, eq. (3.66)."""
    return np.ones(M), hilbert_kernel(M)


def realify(block: np.ndarray) -> np.ndarray:
    """``3x3`` Hermitian -> ``6x6`` real symmetric, PSD-equivalent."""
    return np.block([[block.real, -block.imag], [block.imag, block.real]])


def gram_block(S: complex, cF: complex, rho: float) -> np.ndarray:
    """The paper's (3.68) matrix, as a complex ``3x3``."""
    return np.array([[1.0, S, cF],
                     [np.conj(S), 1.0, np.conj(cF)],
                     [np.conj(cF), cF, rho]], dtype=complex)
