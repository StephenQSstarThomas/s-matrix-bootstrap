"""Discrete Cauchy transforms on the collocation grid -- (2.7), (3.66), (3.67).

Let ``g(nu) = (1/pi) int_4^inf sigma(x) / (x - nu) dx``.  Because ``g`` is
analytic in the nu-plane cut along ``[4, inf)`` and real below the cut, the
conformal map (3.58) turns it into a function ``G(z)`` analytic on the unit
disk with real Taylor coefficients.  On the upper half circle
``z = exp(i phi)``:

    Im G(exp(i phi_i)) = sigma(s_i),
    Re G(exp(i phi_i)) = G(0) + (K sigma)_i ,

where ``K`` is the paper's kernel (3.67) and ``G(0) = g(nu = 0)``.  We *derive*
``K`` here rather than quoting it: writing ``G = sum_n a_n z^n`` with real
``a_n`` makes ``Im G`` odd and ``Re G`` even in ``phi``, so the odd extension of
``sigma`` to the full circle lives on the 2M-point staggered grid
``theta_m = (pi/M)(m - 1/2)``, and the conjugate-function operator there is the
classical odd-offset cotangent kernel.  Collecting the node ``j`` (at
``theta = +phi_j``, value ``+sigma_j``) with its mirror ``1 - j`` (at
``theta = -phi_j``, value ``-sigma_j``) reproduces (3.67) exactly -- see
:func:`hilbert_kernel` and the equality test in ``tests/sdp``.

Off the cut we use the midpoint rule in ``phi`` induced by (3.61), which is the
same rule the paper prescribes for the FESR in (3.72).
"""
from __future__ import annotations

import numpy as np

from .grid import dsdphi, phi_nodes, s_nodes


def hilbert_kernel(M: int) -> np.ndarray:
    """Paper eq. (3.67):  ``K_ij = Kt_{i+j-2M-1} - Kt_{i-j}`` with
    ``Kt_m = 0`` for ``m = 0`` and ``(1-(-1)^m)/(2M) cot(m pi / 2M)`` otherwise.
    """
    i = np.arange(1, M + 1)[:, None]
    j = np.arange(1, M + 1)[None, :]
    return _k_tilde(i + j - 2 * M - 1, M) - _k_tilde(i - j, M)


def _k_tilde(m: np.ndarray, M: int) -> np.ndarray:
    m = np.asarray(m)
    out = np.zeros(m.shape, dtype=float)
    nz = m != 0
    mm = m[nz]
    out[nz] = (1.0 - (-1.0) ** mm) / (2.0 * M) / np.tan(mm * np.pi / (2.0 * M))
    return out


def hilbert_kernel_from_odd_extension(M: int) -> np.ndarray:
    """Independent construction of (3.67) from the 2M-point staggered grid.

    The conjugate-function operator on N = 2M equidistant points keeps only odd
    index offsets, with weight ``(2/N) cot(offset * pi / N)``.  Applying it to
    the odd extension of ``sigma`` (node ``j`` at ``+phi_j``, mirror node
    ``1-j`` at ``-phi_j`` carrying ``-sigma_j``) must reproduce
    :func:`hilbert_kernel`.
    """
    N = 2 * M

    def d(off):
        off = np.asarray(off)
        out = np.zeros(off.shape, dtype=float)
        odd = (off % 2) != 0
        out[odd] = (2.0 / N) / np.tan(off[odd] * np.pi / N)
        return out

    i = np.arange(1, M + 1)[:, None]
    j = np.arange(1, M + 1)[None, :]
    return d(j - i) - d((1 - j) - i)


def cauchy_offcut_row(M: int, nu: float) -> np.ndarray:
    """Row ``c`` with ``(1/pi) int sigma(x)/(x-nu) dx ~ c . sigma`` for ``nu`` off the cut.

    Midpoint rule in phi:  ``(1/pi)(pi/M) sum_j (ds/dphi)_j sigma_j/(s_j - nu)``.
    """
    s, w = s_nodes(M), dsdphi(M)
    if np.any(np.abs(s - nu) < 1e-12):
        raise ValueError(f"nu={nu} collides with a collocation node")
    return (w / (s - nu)) / M


def cauchy_oncut_real(M: int) -> np.ndarray:
    """Matrix ``A`` with ``Re g(s_k + i0) = (A sigma)_k``.

    ``A = K + 1 c(0)^T``: the conjugate-function part (3.67) plus the constant
    ``G(0) = g(nu=0)`` evaluated with the same midpoint rule.
    """
    return hilbert_kernel(M) + cauchy_offcut_row(M, 0.0)[None, :]


def cauchy_oncut_pv_midpoint(M: int) -> np.ndarray:
    """Punctured midpoint rule for the principal value (comparison only).

    Not used by the model: it is the low-order alternative to
    :func:`cauchy_oncut_real`, retained so the two discretisations of the same
    principal value can be compared in the self-check.
    """
    s, w = s_nodes(M), dsdphi(M)
    with np.errstate(divide="ignore"):
        out = (w[None, :] / (s[None, :] - s[:, None])) / M
    np.fill_diagonal(out, 0.0)
    return out


def disk_series_coeffs(M: int) -> np.ndarray:
    """Matrix ``T`` with ``a_n = (T sigma)_n``, ``n = 1..M``, from ``sigma_j = Im G(e^{i phi_j})``.

    Solves ``sigma_j = sum_{n=1}^{M} a_n sin(n phi_j)``; used only by self-checks
    (it gives an independent route to ``Re G`` and to ``G`` at interior points).
    """
    p = phi_nodes(M)
    S = np.sin(np.outer(p, np.arange(1, M + 1)))
    return np.linalg.inv(S)
