"""Legendre functions of the second kind Q_ell(z) for real |z| > 1.

Needed for the closed form of the crossed-channel angular integral derived in
:mod:`smatrix_bootstrap.sdp.projector`:

    int_{-1}^{1} P_ell(mu) / (z - mu) dmu = 2 Q_ell(z)        (Neumann formula)

Q_ell is the *minimal* solution of the Legendre recursion, so the upward
recursion is unstable for large z (where Q_ell ~ z^{-ell-1}).  We therefore use
Miller's backward recursion carried out in mpmath, whose unbounded exponent
range also copes with z ~ 1e8 (which occurs at the first collocation node,
where s - 4 ~ 1e-3).  Results are cached on disk.
"""
from __future__ import annotations

import os
import numpy as np
from mpmath import mp, mpf, atanh

_ELL_SAFETY = 25          # guard digits kept after the forward-recursion loss


def _required_dps(z_abs: float, ell_max: int) -> int:
    """Working precision for a stable *forward* Legendre recursion.

    Q_ell is the minimal solution, P_ell the dominant one, and
    P_ell / Q_ell ~ exp((2 ell + 1) xi) with z = cosh(xi).  A forward recursion
    therefore contaminates Q_ell with the dominant solution at relative size
    ~ eps * exp((2 ell + 1) xi); we simply buy back those digits up front.
    """
    xi = float(np.arccosh(max(z_abs, 1.0)))
    lost = (2 * ell_max + 1) * xi / np.log(10.0)
    return int(_ELL_SAFETY + np.ceil(lost)) + 5


def _q_column_mp(z: mpf, ell_max: int) -> list[mpf]:
    """Q_0..Q_ell_max at a single real z with |z| > 1 by upward recursion.

    Q_0(z) = atanh(1/z)  [= (1/2) log((z+1)/(z-1)), the correct branch for
    z < -1 as well, since there (z+1)/(z-1) > 0],
    Q_1(z) = z Q_0(z) - 1,
    (ell+1) Q_{ell+1} = (2 ell + 1) z Q_ell - ell Q_{ell-1}.
    """
    q0 = atanh(1 / z)
    if ell_max == 0:
        return [q0]
    out = [q0, z * q0 - 1]
    for ell in range(1, ell_max):
        out.append(((2 * ell + 1) * z * out[ell] - ell * out[ell - 1]) / (ell + 1))
    return out


def q_table(z, ell_max: int):
    """``(ell_max+1, len(z))`` array of Q_ell(z_k) as float64.

    Valid for real z with |z| > 1 on either side of the cut [-1, 1].
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    if np.any(np.abs(z) <= 1.0):
        raise ValueError("Q_ell requires |z| > 1 on the collocation grid")
    old = mp.dps
    try:
        out = np.empty((ell_max + 1, z.size), dtype=float)
        for k, zk in enumerate(z):
            mp.dps = _required_dps(abs(float(zk)), ell_max)
            col = _q_column_mp(mpf(float(zk)), ell_max)
            for ell, v in enumerate(col):
                out[ell, k] = float(v)
    finally:
        mp.dps = old
    return out


class QCache:
    """Disk-backed cache of ``q_table`` keyed by (rounded z vector, ell_max)."""

    def __init__(self, path: str | None = None) -> None:
        self.path = path or os.environ.get("SDP_CACHE", "/tmp")
        self._mem: dict[tuple, np.ndarray] = {}

    def get(self, z: np.ndarray, ell_max: int) -> np.ndarray:
        key = (ell_max, hash(np.asarray(z, dtype=float).tobytes()))
        hit = self._mem.get(key)
        if hit is not None:
            return hit
        fn = os.path.join(self.path, f"q_{ell_max}_{key[1] & 0xFFFFFFFFFFFF:012x}.npy")
        if os.path.exists(fn):
            val = np.load(fn)
        else:
            val = q_table(z, ell_max)
            os.makedirs(self.path, exist_ok=True)
            np.save(fn, val)
        self._mem[key] = val
        return val
