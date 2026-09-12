"""Task 5a.10: compare the from-scratch operators against independent targets.

Nothing here is an implementation source.  The repository's historical
``kernels``/``operators``/``model`` modules (Arb arithmetic, a completely
separate code path written for the Newton mainline) are loaded *only* here, so
that the new rows can be diffed against them row by row.

Packing note (the one documented convention difference).  The historical row is
indexed by ``kernels.density_labels`` -- ``T0``, ``sigma1_i``, ``sigma2_i``,
``rho1_{ij}`` (full M x M), ``rho2_{ij}`` for ``i <= j`` -- which is the same
ordering as :class:`smatrix_bootstrap.sdp.projector.Layout`.  It then applies
``kernels.density_row_to_cflat``, halving every off-diagonal ``rho2``
coefficient, because its free variable is "C_flat" = 2 rho2_{ij} for i != j
(see ``kernels.coefficient_blocks``, which divides by two on the way back).
Our packed variable is rho2_{ij} itself, so our coefficient is the sum over the
symmetric pair.  The two describe the same bilinear form; comparison therefore
applies the same halving to the new row.
"""
from __future__ import annotations

import numpy as np


def old_row(M: int, L: int, isospin: int, ell: int, s: float, node: int | None, bits: int = 384):
    """One partial-wave row from the historical Arb implementation."""
    src = _src(M, L, bits)
    if node is not None:
        s = float(src.x[node].str(25, radius=False))   # its own float64 node value
    return np.array(_to_complex(src.row(s, ell, isospin, node=node)))


def to_cflat(row: np.ndarray, lay) -> np.ndarray:
    """Convert a Layout-packed row to the historical C_flat convention."""
    out = row.copy()
    i, j = lay.triu
    blk = out[lay.r2].copy()
    blk[i != j] *= 0.5
    out[lay.r2] = blk
    return out


_SRC: dict = {}


def _src(M, L, bits):
    if (M, L, bits) not in _SRC:
        from ..kernels import PVSourceRows
        _SRC[(M, L, bits)] = PVSourceRows(M=M, L=L, bits=bits, subtracted=False)
    return _SRC[(M, L, bits)]


def _to_complex(row):
    from flint import acb
    out = []
    for v in row:
        a = acb(v)
        out.append(complex(float(a.real.str(20, radius=False)),
                           float(a.imag.str(20, radius=False))))
    return out


def compare_rows(new_op, M: int, L: int, cases) -> list[dict]:
    """Row-by-row diff for a list of ``(isospin, ell, s, node)`` cases."""
    out = []
    for isospin, ell, s, node in cases:
        new = new_op.rows(isospin, ell, s, node)
        new_c = new[0] + 1j * new[1]
        new_c = to_cflat(new_c, new_op.lay)
        old = old_row(M, L, isospin, ell, s, node)
        scale = max(np.abs(old).max(), np.abs(new_c).max(), 1e-300)
        out.append({"isospin": isospin, "ell": ell, "s": float(s), "node": node,
                    "max_abs_diff": float(np.abs(new_c - old).max()),
                    "row_scale": float(scale),
                    "max_rel_diff": float(np.abs(new_c - old).max() / scale)})
    return out


def compare_hilbert_kernel(M: int) -> float:
    """New (3.67) vs the historical ``kernels.pv_matrix`` (Arb)."""
    from ..kernels import pv_matrix
    from .hilbert import hilbert_kernel
    old = pv_matrix(M)
    o = np.array([[float(old[i, j].str(20, radius=False)) for j in range(M)] for i in range(M)])
    return float(np.abs(hilbert_kernel(M) - o).max())


def compare_kinematic_squares(s_values) -> float:
    """New (2.33)^2 vs the historical ``model.current_kinematic_squares``."""
    from ..model import current_kinematic_squares
    from .formfactor import kinematic_factor
    worst = 0.0
    for s in s_values:
        old = current_kinematic_squares(float(s))
        for ell in (0, 1):
            o = float(old[ell].str(25, radius=False))
            n = float(kinematic_factor(ell, np.array([float(s)]))[0] ** 2)
            worst = max(worst, abs(n / o - 1.0))
    return worst


def compare_grid(M: int) -> dict:
    """New (3.60)-(3.61) grid and weights vs ``operators.midpoint_grid`` (Arb)."""
    from ..operators import midpoint_grid
    from .grid import dsdphi, s_nodes
    x, w = midpoint_grid(M)
    xo = np.array([float(v.str(25, radius=False)) for v in x])
    wo = np.array([float(v.str(25, radius=False)) for v in w])
    return {"nodes_max_rel": float(np.abs(s_nodes(M) / xo - 1).max()),
            "weights_max_rel": float(np.abs((dsdphi(M) / M) / wo - 1).max())}
