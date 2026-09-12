"""Self-check 5a.10: the from-scratch operators against independent targets.

The repository's historical Arb implementation (``kernels``, ``operators``,
``model``) is a completely separate code path written for the Newton mainline.
It is used here as a *comparison target only*; no part of the sdp subpackage
imports it outside :mod:`smatrix_bootstrap.sdp.crosscheck`.
"""
import numpy as np
import pytest

from smatrix_bootstrap.sdp import crosscheck as X
from smatrix_bootstrap.sdp import grid
from smatrix_bootstrap.sdp.projector import PartialWaveOperator

M, L = 50, 10
S = grid.s_nodes(M)


def test_no_legacy_import_in_the_subpackage():
    """The rule of task 1b: only crosscheck.py may reach for the old modules."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[2] / "src" / "smatrix_bootstrap" / "sdp"
    banned = ("linear", "scattering", "merit", "ir", "gauge", "quotient", "conic",
              "kernels", "operators", "model", "certificates", "basis", "analytic")
    for f in root.glob("*.py"):
        if f.name == "crosscheck.py":
            continue
        text = f.read_text()
        for name in banned:
            assert f"from ..{name}" not in text and f"import {name}" not in text, (f.name, name)


def test_grid_matches_legacy():
    r = X.compare_grid(M)
    assert r["nodes_max_rel"] < 1e-12
    assert r["weights_max_rel"] < 1e-12


def test_hilbert_kernel_matches_legacy():
    """(3.67) built here vs ``kernels.pv_matrix`` built in Arb."""
    assert X.compare_hilbert_kernel(M) < 1e-13


def test_kinematic_squares_match_legacy():
    """(2.33)^2 vs ``model.current_kinematic_squares``."""
    assert X.compare_kinematic_squares(list(S[:6]) + [9.0, 73.0, 500.0]) < 1e-11


CASES = [(0, 0, 3.0, None), (1, 1, 3.0, None), (2, 0, 1.5, None), (0, 2, 0.5, None),
         (2, 4, 2.0, None), (0, 0, S[10], 10), (1, 1, S[10], 10), (2, 0, S[25], 25),
         (0, 4, S[30], 30), (1, 3, S[5], 5), (0, 0, S[44], 44), (1, 9, S[20], 20),
         (2, 8, S[40], 40), (0, 18, S[35], 35), (1, 19, S[48], 48)]


def test_partial_wave_rows_match_legacy():
    """Every coefficient of 15 rows, after the one documented packing change.

    The historical row uses ``C_flat`` = 2 rho2_{ij} off the diagonal; ours uses
    rho2_{ij}.  With that conversion the two independent implementations agree
    to ~1e-13, which is the accuracy of the legacy side's Arb -> float64 export.
    """
    rep = X.compare_rows(PartialWaveOperator(M), M, L, CASES)
    assert max(r["max_rel_diff"] for r in rep) < 1e-11


def test_rho2_convention_is_the_only_difference():
    """Without the conversion the rho2 off-diagonal is off by exactly 2, and
    nothing else differs -- the difference is explained, not tolerated."""
    op = PartialWaveOperator(M)
    lay = op.lay
    new = op.rows(0, 0, 3.0, None)
    new_c = new[0] + 1j * new[1]
    old = X.old_row(M, L, 0, 0, 3.0, None)
    i, j = lay.triu
    ratio = np.abs(new_c[lay.r2] / old[lay.r2])
    assert np.abs(ratio[i == j] - 1.0).max() < 1e-10
    assert np.abs(ratio[i != j] - 2.0).max() < 1e-10
    for sl in (slice(0, 1), lay.s1, lay.s2, lay.r1):
        m = np.abs(old[sl]) > 0
        assert np.abs(np.abs(new_c[sl][m] / old[sl][m]) - 1.0).max() < 1e-10
