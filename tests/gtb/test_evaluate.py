"""Amplitude evaluator at arbitrary real or complex s (Phase B, Task B2)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from flint import acb, arb

KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")
needs = pytest.mark.skipif(not (KS10 / "manifest.json").exists(), reason="M=10 export missing")


@pytest.fixture(scope="module")
def grid():
    from smatrix_bootstrap.gtb.arbkernels import Grid
    return Grid(10, prec_bits=256)


def _close(a, b, tol):
    return abs(complex(a) - complex(b)) <= tol * max(1.0, abs(complex(b)))


def test_interpolant_reproduces_the_conjugate_kernel_at_the_nodes(grid):
    from smatrix_bootstrap.gtb.evaluate import kappa
    for k in (0, 4, 9):
        kap = kappa(grid, ("node", k), "interp")
        for m in range(grid.M):
            want = complex(float(grid.K[k][m].mid()), 1.0 if m == k else 0.0)
            assert _close(kap[m].mid(), want, 1e-30), (k, m)


@needs
def test_node_rows_equal_the_export(grid):
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.evaluate import f_row
    ks = kernels.load(KS10)
    M = grid.M
    for (I, ell, wave_row, k) in ((1, 1, 2 * 10 + 0, 7), (0, 2, 1, 3), (2, 0, 10, 0)):
        n = wave_row * M + k
        def exact_row(a, b):
            flat_a, flat_b = ks.exact(a), ks.exact(b)
            ca = ks.manifest["matrices"][a]["shape"][1]
            cb = ks.manifest["matrices"][b]["shape"][1]
            return [arb(t) for t in flat_a[n * ca:(n + 1) * ca] + flat_b[n * cb:(n + 1) * cb]]
        re, im = exact_row("Re_sig", "Re_rho"), exact_row("Im_sig", "Im_rho")
        got = f_row(grid, I, ell, ("node", k), "interp")
        err = max(float(abs(g * grid.omega[k] - acb(r, i)).mid()) / (1 + float(abs(acb(r, i)).mid()))
                  for g, r, i in zip(got, re, im))
        assert err < 1e-35, (I, ell, k, err)


@needs
def test_midpoint_prescription_at_a_real_point_below_threshold_equals_the_chiral_rows(grid):
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.evaluate import f_row
    ks = kernels.load(KS10)
    got = f_row(grid, 0, 0, ("point", acb(arb(1) / 2)), "midpoint")       # fchi row 0: I=0, l=0, s*=1/2
    want = np.r_[ks.fchi_sig[0], ks.fchi_rho[0]]
    assert max(abs(complex(g.mid()) - w) / (1 + abs(w)) for g, w in zip(got, want)) < 1e-30


@pytest.mark.parametrize("presc", ["interp", "midpoint"])
def test_schwarz_reflection_off_the_cut(grid, presc):
    from smatrix_bootstrap.gtb.evaluate import f_row
    s = acb(30, -5)
    a = f_row(grid, 1, 1, ("point", s), presc)
    b = f_row(grid, 1, 1, ("point", s.conjugate()), presc)
    assert max(abs(complex(x.mid()) - complex(y.mid()).conjugate()) for x, y in zip(a, b)) < 1e-30


def test_constant_amplitude_everywhere(grid):
    from smatrix_bootstrap.gtb.evaluate import f_row
    for I, c in ((0, 2.5), (2, 1.0), (1, 0.0)):
        for where in (("point", acb(25, -3)), ("cut", arb(17)), ("node", 2)):
            ell = 1 if I == 1 else 0
            assert _close(f_row(grid, I, ell, where, "interp")[0].mid(), c if ell == 0 else 0.0, 1e-30)
