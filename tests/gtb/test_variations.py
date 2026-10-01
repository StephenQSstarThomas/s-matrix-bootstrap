"""App. B variations (Phase B, Task B6): the program follows the kernel export (L, N_B, n_SR) and spec overrides."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")
needs_ks10 = pytest.mark.skipif(not (KS10 / "manifest.json").exists(), reason="M=10 kernel export missing")


@needs_ks10
def test_fesr_blocks_follow_the_number_of_exported_moments():
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks = kernels.load(KS10)

    class Ks4:                                   # an n_SR = 4 export: one more moment row and target per channel
        def __getattr__(self, name):
            if name in ("intS0", "intP1", "intD0"):
                a = getattr(ks, name)
                return np.vstack([a, 2 * a[-1:]])
            if name == "SVZSR":
                s = ks.SVZSR.reshape(3, 3)
                return np.hstack([s, 7 * s[:, -1:]]).ravel()
            return getattr(ks, name)
    k4 = Ks4()
    p = Program(GtbSpec(M=10), k4)
    fesr = [b for b in p.blocks() if b.name.startswith("fesr:")]
    assert [b.size for b in fesr] == [5, 5, 5]
    g = k4.intD0 @ ks.BSpline3[2 * p.n0:3 * p.n0]
    want = (p.c_form(g[3]) - p.const(k4.SVZSR[11])) / p.num(p.spec.eps_sr[2])
    assert np.allclose(fesr[2].forms[(0, 4)], want)


@needs_ks10
def test_program_rejects_a_spec_that_does_not_match_the_export():
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks = kernels.load(KS10)
    for bad in (GtbSpec(M=10, L=24), GtbSpec(M=10, n_b=46), GtbSpec(M=12)):
        with pytest.raises(ValueError):
            Program(bad, ks)


@needs_ks10
def test_cli_spec_follows_the_export_and_the_overrides():
    import argparse
    from smatrix_bootstrap.gtb.cli import spec_from_args
    a = argparse.Namespace(kernels=str(KS10), M=50, ni=44, objective="F1min",
                           set=["eps_chi=0.001", "r=0.02,2,6", "reg_bound=500", "f0_max=0.145743"])
    sp = spec_from_args(a)
    assert (sp.M, sp.L, sp.n_b, sp.f0_ni, sp.objective) == (10, 20, 40, 44, "F1min")
    assert sp.eps_chi == 1e-3 and sp.r == (0.02, 2.0, 6.0) and sp.reg_bound == 500.0 and sp.f0_max == 0.145743
    with pytest.raises(ValueError):
        spec_from_args(argparse.Namespace(**{**vars(a), "set": ["no_such_field=1"]}))
