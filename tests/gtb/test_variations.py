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


def test_sdpb_arguments_for_a_cold_and_a_hot_start():
    from smatrix_bootstrap.gtb.pmpwriter import sdpb_args
    cold = sdpb_args(precision=128, timeout=60, duality_gap=1e-6, primal_error=1e-10, dual_error=1e-10,
                     initial_scale=1e4, max_iterations=500, checkpoint_interval=300, write_solution="x,y")
    hot = sdpb_args(precision=128, timeout=60, duality_gap=1e-6, primal_error=1e-10, dual_error=1e-10,
                    initial_scale=1e4, max_iterations=500, checkpoint_interval=300, write_solution="x,y",
                    initial_checkpoint="init.ck")
    assert "--initialMatrixScalePrimal=10000.0" in cold and "-i" not in cold
    assert hot[hot.index("-i") + 1] == "init.ck" and hot[hot.index("-c") + 1] == "sdp.ck"
    assert not any(a.startswith("--initialMatrixScale") for a in hot)


@needs_ks10
def test_hot_start_refuses_a_seed_with_another_structure(tmp_path):
    import json
    from smatrix_bootstrap.gtb.solve import run_leaf
    from smatrix_bootstrap.gtb.spec import GtbSpec
    seed = tmp_path / "seed"
    (seed / "sdp.ck").mkdir(parents=True)
    (seed / "report.json").write_text(json.dumps({"n_vars": 7, "reduce_basis": False, "regulariser": "linf-relax"}))
    with pytest.raises(ValueError, match="hot start"):
        run_leaf(GtbSpec(M=10), KS10, tmp_path / "leaf", hot_start=seed)


@needs_ks10
def test_box_inner_regulariser_uses_per_value_radii_inside_the_l4_ball(tmp_path):
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks, R, n = kernels.load(KS10), 100.0 * 155, 155
    b = np.linspace(1.0, 2.0, n)
    b *= R / (b ** 4).sum() ** 0.25                              # on the l4 sphere: allowed
    np.save(tmp_path / "b.npy", b)
    p = Program(GtbSpec(M=10), ks, reg=f"box-inner:{tmp_path / 'b.npy'}")
    hi = [blk for blk in p.blocks() if blk.name.startswith("reg:") and blk.name.endswith(":hi")]
    assert len(hi) == n and np.isclose(hi[7].forms[(0, 0)][1 + p.groups["rho"].start - 1 + 7], -1 / b[7])
    np.save(tmp_path / "big.npy", b * 1.001)                     # outside the l4 ball: refused
    with pytest.raises(ValueError):
        Program(GtbSpec(M=10), ks, reg=f"box-inner:{tmp_path / 'big.npy'}")


@needs_ks10
def test_l4_cuts_are_supporting_hyperplanes(tmp_path):
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks, R, n = kernels.load(KS10), 100.0 * 155, 155
    rb = np.sin(np.arange(n) + 1.0)
    rb *= R / (rb ** 4).sum() ** 0.25
    np.save(tmp_path / "cuts.npy", rb[None])
    p = Program(GtbSpec(M=10), ks, reg=f"linf-relax+cuts:{tmp_path / 'cuts.npy'}")
    cut = [blk for blk in p.blocks() if blk.name.startswith("reg:cut:")]
    assert len(cut) == 1
    f = np.asarray(cut[0].forms[(0, 0)])
    v = np.zeros(p.n_raw); v[p.groups["rho"]] = rb
    assert abs(f @ np.r_[1.0, p.unlift(v)]) < 1e-9               # rho on the sphere at the cut point: equality
    assert abs(f[0] - 1.0) < 1e-12                               # rho = 0 (and every other amplitude 0): value 1
    np.save(tmp_path / "off.npy", 1.01 * rb[None])
    with pytest.raises(ValueError):
        Program(GtbSpec(M=10), ks, reg=f"linf-relax+cuts:{tmp_path / 'off.npy'}")


@needs_ks10
def test_hot_started_reduced_solve_reuses_the_seed_basis(tmp_path):
    import json
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.solve import leaf_basis
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks = kernels.load(KS10)
    fresh, _ = leaf_basis(ks, GtbSpec(M=10), True, 1e-12, None)
    seed = tmp_path / "seed"
    seed.mkdir()
    np.save(seed / "basis.npy", -fresh)                        # same span, other signs (as another BLAS gives)
    (seed / "report.json").write_text(json.dumps({"basis": {"rank": fresh.shape[1], "from": "test"}}))
    got, info = leaf_basis(ks, GtbSpec(M=10), True, 1e-12, seed)
    assert np.array_equal(got, -fresh) and info["reused_from"] == str(seed)
    assert leaf_basis(ks, GtbSpec(M=10), False, 1e-12, seed) == (None, None)


LEAF10 = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10")


@needs_ks10
def test_midpoint_export_rows_and_scan_agree_with_the_between_node_scan(tmp_path):
    from smatrix_bootstrap.gtb import midpoints
    from smatrix_bootstrap.gtb.arbkernels import Grid
    from smatrix_bootstrap.gtb.leafanalysis import _load, between_nodes
    midpoints.build(KS10, tmp_path / "mid", l_max=4)
    mid = midpoints.load(tmp_path / "mid")
    assert mid["Re"].shape == (3 * 2 * 9, 176) and mid["meta"].shape == (54, 3)
    _, _, prog, v = _load(LEAF10, KS10)
    amp = v[prog.groups["amp"]]
    S = midpoints.scan(mid, amp)
    ref = between_nodes(Grid(10, prec_bits=192), amp, 10)
    for w, (I, ell) in {"S0": (0, 0), "P1": (1, 1), "D2": (2, 2)}.items():
        rows = [n for n, (i, l, k) in enumerate(mid["meta"]) if (i, l) == (I, ell)]
        assert len(rows) == 9 and np.allclose(S[rows][:len(ref[w]["abs_S"])], ref[w]["abs_S"], atol=1e-9)


@needs_ks10
def test_program_adds_midpoint_unitarity_blocks(tmp_path):
    from smatrix_bootstrap.gtb import kernels, midpoints
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    midpoints.build(KS10, tmp_path / "mid", l_max=2)
    mid = midpoints.load(tmp_path / "mid", rows=[3, 7])
    p = Program(GtbSpec(M=10), kernels.load(KS10), extra_unitarity=mid)
    blk = [b for b in p.blocks() if b.name.startswith("mid:")]
    assert len(blk) == 2 and all(b.size == 3 for b in blk)
    yfree = np.random.default_rng(1).normal(scale=1e-3, size=p.n_vars - 1)
    y, amp = np.r_[1.0, yfree], p.lift(yfree)[p.groups["amp"]]
    h = mid["Re"][1] @ amp + 1j * (mid["Im"][1] @ amp)
    lam = mid["Lam"][1]
    assert np.isclose(np.asarray(blk[1].forms[(0, 0)]) @ y, 2 * h.imag / lam ** 2, rtol=1e-12, atol=1e-14)
    assert np.isclose(np.asarray(blk[1].forms[(0, 1)]) @ y, h.real / lam, rtol=1e-12, atol=1e-14)
