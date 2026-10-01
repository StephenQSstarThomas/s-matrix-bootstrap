"""Verification layer: V1 (Arb primal certificate) and V2 (dual bound) on small problems with known answers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from flint import arb

KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")
LEAF10 = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10")


def test_verdict_of_a_ball():
    from smatrix_bootstrap.gtb.certify import verdict
    assert verdict(arb("1e-30", "1e-40")) == "certified_pass"
    assert verdict(arb("-1e-30", "1e-40")) == "certified_fail"
    assert verdict(arb("1e-40", "1e-30")) == "inconclusive"


def test_exact_ball_from_a_40_digit_string_contains_the_value():
    from smatrix_bootstrap.gtb.certify import ball
    b = ball("1234567890123456789012345678901234567890e-40")
    assert b.contains(arb("0.1234567890123456789012345678901234567890"))
    assert b.contains(arb("0.12345678901234567890123456789012345678905"))   # half an ulp away
    assert ball("0") == 0


needs_leaf10 = pytest.mark.skipif(not (LEAF10 / "solution.npz").exists(), reason="M=10 leaf missing")


@needs_leaf10
def test_primal_certifies_the_m10_leaf():
    import json
    from smatrix_bootstrap.gtb.certify import primal
    rep = primal(LEAF10, KS10)
    for fam, r in rep["families"].items():
        assert r["certified_fail"] == 0 and r["inconclusive"] == 0, (fam, r)
    sd = json.loads((LEAF10 / "report.json").read_text())["sdpb_result"]["primalObjective"]
    f1 = rep["objective_F1"]
    assert abs(float(f1["mid"]) - sd) < 1e-6 * max(1.0, abs(sd)) and f1["radius"] < 1e-25


@needs_leaf10
def test_primal_flags_a_violated_gram_block():
    from smatrix_bootstrap.gtb.certify import primal
    rep = primal(LEAF10, KS10, perturb={"ImF": (10 + 0, 50.0)})        # Im F of P1 at node 0 (M = 10) set huge
    assert rep["families"]["gram"]["certified_fail"] >= 1


TOY = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/pmp_roundtrip")


@pytest.mark.skipif(not (TOY / "out/X_matrix_0.txt").exists(), reason="toy SDPB solution missing")
def test_dual_bound_on_the_toy_problem_is_rigorous_and_tight():
    """max y1 s.t. [[1,y1],[y1,1]] >= 0, 0 <= y2 <= 2: optimum 1; the certified bound must be in [1, 1 + 1e-20]."""
    from smatrix_bootstrap.gtb.certify import dual_bound, read_dual_matrices
    from smatrix_bootstrap.gtb.program import Block
    blocks = [Block("disk", 2, {(0, 0): ["1", "0", "0"], (0, 1): ["0", "1", "0"], (1, 1): ["1", "0", "0"]}),
              Block("y2_lo", 1, {(0, 0): ["0", "0", "1"]}),
              Block("y2_hi", 1, {(0, 0): ["2", "0", "-1"]})]
    Z = read_dual_matrices(TOY / "out", [b.size for b in blocks])
    res = dual_bound(blocks, ["0", "1", "0"], Z, tau=[4, 2, 2], ybar=[2, 2])
    ub = res["bound"]
    assert ub >= 1 and ub < 1 + arb("1e-20")
    assert res["residual_max"] < 1e-40


@pytest.mark.skipif(not (KS10 / "manifest.json").exists(), reason="M=10 kernel export missing")
def test_exact_program_blocks_agree_with_the_float_program():
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.certify import ExactKernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    pf = Program(GtbSpec(M=10), kernels.load(KS10))
    pe = Program(GtbSpec(M=10), ExactKernels(kernels.load(KS10)))
    bf, be = pf.blocks(), pe.blocks()
    assert len(bf) == len(be) and pf.n_vars == pe.n_vars
    worst = 0.0
    for x, y in zip(bf, be):
        assert x.name == y.name and set(x.forms) == set(y.forms)
        for key in x.forms:
            ex = np.array([float(v.mid()) for v in y.forms[key]])
            worst = max(worst, float(np.max(np.abs(ex - x.forms[key]) / (1 + np.abs(ex)))))
    assert worst < 1e-12


@needs_leaf10
def test_dual_bound_brackets_the_m10_leaf_from_above():
    from smatrix_bootstrap.gtb.certify import certify_leaf
    res = certify_leaf(LEAF10, KS10)
    lo, hi = res["lower"], res["upper"]           # V1 objective ball (feasible point) and V2 bound
    assert lo["certified"] and hi["certified"]
    assert float(hi["value"]) >= float(lo["value"]) and float(hi["value"]) - float(lo["value"]) < 1e-5


@needs_leaf10
def test_dual_residual_is_driven_below_1e_30():
    from smatrix_bootstrap.gtb.certify import certify_leaf
    res = certify_leaf(LEAF10, KS10)
    assert res["upper"]["residual_history"][-1] < 1e-30


@needs_leaf10
def test_cli_certify_writes_a_certificate():
    import json, subprocess, sys, os
    out = subprocess.run([sys.executable, "-m", "smatrix_bootstrap.run", "gtb", "certify", "--leaf", str(LEAF10),
                          "--kernels", str(KS10)], capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"})
    assert out.returncode == 0, out.stderr[-2000:]
    cert = json.loads((LEAF10 / "certificate.json").read_text())
    assert cert["lower"]["certified"] and cert["upper"]["certified"]
    assert set(cert) >= {"lower", "upper", "regulariser", "l4_feasible", "families"}


@pytest.mark.skipif(not Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10_watson/solution.npz").exists(),
                    reason="M=10 Watson leaf missing")
def test_primal_certifies_the_m10_watson_leaf():
    from smatrix_bootstrap.gtb.certify import primal
    rep = primal(Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10_watson"), KS10)
    for fam, r in rep["families"].items():
        assert r["certified_fail"] == 0 and r["inconclusive"] == 0, (fam, r)
