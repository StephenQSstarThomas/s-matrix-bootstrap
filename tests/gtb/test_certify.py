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


@needs_leaf10
def test_primal_reports_a_custom_objective_when_the_leaf_has_one(tmp_path):
    import shutil
    from smatrix_bootstrap.gtb.certify import primal
    leaf = tmp_path / "leaf"
    shutil.copytree(LEAF10, leaf, ignore=shutil.ignore_patterns("pmp.json", "sdp", "sdp.ck", "out"))
    v = np.load(leaf / "solution.npz")["v"]
    row = np.zeros(v.size); row[5] = 2.0; row[-1] = -1.0
    np.save(leaf / "objective_row.npy", row)
    rep = primal(leaf, KS10)
    assert abs(float(rep["objective_custom"]["mid"]) - row @ v) < 1e-9 * max(1, abs(row @ v))


HOT = Path("/playpen1/shiqiu/gtb2403_runs/scratch_hot")


@pytest.mark.skipif(not (HOT / "cold46_echi204" / "report.json").exists(), reason="M=10 eps_chi pair missing")
def test_shadow_price_of_eps_chi_matches_two_independent_solves():
    import json
    from smatrix_bootstrap.gtb.sensitivity import shadow_prices
    f1 = lambda d: json.loads((HOT / d / "report.json").read_text())["sdpb_result"]["dualObjective"]
    fd = (f1("cold46_echi204") - f1("cold46")) / 0.04e-3                    # two SDPB solves, eps_chi +2 %
    sp = shadow_prices(HOT / "cold46", KS10, which={"eps_chi"})["eps_chi"]["dopt_dtheta"]
    assert fd > 1 and abs(sp - fd) < 0.05 * fd, (sp, fd)


def _relabelled_leaf(tmp_path, reg):
    import json, shutil
    d = tmp_path / "leaf"
    d.mkdir(parents=True)
    shutil.copy(LEAF10 / "solution.npz", d)
    rep = json.loads((LEAF10 / "report.json").read_text())
    rep["regulariser"] = reg
    (d / "report.json").write_text(json.dumps(rep))
    return d


@needs_leaf10
def test_primal_checks_l4_cuts_and_box_radii(tmp_path):
    from smatrix_bootstrap.gtb.certify import primal
    v = np.load(LEAF10 / "solution.npz")["v"]
    rho, R = v[1 + 20:176], 100.0 * 155                    # M = 10: rho after T0, sigma1, sigma2; R = Mreg Mrho
    assert (rho ** 4).sum() ** 0.25 > R                     # the M = 10 relax leaf is outside the l4 ball
    np.save(tmp_path / "cut.npy", (R * rho / (rho ** 4).sum() ** 0.25)[None])   # its own cut excludes it
    rep = primal(_relabelled_leaf(tmp_path, f"linf-relax+cuts:{tmp_path / 'cut.npy'}"), KS10)
    assert rep["families"]["regulariser"]["certified_fail"] == 1
    np.save(tmp_path / "b.npy", np.full(rho.size, R / rho.size ** 0.25))
    rep = primal(_relabelled_leaf(tmp_path / "x", f"box-inner:{tmp_path / 'b.npy'}"), KS10)
    assert rep["families"]["regulariser"]["certified_fail"] == int((np.abs(rho) > R / rho.size ** 0.25).sum())


@pytest.mark.skipif(not (HOT / "infeas_ni60" / "out").exists(), reason="M=10 infeasible run missing")
def test_infeasibility_is_certified_beyond_the_tip_and_not_inside():
    from smatrix_bootstrap.gtb.dual import certify_infeasible
    bad = certify_infeasible(HOT / "infeas_ni60", KS10)          # F0 = 1.2 x 0.163108, beyond the allowed region
    assert bad["certified_infeasible"] and float(bad["bound"]) < 0
    ok = certify_infeasible(HOT / "cold46", KS10)                # a feasible problem cannot be proven infeasible
    assert not ok["certified_infeasible"]


KS10_L24 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10_L24")


@needs_leaf10
@pytest.mark.skipif(not (KS10_L24 / "manifest.json").exists(), reason="M=10 L=24 export missing")
def test_unitarity_bound_beyond_the_imposed_waves_matches_the_float_scan():
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.certify import unitarity_beyond
    from smatrix_bootstrap.gtb.leafanalysis import omitted_waves
    v = np.load(LEAF10 / "solution.npz")["v"]
    exact = unitarity_beyond(LEAF10, KS10, KS10_L24)
    flt = max(omitted_waves(v[:176], kernels.load(KS10_L24), 20)["eta_max"].values())
    assert exact["L_from"] == 20 and exact["L_to"] == 24 and exact["rows"] == 3 * 2 * 10
    assert abs(exact["max_abs_S_upper"] - flt) < 1e-9 and exact["max_abs_S_upper"] >= flt
