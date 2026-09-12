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


def _uv_control(*, rho_value=3, bits=384, eps_ff=6e-5, frozen=True):
    """S=0, F=1 gives Gram determinant k²(rho_hat-2), independently."""
    from flint import acb
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    verifier = ArbAudit(50, 1, bits)
    return verifier._uv_audit({(ell, i): acb(0) for ell in (0, 1) for i in range(50)},
                             np.zeros((2, 50)), [[rho_value] * 50] * 2,
                             "SR-a", eps_ff, grid.M_Q, frozen_at_s0=frozen)


def test_uv_verifier_gram_positive_control_and_actual_uv_negative_controls():
    """The former UV entry crashes; a PSD Gram alone cannot pass FESR or FF."""
    report = _uv_control()
    assert report["gram_ok"]
    assert report["gram_counts"] == {"certified_pass": 700, "certified_fail": 0,
                                      "inconclusive": 0}
    assert {r["minor"] for r in report["gram"]} == {"1", "2", "3", "12", "13", "23", "123"}
    assert not report["fesr_ok"] and not report["form_factor_ok"]
    assert report["form_factor_counts"]["certified_fail"] == 14
    assert report["n_ff_bounds"] == 14 and report["n_gram_blocks"] == 100


@pytest.mark.parametrize("value,verdict", [("1e-400", "certified_pass"),
                                          ("-1e-400", "certified_fail"),
                                          ("0 +/- 1e-400", "inconclusive"),
                                          ("0", "certified_pass")])
def test_uv_interval_verdict_preserves_tiny_sign_and_zero_overlap(value, verdict):
    """Underflow or 'no proved failure' must never turn a failure into a pass."""
    from flint import arb
    from smatrix_bootstrap.sdp.arbaudit import _constraint_row
    row = _constraint_row(arb(value))
    assert row["verdict"] == verdict
    if verdict == "certified_fail":
        assert arb(row["slack"]["upper"]) < 0
    elif verdict == "inconclusive":
        assert arb(row["slack"]["lower"]) < 0 < arb(row["slack"]["upper"])


@pytest.mark.parametrize("side,verdict", [(1, "certified_pass"),
                                         (-1, "certified_fail"),
                                         (0, "inconclusive")])
def test_uv_gram_near_singular_determinant_controls(side, verdict):
    from flint import arb, ctx
    ctx.prec = 512
    rho = arb(2) + side * arb("1e-100") if side else arb(2, "1e-100")
    report = _uv_control(rho_value=rho, bits=512)
    determinants = [r for r in report["gram"] if r["minor"] == "123"]
    assert len(determinants) == 100
    assert {r["verdict"] for r in determinants} == {verdict}
    assert report["gram_ok"] == (side == 1)


def test_uv_frozen_factor_is_explicit_and_changes_only_ff_checks():
    """With F=1 and eps=.1, frozen k(s0) passes; k(s_i) fails at high P1 s."""
    frozen, per_node = _uv_control(eps_ff=.1), _uv_control(eps_ff=.1, frozen=False)
    assert frozen["form_factor_ok"] and not per_node["form_factor_ok"]
    assert frozen["gram"] == per_node["gram"] and frozen["fesr"] == per_node["fesr"]
    assert frozen["uv_contract"]["m_q_exact"] == "113/2800"


@pytest.mark.parametrize("caliber", ["SR-a", "SR-b", "SR-c"])
def test_uv_exact_printed_targets_and_raw_moments_against_independent_arithmetic(caliber):
    """Catch copied float targets, a missing pi, a wrong k² or a wrong cutoff."""
    from fractions import Fraction as Q
    import mpmath as mp
    from flint import acb, arb
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    verifier = ArbAudit(50, 1)
    report = verifier._uv_audit({(e, i): acb(0) for e in (0, 1) for i in range(50)},
                               [[0] * 50] * 2, [[3] * 50] * 2, caliber, 6e-5, grid.M_Q)
    # Simplified by hand from the four printed equations; no production helper.
    exact_targets = [Q("4.4187e-7") * Q(3600, 49) ** 2,
                     Q("2.82014e-7") * Q(3600, 49) ** 3,
                     Q("5.75484e-5") * Q(3600, 49),
                     Q("2.69948e-5") * Q(3600, 49) ** 2]
    with mp.workdps(150):
        for row, target, ell in zip(report["fesr"], exact_targets, (0, 0, 1, 1)):
            tol = Q(".002") if caliber == "SR-a" else target * Q(".1" if caliber == "SR-b" else ".2")
            for field, value in (("target", target), ("tolerance", tol)):
                assert Q(row[field]["lower"]) <= value <= Q(row[field]["upper"])
                assert Q(row[field]["upper"]) - Q(row[field]["lower"]) < Q("1e-100")
            moment = mp.mpf(0)
            for i in range(43):
                phi = mp.pi * (mp.mpf(i) + mp.mpf(".5")) / 50
                s = 8 / (1 + mp.cos(phi))
                weight = mp.pi / 50 * 8 * mp.sin(phi) / (1 + mp.cos(phi)) ** 2
                # Original unsquared curF prefactors, then square independently.
                k = (mp.sqrt(6 * mp.pi) / (16 * mp.pi ** 3) * s ** (-mp.mpf(1) / 4)
                     * ((s - 4) / 4) ** (mp.mpf(1) / 4)) if ell == 0 else (
                     mp.sqrt(4 * mp.pi / 3) / (8 * mp.pi ** 3) * s ** (-mp.mpf(1) / 4)
                     * ((s - 4) / 4) ** (mp.mpf(3) / 4))
                moment += weight * s ** row["n"] * 3 * k ** 2
            interval = row["moment_interval"]
            assert arb(interval["lower"]) < arb(str(moment)) < arb(interval["upper"])


@pytest.mark.parametrize("side,verdict", [(-1, "certified_pass"),
                                         (1, "certified_fail"),
                                         (0, "inconclusive")])
def test_uv_fesr_boundary_uses_whole_ball(side, verdict):
    """An abs_lower below tolerance is insufficient when the residual overlaps its edge."""
    from fractions import Fraction as Q
    from flint import acb, arb
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    verifier = ArbAudit(50, 1, 512)
    target = arb(str(Q("4.4187e-7") * Q(3600, 49) ** 2))
    desired = target + arb(".002")
    desired += side * arb("1e-100") if side else arb(0, "1e-100")
    s = verifier.x[0]
    k = (6 * arb.pi()).sqrt() / (16 * arb.pi() ** 3) * ((s - 4) / (4 * s)).root(4)
    rho_hat = [[arb(0) for _ in range(50)] for _ in range(2)]
    rho_hat[0][0] = desired / (arb.pi() * verifier.w[0] * k ** 2)
    report = verifier._uv_audit({(e, i): acb(0) for e in (0, 1) for i in range(50)},
                               [[0] * 50] * 2, rho_hat, "SR-a", 6e-5, grid.M_Q)
    assert report["fesr"][0]["verdict"] == verdict


def test_uv_audit_forwards_frozen_setting_and_custom_mass_and_epsilon():
    from fractions import Fraction as Q
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    verifier = ArbAudit(8, 1)
    args = (np.zeros(grid.n_amplitude_vars(8)), np.zeros((2, 8)), np.full((2, 8), 3.))
    frozen = verifier.audit(*args, chi_caliber=None, eps_ff=.005)
    per_node = verifier.audit(*args, chi_caliber=None, eps_ff=.005, ff_frozen_at_s0=False)
    assert all(r["verdict"] == "certified_pass" for r in frozen["form_factor"] if r["ell"] == 1)
    assert any(r["verdict"] == "certified_fail" for r in per_node["form_factor"] if r["ell"] == 1)
    custom = verifier.audit(*args, chi_caliber=None, eps_ff=.25, m_q=.125)
    for row in custom["form_factor"]:
        target = Q(1, 128) if row["ell"] == 0 else Q(1, 8)
        assert Q(row["cap"]["lower"]) <= target <= Q(row["cap"]["upper"])
