"""The 2103.11484 M-regulariser as SDPB blocks: layout, basis rows, verification."""
import json

import numpy as np
import pytest

from smatrix_bootstrap.sdp.pmp import Pmp
from smatrix_bootstrap.sdp.spec import ModelSpec


def _unit_blocks(path):
    blocks = json.load(open(path))["PositiveMatrixWithPrefactorArray"]
    return [b for b in blocks if len(b["polynomials"]) == 1]


def _coeffs(block):
    return [float(v[0]) for v in block["polynomials"][0][0]]


def test_linf_regulariser_emits_two_unit_blocks_per_density_value(tmp_path):
    w = Pmp(ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=7.5))
    info = w.write(str(tmp_path / "pmp.json"))
    lay = w.ops.lay
    n_rho = (lay.r1.stop - lay.r1.start) + (lay.r2.stop - lay.r2.start)
    ones = _unit_blocks(tmp_path / "pmp.json")
    assert info["block_sizes"][1] == 2 * n_rho == len(ones)
    upper, lower = _coeffs(ones[0]), _coeffs(ones[1])
    assert upper[0] == 1.0 and upper[1 + lay.r1.start] == -1.0 / 7.5
    assert lower[0] == 1.0 and lower[1 + lay.r1.start] == 1.0 / 7.5
    assert sum(x != 0 for x in upper) == 2 and sum(x != 0 for x in lower) == 2
    last = _coeffs(ones[-1])
    assert last[0] == 1.0 and last[1 + lay.r2.stop - 1] == 1.0 / 7.5


def test_linf_regulariser_uses_exact_basis_rows_when_reduced(tmp_path):
    w = Pmp(ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=3.0, reduce_basis=True))
    w.write(str(tmp_path / "pmp.json"))
    lay = w.ops.lay
    row = _coeffs(_unit_blocks(tmp_path / "pmp.json")[0])
    assert row[0] == 1.0
    np.testing.assert_array_equal(np.array(row[1:1 + w.n_a]), -w.basis[lay.r1.start] / 3.0)


def test_precise_assembly_prints_the_same_regulariser_blocks(tmp_path):
    spec = ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=2.0, reduce_basis=True, operator_dps=30)
    w = Pmp(spec, digits=20)
    info = w.write(str(tmp_path / "pmp.json"))
    lay = w.ops.lay
    n_rho = (lay.r1.stop - lay.r1.start) + (lay.r2.stop - lay.r2.start)
    assert info["block_sizes"][1] == 2 * n_rho
    row = _coeffs(_unit_blocks(tmp_path / "pmp.json")[0])
    assert row[0] == 1.0
    np.testing.assert_array_equal(np.array(row[1:1 + w.n_a]), -w.basis[lay.r1.start] / 2.0)


def test_float_verification_reports_norms_and_flags_a_violation():
    w = Pmp(ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=2.0))
    c = np.zeros(w.ops.lay.n)
    c[w.ops.lay.r1.start] = 2.5
    bad = w.verify({"a": c, "c": c})
    assert bad["rho_linf"] == 2.5
    assert bad["regulariser"]["active"] and not bad["constraint_checks"]["regulariser"]
    assert not bad["primal_feasible"]
    z = np.zeros(w.ops.lay.n)
    good = w.verify({"a": z, "c": z})
    assert good["constraint_checks"]["regulariser"] and not good["regulariser"]["active"]
    assert good["primal_feasible"]


def test_precise_verification_flags_a_violation_from_the_full_y_text():
    w = Pmp(ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=2.0, operator_dps=30), digits=20)
    n = w.n_vars - 1
    y = ["0"] * n
    y[w.ops.lay.r1.start] = "3"
    sol = {"y_text": y, "a": np.zeros(n), "c": np.zeros(w.ops.lay.n)}
    sol["c"][w.ops.lay.r1.start] = 3.0
    rep = w.verify(sol)
    assert rep["density_norms"]["rho_linf"] == 3.0
    assert not rep["constraint_checks"]["regulariser"] and not rep["primal_feasible"]


def test_regulariser_spec_validation():
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, reg_norm="linf"))
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, reg_norm="l4", reg_bound=1.0))
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, reg_bound=1.0))
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, reg_norm="linf", reg_bound=-1.0))


def test_cli_defaults_follow_the_authors_conventions():
    from smatrix_bootstrap.sdp import cli
    import argparse
    src = open(cli.__file__).read()
    assert 'default="chi-b"' in src
    assert "default='mixed-pv'" in src
    assert "no density regulator" not in src
