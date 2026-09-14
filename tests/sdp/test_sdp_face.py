"""Degenerate-face diagnostic: the near-optimal slab block and node functionals as the PMP objective."""
import json

import numpy as np
import pytest

from smatrix_bootstrap.sdp import sdpb
from smatrix_bootstrap.sdp.pmp import Pmp
from smatrix_bootstrap.sdp.spec import ModelSpec

FACE = {"direction": [0.0, 1.0], "value": -0.0049, "margin": 2e-6}
FUNC = {"kind": "ImKH", "wave": "P1", "node": 2, "sense": "max"}


def _blocks(path):
    return json.load(open(path))["PositiveMatrixWithPrefactorArray"]


def _coeffs(block):
    return [float(v[0]) for v in block["polynomials"][0][0]]


def test_face_block_follows_the_section_and_the_objective_is_the_node_row(tmp_path):
    w = Pmp(ModelSpec(M=4, L=1), fix_f00=0.07, face=FACE, functional=FUNC)
    w.write(str(tmp_path / "pmp.json"))
    doc = json.load(open(tmp_path / "pmp.json"))
    obj = [float(v) for v in doc["objective"]]
    re_row, im_row = w.gram[1]
    assert obj[0] == 0.0
    np.testing.assert_allclose(obj[1:1 + w.n_a], im_row[2])
    blocks = doc["PositiveMatrixWithPrefactorArray"]
    section = _coeffs(blocks[-2])
    assert section[0] == pytest.approx(-0.07)
    np.testing.assert_allclose(section[1:1 + w.n_a], w.f00)
    face = _coeffs(blocks[-1])
    assert face[0] == pytest.approx(-(FACE["value"] - FACE["margin"]))
    np.testing.assert_allclose(face[1:1 + w.n_a], w.f11)


def test_functional_sense_and_kinds():
    w = Pmp(ModelSpec(M=4, L=1), functional=dict(FUNC, sense="min"))
    np.testing.assert_allclose(w.objective()[1:1 + w.n_a], -w.gram[1][1][2])
    w = Pmp(ModelSpec(M=4, L=1), functional=dict(FUNC, kind="ImS", wave="S0", node=3))
    np.testing.assert_allclose(w.objective()[1:1 + w.n_a], w.gram[0][0][3])
    w = Pmp(ModelSpec(M=4, L=1, uv=True), functional={"kind": "ImF", "wave": "P1", "node": 1, "sense": "max"})
    obj = w.objective()
    assert obj[w.i_ImF + 4 + 1] == 1.0 and np.count_nonzero(obj) == 1
    w = Pmp(ModelSpec(M=4, L=1, uv=True), functional={"kind": "rho", "wave": "S0", "node": 0, "sense": "min"})
    obj = w.objective()
    assert obj[w.i_rho] == -1.0 and np.count_nonzero(obj) == 1


@pytest.mark.parametrize("bad", [
    dict(FUNC, kind="ReS"), dict(FUNC, wave="S2"), dict(FUNC, node=4), dict(FUNC, node=-1),
    dict(FUNC, node=True), dict(FUNC, sense="up"), {"kind": "ImF", "wave": "P1", "node": 1, "sense": "max"}])
def test_functional_validation(bad):
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1), functional=bad)


@pytest.mark.parametrize("bad", [
    dict(FACE, margin=0.0), dict(FACE, margin=-1e-6), dict(FACE, direction=[0.0, 0.0]),
    dict(FACE, direction=[1.0, float("nan")]), dict(FACE, value=float("inf"))])
def test_face_validation(bad):
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1), face=bad)


def test_float_verification_reports_the_functional_and_checks_the_face():
    w = Pmp(ModelSpec(M=4, L=1), face=FACE, functional=FUNC)
    z = np.zeros(w.ops.lay.n)
    rep = w.verify({"a": z, "c": z})
    assert rep["functional"]["value"] == 0.0 and rep["objective_recomputed"] == 0.0
    assert rep["face"]["slack"] == pytest.approx(0.0049 + 2e-6) and rep["constraint_checks"]["face"]
    w = Pmp(ModelSpec(M=4, L=1), face=dict(FACE, value=0.01), functional=dict(FUNC, sense="min"))
    rep = w.verify({"a": z, "c": z})
    assert not rep["constraint_checks"]["face"] and not rep["primal_feasible"]


def test_precise_assembly_uses_the_arb_node_rows_and_verifies_the_face(tmp_path):
    spec = ModelSpec(M=4, L=1, operator_dps=30)
    w = Pmp(spec, face=dict(FACE, value=0.01), functional=FUNC, digits=20)
    w.write(str(tmp_path / "pmp.json"))
    doc = json.load(open(tmp_path / "pmp.json"))
    obj = np.array([float(v) for v in doc["objective"]])
    rows = w.precise.rows(1, 1, node=2) * w.precise.kap[2]
    expected = np.array([float(v) for v in w.precise.projected(rows, None)[1]])
    np.testing.assert_allclose(obj[1:1 + w.n_a], expected, rtol=1e-15, atol=1e-300)
    face = _coeffs(doc["PositiveMatrixWithPrefactorArray"][-1])
    assert face[0] == pytest.approx(-(0.01 - 2e-6))
    y = ["0"] * (w.n_vars - 1)
    rep = w.verify({"y_text": y, "a": np.zeros(w.n_vars - 1), "c": np.zeros(w.ops.lay.n)})
    assert rep["functional"]["value"] == 0.0 and rep["objective_recomputed"] == 0.0
    assert not rep["constraint_checks"]["face"] and not rep["primal_feasible"]


def _source(tmp_path):
    source = {"solver": "SDPB", "status": "numerically_accepted", "accepted": True,
              "verification": {"primal_feasible": True, "f00_3": .0733, "f11_3": -.0046, "objective_recomputed": -.0046},
              "convergence": {"solver_optimal": True}, "spec": ModelSpec().__dict__,
              "settings": sdpb.Settings().__dict__, "direction": [0., 1.], "fix_f00": .0733}
    path = tmp_path / "report.json"
    path.write_text(json.dumps(source))
    return path


def test_face_support_keeps_the_source_direction_and_section(tmp_path, monkeypatch):
    path, seen = _source(tmp_path), []
    monkeypatch.setattr(sdpb, "run_once", lambda *a, **k: seen.append((a, k)) or ({"accepted": False}, None, None))
    sdpb.support_from_saved(path, tmp_path / "out", face_margin=2e-6, functional=FUNC)
    args, kw = seen[0]
    assert list(args[2]) == [0., 1.] and args[3] == .0733
    assert kw["face"] == {"direction": [0., 1.], "value": -.0046, "margin": 2e-6}
    assert kw["functional"] == FUNC and kw["prepared"] == path


def test_functional_only_support_optimises_over_the_whole_feasible_set(tmp_path, monkeypatch):
    path, seen = _source(tmp_path), []
    monkeypatch.setattr(sdpb, "run_once", lambda *a, **k: seen.append((a, k)) or ({"accepted": False}, None, None))
    sdpb.support_from_saved(path, tmp_path / "out", functional=FUNC)
    args, kw = seen[0]
    assert list(args[2]) == [0., 1.] and args[3] == .0733 and kw["face"] is None and kw["functional"] == FUNC


@pytest.mark.parametrize("kw", [
    dict(face_margin=2e-6), dict(face_margin=0., functional=FUNC), dict(functional=FUNC, point="ref"),
    dict(face_margin=2e-6, functional=FUNC, point="tip"), dict(face_margin=2e-6, functional=FUNC, direction=[1., 0.]),
    dict(face_margin=2e-6, functional=FUNC, warm_start=True)])
def test_face_support_rejects_conflicting_requests(tmp_path, monkeypatch, kw):
    path = _source(tmp_path)
    monkeypatch.setattr(sdpb, "run_once", lambda *a, **k: ({"accepted": False}, None, None))
    with pytest.raises(ValueError):
        sdpb.support_from_saved(path, tmp_path / "out", **kw)


def test_run_once_refuses_exact_input_or_warm_start_with_a_face(tmp_path):
    with pytest.raises(ValueError):
        sdpb.run_once(ModelSpec(), tmp_path / "x", prepared=tmp_path / "r.json", exact_input=True, face=FACE, functional=FUNC)
    with pytest.raises(ValueError):
        sdpb.run_once(ModelSpec(), tmp_path / "x", prepared=tmp_path / "r.json", warm_start=True, functional=FUNC)


@pytest.mark.parametrize("argv", [
    ["support", "--source-report", "r", "--out", "o", "--face-margin", "2e-6"],
    ["support", "--source-report", "r", "--out", "o", "--functional", "ImKH", "P1", "38", "max", "--direction", "1", "0"],
    ["support", "--source-report", "r", "--out", "o", "--face-margin", "2e-6", "--functional", "ImKH", "P1", "x", "max"],
    ["support", "--source-report", "r", "--out", "o", "--face-margin", "2e-6", "--functional", "ImKH", "P1", "38", "max", "--point", "tip"],
    ["refine", "--source-report", "r", "--out", "o", "--face-margin", "2e-6", "--functional", "ImKH", "P1", "38", "max"]])
def test_cli_face_arguments_are_validated(argv):
    from smatrix_bootstrap.sdp.__main__ import main
    with pytest.raises(SystemExit):
        main(argv)


def test_cli_face_arguments_reach_the_driver(tmp_path, monkeypatch):
    from smatrix_bootstrap.sdp import __main__ as entry
    seen = []
    monkeypatch.setattr(sdpb, "support_from_saved", lambda *a, **k: seen.append((a, k)) or {"accepted": True})
    assert entry.main(["support", "--source-report", "r", "--out", "o", "--face-margin", "2e-6",
                       "--functional", "ImKH", "P1", "38", "max"]) == 0
    assert seen[0][1]["face_margin"] == 2e-6
    assert seen[0][1]["functional"] == {"kind": "ImKH", "wave": "P1", "node": 38, "sense": "max"}


def test_cli_functional_without_face_reaches_the_driver(tmp_path, monkeypatch):
    from smatrix_bootstrap.sdp import __main__ as entry
    seen = []
    monkeypatch.setattr(sdpb, "support_from_saved", lambda *a, **k: seen.append((a, k)) or {"accepted": True})
    assert entry.main(["support", "--source-report", "r", "--out", "o", "--functional", "SRmom", "S0", "0", "min"]) == 0
    assert "face_margin" not in seen[0][1]
    assert seen[0][1]["functional"] == {"kind": "SRmom", "wave": "S0", "node": 0, "sense": "min"}


# ---------------------------------------------------------------- moment-range diagnostic (SRmom, sr_free)
from smatrix_bootstrap.sdp import constraints as C
from smatrix_bootstrap.sdp import formfactor as FFM


def test_srmom_objective_is_the_fesr_moment_row():
    spec = ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-b")
    w = Pmp(spec, functional={"kind": "SRmom", "wave": "S0", "node": 0, "sense": "min"})
    obj = w.objective()
    mr = C.moment_row(4, 0) * FFM.gram_scale(0, w.ops.s) ** 2
    np.testing.assert_allclose(obj[w.i_rho:w.i_rho + 4], -mr)
    assert np.count_nonzero(obj[:w.i_rho]) == 0 and np.count_nonzero(obj[w.i_rho + 4:]) == 0
    w = Pmp(spec, functional={"kind": "SRmom", "wave": "P1", "node": -1, "sense": "max"})
    obj = w.objective()
    np.testing.assert_allclose(obj[w.i_rho + 4:w.i_rho + 8], C.moment_row(4, -1) * FFM.gram_scale(1, w.ops.s) ** 2)


@pytest.mark.parametrize("bad", [
    {"kind": "SRmom", "wave": "S0", "node": 2, "sense": "min"},
    {"kind": "SRmom", "wave": "P1", "node": 1, "sense": "min"}])
def test_srmom_validation(bad):
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, uv=True), functional=bad)
    with pytest.raises(ValueError):
        Pmp(ModelSpec(M=4, L=1, uv=True, uv_parts=("gram", "ff")), functional={"kind": "SRmom", "wave": "S0", "node": 0, "sense": "min"})


def test_sr_free_drops_exactly_that_moment_box(tmp_path):
    full = Pmp(ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-b"))
    freed = Pmp(ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-b", sr_free=(("S0", 0),)))
    a = full.write(str(tmp_path / "a.json")); b = freed.write(str(tmp_path / "b.json"))
    assert a["n_blocks"] - b["n_blocks"] == 2 and a["block_sizes"][1] - b["block_sizes"][1] == 2
    full_d = Pmp(ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-d")).write(str(tmp_path / "d.json"))
    c = Pmp(ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-d", sr_free=(("P1", -1),))).write(str(tmp_path / "c.json"))
    # SR-d keeps one arrow per wave; the P1 arrow now spans one residual (2x2) instead of two (3x3)
    assert c["n_blocks"] == full_d["n_blocks"]
    assert c["block_sizes"].get(2, 0) == full_d["block_sizes"].get(2, 0) + 1
    assert c["block_sizes"][3] == full_d["block_sizes"][3] - 1


def test_float_verification_ignores_a_freed_moment_and_reports_the_moment():
    spec = ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-b", sr_free=(("S0", 0),))
    w = Pmp(spec, functional={"kind": "SRmom", "wave": "S0", "node": 0, "sense": "min"})
    z = np.zeros(w.ops.lay.n)
    sol = {"a": z, "c": z, "ImF": np.zeros((2, 4)), "rho_hat": np.zeros((2, 4))}
    sol["rho_hat"][0] = 7.0
    rep = w.verify(sol)
    mr = C.moment_row(4, 0) * FFM.gram_scale(0, w.ops.s) ** 2
    assert rep["functional"]["value"] == pytest.approx(float(mr @ sol["rho_hat"][0]))
    assert rep["objective_recomputed"] == pytest.approx(-rep["functional"]["value"])
    rows = {(r["wave"], r["n"]): r for r in rep["fesr"]["rows"]}
    assert rows[("S0", 0)]["free"] and not rows[("S0", 1)]["free"]
    assert rep["fesr"]["max_violation_over_tolerance"] == pytest.approx(
        max(r["violation"] / r["tolerance"] for k, r in rows.items() if k != ("S0", 0)))


def test_precise_verification_marks_the_freed_row_not_imposed():
    spec = ModelSpec(M=4, L=1, uv=True, sr_caliber="SR-b", sr_free=(("S0", 0),), operator_dps=30)
    w = Pmp(spec, functional={"kind": "SRmom", "wave": "S0", "node": 0, "sense": "max"}, digits=20)
    y = ["0"] * (w.n_vars - 1)
    rep = w.verify({"y_text": y, "a": np.zeros(w.n_vars - 1), "c": np.zeros(w.ops.lay.n)})
    rows = {(r["wave"], r["n"]): r for r in rep["fesr"]["rows"]}
    assert rows[("S0", 0)]["verdict"] == "not_imposed" and rows[("S0", 0)]["free"]
    assert rep["functional"]["value"] == 0.0 and rep["objective_recomputed"] == 0.0
    # the freed row's violation must not enter the FESR check; the other rows do (zero rho violates them)
    assert not rep["constraint_checks"]["fesr"]


def test_solve_cli_parses_sr_free_and_functional(tmp_path, monkeypatch):
    from smatrix_bootstrap.sdp import cli
    seen = []
    monkeypatch.setattr(cli, "run_once", lambda spec, *a, **k: seen.append((spec, k)) or ({"accepted": True, "verification": {"f00_3": 0.}}, None, None))
    rc = cli.main(["--workdir", str(tmp_path / "w"), "--M", "4", "--L", "1", "--uv", "--sr", "SR-b", "--sr-free", "S0", "0",
                   "--functional", "SRmom", "S0", "0", "min", "--skip-mma-audit", "--operator-dps", "17"])
    assert rc == 0 and seen[0][0].sr_free == (("S0", 0),)
    assert seen[0][1]["functional"] == {"kind": "SRmom", "wave": "S0", "node": 0, "sense": "min"}
    with pytest.raises(SystemExit):
        cli.main(["--workdir", str(tmp_path / "v"), "--M", "4", "--L", "1", "--sr-free", "S0", "0", "--skip-mma-audit"])
    with pytest.raises(SystemExit):
        cli.main(["--workdir", str(tmp_path / "u"), "--M", "4", "--L", "1", "--uv", "--functional", "SRmom", "S0", "0", "min",
                  "--points", "tip", "--skip-mma-audit"])
