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
    import ast,pathlib
    root = pathlib.Path(__file__).resolve().parents[2] / "src" / "smatrix_bootstrap" / "sdp"
    assert list(root.glob('*.py')), 'The import guard must inspect real source files'
    banned = ("linear", "scattering", "merit", "ir", "gauge", "quotient", "conic",
              "kernels", "operators", "model", "certificates", "basis", "analytic")
    for f in root.glob("*.py"):
        if f.name == "crosscheck.py":
            continue
        for node in ast.walk(ast.parse(f.read_text())):
            imports=[]
            if isinstance(node,ast.Import):imports=[v.name for v in node.names]
            if isinstance(node,ast.ImportFrom) and node.level in (0,2):
                imports=[node.module] if node.module else [v.name for v in node.names]
            roots={v.removeprefix('smatrix_bootstrap.').split('.')[0] for v in imports}
            assert not roots.intersection(banned),(f.name,roots.intersection(banned))
def test_cli_default_follows_source_tolerances(tmp_path, monkeypatch):
    from smatrix_bootstrap.sdp import cli
    seen = []
    def capture(spec, *args, **kwargs):
        seen.append(spec)
        return {'accepted': False}, None, None
    monkeypatch.setattr(cli, 'run_once', capture)
    cli.main(['--workdir', str(tmp_path), '--skip-mma-audit', '--chiral', '--uv'])
    assert (seen[0].chi_caliber, seen[0].sr_caliber, seen[0].ff_frozen_at_s0) == ('chi-c', 'SR-a', False)
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
def _uv_control(*, rho_value=3, bits=384, eps_ff=6e-5, frozen=True, S=0):
    """S=0, F=1 gives Gram determinant k²(rho_hat-2), independently."""
    from flint import acb
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    verifier = ArbAudit(50, 1, bits)
    return verifier._uv_audit({(ell, i): acb(S) for ell in (0, 1) for i in range(50)},
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
    from flint import acb, arb
    zero = _uv_control(S=acb(arb('1 +/- 1e-50'), arb('0 +/- 1e-50')))
    assert zero['gram_counts']['inconclusive'] > 0
    assert zero['gram_counts']['certified_fail'] == 0
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
@pytest.mark.parametrize("part", ["gram", "fesr", "ff"])
def test_sdpb_rejects_zero_current_with_nonzero_charge(part):
    """F(0)=1, rho=0 cannot satisfy Gram/FESR, and F=1 violates the FF cap."""
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.problem import ModelSpec
    w = Pmp(ModelSpec(M=8, L=2, chiral=True, uv=True, uv_parts=(part,)))
    sol = {"a": np.zeros(w.n_a), "c": np.zeros(w.ops.lay.n),
           "ImF": np.zeros((2, 8)), "rho_hat": np.zeros((2, 8))}
    report = w.verify(sol)
    assert report["unitarity"]["feasible"]
    assert not report["certified"]
    assert not report["primal_feasible"]
    assert not report["constraint_checks"][part]
def test_sdpb_gram_linear_functions_survive_disk_reduction():
    """A Gram replaces a disk, so removing that disk must retain its amplitude rows."""
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.problem import ModelSpec
    w = Pmp(ModelSpec(M=12, L=3, uv=True, reduce_basis=True,
                      disk_mask=np.zeros(3 * 3 * 12, dtype=bool)))
    for pair in w.ops.gram_rows.values():
        for rows in pair:
            rows = rows / np.maximum(np.max(abs(rows), axis=1, keepdims=True), 1e-300)
            assert np.max(abs(rows - (rows @ w.basis) @ w.basis.T)) < 1e-10
def test_sdpb_feasible_point_alone_is_not_an_optimum():
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.problem import ModelSpec
    w = Pmp(ModelSpec(M=8, L=2))
    report = w.verify({"a": np.zeros(w.n_a), "c": np.zeros(w.ops.lay.n)})
    assert report["primal_feasible"]
    assert not report["certified"]
def test_sdpb_import_has_no_cvxpy_or_legacy_solver_dependency():
    import subprocess
    import sys
    code = '''
import sys
class BlockSolver:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'cvxpy', 'mosek', 'clarabel'}:
            raise ImportError('forbidden solver dependency: ' + fullname)
sys.meta_path.insert(0, BlockSolver())
from smatrix_bootstrap.sdp.pmp import Pmp
from smatrix_bootstrap.sdp.spec import ModelSpec
Pmp(ModelSpec(M=4, L=1))
'''
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
@pytest.mark.parametrize("reason,p,d,pe,expected", [
    ("found primal-dual optimal solution", .08, .0799999, 1e-12, True),
    ("maxIterations exceeded", .08, .0799999, 1e-12, False),
    ("found primal-dual optimal solution", .08, .0799999, 1., False),
    ("found primal-dual optimal solution", .08, .01, 1e-12, False)])
def test_sdpb_status_requires_residuals_gap_and_terminal_reason(reason, p, d, pe, expected):
    from smatrix_bootstrap.sdp.sdpb import convergence_report
    data = dict(terminateReason=reason, primalObjective=p, dualObjective=d,
                primalError=pe, dualError=1e-12, dualityGap=abs(p-d))
    report = convergence_report(data, returncode=0, objective=d)
    assert report["solver_optimal"] is expected
    assert not report["independent_dual_certified"]
def test_sdpb_solution_objective_is_dual_not_primal():
    from smatrix_bootstrap.sdp.sdpb import convergence_report
    data = dict(terminateReason="found primal-dual optimal solution",
                primalObjective=.0800005, dualObjective=.08,
                primalError=1e-12, dualError=1e-12, dualityGap=5e-7)
    assert convergence_report(data, returncode=0, objective=.08)["solver_optimal"]
    assert not convergence_report(data, returncode=0, objective=.1)["solver_optimal"]
def test_sdpb_preprocess_failure_is_retained_by_generation(tmp_path, monkeypatch):
    from smatrix_bootstrap.sdp import sdpb
    from smatrix_bootstrap.sdp.spec import ModelSpec
    failed = {"accepted": False, "status": "preprocess_failed"}
    monkeypatch.setattr(sdpb, "run_once", lambda *args,**kwargs: (failed, None, None))
    assert sdpb.generate(ModelSpec(M=2, L=1), tmp_path) == failed
@pytest.mark.parametrize("args", [["--M", "2", "--out", "/unused", "solve"], ["audit-uv", "--out", "/unused"]])
def test_sdp_cli_cannot_reach_other_solvers(args, monkeypatch):
    from smatrix_bootstrap.sdp.__main__ import main
    import smatrix_bootstrap.sdp.runner as runner
    def forbidden(*args, **kwargs):
        pytest.fail("SDPB CLI dispatched another solver")
    monkeypatch.setattr(runner, "run_job", forbidden)
    with pytest.raises(SystemExit):
        main(args)
def test_sdpb_figure_loader_uses_complete_numerical_acceptance(tmp_path):
    import json
    from smatrix_bootstrap.sdp.figures import load_reports, boundary_points, c1_table
    doc = {"solver": "SDPB", "spec": {"M": 8}, "accepted": True,
           "direction": [1., 0.], "verification": {"f00_3": 1., "f11_3": .1,
               "objective_recomputed": 1., "primal_feasible": True,
               "unitarity": {"feasible": True}},
           "convergence": {"solver_optimal": True}}
    (tmp_path / "report.json").write_text(json.dumps(doc))
    rows = load_reports(str(tmp_path))
    assert len(rows) == 1
    assert boundary_points(rows).shape == (1, 2)
    assert c1_table(rows * 24)["n_directions"] == 1
    assert not c1_table(rows * 24)["pass"]
    doc["accepted"] = False
    (tmp_path / "report.json").write_text(json.dumps(doc))
    assert boundary_points(load_reports(str(tmp_path))).shape == (0, 2)
def test_sdpb_timeout_reaps_only_its_own_process(tmp_path):
    import subprocess
    import sys
    from smatrix_bootstrap.sdp.sdpb import wait_process, stop_owned_process
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"],
                            start_new_session=True)
    try:
        assert wait_process(proc, .01) == 124
    finally:
        with (tmp_path / "log").open("w") as log:
            stop_owned_process(proc, tmp_path / "absent.cid", log)
    assert proc.poll() is not None
def test_reuse_rejects_unknown_operator_producer(tmp_path):
    import json
    from smatrix_bootstrap.sdp.pmp import Pmp
    path = tmp_path/'report.json'
    path.write_text(json.dumps({'solver':'SDPB','spec':{'operator_dps':40}}))
    with pytest.raises(ValueError,match='Source operator changed'):
        Pmp.from_saved(path,(1,0))
@pytest.mark.parametrize('bad,source_L,target_L',[(b,2,2) for b in
    (None,'hash','M','family','layout','packing','producer','finite','unaccepted','no_reduce')]+[(None,12,10),(None,12,13),(None,10,12)])
def test_saved_basis_transfer_is_exact_and_rejects_wrong_identity(tmp_path,bad,source_L,target_L,monkeypatch):
    import hashlib,json
    from dataclasses import asdict,replace
    from pathlib import Path
    from smatrix_bootstrap.sdp import assembly
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.spec import ModelSpec
    spec=ModelSpec(M=4,L=source_L,uv=True,reduce_basis=True,operator_dps=40,scattering_prescription='sine-cardinal')
    source=Pmp(spec,digits=30);basis=source.basis.copy();np.save(tmp_path/'basis.npy',basis)
    record={'solver':'SDPB','status':'numerically_accepted','accepted':True,'spec':asdict(spec),
            'verification':{'primal_feasible':True,'source_audit':{'scattering_prescription':'sine-cardinal'}},
            'convergence':{'solver_optimal':True},'pmp':{'n_a':source.n_a,'n_vars':source.n_vars,'reduce_basis':True},
            'basis':{'sha256':hashlib.sha256((tmp_path/'basis.npy').read_bytes()).hexdigest(),'tolerance':spec.basis_tol},
            'source_sha256':{'src/smatrix_bootstrap/sdp/'+n:hashlib.sha256(Path(assembly.__file__).with_name(n).read_bytes()).hexdigest()
                             for n in ('projector.py','grid.py')}}
    if bad=='hash':record['basis']['sha256']='0'*64
    if bad=='M':record['spec']['M']=5
    if bad=='family':record['spec']['scattering_prescription']='mixed-pv'
    if bad=='layout':record['pmp']['n_a']+=1
    if bad=='packing':record['basis']['packing']='wrong-packing'
    if bad=='producer':record['source_sha256']['src/smatrix_bootstrap/sdp/projector.py']='0'*64
    if bad=='finite':
        basis[0,0]=np.nan;np.save(tmp_path/'basis.npy',basis)
        record['basis']['sha256']=hashlib.sha256((tmp_path/'basis.npy').read_bytes()).hexdigest()
    if bad=='unaccepted':record['accepted']=False
    path=tmp_path/'report.json';path.write_text(json.dumps(record))
    target=replace(spec,L=target_L,uv=False,cone_scaling='none',basis_tol=.5,reduce_basis=bad!='no_reduce')
    if bad:
        with pytest.raises(ValueError):Pmp(target,digits=30,basis_source_report=path)
    else:
        w=Pmp(target,digits=30,basis_source_report=path)
        assert np.array_equal(w.basis,basis) and w.n_a==source.n_a and w.i_ImF is None
        assert w.basis_source['source_basis_sha256']==record['basis']['sha256']
        assert w.basis_source['feasibility_transferred'] is False and not w.basis_source['target_svd_tolerance_applied']
        assert (w.basis_source['source_L'],w.basis_source['target_L'],w.basis_source['cross_L'])==(source_L,target_L,source_L!=target_L)
        from smatrix_bootstrap.sdp import figures,sdpb
        loaded=figures.load_reports(str(tmp_path))[0]
        assert loaded['basis']['sha256']==record['basis']['sha256'] and figures.accepted_support(loaded)
        monkeypatch.setattr(sdpb,'execute',lambda *args,**kwargs:{'returncode':99})
        fresh,_,_=sdpb.run_once(target,tmp_path/'fresh',settings=sdpb.Settings(digits=30),basis_source_report=path)
        assert fresh['basis']['sha256']==record['basis']['sha256'] and not fresh['spec']['uv']
        assert fresh['basis']['source']['constraints_reassembled'] and 'reused_source_report' not in fresh
