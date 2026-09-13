"""Self-checks 5a.6 (S matrix / unitarity cone), 5a.7 (form factors and Gram),
5a.8 (FESR quadrature and its cutoff bias) and 5a.9 (FESR target recomputation).
"""
import numpy as np
import pytest
from mpmath import mp, mpf, quad
from smatrix_bootstrap.sdp import constraints as C
from smatrix_bootstrap.sdp import formfactor as FF
from smatrix_bootstrap.sdp import grid
M = 50
S = grid.s_nodes(M)
# --------------------------------------------------------------- 5a.6
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_unitarity_cone_equivalent_to_eta_le_one(seed):
    """|S| <= 1  <=>  |h|^2 <= 2 Im h  with h = kappa f, from (2.11)."""
    rng = np.random.default_rng(seed)
    for _ in range(400):
        eta, delta = rng.uniform(0.0, 1.6), rng.uniform(-np.pi, np.pi)
        s = float(rng.uniform(4.05, 500.0))
        Smat = eta * np.exp(2j * delta)
        f = (Smat - 1.0) / (1j * grid.kappa(s))     # invert (2.11)
        h = grid.kappa(s) * f
        assert abs(Smat - (1 + 1j * grid.kappa(s) * f)) < 1e-12
        assert (abs(h) ** 2 <= 2 * h.imag + 1e-12) == (eta <= 1 + 1e-12)
def test_lambda_rescaling_preserves_the_cone():
    """|h|^2 <= 2 Im h  <=>  |h/L|^2 <= 2 (Im h)/L^2 for any L > 0 (multiply by L^2)."""
    rng = np.random.default_rng(7)
    for _ in range(2000):
        h = complex(rng.normal(), rng.normal()) * rng.choice([1e-6, 1.0, 1e3])
        L = float(np.exp(rng.uniform(-18, 2)))
        a = abs(h) ** 2 <= 2 * h.imag
        b = abs(h / L) ** 2 <= 2 * h.imag / L ** 2
        assert a == b
# --------------------------------------------------------------- 5a.7
def test_centrifugal_scaling_does_not_clip_high_spin_threshold_rows():
    from smatrix_bootstrap.sdp.assembly import Operators
    ops = Operators(20, 6)
    ops.set_cone_scaling("centrifugal")
    assert np.min(ops.lambda_analytic) < 1e-16
    np.testing.assert_array_equal(ops.row_scale, ops.lambda_analytic)
@pytest.mark.parametrize("constant,verdict", [(0., "certified_pass"), (1., "certified_fail")])
def test_independent_arb_scattering_requires_every_slack_ball(constant, verdict):
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    from smatrix_bootstrap.sdp.grid import n_amplitude_vars
    c = np.zeros(n_amplitude_vars(4)); c[0] = constant
    report = ArbAudit(4, 1).audit(c, chi_caliber="chi-b")
    assert report["unitarity_verdict"] == verdict
    assert report["chiral_verdict"] == verdict
    assert sum(report["unitarity_counts"].values()) == 12
@pytest.mark.parametrize("late,eta,accepted", [(True, 1., False), (False, .95, True), (False, .5, False)])
def test_rho_claim_uses_registered_low_energy_window(late, eta, accepted):
    from smatrix_bootstrap.sdp.accuracy import rho_diagnostic
    ph = {"E_GeV": [.3, .6, .78, .82, .9, 1.2, 1.5],
          "delta_deg": [0, 10, 15, 20, 22, 25, 90] if late else [0, 20, 70, 90, 110, 140, 150],
          "eta": [eta] * 7}
    report = rho_diagnostic(ph)
    assert report["single_point_C6_pass"] is accepted
    if late:
        assert report["low_energy_crossing_GeV"] is None
def test_raw_pmp_dual_replay_checks_psd_and_stationarity(tmp_path):
    import json
    from smatrix_bootstrap.sdp.dual import replay_dual
    pmp = {"objective": ["0", "1"], "normalization": ["1", "0"],
           "PositiveMatrixWithPrefactorArray": [{"polynomials": [[[["2"], ["-1"]]]]}]}
    (tmp_path / "pmp.json").write_text(json.dumps(pmp))
    (tmp_path / "out").mkdir()
    (tmp_path / "out/x_0.txt").write_text("1 1\n1\n")
    good = replay_dual(tmp_path)
    assert good["emitted_pmp_bound_certified"]
    assert good["upper_bound_midpoint"] == 2
    (tmp_path / "out/x_0.txt").write_text("1 1\n0.9\n")
    assert not replay_dual(tmp_path)["emitted_pmp_bound_certified"]
def test_arb_amplitude_does_not_round_small_density_to_zero():
    from flint import arb
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    from smatrix_bootstrap.sdp.grid import n_amplitude_vars
    c = [arb(0) for _ in range(n_amplitude_vars(4))]
    c[0] = arb("1e-400")
    report = ArbAudit(4, 1).audit(c, chi_caliber=None)
    assert report["unitarity_counts"]["certified_fail"] == 8
    assert not report["unitarity_no_certified_violation"]
def test_phase_claim_rejects_large_rms_even_with_correct_endpoint():
    from smatrix_bootstrap.sdp.claims import c7
    series = lambda v: {'E_GeV':[.28,1.196],'delta_deg':[v,v]}
    points = {name:{'S0':series(95.),'S2':series(-25.)} for name in ('tip','ref','mid')}
    assert c7(points)['verdict'] == 'FAIL'
def test_kinematic_factor_squares_match_paper():
    """(2.33) squared reproduces (2.46)-(2.47)."""
    s = np.concatenate([S, np.linspace(4.01, 900.0, 57)])
    for ell in (0, 1):
        got = FF.kinematic_factor(ell, s) ** 2
        ref = FF.kinematic_square_reference(ell, s)
        assert np.abs(got / ref - 1).max() < 1e-13
def test_gram_realification_is_psd_equivalent():
    """3x3 Hermitian >= 0  <=>  [[Re,-Im],[Im,Re]] >= 0."""
    rng = np.random.default_rng(11)
    agree = 0
    for _ in range(3000):
        eta, d = rng.uniform(0.5, 1.3), rng.uniform(-np.pi, np.pi)
        Smat = eta * np.exp(2j * d)
        cF = complex(rng.normal(), rng.normal()) * 0.4
        rho = float(rng.uniform(0.0, 2.0))
        B = FF.gram_block(Smat, cF, rho)
        identity = (1-abs(Smat)**2)*(rho-abs(cF)**2)-abs(cF-Smat*cF.conjugate())**2
        assert abs(np.linalg.det(B)-identity)<1e-12
        h_psd = np.linalg.eigvalsh(B).min() >= -1e-12
        r_psd = np.linalg.eigvalsh(FF.realify(B)).min() >= -1e-12
        assert h_psd == r_psd
        agree += h_psd
    assert 0 < agree < 3000          # the test actually explores both sides
def test_gram_congruence_leaves_feasibility_alone():
    """V = diag(1,1,1/g) congruence: cF -> cF/g, rho -> rho/g^2."""
    rng = np.random.default_rng(13)
    for _ in range(2000):
        Smat = rng.uniform(0.3, 1.0) * np.exp(2j * rng.uniform(-np.pi, np.pi))
        cF = complex(rng.normal(), rng.normal()) * 1e-2
        rho = float(rng.uniform(0.0, 1e-3))
        g = float(rng.uniform(1e-3, 1.0))
        a = np.linalg.eigvalsh(FF.gram_block(Smat, cF, rho)).min() >= -1e-14
        b = np.linalg.eigvalsh(FF.gram_block(Smat, cF / g, rho / g ** 2)).min() >= -1e-14
        assert a == b
def test_real_part_operator_is_the_paper_kernel():
    """(3.66): ReF = 1 + K ImF, with K from (3.67)."""
    const, K = FF.real_part_operator(M)
    assert np.allclose(const, 1.0)
    phi = grid.phi_nodes(M)
    # F(s) = 1 + z(s) is analytic on the disk with F(0) = 1, so ImF = sin phi
    assert np.abs(const + K @ np.sin(phi) - (1.0 + np.cos(phi))).max() < 1e-10
# --------------------------------------------------------------- 5a.8
def test_fesr_quadrature_and_cutoff_bias():
    """(3.72) against the exact integral for smooth densities.
    Two effects are separated.  (i) The midpoint rule itself is accurate to a
    few 0.1 per mille.  (ii) The *range* it actually covers is not [4, s0]: the
    sum keeps every node with s_i <= s0, and the quadrature cell of the last
    kept node (i = 43) ends at phi = 43 pi / M, i.e. at s = 84.06, overshooting
    s0 = 73.47 by 14.4%.  That is a systematic bias of the paper's prescription,
    biggest for the highest moment; it is reported, not silently absorbed.
    """
    mp.dps = 30
    n0 = int(C.below_s0(M).sum())
    s_edge = 8.0 / (1.0 + np.cos(np.pi * n0 / M))
    assert s_edge == pytest.approx(84.0575, abs=1e-3)
    densities = {"(1+x/20)^-2": lambda x: 1.0 / (1.0 + x / 20.0) ** 2,
                 "exp(-x/30)": lambda x: mp.e ** (-x / 30) if isinstance(x, mpf) else np.exp(-x / 30),
                 "x^-1.5": lambda x: x ** -1.5}
    report = {}
    for name, rho in densities.items():
        for n in (-1, 0, 1):
            disc = float(C.moment_row(M, n) @ rho(S))
            e_edge = float(quad(lambda x: rho(x) * x ** n, [4, mpf(float(s_edge))]))
            e_s0 = float(quad(lambda x: rho(x) * x ** n, [4, mpf(float(grid.S0))]))
            report[(name, n)] = (disc / e_edge - 1.0, disc / e_s0 - 1.0)
            assert abs(disc / e_edge - 1.0) < 5e-3          # the rule is fine
    # ... but the cutoff bias is large and grows with the moment
    assert all(abs(v[1]) > 5e-3 for k, v in report.items() if k[1] == 1)
    assert max(abs(v[1]) for v in report.values()) > 0.05
def test_fesr_node_cutoff_geometry():
    below = C.below_s0(M)
    assert below.sum() == 43
    assert S[42] < grid.S0 < S[43]
# --------------------------------------------------------------- 5a.9
def test_fesr_targets_recomputed_from_svz():
    """Printed (2.56) vs an independent evaluation of (2.50).
    The S0 targets pin the quark-mass convention: they agree with
    m_q = sqrt((m_u^2+m_d^2)/2) to <1%, and disagree with the arithmetic mean by
    ~9%.  P1 carries no quark mass and agrees to 0.1%.
    """
    a = C.fesr_target_audit()
    for r in a["rows"]:
        assert abs(r["ratio_rms"] - 1) < 0.01
        if r["wave"] == "S0":
            assert abs(r["ratio_mean"] - 1) > 0.08
        else:
            assert abs(r["ratio_mean"] - 1) < 0.01
def test_ff_asymptotic_bound_counts():
    idx, bounds, kin = C.ff_asymptotic_bounds(M)
    assert idx.size == 7                    # 50 - 43 nodes above s0
    assert bounds[0] == pytest.approx(np.sqrt(2 * grid.M_Q ** 2 * C.EPS_FF))
    assert bounds[1] == pytest.approx(np.sqrt(0.5 * C.EPS_FF))
def test_ff_factor_reading_matters_only_for_P1():
    """(3.75) with the (2.33) factor frozen at s0 -- the paper's literal text --
    against the per-node reading.  k_0 is flat across the nodes above s0, k_1
    rises by a factor 15, so the two readings differ only for P1."""
    _, _, frozen = C.ff_asymptotic_bounds(M, frozen_at_s0=True)
    _, _, pernode = C.ff_asymptotic_bounds(M, frozen_at_s0=False)
    r0 = pernode[0] / frozen[0]
    r1 = pernode[1] / frozen[1]
    assert np.allclose(frozen[0], frozen[0][0]) and np.allclose(frozen[1], frozen[1][0])
    assert r0.max() < 1.02                  # S0: the factor is flat
    assert r1.max() > 10.0                  # P1: it is not
def test_chiral_reference_point():
    x, y = C.chiral_reference_point()
    assert x == pytest.approx(0.073321, abs=1e-6)
    assert y == pytest.approx(-x / 15.0)
# --------------------------------------------------------------- P1 functional
def test_lambda_row_reproduces_the_weinberg_value():
    """Below (2.14): lambda = (pi/4) T_{33,33}(4/3,4/3,4/3) = 3 pi/4 A(4/3,4/3,4/3).
    On a configuration where the amplitude is the constant T0, the row must give
    (3 pi/4) T0; and the Weinberg amplitude (2.12) must give the paper's
    lambda = m_pi^2 / (32 pi f_pi^2) = 0.023.
    """
    from smatrix_bootstrap.sdp.projector import Layout, PartialWaveOperator
    op = PartialWaveOperator(M)
    c = np.zeros(Layout(M).n)
    c[0] = 1.0
    assert op.lambda_row() @ c == pytest.approx(0.75 * np.pi)
    # Weinberg: A = (s-1)/(8 pi^2 f_pi^2), so A(4/3) = (1/3)/(8 pi^2 f_pi^2)
    a_weinberg = (4.0 / 3.0 - 1.0) / (8 * np.pi ** 2 * grid.F_PI ** 2)
    assert 0.75 * np.pi * a_weinberg == pytest.approx(
        1.0 / (32 * np.pi * grid.F_PI ** 2))
    assert 0.75 * np.pi * a_weinberg == pytest.approx(0.023, abs=5e-4)
def test_lambda_row_is_crossing_symmetric_by_construction():
    """All three Mandelstam slots sit at 4/3, so sigma2 carries twice the weight
    of sigma1 and rho1 twice that of rho2's folded diagonal."""
    from smatrix_bootstrap.sdp.projector import Layout, PartialWaveOperator
    op, lay = PartialWaveOperator(M), Layout(M)
    row = op.lambda_row()
    assert np.allclose(row[lay.s2], 2.0 * row[lay.s1])
# --------------------------------------------------------------- Arb audit
def test_arb_det3_matches_numpy():
    """arbaudit._det3 is the determinant of the Hermitian (3.68) block."""
    from flint import acb, arb, ctx
    from smatrix_bootstrap.sdp.arbaudit import _det3
    ctx.prec = 256
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(200):
        S = complex(rng.normal(), rng.normal()) * 0.6
        F = complex(rng.normal(), rng.normal()) * 0.5
        rho = float(rng.uniform(0.0, 2.0))
        ref = np.linalg.det(FF.gram_block(S, F, rho)).real
        got = float(_det3(acb(S.real, S.imag), acb(F.real, F.imag),
                          arb(rho)).str(25, radius=False))
        worst = max(worst, abs(got - ref) / max(abs(ref), 1e-12))
    assert worst < 1e-10
def test_arb_legendre_q_matches_the_float_implementation():
    """The Arb Q_ell used by the audit agrees with the float64 one."""
    from flint import arb, ctx
    from smatrix_bootstrap.sdp.arbaudit import _legendre_q
    from smatrix_bootstrap.sdp.legendreq import q_table
    ctx.prec = 384
    for z in (1.5, 3.0, 40.0, 1e4, -2.5, -700.0):
        col = q_table(np.array([z]), 12)[:, 0]
        for ell in (0, 1, 2, 5, 9, 12):
            got = float(_legendre_q(ell, arb(z)).str(25, radius=False))
            assert abs(got - col[ell]) <= 1e-11 * max(abs(col[ell]), 1e-300)
def test_phase_reconstructs_s_and_stops_at_zero():
    """A pi/2 patch changes S's sign; a zero has no threshold-connected phase."""
    from types import SimpleNamespace
    from smatrix_bootstrap.sdp.observables import phase_shift
    for S in (np.exp(1j * np.deg2rad([-170., -150., -120.])),
              np.exp(1j * np.deg2rad([0., 140., 220.])),
              np.array([1., 1e-15 * np.exp(0.7j), np.exp(1.1j)]),
              np.array([1., 0., -1.])):
        ops = SimpleNamespace(M=3, s=np.array([4.01, 4.024775834628433, 16.]), index=[(0, 0)],
                              h_re=S.imag[:, None], h_im=(1-S.real)[:, None])
        ph = phase_shift(ops, np.ones(1), "S0")
        if np.any(S == 0):
            assert np.isnan(ph["delta_deg"][1:]).all()
            assert ph["first_zero_node"] == 1
        else:
            recovered = ph["eta"] * np.exp(2j * np.deg2rad(ph["delta_deg"]))
            np.testing.assert_allclose(recovered, S, atol=5e-15)
            if S[0].imag < 0:
                assert abs(ph["delta_deg"][0] + 85.) < 1e-12
def test_compiled_gram_is_the_unitarily_equivalent_real_three_by_three():
    """The conjugate-paired form permits a real 3x3 cone, without doubling it."""
    from smatrix_bootstrap.sdp.problem import _gram_psd
    U = np.array([[1, 1j, 0], [1, -1j, 0], [0, 0, np.sqrt(2)]]) / np.sqrt(2)
    rng = np.random.default_rng(941)
    for _ in range(16):
        S, F = rng.normal(size=2) + 1j * rng.normal(size=2)
        rho = float(rng.normal())
        cone = _gram_psd(S.real, S.imag, F.real, F.imag, rho)
        assert cone.shape == (3, 3)
        expected = U.conj().T @ FF.gram_block(S, F, rho) @ U
        np.testing.assert_allclose(cone.args[0].value, expected, atol=2e-15)
@pytest.mark.parametrize('claim', ['C6', 'C8'])
def test_missing_physics_evidence_cannot_pass(claim):
    from smatrix_bootstrap.sdp import claims
    point = {'P1': {'crossing_90_GeV': .82, 'min_eta_below_1p2GeV': 1.},
             'S0': {'E_GeV': [.9, 1.1], 'delta_deg': [95., 95.]}}
    if claim == 'C6':
        result = claims.c6({'tip': point})
    else:
        runs = {(50, 8): {'tip': point}, (50, 10): {'tip': point},
                (50, 12): {'tip': point}, (45, 10): {'tip': point},
                (60, 10): {'tip': dict(point, P1=dict(point['P1'], crossing_90_GeV=.77))}}
        result = claims.c8(runs)
    assert result['verdict'] != 'pass'
def test_reduced_basis_identity_separates_figure_and_claim_grouping():
    import copy
    from smatrix_bootstrap.sdp import claims,figures
    from smatrix_bootstrap.sdp.spec import ModelSpec
    spec=ModelSpec(chiral=True,reduce_basis=True).__dict__
    a={'solver':'SDPB','accepted':True,'verification':{'primal_feasible':True},
       'convergence':{'solver_optimal':True},'spec':spec,'basis':{'sha256':'a'*64},
       'job':'tip','fix_f00':None,'result':{'direction':[1.,0.],'f00_3':.0826,'f11_3':-.001}}
    b=copy.deepcopy(a);b['basis']['sha256']='b'*64
    assert figures.model_key(a)!=figures.model_key(b)
    unknown=copy.deepcopy(a);unknown.pop('basis')
    assert not figures.accepted_support(unknown)
    with pytest.raises(ValueError):figures.model_key(unknown)
    assert claims.c2([a,b])['verdict']=='not run' and 'contract' in claims.c2([a,b])['evidence']
    by_eps={e:copy.deepcopy(a) for e in (.002,.004,.006)}
    for e,r in by_eps.items():r['spec']['eps_chi']=e
    by_eps[.004]['basis']['sha256']='b'*64
    assert claims.c3(by_eps)['verdict']=='not run'
    records=[]
    for uv in (False,True):
        tip=copy.deepcopy(a);tip['spec']['uv']=uv;tip['result']['f00_3']=.0811 if uv else .0826
        tip['basis']['sha256']=('b' if uv else 'a')*64;records.append(tip)
        for name,sign,y in [('section_hi',1.,-.001-(.00023 if uv else 0)),('section_lo',-1.,-.002+(.00005 if uv else 0))]:
            r=copy.deepcopy(tip);r.update(job=name,fix_f00=C.chiral_reference_point()[0])
            r['result'].update(direction=[0.,sign],f11_3=y);records.append(r)
    assert claims.c5(records)['verdict']=='not run'
@pytest.mark.parametrize('override,wanted',[(None,'1'),('3','3')])
def test_public_selfcheck_bounds_blas_threads(monkeypatch,override,wanted):
    import subprocess
    from smatrix_bootstrap.sdp.__main__ import main
    seen=[]
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):monkeypatch.setenv(name,'64')
    if override is None:monkeypatch.delenv('SDP_TEST_THREADS',raising=False)
    else:monkeypatch.setenv('SDP_TEST_THREADS',override)
    monkeypatch.setattr(subprocess,'call',lambda command,**kw:seen.append((command,kw.get('env',{}))) or 0)
    assert main(['selfcheck'])==0 and 'pytest' in seen[0][0]
    assert all(seen[0][1].get(name)==wanted for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'))
@pytest.mark.parametrize('direction,fixed,point,valid',[
    ([0,-1],.07,None,True),([2,1],None,None,True),([0,0],None,None,False),
    ([float('nan'),1],None,None,False),([1,0],float('inf'),None,False),
    ([0,-1],.07,'ref',False),(None,.07,None,False),(None,None,None,True)])
def test_saved_support_direction_contract(tmp_path,monkeypatch,direction,fixed,point,valid):
    import json
    from smatrix_bootstrap.sdp import sdpb
    from smatrix_bootstrap.sdp.spec import ModelSpec
    source={'solver':'SDPB','status':'numerically_accepted','accepted':True,
            'verification':{'primal_feasible':True,'f00_3':.2},'convergence':{'solver_optimal':True},
            'spec':ModelSpec().__dict__,'settings':sdpb.Settings().__dict__,'direction':[1.,0.],'fix_f00':None}
    path=tmp_path/'report.json';path.write_text(json.dumps(source));seen=[]
    monkeypatch.setattr(sdpb,'run_once',lambda *a,**k:seen.append((a,k)) or ({'accepted':False},None,None))
    if not valid:
        with pytest.raises(ValueError):sdpb.support_from_saved(path,tmp_path/'out',point,direction=direction,fix_f00=fixed)
    else:
        sdpb.support_from_saved(path,tmp_path/'out',point,direction=direction,fix_f00=fixed)
        assert list(seen[0][0][2])==(direction if direction is not None else [0.,1.])
        assert seen[0][0][3]==(fixed if direction is not None else C.chiral_reference_point()[0])
        assert seen[0][1]['prepared']==path and 'exact_input' not in seen[0][1]
