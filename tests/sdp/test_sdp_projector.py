"""Self-checks 5a.3 (crossed-channel Legendre-Q closed forms), 5a.4 (subthreshold
rows) and 5a.5 (isospin / parity selection rule).

The reference side shares no code with the object under test: the angular
integrals are redone with mpmath quadrature, and the whole partial wave is
redone by evaluating A(nu1,nu2,nu3) pointwise in mu and integrating with
Gauss-Legendre.
"""
import numpy as np
import pytest
from mpmath import mp, mpf, legendre, quad
from smatrix_bootstrap.sdp import grid
from smatrix_bootstrap.sdp.projector import Layout, PartialWaveOperator, ells_for
M = 50
OP = PartialWaveOperator(M)
S = grid.s_nodes(M)
OMEGA = grid.dsdphi(M) / M
LAY = Layout(M)
@pytest.mark.parametrize('tol',[1e-12,1e-3])
def test_basis_is_identical_for_equivalent_disks_grams_and_scalings(tol):
    from smatrix_bootstrap.sdp.assembly import Operators
    ops=Operators(4,2,'sine-cardinal',40)
    mask=np.ones(24,dtype=bool)
    ops.set_cone_scaling('none');reference=ops.build_basis(tol,mask=mask).copy()
    if tol==1e-3:assert reference.shape[1]<ops.lay.n
    for wave in ((0,0),(1,1)):
        start=ops.index.index(wave)*4;mask[start:start+4]=False
    for scaling in ('none','rownorm','centrifugal'):
        ops.set_cone_scaling(scaling)
        assert np.array_equal(reference,ops.build_basis(tol,mask=mask,include_gram=True))
        assert np.array_equal(reference,ops.build_basis(tol,mask=np.ones(24,dtype=bool)))
def test_resolution_L12_contains_twelve_waves_in_each_isospin():
    assert ells_for(0, 12) == tuple(range(0, 24, 2))
    assert ells_for(1, 12) == tuple(range(1, 24, 2))
    assert ells_for(2, 12) == tuple(range(0, 24, 2))
    row = PartialWaveOperator(4).rows(1, 23, 40., 1)
    assert np.all(np.isfinite(row))
    assert np.max(abs(row)) > 0
@pytest.mark.parametrize("I,ell,s,node", [(0,0,None,2), (1,3,None,4), (2,0,"1.5",None)])
def test_precise_rows_against_independent_arb_amplitude(I, ell, s, node):
    from flint import arb, arb_mat
    from smatrix_bootstrap.sdp.precision import PrecisionRows
    from smatrix_bootstrap.sdp.arbaudit import ArbAudit
    p = PrecisionRows(6, 100)
    audit = ArbAudit(6, 3, bits=448)
    c = np.random.default_rng(987).normal(size=Layout(6).n) * 1e12
    v = Layout(6).unpack(c)
    dens = {"T0": arb(float(v['T0'])),
            **{k: [arb(float(x)) for x in v[k]] for k in ('sigma1','sigma2')},
            **{k: [[arb(float(x)) for x in row] for row in v[k]] for k in ('rho1','rho2')}}
    rows = p.rows(I, ell, point=s, node=node)
    value = arb_mat(rows.tolist()) * arb_mat([[arb(float(x))] for x in c])
    reference = audit.partial_wave(dens, I, ell, node or 0, point=s)
    assert abs(value[0,0]-reference.real) < arb('1e-85')
    assert abs(value[1,0]-reference.imag) < arb('1e-85')
def test_precise_pmp_exports_genuine_precision_and_degree_zero(tmp_path):
    import json
    from flint import arb
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.spec import ModelSpec
    w = Pmp(ModelSpec(M=4,L=2,chiral=True,uv=True,reduce_basis=True,operator_dps=50), digits=40)
    w.write(str(tmp_path/'pmp.json'))
    doc = json.loads((tmp_path/'pmp.json').read_text())
    blocks = doc['PositiveMatrixWithPrefactorArray']
    gram = blocks[int(w.keep.sum())+1]['polynomials']
    # The original real Gram has sqrt(2) F_re. A promoted double fails by ~1e-16.
    assert abs(arb(gram[0][2][0][0])-arb(2).sqrt()) < arb('1e-38')
    assert all(len(poly)==1 for block in blocks for row in block['polynomials']
               for cell in row for poly in cell)
def test_sine_completion_has_native_cardinality_and_unsubtracted_dispersion():
    from flint import arb, acb
    from smatrix_bootstrap.sdp.sine import SineFamily
    p = SineFamily(6,2,dps=40)
    assert all(x == 0 for x in p.cardinal(arb(-1)))
    for k in range(6):
        phi = arb.pi()*(2*k+1)/12
        row = p.cardinal(acb(phi.cos(),phi.sin()))
        for j,g in enumerate(row):
            assert abs(g.real-p.on[k,j]) < arb('1e-38')
            assert abs(g.imag-int(k==j)) < arb('1e-38')
    mp.dps=50; nu=mpf(3)
    ref=quad(lambda t: (4*mp.cos(t)**2-1)*(1-mp.cos(t))/(1-nu*mp.cos(t/2)**2/4)/mp.pi,[0,mp.pi])
    c=[(arb.pi()*3*(2*j+1)/12).sin() for j in range(6)]
    assert abs(sum((g*v for g,v in zip(p.cardinal_at(3),c)),arb(0))-arb(str(ref))) < arb('1e-38')
def t_of_mu(s, mu):
    return -(s - 4.0) * (1.0 - mu) / 2.0
def u_of_mu(s, mu):
    return -(s - 4.0) * (1.0 + mu) / 2.0
# ------------------------------------------------------- 5a.3: angular kernels
@pytest.mark.parametrize("ell,s,x", [(0, 3.0, 4.5), (1, 3.0, 30.0), (2, 0.5, 9.0),
                                     (3, 12.0, 4.2), (5, 7.756, 100.0),
                                     (10, 200.0, 7.0), (19, 40.0, 60.0), (4, 1.5, 1000.0),
                                     (23, 40.0, 60.0)])
def test_Qc_closed_form(ell, s, x):
    """Qc_ell(x) = (2/a) Q_ell(1 + x/a) reproduces the direct mu integral."""
    mp.dps = 50
    ref = float(quad(lambda m: legendre(ell, m) / (mpf(x) - t_of_mu(mpf(s), m)), [-1, 1]))
    a = 0.5 * (s - 4.0)
    got = (2.0 / a) * OP.qc.get(np.array([1.0 + x / a]), max(19, ell))[ell, 0]
    assert abs(got - ref) <= 1e-11 * abs(ref)
@pytest.mark.parametrize("ell,s,x,y", [(0, 3.0, 4.5, 9.0), (1, 12.0, 30.0, 4.2),
                                       (2, 0.5, 9.0, 60.0), (3, 7.756, 4.2, 4.2),
                                       (6, 200.0, 7.0, 1000.0), (19, 40.0, 60.0, 5.0),
                                       (7, 1.0, 12.0, 12.0), (12, 500.0, 4.1, 900.0)])
def test_Dc_closed_form(ell, s, x, y):
    """Dc_ell(x,y) = [Qc(x) + (-1)^ell Qc(y)] / (s-4+x+y) reproduces the direct integral."""
    mp.dps = 50
    ref = float(quad(lambda m: legendre(ell, m)
                     / ((mpf(x) - t_of_mu(mpf(s), m)) * (mpf(y) - u_of_mu(mpf(s), m))),
                     [-1, 1]))
    a = 0.5 * (s - 4.0)
    qx = (2.0 / a) * OP.qc.get(np.array([1.0 + x / a]), 19)[ell, 0]
    qy = (2.0 / a) * OP.qc.get(np.array([1.0 + y / a]), 19)[ell, 0]
    got = (qx + (-1.0) ** ell * qy) / (s - 4.0 + x + y)
    assert abs(got - ref) <= 1e-11 * abs(ref)
# ----------------------------------------- independent amplitude / projection
def _cw(nu, node=None):
    """Cauchy weights: C[f](nu) = cw . f, as a complex (..., M) array."""
    if node is not None:
        w = OP.A_on[node].astype(complex)
        w[node] += 1j
        return w
    nu = np.atleast_1d(np.asarray(nu, dtype=float))
    return (OMEGA[None, :] / (S[None, :] - nu[:, None])).astype(complex)
def _amp(c, w1, w2, w3):
    """A(nu1,nu2,nu3) = T0 + C[s1](n1) + C[s2](n2) + C[s2](n3)
                        + D[r1](n1,n2) + D[r1](n1,n3) + D[r2](n2,n3),
    with D[g](a,b) = cw(a)^T g cw(b).  Broadcasts over a leading mu axis."""
    v = LAY.unpack(c)
    s1, s2, r1, r2 = v["sigma1"], v["sigma2"], v["rho1"], v["rho2"]
    out = complex(v["T0"]) + w1 @ s1 + w2 @ s2 + w3 @ s2
    out = out + np.einsum("...i,ij,...j->...", w1, r1, w2)
    out = out + np.einsum("...i,ij,...j->...", w1, r1, w3)
    out = out + np.einsum("...i,ij,...j->...", w2, r2, w3)
    return out
def _f_reference(c, isospin, ell, s, node=None, nmu=2000):
    """f^I_ell(s) = (1/4) int P_ell(mu) T^I dmu by Gauss-Legendre, no Legendre-Q.

    Also returns the L1 norm of the integrand.  For high ell the projection is
    many orders of magnitude below that norm, i.e. the Gauss-Legendre sum
    cancels heavily; the achievable float64 accuracy of the *reference* is
    therefore set by the L1 norm, not by the answer, and that is the scale the
    agreement must be judged against.
    """
    mu, wq = np.polynomial.legendre.leggauss(nmu)
    t, u = t_of_mu(s, mu), u_of_mu(s, mu)
    ws = _cw(s, node)
    ws = np.broadcast_to(ws, (nmu, M)) if ws.ndim == 1 else np.broadcast_to(ws, (nmu, M))
    wt, wu = _cw(t), _cw(u)
    A_stu = _amp(c, ws, wt, wu)
    A_tsu = _amp(c, wt, ws, wu)
    A_uts = _amp(c, wu, wt, ws)
    T = {0: 3 * A_stu + A_tsu + A_uts, 1: A_tsu - A_uts, 2: A_tsu + A_uts}[isospin]
    # Bound on the reference's own float64 roundoff: the isospin combination
    # (2.5) and the ell projection both cancel, so the achievable accuracy is
    # set by the largest intermediate, not by the answer.
    Tabs = {0: 3 * np.abs(A_stu) + np.abs(A_tsu) + np.abs(A_uts)}.get(
        isospin, np.abs(A_tsu) + np.abs(A_uts))
    P = np.polynomial.legendre.Legendre.basis(ell)(mu)
    return 0.25 * np.sum(wq * P * T), 0.25 * np.sum(wq * np.abs(P) * Tabs)
@pytest.fixture(scope="module")
def sample_c():
    rng = np.random.default_rng(20260912)
    c = rng.standard_normal(LAY.n) * 1e-3
    r2 = c[LAY.r2]                       # already stored packed-symmetric
    c[LAY.r2] = r2
    return c
SUBTHRESHOLD = [0.5, 1.0, 1.5, 2.0, 3.0]
@pytest.mark.parametrize("isospin,ell,s", [(I, e, s) for s in SUBTHRESHOLD
                                           for I in (0, 1, 2) for e in ells_for(I, 3)])
def test_subthreshold_rows(sample_c, isospin, ell, s):
    """5a.4: below threshold the partial waves are real and match the direct integral."""
    got = OP.rows(isospin, ell, s) @ sample_c
    ref, l1 = _f_reference(sample_c, isospin, ell, s)
    assert abs(ref.imag) < 1e-13 * l1
    assert abs(got[1]) < 1e-18
    assert abs(got[0] - ref.real) <= 1e-11 * l1
@pytest.mark.parametrize("isospin,ell,node", [(I, e, n) for n in (0, 15, 30, 44)
                                              for I in (0, 1, 2) for e in ells_for(I, 3)])
def test_physical_region_rows(sample_c, isospin, ell, node):
    """5a.3/5a.6: on the cut, both real and imaginary parts match the direct integral."""
    s = S[node]
    got = OP.rows(isospin, ell, s, node) @ sample_c
    ref, l1 = _f_reference(sample_c, isospin, ell, s, node)
    assert abs(got[0] - ref.real) <= 1e-11 * l1
    assert abs(got[1] - ref.imag) <= 1e-11 * l1
def test_selection_rule_is_an_identity():
    """5a.5: the forbidden isospin/parity projections vanish to machine precision."""
    worst = 0.0
    for ell in range(20):
        for s, node in [(3.0, None), (0.5, None), (S[0], 0), (S[25], 25), (S[49], 49)]:
            worst = max(worst, OP.parity_residual(ell, s, node))
    assert worst < 1e-13
def test_crossing_symmetry_t_u(sample_c):
    """A(s,t,u) = A(s,u,t) at the implementation level (2.4)."""
    for s, t in [(3.0, -1.2), (0.5, -0.3), (12.0, -5.0)]:
        u = 4.0 - s - t
        a1 = _amp(sample_c, _cw(s)[0], _cw(t)[0], _cw(u)[0])
        a2 = _amp(sample_c, _cw(s)[0], _cw(u)[0], _cw(t)[0])
        assert abs(a1 - a2) < 1e-14 * max(abs(a1), 1e-12)
# Independent grid and Hilbert-transform controls, retained intact.
"""Self-checks 5a.1 (z map / grid) and 5a.2 (discrete Cauchy integral).

Independent mathematics only: nothing here calls the model being tested twice.
"""
import numpy as np
import pytest
from smatrix_bootstrap.sdp import grid, hilbert
# ------------------------------------------------------------------ 5a.1
def test_z_map_sends_nodes_to_unit_circle():
    """(3.58)-(3.60): z(s_i + i0) = exp(i phi_i)."""
    for M in (20, 50, 60):
        phi = grid.phi_nodes(M)
        z = grid.z_map(grid.s_nodes(M) + 1e-300j)
        assert np.abs(z - np.exp(1j * phi)).max() < 1e-12
def test_dsdphi_matches_numerical_derivative():
    """(3.71) against a central difference of nu(phi) = 8/(1+cos phi)."""
    for M in (20, 50):
        phi, h = grid.phi_nodes(M), 1e-6
        num = (8.0 / (1 + np.cos(phi + h)) - 8.0 / (1 + np.cos(phi - h))) / (2 * h)
        assert np.abs(num / grid.dsdphi(M) - 1).max() < 1e-8
# ------------------------------------------------------------------ 5a.2
def test_paper_kernel_equals_odd_extension_derivation():
    """(3.67) is the odd-offset cotangent conjugation on the 2M staggered grid."""
    for M in (8, 17, 50, 60):
        assert np.abs(hilbert.hilbert_kernel(M)
                      - hilbert.hilbert_kernel_from_odd_extension(M)).max() < 1e-13
def test_kernel_maps_sin_to_cos():
    """The conjugate function sends sin(n phi) -> cos(n phi) for n < M."""
    M = 50
    K, phi = hilbert.hilbert_kernel(M), grid.phi_nodes(M)
    for n in range(1, M):
        assert np.abs(K @ np.sin(n * phi) - np.cos(n * phi)).max() < 1e-10
def _reference(a, nu):
    """g(nu) = G(z(nu)) with G(z) = (1+z)/(a-z): analytic on the disk, real
    coefficients, G(-1) = 0, so g is the unsubtracted Cauchy transform of
    sigma(phi) = Im G(exp(i phi))."""
    z = grid.z_map(nu)
    return (1.0 + z) / (a - z)
def _sigma(a, M):
    phi = grid.phi_nodes(M)
    return (a + 1) * np.sin(phi) / (a * a - 2 * a * np.cos(phi) + 1)
@pytest.mark.parametrize("a", [1.2, 2.0, 5.0])
def test_cauchy_on_cut_real_part(a):
    """Re g(s_k + i0) = G(0) + (K sigma)_k, and convergence with M."""
    errs = {}
    for M in (25, 50, 100):
        got = hilbert.cauchy_oncut_real(M) @ _sigma(a, M)
        exact = _reference(a, grid.s_nodes(M) + 1e-300j).real
        errs[M] = np.abs(got - exact).max()
    assert errs[50] < 3e-3
    if errs[25] > 1e-12:                          # otherwise already at roundoff
        assert errs[100] < errs[50] < errs[25]
@pytest.mark.parametrize("a", [1.2, 2.0, 5.0])
def test_cauchy_off_cut(a):
    """Midpoint-in-phi quadrature of the same transform below the cut."""
    for nu in (-40.0, -5.0, 0.0, 1.5, 3.0, 3.9):
        got = hilbert.cauchy_offcut_row(50, nu) @ _sigma(a, 50)
        exact = _reference(a, complex(nu)).real
        assert abs(got - exact) < 5e-3 * max(abs(exact), 1.0)
def test_disk_series_reproduces_conjugate():
    """Independent route: interpolate sigma in sin(n phi), read off cos(n phi)."""
    M, a = 50, 2.0
    sig = _sigma(a, M)
    coef = hilbert.disk_series_coeffs(M) @ sig
    phi = grid.phi_nodes(M)
    re_series = np.cos(np.outer(phi, np.arange(1, M + 1))) @ coef
    assert np.abs(re_series - hilbert.hilbert_kernel(M) @ sig).max() < 1e-9
def test_infinity_constant_includes_nyquist_and_equals_midpoint():
    """A missing half weight at n=M changes G(-1) and the on-cut operator."""
    for M in (3, 20, 30, 50, 100):
        p = grid.phi_nodes(M)
        n = np.arange(1, M + 1)
        sine = np.sin(np.outer(p, n))
        a0 = hilbert.infinity_constant_row(M)
        # G(z)=(-1)^(n+1)+z^n vanishes at infinity z=-1.
        np.testing.assert_allclose(a0 @ sine, (-1.) ** (n + 1), atol=2e-13)
        np.testing.assert_allclose(a0, hilbert.cauchy_offcut_row(M, 0),
                                   rtol=2e-12, atol=2e-14)
        expected = (-1.) ** (n + 1) + np.cos(np.outer(p, n))
        np.testing.assert_allclose(hilbert.cauchy_oncut_real(M, "infinity") @ sine,
                                   expected, atol=3e-13)
def test_cot_kernel_beats_punctured_midpoint_pv():
    """Both discretise the same principal value, but only the (3.67) kernel is
    spectrally accurate: the punctured midpoint rule decays merely like 1/M.

    This is *why* the task prescribes the cot kernel on the cut; the numbers are
    reported in the reproduction report.
    """
    a = 2.0
    exact = lambda M: _reference(a, grid.s_nodes(M) + 1e-300j).real
    e_cot, e_mid = {}, {}
    for M in (50, 100, 200):
        sig = _sigma(a, M)
        e_cot[M] = np.abs(hilbert.cauchy_oncut_real(M) @ sig - exact(M)).max()
        e_mid[M] = np.abs(hilbert.cauchy_oncut_pv_midpoint(M) @ sig - exact(M)).max()
    assert e_cot[50] < 1e-12                       # spectral
    assert e_mid[50] > 1e-3                        # first order
    assert 1.6 < e_mid[100] / e_mid[200] < 2.4     # ~1/M decay
def test_physical_offnode_integral_encloses_imaginary_tail():
    from flint import arb
    from smatrix_bootstrap.sdp.sine import SineFamily
    f = SineFamily(4, 2, 40)
    c = [arb(0) for _ in range(f.lay.n)]; c[f.lay.s1.start+1] = 1
    value, _ = f.wave(c, 0, 0, point='8', order=2)
    # The crossed legs are real; Im f00=3 Im g_1(i)/2 exactly.
    expected = arb(3)/4*((3*arb.pi()/8).sin()+(arb.pi()/8).sin())
    assert value.imag.contains(expected)
@pytest.mark.parametrize('distances,errors,expected',[
    ((3.,1.),(.1,.1),'lower'),((1.,3.),(.1,.1),'upper'),((1.,1.1),(.2,.1),None),((1.,1.),(0.,0.),None)])
def test_ir_endpoint_choice_respects_support_gap(distances,errors,expected):
    from smatrix_bootstrap.sdp.ir_selection import distance_order
    assert distance_order(distances,errors)['selected']==expected
@pytest.mark.parametrize("kind,truncate", [(kind,cut) for kind in ("IR","UV") for cut in (False,True)])
def test_paper_comparison_requires_its_published_endpoint(kind,truncate):
    from smatrix_bootstrap.sdp.delivery import c4_source,c7_source
    from smatrix_bootstrap.sdp.figures import load_csv
    from smatrix_bootstrap.sdp.claims import PAIRING
    groups={"ref":"ir_magenta"} if kind=="IR" else PAIRING
    waves=("S0","S2","P1") if kind=="IR" else ("S0","S2")
    curves={}
    for point,group in groups.items():
        curves[point]={}
        for wave in waves:
            prefix="figure7" if kind=="IR" else "figure10"
            ref=load_csv(f"{prefix}_{wave.lower()}_phases.csv");mask=ref["group"]==group
            E,d=ref["energy_gev"][mask],ref["phase_deg"][mask];order=np.argsort(E)
            if truncate:order=order[:-1]
            curves[point][wave]={"E_GeV":E[order],"delta_deg":d[order]}
    verdict=c4_source(curves["ref"]) if kind=="IR" else c7_source(curves)
    assert verdict["verdict"]==("not run" if truncate else "pass")
@pytest.mark.parametrize('ours,ref,passed',[([1e-18],[0.],True),([1e-9],[0.],False),
    ([2e-40],[1e-40],False),([1+1e-11],[1.],True),([1+1e-9],[1.],False),
    ([np.nan],[0.],None),([0.],[np.nan],None),([0.,0.],[0.],None)])
def test_mma_exact_zero_roundoff_does_not_relax_nonzero_relative_check(ours,ref,passed):
    from smatrix_bootstrap.sdp.mma import _compare_values
    if passed is None:
        with pytest.raises(ValueError):_compare_values(ours,ref)
    else:
        result=_compare_values(ours,ref)
        assert result['passed'] is passed
        assert result['nonzero_relative_tolerance']==1e-10
