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


def t_of_mu(s, mu):
    return -(s - 4.0) * (1.0 - mu) / 2.0


def u_of_mu(s, mu):
    return -(s - 4.0) * (1.0 + mu) / 2.0


# ------------------------------------------------------- 5a.3: angular kernels
@pytest.mark.parametrize("ell,s,x", [(0, 3.0, 4.5), (1, 3.0, 30.0), (2, 0.5, 9.0),
                                     (3, 12.0, 4.2), (5, 7.756, 100.0),
                                     (10, 200.0, 7.0), (19, 40.0, 60.0), (4, 1.5, 1000.0)])
def test_Qc_closed_form(ell, s, x):
    """Qc_ell(x) = (2/a) Q_ell(1 + x/a) reproduces the direct mu integral."""
    mp.dps = 50
    ref = float(quad(lambda m: legendre(ell, m) / (mpf(x) - t_of_mu(mpf(s), m)), [-1, 1]))
    a = 0.5 * (s - 4.0)
    got = (2.0 / a) * OP.qc.get(np.array([1.0 + x / a]), 19)[ell, 0]
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
