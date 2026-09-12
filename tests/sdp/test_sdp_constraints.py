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
    idx, bounds = C.ff_asymptotic_bounds(M)
    assert idx.size == 7                    # 50 - 43 nodes above s0
    assert bounds[0] == pytest.approx(np.sqrt(2 * grid.M_Q ** 2 * C.EPS_FF))
    assert bounds[1] == pytest.approx(np.sqrt(0.5 * C.EPS_FF))


def test_chiral_reference_point():
    x, y = C.chiral_reference_point()
    assert x == pytest.approx(0.073321, abs=1e-6)
    assert y == pytest.approx(-x / 15.0)
