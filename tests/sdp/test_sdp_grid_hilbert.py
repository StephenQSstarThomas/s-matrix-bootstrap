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
