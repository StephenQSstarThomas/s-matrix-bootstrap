"""Analysis layer (Phase B, Task B3): checks against the paper's own algebra and against synthetic amplitudes."""
from __future__ import annotations

import numpy as np
import pytest


def test_lbar_algebra_reproduces_table_3_from_tables_1_and_2():
    from smatrix_bootstrap.gtb.analysis import lbar
    # paper Table 1 (a02, a22 in units 1e-4 m_pi^-5) and Table 2 (radii in fm^2), blue and red
    blue = lbar(a02=14.6e-4, a22=3.41e-4, r2S_fm2=0.64, r2V_fm2=0.394)
    red = lbar(a02=14.9e-4, a22=2.53e-4, r2S_fm2=0.55, r2V_fm2=0.391)
    for got, want in ((blue, (0.96, 4.5, 4.8, 14.6)), (red, (-0.10, 4.8, 4.2, 14.5))):
        assert np.allclose([got["l1"], got["l2"], got["l4"], got["l6"]], want, atol=0.12)


def test_threshold_fit_recovers_known_parameters():
    from smatrix_bootstrap.gtb.analysis import threshold_fit
    for ell, a, b in ((0, 0.19, 0.28), (1, 0.031, 0.007), (2, 1.5e-3, 2e-4)):
        k2 = np.linspace(0.02, 0.25, 12)
        ref = k2 ** ell * (a + b * k2 + 0.01 * k2 ** 2)
        fa, fb = threshold_fit(k2, ref, ell)
        assert abs(fa - a) < 1e-12 and abs(fb - b) < 1e-10


def test_phase_convention_and_crossing_on_a_breit_wigner():
    from smatrix_bootstrap.gtb.analysis import crossing, phases
    E = np.linspace(0.3, 1.5, 200)
    m, g = 0.8, 0.15
    delta = np.arctan2(m * g, m ** 2 - E ** 2)              # elastic BW, rises through 90 deg at E = m
    S = np.exp(2j * delta)
    d, eta = phases(S)
    assert np.allclose(d, np.degrees(delta)) and np.allclose(eta, 1)
    assert abs(crossing(E, d, 90.0) - m) < 1e-4


def test_radius_units():
    from smatrix_bootstrap.gtb.analysis import radius_fm2
    # <r^2> = 6 F'(0) in m_pi^-2; 1 m_pi^-2 = (197.327/139.57)^2 fm^2
    assert abs(radius_fm2(1.0 / 6.0) - (197.327 / 139.57) ** 2) < 1e-12


MEV = 1 / 139.57                      # MeV -> m_pi units


@pytest.mark.parametrize("form", ["bw", "gs"])
def test_rho_fits_recover_known_mass_and_width(form):
    from smatrix_bootstrap.gtb.analysis import FF_FORMS, fit_rho
    m, g = 838 * MEV, 114 * MEV
    s = np.linspace((500 * MEV) ** 2, (1100 * MEV) ** 2, 30)
    data = np.abs(FF_FORMS[form](s, m, g))
    fm, fg = fit_rho(s, data, form, guess=(800 * MEV, 150 * MEV))
    assert abs(fm - m) < 1e-8 and abs(fg - g) < 1e-8


def test_two_bw_fit_recovers_known_parameters():
    from smatrix_bootstrap.gtb.analysis import fit_two_bw, two_bw
    p = (824 * MEV, 137 * MEV, 1635 * MEV, 248 * MEV, 1.29, 0.27)
    s = np.linspace((500 * MEV) ** 2, (1830 * MEV) ** 2, 40)
    got = fit_two_bw(s, np.abs(two_bw(s, *p)) ** 2, guess=(800 * MEV, 150 * MEV, 1600 * MEV, 300 * MEV, 1.2, 0.2))
    assert np.allclose(got, p, rtol=1e-7)


def test_g_rho_pi_pi_from_pdg_mass_and_width_is_6():
    from smatrix_bootstrap.gtb.analysis import g_rho_pi_pi
    assert abs(g_rho_pi_pi(775 * MEV, 149 * MEV) - 6.0) < 0.05
