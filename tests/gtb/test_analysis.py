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


def test_g_f2_pi_pi_reproduces_the_paper_from_table_7():
    from smatrix_bootstrap.gtb.analysis import g_f2_pi_pi
    assert abs(g_f2_pi_pi(1366 * MEV, 164 * MEV) - 8.5) < 0.05
    assert abs(g_f2_pi_pi(1388 * MEV, 188 * MEV) - 8.9) < 0.05


from pathlib import Path

LEAF10W = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10_watson")
KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")


@pytest.mark.skipif(not (LEAF10W / "solution.npz").exists(), reason="M=10 Watson leaf missing")
def test_analyze_runs_end_to_end_on_the_m10_watson_leaf():
    from smatrix_bootstrap.gtb.leafanalysis import analyze
    r = analyze(LEAF10W, KS10, pole=True)
    for w in ("S0", "D0", "S2", "D2", "P1", "F1"):
        assert len(r["phases"][w]["delta_deg"]) == r["n0"] and max(r["phases"][w]["eta"]) <= 1 + 1e-6
    assert set(r["threshold"]) == {"a00", "b00", "a20", "b20", "a11", "b11", "a02", "a22"}
    assert all(np.isfinite(v) for v in r["threshold"].values())
    assert set(r["radii_fm2"]) == {"S", "V", "T"} and set(r["lbar"]) == {"l1", "l2", "l4", "l6"}
    if r["pole"]["interp"] is not None:
        assert r["pole"]["interp"]["sqrt_s_MeV"][1] < 0


def test_gtbii_reproduces_table_9_from_table_3():
    from smatrix_bootstrap.gtb.analysis import gtbii
    blue = gtbii({"l1": 0.96, "l2": 4.5, "l4": 4.8, "l6": 14.6})
    red = gtbii({"l1": -0.10, "l2": 4.8, "l4": 4.2, "l6": 14.5})
    want_b = {"a00": 0.21, "a20": -0.043, "b00": 0.26, "b20": -0.080, "a11": 36.3e-3, "b11": 2.62e-3, "a02": 14.8e-4, "a22": 3.4e-4}
    want_r = {"a00": 0.21, "a20": -0.043, "b00": 0.26, "b20": -0.079, "a11": 36.7e-3, "b11": 3.43e-3, "a02": 15.0e-4, "a22": 2.5e-4}
    for got, want in ((blue, want_b), (red, want_r)):
        for k, w in want.items():
            assert abs(got[k] - w) <= 0.06 * abs(w), (k, got[k], w)


@pytest.mark.skipif(not (LEAF10W / "solution.npz").exists(), reason="M=10 Watson leaf missing")
def test_cli_analyze_writes_analysis_json():
    import json, os, subprocess, sys
    out = subprocess.run([sys.executable, "-m", "smatrix_bootstrap.run", "gtb", "analyze", "--leaf", str(LEAF10W),
                          "--kernels", str(KS10)], capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"})
    assert out.returncode == 0, out.stderr[-2000:]
    an = json.loads((LEAF10W / "analysis.json").read_text())
    assert {"phases", "threshold", "fits", "forward", "pole"} <= set(an)


@pytest.mark.skipif(not (LEAF10W / "solution.npz").exists(), reason="M=10 Watson leaf missing")
def test_lec_rows_reproduce_the_leaf_analysis():
    import json
    from smatrix_bootstrap.gtb.lec import lec_rows
    from smatrix_bootstrap.gtb.leafanalysis import _load, analyze
    rep, ks, prog, v = _load(LEAF10W, KS10)
    rows = lec_rows(prog)
    an = analyze(LEAF10W, KS10, pole=False)
    for key in ("a00", "a20", "a11", "b00", "a02", "a22"):
        assert abs(rows[key] @ v - an["threshold"][key]) <= 1e-10 * max(1, abs(an["threshold"][key])), key
    for key, ch in (("r2S", "S"), ("r2V", "V"), ("r2T", "T")):
        assert abs(rows[key] @ v - an["radii_fm2"][ch]) <= 1e-12 * max(1, abs(an["radii_fm2"][ch])), key


@pytest.mark.skipif(not (LEAF10W / "solution.npz").exists(), reason="M=10 Watson leaf missing")
def test_watson_residual_is_zero_at_own_targets_and_sees_a_shifted_target():
    from smatrix_bootstrap.gtb.leafanalysis import _load
    from smatrix_bootstrap.gtb.watson import residual, targets_from_solution
    _, _, prog, v = _load(LEAF10W, KS10)
    t = targets_from_solution(prog, v)
    r = residual(prog, v, t)
    assert all(r[w] == 0.0 for w in ("S2", "D2", "F1"))
    t["F1"][1][2] += 0.125
    assert abs(residual(prog, v, t)["F1"] - 0.125) < 1e-12


KS10_L24 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10_L24")


@pytest.mark.skipif(not (KS10_L24 / "manifest.json").exists() or not (LEAF10W / "solution.npz").exists(),
                    reason="M=10 L=24 export or leaf missing")
def test_omitted_waves_come_from_a_larger_export_of_the_same_amplitude():
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.leafanalysis import _load, omitted_waves
    _, _, prog, v = _load(LEAF10W, KS10)
    amp = v[prog.groups["amp"]]
    big = kernels.load(KS10_L24)
    res = omitted_waves(amp, big, 20, include_kept=True)
    S = 1 + 1j * (prog.Rre @ amp + 1j * (prog.Rim @ amp))
    assert abs(res["kept"]["P1"][3] - abs(S[prog.row_index("P1", 3)])) < 1e-12      # same rows, same S
    assert set(res["eta_max"]) == {"I0_l20", "I0_l22", "I2_l20", "I2_l22", "I1_l21", "I1_l23"}
    assert all(np.isfinite(x) for x in res["eta_max"].values())
