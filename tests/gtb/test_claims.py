"""Pre-registered claim rules (docs/gtb2403/PREREGISTRATION.md)."""
from __future__ import annotations


def test_g3_rule_on_brackets():
    from smatrix_bootstrap.gtb.claims import c1_verdict
    p = -0.028770
    assert c1_verdict(p - 2e-5, p + 3e-5, p) == "PASS"                 # bracket inside the band
    assert c1_verdict(p + 2e-4, p + 3e-4, p) == "FAIL"                 # wholly outside
    assert c1_verdict(p - 3e-4, p - 1.5e-4, p) == "FAIL"
    assert c1_verdict(p - 2e-4, p + 2e-4, p) == "INCONCLUSIVE"        # straddles, wide
    assert c1_verdict(p + 0.9e-4, p + 1.1e-4, p) == "INCONCLUSIVE"     # straddles the edge, narrow
    try:
        c1_verdict(p + 1e-4, p, p)
        raise AssertionError("an inverted bracket must be rejected")
    except ValueError:
        pass


import copy
from pathlib import Path

import pytest

PAPER = Path("/playpen1/shiqiu/gtb2403_runs/paper_data")
needs_paper = pytest.mark.skipif(not (PAPER / "index.json").exists(), reason="paper extraction missing")


@needs_paper
def test_phase_rules_accept_the_papers_own_curves():
    from smatrix_bootstrap.gtb.claims import paper_as_analysis, phase_claims
    an = {ni: paper_as_analysis(ni) for ni in (44, 46, 48)}
    res = phase_claims(an)
    for cid in ("C4", "C5", "C6", "C7", "C12"):
        assert res[cid]["verdict"] == "PASS", (cid, res[cid])


@needs_paper
def test_phase_rules_reject_a_shifted_p1_wave():
    import numpy as np
    from smatrix_bootstrap.gtb.claims import paper_as_analysis, phase_claims
    an = {ni: paper_as_analysis(ni) for ni in (44, 46, 48)}
    bad = copy.deepcopy(an)
    for ni in bad:
        bad[ni]["phases"]["P1"]["delta_deg"] = (np.array(bad[ni]["phases"]["P1"]["delta_deg"]) + 15).tolist()
    assert phase_claims(bad)["C4"]["verdict"] == "FAIL"


def _paper_numbers_as_analysis(colour):
    from smatrix_bootstrap.gtb.claims import TABLES
    t = {k: v[colour] for k, v in TABLES.items()}
    return {"threshold": {k: t[k] for k in ("a00", "a20", "b00", "b20", "a11", "b11", "a02", "a22")},
            "radii_fm2": {"S": t["r2S"], "V": t["r2V"], "T": t["r2T"]},
            "lbar": {k: t[k] for k in ("l1", "l2", "l4", "l6")},
            "pole": {"interp": {"sqrt_s_MeV": [t["pole_re"], -t["pole_im"]]}},
            "fits": {"rho_gs": {"m_MeV": t["m_rho"], "Gamma_MeV": t["G_rho"], "g_rho_pi_pi": t["g_rho"]},
                     "two_bw": {"m1_MeV": t["m1"], "G1_MeV": t["G1"], "m2_MeV": t["m2"], "G2_MeV": t["G2"],
                                "A1": t["A1"], "A2": t["A2"]},
                     "f2_bw": {"m_MeV": t["m_f2"], "Gamma_MeV": t["G_f2"], "g_f2_pi_pi": t["g_f2"]},
                     "four_pi_over_g_rho2": t["4pi_grho2"]}}


def test_table_rules_accept_the_papers_numbers_and_reject_a_wrong_one():
    from smatrix_bootstrap.gtb.claims import table_claims
    blue, red = _paper_numbers_as_analysis(0), _paper_numbers_as_analysis(1)
    res = table_claims(blue, red)
    assert all(v["verdict"] == "PASS" for k, v in res.items() if k.startswith("C")), res
    red["threshold"]["a00"] = 0.30
    assert table_claims(blue, red)["C8"]["verdict"] == "FAIL"


@needs_paper
def test_forward_and_watson_rules_accept_the_papers_own_curves():
    from smatrix_bootstrap.gtb.claims import c20, c21, paper_forward_as_analysis, paper_stage_as_analysis
    assert c20(paper_forward_as_analysis())["verdict"] == "PASS"
    before = {ni: paper_stage_as_analysis(ni, "before") for ni in (44, 46, 48)}
    after = {ni: paper_stage_as_analysis(ni, "after") for ni in (44, 46, 48)}
    assert c21(before, after)["verdict"] == "PASS"
