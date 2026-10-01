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
