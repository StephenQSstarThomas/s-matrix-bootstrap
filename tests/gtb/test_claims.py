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
