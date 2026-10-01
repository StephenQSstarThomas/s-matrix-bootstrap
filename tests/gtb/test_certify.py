"""Verification layer: V1 (Arb primal certificate) and V2 (dual bound) on small problems with known answers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from flint import arb

KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")
LEAF10 = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests/leaf_M10")


def test_verdict_of_a_ball():
    from smatrix_bootstrap.gtb.certify import verdict
    assert verdict(arb("1e-30", "1e-40")) == "certified_pass"
    assert verdict(arb("-1e-30", "1e-40")) == "certified_fail"
    assert verdict(arb("1e-40", "1e-30")) == "inconclusive"


def test_exact_ball_from_a_40_digit_string_contains_the_value():
    from smatrix_bootstrap.gtb.certify import ball
    b = ball("1234567890123456789012345678901234567890e-40")
    assert b.contains(arb("0.1234567890123456789012345678901234567890"))
    assert b.contains(arb("0.12345678901234567890123456789012345678905"))   # half an ulp away
    assert ball("0") == 0


needs_leaf10 = pytest.mark.skipif(not (LEAF10 / "solution.npz").exists(), reason="M=10 leaf missing")


@needs_leaf10
def test_primal_certifies_the_m10_leaf():
    import json
    from smatrix_bootstrap.gtb.certify import primal
    rep = primal(LEAF10, KS10)
    for fam, r in rep["families"].items():
        assert r["certified_fail"] == 0 and r["inconclusive"] == 0, (fam, r)
    sd = json.loads((LEAF10 / "report.json").read_text())["sdpb_result"]["primalObjective"]
    f1 = rep["objective_F1"]
    assert abs(float(f1["mid"]) - sd) < 1e-6 * max(1.0, abs(sd)) and f1["radius"] < 1e-25


@needs_leaf10
def test_primal_flags_a_violated_gram_block():
    from smatrix_bootstrap.gtb.certify import primal
    rep = primal(LEAF10, KS10, perturb={"ImF": (10 + 0, 50.0)})        # Im F of P1 at node 0 (M = 10) set huge
    assert rep["families"]["gram"]["certified_fail"] >= 1
