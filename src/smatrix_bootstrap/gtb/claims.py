"""Pre-registered claim verdicts (docs/gtb2403/PREREGISTRATION.md, amendment 1)."""
from __future__ import annotations

C1_PAPER = {44: -0.027229, 46: -0.028770, 48: -0.030397}       # Fig. 3 vector read-off (+-3e-5)
C1_TOL = 1e-4


def c1_verdict(lower: float, upper: float, paper: float, tol: float = C1_TOL) -> str:
    """Certified bracket [lower, upper] of F1max against the band paper +- tol: PASS only if inside, FAIL only if
    disjoint, otherwise INCONCLUSIVE (a verdict needs a bracket that decides the inequality)."""
    if lower > upper:
        raise ValueError("inverted bracket")
    lo, hi = paper - tol, paper + tol
    if lo <= lower and upper <= hi:
        return "PASS"
    if upper < lo or lower > hi:
        return "FAIL"
    return "INCONCLUSIVE"
