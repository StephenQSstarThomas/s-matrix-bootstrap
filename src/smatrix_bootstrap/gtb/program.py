"""The 2403 conic program as SDPB blocks.

A block is a symmetric matrix of affine forms in the SDPB variables [y0, y1, ...] (y0 is pinned to 1 by the
normalisation and carries the constants).  ``forms[(r, c)]`` with r <= c is the coefficient sequence of entry
(r, c); missing entries are identically zero.  The constraint is  M(y) = sum_n y_n M^(n)  >= 0.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Block:
    name: str
    size: int
    forms: dict = field(default_factory=dict)

    def entry(self, r: int, c: int):
        return self.forms.get((min(r, c), max(r, c)))
