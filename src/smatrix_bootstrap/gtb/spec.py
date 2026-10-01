"""Every parameter and reading switch of one 2403 finite-model calculation.

Defaults are the authors' code (c840933) mainline, see PLAN_SDPB_2403_ZH.md section 2;
text readings are controls selected through ``reading`` and the individual fields.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class GtbSpec:
    M: int = 50
    L: int = 20
    nu0: int = -20
    n_b: int = 40
    eps_chi: float = 2e-3
    eps_sr: tuple = (1e-7, 6e-6, 5e-6)
    r: tuple = (0.05, 2.0, 6.0)
    fp1asym: float = 6.87076
    fd0asym: float = 20.6123
    reg_norm: str = "l4"
    reg_bound: float = 1e2
    reg_divide_by_mrho: bool = True
    chi_second_norm_waves: tuple = ("D0", "D2", "F1")
    watson_factor: int = 2
    f0_ni: int = 46
    f0_max: float = 0.163108
    objective: str = "F1max"
    reading: str = "code"

    def key(self) -> str:
        """First 16 hex digits of the sha256 of the fields, serialised with sorted keys."""
        text = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode()).hexdigest()[:16]
