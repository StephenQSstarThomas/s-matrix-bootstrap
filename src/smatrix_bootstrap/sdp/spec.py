"""Frozen finite-problem inputs, independent of every solver package."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

from . import constraints as C
from .grid import M_Q


@dataclass
class ModelSpec:
    M: int = 50
    L: int = 10
    chiral: bool = False
    chi_caliber: str = "chi-b"          # chi-a | chi-b | chi-c
    eps_chi: float = C.EPS_CHI_MAIN
    uv: bool = False                    # Gram + FESR + form-factor asymptotics
    uv_parts: tuple = ("gram", "fesr", "ff")   # for diagnosing infeasibility
    sr_caliber: str = "SR-b"            # SR-a | SR-b | SR-c
    eps_ff: float = C.EPS_FF
    ff_frozen_at_s0: bool = True        # historical replay; CLI uses the per-node (3.75) factors
    m_q: float = M_Q
    B: float | None = None              # legacy l2 arrow on (rho1, rho2); retired, see reg_*
    B_norm: str = "l2"                  # l2 (stronger, cheap) | l4 (pre-registered)
    # M-regularisation of He-Kruczenski 2103.11484 section 3: a cap on the
    # double-spectral-density node values, which unitarity at finitely many
    # nodes leaves unconstrained.  'linf' imposes |rho_{a,ij}| <= reg_bound as
    # two 1x1 SDPB blocks per value; the bound is fixed by the plateau rule of
    # that paper's section 3.3, never by any output curve.
    reg_norm: str | None = None         # None | 'linf'
    reg_bound: float | None = None      # Mreg for reg_norm
    # M50 spans: none 85.8, unclipped centrifugal 6.24, rownorm 1.25 decades.
    # The old 70-decade result clipped Lambda² at 1e-16; that bug is removed.
    cone_scaling: str = "rownorm"       # none | centrifugal | rownorm (all exact)
    sparsify: float = 0.0               # zero entries below this fraction of their row scale
    reduce_basis: bool = False          # SVD approximation to span(imposed rows)
    basis_tol: float = 1e-12
    operator_dps: int = 17               # CLI uses 40; 17 preserves historical matrix replay only
    scattering_prescription: str = 'mixed-pv'  # CLI uses sine-cardinal; historical loads remain explicit
    tag: str = ""
    disk_mask: object = None             # bool array over (wave, node): which disks to impose

    def key(self) -> str:
        return hashlib.sha256(repr(sorted(self.__dict__.items())).encode()).hexdigest()[:16]
