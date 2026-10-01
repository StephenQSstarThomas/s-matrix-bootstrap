"""Encodings of the authors' double-spectral-density regulariser  ||rho / Mrho||_4 <= Mreg  (code reading).

R = Mreg * Mrho (= 1e2 * 3775 = 3.775e5) in the code reading, R = Mreg in the text reading.  With n = 3775 values,
    {|rho_i| <= R / n^(1/4)}  is inside  {||rho||_4 <= R}  which is inside  {|rho_i| <= R},
so 'linf-inner' and 'linf-relax' bracket the exact optimum from below and above (neither is the authors'
constraint; every number they produce is labelled with its encoding).  'none' drops the regulariser.
"""
from __future__ import annotations

import numpy as np


class _Linf:
    def __init__(self, prog, bound, label):
        self.prog, self.bound, self.label, self.n_aux = prog, float(bound), label, 0

    def rho_forms(self):
        p = self.prog
        rho = np.arange(p.groups["rho"].start, p.groups["rho"].stop)
        out = np.zeros((rho.size, p.n_vars))
        out[:, 0] = p.a0[rho]
        out[:, 1:1 + p.n_af] = p.Ea[rho]
        return out

    def blocks(self):
        from .program import Block
        one, out = self.prog.const(1.0), []
        for i, f in enumerate(self.rho_forms()):
            out.append(Block(f"reg:{self.label}:{i}:lo", 1, {(0, 0): one + f / self.bound}))
            out.append(Block(f"reg:{self.label}:{i}:hi", 1, {(0, 0): one - f / self.bound}))
        return out


class _NoReg:
    n_aux = 0

    def blocks(self):
        return []


def l4_bound(spec, n_rho: int) -> float:
    return spec.reg_bound * (n_rho if spec.reg_divide_by_mrho else 1)


def make(kind: str, spec, prog):
    n_rho = prog.groups["rho"].stop - prog.groups["rho"].start
    R = l4_bound(spec, n_rho)
    if kind == "linf-relax":
        return _Linf(prog, R, "linf-relax")
    if kind == "linf-inner":
        return _Linf(prog, R / n_rho ** 0.25, "linf-inner")
    if kind == "none":
        return _NoReg()
    raise ValueError(f"unknown regulariser encoding {kind!r}")
