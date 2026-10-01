"""Encodings of the authors' double-spectral-density regulariser  ||rho / Mrho||_4 <= Mreg  (code reading).

R = Mreg * Mrho (= 1e2 * 3775 = 3.775e5) in the code reading, R = Mreg in the text reading.  With n = 3775 values,
    {|rho_i| <= R / n^(1/4)}  is inside  {||rho||_4 <= R}  which is inside  {|rho_i| <= R},
so 'linf-inner' and 'linf-relax' bracket the exact optimum from below and above (neither is the authors'
constraint; every number they produce is labelled with its encoding).  'none' drops the regulariser.
Sharper, equally cheap brackets: 'box-inner:FILE' takes per-value radii b_i with sum b_i^4 <= R^4 (any such box is
inside the l4 ball); 'linf-relax+cuts:FILE' adds supporting hyperplanes sum_i (rb_i/R)^3 (rho_i/R) <= 1 at points rb
with ||rb||_4 = R (valid on the ball by Hoelder, |sum rb^3 rho| <= ||rb^3||_{4/3} ||rho||_4 = R^3 R).
"""
from __future__ import annotations

import numpy as np


class _Linf:
    def __init__(self, prog, bound, label, cuts=None):
        self.prog, self.label, self.n_aux, self.cuts = prog, label, 0, cuts
        n = prog.groups["rho"].stop - prog.groups["rho"].start
        self.bound = np.broadcast_to(np.asarray(bound, dtype=float), (n,))

    def rho_forms(self):
        p = self.prog
        rho = np.arange(p.groups["rho"].start, p.groups["rho"].stop)
        if p.basis is None:                        # rho_i is the free variable y_{rho_i} (T0 is the pivot)
            return [p.unit(int(i)) for i in rho]
        out = np.zeros((rho.size, p.n_vars))
        out[:, 0] = p.a0[rho]
        out[:, 1:1 + p.n_af] = p.Ea[rho]
        return out

    def blocks(self):
        from .program import Block
        one, out = self.prog.const(1.0), []
        forms = self.rho_forms()
        for i, f in enumerate(forms):
            out.append(Block(f"reg:{self.label}:{i}:lo", 1, {(0, 0): one + f / self.bound[i]}))
            out.append(Block(f"reg:{self.label}:{i}:hi", 1, {(0, 0): one - f / self.bound[i]}))
        for k, g in enumerate(self.cuts if self.cuts is not None else []):
            out.append(Block(f"reg:cut:{k}", 1, {(0, 0): one - np.asarray(g) @ np.asarray(forms)}))
        return out


class _L4Exact(_Linf):
    """||rho||_4 <= T exactly:  u_i >= (rho_i/T)^2,  w_i >= u_i^2,  sum w_i <= 1  (2 auxiliaries per value)."""

    def __init__(self, prog, bound):
        super().__init__(prog, bound, "l4-exact")
        self.n_rho = prog.groups["rho"].stop - prog.groups["rho"].start
        self.n_aux = 2 * self.n_rho

    def blocks(self):
        from .program import Block
        p, T, n = self.prog, float(self.bound[0]), self.n_rho
        one, out = p.const(1.0), []
        last = p.const(1.0)
        for i, f in enumerate(self.rho_forms()):
            u, w = p.unit(p.o_aux + i), p.unit(p.o_aux + n + i)
            out.append(Block(f"reg:l4:{i}:u", 2, {(0, 0): one, (0, 1): f / T, (1, 1): u}))
            out.append(Block(f"reg:l4:{i}:w", 2, {(0, 0): one, (0, 1): u, (1, 1): w}))
            last = last - w
        out.append(Block("reg:l4:sum", 1, {(0, 0): last}))
        return out


class _NoReg:
    n_aux = 0

    def blocks(self):
        return []


def l4_bound(spec, n_rho: int) -> float:
    return spec.reg_bound * (n_rho if spec.reg_divide_by_mrho else 1)


def linf_data(kind: str, R: float, n_rho: int):
    """(radii, cut rows) of an l-infinity-type encoding as floats: |rho_i| <= radii_i and 1 - g_k . rho >= 0."""
    if kind == "linf-relax":
        return np.full(n_rho, R), None
    if kind == "linf-inner":
        return np.full(n_rho, R / n_rho ** 0.25), None
    if kind.startswith("box-inner:"):
        b = np.load(kind.split(":", 1)[1])
        if b.shape != (n_rho,) or np.any(b <= 0) or (b ** 4).sum() > R ** 4 * (1 + 1e-12):
            raise ValueError("box-inner radii must be positive with sum b^4 <= R^4")
        return b, None
    if kind.startswith("linf-relax+cuts:"):
        rb = np.atleast_2d(np.load(kind.split(":", 1)[1]))
        if rb.shape[1] != n_rho or np.any(np.abs((rb ** 4).sum(axis=1) ** 0.25 / R - 1) > 1e-9):
            raise ValueError("cut points must lie on the l4 sphere ||rb||_4 = R")
        return np.full(n_rho, R), (rb / R) ** 3 / R
    return None, None


def make(kind: str, spec, prog):
    n_rho = prog.groups["rho"].stop - prog.groups["rho"].start
    R = l4_bound(spec, n_rho)
    radii, cuts = linf_data(kind, R, n_rho)
    if radii is not None:
        return _Linf(prog, radii, kind.split(":")[0].replace("+cuts", ""), cuts=cuts)
    if kind == "l4-exact":
        return _L4Exact(prog, R)
    if kind == "none":
        return _NoReg()
    raise ValueError(f"unknown regulariser encoding {kind!r}")
