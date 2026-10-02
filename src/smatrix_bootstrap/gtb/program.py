"""The 2403 conic program as SDPB blocks (authors' code mainline, PLAN_SDPB_2403_ZH.md section 2).

A block is a symmetric matrix of affine forms in the SDPB variables [y0, y1, ...] (y0 is pinned to 1 by the
normalisation and carries the constants).  ``forms[(r, c)]`` with r <= c is the coefficient vector of entry (r, c);
missing entries are identically zero.  The constraint is  M(y) = sum_n y_n M^(n)  >= 0.

Raw variables v = [amplitude (T0, sigma1, sigma2, rho1, rho2) | Im F (S0, P1, D0 at all M nodes) | B-spline
coefficients (S0, P1, D0, 40 each)].  Equalities are eliminated, never written as two inequalities:
c_1 = 1 per channel, and F0 = 2 (f00(3) + 5 f02(3)) = F0fix through the T0 pivot (full space) or the largest
component in a reduced basis.  Every inequality block is scaled to O(1) entries by a positive congruence.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

CHANNELS = ("S0", "P1", "D0")                 # Gram / form factor / B-spline order (authors' code)
WAVES = {"S0": (0, 0), "D0": (0, 1), "S2": (1, 0), "D2": (1, 1), "P1": (2, 0), "F1": (2, 1)}  # (isospin block, wave)
CHI_POINTS = (0.5, 1.0, 1.5, 2.0)
WATSON_WAVES = ("S0", "P1", "D0", "S2", "D2", "F1")    # the six waves of the authors' Watson functional


@dataclass
class Block:
    name: str
    size: int
    forms: dict = field(default_factory=dict)

    def entry(self, r: int, c: int):
        return self.forms.get((min(r, c), max(r, c)))


class Program:
    def __init__(self, spec, ks, reg: str = "linf-relax", basis=None, watson_targets=None, objective_row=None,
                 extra_unitarity=None, psd_margin: float = 0.0):
        from . import regulariser
        M, nw, nb = spec.M, spec.L // 2, spec.n_b
        self.spec, self.ks, self.M, self.nw, self.nb = spec, ks, M, nw, nb
        self.exact = bool(getattr(ks, "exact_mode", False))   # Arb balls in object arrays (certification)
        man = ks.manifest
        if (int(man["M"]), int(man["L"]), man["matrices"]["BSpline3"]["shape"][1]) != (M, spec.L, 3 * nb):
            raise ValueError(f"spec (M, L, n_b) = {(M, spec.L, nb)} does not match the export "
                             f"(M, L, 3 n_b) = {(man['M'], man['L'], man['matrices']['BSpline3']['shape'][1])}")
        self.n0 = int(ks.manifest["n0"])
        n_amp = 1 + 2 * M + M * M + M * (M + 1) // 2
        self.groups = {"amp": slice(0, n_amp), "sig": slice(0, 2 * M + 1), "rho": slice(2 * M + 1, n_amp),
                       "ImF": slice(n_amp, n_amp + 3 * M), "c": slice(n_amp + 3 * M, n_amp + 3 * M + 3 * nb)}
        self.n_raw = n_amp + 3 * M + 3 * nb
        self.Rre = np.hstack([ks.Re_sig, ks.Re_rho])
        self.Rim = np.hstack([ks.Im_sig, ks.Im_rho])
        self.Fchi = np.hstack([ks.fchi_sig, ks.fchi_rho])
        fF = np.hstack([ks.fF_sig, ks.fF_rho])
        self.F0row = 2 * (fF[0] + 5 * fF[1])
        self.F1row = 2 * (3 * fF[2 * nw] + 7 * fF[2 * nw + 1])
        self.f0fix = self.num(spec.f0_ni) / 50 * self.num(spec.f0_max)
        self.watson = spec.objective == "watson"
        if self.watson and watson_targets is None:
            raise ValueError("the Watson objective needs targets (watson.targets_from_solution)")
        self.watson_targets = watson_targets
        self.objective_row = objective_row
        self.extra_unitarity = extra_unitarity            # midpoint rows (midpoints.load), WELLPOSED W3
        self.psd_margin = psd_margin                      # solve M - delta I >= 0: a strictly feasible point for V1
        if spec.objective.startswith("custom") and objective_row is None:
            raise ValueError("a custom objective needs objective_row over the raw variables")
        self._amplitude_map(n_amp, basis)
        self.c_first = np.array([0, nb, 2 * nb])
        # coefficients of B-splines with no node in their support enter no constraint: fixed to 0, removed
        self.c_dropped = np.setdiff1d(np.flatnonzero(~np.any(ks.BSpline3 != 0, axis=0)), self.c_first)
        self.c_free = np.setdiff1d(np.arange(3 * nb), np.r_[self.c_first, self.c_dropped])
        self.o_imf = 1 + self.n_af
        self.o_c = self.o_imf + 3 * M
        self.o_aux = self.o_c + len(self.c_free)
        self.reg = regulariser.make(reg, spec, self)
        self.n_vars = self.o_aux + self.reg.n_aux
        self._blocks = None

    # ------------------------------------------------------------------ numbers
    def num(self, x):
        """A spec constant: the float itself, or (exact mode) the Arb value of its decimal repr."""
        if not self.exact:
            return float(x)
        from flint import arb
        return arb(repr(x)) if isinstance(x, float) else arb(x)

    def zeros(self, shape):
        if not self.exact:
            return np.zeros(shape)
        from flint import arb
        return np.full(shape, arb(0), dtype=object)

    # ------------------------------------------------------------------ parametrisation
    def _amplitude_map(self, n_amp, basis):
        if self.watson:                                        # F0 is released in the Watson step
            self.pivot, self.a0 = None, self.zeros(n_amp)
            if basis is None:
                self.basis, self.n_af = None, n_amp
                self.Ea = None if self.exact else np.eye(n_amp)
            else:                                              # reduced basis: amp = V z, every z free
                self.basis, self.Ea, self.n_af = np.asarray(basis), np.asarray(basis), basis.shape[1]
            return
        if basis is None:
            p = 0                                              # T0: F0 row coefficient 5
            free = np.arange(1, n_amp)
            self.free_amp, self.coef_p = free, -self.F0row[free] / self.F0row[p]
            a0 = self.zeros(n_amp); a0[p] = self.f0fix / self.F0row[p]
            self.basis, self.pivot = None, ("amplitude", p)
            self.n_af = n_amp - 1
            if self.exact:
                self.a0, self.Ea = a0, None
                return
            Ea = np.zeros((n_amp, n_amp - 1))
            Ea[free, np.arange(n_amp - 1)] = 1.0
            Ea[p] = self.coef_p
        else:
            V = np.asarray(basis)
            g = self.F0row @ V
            r = int(np.argmax(np.abs(g)))
            k = V.shape[1]
            free = np.setdiff1d(np.arange(k), [r])
            Ez = np.zeros((k, k - 1)); Ez[free, np.arange(k - 1)] = 1.0; Ez[r] = -g[free] / g[r]
            z0 = np.zeros(k); z0[r] = self.f0fix / g[r]
            Ea, a0 = V @ Ez, V @ z0
            self.basis, self.pivot = V, ("basis", r)
        self.a0, self.Ea, self.n_af = a0, Ea, Ea.shape[1]

    def lift(self, y):
        """Raw variables v from SDPB's free variables y (without y0)."""
        y = np.asarray(y, dtype=float)
        v = np.zeros(self.n_raw)
        v[self.groups["amp"]] = self.a0 + self.Ea @ y[:self.n_af]
        v[self.groups["ImF"]] = y[self.n_af:self.n_af + 3 * self.M]
        c = np.zeros(3 * self.nb); c[self.c_first] = 1.0
        c[self.c_free] = y[self.o_c - 1:self.o_aux - 1]
        v[self.groups["c"]] = c
        return v

    def unlift(self, v):
        """Free variables y for a raw v that satisfies the eliminated equalities (aux variables set to 0)."""
        v = np.asarray(v, dtype=float)
        amp = v[self.groups["amp"]]
        if self.pivot is None:
            ya = amp.copy() if self.basis is None else self.basis.T @ amp
        elif self.basis is None:
            ya = np.delete(amp, self.pivot[1])
        else:
            ya = np.delete(self.basis.T @ amp, self.pivot[1])
        return np.r_[ya, v[self.groups["ImF"]], v[self.groups["c"]][self.c_free], np.zeros(self.reg.n_aux)]

    # ------------------------------------------------------------------ affine forms over [y0, y]
    def _empty(self):
        return self.zeros(self.n_vars)

    def amp_forms(self, rows):
        """Forms of linear functionals of the amplitude (rows: k x n_amp) -> k x n_vars."""
        rows = np.atleast_2d(rows)
        out = self.zeros((rows.shape[0], self.n_vars))
        if self.pivot is None:                     # no equality eliminated: the amplitude variables are free
            out[:, 1:1 + self.n_af] = rows if self.basis is None else rows @ self.Ea
        elif self.basis is None:                   # rank-one pivot update: identity on the free variables
            p = self.pivot[1]
            out[:, 0] = rows[:, p] * self.a0[p]
            out[:, 1:1 + self.n_af] = rows[:, self.free_amp] + rows[:, [p]] * self.coef_p[None, :]
        else:
            out[:, 0] = rows @ self.a0
            out[:, 1:1 + self.n_af] = rows @ self.Ea
        return out

    def unit(self, k):
        """Form of the free variable y_k (k >= 1 in [y0, y...] numbering)."""
        f = self._empty(); f[k] = self.num(1)
        return f

    def const(self, value):
        f = self._empty(); f[0] = self.num(value) if isinstance(value, (int, float)) else value
        return f

    def ff_forms(self, ch, j):
        """(Re F, Im F) of channel ch at node j:  Re F = 1 + (K + K0) Im F."""
        M, k = self.M, CHANNELS.index(ch)
        re = self.const(1.0)
        re[self.o_imf + k * M:self.o_imf + (k + 1) * M] = self.ks.KFF3[k * M + j, k * M:(k + 1) * M]
        im = self.unit(self.o_imf + k * M + j)
        return re, im

    def c_form(self, row):
        """Form of a linear functional of the 120 B-spline coefficients."""
        f = self._empty()
        f[0] = row[self.c_first].sum()
        f[self.o_c:self.o_aux] = row[self.c_free]
        return f

    def row_index(self, wave, j):
        b, w = WAVES[wave]
        return (b * self.nw + w) * self.M + j

    def chi_index(self, wave, a):
        b, w = WAVES[wave]
        return (b * self.nw + w) * len(CHI_POINTS) + a

    # ------------------------------------------------------------------ blocks
    def blocks(self):
        if self._blocks is None:
            self._blocks = (self._unitarity() + self._midpoints() + self._gram() + self._chiral() + self._fesr() + self._ff()
                            + self._cpos() + self.reg.blocks())
            if self.psd_margin:
                for blk in self._blocks:
                    for (r, c), f in list(blk.forms.items()):
                        if r == c:
                            g = f.copy()
                            g[0] = g[0] - self.num(self.psd_margin)
                            blk.forms[(r, c)] = g
        return self._blocks

    def _gram_rows(self):
        return {self.row_index(ch, j) for ch in CHANNELS for j in range(self.n0)}

    def _scaled_h(self):
        lam = self.ks.Lam[:, None]
        return self.amp_forms(self.Rre / lam), self.amp_forms(self.Rim / lam), self.amp_forms(self.Rim / lam ** 2)

    def _unitarity(self):
        hre, him, hhat = self._scaled_h()
        self._h = (hre, him, hhat)
        one, gram = self.const(1.0), self._gram_rows()
        return [Block(f"uni:{n}", 3, {(0, 0): 2 * hhat[n], (0, 1): hre[n], (0, 2): him[n], (1, 1): one, (2, 2): one})
                for n in range(self.Rre.shape[0]) if n not in gram]

    def _midpoints(self):
        """The node cone 2 Im h >= |h|^2 (Lambda-rescaled) at selected phi-midpoints (interpolant rows)."""
        m = self.extra_unitarity
        if m is None:
            return []
        lam = m["Lam"][:, None]
        hre, him, hhat = self.amp_forms(m["Re"] / lam), self.amp_forms(m["Im"] / lam), self.amp_forms(m["Im"] / lam ** 2)
        one, idx = self.const(1.0), m.get("rows", range(len(lam)))
        return [Block(f"mid:{int(r)}", 3, {(0, 0): 2 * hhat[n], (0, 1): hre[n], (0, 2): him[n], (1, 1): one, (2, 2): one})
                for n, r in enumerate(idx)]

    def _gram(self):
        hre, him, hhat = self._h
        one, out = self.const(1.0), []
        for k, ch in enumerate(CHANNELS):
            for j in range(self.n0):
                n = self.row_index(ch, j)
                fre, fim = self.ff_forms(ch, j)
                imfhat2 = 2 * fim / self.ks.Lambda3[k * self.n0 + j]
                rhat = self.c_form(self.ks.BSpline3[k * self.n0 + j])
                f = {(0, 0): one, (0, 1): hre[n], (0, 2): fre, (1, 1): 2 * hhat[n], (1, 2): imfhat2, (2, 2): rhat,
                     (0, 4): -him[n], (0, 5): -fim, (1, 3): him[n], (2, 3): fim,
                     (3, 3): one, (3, 4): hre[n], (3, 5): fre, (4, 4): 2 * hhat[n], (4, 5): imfhat2, (5, 5): rhat}
                out.append(Block(f"gram:{ch}:{j}", 6, f))
        return out

    @staticmethod
    def arrow(name, forms, eps, one):
        """||r|| <= eps  as  [[1, r^T/eps], [r/eps, I]] >= 0."""
        f = {(0, 0): one}
        for k, r in enumerate(forms, start=1):
            f[(0, k)] = r / eps
            f[(k, k)] = one
        return Block(name, len(forms) + 1, f)

    def _chiral(self):
        F, eps, one = self.Fchi, self.num(self.spec.eps_chi), self.const(1.0)
        r01 = [F[self.chi_index("S0", a)] - 3 * (2 * s - 1) / (s - 4) * F[self.chi_index("P1", a)]
               for a, s in enumerate(CHI_POINTS)]
        r21 = [F[self.chi_index("S2", a)] - 3 * (2 - s) / (s - 4) * F[self.chi_index("P1", a)]
               for a, s in enumerate(CHI_POINTS)]
        hi = [F[self.chi_index(w, a)] for w in self.spec.chi_second_norm_waves for a in range(len(CHI_POINTS))]
        return [self.arrow("chi:ratio", self.amp_forms(np.array(r01 + r21)), eps, one),
                self.arrow("chi:higher", self.amp_forms(np.array(hi)), eps, one)]

    def _fesr(self):
        ks, n0, out = self.ks, self.n0, []
        for k, (ch, ints) in enumerate(zip(CHANNELS, (ks.intS0, ks.intP1, ks.intD0))):
            nsr = ints.shape[0]                                    # n_SR moments per channel (3 in the paper)
            g = ints @ ks.BSpline3[k * n0:(k + 1) * n0]            # n_SR x 120 (block diagonal in channels)
            w = [self.c_form(g[i]) - self.const(ks.SVZSR[nsr * k + i]) for i in range(nsr)]
            out.append(self.arrow(f"fesr:{ch}", w, self.num(self.spec.eps_sr[k]), self.const(1.0)))
        return out

    def _ff(self):
        sp, x, one, out, n = self.spec, self.ks.x, self.const(1.0), [], self.num
        caps = {"S0": lambda j: n(sp.r[0]), "P1": lambda j: n(sp.r[1]) * n(sp.fp1asym) / x[j],
                "D0": lambda j: n(sp.r[2]) * n(sp.fd0asym) / x[j]}
        for ch in CHANNELS:
            for j in range(self.n0, self.M):
                fre, fim = self.ff_forms(ch, j)
                cap = caps[ch](j)
                out.append(Block(f"ff:{ch}:{j}", 3, {(0, 0): one, (0, 1): fre / cap, (0, 2): fim / cap,
                                                     (1, 1): one, (2, 2): one}))
        return out

    def _cpos(self):
        out = []
        for k, i in enumerate(self.c_free):
            out.append(Block(f"cpos:{i}", 1, {(0, 0): self.unit(self.o_c + k)}))
        return out

    def objective(self):
        if self.spec.objective.startswith("custom"):
            row = np.asarray(self.objective_row, dtype=float) * (-1 if self.spec.objective == "custom_min" else 1)
            f = self.amp_forms(row[self.groups["amp"]])[0]
            f[self.o_imf:self.o_imf + 3 * self.M] += row[self.groups["ImF"]]
            return f + self.c_form(row[self.groups["c"]])
        if not self.watson:
            sign = -1 if self.spec.objective == "F1min" else 1
            return self.amp_forms(sign * self.F1row)[0]
        lam = self.ks.Lam
        row = self.zeros(self.Rre.shape[1])
        for w in WATSON_WAVES:
            tre, tim = self.watson_targets[w]
            for j in range(self.n0):
                n = self.row_index(w, j)
                row = row + (self.num(float(tre[j])) * self.Rre[n] + self.num(float(tim[j])) * self.Rim[n]) / lam[n] \
                    - self.Rim[n] / lam[n] ** 2
        return self.amp_forms(row)[0]
