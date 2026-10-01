"""Low-energy parameters as linear functionals of the raw variables (for [min, max] over the feasible set, B7).

a, b: the pre-registered threshold fit Re f/k^(2l) = a + b k^2 + c k^4 on 12 k^2 in (0, 0.25] is linear in Re f, which is
linear in the amplitude; the fit weights times the interpolant rows give one row each.  <r^2>: 6 F'(0) with the
midpoint F'(0) = sum_m (Delta/pi) x'_m Im F_m / x_m^2, linear in Im F.  Rows are over the raw vector v of program.py.
"""
from __future__ import annotations

import numpy as np

from .analysis import HBARC, M_PI_MEV
from .leafanalysis import WAVE_IL


def lec_rows(prog) -> dict:
    from flint import arb
    from .arbkernels import Grid
    from .evaluate import f_row
    grid = Grid(prog.M, prec_bits=192)
    k2 = np.linspace(0.25 / 12, 0.25, 12)
    out = {}
    for w, key in (("S0", "00"), ("S2", "20"), ("P1", "11"), ("D0", "02"), ("D2", "22")):
        I, ell = WAVE_IL[w]
        F = np.array([[float(c.real.mid()) for c in f_row(grid, I, ell, ("cut", arb(repr(float(4 + 4 * q)))), "interp")]
                      for q in k2])                                    # 12 x n_amp: Re f rows
        X = np.vstack([np.ones_like(k2), k2, k2 ** 2]).T
        W = np.linalg.pinv(X) / k2[None, :] ** ell                     # (a, b, c) = W @ Re f
        for name, wrow in (("a", W[0]), ("b", W[1])):
            row = np.zeros(prog.n_raw)
            row[prog.groups["amp"]] = wrow @ F
            out[name + key] = row
    M, d, ks = prog.M, np.pi / prog.M, prog.ks
    for k, key in enumerate(("r2S", "r2V", "r2T")):
        row = np.zeros(prog.n_raw)
        sl = slice(prog.groups["ImF"].start + k * M, prog.groups["ImF"].start + (k + 1) * M)
        row[sl] = 6 * d / np.pi * ks.xp / ks.x ** 2 * (HBARC / M_PI_MEV) ** 2
        out[key] = row
    return out
