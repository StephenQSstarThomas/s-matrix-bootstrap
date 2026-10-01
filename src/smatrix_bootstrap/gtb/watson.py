"""Targets of the authors' Watson unitarisation step (2403 code c840933, second CVX block; DERIVATION section 6).

The step re-imposes every constraint, releases F0 and maximises  sum_n [Re t_n Re h~_n + Im t_n Im h~_n - Im h^_n]
over S0, P1, D0, S2, D2, F1 at the nodes below s0.  Its fixed point is h~ = t (saturated unitarity).  Targets:
S2, D2, F1 -> the previous h~;  S0, P1, D0 -> 2 Im F^ F / |F|^2 with Im F^ = Im F / Lambda_l (factor 2: the code,
consistent with the rank-one saturation of the rescaled Gram matrix; the text's (8.76) drops it).
"""
from __future__ import annotations

import numpy as np

from .program import CHANNELS, WATSON_WAVES


def targets_from_solution(prog, v) -> dict:
    """{wave: (t_re[n0], t_im[n0])} from a lifted float solution v of a previous leaf."""
    ks, M, n0 = prog.ks, prog.M, prog.n0
    amp, imf = v[prog.groups["amp"]], v[prog.groups["ImF"]]
    hre, him = (prog.Rre @ amp) / ks.Lam, (prog.Rim @ amp) / ks.Lam
    out = {}
    for w in WATSON_WAVES:
        if w in CHANNELS:
            k = CHANNELS.index(w)
            F = np.array([1 + ks.KFF3[k * M + j, k * M:(k + 1) * M] @ imf[k * M:(k + 1) * M] + 1j * imf[k * M + j]
                          for j in range(n0)])
            imfhat = imf[k * M:k * M + n0] / ks.Lambda3[k * n0:(k + 1) * n0]
            t = 2 * imfhat * F / np.abs(F) ** 2
            out[w] = (t.real.copy(), t.imag.copy())
        else:
            rows = [prog.row_index(w, j) for j in range(n0)]
            out[w] = (hre[rows].copy(), him[rows].copy())
    return out
