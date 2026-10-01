"""Verification of 2403 leaves.

V1 (primal, rigorous at the collocation nodes): every constraint is recomputed in Arb from the exact kernel export
(each 40-digit entry becomes a ball of radius one unit in its last digit) at a point lifted from SDPB's printed free
variables, with the eliminated equalities solved exactly in Arb.  A strictly positive lower bound of every slack is a
proof that this point is feasible for the finite problem, hence its objective is a rigorous lower bound of the max.
V2 (dual bound) is added below; it gives the rigorous upper bound.
"""
from __future__ import annotations

from flint import arb


def ball(text: str) -> arb:
    """Decimal ``<int>e<exp>`` (or ``0``) as a ball containing the exact value it was rounded from."""
    if text == "0":
        return arb(0)
    mant, _, exp = text.partition("e")
    e = int(exp) if exp else 0
    return arb(text) + arb(0, arb(10) ** e)


def verdict(x: arb) -> str:
    if x > 0:
        return "certified_pass"
    if x < 0:
        return "certified_fail"
    return "inconclusive"


# ---------------------------------------------------------------------------------------------- V1 primal
import json
from pathlib import Path

import numpy as np
from flint import acb, arb_mat, ctx

from .program import CHANNELS, CHI_POINTS, WAVES


def _vec(ks, name):
    return [ball(t) for t in ks.exact(name)]


def _mat(ks, name, cols=None):
    shape = ks.manifest["matrices"][name]["shape"]
    flat = ks.exact(name)
    r, c = shape
    return [[ball(flat[i * c + j]) for j in range(c)] for i in range(r)]


def _dot(row, v):
    return arb.dot_sum(row, v) if hasattr(arb, "dot_sum") else sum((a * b for a, b in zip(row, v)), arb(0))


class _Family:
    def __init__(self):
        self.counts = {"certified_pass": 0, "certified_fail": 0, "inconclusive": 0}
        self.min_lower, self.where = None, None

    def add(self, slack, where):
        self.counts[verdict(slack)] += 1
        lo = float(slack.lower()) if slack.is_finite() else float("-inf")
        if self.min_lower is None or lo < self.min_lower:
            self.min_lower, self.where = lo, where

    def report(self):
        return {**self.counts, "min_lower_slack": self.min_lower, "at": self.where}


def _lift(prog_info, ks, y_text, basis, rows, perturb):
    """Raw variables as Arb balls; the F0 equality and c_1 = 1 hold exactly (pivot solved in Arb)."""
    from .spec import GtbSpec
    spec, M, nb = prog_info["spec"], prog_info["M"], prog_info["nb"]
    y = [arb(str(t)) for t in y_text]
    n_amp = 1 + 2 * M + M * M + M * (M + 1) // 2
    F0row = rows["F0"]
    f0fix = arb(spec["f0_ni"]) / 50 * arb(repr(spec["f0_max"]))
    if prog_info["pivot"] is None:            # Watson leaf: F0 released, no equality to solve
        if basis is None:
            amp, n_af = y[:n_amp], n_amp
        else:
            V = [[arb(repr(float(x))) for x in row] for row in basis]
            k = len(V[0])
            amp, n_af = [_dot(V[i], y[:k]) for i in range(n_amp)], k
    elif basis is None:
        rest = y[:n_amp - 1]
        t0 = (f0fix - _dot(F0row[1:], rest)) / F0row[0]
        amp = [t0] + rest
        n_af = n_amp - 1
    else:
        V = [[arb(repr(float(x))) for x in row] for row in basis]
        k = len(V[0])
        g = [_dot([V[i][j] for i in range(n_amp)], F0row) for j in range(k)]
        r = int(prog_info["pivot"][1])
        free = [j for j in range(k) if j != r]
        z = [arb(0)] * k
        for jj, j in enumerate(free):
            z[j] = y[jj]
        z[r] = (f0fix - sum((g[j] * z[j] for j in free), arb(0))) / g[r]
        amp = [_dot(V[i], z) for i in range(n_amp)]
        n_af = k - 1
    imf = y[n_af:n_af + 3 * M]
    c = [arb(0)] * (3 * nb)
    empty = prog_info["empty_c"]           # exactly-zero B-spline columns, recomputed from the exact export
    free_c = [i for i in range(3 * nb) if i not in (0, nb, 2 * nb) and i not in empty]
    for i in (0, nb, 2 * nb):
        c[i] = arb(1)
    for k_, i in enumerate(free_c):
        c[i] = y[n_af + 3 * M + k_]
    for name, (idx, val) in (perturb or {}).items():
        {"ImF": imf, "amp": amp, "c": c}[name][idx] = arb(val)
    return amp, imf, c


def primal(leaf_dir, kernels_dir, perturb=None) -> dict:
    """V1 for one leaf: per-family verdicts, the minimum certified slack, and the objective as a ball."""
    from .kernels import load
    ctx.prec = 256
    leaf, ks = Path(leaf_dir), load(kernels_dir)
    rep = json.loads((leaf / "report.json").read_text())
    spec = rep["spec"]
    M, nb, nw, n0 = spec["M"], spec["n_b"], spec["L"] // 2, int(ks.manifest["n0"])
    y_text = list(np.load(leaf / "solution.npz")["y_text"])
    basis = np.load(leaf / "basis.npy") if (leaf / "basis.npy").exists() else None
    cat = lambda a, b: [ra + rb for ra, rb in zip(_mat(ks, a), _mat(ks, b))]
    Rre, Rim, Fchi, fF = cat("Re_sig", "Re_rho"), cat("Im_sig", "Im_rho"), cat("fchi_sig", "fchi_rho"), cat("fF_sig", "fF_rho")
    F0row = [2 * (a + 5 * b) for a, b in zip(fF[0], fF[1])]
    F1row = [2 * (3 * a + 7 * b) for a, b in zip(fF[2 * nw], fF[2 * nw + 1])]
    kff, bs, svz = _mat(ks, "KFF3"), _mat(ks, "BSpline3"), _vec(ks, "SVZSR")
    empty_c = {i for i in range(3 * nb) if all(r[i] == 0 for r in bs)}
    info = {"spec": spec, "M": M, "nb": nb, "pivot": rep["pivot"], "empty_c": empty_c}
    amp, imf, c = _lift(info, ks, y_text, basis, {"F0": F0row}, perturb)
    lam, lam3, x = _vec(ks, "Lam"), _vec(ks, "Lambda3"), _vec(ks, "x")
    hre = [_dot(r, amp) for r in Rre]
    him = [_dot(r, amp) for r in Rim]
    fam = {k: _Family() for k in ("unitarity", "gram", "chiral", "fesr", "ff", "c_positive", "regulariser")}
    row = lambda w, j: (WAVES[w][0] * nw + WAVES[w][1]) * M + j
    gram_rows = {row(ch, j) for ch in CHANNELS for j in range(n0)}
    for n in range(len(hre)):
        if n not in gram_rows:
            fam["unitarity"].add((2 * him[n] - hre[n] ** 2 - him[n] ** 2) / lam[n] ** 2, f"row {n}")
    F = []
    for k, ch in enumerate(CHANNELS):
        F.append([acb(1 + _dot(kff[k * M + j][k * M:(k + 1) * M], imf[k * M:(k + 1) * M]), imf[k * M + j])
                  for j in range(M)])
    for k, ch in enumerate(CHANNELS):
        for j in range(n0):
            n = row(ch, j)
            ht = acb(hre[n], him[n]) / lam[n]
            hh = 2 * him[n] / lam[n] ** 2
            e = 2 * F[k][j].imag / lam3[k * n0 + j]
            rh = _dot(bs[k * n0 + j], c)
            f = F[k][j]
            m2 = hh - abs(ht) ** 2
            det = (hh * rh - e * e) - ht * (ht.conjugate() * rh - e * f.conjugate()) \
                + f * (ht.conjugate() * e - hh * f.conjugate())
            fam["gram"].add(m2, f"{ch}:{j}:m2")
            fam["gram"].add(det.real, f"{ch}:{j}:det")
    chi = lambda w, a: (WAVES[w][0] * nw + WAVES[w][1]) * len(CHI_POINTS) + a
    fv = [_dot(r, amp) for r in Fchi]
    r1 = [fv[chi("S0", a)] - 3 * (2 * arb(repr(s)) - 1) / (arb(repr(s)) - 4) * fv[chi("P1", a)] for a, s in enumerate(CHI_POINTS)]
    r1 += [fv[chi("S2", a)] - 3 * (2 - arb(repr(s))) / (arb(repr(s)) - 4) * fv[chi("P1", a)] for a, s in enumerate(CHI_POINTS)]
    r2 = [fv[chi(w, a)] for w in spec["chi_second_norm_waves"] for a in range(len(CHI_POINTS))]
    eps = arb(repr(spec["eps_chi"]))
    fam["chiral"].add(eps ** 2 - sum((v * v for v in r1), arb(0)), "ratio")
    fam["chiral"].add(eps ** 2 - sum((v * v for v in r2), arb(0)), "higher")
    for k, (ch, nm) in enumerate(zip(CHANNELS, ("intS0", "intP1", "intD0"))):
        ints = _mat(ks, nm)
        rhat = [_dot(bs[k * n0 + j], c) for j in range(n0)]
        w = [_dot(ints[i], rhat) - svz[3 * k + i] for i in range(3)]
        e = arb(repr(spec["eps_sr"][k]))
        fam["fesr"].add(e ** 2 - sum((v * v for v in w), arb(0)), ch)
    r = [arb(repr(v)) for v in spec["r"]]
    caps = {"S0": lambda j: r[0], "P1": lambda j: r[1] * arb(repr(spec["fp1asym"])) / x[j],
            "D0": lambda j: r[2] * arb(repr(spec["fd0asym"])) / x[j]}
    for k, ch in enumerate(CHANNELS):
        for j in range(n0, M):
            fam["ff"].add(caps[ch](j) ** 2 - abs(F[k][j]) ** 2, f"{ch}:{j}")
    for i in range(3 * nb):
        if i not in (0, nb, 2 * nb) and i not in empty_c:
            fam["c_positive"].add(c[i], f"c{i}")
    rho = amp[1 + 2 * M:]
    R = arb(repr(spec["reg_bound"])) * (len(rho) if spec["reg_divide_by_mrho"] else 1)
    l4 = sum((v ** 4 for v in rho), arb(0))
    if rep["regulariser"].startswith("linf"):
        bound = R if rep["regulariser"] == "linf-relax" else R / arb(len(rho)).root(4)
        for i, v in enumerate(rho):
            fam["regulariser"].add(bound - abs(v), f"rho{i}")
    l4_verdict = verdict(R ** 4 - l4)
    f1 = _dot(F1row, amp)
    return {"families": {k: v.report() for k, v in fam.items()},
            "objective_F1": {"mid": f1.mid().str(30, radius=False), "radius": float(f1.rad())},
            "l4_feasible": l4_verdict, "rho_l4": float(l4.root(4).mid()), "l4_bound": float(R.mid()),
            "F0_check": float((_dot(F0row, amp)).mid())}


# V2 lives in dual.py (size limit); re-exported here.
from .dual import ExactKernels, certify_leaf, dual_bound, psd_shift, read_dual_matrices  # noqa: E402,F401
