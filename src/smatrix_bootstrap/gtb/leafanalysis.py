"""Everything the 2403 claims need from one leaf (PREREGISTRATION amendment 1).  Units m_pi = 1 unless 'MeV'/'GeV'."""
from __future__ import annotations

import json
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

from . import analysis as A

WAVE_IL = {"S0": (0, 0), "D0": (0, 2), "S2": (2, 0), "D2": (2, 2), "P1": (1, 1), "F1": (1, 3)}
KAPPA = {"S0": lambda s: 3 * np.sqrt(s - 4) / (256 * np.pi ** 5 * np.sqrt(s)),
         "P1": lambda s: (s - 4) ** 1.5 / (384 * np.pi ** 5 * np.sqrt(s)),
         "D0": lambda s: (s - 4) ** 2.5 / (2560 * np.pi ** 5 * np.sqrt(s))}
GEV = A.M_PI_MEV / 1000


def _load(leaf, kernels_dir):
    from .kernels import load
    from .program import Program
    from .spec import GtbSpec
    rep = json.loads((Path(leaf) / "report.json").read_text())
    sp = {f.name: (tuple(rep["spec"][f.name]) if isinstance(rep["spec"][f.name], list) else rep["spec"][f.name])
          for f in fields(GtbSpec)}
    spec = replace(GtbSpec(**sp), objective="F1max")          # analysis needs the rows, not the objective
    ks = load(kernels_dir)
    return rep, ks, Program(spec, ks, reg="none"), np.load(Path(leaf) / "solution.npz")["v"]


def _form_factors(prog, v):
    ks, M = prog.ks, prog.M
    imf = v[prog.groups["ImF"]]
    out = {}
    for k, ch in enumerate(("S0", "P1", "D0")):
        im = imf[k * M:(k + 1) * M]
        out[ch] = 1 + ks.KFF3[k * M:(k + 1) * M, k * M:(k + 1) * M] @ im + 1j * im
    return out


def _threshold(prog, amp, grid):
    """a, b from the interpolant on 12 k^2 in (0, 0.25] (mainline) and from the first 4 nodes (control)."""
    from flint import arb
    from .evaluate import f_row, value
    k2 = np.linspace(0.25 / 12, 0.25, 12)
    out, ctrl = {}, {}
    x = prog.ks.x
    for w, key in (("S0", "00"), ("S2", "20"), ("P1", "11"), ("D0", "02"), ("D2", "22")):
        I, ell = WAVE_IL[w]
        ref = [float(value(f_row(grid, I, ell, ("cut", arb(repr(float(4 + 4 * q)))), "interp"), amp).real.mid()) for q in k2]
        a, b = A.threshold_fit(k2, np.array(ref), ell)
        out["a" + key] = a
        if ell < 2:
            out["b" + key] = b
        rows = [prog.row_index(w, j) for j in range(4)]
        kn = (x[:4] - 4) / 4
        ctrl["a" + key] = A.threshold_fit(kn, (prog.Rre[rows] @ amp) / (np.pi * np.sqrt((x[:4] - 4) / x[:4])), ell)[0]
    return out, ctrl


def threshold_below(grid, amp, n=12, lo=2.0, hi=3.8) -> dict:
    """Threshold parameters from the interpolant below threshold (PREREGISTRATION amendment 8): f_l(s) on n points
    s in [lo, hi] (z inside the disk, where the node-alternating artefacts of the discrete K decay as |z|^M), fit of
    f/k^(2l) = a + b k^2 + c k^4 with k^2 = (s-4)/4 < 0, read at k^2 = 0."""
    from flint import arb
    from .evaluate import f_row, value
    s_pts = np.linspace(lo, hi, n)
    out = {}
    for w, key in (("S0", "00"), ("S2", "20"), ("P1", "11"), ("D0", "02"), ("D2", "22")):
        I, ell = WAVE_IL[w]
        f = [float(value(f_row(grid, I, ell, ("point", arb(repr(float(s)))), "interp"), amp).real.mid()) for s in s_pts]
        a, b = A.threshold_fit((s_pts - 4) / 4, np.array(f), ell)
        out["a" + key] = a
        if ell < 2:
            out["b" + key] = b
    return out


def s_abs(grid, amp, wave, where, prescription="interp") -> float:
    """|S^I_l| = |1 + i pi sqrt((s-4)/s) f^I_l(s)| at a node or a real point on the cut."""
    from flint import acb
    from .evaluate import _s_and_z, f_row, value
    I, ell = WAVE_IL[wave]
    s = _s_and_z(grid, where)[0]
    f = value(f_row(grid, I, ell, where, prescription), amp)
    return float(abs(1 + acb(0, 1) * acb.pi() * ((s - 4) / s).sqrt() * f).mid())


def between_nodes(grid, amp, n0, prescription="interp") -> dict:
    """Unitarity between the nodes (Phase B, Task B7): |S| at phi = (k+1) pi/M, halfway in angle between nodes k
    and k+1, for the six waves below s0.  The program imposes |S| <= 1 only at the nodes."""
    from flint import arb
    M, out = grid.M, {}
    pts = [arb(grid.nu0) + (8 - 2 * arb(grid.nu0)) / (1 + (arb.pi() * (k + 1) / M).cos()) for k in range(n0 - 1)]
    for w in WAVE_IL:
        a = [s_abs(grid, amp, w, ("cut", s), prescription) for s in pts]
        out[w] = {"s": [float(s.mid()) for s in pts], "abs_S": a, "max_abs_S": max(a)}
    return out


def omitted_waves(amp, big, L_leaf, include_kept=False) -> dict:
    """|S| = |1 + i h| of the partial waves l >= L_leaf at every node, from a larger-L export of the same kernels
    (the amplitude variables depend on M only).  Unitarity is imposed only for l < L_leaf (2309 rule: eta <= 1.02)."""
    M, nw = int(big.manifest["M"]), int(big.manifest["L"]) // 2
    re, im = np.hstack([big.Re_sig, big.Re_rho]), np.hstack([big.Im_sig, big.Im_rho])
    out, kept = {"eta_max": {}, "eta": {}}, {}
    names = {(0, 0): "S0", (0, 2): "D0", (2, 0): "S2", (2, 2): "D2", (1, 1): "P1", (1, 3): "F1"}
    for b, I in enumerate((0, 2, 1)):
        for w in range(nw):
            ell = 2 * w + (1 if I == 1 else 0)
            rows = slice((b * nw + w) * M, (b * nw + w + 1) * M)
            eta = np.abs(1 + 1j * (re[rows] @ amp + 1j * (im[rows] @ amp)))
            if ell >= L_leaf:
                out["eta_max"][f"I{I}_l{ell}"] = float(eta.max())
                out["eta"][f"I{I}_l{ell}"] = eta.tolist()
            elif (I, ell) in names:
                kept[names[(I, ell)]] = eta.tolist()
    if include_kept:
        out["kept"] = kept
    return out


def analyze(leaf, kernels_dir, pole=True, omitted_kernels=None) -> dict:
    from .arbkernels import Grid
    rep, ks, prog, v = _load(leaf, kernels_dir)
    M, n0, x = prog.M, prog.n0, ks.x
    amp = v[prog.groups["amp"]]
    E = np.sqrt(x) * GEV
    h = prog.Rre @ amp + 1j * (prog.Rim @ amp)
    res = {"leaf": str(leaf), "n0": n0, "E_GeV": E[:n0].tolist(), "phases": {}, "crossings_GeV": {}}
    for w in WAVE_IL:
        S = 1 + 1j * h[[prog.row_index(w, j) for j in range(n0)]]
        d, eta = A.phases(S)
        res["phases"][w] = {"delta_deg": d.tolist(), "eta": eta.tolist()}
    for w in ("P1", "D0", "S0"):
        res["crossings_GeV"][w] = A.crossing(E[:n0], np.array(res["phases"][w]["delta_deg"]), 90.0)
    F = _form_factors(prog, v)
    res["FF_abs"] = {ch: np.abs(f).tolist() for ch, f in F.items()}
    d = np.pi / M
    imf = v[prog.groups["ImF"]]
    fp = {ch: float(np.sum(d / np.pi * ks.xp * imf[k * M:(k + 1) * M] / x ** 2))
          for k, ch in enumerate(("S0", "P1", "D0"))}
    res["radii_fm2"] = {"S": A.radius_fm2(fp["S0"]), "V": A.radius_fm2(fp["P1"]), "T": A.radius_fm2(fp["D0"])}
    grid = Grid(M, prec_bits=192)
    res["threshold"], res["threshold_control"] = _threshold(prog, amp, grid)
    res["threshold_below"] = threshold_below(grid, amp)
    res["between_nodes"] = between_nodes(grid, amp, n0)
    from .watson import residual, targets_from_solution
    res["watson_self_residual"] = residual(prog, v, targets_from_solution(prog, v))
    if (Path(leaf) / "watson_targets.npz").exists():
        z = np.load(Path(leaf) / "watson_targets.npz")
        seed = {w: (z[f"{w}_re"], z[f"{w}_im"]) for w in {k.rsplit("_", 1)[0] for k in z.files}}
        res["watson_seed_residual"] = residual(prog, v, seed)
    res["lbar"] = A.lbar(res["threshold"]["a02"], res["threshold"]["a22"], res["radii_fm2"]["S"], res["radii_fm2"]["V"])
    s = x[:n0]
    win = lambda lo, hi: (E[:n0] >= lo) & (E[:n0] <= hi)
    fits = {}
    w1 = win(0.5, 1.1)
    for form in ("bw", "gs"):
        m, g = A.fit_rho(s[w1], np.abs(F["P1"][:n0][w1]), form, guess=(0.80 / GEV, 0.15 / GEV))
        fits[f"rho_{form}"] = {"m_MeV": m * A.M_PI_MEV, "Gamma_MeV": g * A.M_PI_MEV, "g_rho_pi_pi": A.g_rho_pi_pi(m, g)}
    w2 = win(0.5, 1.83)
    p = A.fit_two_bw(s[w2], KAPPA["P1"](s[w2]) * np.abs(F["P1"][:n0][w2]) ** 2,
                     guess=(0.82 / GEV, 0.14 / GEV, 1.6 / GEV, 0.25 / GEV, 1.3, 0.27))
    fits["two_bw"] = dict(zip(("m1_MeV", "G1_MeV", "m2_MeV", "G2_MeV"), [q * A.M_PI_MEV for q in p[:4]]),
                          A1=p[4], A2=p[5])
    w3 = win(1.0, 1.83)
    m, g = A.fit_rho(s[w3], np.abs(F["D0"][:n0][w3]), "bw", guess=(1.37 / GEV, 0.17 / GEV))
    fits["f2_bw"] = {"m_MeV": m * A.M_PI_MEV, "Gamma_MeV": g * A.M_PI_MEV, "g_f2_pi_pi": A.g_f2_pi_pi(m, g)}
    w4 = win(0.6, 1.0)
    mr = fits["rho_gs"]["m_MeV"] / A.M_PI_MEV
    integ = np.sum((d * ks.xp[:n0] * KAPPA["P1"](s) * np.abs(F["P1"][:n0]) ** 2 / s)[w4])
    fits["four_pi_over_g_rho2"] = float(2 * (2 * np.pi) ** 4 * integ / mr ** 2)
    res["fits"] = fits
    f_all = {I: {} for I in (0, 1, 2)}
    for (b, nw) in ((0, 0), (1, 2), (2, 1)):
        for wi in range(prog.nw):
            ell = 2 * wi + (1 if nw == 1 else 0)
            rows = [(b * prog.nw + wi) * M + j for j in range(n0)]
            f_all[nw][ell] = h[rows] / (np.pi * np.sqrt((s - 4) / s))
    res["forward"] = {k: {"re": v_.real.tolist(), "im": v_.imag.tolist()} for k, v_ in A.forward(f_all).items()}
    if omitted_kernels:
        from .kernels import load
        res["omitted_waves"] = omitted_waves(amp, load(omitted_kernels), prog.spec.L)
    res["pole"] = {}
    if pole:
        from .pole import find_zero, s_matrix_p1, sqrt_s_mev
        for presc in ("interp", "midpoint"):
            try:
                z = find_zero(s_matrix_p1(grid, amp, presc), (0.83 / GEV - 0.04j / GEV) ** 2)
                r = sqrt_s_mev(z)
                if not (r.imag < 0 and abs(r - (830 - 40j)) < 400):
                    raise ValueError(f"zero at {r:.1f} MeV is outside the accepted window")
                res["pole"][presc] = {"sqrt_s_MeV": [float(r.real), float(r.imag)]}
            except (RuntimeError, ZeroDivisionError, ValueError) as exc:
                res["pole"][presc] = None
                res.setdefault("pole_errors", {})[presc] = str(exc)
    return res
