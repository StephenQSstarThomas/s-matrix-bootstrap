"""V2 retro for 2309 leaves: rigorous upper bound from SDPB's multipliers on a full-space degree-0 PMP.

Usage:
  python scripts/gtb/v2_2309.py p1 LEAF OUT_JSON [--paper X]          full-space leaf, its own multipliers (P0+P1)
  python scripts/gtb/v2_2309.py p2 FULL_LEAF Z_LEAF OUT_JSON          Z of a reduced leaf on FULL_LEAF's blocks
  python scripts/gtb/v2_2309.py p3 OUT_JSON [LEAF ...]                UV exclusion leaves (table P3 below)
  [--refine K]  Gram-repair iterations of dual_bound (default 8 for p1/p2; p3 uses 0, then 1 if needed)

Weak duality (gtb/dual.py:dual_bound): for every feasible y,
  b.y <= b0 + sum_j <Z_j, C_j> + sum_n |r_n| ybar_n + sum_j max(0, -lambda_min Z_j) tau_j.
The certified problem is an emitted 30-digit PMP, parsed exactly into Arb balls: pmp.json of a full-space leaf
(p1/p2), or for p3 the full-space PMP regenerated in memory by the unchanged generator sdp/pmp.py:Pmp._blocks with
the leaf's spec and reduce_basis=False, checked block by block against the leaf's own reduced pmp.json (non-amplitude
coefficients equal as strings, amplitude rows equal after projection through basis.npy).  Z_j from out/x_j.txt
(column-major upper triangle, off-diagonals doubled; checked in P0).
UV relaxation (p3): rho_hat at nodes above s0 appears only in its Gram (3,3) entry and rho_hat >= 0, so it has no
a-priori bound.  Those Gram blocks are replaced by their leading 2x2 (|S| <= 1) and the variables and their 1x1
blocks are dropped.  The relaxed feasible set contains the stated one, so its upper bound is valid for the stated
problem; exclusions 'max X < c' need this upper bound, 'min X > c' is certified as max(-X) < -c.
A-priori bounds (B8 memo section 3), all Arb upper ends:
  rho:   single-variable 1x1 blocks (1 -/+ rho/Mreg >= 0) -> |rho_i| <= Mreg; rho_hat >= 0 likewise one-sided;
  sigma: S-wave disks [[A,P,B],[P,1,0],[B,0,1]]: A.y >= (B.y)^2 + (P.y)^2, B ~ mu A  ->  0 <= A.y <= Xbar;
         S0 Gram blocks: (0,0) = 2 y0 - D1 exactly (string check), so 0 <= Im h = D1.y <= 2;
  T0:    |P.y| <= sqrt(Xbar) on disks, |Re h| <= 1 on Gram blocks;
  rho_hat (nodes <= s0): FESR 1x1 boxes  t + d - sum_i w_i rho_hat_i >= 0  with w_i > 0  ->  rho_hat_i <= (t+d)/w_i;
  ImF:   Gram rows (1,2): beta^2 ImF^2 <= (Im h)(gamma rho_hat) <= 2 gamma rho_hat_bar;  FF blocks: |beta ImF| <= b.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from flint import arb, ctx

from smatrix_bootstrap.gtb.dual import dual_bound
from smatrix_bootstrap.sdp.dual import _stream_pmp

ZERO = arb(0)
PREC = 256
ZSTR = ("0", "0.0", "-0", "-0.0")
R = "/playpen1/shiqiu/sdpb_regularised_20260913/"
T_S0N0 = arb("3.09e-8") * (arb("27.38") / 2 + arb(".61")) * (arb(3600) / 49) ** 2    # precision.py:91-94, n=0
P3 = [  # (leaf, statement, threshold c of 'max objective < c')
    ("uv_SRa_tip/tip", "UV +x end f00(3) < paper 0.0811249 (C5 FAIL is not a solver artifact)", arb("0.0811249")),
    ("srmom_SRb_free_S0n0_min/support", "SR-1: min S0 n=0 moment > 1.1 T, i.e. the SR-b +/-10% box is empty "
     "(objective is -moment, so max(-m) < -1.1 T)", -arb("1.1") * T_S0N0),
    ("face_SRa_tip_ImKH_P1_38_max", "F1 (elastic BW): max on tip face of 1-Re S_P1(0.792 GeV) < 1.5", arb("1.5")),
    ("face_SRa_tip_ImKH_P1_39_max", "F3': max on tip face of 1-Re S_P1(0.864 GeV) < 1", arb(1)),
    ("face_SRa_tip_ImS_P1_38_max", "F6: max on tip face of Im S_P1(0.792 GeV) < 0", arb(0)),
    ("face_SRa_xrefhi_ImKH_P1_38_max", "F1': max on x_ref upper face of 1-Re S_P1(0.792 GeV) < 1", arb(1)),
    ("face_SRa_xrefhi_ImKH_P1_39_max", "F3': max on x_ref upper face of 1-Re S_P1(0.864 GeV) < 1", arb(1)),
    ("face_SRa_ffs0_tip_ImKH_P1_38_max", "F1': max on frozen-factor tip face of 1-Re S_P1(0.792 GeV) < 1", arb(1)),
    ("face_SRa_ffs0_tip_ImKH_P1_39_max", "F3': max on frozen-factor tip face of 1-Re S_P1(0.864 GeV) < 1", arb(1)),
]


class Block:
    """Same interface as gtb.program.Block: forms[(r, c)] (r <= c) is the coefficient vector over [y0, y1, ...]."""

    def __init__(self, size, forms):
        self.size, self.forms = size, forms

    def entry(self, r, c):
        return self.forms.get((min(r, c), max(r, c)))


def _neg(s):
    return s[1:] if s.startswith("-") else "-" + s


def ingest(m, cells, n):
    """cells[(r, c)], r <= c: n decimal strings each.  Returns Block, sparse {(r,c): (idx, vals)} and a kind tag."""
    forms, sp = {}, {}
    for key, strs in cells.items():
        if len(strs) != n:
            raise ValueError("entry length differs from the number of variables")
        idx, vals = [], []
        for k, s in enumerate(strs):
            if s in ZSTR:
                continue
            v = arb(s)
            if v != 0:
                idx.append(k)
                vals.append(v)
        if idx:
            va = np.empty(len(vals), dtype=object)
            va[:] = vals
            f = np.full(n, ZERO, dtype=object)
            f[idx] = va
            forms[key] = f
            sp[key] = (np.array(idx, dtype=int), va)
    return Block(m, forms), sp, classify(m, sp, cells)


def _single(sp, key):
    """(var, coef, const) if entry key is const + coef * y_var with one variable (const may be zero)."""
    if key not in sp:
        return None
    idx, vals = sp[key]
    var = [(k, v) for k, v in zip(idx, vals) if k != 0]
    const = vals[0] if idx[0] == 0 else ZERO
    return (int(var[0][0]), var[0][1], const) if len(var) == 1 else None


def classify(m, sp, cells):
    if m == 1:
        return "1x1"
    if m != 3:
        return "other"
    y0one = all(k in sp and list(sp[k][0]) == [0] and sp[k][1][0] == 1 for k in ((1, 1), (2, 2)))
    if y0one and (1, 2) not in sp and set(sp) == {(0, 0), (0, 1), (0, 2), (1, 1), (2, 2)}:
        return "disk"
    d0, d1 = cells[(0, 0)], cells[(1, 1)]
    if (d1[0] in ZSTR and arb(d0[0]) == 2 and len(sp.get((1, 1), ((), ()))[0]) > 1
            and all((a in ZSTR and b in ZSTR) or a == _neg(b) for a, b in zip(d0[1:], d1[1:]))):
        g = _single(sp, (2, 2))
        if g is not None and g[2] == 0 and g[1] == 1:
            return "gram"                            # [[2-Im h, Re h, .], [Re h, Im h, .], [., ., rho_hat]]
    diag = [sp.get(k) for k in ((0, 0), (1, 1), (2, 2))]
    if (all(d is not None and list(d[0]) == [0] for d in diag) and cells[(0, 0)] == cells[(1, 1)] == cells[(2, 2)]
            and (1, 2) not in sp and _single(sp, (0, 2)) is not None):
        return "ff"                                  # [[b, cF_re, cF_im], [cF_re, b, 0], [cF_im, 0, b]]
    return "other"


def load_pmp(path):
    """Parse pmp.json into (objective, blocks, sparse, kinds)."""
    ctx.prec = PREC
    stream = _stream_pmp(path)
    head = next(stream)
    n = len(head["objective"])
    norm = [arb(v) for v in head["normalization"]]
    if not (norm[0] == 1 and all(v == 0 for v in norm[1:])):
        raise ValueError("expected the y0 = 1 normalisation")
    blocks, sparse, kinds = [], [], []
    for blk in stream:
        P = blk["polynomials"]
        if any(len(p) != 1 for row in P for cell in row for p in cell):
            raise ValueError("only degree-0 PMPs are supported")
        cells = {(r, c): [p[0] for p in P[r][c]] for c in range(len(P)) for r in range(c + 1)}
        b, sp, kind = ingest(len(P), cells, n)
        blocks.append(b), sparse.append(sp), kinds.append(kind)
    return [arb(v) for v in head["objective"]], blocks, sparse, kinds


def full_space_from_reduced(leaf):
    """Regenerate the leaf's full-space PMP in memory (unchanged generator, reduce_basis=False) and check it block by
    block against the leaf's reduced pmp.json.  Returns (objective, blocks, sparse, kinds, identity report)."""
    from smatrix_bootstrap.sdp.pmp import Pmp
    from smatrix_bootstrap.sdp.spec import ModelSpec
    leaf = Path(leaf)
    rep = json.loads((leaf / "report.json").read_text())
    s = dict(rep["spec"], reduce_basis=False)
    s["uv_parts"] = tuple(s["uv_parts"])
    s["sr_free"] = tuple(tuple(x) for x in s.get("sr_free") or ())
    w = Pmp(ModelSpec(**s), direction=tuple(rep["direction"]), fix_f00=rep.get("fix_f00"),
            digits=rep["pmp"]["digits"], face=rep.get("face"), functional=rep.get("functional"))
    ctx.threads = 4
    bits, n, na = w.precise.bits, w.n_vars, w.n_a
    V = np.load(leaf / "basis.npy")
    red = _stream_pmp(leaf / "pmp.json")
    rhead = next(red)
    nr = len(rhead["objective"])
    worst = [0.0]

    def same(fs, rs):
        if fs[0] != rs[0] or fs[1 + na:] != rs[nr - (n - 1 - na):]:
            raise ValueError("non-amplitude coefficients differ from the reduced leaf")
        fa = np.array([float(x) for x in fs[1:1 + na]])
        ra = np.array([float(x) for x in rs[1:nr - (n - 1 - na)]])
        if fa.any() or ra.any():
            dev = np.abs(fa @ V - ra).max() / max((np.abs(fa) @ np.abs(V)).max(), 1e-300)
            worst[0] = max(worst[0], float(dev))
            if dev > 1e-9:
                raise ValueError(f"amplitude row differs from the reduced leaf after projection ({dev:.2e})")
    ctx.prec = bits
    obj = w._num_vec(w.objective())
    same(obj, rhead["objective"])
    blocks, sparse, kinds = [], [], []
    for blk in w._blocks():
        m = len(blk)
        cells = {(r, c): w._num_vec(blk[r][c]) for c in range(m) for r in range(c + 1)}
        rb = next(red)["polynomials"]
        if len(rb) != m:
            raise ValueError(f"block {len(blocks)}: size {m} vs reduced {len(rb)}")
        for (r, c), fs in cells.items():
            same(fs, [p[0] for p in rb[r][c]])
        ctx.prec = PREC
        b, sp, kind = ingest(m, cells, n)
        blocks.append(b), sparse.append(sp), kinds.append(kind)
        ctx.prec = bits
    if next(red, None) is not None:
        raise ValueError("reduced leaf has more blocks")
    ctx.prec = PREC
    return [arb(v) for v in obj], blocks, sparse, kinds, {
        "n_vars_full": n, "n_vars_reduced": nr, "n_blocks": len(blocks), "nonamplitude_strings_equal": True,
        "max_rel_amplitude_projection_deviation": worst[0], "generator": "sdp/pmp.py Pmp._blocks, reduce_basis=False"}


def relax_rho_hat(objective, blocks, sparse, kinds, Z):
    """Drop variables seen only in one Gram (2,2) entry and single-variable 1x1 blocks; shrink that Gram to 2x2."""
    seen = defaultdict(list)
    for j, sp in enumerate(sparse):
        for key, (idx, _) in sp.items():
            for k in idx:
                seen[int(k)].append((j, key))
    drop, shrink, kill = [], [], set()
    for v, locs in seen.items():
        g = [j for j, key in locs if kinds[j] == "gram" and key == (2, 2)]
        ones = [j for j, key in locs if kinds[j] == "1x1" and _single(sparse[j], (0, 0)) is not None]
        if v and len(g) == 1 and len(locs) == 1 + len(ones) and objective[v] == 0:
            drop.append(v), shrink.append(g[0]), kill.update(ones)
    for j in shrink:
        keep = ((0, 0), (0, 1), (1, 1))
        blocks[j] = Block(2, {k: f for k, f in blocks[j].forms.items() if k in keep})
        sparse[j] = {k: e for k, e in sparse[j].items() if k in keep}
        kinds[j], Z[j] = "gram2", [row[:2] for row in Z[j][:2]]
    alive = [j for j in range(len(blocks)) if j not in kill]
    blocks, sparse, kinds, Z = ([x[j] for j in alive] for x in (blocks, sparse, kinds, Z))
    n = len(objective)
    keepv = np.array([k for k in range(n) if k not in set(drop)])
    remap = -np.ones(n, dtype=int)
    remap[keepv] = np.arange(len(keepv))
    for b, sp in zip(blocks, sparse):
        for key in list(sp):
            idx, vals = sp[key]
            if (remap[idx] < 0).any():
                raise ValueError("a dropped variable is still referenced")
            sp[key] = (remap[idx], vals)
            b.forms[key] = b.forms[key][keepv]
    return [objective[k] for k in keepv], blocks, sparse, kinds, Z, sorted(drop)


def read_x(outdir, sizes, mapping="colmajor_halved"):
    """Z_j from SDPB's out/x_j.txt (degree 0: one entry per upper-triangle element)."""
    Z = []
    for j, m in enumerate(sizes):
        tok = Path(outdir, f"x_{j}.txt").read_text().split()
        cnt = int(tok[0])
        if cnt != m * (m + 1) // 2 or int(tok[1]) != 1 or len(tok) != 2 + cnt:
            raise ValueError(f"x_{j}.txt: {tok[:2]} does not match a {m}x{m} degree-0 block")
        order = ([(r, c) for c in range(m) for r in range(c + 1)] if mapping.startswith("colmajor")
                 else [(r, c) for r in range(m) for c in range(r, m)])
        Zj = [[ZERO] * m for _ in range(m)]
        for (r, c), s in zip(order, tok[2:]):
            x = arb(s)
            Zj[r][c] = Zj[c][r] = x / 2 if (r != c and mapping.endswith("_halved")) else x
        Z.append(Zj)
    if Path(outdir, f"x_{len(sizes)}.txt").exists():
        raise ValueError("more x_j files than PMP blocks")
    return Z


def stationarity(objective, sparse, Z):
    """r_n = b_n + sum_j <Z_j, A_jn> in Arb; r_0 = b0 + sum_j <Z_j, C_j> is the dual objective."""
    r = np.array(objective, dtype=object)
    for sp, z in zip(sparse, Z):
        for (rr, cc), (idx, vals) in sp.items():
            r[idx] = r[idx] + (z[rr][cc] if rr == cc else 2 * z[rr][cc]) * vals
    return r


def _up(x):
    return arb(x.upper())


def _absup(x):
    return arb(abs(x).upper())


def a_priori_bounds(sparse, kinds, n, M):
    """ybar[n-1] for every free variable (amplitude layout of sdp/projector.py first: T0, sigma1, sigma2, rho)."""
    i_T0, s1, s2, rho0, amp = 1, 2, 2 + M, 2 + 2 * M, 2 + 2 * M + M * M + M * (M + 1) // 2
    lo, hi = [None] * n, [None] * n
    pick = {True: lambda a, b: b if a is None else max(a, b, key=lambda t: float(t.lower())),
            False: lambda a, b: b if a is None else min(a, b, key=lambda t: float(t.upper()))}
    for sp, kind in zip(sparse, kinds):              # (B0) single-variable 1x1 blocks: a0 + a_k y_k >= 0
        one = _single(sp, (0, 0)) if kind == "1x1" else None
        if one:
            k, a, a0 = one
            if a > 0:
                lo[k] = pick[True](lo[k], -a0 / a)
            elif a < 0:
                hi[k] = pick[False](hi[k], a0 / (-a))
    ybar = [None] * n
    for k in range(rho0, amp):
        if lo[k] is None or hi[k] is None:
            raise ValueError(f"no two-sided 1x1 bound on rho variable {k}")
        ybar[k] = max(_absup(lo[k]), _absup(hi[k]), key=float)

    def rest(d, skip, extra=lambda k: ybar[k]):
        return sum((_absup(v) * (arb(1) if k == 0 else extra(k)) for k, v in d.items() if k not in skip), ZERO)

    def xbar(j, key, extra):
        """Disk j: 0 <= A.y <= Xbar from A.y >= (B.y)^2, B.y = mu* A.y + eps, mu* = B_key/A_key exact, |eps| <= eta."""
        Ad, Bd = dict(zip(*sparse[j][(0, 0)])), dict(zip(*sparse[j][(0, 2)]))
        if set(Ad) != set(Bd) or 0 in Ad or i_T0 in Ad:
            raise ValueError("Im rows must share support and carry no constant or T0")
        mu = Bd[key] / Ad[key]
        if not mu > 0:
            raise ValueError("non-positive Im-row ratio")
        eta = sum((_absup(Bd[k] - mu * Ad[k]) * extra(k) for k in Ad if k != key), ZERO)
        mlo = arb(mu.lower())
        return max(_up((1 + 2 * mlo * eta) / mlo ** 2), _up(eta / mlo), key=float), Ad

    sig = lambda j, key: {int(v) for v in sparse[j][key][0] if s1 <= v < rho0}
    i2 = {next(iter(sig(j, (0, 0)))) - s2: j for j, kd in enumerate(kinds) if kd == "disk" and len(sig(j, (0, 0))) == 1}
    i0 = {min(sig(j, (0, 0))) - s1: j for j, kd in enumerate(kinds) if kd == "disk" and len(sig(j, (0, 0))) == 2}
    g0 = {min(sig(j, (1, 1))) - s1: j for j, kd in enumerate(kinds) if kd in ("gram", "gram2") and len(sig(j, (1, 1))) == 2}
    cands = []                                       # (Xbar of |Re part|^2, j, Re-entry key)
    for k in range(M):                               # (B2) I=2, l=0 disk: A = c sigma2_k + rho terms
        j = i2[k]
        if sig(j, (0, 0)) != {s2 + k}:
            raise ValueError(f"block {j} does not isolate sigma2_{k}")
        Xb, Ad = xbar(j, s2 + k, lambda v: ybar[v] if v >= rho0 else None)
        ybar[s2 + k] = _up((Xb + rest(Ad, {s2 + k})) / abs(Ad[s2 + k]))
        cands.append((Xb, j, (0, 1)))
    for k in range(M):                               # (B2) I=0, l=0: disk (Xbar) or S0 Gram (0 <= Im h <= 2)
        if k in i0:
            j = i0[k]
            Xb, Ad = xbar(j, s1 + k, lambda v: ybar[v])
            cands.append((Xb, j, (0, 1)))
        else:
            j, Xb = g0[k], arb(2)
            Ad = dict(zip(*sparse[j][(1, 1)]))
            cands.append((arb(1), j, (0, 1)))        # (Re h)^2 <= Im h (2 - Im h) <= 1
        if sig(j, (0, 0) if k in i0 else (1, 1)) != {s1 + k, s2 + k} or 0 in Ad or i_T0 in Ad:
            raise ValueError(f"block {j} does not isolate sigma1_{k}, sigma2_{k}")
        ybar[s1 + k] = _up((Xb + rest(Ad, {s1 + k})) / abs(Ad[s1 + k]))
    best = None                                      # (B3) T0 from |Re part| <= sqrt(Xbar)
    for Xb, j, key in cands:
        Pd = dict(zip(*sparse[j][key]))
        if i_T0 in Pd:
            t = _up((Xb.sqrt() + rest(Pd, {i_T0})) / abs(Pd[i_T0]))
            best = (t, j) if best is None or float(t) < float(best[0]) else best
    ybar[i_T0] = best[0]
    for sp, kind in zip(sparse, kinds):              # (B5) t + d - sum w_i y_i >= 0, w_i > 0, y_i >= 0
        if kind != "1x1" or _single(sp, (0, 0)) is not None:
            continue
        idx, vals = sp[(0, 0)]
        a0 = vals[0] if idx[0] == 0 else ZERO
        var = [(int(k), v) for k, v in zip(idx, vals) if k != 0]
        if all(v < 0 and k >= amp and lo[k] is not None and lo[k] >= 0 for k, v in var):
            for k, v in var:
                hi[k] = pick[False](hi[k], a0 / (-v))
    for k in range(amp, n):
        if lo[k] is not None and hi[k] is not None and lo[k] >= 0:
            ybar[k] = _absup(hi[k])
    rho_hat = [k for k in range(amp, n) if ybar[k] is not None]
    for sp, kind in zip(sparse, kinds):              # (B6) Im F from Gram rows (1,2) and from the FF caps
        if kind == "gram":
            f, g = _single(sp, (1, 2)), _single(sp, (2, 2))
            if f and g and f[2] == 0 and ybar[g[0]] is not None:
                b = _up((2 * _absup(g[1]) * ybar[g[0]]).sqrt() / abs(f[1]))
                ybar[f[0]] = b if ybar[f[0]] is None else min(ybar[f[0]], b, key=float)
        elif kind == "ff":
            f, b0 = _single(sp, (0, 2)), sp[(0, 0)][1][0]
            if f[2] == 0:
                b = _up(abs(b0) / abs(f[1]))
                ybar[f[0]] = b if ybar[f[0]] is None else min(ybar[f[0]], b, key=float)
    if any(v is None for v in ybar[1:]):
        raise ValueError(f"unbounded free variables: {[k for k in range(1, n) if ybar[k] is None][:20]}")
    rng = lambda ks: [min(float(ybar[k]) for k in ks), max(float(ybar[k]) for k in ks)] if ks else None
    imf = [k for k in range(amp, n) if k not in rho_hat]
    return ybar[1:], {"T0_block": best[1], "T0": float(best[0]), "sigma1_range": rng(range(s1, s2)),
                      "sigma2_range": rng(range(s2, rho0)), "rho": rng(range(rho0, amp)),
                      "rho_hat_range": rng(rho_hat), "ImF_range": rng(imf), "n_rho_hat": len(rho_hat), "n_ImF": len(imf)}


def trace_bounds(sparse, ybar):
    """tau_j >= tr M_j(y): sum over diagonal entries of |const| + sum_n |coef_n| ybar_n (Arb upper ends)."""
    return [_up(sum((_absup(v) * (arb(1) if k == 0 else ybar[k - 1]) for (r, c), (idx, vals) in sp.items() if r == c
                     for k, v in zip(idx, vals)), ZERO)) for sp in sparse]


def upper_str(x, digits=25):
    """A decimal string >= the upper end of x."""
    u = arb(x.upper())
    s = (u + abs(u) * arb(10) ** (3 - digits) + arb(10) ** (-60)).upper().str(digits, radius=False)
    if not arb(s) >= u:
        raise ArithmeticError("upper_str rounding")
    return s


def sdpb_out(leaf):
    txt = Path(leaf, "out", "out.txt").read_text()
    return {k.strip(): v.strip().rstrip(";").strip('"') for k, v in (ln.split("=", 1) for ln in txt.splitlines() if "=" in ln)}


def p0(objective, sparse, sizes, outdir, solver):
    """Mapping check: Arb stationarity residual and dual objective vs SDPB's primalObjective."""
    res = {}
    for mapping in ("colmajor_halved", "colmajor_unhalved", "rowmajor_halved"):
        r = stationarity(objective, sparse, read_x(outdir, sizes, mapping))
        res[mapping] = {"dual_objective": r[0].mid().str(40, radius=False),
                        "max_abs_residual": max(float(abs(x).upper()) for x in r[1:]),
                        "minus_primalObjective": float((r[0] - arb(solver["primalObjective"])).mid())}
    ok = res["colmajor_halved"]
    passed = ok["max_abs_residual"] < 1e-30 and abs(ok["minus_primalObjective"]) < 1e-30
    return {"mappings": res, "accepted_mapping": "colmajor_halved" if passed else None, "passed": passed,
            "sdpb": {k: solver.get(k) for k in ("primalObjective", "dualObjective", "dualityGap", "primalError", "dualError")}}


def bound(blocks, sparse, kinds, objective, Z, M, refine):
    ybar, binfo = a_priori_bounds(sparse, kinds, len(objective), M)
    tau = trace_bounds(sparse, ybar)
    t1 = time.time()
    res = dual_bound(blocks, objective, Z, tau, ybar, refine=refine)      # refine=0: SDPB's Z, no Gram repair
    ub = res["bound"]
    return ub, {"a_priori_bounds": binfo, "tau_max": float(max(tau, key=float)), "gram_repair_iterations": refine,
                "certified": bool(ub.is_finite()), "bound_upper": upper_str(ub) if ub.is_finite() else "inf",
                "bound_mid": float(ub.mid()) if ub.is_finite() else None, "residual_history": res.get("residual_history"),
                "max_psd_shift": res.get("max_shift"), "failed_blocks": res.get("failed_blocks"),
                "dual_bound_seconds": time.time() - t1}


def certify(full_leaf, z_leaf, out_json, paper=None, refine=8):
    """P1 (z_leaf == full_leaf, with P0) and P2 (reduced z_leaf on a full-space leaf with the same block list)."""
    t0 = time.time()
    full_leaf, z_leaf = Path(full_leaf), Path(z_leaf)
    spec = json.loads((full_leaf / "report.json").read_text())["spec"]
    if spec.get("uv") or spec.get("reduce_basis"):
        raise ValueError("p1/p2 need a full-space chiral-only leaf for the blocks; use p3 for UV leaves")
    objective, blocks, sparse, kinds = load_pmp(full_leaf / "pmp.json")
    sizes = [b.size for b in blocks]
    print(f"parsed {len(blocks)} blocks, n={len(objective)} in {time.time() - t0:.0f}s", flush=True)
    out = {"leaf_blocks": str(full_leaf), "leaf_Z": str(z_leaf), "scope": "emitted full-space degree-0 PMP "
           "(30-digit decimals parsed exactly into Arb balls); chiral-only leaf, so no rho_hat variables and no "
           "relaxation is used", "n_vars": len(objective), "n_blocks": len(blocks)}
    solver_z = sdpb_out(z_leaf)
    if z_leaf == full_leaf:
        out["P0"] = p0(objective, sparse, sizes, z_leaf / "out", solver_z)
        if not out["P0"]["passed"]:
            Path(out_json).write_text(json.dumps(out, indent=1))
            raise SystemExit("P0 failed: x_j -> Z mapping not confirmed")
    Z = read_x(z_leaf / "out", sizes)
    if z_leaf != full_leaf:                          # P2: raw full-space stationarity of a foreign Z
        r = stationarity(objective, sparse, Z)
        out["raw_foreign_Z"] = {"dual_objective": float(r[0].mid()), "Z_leaf_primalObjective": solver_z["primalObjective"],
                                "max_abs_residual": max(float(abs(x).upper()) for x in r[1:]),
                                "sum_abs_residual": sum(float(abs(x).upper()) for x in r[1:])}
    ub, info = bound(blocks, sparse, kinds, objective, Z, spec["M"], refine)
    out.update(info, wall_seconds=time.time() - t0)
    print(f"a-priori bounds {info['a_priori_bounds']}", flush=True)
    s_full = sdpb_out(full_leaf)
    if ub.is_finite():
        bu = float(arb(out["bound_upper"]))
        out["comparison"] = {"sdpb_primalObjective": s_full["primalObjective"], "sdpb_dualObjective": s_full["dualObjective"],
                             "bound_minus_primalObjective": (arb(out["bound_upper"]) - arb(s_full["primalObjective"])).mid().str(5, radius=False),
                             "bound_minus_dualObjective": bu - float(s_full["dualObjective"])}
        if paper is not None:
            out["comparison"].update(paper=paper, bound_minus_paper=bu - paper, rel_to_paper=(bu - paper) / paper)
    Path(out_json).write_text(json.dumps(out, indent=1))
    return out


def certify_uv(leaf, statement, c):
    """P3: full-space relaxed UV problem, the reduced leaf's own multipliers, verdict on 'max objective < c'."""
    t0 = time.time()
    leaf = Path(R + leaf)
    objective, blocks, sparse, kinds, ident = full_space_from_reduced(leaf)
    print(f"{leaf.name}: full-space PMP regenerated and matched in {time.time() - t0:.0f}s", flush=True)
    Z = read_x(leaf / "out", [b.size for b in blocks])
    objective, blocks, sparse, kinds, Z, drop = relax_rho_hat(objective, blocks, sparse, kinds, Z)
    M = json.loads((leaf / "report.json").read_text())["spec"]["M"]
    out = {"leaf": str(leaf), "statement": statement, "threshold": c.mid().str(20, radius=False), "identity": ident,
           "relaxation": {"dropped_variables": len(drop), "dropped_index_range": [drop[0], drop[-1]] if drop else None,
                          "scope": "Gram blocks of the rho_hat variables above s0 replaced by their leading 2x2 (|S|<=1); "
                                   "this relaxed set contains the stated feasible set, so the bound is an upper bound for "
                                   "the stated full-space problem as well"},
           "sdpb_reduced_leaf": {k: sdpb_out(leaf).get(k) for k in ("primalObjective", "dualObjective")}, "runs": []}
    best = None
    for refine in (0, 1):
        ub, info = bound(blocks, sparse, kinds, objective, Z, M, refine)
        out["runs"].append(info)
        print(f"  refine={refine}: bound {info['bound_upper']} shift {info['max_psd_shift']}", flush=True)
        if ub.is_finite() and (best is None or float(ub.upper()) < float(best.upper())):
            best = ub
        if best is not None and best < c:
            break
    out.update({"certified_bound_upper": upper_str(best) if best is not None else "inf",
                "margin": float((c - arb(best.upper())).mid()) if best is not None else None,
                "verdict": "proven" if best is not None and best < c else "numerical-only",
                "seconds": time.time() - t0})
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] not in ("p1", "p2", "p3"):
        raise SystemExit(__doc__)
    paper = float(a[a.index("--paper") + 1]) if "--paper" in a else None
    refine = int(a[a.index("--refine") + 1]) if "--refine" in a else 8
    if a[0] == "p3":
        names = [x for x in a[2:] if not x.startswith("--")] or [p[0] for p in P3]
        results, path = [], Path(a[1])
        for leaf, statement, c in P3:
            if leaf in names:
                results.append(certify_uv(leaf, statement, c))
                path.write_text(json.dumps({"recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                            "leaves": results}, indent=1))
        print(json.dumps([{k: r[k] for k in ("leaf", "certified_bound_upper", "margin", "verdict")} for r in results], indent=1))
    else:
        r = certify(a[1], a[1] if a[0] == "p1" else a[2], a[2] if a[0] == "p1" else a[3], paper, refine)
        print(json.dumps({k: r.get(k) for k in ("bound_upper", "certified", "max_psd_shift", "wall_seconds", "comparison")}, indent=1))
