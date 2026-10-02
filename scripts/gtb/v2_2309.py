"""V2 retro for 2309 leaves: rigorous upper bound from SDPB's multipliers on the emitted degree-0 PMP.

Usage:
  python scripts/gtb/v2_2309.py p1 LEAF OUT_JSON [--paper 0.0825728]
  python scripts/gtb/v2_2309.py p2 FULL_LEAF Z_LEAF OUT_JSON       (Z of another leaf on FULL_LEAF's blocks)
  [--refine K]  Gram-repair iterations of dual_bound (default 8; 0 keeps SDPB's Z and pays sum |r_n| ybar_n)

Weak duality (gtb/dual.py:dual_bound): for every feasible y,
  b.y <= b0 + sum_j <Z_j, C_j> + sum_n |r_n| ybar_n + sum_j max(0, -lambda_min Z_j) tau_j.
Blocks are parsed from pmp.json (30-digit decimal strings -> Arb balls containing the exact decimals), so the
certificate is for the emitted PMP.  Z_j comes from out/x_j.txt (no X_matrix files exist for 2309 leaves); the
mapping is checked first (P0) against SDPB's primal objective and the dual equalities.
A-priori bounds (B8 memo section 3), all from the PMP rows in Arb:
  rho:   single-variable 1x1 blocks (regulariser 1 -/+ rho/Mreg >= 0)  ->  |rho_i| <= Mreg;
  sigma: S-wave disks [[A,P,B],[P,1,0],[B,0,1]] >= 0 give A.y >= (B.y)^2 + (P.y)^2 with B ~ mu A, hence
         0 <= A.y <= Xbar; A carries only sigma_{1,k}, sigma_{2,k} and rho (nodal prescription);
  T0:    |P.y| <= sqrt(Xbar) on the same blocks.
Full-space chiral-only leaves only.  UV leaves (rho_hat above s0 unbounded; relaxation to the 2x2 |S|<=1 minor)
are refused here.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, ctx

from smatrix_bootstrap.gtb.dual import dual_bound
from smatrix_bootstrap.sdp.dual import _stream_pmp

ZERO = arb(0)
PREC = 256


class Block:
    """Same interface as gtb.program.Block: forms[(r, c)] (r <= c) is the coefficient vector over [y0, y1, ...]."""

    def __init__(self, size, forms):
        self.size, self.forms = size, forms

    def entry(self, r, c):
        return self.forms.get((min(r, c), max(r, c)))


def load_pmp(path):
    """Parse pmp.json; returns objective (arb list), Blocks (dense object forms) and sparse forms {(r,c): (idx, vals)}."""
    ctx.prec = PREC
    stream = _stream_pmp(path)
    head = next(stream)
    n = len(head["objective"])
    norm = [arb(v) for v in head["normalization"]]
    if not (norm[0] == 1 and all(v == 0 for v in norm[1:])):
        raise ValueError("expected the y0 = 1 normalisation")
    objective = [arb(v) for v in head["objective"]]
    blocks, sparse = [], []
    for blk in stream:
        P = blk["polynomials"]
        m = len(P)
        forms, sp = {}, {}
        for c in range(m):
            for r in range(c + 1):
                cell = P[r][c]
                if len(cell) != n or any(len(p) != 1 for p in cell):
                    raise ValueError("only degree-0 PMPs with n coefficients per entry are supported")
                idx, vals = [], []
                for k, p in enumerate(cell):
                    if p[0] in ("0", "0.0", "-0", "-0.0"):
                        continue
                    v = arb(p[0])
                    if v != 0:
                        idx.append(k)
                        vals.append(v)
                if not idx:
                    continue
                va = np.empty(len(vals), dtype=object)
                va[:] = vals
                f = np.full(n, ZERO, dtype=object)
                f[idx] = va
                forms[(r, c)] = f
                sp[(r, c)] = (np.array(idx, dtype=int), va)
        blocks.append(Block(m, forms))
        sparse.append(sp)
    return objective, blocks, sparse


def read_x(outdir, sizes, mapping="colmajor_halved"):
    """Z_j from SDPB's out/x_j.txt (degree 0: one entry per upper-triangle element)."""
    Z = []
    for j, m in enumerate(sizes):
        tok = Path(outdir, f"x_{j}.txt").read_text().split()
        cnt = int(tok[0])
        if cnt != m * (m + 1) // 2 or int(tok[1]) != 1 or len(tok) != 2 + cnt:
            raise ValueError(f"x_{j}.txt: {tok[:2]} does not match a {m}x{m} degree-0 block")
        vals = tok[2:]
        order = ([(r, c) for c in range(m) for r in range(c + 1)] if mapping.startswith("colmajor")
                 else [(r, c) for r in range(m) for c in range(r, m)])
        Zj = [[ZERO] * m for _ in range(m)]
        for (r, c), s in zip(order, vals):
            x = arb(s)
            v = x / 2 if (r != c and mapping.endswith("_halved")) else x
            Zj[r][c] = Zj[c][r] = v
        Z.append(Zj)
    if Path(outdir, f"x_{len(sizes)}.txt").exists():
        raise ValueError("more x_j files than PMP blocks")
    return Z


def stationarity(objective, sparse, Z):
    """r_n = b_n + sum_j <Z_j, A_jn> in Arb; r_0 = b0 + sum_j <Z_j, C_j> is the dual objective."""
    r = np.array(objective, dtype=object)
    for sp, z in zip(sparse, Z):
        for (rr, cc), (idx, vals) in sp.items():
            w = z[rr][cc] if rr == cc else 2 * z[rr][cc]
            r[idx] = r[idx] + w * vals
    return r


def _up(x):
    return arb(x.upper())


def _absup(x):
    return arb(abs(x).upper())


def a_priori_bounds(sparse, n, M, L):
    """ybar[n-1] for every free variable n >= 1 of the full-space chiral-only problem (amplitude layout of
    sdp/projector.py: T0, sigma1[M], sigma2[M], rho1[M*M], rho2[M(M+1)/2])."""
    i_T0, s1, s2, rho0 = 1, 2, 2 + M, 2 + 2 * M
    if n != rho0 + M * M + M * (M + 1) // 2:
        raise ValueError(f"n_vars={n} is not the full-space chiral-only layout for M={M}")
    lo, hi = [None] * n, [None] * n
    for sp in sparse:                                  # (B0) single-variable 1x1 blocks: a0 + a_k y_k >= 0
        if set(sp) != {(0, 0)}:
            continue
        idx, vals = sp[(0, 0)]
        a0 = vals[0] if idx[0] == 0 else ZERO
        var = [(k, v) for k, v in zip(idx, vals) if k != 0]
        if len(var) != 1:
            continue
        k, a = var[0]
        if a > 0:
            b = -a0 / a
            lo[k] = b if lo[k] is None else max(lo[k], b, key=lambda t: float(t.lower()))
        elif a < 0:
            b = a0 / (-a)
            hi[k] = b if hi[k] is None else min(hi[k], b, key=lambda t: float(t.upper()))
    ybar = [None] * n
    for k in range(rho0, n):
        if lo[k] is None or hi[k] is None:
            raise ValueError(f"no two-sided 1x1 bound on rho variable {k}")
        ybar[k] = max(_absup(lo[k]), _absup(hi[k]), key=float)
    rho_bar = ybar[rho0:]

    def disk(j):                                       # [[A, P, B], [P, 1, 0], [B, 0, 1]] with exact 1 and 0
        sp = sparse[j]
        if set(sp) != {(0, 0), (0, 1), (0, 2), (1, 1), (2, 2)}:
            raise ValueError(f"block {j} is not a unitarity disk")
        for d in ((1, 1), (2, 2)):
            idx, vals = sp[d]
            if list(idx) != [0] or not vals[0] == 1:     # exact 1 (Arb == needs exact balls)
                raise ValueError(f"block {j}: diagonal {d} is not y0")
        return sp[(0, 0)], sp[(0, 1)], sp[(0, 2)]

    def xbar(A, B, key, extra):
        """0 <= A.y <= Xbar, using B.y = mu* A.y + eps with mu* = B_key/A_key exact and |eps| <= eta."""
        Ad, Bd = dict(zip(*A)), dict(zip(*B))
        if set(Ad) != set(Bd) or 0 in Ad or i_T0 in Ad:
            raise ValueError("Im rows must share support and carry no constant or T0")
        mu = Bd[key] / Ad[key]
        if not mu > 0:
            raise ValueError("non-positive Im-row ratio")
        eta = sum((_absup(Bd[k] - mu * Ad[k]) * extra(k) for k in Ad if k != key), ZERO)
        mlo = arb(mu.lower())
        return max(_up((1 + 2 * mlo * eta) / mlo ** 2), _up(eta / mlo), key=float), Ad

    def rho_only(k):
        if k < rho0:
            raise ValueError(f"unexpected variable {k} in an Im row")
        return rho_bar[k - rho0]

    idx_of = [(I, e) for I in (0, 1, 2) for e in range(I % 2, 2 * L, 2)]
    xbars = {}
    for k in range(M):                                 # (B2) I=2, l=0: A = c sigma2_k + rho terms
        j = idx_of.index((2, 0)) * M + k
        A, P, B = disk(j)
        Xb, Ad = xbar(A, B, s2 + k, rho_only)
        if {v for v in Ad if v < rho0} != {s2 + k}:
            raise ValueError(f"block {j}: I=2 Im row does not isolate sigma2_{k}")
        gam = sum((_absup(v) * rho_bar[v_ - rho0] for v_, v in Ad.items() if v_ >= rho0), ZERO)
        ybar[s2 + k] = _up((Xb + gam) / abs(Ad[s2 + k]))     # |A_s sigma| <= |A.y| + gam <= Xbar + gam
        xbars[j] = Xb
    for k in range(M):                                 # (B2) I=0, l=0: A = 3/2 c sigma1_k + c sigma2_k + rho terms
        j = idx_of.index((0, 0)) * M + k
        A, P, B = disk(j)
        Xb, Ad = xbar(A, B, s1 + k, lambda v: ybar[s2 + k] if v == s2 + k else rho_only(v))
        if {v for v in Ad if v < rho0} != {s1 + k, s2 + k}:
            raise ValueError(f"block {j}: I=0 Im row does not isolate sigma1_{k}, sigma2_{k}")
        gam = _absup(Ad[s2 + k]) * ybar[s2 + k] + sum((_absup(v) * rho_bar[v_ - rho0]
                                                       for v_, v in Ad.items() if v_ >= rho0), ZERO)
        ybar[s1 + k] = _up((Xb + gam) / abs(Ad[s1 + k]))
        xbars[j] = Xb
    best = None                                        # (B3) T0 from |P.y| <= sqrt(Xbar)
    for j, Xb in xbars.items():
        Pd = dict(zip(*disk(j)[1]))
        if i_T0 not in Pd:
            continue
        rest = sum((_absup(v) * (arb(1) if k == 0 else ybar[k]) for k, v in Pd.items() if k != i_T0), ZERO)
        t = _up((Xb.sqrt() + rest) / abs(Pd[i_T0]))     # ball division: rigorous upper end
        if best is None or float(t) < float(best[0]):
            best = (t, j)
    ybar[i_T0] = best[0]
    if any(v is None for v in ybar[1:]):
        raise ValueError("unbounded free variable")
    return ybar[1:], {"T0_block": best[1], "T0": float(best[0]),
                      "sigma1_range": [min(float(v) for v in ybar[s1:s1 + M]), max(float(v) for v in ybar[s1:s1 + M])],
                      "sigma2_range": [min(float(v) for v in ybar[s2:s2 + M]), max(float(v) for v in ybar[s2:s2 + M])],
                      "rho": [min(float(v) for v in rho_bar), max(float(v) for v in rho_bar)]}


def trace_bounds(sparse, ybar):
    """tau_j >= tr M_j(y): sum over diagonal entries of |const| + sum_n |coef_n| ybar_n (Arb upper ends)."""
    out = []
    for sp in sparse:
        t = ZERO
        for (r, c), (idx, vals) in sp.items():
            if r == c:
                t += sum((_absup(v) * (arb(1) if k == 0 else ybar[k - 1]) for k, v in zip(idx, vals)), ZERO)
        out.append(_up(t))
    return out


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
        Z = read_x(outdir, sizes, mapping)
        r = stationarity(objective, sparse, Z)
        res[mapping] = {"dual_objective": r[0].mid().str(40, radius=False),
                        "max_abs_residual": max(float(abs(x).upper()) for x in r[1:]),
                        "minus_primalObjective": float((r[0] - arb(solver["primalObjective"])).mid())}
    ok = res["colmajor_halved"]
    passed = ok["max_abs_residual"] < 1e-30 and abs(ok["minus_primalObjective"]) < 1e-30
    return {"mappings": res, "accepted_mapping": "colmajor_halved" if passed else None, "passed": passed,
            "sdpb": {k: solver.get(k) for k in ("primalObjective", "dualObjective", "dualityGap", "primalError", "dualError")}}


def certify(full_leaf, z_leaf, out_json, paper=None, refine=8):
    t0 = time.time()
    full_leaf, z_leaf = Path(full_leaf), Path(z_leaf)
    rep = json.loads((full_leaf / "report.json").read_text())
    spec = rep["spec"]
    if spec.get("uv"):
        raise NotImplementedError("UV leaf: rho_hat above s0 has no a-priori bound; the 2x2 |S|<=1 relaxation is not "
                                  "implemented in this script")
    if spec.get("reduce_basis"):
        raise ValueError("blocks must come from a full-space leaf (a reduced PMP is a restriction of the problem)")
    objective, blocks, sparse = load_pmp(full_leaf / "pmp.json")
    sizes = [b.size for b in blocks]
    t_parse = time.time() - t0
    print(f"parsed {len(blocks)} blocks, n={len(objective)} in {t_parse:.0f}s", flush=True)
    out = {"leaf_blocks": str(full_leaf), "leaf_Z": str(z_leaf), "scope": "emitted full-space degree-0 PMP "
           "(30-digit decimals parsed exactly into Arb balls); chiral-only leaf, so no rho_hat variables and no "
           "relaxation is used", "n_vars": len(objective), "n_blocks": len(blocks), "parse_seconds": t_parse}
    solver_z = sdpb_out(z_leaf)
    if z_leaf == full_leaf:
        out["P0"] = p0(objective, sparse, sizes, z_leaf / "out", solver_z)
        if not out["P0"]["passed"]:
            Path(out_json).write_text(json.dumps(out, indent=1))
            raise SystemExit("P0 failed: x_j -> Z mapping not confirmed")
    Z = read_x(z_leaf / "out", sizes, "colmajor_halved")
    if z_leaf != full_leaf:                            # P2: raw full-space stationarity of a foreign Z
        r = stationarity(objective, sparse, Z)
        out["raw_foreign_Z"] = {"dual_objective": float(r[0].mid()), "Z_leaf_primalObjective": solver_z["primalObjective"],
                                "max_abs_residual": max(float(abs(x).upper()) for x in r[1:]),
                                "sum_abs_residual": sum(float(abs(x).upper()) for x in r[1:])}
    print(f"Z read; P0={out.get('P0', {}).get('passed')} at {time.time() - t0:.0f}s", flush=True)
    ybar, binfo = a_priori_bounds(sparse, len(objective), spec["M"], spec["L"])
    tau = trace_bounds(sparse, ybar)
    print(f"a-priori bounds {binfo} at {time.time() - t0:.0f}s", flush=True)
    t1 = time.time()
    res = dual_bound(blocks, objective, Z, tau, ybar, refine=refine)      # refine=0: no Gram repair of Z
    out["gram_repair_iterations"] = refine
    ub = res["bound"]
    out.update({"a_priori_bounds": binfo, "tau_max": float(max(tau, key=float)),
                "certified": bool(ub.is_finite()), "bound_upper": upper_str(ub) if ub.is_finite() else "inf",
                "bound_mid": float(ub.mid()) if ub.is_finite() else None,
                "residual_history": res.get("residual_history"), "residual_max_after_repair": res.get("residual_max"),
                "max_psd_shift": res.get("max_shift"), "failed_blocks": res.get("failed_blocks"),
                "dual_bound_seconds": time.time() - t1, "wall_seconds": time.time() - t0})
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


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] not in ("p1", "p2"):
        raise SystemExit(__doc__)
    paper = float(a[a.index("--paper") + 1]) if "--paper" in a else None
    refine = int(a[a.index("--refine") + 1]) if "--refine" in a else 8
    if a[0] == "p1":
        r = certify(a[1], a[1], a[2], paper, refine)
    else:
        r = certify(a[1], a[2], a[3], paper, refine)
    print(json.dumps({k: r.get(k) for k in ("bound_upper", "certified", "max_psd_shift", "wall_seconds", "comparison")}, indent=1))
