"""V2: rigorous upper bound on the optimum of a 2403 leaf from SDPB's multipliers (VERIFICATION_AUDIT_ZH.md section 2).

Weak duality with a corrected multiplier and a-priori bounds (DERIVATION_ZH.md section 12):
    b.y <= b0 + sum_j <Z_j, C_j> + sum_n |r_n| ybar_n + sum_j max(0, -lambda_min(Z_j)) tau_j .
The certified problem is the intended one: block forms are rebuilt in Arb from the exact kernel export
(ExactKernels), never from the float64 PMP.  Valid for the full-space problem only (no basis reduction).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from flint import arb, ctx

def read_dual_matrices(outdir, sizes):
    """SDPB's X_matrix_{2j}.txt is the multiplier Z_j of PMP block j (degree 0: the odd parity block is empty)."""
    out = []
    for j, m in enumerate(sizes):
        lines = Path(outdir, f"X_matrix_{2 * j}.txt").read_text().split()
        r, c = int(lines[0]), int(lines[1])
        if (r, c) != (m, m):
            raise ValueError(f"block {j}: X is {r}x{c}, expected {m}x{m}")
        out.append([[lines[2 + a * m + b] for b in range(m)] for a in range(m)])
    return out


def _to_arb(v):
    if isinstance(v, arb):
        return v
    if isinstance(v, str):
        return arb(v)
    return arb(float(v))


def psd_shift(A):
    """Smallest mu in {0, 10^-k scale} with A + mu I certified positive definite by an Arb LDL^T; lambda_min >= -mu."""
    m = len(A)
    scale = max((abs(float(A[i][j].mid())) for i in range(m) for j in range(m)), default=1.0) or 1.0
    for mu in [arb(0)] + [arb(10) ** (-k) * arb(scale) for k in range(70, -3, -2)]:
        L, d, ok = [[arb(0)] * m for _ in range(m)], [arb(0)] * m, True
        for k in range(m):
            d[k] = A[k][k] + mu - sum((L[k][i] ** 2 * d[i] for i in range(k)), arb(0))
            if not d[k] > 0:
                ok = False
                break
            for j in range(k + 1, m):
                L[j][k] = (A[j][k] - sum((L[j][i] * L[k][i] * d[i] for i in range(k)), arb(0))) / d[k]
        if ok:
            return mu
    return None


def dual_bound(blocks, objective, Z, tau, ybar, refine: int = 8) -> dict:
    """Rigorous upper bound of max b.y s.t. sum_n y_n M_j^(n) >= 0 (y0 = 1) from approximate multipliers Z_j.

    b.y <= b0 + sum_j <Z_j, C_j> + sum_n |r_n| ybar_n + sum_j max(0, -lambda_min(Z_j)) tau_j, where
    r_n = b_n + sum_j <Z_j, A_jn>, |y_n| <= ybar_n and tr M_j(y) <= tau_j on the feasible set.  Z is first corrected
    in the span of the A_jn (float64 Gram solve, exact Arb re-evaluation, iterated) so that r is tiny.
    """
    ctx.prec = 256
    n = len(objective)
    b = [_to_arb(v) for v in objective]
    ents = []                       # (block, r, c, weight, nonzero indices, exact values there, float coefficients)
    for j, blk in enumerate(blocks):
        for (r, c), f in blk.forms.items():
            ex = np.array([_to_arb(v) for v in f], dtype=object)
            nz = np.array([k for k in range(n) if ex[k] != 0], dtype=int)
            fl = np.zeros(n)
            fl[nz] = [float(x.mid()) for x in ex[nz]]
            ents.append((j, r, c, 1 if r == c else 2, nz, ex[nz], fl))
    Zc = [[[_to_arb(x) for x in row] for row in z] for z in Z]

    def residual():
        r = np.array(b, dtype=object)
        for j, rr, cc, w, nz, vals, _ in ents:
            if nz.size:
                r[nz] = r[nz] + (w * Zc[j][rr][cc]) * vals
        return list(r)

    C = np.array([e[6][1:] for e in ents])
    G = C.T @ (np.array([e[3] for e in ents], dtype=float)[:, None] * C)
    dg = np.sqrt(np.diag(G))
    D = 1.0 / np.where(dg > 0, dg, 1.0)                 # Jacobi scaling: cond(G) ~ 1e13 at M = 10
    Gs = D[:, None] * G * D[None, :]
    history = []
    for it in range(refine + 1):
        r = residual()
        history.append(max(abs(float(x.mid())) for x in r[1:]))
        if it == refine or history[-1] < 1e-45:
            break
        rhs = -D * np.array([float(x.mid()) for x in r[1:]])
        try:
            lam = D * np.linalg.solve(Gs, rhs)
        except np.linalg.LinAlgError:
            lam = D * np.linalg.lstsq(Gs, rhs, rcond=1e-18)[0]
        for j, rr, cc, w, _, _, fl in ents:
            d = arb(float(fl[1:] @ lam))
            Zc[j][rr][cc] += d
            if rr != cc:
                Zc[j][cc][rr] += d
    mus = [psd_shift(z) for z in Zc]
    if any(mu is None for mu in mus):
        return {"bound": arb("inf"), "residual_history": history, "failed_blocks": [i for i, m in enumerate(mus) if m is None]}
    bound = r[0] + sum((abs(r[k]) * _to_arb(ybar[k - 1]) for k in range(1, n)), arb(0)) \
        + sum((mu * _to_arb(t) for mu, t in zip(mus, tau)), arb(0))
    return {"bound": bound, "residual_max": max(abs(float(x.mid())) + float(x.rad()) for x in r[1:]),
            "residual_history": history, "max_shift": max(float(m.mid()) for m in mus),
            "objective_constant_part": r[0]}


class ExactKernels:
    """A kernel export whose arrays are object arrays of Arb balls (1-ulp radius at the 40th digit)."""

    exact_mode = True

    def __init__(self, ks):
        self._ks, self.manifest, self._cache = ks, ks.manifest, {}

    def __getattr__(self, name):
        if name.startswith("_") or name not in self.manifest.get("matrices", {}):
            raise AttributeError(name)
        if name not in self._cache:
            from .certify import ball
            shape = tuple(self.manifest["matrices"][name]["shape"])
            arr = np.empty(int(np.prod(shape)), dtype=object)
            arr[:] = [ball(t) for t in self._ks.exact(name)]
            self._cache[name] = arr.reshape(shape)
        return self._cache[name]


def _abs_upper(x):
    return arb(abs(x).upper()) if isinstance(x, arb) else arb(abs(float(x)))


def free_variable_bounds(prog) -> list:
    """|y_n| <= ybar_n for every free variable of a full-space exact program (DERIVATION section 12, B0-B6)."""
    ks, M, n0, sp = prog.ks, prog.M, prog.n0, prog.spec
    n_rho = prog.groups["rho"].stop - prog.groups["rho"].start
    R = getattr(prog.reg, "bound", None)
    if R is None:
        raise ValueError("V2 needs a regularised program (a-priori bound on rho)")
    rb = [arb(repr(float(r))) for r in R]              # |rho_i| <= radius_i (uniform for linf-relax)
    ybar = [arb(0)] * (prog.n_vars - 1)
    rho0 = prog.groups["rho"].start
    for i in range(n_rho):
        ybar[rho0 + i - 1] = rb[i]
    if all(float(r) == float(R[0]) for r in R):
        beta = [sum((_abs_upper(v) for v in row), arb(0)) * rb[0] for row in ks.Im_rho]   # |Im_rho . rho|
    else:
        beta = [sum((_abs_upper(v) * r for v, r in zip(row, rb)), arb(0)) for row in ks.Im_rho]
    row = prog.row_index
    s2, s1 = [None] * M, [None] * M
    for m in range(M):                     # Im h^{I=2}_0 = c2 sigma2_m + rho terms, 0 <= Im h <= 2 (B2)
        n = row("S2", m)
        c2 = ks.Im_sig[n][1 + M + m]
        s2[m] = (2 + beta[n]) / abs(c2)
        n = row("S0", m)
        c1, c2b = ks.Im_sig[n][1 + m], ks.Im_sig[n][1 + M + m]
        s1[m] = (2 + abs(c2b) * s2[m] + beta[n]) / abs(c1)
        ybar[m] = s1[m]                    # sigma1_m is free variable index m (amp index 1 + m)
        ybar[M + m] = s2[m]
    # B-spline coefficients from the FESR balls (B5) and Im F from the Gram blocks / caps (B6)
    nb, chans = prog.nb, ("S0", "P1", "D0")
    cbar = [arb(0)] * (3 * nb)
    for k, nm in enumerate(("intS0", "intP1", "intD0")):
        mom = getattr(ks, nm) @ ks.BSpline3[k * n0:(k + 1) * n0]           # n_SR x 120, nonnegative entries
        nsr = len(mom)
        assert 3 * nsr == len(ks.SVZSR), "moment rows and SVZ targets disagree"
        eps = prog.num(sp.eps_sr[k])
        for i in range(k * nb, (k + 1) * nb):
            cands = [(ks.SVZSR[nsr * k + q] + eps) / mom[q][i] for q in range(nsr) if mom[q][i] > 0]
            cbar[i] = min(cands, key=lambda z: float(z.upper())) if cands else arb(0)
    for j, i in enumerate(prog.c_free):
        ybar[prog.o_c - 1 + j] = cbar[i]
    caps = {0: lambda j: prog.num(sp.r[0]), 1: lambda j: prog.num(sp.r[1]) * prog.num(sp.fp1asym) / ks.x[j],
            2: lambda j: prog.num(sp.r[2]) * prog.num(sp.fd0asym) / ks.x[j]}
    for k in range(3):
        for j in range(M):
            if j < n0:
                rhat = sum((ks.BSpline3[k * n0 + j][i] * (arb(1) if i in (0, nb, 2 * nb) else cbar[i])
                            for i in range(k * nb, (k + 1) * nb)), arb(0))
                bnd = rhat.sqrt()          # |F|^2 <= rho_hat from the (1,3) principal minor of the Gram block
            else:
                bnd = caps[k](j)
            ybar[prog.o_imf - 1 + k * M + j] = bnd
    for a in range(prog.reg.n_aux):        # l4-exact auxiliaries u, w in [0, 1]
        ybar[prog.o_aux - 1 + a] = arb(1)
    return ybar


def trace_bounds(blocks, ybar) -> list:
    """tr M_j(y) <= sum over diagonal entries of |const| + sum_n |coef_n| ybar_n."""
    out = []
    for blk in blocks:
        t = arb(0)
        for r in range(blk.size):
            f = blk.entry(r, r)
            if f is None:
                continue
            t += _abs_upper(f[0]) + sum((_abs_upper(f[k]) * ybar[k - 1] for k in range(1, len(f)) if f[k] != 0), arb(0))
        out.append(t)
    return out


def certify_leaf(leaf_dir, kernels_dir) -> dict:
    """V1 lower bound (certified feasible point) and V2 upper bound for one full-space leaf."""
    import json
    from dataclasses import fields
    from .certify import primal
    from .kernels import load
    from .program import Program
    from .spec import GtbSpec
    leaf = Path(leaf_dir)
    rep = json.loads((leaf / "report.json").read_text())
    full_space = not rep.get("reduce_basis")
    spec = GtbSpec(**{f.name: (tuple(rep["spec"][f.name]) if isinstance(rep["spec"][f.name], list) else rep["spec"][f.name])
                      for f in fields(GtbSpec)})
    p1 = primal(leaf, kernels_dir)
    ok1 = all(v["certified_fail"] == 0 and v["inconclusive"] == 0 for v in p1["families"].values())
    lower = {"value": p1["objective_F1"]["mid"], "radius": p1["objective_F1"]["radius"], "certified": ok1}
    upper = {"value": None, "certified": False, "reason": "reduced basis: V2 needs the full-space problem"}
    if full_space:
        targets = None
        if (leaf / "watson_targets.npz").exists():
            z = np.load(leaf / "watson_targets.npz")
            targets = {w: (z[f"{w}_re"], z[f"{w}_im"]) for w in ("S0", "P1", "D0", "S2", "D2", "F1")}
        orow = np.load(leaf / "objective_row.npy") if (leaf / "objective_row.npy").exists() else None
        prog = Program(spec, ExactKernels(load(kernels_dir)), reg=rep["regulariser"], watson_targets=targets,
                       objective_row=orow)
        blocks = prog.blocks()
        Z = read_dual_matrices(leaf / "out", [b.size for b in blocks])
        ybar = free_variable_bounds(prog)
        res = dual_bound(blocks, list(prog.objective()), Z, trace_bounds(blocks, ybar), ybar)
        ub = res["bound"]
        upper = {"value": ub.upper().str(30, radius=False) if ub.is_finite() else "inf",
                 "certified": bool(ub.is_finite()), "residual_history": res.get("residual_history"),
                 "max_shift": res.get("max_shift")}
    return {"lower": lower, "upper": upper, "regulariser": rep["regulariser"], "l4_feasible": p1["l4_feasible"],
            "families": p1["families"], "rho_l4": p1["rho_l4"], "l4_bound": p1["l4_bound"],
            "objective": rep["spec"]["objective"], "f0_ni": rep["spec"]["f0_ni"]}


def certify_infeasible(leaf_dir, kernels_dir) -> dict:
    """V2 infeasibility certificate from SDPB's multipliers of a (diverging) leaf.

    For every feasible y:  0 <= sum_j <Z_j, M_j(y)> = sum_j <Z_j, C_j> + sum_n y_n r_n,  r_n = sum_j <Z_j, A_jn>.
    dual_bound with a zero objective returns  sum<Z, C> + sum |r_n| ybar_n + PSD shifts  (Z repaired so that r is
    tiny); a value below 0 proves that no feasible y exists.  A ray has no scale: Z is normalised to max entry 1.
    """
    import json
    from dataclasses import fields
    from .kernels import load
    from .program import Program
    from .spec import GtbSpec
    leaf = Path(leaf_dir)
    rep = json.loads((leaf / "report.json").read_text())
    spec = GtbSpec(**{f.name: (tuple(rep["spec"][f.name]) if isinstance(rep["spec"][f.name], list) else rep["spec"][f.name])
                      for f in fields(GtbSpec)})
    prog = Program(spec, ExactKernels(load(kernels_dir)), reg=rep["regulariser"])
    blocks = prog.blocks()
    Z = read_dual_matrices(leaf / "out", [b.size for b in blocks])
    scale = max(abs(float(x)) for z in Z for row in z for x in row) or 1.0
    Z = [[[arb(x) / arb(repr(scale)) for x in row] for row in z] for z in Z]
    ybar = free_variable_bounds(prog)
    res = dual_bound(blocks, [0] * prog.n_vars, Z, trace_bounds(blocks, ybar), ybar)
    b = res["bound"]
    return {"certified_infeasible": bool(b.is_finite() and b < 0),
            "bound": b.upper().str(20, radius=False) if b.is_finite() else "inf",
            "residual_history": res.get("residual_history"), "max_shift": res.get("max_shift")}
