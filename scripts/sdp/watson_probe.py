"""S4b: measure the Watson-saturation degeneracy of the UV sector (blocker B2).

The real 3x3 form of (3.68) built in ``problem._gram_psd`` is

    [[1+S_re,   S_im,     sqrt2 F_re],
     [S_im,     1-S_re,   sqrt2 F_im],
     [sqrt2 F_re, sqrt2 F_im, rho]]  >= 0

whose leading 2x2 minor is exactly ``1 - |S|^2``.  At unitarity saturation the
2x2 is singular with null vector ``(-S_im, 1+S_re)``, so PSD forces

    W := -S_im * F_re + (1+S_re) * F_im  ->  0,

i.e. ``arg F = delta`` -- Watson's theorem.  In a conic program that is an
*equality hidden inside a PSD constraint*: where the paper expects saturation
the feasible set has no interior and Slater fails.

This script does not decide whether the paper's UV problem is strictly feasible
-- that is the physics question (TASK_PHYSICS_CALIBERS_ZH.md Q0).  It measures
what is measurable: can a solver close the gap, and how tightly is W locked.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import cvxpy as cp
import numpy as np

from smatrix_bootstrap.sdp import constraints as C
from smatrix_bootstrap.sdp import formfactor as FFM
from smatrix_bootstrap.sdp.problem import Model, ModelSpec


def _solve(model, direction, solver, feasibility, max_iter, time_limit):
    if feasibility:
        prob = cp.Problem(cp.Maximize(cp.Constant(0.0)), model.constraints)
    else:
        model.direction.value = np.asarray(direction, dtype=float)
        prob = model.finalize(objective="plane")
    kw = {"verbose": False}
    if solver == "MOSEK":
        kw["mosek_params"] = {
            "MSK_IPAR_INTPNT_MAX_ITERATIONS": max_iter,
            "MSK_DPAR_OPTIMIZER_MAX_TIME": time_limit,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-10,
        }
    else:
        kw.update(max_iter=max_iter, time_limit=time_limit)
    t0 = time.time()
    try:
        prob.solve(solver=solver, **kw)
        err = None
    except Exception as exc:                    # solver blow-up is a result here
        err = f"{type(exc).__name__}: {exc}"
    return prob, time.time() - t0, err


def _watson(model) -> dict:
    """|S|, the Watson residual W and rho on the returned solution, per wave."""
    out = {}
    if model.ImF is None or model.ImF.value is None:
        return out
    M = model.spec.M
    K = FFM.hilbert_kernel(M)
    a = model.a.value
    for ell, I in ((0, 0), (1, 1)):
        g = FFM.gram_scale(ell, model.ops.s)
        kin = FFM.kinematic_factor(ell, model.ops.s)
        re_row, im_row = (model.ops.gram_rows[ell] if model.basis is None else
                          tuple(x @ model.basis for x in model.ops.gram_rows[ell]))
        S_re = 1.0 - im_row @ a
        S_im = re_row @ a
        ImF = model.ImF.value[ell]
        ReF = 1.0 + K @ ImF
        F_re, F_im = (kin / g) * ReF, (kin / g) * ImF
        W = -S_im * F_re + (1.0 + S_re) * F_im
        absS = np.hypot(S_re, S_im)
        scale = np.maximum(np.hypot(F_re, F_im) * (1.0 + absS), 1e-300)
        out["S0" if ell == 0 else "P1"] = {
            "max_absS": float(absS.max()),
            "n_nodes_absS_gt_0.999": int((absS > 0.999).sum()),
            "n_nodes_absS_gt_0.99": int((absS > 0.99).sum()),
            "min_1_minus_absS": float((1.0 - absS).min()),
            "watson_residual_rel_max": float((np.abs(W) / scale).max()),
            "watson_residual_rel_at_most_saturated": float(
                (np.abs(W) / scale)[int(np.argmax(absS))]),
            "rho_hat_min": float(model.rho_hat.value[ell].min()),
        }
    return out


def run(M, L, parts, chiral, solver, feasibility, max_iter, time_limit, B):
    spec = ModelSpec(M=M, L=L, chiral=chiral, uv=True, uv_parts=tuple(parts),
                     B=B, B_norm="l4", cone_scaling="rownorm",
                     tag=f"watson-{'+'.join(parts)}")
    t0 = time.time()
    model = Model(spec)
    build = time.time() - t0
    prob, secs, err = _solve(model, (1.0, 0.0), solver, feasibility,
                             max_iter, time_limit)
    rec = {
        "M": M, "L": L, "uv_parts": list(parts), "chiral": chiral,
        "solver": solver, "mode": "feasibility" if feasibility else "extremal",
        "build_seconds": round(build, 2), "solve_seconds": round(secs, 2),
        "status": None if err else prob.status, "error": err,
        "objective": None if err or prob.value is None else float(prob.value),
        "n_unitarity_cones": model.n_unitarity_cones,
        "n_constraints": len(model.constraints),
    }
    st = getattr(prob, "solver_stats", None)
    if st is not None:
        rec["iterations"] = getattr(st, "num_iters", None)
    if not err and model.a.value is not None:
        rec["f00_3"] = float(model.f00.value)
        rec["f11_3"] = float(model.f11.value)
        rec["watson"] = _watson(model)
    return rec


def main(argv=None) -> int:
    p = argparse.ArgumentParser("watson_probe")
    p.add_argument("--M", type=int, nargs="+", default=[20, 30])
    p.add_argument("--L", type=int, default=8)
    p.add_argument("--solver", nargs="+", default=["MOSEK", "CLARABEL"])
    p.add_argument("--parts", nargs="+", default=["gram", "gram+fesr", "gram+fesr+ff"])
    p.add_argument("--chiral", action="store_true", default=True)
    p.add_argument("--no-chiral", dest="chiral", action="store_false")
    p.add_argument("--B", type=float, default=377500.0)
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--time-limit", type=float, default=1800.0)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    recs = []
    for M in a.M:
        for parts in a.parts:
            for solver in a.solver:
                for feas in (True, False):
                    rec = run(M, a.L, parts.split("+"), a.chiral, solver, feas,
                              a.max_iter, a.time_limit, a.B)
                    recs.append(rec)
                    print(json.dumps({k: rec[k] for k in
                                      ("M", "uv_parts", "solver", "mode", "status",
                                       "objective", "solve_seconds")}), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump({"records": recs}, fh, indent=2)
    print(f"\nwrote {a.out}  ({len(recs)} records)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
