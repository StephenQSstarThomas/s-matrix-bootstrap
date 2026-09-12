#!/usr/bin/env python
"""Smallest attainable eps^FF of (3.75), by constraint generation.

The paper fixes eps^FF = 6e-5 and says only that it should "reduce the allowed
space without making the problem unfeasible".  This script measures the boundary
of that statement directly: minimise the multiplier t with

    |cF_ell(s_i)| <= t * bound_ell   for the nodes above s0,

subject to the chiral constraints (3.64), the Gram matrices (3.68) and the FESR
(3.73).  Then the smallest attainable tolerance is eps^FF_min = 6e-5 * t^2, and
t <= 1 means the paper's value is attainable at that resolution.

The unitarity disks are handled by the same constraint generation used
elsewhere, so each t is certified: the relaxation optimum whose returned point
satisfies all 3 L M disks of the full problem.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings

import cvxpy as cp
import numpy as np

warnings.filterwarnings("ignore")

from smatrix_bootstrap.sdp import constraints as C
from smatrix_bootstrap.sdp import formfactor as FFM
from smatrix_bootstrap.sdp.problem import Model, ModelSpec, Operators

EPS_FF = C.EPS_FF


def one(M, L, frozen, time_limit, start_tol=1e-6, max_rounds=8):
    ops = Operators(M, L)
    ops.set_cone_scaling("rownorm")
    nu = ops.nu_measured
    mask = nu > start_tol * nu.max()
    rounds = []
    for _ in range(max_rounds):
        spec = ModelSpec(M=M, L=L, cone_scaling="rownorm", B=3.775e5, B_norm="l4",
                         chiral=True, uv=True, uv_parts=("gram", "fesr"),
                         ff_frozen_at_s0=frozen, disk_mask=mask.copy())
        m = Model(spec)
        K = FFM.hilbert_kernel(M)
        idx, ffb, ffk = C.ff_asymptotic_bounds(M, frozen_at_s0=frozen)
        t = cp.Variable(nonneg=True)
        extra = []
        for ell in (0, 1):
            g = FFM.gram_scale(ell, ops.s)
            kin = FFM.kinematic_factor(ell, ops.s)
            ReF = 1.0 + K @ m.ImF[ell]
            for n, i in enumerate(idx):
                extra.append(cp.SOC(t * (ffb[ell] / ffk[ell][n]),
                                    cp.hstack([kin[i] / g[i] * ReF[i],
                                               kin[i] / g[i] * m.ImF[ell][i]])))
        prob = cp.Problem(cp.Minimize(t), m.constraints + extra)
        try:
            prob.solve(solver="CLARABEL", max_iter=800, time_limit=time_limit,
                       tol_gap_abs=1e-8, tol_gap_rel=1e-8, tol_feas=1e-8)
        except Exception as exc:
            rounds.append({"n_disks": int(mask.sum()), "status": f"error: {exc!s:.60}"})
            return {"M": M, "L": L, "frozen_at_s0": frozen, "status": "failed",
                    "rounds": rounds}
        if t.value is None:
            rounds.append({"n_disks": int(mask.sum()), "status": prob.status})
            return {"M": M, "L": L, "frozen_at_s0": frozen, "status": prob.status,
                    "rounds": rounds}
        c = m.solution()["c"]
        hre, him = ops.h_re @ c, ops.h_im @ c
        mag2 = hre ** 2 + him ** 2
        bad = (mag2 - 2.0 * him) / np.maximum(mag2, 1e-300) > 1e-8
        rounds.append({"n_disks": int(mask.sum()), "status": prob.status,
                       "t": float(t.value), "n_violated": int(bad.sum())})
        if not bad.any():
            return {"M": M, "L": L, "frozen_at_s0": frozen, "status": prob.status,
                    "certified": True, "t": float(t.value),
                    "eps_ff_min": EPS_FF * float(t.value) ** 2,
                    "eps_ff_paper": EPS_FF, "rounds": rounds,
                    "n_disks_imposed": int(mask.sum())}
        mask |= bad
    return {"M": M, "L": L, "frozen_at_s0": frozen, "status": prob.status,
            "certified": False, "t": float(t.value),
            "eps_ff_min": EPS_FF * float(t.value) ** 2, "rounds": rounds}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--grid", default="20:6,25:8,30:8,35:10,40:10,50:10")
    p.add_argument("--time-limit", type=float, default=2400.0)
    p.add_argument("--out", default=None)
    a = p.parse_args()
    rows = []
    for spec in a.grid.split(","):
        M, L = (int(x) for x in spec.split(":"))
        for frozen in (True, False):
            t0 = time.time()
            r = one(M, L, frozen, a.time_limit)
            r["seconds"] = time.time() - t0
            rows.append(r)
            print(json.dumps(r, default=float), flush=True)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump({"paper_eps_ff": EPS_FF, "rows": rows}, fh, indent=1, default=float)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
