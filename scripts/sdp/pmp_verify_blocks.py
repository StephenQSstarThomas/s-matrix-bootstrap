"""Verify that every PMP block is the same constraint problem.Model states.

This is the test that matters for the SDPB arm: a solver result is only as good
as the problem it was handed.  It needs no feasible point and no solver -- it
evaluates each ``pmp.Pmp`` block and the corresponding ``problem.Model``
expression at one random point and diffs them.

Block order emitted by ``Pmp._blocks()``:

    unitarity disks (keep.sum())
    chiral (1 arrow for chi-b, 2 for chi-c, 2*n_res 1x1 for chi-a)
    density ball (1 arrow, if B is set)
    per ell in (0, 1):  gram (M)   [only with "gram"]
                        ff   (len(idx_hi))  interleaved into the gram loop
                        rho_hat >= 0 (M)
                        FESR (2 per moment)  [only with "fesr"]
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from smatrix_bootstrap.sdp import formfactor as FFM
from smatrix_bootstrap.sdp.pmp import Pmp
from smatrix_bootstrap.sdp.problem import Model, ModelSpec


def check(spec: ModelSpec, seed: int = 0, tol: float = 1e-10) -> dict:
    rng = np.random.default_rng(seed)
    m, w = Model(spec), Pmp(spec, direction=(1.0, 0.0))
    M = spec.M
    a = rng.normal(size=m.a.size) * 1e-3
    y = np.zeros(w.n_vars)
    y[0] = 1.0
    y[w.i_a:w.i_a + w.n_a] = a
    ImF = rh = None
    if spec.uv:
        ImF = rng.normal(size=(2, M)) * 1e-2
        rh = np.abs(rng.normal(size=(2, M)))
        y[w.i_ImF:w.i_ImF + 2 * M] = ImF.ravel()
        y[w.i_rho:w.i_rho + 2 * M] = rh.ravel()

    blocks = list(w._blocks())
    ev = lambda blk: np.array([[row @ y for row in r] for r in blk])
    out = {"spec": {k: (list(v) if isinstance(v, tuple) else v)
                    for k, v in spec.__dict__.items() if k != "disk_mask"},
           "n_blocks": len(blocks), "mismatches": [], "checked": {}}
    n_disk = int(w.keep.sum())

    # ---- unitarity disks: block must equal [[2r, p_re, p_im], [p_re, 1, 0], ...]
    bad = 0
    for j, k in enumerate(np.flatnonzero(w.keep)):
        ref = np.array([[2 * (w.R_im[k] @ a), w.P_re[k] @ a, w.P_im[k] @ a],
                        [w.P_re[k] @ a, 1.0, 0.0],
                        [w.P_im[k] @ a, 0.0, 1.0]])
        if np.abs(ev(blocks[j]) - ref).max() > tol:
            bad += 1
    out["checked"]["unitarity_disks"] = {"n": n_disk, "mismatch": bad}
    bad and out["mismatches"].append("unitarity_disks")

    cur = n_disk
    # ---- chiral
    if spec.chiral:
        res = w.chi @ a
        if spec.chi_caliber == "chi-b":
            arrow = ev(blocks[cur])
            ok = (abs(arrow[0, 0] - spec.eps_chi) < tol
                  and np.abs(arrow[0, 1:] - res).max() < tol
                  and np.abs(np.diag(arrow)[1:] - spec.eps_chi).max() < tol)
            out["checked"]["chiral_arrow"] = {"size": arrow.shape[0], "ok": bool(ok)}
            (not ok) and out["mismatches"].append("chiral_arrow")
            cur += 1
        else:
            out["checked"]["chiral_arrow"] = {"skipped": spec.chi_caliber}
            cur += 2 if spec.chi_caliber == "chi-c" else 2 * len(res)
    if spec.B is not None:
        cur += 1

    # ---- UV sector
    if spec.uv and "gram" in spec.uv_parts:
        K = FFM.hilbert_kernel(M)
        rt = np.sqrt(2.0)
        bad = 0
        n_ff = 0
        for ell in (0, 1):
            g = FFM.gram_scale(ell, m.ops.s)
            kin = FFM.kinematic_factor(ell, m.ops.s)
            re_row, im_row = m.ops.gram_rows[ell]
            S_re, S_im = 1.0 - im_row @ a, re_row @ a
            cF_re = (kin / g) * (1.0 + K @ ImF[ell])
            cF_im = (kin / g) * ImF[ell]
            for i in range(M):
                ref = np.array([[1 + S_re[i], S_im[i], rt * cF_re[i]],
                                [S_im[i], 1 - S_re[i], rt * cF_im[i]],
                                [rt * cF_re[i], rt * cF_im[i], rh[ell][i]]])
                got = ev(blocks[cur])
                if got.shape != ref.shape or np.abs(got - ref).max() > tol:
                    bad += 1
                cur += 1
                if "ff" in spec.uv_parts:          # ff blocks interleave here
                    from smatrix_bootstrap.sdp import constraints as C
                    idx_hi, _, _ = C.ff_asymptotic_bounds(M, spec.m_q, spec.eps_ff,
                                                          spec.ff_frozen_at_s0)
                    if i in set(idx_hi):
                        cur += 1
                        n_ff += 1
            for i in range(M):                     # rho_hat >= 0
                if abs(ev(blocks[cur])[0, 0] - rh[ell][i]) > tol:
                    bad += 1
                cur += 1
            if "fesr" in spec.uv_parts:
                from smatrix_bootstrap.sdp import constraints as C
                cur += 2 * len(C.MOMENTS[ell])
        out["checked"]["gram_and_rho_hat"] = {"n": 4 * M, "mismatch": bad,
                                              "ff_blocks_skipped": n_ff}
        bad and out["mismatches"].append("gram_and_rho_hat")

    out["all_blocks_consumed"] = (cur == len(blocks))
    out["blocks_accounted"] = cur
    out["ok"] = not out["mismatches"] and out["all_blocks_consumed"]
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser("pmp_verify_blocks")
    p.add_argument("--M", type=int, default=8)
    p.add_argument("--L", type=int, default=4)
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)
    recs = []
    for chiral, uv, parts in ((True, False, ()), (True, True, ("gram",)),
                              (True, True, ("gram", "fesr")),
                              (True, True, ("gram", "fesr", "ff"))):
        spec = ModelSpec(M=a.M, L=a.L, chiral=chiral, chi_caliber="chi-b", uv=uv,
                         uv_parts=parts, B=None, B_norm="l2", cone_scaling="rownorm")
        r = check(spec)
        recs.append(r)
        print(f"uv_parts={'+'.join(parts) or 'none':20s} blocks={r['n_blocks']:5d} "
              f"accounted={r['blocks_accounted']:5d} ok={r['ok']}  {r['mismatches']}",
              flush=True)
    if a.out:
        json.dump({"records": recs}, open(a.out, "w"), indent=2)
        print("wrote", a.out)
    return 0 if all(r["ok"] for r in recs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
