"""A posteriori verification of a returned solution against the exact constraints.

Every conditioning device used by :mod:`smatrix_bootstrap.sdp.problem` (the
per-row cone rescaling, the sparsification, the basis reduction, the Gram
congruence) is exact in theory but approximate in float64.  Nothing is trusted:
a solution is always replayed through the *unmodified* operators here, and the
worst violation of each constraint family is reported alongside the objective.

Unitarity is measured by the margin of (2.12),

    m = 2 Im h - |h|^2 ,   h = kappa f ,

which must be >= 0; we report it both absolutely and relative to the row scale.
"""
from __future__ import annotations

import numpy as np

from . import constraints as C
from . import formfactor as FFM
from .grid import quad_weights, s_nodes


def unitarity_report(ops, c: np.ndarray) -> dict:
    """Worst violation of |h|^2 <= 2 Im h over all (I, ell, i), unrescaled."""
    hre, him = ops.h_re @ c, ops.h_im @ c
    margin = 2.0 * him - (hre ** 2 + him ** 2)
    # The primary metric is eta = |S| <= 1 of (2.12): it is scale free, so it is
    # not swamped by the centrifugally suppressed rows where both |h|^2 and
    # Im h are ~1e-40 and the margin is meaninglessly small either way.
    eta = np.sqrt(hre ** 2 + (him - 1.0) ** 2)
    scale = np.maximum(np.maximum(hre ** 2 + him ** 2, 2.0 * np.abs(him)), 1e-300)
    rel = margin / scale
    # The discriminating metric: relative violation restricted to the disks that
    # actually carry amplitude.  `max eta - 1` is blind here because the
    # centrifugally suppressed disks have h ~ 0, hence eta == 1 exactly, and they
    # dominate the maximum; `min margin` is blind for the mirror reason.
    mag = np.sqrt(hre ** 2 + him ** 2)
    active = mag > 1e-6
    rel_v = ((mag ** 2 - 2.0 * him)[active] / np.maximum(mag[active] ** 2, 1e-300)
             if active.any() else np.array([-1.0]))
    # A disk with |h| below the active cut can still violate the constraint in
    # absolute terms (Im h slightly negative), which the relative metric would
    # divide away, so feasibility needs both tests.
    abs_v = float((mag ** 2 - 2.0 * him).max())
    nw = len(ops.index)
    return {"max_relative_violation_active": float(rel_v.max()),
            "max_absolute_violation_all": abs_v,
            "n_active_disks": int(active.sum()),
            "max_abs_h": float(mag.max()),
            "feasible": bool(rel_v.max() <= 1e-6 and abs_v <= 1e-8),
            "max_eta_minus_1": float(eta.max() - 1.0),
            "max_eta_wave": _label(ops, int(np.argmax(eta)), nw),
            "n_rows_eta_above_1p1e_minus_8": int((eta > 1 + 1e-8).sum()),
            "min_margin": float(margin.min()),
            "min_margin_wave": _label(ops, int(np.argmin(margin)), nw),
            "min_margin_relative": float(rel.min()),
            "min_margin_relative_wave": _label(ops, int(np.argmin(rel)), nw),
            "n_rows": int(margin.size)}


def _label(ops, k, nw):
    a, node = divmod(k, ops.M)
    I, ell = ops.index[a]
    return {"isospin": I, "ell": ell, "node": node, "s": float(ops.s[node])}


def gram_report(ops, c, ImF, rho_hat) -> dict:
    """Smallest eigenvalue of the *unrescaled* (3.68) block at every node."""
    M = ops.M
    K = FFM.hilbert_kernel(M)
    worst = np.inf
    where = None
    for ell in (0, 1):
        kin = FFM.kinematic_factor(ell, ops.s)
        ReF = 1.0 + K @ ImF[ell]
        cF = kin * (ReF + 1j * ImF[ell])
        rho = rho_hat[ell] * FFM.gram_scale(ell, ops.s) ** 2
        hre, him = ops.gram_rows[ell][0] @ c, ops.gram_rows[ell][1] @ c
        S = (1.0 - him) + 1j * hre
        for i in range(M):
            ev = np.linalg.eigvalsh(FFM.gram_block(S[i], cF[i], rho[i])).min()
            if ev < worst:
                worst, where = float(ev), {"ell": ell, "node": i, "s": float(ops.s[i])}
    return {"min_eigenvalue": worst, "at": where}


def fesr_report(ops, rho_hat, caliber: str) -> dict:
    M = ops.M
    tgt, tol = C.printed_targets(), C.sr_tolerances(caliber)
    rows = []
    for ell, wave in ((0, "S0"), (1, "P1")):
        rho = rho_hat[ell] * FFM.gram_scale(ell, ops.s) ** 2
        for n in C.MOMENTS[ell]:
            m = float(C.moment_row(M, n) @ rho)
            rows.append({"wave": wave, "n": n, "moment": m, "target": tgt[(wave, n)],
                         "residual": m - tgt[(wave, n)], "tolerance": tol[(wave, n)],
                         "violation": max(0.0, abs(m - tgt[(wave, n)]) - tol[(wave, n)])})
    return {"caliber": caliber, "rows": rows,
            "max_violation": max(r["violation"] for r in rows)}


def ff_report(ops, ImF, m_q, eps_ff, frozen_at_s0=True) -> dict:
    M = ops.M
    K = FFM.hilbert_kernel(M)
    idx, bnd, kuse = C.ff_asymptotic_bounds(M, m_q, eps_ff, frozen_at_s0)
    worst, where = 0.0, None
    for ell in (0, 1):
        F = (1.0 + K @ ImF[ell]) + 1j * ImF[ell]
        cF = {i: kuse[ell][n] * F[i] for n, i in enumerate(idx)}
        for i in idx:
            v = abs(cF[i]) - bnd[ell]
            if v > worst:
                worst, where = float(v), {"ell": ell, "node": int(i), "bound": bnd[ell],
                                          "value": float(abs(cF[i]))}
    return {"max_violation": worst, "at": where, "n_constraints": 2 * len(idx)}


def chiral_report(ops, c, caliber: str, eps: float) -> dict:
    r = ops.chi_rows @ c
    if caliber == "chi-a":
        used = float(np.abs(r).max())
    elif caliber == "chi-b":
        used = float(np.linalg.norm(r))
    elif caliber == "chi-c":
        used = float(max(np.linalg.norm(r[0::2]), np.linalg.norm(r[1::2])))
    else:
        raise ValueError(caliber)
    return {"caliber": caliber, "eps": eps, "norm_used": used,
            "violation": max(0.0, used - eps), "residuals": [float(x) for x in r]}


def full_report(model, sol: dict) -> dict:
    ops, spec = model.ops, model.spec
    c = sol["c"]
    out = {"unitarity": unitarity_report(ops, c),
           "f00_3": float(ops.f_proj["f00"] @ c),
           "f11_3": float(ops.f_proj["f11"] @ c),
           "c_norm_inf": float(np.abs(c).max()),
           "c_norm_1": float(np.abs(c).sum()),
           "rho_l4": float(np.sum(np.abs(np.concatenate(
               [c[ops.lay.r1], c[ops.lay.r2]])) ** 4) ** 0.25),
           "rho_l2": float(np.linalg.norm(np.concatenate(
               [c[ops.lay.r1], c[ops.lay.r2]])))}
    if spec.B is not None:
        out["B"] = spec.B
        out["B_norm"] = spec.B_norm
        out["B_active"] = bool(out["rho_l4" if spec.B_norm == "l4" else "rho_l2"]
                               >= 0.9 * spec.B)
    if spec.chiral:
        out["chiral"] = chiral_report(ops, c, spec.chi_caliber, spec.eps_chi)
    if spec.uv and "ImF" in sol:
        out["gram"] = gram_report(ops, c, sol["ImF"], sol["rho_hat"])
        out["fesr"] = fesr_report(ops, sol["rho_hat"], spec.sr_caliber)
        out["form_factor"] = ff_report(ops, sol["ImF"], spec.m_q, spec.eps_ff,
                                       spec.ff_frozen_at_s0)
    return out
