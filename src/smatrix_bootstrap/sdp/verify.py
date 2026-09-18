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


def fesr_report(ops, rho_hat, caliber: str, eps_sr: float = C.EPS_SR) -> dict:
    M = ops.M
    tgt, tol = C.printed_targets(), C.sr_tolerances(caliber, eps_sr)
    rows = []
    for ell, wave in ((0, "S0"), (1, "P1")):
        rho = rho_hat[ell] * FFM.gram_scale(ell, ops.s) ** 2
        wave_rows = []
        for n in C.MOMENTS[ell]:
            m = float(C.moment_row(M, n) @ rho)
            wave_rows.append({"wave": wave, "n": n, "moment": m, "target": tgt[(wave, n)],
                              "residual": m - tgt[(wave, n)], "tolerance": tol[(wave, n)],
                              "violation": max(0.0, abs(m - tgt[(wave, n)]) - tol[(wave, n)])})
        if caliber == "SR-d":
            # per-wave L2 packaging: one ball of radius eps_SR over both moments
            norm = float(np.sqrt(sum(r["residual"] ** 2 for r in wave_rows)))
            for r in wave_rows:
                r["wave_l2_residual"] = norm
                r["violation"] = max(0.0, norm - eps_sr)
        rows += wave_rows
    return {"caliber": caliber, "rows": rows,
            "packaging": "per-wave L2 ball" if caliber == "SR-d" else "per-moment box",
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
    rho = np.concatenate([c[ops.lay.r1], c[ops.lay.r2]])
    out["rho_linf"] = float(np.abs(rho).max())
    if spec.reg_norm is not None:
        used = out["rho_" + spec.reg_norm]
        out["regulariser"] = {"norm": spec.reg_norm, "bound": spec.reg_bound, "norm_used": used,
                              "rho_linf": out["rho_linf"],
                              "active_fraction": (float(np.mean(np.abs(rho) >= 0.99 * spec.reg_bound))
                                                  if spec.reg_norm == "linf" else None),
                              "active": bool(used >= 0.99 * spec.reg_bound)}
    if spec.B is not None:
        out["B"] = spec.B
        out["B_norm"] = spec.B_norm
        out["B_active"] = bool(out["rho_l4" if spec.B_norm == "l4" else "rho_l2"]
                               >= 0.9 * spec.B)
    if spec.chiral:
        out["chiral"] = chiral_report(ops, c, spec.chi_caliber, spec.eps_chi)
    if spec.uv and "ImF" in sol:
        out["gram"] = gram_report(ops, c, sol["ImF"], sol["rho_hat"])
        out["fesr"] = fesr_report(ops, sol["rho_hat"], spec.sr_caliber, spec.eps_sr)
        out["form_factor"] = ff_report(ops, sol["ImF"], spec.m_q, spec.eps_ff,
                                       spec.ff_frozen_at_s0)
    return out


def solution_report(model, sol, tolerance=1e-8):
    """Numerical original-coordinate feasibility; never an optimality certificate.

    Residuals for small chiral, moment and FF constraints are divided by their
    own tolerance/cap before acceptance. The Gram is also checked after its
    invertible charge-factor congruence to avoid hiding tiny spectral errors.
    """
    spec, ops, c = model.spec, model.ops, sol["c"]
    finite = all(np.all(np.isfinite(sol.get(key, np.nan))) for key in
                 (("c", "ImF", "rho_hat") if spec.uv else ("c",)))
    out = {"certified": False, "primal_feasible": False,
           "scope": "numerical native-node feasibility; no rigorous dual certificate",
           "constraint_checks": {"finite": bool(finite)}}
    if not finite:
        return out
    out.update(full_report(model, sol))
    checks = out["constraint_checks"]
    checks["unitarity"] = out["unitarity"]["feasible"]
    if spec.chiral:
        checks["chiral"] = out["chiral"]["violation"] <= tolerance * spec.eps_chi
    if spec.B is not None:
        checks["density"] = out["rho_" + spec.B_norm] <= spec.B * (1 + tolerance)
    if spec.reg_norm is not None:
        checks["regulariser"] = out["rho_" + spec.reg_norm] <= spec.reg_bound * (1 + tolerance)
    if spec.uv:
        checks["rho_nonnegative"] = bool(np.min(sol["rho_hat"]) >= -tolerance)
        if "gram" in spec.uv_parts:
            worst = np.inf
            K = FFM.hilbert_kernel(spec.M)
            for ell in (0, 1):
                re, im = ops.gram_rows[ell]
                S = 1 - im @ c + 1j * (re @ c)
                F = 1 + K @ sol["ImF"][ell] + 1j * sol["ImF"][ell]
                for i in range(spec.M):
                    G = FFM.gram_block(S[i], F[i], sol["rho_hat"][ell, i])
                    scale = np.sqrt(np.maximum(np.abs(np.diag(G)), 1.))
                    worst = min(worst, float(np.linalg.eigvalsh(
                        G / scale[:, None] / scale[None, :]).min()))
            out["gram"]["min_equilibrated_eigenvalue"] = worst
            checks["gram"] = worst >= -tolerance
        if "fesr" in spec.uv_parts:
            free = {(str(w), int(n)) for w, n in (spec.sr_free or ())}
            for row in out["fesr"]["rows"]:
                row["free"] = (row["wave"], row["n"]) in free
            imposed = [row for row in out["fesr"]["rows"] if not row["free"]]
            v = max((row["violation"] / row["tolerance"] for row in imposed), default=0.0)
            out["fesr"]["max_violation_over_tolerance"] = v
            checks["fesr"] = v <= tolerance
        if "ff" in spec.uv_parts:
            idx, cap, kin = C.ff_asymptotic_bounds(spec.M, spec.m_q, spec.eps_ff,
                                                  spec.ff_frozen_at_s0, spec.eps_ff_s0, spec.eps_ff_p1)
            K = FFM.hilbert_kernel(spec.M)
            v = max(float(np.max(np.abs(kin[e] * (1 + K @ sol["ImF"][e]
                    + 1j * sol["ImF"][e])[idx]) / cap[e] - 1)) for e in (0, 1))
            out["form_factor"]["max_relative_violation"] = v
            checks["ff"] = v <= tolerance
    if model.fix_f00 is not None:
        residual = abs(out["f00_3"] - model.fix_f00)
        out["section_residual"] = residual
        checks["section"] = residual <= tolerance * max(1, abs(model.fix_f00))
    if getattr(model, "functional", None) is not None:
        f = model.functional
        if f["kind"] == "watson":
            from .observables import _wave_h
            from .watson import weight_of
            raw = 0.0
            for wave in f["waves"]:
                h = _wave_h(ops, c, wave)
                for i, ((tr, ti), k) in enumerate(zip(f["targets"][wave], f["nodes"])):
                    raw += weight_of(f, wave, i) * (float(tr) * float(h[k].real) + (float(ti) - 1.0) * float(h[k].imag))
        elif f["kind"] in ("ImKH", "ImS"):
            from .observables import _wave_h
            h = _wave_h(ops, c, f["wave"])[f["node"]]
            raw = float(h.imag) if f["kind"] == "ImKH" else float(h.real)
        elif f["kind"] == "SRmom":
            raw = float(C.moment_row(spec.M, f["node"]) * FFM.gram_scale(f["ell"], ops.s) ** 2
                        @ sol["rho_hat"][f["ell"]])
        else:
            raw = float(sol["ImF" if f["kind"] == "ImF" else "rho_hat"][f["ell"]][f["node"]])
        out["functional"] = dict(f, value=raw)
        out["objective_recomputed"] = raw if f["sense"] == "max" else -raw
    else:
        out["objective_recomputed"] = float(np.dot(model.direction,
                                                   [out["f00_3"], out["f11_3"]]))
    if getattr(model, "face", None) is not None:
        d0, d1 = model.face["direction"]; floor = model.face["value"] - model.face["margin"]
        here = d0 * out["f00_3"] + d1 * out["f11_3"]
        out["face"] = dict(model.face, floor=floor, value_here=here, slack=here - floor)
        checks["face"] = here >= floor - tolerance * max(1.0, abs(floor))
    out["primal_feasible"] = bool(all(checks.values()))
    out["verification_tolerance"] = tolerance
    return out
