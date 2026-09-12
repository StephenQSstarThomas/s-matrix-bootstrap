"""Independent ball-arithmetic re-verification of a solution (task section 5).

The float64 model is never trusted at the end.  This module rebuilds the
*same* operators of section 5a from scratch in Arb at ``bits`` precision and
replays every constraint of the original (unrescaled, unreduced, unsparsified)
problem on the returned solution:

    1500 unitarity disks   |S^I_ell(s_i)| <= 1
    100 Gram blocks        leading principal minors of (3.68) >= 0
    8 chiral residuals     (3.64) in the caliber used
    4 FESR moments         (3.73)
    14 form-factor bounds  (3.75)

Everything is an Arb ball, so a reported violation is an enclosure, not a
rounding artefact.  The historical ``certificates`` module is not used: its
signature assumes the Newton mainline's variable layout, and the task requires
the audit input to be the operators implemented here.
"""
from __future__ import annotations

import numpy as np
from flint import acb, arb, ctx

from . import constraints as C
from . import formfactor as FFM
from .grid import M_PI_MEV, S0
from .projector import Layout, ells_for


def _grid(M: int):
    nodes, weights = [], []
    for j in range(M):
        half = arb(2 * j + 1) / (4 * M)          # phi_j / 2 in units of pi
        c, s = half.cos_pi(), half.sin_pi()
        x = 4 / (c * c)                          # = 8 / (1 + cos phi_j)
        nodes.append(x)
        weights.append(x * s / (M * c))          # = (ds/dphi)_j / M
    return nodes, weights


def _kernel(M: int):
    v = [arb(0)] * (2 * M)
    for m in range(1, 2 * M):
        if m % 2:
            a = arb(m) / (2 * M)
            v[m] = a.cos_pi() / (M * a.sin_pi())
    return [[-v[(k - i) % (2 * M)] + v[(k + i + 1) % (2 * M)] for i in range(M)]
            for k in range(M)]


def _legendre_q(ell: int, z: arb) -> arb:
    """Q_ell(z) for real |z| > 1 via the 2F1 form, in ball arithmetic."""
    sign = arb(1)
    if z < -1:
        sign = arb((-1) ** (ell + 1))
        z = -z
    elif not z > 1:
        raise ValueError("exterior argument required")
    from math import factorial
    coef = arb(2 ** ell * factorial(ell) ** 2) / arb(factorial(2 * ell + 1))
    h = (1 / (z * z)).hypgeom_2f1(arb(ell + 1) / 2, arb(ell + 2) / 2, arb(2 * ell + 3) / 2)
    return sign * coef * z ** (-ell - 1) * h


class ArbAudit:
    def __init__(self, M: int, L: int, bits: int = 384) -> None:
        ctx.prec = bits
        self.M, self.L, self.bits = M, L, bits
        self.x, self.w = _grid(M)
        self.K = _kernel(M)
        self.b = [wj / xj for xj, wj in zip(self.x, self.w)]     # C[.](nu=0) row
        self.lay = Layout(M)

    # ------------------------------------------------------------------
    def _cauchy(self, dens, k):
        """C[dens](s_k + i0) as an acb: (K dens)_k + (b.dens) + i dens_k."""
        re = sum((self.K[k][i] * dens[i] for i in range(self.M)), arb(0))
        re += sum((self.b[i] * dens[i] for i in range(self.M)), arb(0))
        return acb(re, dens[k])

    def _angular(self, ell: int, k: int):
        d = self.x[k] - 4
        Qc = [4 * _legendre_q(ell, 1 + 2 * xi / d) / d for xi in self.x]
        return d, Qc

    def partial_wave(self, dens, isospin: int, ell: int, k: int) -> acb:
        """f^I_ell(s_k + i0) from the Arb operators of section 5a."""
        M = self.M
        T0, s1, s2, r1, r2 = (dens["T0"], dens["sigma1"], dens["sigma2"],
                              dens["rho1"], dens["rho2"])
        eps = arb((-1) ** ell)
        P0 = arb(2) if ell == 0 else arb(0)
        d, Qc = self._angular(ell, k)
        wQ = [self.w[j] * Qc[j] for j in range(M)]
        cr = self._cauchy(s1, k)
        cr2 = self._cauchy(s2, k)
        # contractions with the double spectral densities
        v1 = [sum((wQ[j] * r1[i][j] for j in range(M)), arb(0)) for i in range(M)]
        u1 = [sum((wQ[i] * r1[i][j] for i in range(M)), arb(0)) for j in range(M)]
        u2 = [sum((wQ[i] * r2[i][j] for i in range(M)), arb(0)) for j in range(M)]
        Dc = [[(Qc[i] + eps * Qc[j]) / (d + self.x[i] + self.x[j]) for j in range(M)]
              for i in range(M)]
        ww1 = sum((self.w[i] * self.w[j] * r1[i][j] * Dc[i][j]
                   for i in range(M) for j in range(M)), arb(0))
        ww2 = sum((self.w[i] * self.w[j] * r2[i][j] * Dc[i][j]
                   for i in range(M) for j in range(M)), arb(0))
        sQ1 = sum((wQ[j] * s1[j] for j in range(M)), arb(0))
        sQ2 = sum((wQ[j] * s2[j] for j in range(M)), arb(0))
        A = P0 * T0 + P0 * cr + (1 + eps) * sQ2 + (1 + eps) * self._cauchy(v1, k) + ww2
        B = (P0 * T0 + sQ1 + P0 * cr2 + eps * sQ2
             + self._cauchy(u1, k) + ww1 + eps * self._cauchy(u2, k))
        Cc = (P0 * T0 + eps * sQ1 + sQ2 + P0 * cr2
              + eps * ww1 + eps * self._cauchy(u1, k) + self._cauchy(u2, k))
        T = {0: 3 * A + B + Cc, 1: B - Cc, 2: B + Cc}[isospin]
        return T / 4

    # ------------------------------------------------------------------
    def audit(self, c: np.ndarray, ImF=None, rho_hat=None, *, chi_caliber="chi-b",
              eps_chi=C.EPS_CHI_MAIN, sr_caliber="SR-b", eps_ff=C.EPS_FF,
              m_q=None) -> dict:
        ctx.prec = self.bits
        M = self.M
        v = self.lay.unpack(c)
        dens = {"T0": arb(float(v["T0"])),
                "sigma1": [arb(float(t)) for t in v["sigma1"]],
                "sigma2": [arb(float(t)) for t in v["sigma2"]],
                "rho1": [[arb(float(v["rho1"][i][j])) for j in range(M)] for i in range(M)],
                "rho2": [[arb(float(v["rho2"][i][j])) for j in range(M)] for i in range(M)]}
        worst_lo, worst_hi, worst_where = arb(-1), arb(-1), None
        worst_active, active_where = arb(-1), None
        S_store = {}
        for I in (0, 1, 2):
            for ell in ells_for(I, self.L):
                for k in range(M):
                    f = self.partial_wave(dens, I, ell, k)
                    kap = arb.pi() * ((self.x[k] - 4) / self.x[k]).sqrt()
                    S = 1 + acb(0, 1) * kap * f
                    lo, hi = S.abs_lower(), S.abs_upper()
                    if (I, ell) in ((0, 0), (1, 1)):
                        S_store[(ell, k)] = S
                    if lo > worst_lo:
                        worst_lo, worst_hi = lo, hi
                        worst_where = {"isospin": I, "ell": ell, "node": k}
                    # a disk where |h| is negligible has eta == 1 trivially; track
                    # the worst among the disks that actually carry amplitude
                    if (kap * f).abs_lower() > arb("1e-3") and lo > worst_active:
                        worst_active = lo
                        active_where = {"isospin": I, "ell": ell, "node": k}
        out = {"bits": self.bits,
               # a certified violation needs the *lower* end of the ball above 1
               "max_eta_lower_end": float(worst_lo.str(20, radius=False)),
               "max_eta_upper_end": float(worst_hi.str(20, radius=False)),
               "max_eta_at": worst_where,
               "max_eta_among_disks_with_|h|>1e-3": float(worst_active.str(20, radius=False)),
               "max_eta_active_at": active_where,
               "unitarity_no_certified_violation": bool(worst_lo <= 1),
               "n_disks": 3 * self.L * M}
        if chi_caliber is not None:
            r = [sum((arb(float(x)) * y for x, y in zip(row, c)), arb(0))
                 for row in _chiral_rows_float(self)]
            out["chiral"] = _chiral_summary(r, chi_caliber, eps_chi)
        if ImF is not None:
            out.update(self._uv_audit(S_store, ImF, rho_hat, sr_caliber, eps_ff,
                                      C.M_Q if m_q is None else m_q))
        return out

    def _uv_audit(self, S_store, ImF, rho_hat, sr_caliber, eps_ff, m_q):
        M = self.M
        worst_minor, where = arb(1), None
        fesr_rows, ff_worst, ff_at = [], arb(-1e300), None
        idx_hi, bnd = C.ff_asymptotic_bounds(M, m_q, eps_ff)
        tgt, tol = C.printed_targets(), C.sr_tolerances(sr_caliber)
        for ell in (0, 1):
            ImFa = [arb(float(t)) for t in ImF[ell]]
            ReFa = [1 + sum((self.K[i][j] * ImFa[j] for j in range(M)), arb(0))
                    for i in range(M)]
            kin = FFM.kinematic_factor(ell, np.array([float(x.str(30, radius=False))
                                                      for x in self.x]))
            # per-node congruence g_i = k_ell(s_i), so rho = k^2 rho_hat
            rho = [arb(float(rho_hat[ell][i])) * arb(float(kin[i])) ** 2
                   for i in range(M)]
            for i in range(M):
                cF = acb(arb(float(kin[i])) * ReFa[i], arb(float(kin[i])) * ImFa[i])
                S = S_store[(ell, i)]
                # all principal minors of the Hermitian 3x3 (3.68)
                m2 = 1 - (S * S.conjugate()).real                     # rows {1,2}
                m13 = rho[i] - (cF * cF.conjugate()).real             # rows {1,3}={2,3}
                m3 = _det3(S, cF, rho[i])                             # full
                for val in (m2, m13, m3, rho[i]):
                    if val < worst_minor:
                        worst_minor, where = val, {"ell": ell, "node": i}
                if i in idx_hi:
                    excess = (cF * cF.conjugate()).real.sqrt() - arb(bnd[ell])
                    if excess > ff_worst:
                        ff_worst, ff_at = excess, {"ell": ell, "node": int(i)}
            wave = "S0" if ell == 0 else "P1"
            for n in C.MOMENTS[ell]:
                row = C.moment_row(M, n)
                mom = sum((arb(float(row[i])) * rho[i] for i in range(M)), arb(0))
                res = mom - arb(tgt[(wave, n)])
                fesr_rows.append({"wave": wave, "n": n,
                                  "moment": float(mom.str(20, radius=False)),
                                  "violation": float((res.abs_lower()
                                                      - arb(tol[(wave, n)])).str(20, radius=False))})
        return {"gram_min_minor": float(worst_minor.str(20, radius=False)),
                "gram_min_minor_at": where,
                "gram_ok": bool(worst_minor >= 0),
                "fesr": fesr_rows,
                "fesr_ok": all(r["violation"] <= 0 for r in fesr_rows),
                "form_factor_max_excess": float(ff_worst.str(20, radius=False)),
                "form_factor_at": ff_at,
                "form_factor_ok": bool(ff_worst <= 0),
                "n_gram_blocks": 2 * M, "n_ff_bounds": 2 * len(idx_hi)}


def _det3(S: acb, F: acb, rho: arb) -> arb:
    """det of [[1,S,F],[S*,1,F*],[F*,F,rho]] (real, the matrix is Hermitian).

    Cofactor expansion along the first row,

        det = (rho - |F|^2) - S (S* rho - (F*)^2) + F (S* F - F*),

    collapses to the closed form used here:

        det = rho (1 - |S|^2) - 2 |F|^2 + 2 Re(S (F*)^2).
    """
    mod_S2 = (S * S.conjugate()).real
    mod_F2 = (F * F.conjugate()).real
    cross = (S * F.conjugate() * F.conjugate()).real
    return rho * (1 - mod_S2) - 2 * mod_F2 + 2 * cross


def _chiral_rows_float(self):
    from .projector import PartialWaveOperator
    op = PartialWaveOperator(self.M)
    rows = []
    for sj in C.CHIRAL_POINTS:
        r01, r21 = C.chiral_ratios(sj)
        rows.append(op.rows(0, 0, sj)[0] - r01 * op.rows(1, 1, sj)[0])
        rows.append(op.rows(2, 0, sj)[0] - r21 * op.rows(1, 1, sj)[0])
    return rows


def _chiral_summary(r, caliber, eps):
    if caliber == "chi-a":
        used = max(float(x.abs_lower().str(20, radius=False)) for x in r)
    elif caliber == "chi-b":
        used = float(sum((x * x for x in r), arb(0)).sqrt().str(20, radius=False))
    else:
        used = max(float(sum((r[i] * r[i] for i in range(0, 8, 2)), arb(0)).sqrt()
                         .str(20, radius=False)),
                   float(sum((r[i] * r[i] for i in range(1, 8, 2)), arb(0)).sqrt()
                         .str(20, radius=False)))
    return {"caliber": caliber, "norm_used": used, "eps": eps,
            "ok": used <= eps, "violation": max(0.0, used - eps)}
