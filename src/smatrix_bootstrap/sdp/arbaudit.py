"""Independent Arb audit of the declared finite collocation model.

UV checks rebuild all seven principal minors per current Gram block, four raw
FESR moments and the high-node squared FF bounds from the paper. Each UV slack
has a certified pass/fail/inconclusive verdict and exact rational endpoints.
These are finite-node checks, not continuum or optimization certificates.
The older scattering summary and float chiral rows remain outside this A4
repair; their absence-of-violation flags are not rigorous feasibility claims.
"""
from __future__ import annotations

import numpy as np
from flint import acb, arb, ctx

from . import constraints as C
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
              m_q=None, ff_frozen_at_s0=True) -> dict:
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
                                      C.M_Q if m_q is None else m_q,
                                      frozen_at_s0=ff_frozen_at_s0))
        return out

    def _uv_audit(self, S_store, ImF, rho_hat, sr_caliber, eps_ff, m_q,
                  frozen_at_s0=True):
        """Replay the original UV inequalities without float operator exports."""
        ctx.prec = self.bits
        M, s0, pi = self.M, arb(3600) / 49, arb.pi()
        if sr_caliber not in ("SR-a", "SR-b", "SR-c"):
            raise ValueError(sr_caliber)
        idx_lo, idx_hi = [], []
        for i, x in enumerate(self.x):
            if x <= s0:
                idx_lo.append(i)
            elif x > s0:
                idx_hi.append(i)
            else:
                raise ValueError("Arb precision cannot resolve the hard cutoff")
        mq_default = m_q is None or m_q == C.M_Q
        mq = arb(113) / 2800 if mq_default else arb(str(m_q))
        eps = arb(str(eps_ff))
        gram, fesr, ff, gram_values, ff_values = [], [], [], [], []
        for ell in (0, 1):
            ImFa = [arb(t) if isinstance(t, arb) else arb(float(t)) for t in ImF[ell]]
            ReFa = [1 + sum((self.K[i][j] * ImFa[j] for j in range(M)), arb(0))
                    for i in range(M)]
            rhat = [arb(t) if isinstance(t, arb) else arb(float(t)) for t in rho_hat[ell]]
            k2 = [_kinematic_square(ell, x) for x in self.x]
            rho = [k2[i] * rhat[i] for i in range(M)]
            cap = 2 * mq * mq * eps if ell == 0 else eps / 2
            for i in range(M):
                F, S = acb(ReFa[i], ImFa[i]), S_store[(ell, i)]
                modF2 = ReFa[i] ** 2 + ImFa[i] ** 2
                m12 = 1 - S.real ** 2 - S.imag ** 2
                m13 = k2[i] * (rhat[i] - modF2)
                # Positive k² congruence gives original (unscaled) minors.
                minors = (arb(1), arb(1), rho[i], m12, m13, m13,
                          k2[i] * _det3(S, F, rhat[i]))
                for label, value in zip(("1", "2", "3", "12", "13", "23", "123"), minors):
                    gram.append(_constraint_row(value, ell=ell, node=i, minor=label))
                    gram_values.append(value)
                if i in idx_hi:
                    ffk2 = _kinematic_square(ell, s0) if frozen_at_s0 else k2[i]
                    used = ffk2 * modF2
                    ff.append(_constraint_row(cap - used, ell=ell, node=i,
                                               used_squared=_enclosure(used), cap=_enclosure(cap)))
                    ff_values.append(cap - used)
            wave = "S0" if ell == 0 else "P1"
            for n in ((0, 1) if ell == 0 else (-1, 0)):
                target = _printed_target(ell, n, s0)
                tol = (arb(".002") if sr_caliber == "SR-a" else
                       arb(".1" if sr_caliber == "SR-b" else ".2") * abs(target))
                mom = pi * sum((self.w[i] * self.x[i] ** n * rho[i] for i in idx_lo), arb(0))
                residual = mom - target
                slack = tol - abs(residual)
                fesr.append(_constraint_row(slack, wave=wave, n=n, moment=float(mom.mid()),
                            moment_interval=_enclosure(mom), target=_enclosure(target),
                            tolerance=_enclosure(tol), residual=_enclosure(residual),
                            violation=float((-slack).mid())))
        out = {"gram": gram, "fesr": fesr, "form_factor": ff,
               "n_gram_blocks": 2 * M, "n_gram_minors": 14 * M,
               "n_ff_bounds": len(ff), "bits": self.bits,
               "uv_scope": "Current constraints at native nodes; no scattering, chiral, or support certificate.",
               "uv_contract": {"cutoff": "3600/49", "cutoff_rule": "hard-midpoint",
                   "moment_units": "m_pi=1, raw integral rho(s) s^n ds",
                   "sr_caliber": sr_caliber, "ff_frozen_at_s0": bool(frozen_at_s0),
                   "eps_ff": _enclosure(eps), "m_q": _enclosure(mq),
                   "m_q_exact": "113/2800" if mq_default else str(m_q),
                   "parameters": "Exact printed decimals; nondefault parameters use their supplied decimal representation.",
                   "current_inputs": "Stored binary64 values are exact dyadics; supplied Arb balls retain uncertainty.",
                   "normalization": "F(0)=1; ReF=1+K ImF; rho=k(s)^2 rho_hat"}}
        for family, rows in (("gram", gram), ("fesr", fesr), ("form_factor", ff)):
            counts = {v: sum(r["verdict"] == v for r in rows) for v in
                      ("certified_pass", "certified_fail", "inconclusive")}
            verdict = ("certified_fail" if counts["certified_fail"] else
                       "inconclusive" if counts["inconclusive"] else "certified_pass")
            out.update({family + "_counts": counts, family + "_verdict": verdict,
                        family + "_ok": verdict == "certified_pass"})
        gi = min(range(len(gram)), key=lambda j: gram_values[j].lower())
        out.update(gram_min_minor=float(gram_values[gi].mid()), gram_min_minor_at=gram[gi])
        if ff:
            fi = min(range(len(ff)), key=lambda j: ff_values[j].lower())
            out.update(form_factor_max_excess=float((-ff_values[fi]).mid()), form_factor_at=ff[fi])
        out["form_factor_excess_kind"] = "squared-current excess (display only); verdict uses full slack ball"
        out["uv_ok"] = all(out[f + "_ok"] for f in ("gram", "fesr", "form_factor"))
        return out


def _kinematic_square(ell, s):
    beta = (1 - 4 / s).sqrt()
    return (3 * beta / (256 * arb.pi() ** 5) if ell == 0 else
            (s - 4) * beta / (384 * arb.pi() ** 5))


def _printed_target(ell, n, s0):
    if ell == 0:
        value = arb("3.09e-8") * (arb("27.38") / (n + 2) + (arb(".61") if n == 0 else 0))
    else:
        value = -arb("4.34e-6") * (-arb("13.26") / (n + 2) + (arb(".41") if n == 0 else 0))
    return value * s0 ** (n + 2)


def _enclosure(value):
    """Exact dyadic endpoint strings avoid both underflow and inward rounding."""
    return {"ball": value.str(40), "lower": str(value.lower().fmpq()),
            "upper": str(value.upper().fmpq())}


def _constraint_row(slack, **metadata):
    verdict = "certified_pass" if slack >= 0 else "certified_fail" if slack < 0 else "inconclusive"
    return {**metadata, "slack": _enclosure(slack), "verdict": verdict}


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
