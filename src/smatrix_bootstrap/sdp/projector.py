"""Partial-wave projection of the Mandelstam representation -- (2.5)-(2.10).

Everything here is re-derived from the paper; no repository operator is used.

Amplitude (2.7), with ``C[f](nu) = (1/pi) int_4^inf f(x)/(x-nu) dx`` and
``D[g](nu,lam) = (1/pi^2) int int g(x,y) / ((x-nu)(y-lam))``::

    A(s,t,u) = T0 + C[s1](s) + C[s2](t) + C[s2](u)
                  + D[r1](s,t) + D[r1](s,u) + D[r2](t,u)

with ``r2`` symmetric, so ``A(s,t,u) = A(s,u,t)`` manifestly.

Angular integrals.  With ``a = (s-4)/2``, ``t(mu) = -a(1-mu)``, ``u(mu) = t(-mu)``
and ``zeta(x) = 1 + x/a``, Neumann's formula gives the closed forms

    Qc_ell(x)      := int_{-1}^{1} P_ell(mu) / (x - t(mu)) dmu = (2/a) Q_ell(zeta(x))
    int P_ell(mu) / (x - u(mu)) dmu                            = (-1)^ell Qc_ell(x)
    Dc_ell(x,y)    := int P_ell(mu) / ((x-t)(y-u)) dmu
                    = [ Qc_ell(x) + (-1)^ell Qc_ell(y) ] / (s - 4 + x + y)

the last by the partial fraction ``1/((zx-mu)(zy+mu)) =
[1/(zx-mu) + 1/(zy+mu)] / (zx+zy)``.  The denominator ``s-4+x+y`` never vanishes
on the grid: for ``s>4`` all terms are positive, and for ``0<s<4`` we have
``x+y >= 8 > 4-s``.

Isospin (2.5) and the projection (2.9)::

    T^0 = 3 A(s,t,u) + A(t,s,u) + A(u,t,s),
    T^1 = A(t,s,u) - A(u,t,s),
    T^2 = A(t,s,u) + A(u,t,s),      f^I_ell(s) = (1/4) int_{-1}^{1} P_ell T^I dmu.

Every term of ``A(s,t,u)`` and of ``A(t,s,u) + A(u,t,s)`` carries a factor
``(1 + (-1)^ell)`` or ``delta_{ell,0}``, and every term of
``A(t,s,u) - A(u,t,s)`` a factor ``(1 - (-1)^ell)``.  The even/odd selection
rule of (2.9) is therefore an identity of the assembled operator, which the
self-checks verify numerically rather than assume.
"""
from __future__ import annotations

import numpy as np

from .grid import dsdphi, n_amplitude_vars, s_nodes, sym_pack_index
from .hilbert import cauchy_offcut_row, cauchy_oncut_real
from .legendreq import QCache

def ells_for(isospin: int, L: int) -> tuple[int, ...]:
    """First ``L`` allowed angular momenta for the given isospin (2.9)."""
    if isospin not in (0, 1, 2) or L < 1:
        raise ValueError("isospin must be 0,1,2 and L positive")
    return tuple(range(isospin % 2, 2 * L, 2))


class Layout:
    """Packing of the amplitude variables (3.62): T0, sigma_{a,i}, rho_{a,ij}."""

    def __init__(self, M: int) -> None:
        self.M = M
        self.n = n_amplitude_vars(M)
        self.i_T0 = 0
        self.s1 = slice(1, 1 + M)
        self.s2 = slice(1 + M, 1 + 2 * M)
        self.r1 = slice(1 + 2 * M, 1 + 2 * M + M * M)
        self.r2 = slice(1 + 2 * M + M * M, self.n)
        self.triu = sym_pack_index(M)

    def fold_sym(self, C: np.ndarray) -> np.ndarray:
        """Fold an ``M x M`` coefficient matrix onto the packed symmetric rho_2."""
        S = C + C.T
        np.fill_diagonal(S, np.diag(C))
        return S[self.triu]

    def unpack(self, c: np.ndarray) -> dict:
        M = self.M
        r2 = np.zeros((M, M))
        r2[self.triu] = c[self.r2]
        r2 = r2 + r2.T - np.diag(np.diag(r2))
        return {"T0": c[self.i_T0], "sigma1": c[self.s1], "sigma2": c[self.s2],
                "rho1": c[self.r1].reshape(M, M), "rho2": r2}


class PartialWaveOperator:
    """Builds the real matrix rows that give ``f^I_ell(s)`` from the variables."""

    def __init__(self, M: int, cache: QCache | None = None) -> None:
        self.M = M
        self.lay = Layout(M)
        self.s = s_nodes(M)
        self.omega = dsdphi(M) / M          # (1/pi) int dx  ->  sum_j omega_j
        self.A_on = cauchy_oncut_real(M)    # Re C[f](s_k + i0)
        self.qc = cache or QCache()

    # ------------------------------------------------------------------
    def _cauchy_rows(self, s: float, node: int | None):
        """(re, im) length-M rows with ``C[f](s) = re.f + i im.f``."""
        if node is None:
            return cauchy_offcut_row(self.M, s), np.zeros(self.M)
        im = np.zeros(self.M)
        im[node] = 1.0
        return self.A_on[node].copy(), im

    def _angular(self, s: float, ell: int):
        """``Qc_ell`` on the grid and ``Dc_ell`` on the grid squared."""
        a = 0.5 * (s - 4.0)
        zeta = 1.0 + self.s / a
        q = self.qc.get(zeta, max(19, ell))[ell]
        Qc = (2.0 / a) * q
        sgn = 1.0 if ell % 2 == 0 else -1.0
        denom = (s - 4.0) + self.s[:, None] + self.s[None, :]
        Dc = (Qc[:, None] + sgn * Qc[None, :]) / denom
        return Qc, Dc, sgn

    def _abc(self, s: float, ell: int, node: int | None):
        """Coefficient vectors of ``int P_ell X dmu`` for X = A(s,t,u), A(t,s,u), A(u,t,s).

        Returns a ``(3, 2, n)`` array indexed by [which, (re,im), variable].
        """
        M, lay, w = self.M, self.lay, self.omega
        Qc, Dc, sgn = self._angular(s, ell)
        P0 = 2.0 if ell == 0 else 0.0
        cr = np.stack(self._cauchy_rows(s, node))          # (2, M)
        wQ = w * Qc                                         # omega_j Qc_j
        wwD = (w[:, None] * w[None, :]) * Dc                # omega_i omega_j Dc_ij

        out = np.zeros((3, 2, lay.n))
        for p in range(2):                                  # p = 0 (Re), 1 (Im)
            c = cr[p]
            # ---- A(s,t,u) = T0 + C[s1](s) + C[s2](t)+C[s2](u) + D[r1](s,t)+D[r1](s,u) + D[r2](t,u)
            out[0, p, lay.i_T0] += P0 if p == 0 else 0.0
            out[0, p, lay.s1] += P0 * c
            out[0, p, lay.s2] += (1.0 + sgn) * wQ if p == 0 else 0.0
            out[0, p, lay.r1] += ((1.0 + sgn) * np.outer(c, wQ)).ravel()
            out[0, p, lay.r2] += lay.fold_sym(wwD) if p == 0 else 0.0
            # ---- A(t,s,u) = T0 + C[s1](t) + C[s2](s) + C[s2](u) + D[r1](t,s)+D[r1](t,u)+D[r2](s,u)
            out[1, p, lay.i_T0] += P0 if p == 0 else 0.0
            out[1, p, lay.s1] += wQ if p == 0 else 0.0
            out[1, p, lay.s2] += P0 * c + (sgn * wQ if p == 0 else 0.0)
            out[1, p, lay.r1] += (np.outer(wQ, c) + (wwD if p == 0 else 0.0)).ravel()
            out[1, p, lay.r2] += lay.fold_sym(sgn * np.outer(c, wQ))
            # ---- A(u,t,s) = T0 + C[s1](u) + C[s2](t) + C[s2](s) + D[r1](u,t)+D[r1](u,s)+D[r2](t,s)
            out[2, p, lay.i_T0] += P0 if p == 0 else 0.0
            out[2, p, lay.s1] += sgn * wQ if p == 0 else 0.0
            out[2, p, lay.s2] += P0 * c + (wQ if p == 0 else 0.0)
            out[2, p, lay.r1] += (sgn * (wwD if p == 0 else 0.0) + sgn * np.outer(wQ, c)).ravel()
            out[2, p, lay.r2] += lay.fold_sym(np.outer(wQ, c))
        return out

    # ------------------------------------------------------------------
    def rows(self, isospin: int, ell: int, s: float, node: int | None = None) -> np.ndarray:
        """``(2, n)`` array: real and imaginary coefficient rows of ``f^I_ell(s)``."""
        abc = self._abc(s, ell, node)
        if isospin == 0:
            T = 3.0 * abc[0] + abc[1] + abc[2]
        elif isospin == 1:
            T = abc[1] - abc[2]
        elif isospin == 2:
            T = abc[1] + abc[2]
        else:
            raise ValueError(isospin)
        return 0.25 * T

    def lambda_row(self) -> np.ndarray:
        """Coefficients of the quartic coupling at the crossing-symmetric point.

        Below (2.14) the paper defines

            lambda = (pi/4) T_{33,33}(4/3, 4/3, 4/3),

        and T_{33,33} = A(s,t,u) + A(t,s,u) + A(u,t,s) collapses to
        3 A(4/3,4/3,4/3) there.  The Weinberg amplitude (2.12) gives
        lambda = m_pi^2 / (32 pi f_pi^2) = 0.023.  The quoted 2.661 is the
        single-scalar benchmark of Paulos et al., not a proved attainable
        maximum for the full O(3) problem.  The neutral S=(S0+2*S2)/3 lies in
        the unit disk whenever both isospin amplitudes do; the converse need
        not hold.  This gives an upper-bound comparison, not an equality gate.

        A(nu1,nu2,nu3) is bilinear in the densities, so the row is assembled
        from the same Cauchy weights used everywhere else.
        """
        M, lay, w = self.M, self.lay, self.omega
        nu = 4.0 / 3.0
        cw = w / (self.s - nu)                      # C[f](nu) = cw . f
        row = np.zeros(lay.n)
        row[lay.i_T0] = 1.0
        row[lay.s1] = cw
        row[lay.s2] = 2.0 * cw
        row[lay.r1] = (2.0 * np.outer(cw, cw)).ravel()
        row[lay.r2] = lay.fold_sym(np.outer(cw, cw))
        return 0.75 * np.pi * row                   # lambda = (3 pi / 4) A

    def parity_residual(self, ell: int, s: float, node: int | None = None) -> float:
        """Size of the projection that the selection rule of (2.9) forbids,
        relative to the size that survives.  Must be ~1e-16."""
        abc = self._abc(s, ell, node)
        keep = abc[1] + abc[2] if ell % 2 == 0 else abc[1] - abc[2]
        drop = abc[1] - abc[2] if ell % 2 == 0 else abc[1] + abc[2]
        drop_A = np.abs(abc[0]).max() if ell % 2 else 0.0
        scale = max(np.abs(keep).max(), np.abs(abc[0]).max(), 1e-300)
        return float(max(np.abs(drop).max(), drop_A) / scale)
