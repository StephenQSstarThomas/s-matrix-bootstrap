"""The 2403 amplitude at arbitrary real or complex s (PREREGISTRATION amendment 1; Phase B Task B2).

f^I_l(s) = sum_v coef_v(s) v over the export's amplitude variables.  The only difference from the node rows of
arbkernels is the s-channel kernel kappa_m(s) of H(s) = (1/pi) int g(x)[1/(x-s) - 1/(x-nu0)] dx ~ sum_m kappa_m g_m:
  'interp'   kappa_m = w_m(z(s)), the disk interpolant with Im w_m(z_k) = delta_mk, Re w_m(z_k) = K_km (mainline);
  'midpoint' kappa_m = (Delta/pi) x'_m [1/(x_m - s) - 1/(x_m - nu0)] (the authors' rule off the cut; control).
On the cut at a node the interpolant gives K + i delta, so the node rows are reproduced exactly.
Points: ("node", k) = x_k + i0; ("cut", s) = s + i0 for real s > 4; ("point", s) = complex s off the cut.
"""
from __future__ import annotations

from flint import acb, arb, ctx, fmpq

T0_COEF = {0: 5, 2: 2, 1: 0}        # 1/2 delta_l0 (3 + 1 + 1 | 1 + 1 | 1 - 1) for I = 0, 2, 1
_COMB = {0: lambda s, t, u: 3 * s + t + u, 2: lambda s, t, u: t + u, 1: lambda s, t, u: t - u}


def legendre_q(ell: int, z: acb) -> acb:
    """Q_l(z) for z off [-1, 1]: sqrt(pi) l! / (Gamma(l+3/2) (2z)^(l+1)) 2F1((l+1)/2, (l+2)/2; l+3/2; 1/z^2)."""
    pre = arb.pi().sqrt() * arb.fac_ui(ell) / arb(fmpq(2 * ell + 3, 2)).gamma()
    return pre / (2 * z) ** (ell + 1) * (1 / (z * z)).hypgeom_2f1(acb(fmpq(ell + 1, 2)), acb(fmpq(ell + 2, 2)),
                                                                   acb(fmpq(2 * ell + 3, 2)))


def _s_and_z(grid, where):
    kind, val = where
    a = arb(4 - grid.nu0).sqrt()
    if kind == "node":
        phi = arb.pi() * (2 * val + 1) / (2 * grid.M)
        return acb(grid.x[val]), acb(grid.xm4[val]), acb(phi.cos(), phi.sin()), val
    if kind == "cut":
        s = arb(val)
        b = (s - 4).sqrt()
        return acb(s), acb(s - 4), acb(a, b) / acb(a, -b), None
    s = acb(val)
    r = (4 - s).sqrt()
    return s, s - 4, (a - r) / (a + r), None


def _w(grid, j, z, node):
    """w_j(z); at z = z_j (node == j) the first term is replaced by its limit i."""
    M = grid.M
    phi = arb.pi() * (2 * j + 1) / (2 * M)
    zj = acb(phi.cos(), phi.sin())
    zi = 1 / zj
    t1 = acb(0, 1) if node == j else acb(0, 1) / (2 * M) * (1 - (z / zj) ** M) * (zj + z) / (zj - z)
    t2 = acb(0, 1) / (2 * M) * (1 - (z / zi) ** M) * (zi + z) / (zi - z)
    return t1 - t2


def kappa(grid, where, prescription="interp"):
    s, _, z, node = _s_and_z(grid, where)
    if prescription == "interp":
        return [_w(grid, j, z, node) for j in range(grid.M)]
    if prescription == "midpoint":
        if node is not None or where[0] == "cut":
            raise ValueError("the midpoint rule is defined off the cut only")
        return [acb(p) / (acb(x) - s) - acb(w0) for p, x, w0 in zip(grid.P, grid.x, grid.W0)]
    raise ValueError(prescription)


def f_row(grid, I: int, ell: int, where, prescription="interp"):
    """Complex coefficients of f^I_l at the point, in the export's variable order."""
    ctx.prec = grid.prec
    M, P, W0 = grid.M, [acb(p) for p in grid.P], [acb(w) for w in grid.W0]
    s, sm4, _, _ = _s_and_z(grid, where)
    kap = kappa(grid, where, prescription)
    eps = 1 if ell % 2 == 0 else -1
    half = acb(fmpq(1, 2)) if ell == 0 else acb(0)
    A = [p * legendre_q(ell, 1 + 2 * acb(x) / sm4) / sm4 for p, x in zip(P, grid.x)]
    Bt = [a - half * w for a, w in zip(A, W0)]
    Bu = [eps * a - half * w for a, w in zip(A, W0)]
    fs = [half * k for k in kap]
    Sv = [a + b for a, b in zip(Bt, Bu)]
    x = [acb(v) for v in grid.x]
    tu = [[(A[m] * P[n] + eps * A[n] * P[m]) / (x[m] + x[n] + sm4) - A[m] * W0[n] - eps * A[n] * W0[m]
           + half * W0[m] * W0[n] for n in range(M)] for m in range(M)]
    rng, pairs = range(M), [(a, b) for a in range(M) for b in range(a, M)]
    sig1 = (fs, Bt, Bu)
    sig2 = (Sv, [a + b for a, b in zip(fs, Bu)], [a + b for a, b in zip(fs, Bt)])
    r1 = ([kap[m] * Sv[n] for m in rng for n in rng],
          [kap[n] * Bt[m] + tu[m][n] for m in rng for n in rng],
          [tu[n][m] + kap[n] * Bu[m] for m in rng for n in rng])
    r2 = ([tu[a][b] + tu[b][a] if a != b else tu[a][a] for a, b in pairs],
          [kap[a] * Bu[b] + kap[b] * Bu[a] if a != b else kap[a] * Bu[a] for a, b in pairs],
          [kap[b] * Bt[a] + kap[a] * Bt[b] if a != b else kap[a] * Bt[a] for a, b in pairs])
    c = _COMB[I]
    t0 = half * T0_COEF[I]
    return [t0] + [c(*v) for grp in (sig1, sig2, r1, r2) for v in zip(*grp)]


def value(row, amp):
    """sum_v coef_v amp_v for a float or Arb amplitude vector."""
    tot = acb(0)
    for c, a in zip(row, amp):
        if a != 0:
            tot += c * (a if isinstance(a, (arb, acb)) else arb(float(a)))
    return tot
