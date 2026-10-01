"""Independent Arb rebuild of the 2403 finite-model kernels (verification levels V0/V1).

Derived from the subtracted Mandelstam representation (units m_pi = 1)

    A(s,t,u) = T0 + H1(s) + H2(t) + H2(u) + Ht1(s,t) + Ht1(s,u) + Ht2(t,u),
    H(nu) = (1/pi) int sigma(x) [1/(x-nu) - 1/(x-nu0)],

projected as f_l = 1/4 int dmu P_l(mu) T^I with 1/4 int P_l/(x-t) = Q_l(1+2x/(s-4))/(s-4)
(u gives (-1)^l) and, for the (t,u) term, 1/((x-t)(y-u)) = [1/(x-t)+1/(y-u)]/(x+y+s-4).
Integrals are midpoint sums in phi (Jacobian x' absorbed), boundary values on the cut
use the discrete conjugate kernel K.  Every number is a ball at ``prec_bits``.
Row/variable order and file names follow the export format of ``kernels.py``.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial, prod

from flint import arb, ctx, fmpq

NU0 = -20
S0 = fmpq(200000, 13957) ** 2                     # (2000 MeV / 139.57 MeV)^2
ALPHA_S, M_U, M_D = fmpq(314588744, 10 ** 9), fmpq(360, 13957), fmpq(650, 13957)
CHI_POINTS = (fmpq(1, 2), fmpq(1), fmpq(3, 2), fmpq(2))
F_POINT = fmpq(3)
N_WAVES = 10                                      # ell <= 19 per isospin
ISO_BLOCKS = (0, 2, 1)                            # row blocks [I=0 | I=2 | I=1]
SLOTS = {0: (3, 1, 1), 2: (0, 1, 1), 1: (0, 1, -1)}   # T^I from A(s,t,u), A(t,s,u), A(u,t,s)
N_CELLS = 37                                      # uniform B-spline cells on [0, x_n0]
Z = arb(0)


def _rat(v) -> fmpq:
    """s* given as fmpq, int, Fraction or text such as "3/2" or "1.5"."""
    f = v if isinstance(v, (int, Fraction)) else Fraction(str(v))
    return v if isinstance(v, fmpq) else fmpq(Fraction(f).numerator, Fraction(f).denominator)


def ells(I: int) -> list[int]:
    return list(range(1, 2 * N_WAVES, 2)) if I == 1 else list(range(0, 2 * N_WAVES, 2))


class Grid:
    """Nodes x_j = nu0 + (8-2nu0)/(1+cos phi_j), phi_j = (j-1/2)pi/M, and node kernels."""

    def __init__(self, M: int, nu0: int = NU0, prec_bits: int = 256):
        ctx.prec = prec_bits
        self.M, self.nu0, self.prec = M, nu0, prec_bits
        a = 4 - nu0                               # x = nu0 + a/cos^2(phi/2)
        self.x, self.xp, self.xm4, self.xmn0 = [], [], [], []
        for j in range(1, M + 1):
            h = arb(fmpq(2 * j - 1, 4 * M))       # phi_j / (2 pi)
            c, s = h.cos_pi(), h.sin_pi()
            self.xmn0.append(a / (c * c))         # x - nu0
            self.x.append(nu0 + a / (c * c))
            self.xm4.append(a * (s / c) ** 2)     # x - 4, without cancellation
            self.xp.append(a * s / c ** 3)        # dx/dphi = (8-2nu0) sin phi/(1+cos phi)^2
        pi = arb.pi()
        self.delta = pi / M
        self.omega = [pi * (xm4 / x).sqrt() for xm4, x in zip(self.xm4, self.x)]
        self.P = [self.delta * xp / pi for xp in self.xp]              # (Delta/pi) x'_m
        self.W0 = [p / d for p, d in zip(self.P, self.xmn0)]           # (Delta/pi) x'/(x-nu0)
        self.K0 = [w * nu0 / x for w, x in zip(self.W0, self.x)]       # F(0) = 1 shift
        below = [bool(x < arb(S0)) for x in self.x]
        assert all(below[j] or bool(self.x[j] > arb(S0)) for j in range(M))
        self.n0 = sum(below)
        c = [Z] * (2 * M)
        for n in range(1, 2 * M, 2):
            c[n] = arb(fmpq(n, 2 * M)).cot_pi() / M
        self.K = [[c[(i + j + 1) % (2 * M)] - c[(i - j) % (2 * M)] for j in range(M)]
                  for i in range(M)]
        self.nvars = 1 + 2 * M + M * M + M * (M + 1) // 2


def legendre_q(ell: int, z: arb) -> arb:
    """Q_l(z) = sqrt(pi) l!/(Gamma(l+3/2)(2z)^(l+1)) 2F1((l+1)/2,(l+2)/2;l+3/2;1/z^2), |z| > 1."""
    if z < -1:
        return (-1) ** (ell + 1) * legendre_q(ell, -z)
    if not z > 1:
        raise ValueError("Q_l needs a real argument with |z| > 1")
    coef = fmpq(factorial(ell), prod(range(2 * ell + 1, 0, -2)))     # sqrt(pi) l!/(Gamma(l+3/2) 2^(l+1))
    h = (1 / (z * z)).hypgeom_2f1(arb(ell + 1) / 2, arb(ell + 2) / 2, arb(2 * ell + 3) / 2,
                                  abc=True, ac=ell % 2 == 0, bc=ell % 2 == 1)
    return arb(coef) * h / z ** (ell + 1)


def lam(ell: int, s) -> arb:
    """Lambda_l(s) = ((sqrt s - 2)/(sqrt s + 2))^(l/2)."""
    s = s if isinstance(s, arb) else arb(_rat(s))
    r = (s.sqrt() - 2) / (s.sqrt() + 2)
    return r.sqrt() ** ell


def svz() -> list[arb]:
    """Identity-operator FESR targets {S0(0,1,2), P1(-1,0,1), D0(-2,-1,0)} from exact inputs."""
    pi = arb.pi()
    a, base = arb(ALPHA_S) / pi, 1 / (2 * pi) ** 4
    s0 = [base * (M_U ** 2 + M_D ** 2) * 3 / (4 * pi * (n + 2)) * (1 + fmpq(13, 3) * a) for n in (0, 1, 2)]
    p1 = [base / (4 * pi * (n + 2)) * (1 + a) for n in (-1, 0, 1)]
    d0 = [base / (4 * pi * (n + 3)) * (fmpq(11, 10) - fmpq(17, 18) * a) for n in (-2, -1, 0)]
    return s0 + p1 + d0


def bspline_u(u: arb) -> list[arb]:
    """The 40 clamped cubic B-splines with knots {0,0,0,0,1,...,36,37,37,37,37} at u (Cox-de Boor)."""
    T = [0, 0, 0] + list(range(N_CELLS + 1)) + [N_CELLS] * 3
    if u < 0 or u > N_CELLS:
        raise ValueError("u outside [0, 37]")
    if u.is_exact() and u == N_CELLS:
        span = 2 + N_CELLS                         # last non-empty span, closed at the right
    else:
        k = int(float(u.mid()))
        k = k if arb(k) <= u else k - 1
        k = k if u < arb(k + 1) else k + 1
        if not (arb(k) <= u and u < arb(k + 1)):
            raise ValueError("knot span is not decided at this precision")
        span = 3 + k
    N = [arb(1) if i == span else Z for i in range(len(T) - 1)]
    for p in (1, 2, 3):
        new = []
        for i in range(len(T) - 1 - p):
            v = Z
            if T[i + p] != T[i] and not N[i].is_zero():
                v += (u - T[i]) / (T[i + p] - T[i]) * N[i]
            if T[i + p + 1] != T[i + 1] and not N[i + 1].is_zero():
                v += (T[i + p + 1] - u) / (T[i + p + 1] - T[i + 1]) * N[i + 1]
            new.append(v)
        N = new
    return N


def bsplines(grid: Grid) -> list[list[arb]]:
    """[B_i(x_j)] for j <= n0 (rows) and i = 1..40 (columns), knots scaled to [0, x_n0]."""
    ctx.prec = grid.prec
    top = grid.n0 - 1                                  # x_n0 is the last knot exactly
    return [bspline_u(arb(N_CELLS) if j == top else N_CELLS * grid.x[j] / grid.x[top])
            for j in range(grid.n0)]


def kappa(ell: int, x: arb, xm4: arb) -> arb:
    """rho = kappa_l(s) rho_hat for the S0, P1, D0 currents."""
    coef = (fmpq(3, 256), fmpq(1, 384), fmpq(1, 2560))[ell]      # times (s-4)^(l+1/2)/(pi^5 sqrt s)
    return coef * xm4 ** ell * xm4.sqrt() / (arb.pi() ** 5 * x.sqrt())


def moments(grid: Grid) -> dict[str, list[list[arb]]]:
    """FESR moment rows kappa_l(x_j) Delta x'_j x_j^n / x_n0^(n+p), j <= n0 (midpoint in phi)."""
    ctx.prec = grid.prec
    xn = grid.x[grid.n0 - 1]
    table = {"intS0": (0, (0, 1, 2), 2), "intP1": (1, (-1, 0, 1), 2), "intD0": (2, (-2, -1, 0), 3)}
    out = {}
    for name, (ell, ns, off) in table.items():
        base = [kappa(ell, grid.x[j], grid.xm4[j]) * grid.delta * grid.xp[j] for j in range(grid.n0)]
        out[name] = [[b * grid.x[j] ** n / xn ** (n + off) for j, b in enumerate(base)] for n in ns]
    return out


def ff_kernel(grid: Grid) -> list[list[arb]]:
    """Re F_j = 1 + sum_m (K_jm + K0_m) Im F_m (F(0) = 1 with the disk centred at nu0)."""
    return [[k + k0 for k, k0 in zip(row, grid.K0)] for row in grid.K]


# --------------------------------------------------------------------------- partial-wave rows
_COMB = {0: lambda s, t, u: 3 * s + t + u, 2: lambda s, t, u: t + u, 1: lambda s, t, u: t - u}


def _point(grid: Grid, kind: str, where):
    """(node j or None, s - 4, s-channel real-part row, h/f scale) for a node or an s* < 4."""
    if kind == "phys":
        j = int(where)
        return j, grid.xm4[j], grid.K[j], grid.omega[j], ("phys", j)
    if kind != "unphys":
        raise ValueError(f"kind must be 'phys' or 'unphys', got {kind!r}")
    s = _rat(where)
    if not 0 < s < 4:
        raise ValueError("unphysical points need 0 < s* < 4")
    sa = arb(s)
    Kr = [p / (x - sa) - w0 for p, x, w0 in zip(grid.P, grid.x, grid.W0)]   # K^u(s*, m)
    return None, sa - 4, Kr, None, ("unphys", s)


def _denominators(grid: Grid, key, sm4) -> list[list[arb]]:
    """1/(x_m + x_n + s - 4), the (t,u) partial-fraction denominator, cached per point."""
    cache = grid.__dict__.setdefault("_dcache", {})
    if key not in cache:
        M, x = grid.M, grid.x
        D = [[Z] * M for _ in range(M)]
        for m in range(M):
            for n in range(m, M):
                D[m][n] = D[n][m] = 1 / (x[m] + x[n] + sm4)
        cache[key] = D
    return cache[key]


def rows_at(grid: Grid, ell: int, kind: str, where, isospins=(0, 1, 2)) -> dict:
    """{I: (re, im)} over the export's variables at one point; h = omega f at a node, f at s*."""
    ctx.prec = grid.prec
    M, P, W0 = grid.M, grid.P, grid.W0
    j, sm4, Kr, scale, key = _point(grid, kind, where)
    eps = 1 if ell % 2 == 0 else -1
    half = arb(fmpq(1, 2)) if ell == 0 else None
    # single functions: A_l(s,x_m) = Delta x'_m Q_l(1 + 2x_m/(s-4)) / (pi (s-4))
    A = [p * legendre_q(ell, 1 + 2 * x / sm4) / sm4 for p, x in zip(P, grid.x)]
    hw = [half * w for w in W0] if half else [Z] * M
    Bt = [a - h for a, h in zip(A, hw)]                      # f^(t): A - 1/2 d_l0 W0
    Bu = [eps * a - h for a, h in zip(A, hw)]                # f^(u)
    fs = [half * k for k in Kr] if half else [Z] * M         # Re f^(s)
    fsi = [half if (half and m == j) else Z for m in range(M)]   # Im f^(s)
    Sv = [a + b for a, b in zip(Bt, Bu)]
    # (t,u): A_m W(s;m,n) + (-1)^l A_n W(s;n,m) + 1/2 d_l0 W00, W(s;a,b) = P_b D_ab - W0_b
    D = _denominators(grid, key, sm4)
    tu = [[D[m][n] * (A[m] * P[n] + eps * A[n] * P[m]) - A[m] * W0[n] - eps * A[n] * W0[m]
           + (half * W0[m] * W0[n] if half else Z) for n in range(M)] for m in range(M)]
    rng, pairs = range(M), [(a, b) for a in range(M) for b in range(a, M)]
    # slot(s), slot(t), slot(u) = contributions of A(s,t,u), A(t,s,u), A(u,t,s)
    sig1 = (fs, Bt, Bu)
    sig2 = (Sv, [a + b for a, b in zip(fs, Bu)], [a + b for a, b in zip(fs, Bt)])
    r1 = ([Kr[m] * Sv[n] for m in rng for n in rng],
          [Kr[n] * Bt[m] + tu[m][n] for m in rng for n in rng],
          [tu[n][m] + Kr[n] * Bu[m] for m in rng for n in rng])
    r2 = ([tu[a][b] + tu[b][a] if a != b else tu[a][a] for a, b in pairs],
          [Kr[a] * Bu[b] + Kr[b] * Bu[a] if a != b else Kr[a] * Bu[a] for a, b in pairs],
          [Kr[b] * Bt[a] + Kr[a] * Bt[b] if a != b else Kr[a] * Bt[a] for a, b in pairs])
    if j is not None:                                        # Im parts on the cut
        isig1 = (fsi, [Z] * M, [Z] * M)
        isig2 = ([Z] * M, fsi, fsi)
        ir1 = ([Sv[n] if m == j else Z for m in rng for n in rng],
               [Bt[m] if n == j else Z for m in rng for n in rng],
               [Bu[m] if n == j else Z for m in rng for n in rng])
        ir2 = ([Z] * len(pairs),
               [(Bu[b] if a == j else Z) + (Bu[a] if b == j and a != b else Z) for a, b in pairs],
               [(Bt[a] if b == j else Z) + (Bt[b] if a == j and a != b else Z) for a, b in pairs])
    out = {}
    for I in isospins:
        c, t0 = _COMB[I], (half * (SLOTS[I][0] + SLOTS[I][1] + SLOTS[I][2]) if half else Z)
        re = [t0] + [c(*v) for grp in (sig1, sig2, r1, r2) for v in zip(*grp)]
        im = [Z] * len(re)
        if j is not None:
            im = [Z] + [c(*v) for grp in (isig1, isig2, ir1, ir2) for v in zip(*grp)]
            re, im = [scale * v for v in re], [scale * v for v in im]
        out[I] = (re, im)
    return out


def row(kind: str, I: int, ell: int, where, grid: Grid):
    """Coefficients of Re and Im of h^I_l(x_j) = omega_j f^I_l (kind 'phys', where = j, 0-based)
    or of f^I_l(s*) (kind 'unphys', where = s*, im all zero) over the 1+2M+M^2+M(M+1)/2 variables
    [T0, sigma1, sigma2, rho1 (s-index slow), rho2 upper triangle (a<b doubled)]."""
    return rows_at(grid, ell, kind, where, (I,))[I]


# --------------------------------------------------------------------------- export
SMALL = ("x", "xp", "K", "K0", "Lam", "Lambda3", "KFF3", "BSpline3", "intS0", "intP1", "intD0", "SVZSR")
PHYS = ("Re_sig", "Im_sig", "Re_rho", "Im_rho")


def matrix_shapes(grid: Grid) -> dict[str, list[int]]:
    M, n0, ns = grid.M, grid.n0, 1 + 2 * grid.M
    nr, R = grid.nvars - ns, 3 * N_WAVES * M
    shapes = {"x": [M], "xp": [M], "K": [M, M], "K0": [M], "Lam": [R], "Lambda3": [3 * n0],
              "KFF3": [3 * M, 3 * M], "BSpline3": [3 * n0, 3 * (N_CELLS + 3)], "SVZSR": [9]}
    shapes.update({k: [3, n0] for k in ("intS0", "intP1", "intD0")})
    for name, nrow in (("", R), ("fchi_", 3 * N_WAVES * 4), ("fF_", 3 * N_WAVES)):
        pre = ("Re_", "Im_") if not name else (name,)
        shapes.update({f"{p}{v}": [nrow, n] for p in pre for v, n in (("sig", ns), ("rho", nr))})
    return shapes


def _small(name: str, grid: Grid) -> list[list[arb]]:
    M, n0 = grid.M, grid.n0
    if name in ("x", "xp", "K0"):
        return [getattr(grid, name)]
    if name == "K":
        return grid.K
    if name == "Lam":
        return [[lam(l, grid.x[j]) for I in ISO_BLOCKS for l in ells(I) for j in range(M)]]
    if name == "Lambda3":
        return [[lam(l, grid.x[j]) for l in (0, 1, 2) for j in range(n0)]]
    if name == "SVZSR":
        return [svz()]
    if name in ("intS0", "intP1", "intD0"):
        return moments(grid)[name]
    blk = ff_kernel(grid) if name == "KFF3" else bsplines(grid)      # I_3 (x) block
    w = len(blk[0])
    return [[Z] * (b * w) + list(r) + [Z] * ((2 - b) * w) for b in range(3) for r in blk]


def matrix_row(name: str, r: int, grid: Grid) -> list[arb]:
    """Row r of export matrix ``name`` (r = 0 gives a 1-D matrix)."""
    ctx.prec = grid.prec
    ns = 1 + 2 * grid.M
    if name in PHYS:
        b, rest = divmod(r, N_WAVES * grid.M)
        (w, j), I = divmod(rest, grid.M), ISO_BLOCKS[b]
        v = row("phys", I, ells(I)[w], j, grid)[name.startswith("Im")]
    elif name.startswith(("fchi", "fF")):
        per = 4 if name.startswith("fchi") else 1
        b, rest = divmod(r, N_WAVES * per)
        (w, p), I = divmod(rest, per), ISO_BLOCKS[b]
        v = row("unphys", I, ells(I)[w], CHI_POINTS[p] if per == 4 else F_POINT, grid)[0]
    else:
        return _small(name, grid)[r]
    return v[:ns] if name.endswith("sig") else v[ns:]


def export(grid: Grid, outdir, digits: int = 40) -> dict:
    """Write every matrix of the export format at ``digits`` digits (generator "arb")."""
    import time
    from .kernels import ExportWriter, fmt
    t0, stats, shapes, ns = time.time(), {"zero_balls": 0}, matrix_shapes(grid), 1 + 2 * grid.M
    w = ExportWriter(outdir, {"M": grid.M, "L": 2 * N_WAVES, "nu0": grid.nu0, "digits": digits,
                              "precision": int(grid.prec * 0.30103), "prec_bits": grid.prec,
                              "generator": "arb"})
    f = lambda vals: [fmt(v, digits, stats) for v in vals]
    for name in SMALL:
        w.write(name, shapes[name], (f(r) for r in _small(name, grid)))
    for name in ("fchi_sig", "fchi_rho", "fF_sig", "fF_rho"):
        w.write(name, shapes[name], (f(matrix_row(name, r, grid)) for r in range(shapes[name][0])))
    hs = {name: w.open(name, shapes[name]) for name in PHYS}
    for b, I in enumerate(ISO_BLOCKS):
        for ell in ells(I):
            for j in range(grid.M):
                re, im = row("phys", I, ell, j, grid)
                for name, v in zip(PHYS, (re[:ns], im[:ns], re[ns:], im[ns:])):
                    hs[name].add(f(v))
    for h in hs.values():
        h.close()
    return w.finish(wall_seconds=round(time.time() - t0, 1), zero_balls=stats["zero_balls"])


def _main(argv=None) -> None:
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Arb export of the 2403 finite-model kernels")
    for flag, default in (("--M", 50), ("--nu0", NU0), ("--prec", 256), ("--digits", 40)):
        ap.add_argument(flag, type=int, default=default)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    man = export(Grid(a.M, a.nu0, a.prec), a.out, a.digits)
    print(json.dumps({k: v for k, v in man.items() if k != "matrices"}))


if __name__ == "__main__":
    _main()
