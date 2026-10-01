"""Independent checks of gtb.arbkernels, the export loader and GtbSpec: hand anchors, mpmath,
exact identities and a direct evaluation of the amplitude that never uses the closed forms."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import random
import shutil
import tempfile
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from flint import arb, ctx, fmpq

from smatrix_bootstrap.gtb import arbkernels as ak
from smatrix_bootstrap.gtb import kernels, spec

SCRATCH = Path(os.environ.get("GTB_SCRATCH", "/playpen1/shiqiu/gtb2403_runs/scratch_a5"))


@pytest.fixture(scope="module")
def g50():
    return ak.Grid(50)


@pytest.fixture(autouse=True)
def _arb_precision():                       # grid-free functions use the global Arb precision
    ctx.prec = 256


def _close(a, ref, rel):
    return bool(abs(a - arb(ref)) <= rel * abs(arb(ref)))


def _mpf(a):
    return mp.mpf(a.mid().str(60, radius=False))


@pytest.fixture
def scratch():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    d = Path(tempfile.mkdtemp(dir=SCRATCH))
    yield d
    shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------------------------------- loader, spec
def _write_synthetic(d: Path, mats: dict) -> None:
    entries = {}
    for name, (shape, lines) in mats.items():
        raw = gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0)
        (d / f"{name}.txt.gz").write_bytes(raw)
        entries[name] = {"shape": shape, "file": f"{name}.txt.gz",
                         "sha256": hashlib.sha256(raw).hexdigest()}
    (d / "manifest.json").write_text(json.dumps({"M": 2, "L": 20, "nu0": -20, "digits": 40,
                                                 "generator": "synthetic", "matrices": entries}))


def test_loader_round_trip(scratch):
    k = ["1234567890123456789012345678901234567890e-39", "0", "-5e-301", "25e0"]
    x = ["4005922736867778268400000000000000000000e-39", "-1e2"]
    _write_synthetic(scratch, {"K": ([2, 2], k), "x": ([2], x)})
    ks = kernels.load(scratch)
    assert ks.manifest["generator"] == "synthetic"
    assert ks.K.dtype == np.float64 and ks.K.shape == (2, 2) and ks.x.shape == (2,)
    assert ks.K[0, 0] == float("1.234567890123456789012345678901234567890")
    assert ks.K[0, 1] == 0.0 and ks.K[1, 0] == -5e-301 and ks.K[1, 1] == 25.0
    assert ks.exact("K") == k and ks.exact("x") == x
    # a tampered file must be rejected
    (scratch / "K.txt.gz").write_bytes(gzip.compress(b"1\n2\n3\n4\n", mtime=0))
    with pytest.raises(ValueError, match="sha256"):
        kernels.load(scratch)


def test_writer_loader_round_trip(scratch):
    w = kernels.ExportWriter(scratch, {"M": 1, "generator": "unit"})
    w.write("v", [3], [["1e0", "0", "-2e-3"]])
    h = w.open("A", [2, 2])
    for r in (["1e0", "2e0"], ["3e0", "4e0"]):
        h.add(r)
    h.close()
    w.finish(wall_seconds=0.0)
    ks = kernels.load(scratch)
    assert ks.v.tolist() == [1.0, 0.0, -0.002] and ks.A.tolist() == [[1, 2], [3, 4]]
    assert ks.manifest["matrices"]["A"]["shape"] == [2, 2]
    with pytest.raises(ValueError):
        kernels.ExportWriter(scratch / "bad", {}).write("v", [3], [["1e0"]])


def test_spec_defaults_and_key():
    s = spec.GtbSpec()
    assert (s.M, s.L, s.nu0, s.n_b, s.reg_norm, s.reading) == (50, 20, -20, 40, "l4", "code")
    assert s.eps_sr == (1e-7, 6e-6, 5e-6) and s.chi_second_norm_waves == ("D0", "D2", "F1")
    assert len(k := s.key()) == 16 and int(k, 16) >= 0 and k == spec.GtbSpec().key()
    assert spec.GtbSpec(f0_ni=44).key() != k


# --------------------------------------------------------------------------- anchors
def test_grid_anchors(g50):
    assert g50.n0 == 39
    assert _close(g50.x[0], "4.0059227368677782684", 2e-20)
    assert _close(g50.x[38], "172.08533939065126099", 2e-20)
    assert _close(g50.x[49], "97256.336691443898924", 2e-20)
    k0 = sorted(float(v.mid()) for v in g50.K0)
    assert abs(k0[0] + 0.0204103) < 5e-8 and abs(k0[-1] + 0.00026181) < 5e-9
    assert _close(ak.legendre_q(5, 1 + 2 * g50.x[10] / g50.xm4[3]), "9.952517107684723245837481e-13", 1e-24)
    assert _close(ak.lam(19, g50.x[0]) ** 2, "6.21643e-66", 1e-5)


def test_svz_targets():
    ref = ("3.11259180672e-7 2.07506120448e-7 1.55629590336e-7 5.61716677265e-5 2.80858338632e-5 "
           "1.87238892422e-5 5.13358746336e-5 2.56679373168e-5 1.71119582112e-5").split()
    got = ak.svz()
    assert len(got) == 9 and all(_close(a, r, 1e-11) for a, r in zip(got, ref))


def test_legendre_q_matches_mpmath():
    """Q_ell(z) = (1/2) int P_ell/(z - mu): mpmath type-3 Q, z > 1 and z < -1 (real part)."""
    rng = random.Random(11)
    mp.mp.dps = 60
    args = [fmpq(1) + fmpq(1, 10 ** k) for k in (1, 3, 5, 8)]
    args += [fmpq(10 ** 4 + rng.randint(1, 10 ** 7), 10 ** 4) for _ in range(16)]
    for z in args:
        for zz in (z, -z):
            for ell in range(20):
                a = ak.legendre_q(ell, arb(zz))
                b = mp.re(mp.legenq(ell, 0, mp.mpf(int(zz.p)) / int(zz.q), type=3))
                assert abs(float((a - arb(str(b))).mid())) <= 1e-45 * abs(float(b)), (ell, zz)


@pytest.mark.parametrize("M", [7, 50])
def test_conjugate_identity(M, g50):
    """sum_j K_ij sin(k phi_j) = cos(k phi_i) for 1 <= k <= M, fails at the band limit M+1."""
    g = g50 if M == 50 else ak.Grid(M)
    def dev(k):
        worst = arb(0)
        for i in range(M):
            lhs = sum((g.K[i][j] * arb(fmpq(k * (2 * j + 1), 2 * M)).sin_pi() for j in range(M)), arb(0))
            worst = max(worst, abs(lhs - arb(fmpq(k * (2 * i + 1), 2 * M)).cos_pi()), key=lambda v: v.mid())
        return worst
    for k in sorted({1, 2, min(10, M), M}):
        assert dev(k) < arb("1e-60"), k
    assert dev(M + 1) > arb("1e-2")


# --------------------------------------------------------------------------- spectral densities
def test_bspline_basis():
    """Clamped cubic basis on 37 uniform cells: partition of unity, closed forms, endpoints."""
    rng = random.Random(5)
    pts = [fmpq(0), fmpq(37), fmpq(1, 2), fmpq(2), fmpq(36)] + [fmpq(rng.randint(0, 37000), 1000) for _ in range(40)]
    for u in pts:
        b = ak.bspline_u(arb(u))
        assert len(b) == 40 and all(v > -arb("1e-70") for v in b)
        assert abs(sum(b, arb(0)) - 1) < arb("1e-70")
    assert abs(ak.bspline_u(arb(0))[0] - 1) < arb("1e-70") and abs(ak.bspline_u(arb(37))[39] - 1) < arb("1e-70")
    assert abs(ak.bspline_u(arb(fmpq(1, 2)))[0] - fmpq(1, 8)) < arb("1e-70")   # B_0 = (1-u)^3
    b2 = ak.bspline_u(arb(2))                                                   # cardinal B_3
    assert abs(b2[3] - fmpq(2, 3)) < arb("1e-70") and abs(b2[2] - fmpq(1, 6)) < arb("1e-70")


def test_bspline_matrix_and_moments(g50):
    B = ak.bsplines(g50)
    assert len(B) == 39 and all(len(r) == 40 for r in B)
    assert all(abs(sum(r, arb(0)) - 1) < arb("1e-60") for r in B)
    assert abs(B[38][39] - 1) < arb("1e-60")            # x_39 is the last knot
    mom = ak.moments(g50)
    assert {k: (len(v), len(v[0])) for k, v in mom.items()} == {k: (3, 39) for k in ("intS0", "intP1", "intD0")}
    # midpoint sums = int_0^{39 Delta} f dphi minus the Euler-Maclaurin midpoint terms,
    # f written analytically in phi (sqrt(x-4) = sqrt(24) tan(phi/2))
    mp.mp.dps = 30
    xn0, h = _mpf(g50.x[38]), mp.pi / 50
    kap = {"intS0": lambda x, r: 3 * r / (256 * mp.pi ** 5 * mp.sqrt(x)),
           "intP1": lambda x, r: r ** 3 / (384 * mp.pi ** 5 * mp.sqrt(x)),
           "intD0": lambda x, r: r ** 5 / (2560 * mp.pi ** 5 * mp.sqrt(x))}
    ns, off = {"intS0": (0, 1, 2), "intP1": (-1, 0, 1), "intD0": (-2, -1, 0)}, {"intS0": 2, "intP1": 2, "intD0": 3}
    for key, rows in mom.items():
        for r, n in enumerate(ns[key]):
            def f(phi):
                c = mp.cos(phi / 2)
                x, xp = -20 + 24 / c ** 2, 24 * mp.sin(phi / 2) / c ** 3
                return kap[key](x, mp.sqrt(24) * mp.tan(phi / 2)) * xp * x ** n / xn0 ** (n + off[key])
            b = 39 * h
            ref = (mp.quad(f, [0, b]) - h ** 2 / 24 * (mp.diff(f, b) - mp.diff(f, 0))
                   + 7 * h ** 4 / 5760 * (mp.diff(f, b, 3) - mp.diff(f, 0, 3)))
            assert abs(_mpf(sum(rows[r], arb(0))) / ref - 1) < 3e-5, (key, n)


# --------------------------------------------------------------------------- partial-wave rows
def _dot(row, v):
    return sum((c * x for c, x in zip(row, v) if not c.is_zero()), arb(0))


def test_constant_amplitude_rows():
    """T0 alone: f^I_l = 1/2 delta_l0 (5, 2, 0)_I T0 at nodes (as h = omega f) and at s*."""
    g = ak.Grid(6)
    e0 = [arb(1)] + [arb(0)] * (g.nvars - 1)
    for I, c in ((0, 5), (2, 2), (1, 0)):
        for ell in range(4):
            want = fmpq(c, 2) if ell == 0 else fmpq(0)
            for j in (0, 3, 5):
                re, im = ak.row("phys", I, ell, j, g)
                assert len(re) == len(im) == g.nvars == 1 + 12 + 36 + 21
                assert abs(_dot(re, e0) - g.omega[j] * want) < arb("1e-70") and _dot(im, e0) == 0
            for s in ("1/2", 1, "3/2", 2, 3):
                re, im = ak.row("unphys", I, ell, s, g)
                assert abs(_dot(re, e0) - want) < arb("1e-70") and all(v.is_zero() for v in im)


def test_functionals_of_constant_amplitude(g50):
    """F0 = 2(f0_0(3) + 5 f0_2(3)) = 5 T0 and F1 = 2(3 f1_1(3) + 7 f1_3(3)) = 0."""
    e0 = [arb(1)] + [arb(0)] * (g50.nvars - 1)
    f = {(I, ell): _dot(ak.row("unphys", I, ell, 3, g50)[0], e0) for I, ell in ((0, 0), (0, 2), (1, 1), (1, 3))}
    assert abs(2 * (f[0, 0] + 5 * f[0, 2]) - 5) < arb("1e-70") and abs(2 * (3 * f[1, 1] + 7 * f[1, 3])) < arb("1e-70")


def test_bose_symmetry_zero_rows(g50):
    """Odd waves of I=0,2 and even waves of I=1 vanish coefficient by coefficient."""
    for kind, where in (("phys", 7), ("phys", 45), ("unphys", "3/2")):
        for I, ell in ((0, 3), (2, 1), (1, 0), (1, 2)):
            ref = max(abs(c.mid()) for c in ak.row(kind, I, ell + 1, where, g50)[0])
            for part in ak.row(kind, I, ell, where, g50):
                assert all(abs(c) < 1e-60 * ref for c in part), (kind, where, I, ell)


def _mp_grid(M):
    phi = [(2 * j + 1) * mp.pi / (2 * M) for j in range(M)]
    x = [-20 + 48 / (1 + mp.cos(p)) for p in phi]
    w = [mp.pi / M * 48 * mp.sin(p) / (1 + mp.cos(p)) ** 2 for p in phi]
    c = [mp.cot(n * mp.pi / (2 * M)) / M if n % 2 else 0 for n in range(2 * M)]
    K = [[c[(i + j + 1) % (2 * M)] - c[(i - j) % (2 * M)] for j in range(M)] for i in range(M)]
    return x, w, K


def _direct_waves(M, data, s, j, ells):
    """1/4 int P_l T^I dmu from A(s,t,u) evaluated by direct sums (no Q, no partial fractions)."""
    x, w, K = _mp_grid(M)
    T0, s1, s2, r1, r2 = data
    def kv(nu):
        return [w[m] / mp.pi * (1 / (x[m] - nu) - 1 / (x[m] + 20)) for m in range(M)]
    def single(d, lab, val):
        if lab == "S" and j is not None:                 # boundary value at s + i0
            return mp.mpc(sum(K[j][m] * d[m] for m in range(M)), d[j])
        return sum(k * dm for k, dm in zip(kv(val[lab]), d))
    def double(r, la, lb, val):
        if "S" in (la, lb) and j is not None:
            if la == "S":                                   # first density index on the cut
                return sum(k * single([rm[n] for rm in r], "S", val) for n, k in enumerate(kv(val[lb])))
            return sum(k * single(r[m], "S", val) for m, k in enumerate(kv(val[la])))
        ka, kb = kv(val[la]), kv(val[lb])
        return sum(ka[m] * kb[n] * r[m][n] for m in range(M) for n in range(M))
    def amp(a, b, c, val):
        return (T0 + single(s1, a, val) + single(s2, b, val) + single(s2, c, val)
                + double(r1, a, b, val) + double(r1, a, c, val) + double(r2, b, c, val))
    memo = {}
    def T(mu):
        if mu not in memo:
            val = {"S": s, "T": -(s - 4) * (1 - mu) / 2, "U": -(s - 4) * (1 + mu) / 2}
            A, B, C = amp("S", "T", "U", val), amp("T", "S", "U", val), amp("U", "T", "S", val)
            memo[mu] = {0: 3 * A + B + C, 2: B + C, 1: B - C}
        return memo[mu]
    return {(I, ell): mp.quad(lambda mu: mp.legendre(ell, mu) * T(mu)[I], [-1, 1]) / 4
            for I in (0, 1, 2) for ell in ells}


@pytest.mark.parametrize("case", ["random", "separable_rho2"])
def test_rows_against_direct_amplitude(case):
    """Closed-form rows (Q_l, (t,u) partial fractions, crossing, isospin, rho2 folding)
    equal the projection of the directly summed amplitude, at nodes and at s* < 4."""
    M, rng = 4, random.Random(3)
    q = (lambda: fmpq(rng.randint(-999, 999), 1000)) if case == "random" else (lambda: fmpq(0))
    T0, s1, s2 = q(), [q() for _ in range(M)], [q() for _ in range(M)]
    r1 = [[q() for _ in range(M)] for _ in range(M)]
    r2 = [[fmpq(0)] * M for _ in range(M)]                  # symmetric; g(x_a) = 1/(a+2)^2
    for a, b in ((a, b) for a in range(M) for b in range(a, M)):
        r2[a][b] = r2[b][a] = q() if case == "random" else fmpq(1, (a + 2) ** 2 * (b + 2) ** 2)
    mp.mp.dps = 40
    tomp = lambda v: mp.mpf(int(v.p)) / int(v.q)
    data = (tomp(T0), [tomp(v) for v in s1], [tomp(v) for v in s2],
            [[tomp(v) for v in r] for r in r1], [[tomp(v) for v in r] for r in r2])
    g, xg = ak.Grid(M), _mp_grid(M)[0]                      # Grid sets the Arb precision
    vec = [arb(v) for v in [T0] + s1 + s2 + [r1[m][n] for m in range(M) for n in range(M)]
           + [r2[a][b] for a in range(M) for b in range(a, M)]]
    points = [("phys", 0, xg[0], 0), ("phys", 2, xg[2], 2), ("phys", 3, xg[3], 3),
              ("unphys", "1/2", mp.mpf(1) / 2, None), ("unphys", 3, mp.mpf(3), None)]
    for kind, where, s, j in points:
        scale = mp.pi * mp.sqrt((s - 4) / s) if kind == "phys" else 1
        direct = _direct_waves(M, data, s, j, range(4))
        floor = max(abs(f) for f in direct.values())            # wrong-parity waves are 0
        for (I, ell), f in direct.items():
            re, im = ak.row(kind, I, ell, where, g)
            size = sum(abs(_mpf(c * v)) for c, v in zip(re + im, vec + vec)) + floor
            for part, ref in ((re, mp.re(scale * f)), (im, mp.im(scale * f))):
                assert abs(_mpf(_dot(part, vec)) - ref) < 1e-25 * size, (case, kind, where, I, ell)


# --------------------------------------------------------------------------- exports
def _decimal(s):
    return fmpq(*Fraction(s).as_integer_ratio())


def _compare(ks, grid, per_matrix, rng):
    """Export decimals vs arbkernels: 1e-30 relative; an exported 0 needs a ball at zero."""
    bad = []
    for name, shape in ak.matrix_shapes(grid).items():
        assert ks.manifest["matrices"][name]["shape"] == shape, name
        rows, cols = (1, shape[0]) if len(shape) == 1 else shape
        exact, cache = ks.exact(name), {}
        cells = ([(r, c) for r in range(rows) for c in range(cols)] if per_matrix is None
                 else [(rng.randrange(rows), rng.randrange(cols)) for _ in range(per_matrix)])
        for r, c in cells:
            if r not in cache:
                cache[r] = ak.matrix_row(name, r, grid)
            want, got = cache[r][c], exact[r * cols + c]
            if got == "0":
                ok = 0 in want or abs(want) < 1e-40 * max(abs(v.mid()) for v in cache[r])
            else:
                ref = arb(_decimal(got))
                ok = bool(abs(want - ref) <= 1e-30 * abs(ref))
            if not ok:
                bad.append((name, r, c, got, want.str(25)))
    return bad


def test_arb_export_layout_and_round_trip(scratch):
    g = ak.Grid(3)
    man = ak.export(g, scratch)
    ks = kernels.load(scratch)
    assert man["generator"] == ks.manifest["generator"] == "arb" and ks.manifest["M"] == 3
    assert g.n0 == 2 and ks.Re_rho.shape == (90, 15) and ks.fchi_sig.shape == (120, 7)
    assert _compare(ks, g, None, None) == []
    # pinned row order: [I=0 | I=2 | I=1] x waves x nodes, s* points fastest in fchi
    re, im = ak.row("phys", 2, 2, 1, g)
    assert ks.Re_sig[30 + 3 + 1].tolist() == [float(v.mid()) for v in re[:7]]
    assert ks.Im_rho[30 + 3 + 1].tolist() == [float(v.mid()) for v in im[7:]]
    assert ks.fchi_rho[80 + 4 + 2].tolist() == [float(v.mid()) for v in ak.row("unphys", 1, 3, "3/2", g)[0][7:]]
    assert ks.fF_sig[1].tolist() == [float(v.mid()) for v in ak.row("unphys", 0, 2, 3, g)[0][:7]]
    assert ks.Lam[34] == float(ak.lam(2, g.x[1]).mid()) and ks.Lam[60 + 3 + 2] == float(ak.lam(3, g.x[2]).mid())
    assert np.allclose(ks.KFF3[3:6, 3:6], ks.K + ks.K0[None, :], rtol=1e-15) and not ks.KFF3[0:3, 3:6].any()


@pytest.mark.skipif(not os.environ.get("GTB_KERNELS"), reason="GTB_KERNELS not set")
def test_export_matches_arbkernels():
    ks = kernels.load(os.environ["GTB_KERNELS"])
    g = ak.Grid(ks.manifest["M"], ks.manifest.get("nu0", -20))
    assert _compare(ks, g, 40, random.Random(2024)) == []
