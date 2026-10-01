"""Conic-program layer of the 2403 reproduction: block algebra, eliminations, PMP round trip through SDPB."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

SCRATCH = Path("/playpen1/shiqiu/gtb2403_runs/scratch_tests")


def _docker_ok():
    return shutil.which("docker") is not None and subprocess.run(
        ["docker", "image", "inspect", "bootstrapcollaboration/sdpb:3.1.0"], capture_output=True).returncode == 0


@pytest.mark.skipif(not _docker_ok(), reason="SDPB docker image unavailable")
def test_pmp_roundtrip_known_optimum():
    """max y1 s.t. [[1, y1], [y1, 1]] >= 0 and 0 <= y2 <= 2: optimum y1 = 1 exactly."""
    from smatrix_bootstrap.gtb.program import Block
    from smatrix_bootstrap.gtb import pmpwriter
    n = 3                                   # [y0, y1, y2], y0 pinned to 1
    blocks = [Block("disk", 2, {(0, 0): ["1", "0", "0"], (0, 1): ["0", "1", "0"], (1, 1): ["1", "0", "0"]}),
              Block("y2_lo", 1, {(0, 0): ["0", "0", "1"]}),
              Block("y2_hi", 1, {(0, 0): ["2", "0", "-1"]})]
    work = SCRATCH / "pmp_roundtrip"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    info = pmpwriter.write(blocks, ["0", "1", "0"], work / "pmp.json", n_vars=n)
    assert info["n_blocks"] == 3 and info["n_vars"] == 3
    res = pmpwriter.solve(work, precision=128, nproc=1, timeout=600, duality_gap=1e-30, primal_error=1e-30, dual_error=1e-30)
    assert res["terminateReason"] == "found primal-dual optimal solution"
    y = pmpwriter.read_y(work / "out", n_vars=n)
    assert abs(float(y["text"][0]) - 1.0) < 1e-20
    assert 0.0 <= float(y["text"][1]) <= 2.0


# ---------------------------------------------------------------- the 2403 program (M = 10 export of our kernels)
import numpy as np

KS10 = Path("/playpen1/shiqiu/gtb2403_runs/kernels/ours_M10")
needs_ks10 = pytest.mark.skipif(not (KS10 / "manifest.json").exists(), reason="M=10 kernel export missing")


def _prog(**kw):
    from smatrix_bootstrap.gtb import kernels
    from smatrix_bootstrap.gtb.program import Program
    from smatrix_bootstrap.gtb.spec import GtbSpec
    ks = kernels.load(KS10)
    return Program(GtbSpec(M=10), ks, **kw), ks


def _matrix(block, yfull):
    m = np.zeros((block.size, block.size))
    for (r, c), f in block.forms.items():
        m[r, c] = m[c, r] = float(np.dot(np.asarray(f, dtype=float), yfull))
    return m


@needs_ks10
def test_lift_satisfies_the_eliminated_equalities():
    prog, ks = _prog()
    rng = np.random.default_rng(1)
    v = prog.lift(rng.normal(size=prog.n_vars - 1))
    c = v[prog.groups["c"]]
    assert np.allclose(c[[0, 40, 80]], 1.0, atol=1e-14)
    f = ks.fF_sig @ v[prog.groups["sig"]] + ks.fF_rho @ v[prog.groups["rho"]]
    assert abs(2 * (f[0] + 5 * f[1]) - prog.f0fix) < 1e-12


@needs_ks10
def test_blocks_equal_the_direct_formulas():
    """Every block evaluated from its forms equals the matrix built here from the kernels and the lifted v."""
    prog, ks = _prog()
    M, n0 = 10, ks.manifest["n0"]
    rng = np.random.default_rng(2)
    y = rng.normal(size=prog.n_vars - 1) * 1e-2
    yfull, v = np.r_[1.0, y], prog.lift(y)
    sig, rho, imf, c = (v[prog.groups[k]] for k in ("sig", "rho", "ImF", "c"))
    hre = (ks.Re_sig @ sig + ks.Re_rho @ rho) / ks.Lam
    him = (ks.Im_sig @ sig + ks.Im_rho @ rho) / ks.Lam
    hhat = (ks.Im_sig @ sig + ks.Im_rho @ rho) / ks.Lam ** 2
    ref = 1 + ks.KFF3 @ imf
    by = {b.name: b for b in prog.blocks()}
    # a unitarity cone not covered by a Gram block: I=2, l=0 (S2), node 3
    n = 10 * M + 3
    want = np.array([[2 * hhat[n], hre[n], him[n]], [hre[n], 1, 0], [him[n], 0, 1]])
    assert np.allclose(_matrix(by[f"uni:{n}"], yfull), want, atol=1e-11)
    # the Gram block of P1 at node 2 (channel 1)
    ch, j, n = 1, 2, 20 * M + 2
    F = complex(ref[ch * M + j], imf[ch * M + j])
    lam = ks.Lambda3[ch * n0 + j]
    rhat = ks.BSpline3[ch * n0 + j] @ c
    B = np.array([[1, complex(hre[n], him[n]), F],
                  [complex(hre[n], -him[n]), 2 * hhat[n], 2 * F.imag / lam],
                  [F.conjugate(), 2 * F.imag / lam, rhat]])
    assert np.allclose(_matrix(by[f"gram:P1:{j}"], yfull), np.block([[B.real, -B.imag], [B.imag, B.real]]), atol=1e-11)
    # FESR ball of D0, scaled by 1/eps
    w = ks.intD0 @ (ks.BSpline3[2 * n0:3 * n0] @ c) - ks.SVZSR[6:9]
    e = 5e-6
    want = np.eye(4); want[0, 1:] = want[1:, 0] = w / e
    assert np.allclose(_matrix(by["fesr:D0"], yfull), want, atol=1e-9)
    # FF cap of P1 at the first node above s0
    j = n0
    X = ks.x[j] * complex(ref[M + j], imf[M + j]) / (2.0 * 6.87076)
    want = np.array([[1, X.real, X.imag], [X.real, 1, 0], [X.imag, 0, 1]])
    assert np.allclose(_matrix(by[f"ff:P1:{j}"], yfull), want, atol=1e-11)


@needs_ks10
def test_chiral_blocks_are_psd_at_an_amplitude_interpolating_weinberg():
    prog, ks = _prog()
    fpi = 92 / 139.57
    s = np.array([0.5, 1.0, 1.5, 2.0])
    target = np.zeros(120)
    target[0:4] = (2 * s - 1) / (16 * np.pi ** 2 * fpi ** 2)        # f00
    target[40:44] = (2 - s) / (16 * np.pi ** 2 * fpi ** 2)          # f20 = (t+u-2)/... projected
    target[80:84] = (s - 4) / (48 * np.pi ** 2 * fpi ** 2)          # f11
    F0row = 2 * (np.hstack([ks.fF_sig[0], ks.fF_rho[0]]) + 5 * np.hstack([ks.fF_sig[1], ks.fF_rho[1]]))
    A = np.vstack([np.hstack([ks.fchi_sig, ks.fchi_rho]), F0row])
    amp = np.linalg.lstsq(A, np.r_[target, prog.f0fix], rcond=None)[0]
    v = np.zeros(prog.n_raw)
    v[prog.groups["amp"]] = amp
    v[prog.groups["c"].start + np.array([0, 40, 80])] = 1.0
    y = prog.unlift(v)
    assert np.allclose(prog.lift(y), v, atol=1e-10)
    yfull = np.r_[1.0, y]
    for name in ("chi:ratio", "chi:higher"):
        blk = next(b for b in prog.blocks() if b.name == name)
        assert np.linalg.eigvalsh(_matrix(blk, yfull)).min() > -1e-9, name


def test_cell_formats_float_arrays_exactly():
    from smatrix_bootstrap.gtb.pmpwriter import _cell
    a = np.array([0.0, 1.5, -2e-300, 0.1])
    text = _cell(a, 4, 40)
    assert text == '[["0"],["1.5"],["-2e-300"],["0.1"]]'
    assert [float(s) for s in __import__("re").findall(r'"([^"]+)"', text)] == a.tolist()


def test_cell_prints_zero_at_either_end():
    from smatrix_bootstrap.gtb.pmpwriter import _cell
    assert _cell(np.array([-0.0, 2.0, 0.0]), 3, 40) == '[["0"],["2.0"],["0"]]'
