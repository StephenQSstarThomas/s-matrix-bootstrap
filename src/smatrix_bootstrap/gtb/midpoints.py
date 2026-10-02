"""Unitarity between the nodes (WELLPOSED_ZH.md W3).

Between the nodes the discrete amplitude is the disk interpolant (evaluate.py, 'interp'), the unique band-limited
function whose boundary values reproduce the conjugate kernel K at the nodes.  At the phi-midpoints
phi = (k+1) pi/M (k = 0..M-2) this module exports the rows of h = pi sqrt((s-4)/s) f^I_l(s) (Re and Im, unrescaled,
as in the node export) with Lambda_l(s) for the same cone 2 Im h >= |h|^2, in the 40-digit export format.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ISO_ORDER = (0, 2, 1)                       # export block order of the node rows: [I=0 | I=2 | I=1]


def _rows(grid, l_max):
    """(I, l, k, s) of every exported row, in export order."""
    from flint import arb
    M, out = grid.M, []
    for I in ISO_ORDER:
        for ell in range(1 if I == 1 else 0, l_max, 2):
            for k in range(M - 1):
                phi = arb.pi() * (k + 1) / M
                out.append((I, ell, k, arb(grid.nu0) + (8 - 2 * arb(grid.nu0)) / (1 + phi.cos())))
    return out


def build(kernels_dir, out_dir, l_max=20, prec_bits=192, digits=40):
    from flint import acb, arb
    from .arbkernels import Grid
    from .evaluate import f_row
    from .kernels import ExportWriter, fmt, load
    man = load(kernels_dir).manifest
    grid = Grid(int(man["M"]), prec_bits=prec_bits)
    rows = _rows(grid, l_max)
    w = ExportWriter(out_dir, {"M": grid.M, "nu0": grid.nu0, "l_max": l_max, "kernels": str(kernels_dir),
                               "generator": "gtb.midpoints", "digits": digits})
    n_amp = len(f_row(grid, 0, 0, ("cut", rows[0][3]), "interp"))
    hre, him = w.open("Re", (len(rows), n_amp)), w.open("Im", (len(rows), n_amp))
    lam, svals = [], []
    for I, ell, k, s in rows:
        om = acb.pi() * ((s - 4) / s).sqrt()
        h = [om * c for c in f_row(grid, I, ell, ("cut", s), "interp")]
        hre.add([fmt(c.real, digits) for c in h])
        him.add([fmt(c.imag, digits) for c in h])
        r = s.sqrt()
        lam.append(fmt(((r - 2) / (r + 2)) ** (arb(ell) / 2), digits))
        svals.append(fmt(s, digits))
    hre.close(), him.close()
    w.write("Lam", (len(rows),), [[v] for v in lam])
    w.write("s", (len(rows),), [[v] for v in svals])
    w.write("meta", (len(rows), 3), [[str(I), str(ell), str(k)] for I, ell, k, _ in rows])
    return w.finish()


def load(mid_dir, rows=None, exact=False) -> dict:
    """{'Re', 'Im', 'Lam', 's', 'meta'} (float arrays, or Arb object arrays with exact=True), optionally a subset."""
    from .kernels import load as kload
    ks = kload(mid_dir)
    if exact:
        from .certify import ball
        get = lambda n: np.array([ball(t) for t in ks.exact(n)], dtype=object).reshape(ks.manifest["matrices"][n]["shape"])
    else:
        get = lambda n: getattr(ks, n)
    out = {n: get(n) for n in ("Re", "Im", "Lam", "s")}
    out["meta"] = np.asarray(ks.meta, dtype=float).astype(int)
    if rows is not None:
        idx = np.asarray(rows, dtype=int)
        out = {n: v[idx] for n, v in out.items()}
        out["rows"] = idx
    return out


def scan(mid, amp) -> np.ndarray:
    """|S| = |1 + i h| at every row for a float amplitude vector."""
    h = mid["Re"] @ amp + 1j * (mid["Im"] @ amp)
    return np.abs(1 + 1j * h)


def violated(mid, amp, tol=1e-6, previous=()):
    """Rows with |S| > 1 + tol, joined with the previous selection (constraint generation, W3)."""
    S = scan(mid, amp)
    return sorted(set(int(i) for i in np.flatnonzero(S > 1 + tol)) | set(int(i) for i in previous))


def save_rows(path, rows, note=""):
    Path(path).write_text(json.dumps({"rows": [int(r) for r in rows], "note": note}) + "\n")
