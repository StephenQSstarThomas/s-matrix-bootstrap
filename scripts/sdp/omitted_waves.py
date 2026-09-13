"""Unitarity of the partial waves that a finite-L solve does not impose.

Usage:  python scripts/sdp/omitted_waves.py LEAF/report.json [--ells N] [--out FILE]

Reconstructs the full 3876-coefficient amplitude from the saved SDPB solution
(basis.npy and out/y.txt) and evaluates, with the float mixed-pv projector, the
partial waves ell = 2L .. 2L+2N-1 (the first N omitted waves per isospin) at all
M nodes.  Reports max eta = |S| per wave with its node, plus the coefficient
scales.  A regularised amplitude should keep these close to or below 1; the
un-regularised SDPB solutions had eta ~ 50.  Diagnostic only: it adds no
constraint and never changes a solve.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from smatrix_bootstrap.sdp.grid import kappa, s_nodes
from smatrix_bootstrap.sdp.projector import Layout, PartialWaveOperator, ells_for


def load_amplitude(leaf: Path) -> tuple[dict, np.ndarray]:
    rec = json.loads(leaf.read_text())
    if rec.get("scattering_prescription", rec["spec"].get("scattering_prescription")) not in (None, "mixed-pv"):
        raise ValueError("float re-projection is defined for the mixed-pv family only")
    y = np.array([float(t) for t in (leaf.parent / "out/y.txt").read_text().split()[2:]])
    n_a = rec["pmp"]["n_a"]
    if rec["spec"]["reduce_basis"]:
        c = np.load(leaf.parent / "basis.npy") @ y[:n_a]
    else:
        c = y[:n_a]
    return rec, c


def omitted_wave_report(rec: dict, c: np.ndarray, n_waves: int = 3) -> dict:
    M, L = rec["spec"]["M"], rec["spec"]["L"]
    op, lay, s = PartialWaveOperator(M), Layout(M), s_nodes(M)
    kap = kappa(s)
    rho = np.r_[c[lay.r1], c[lay.r2]]
    out = {"M": M, "L": L, "imposed_ells": {I: list(ells_for(I, L)) for I in (0, 1, 2)},
           "coefficient_scales": {"T0": float(c[0]), "sigma1_max": float(np.abs(c[lay.s1]).max()),
                                  "sigma2_max": float(np.abs(c[lay.s2]).max()), "rho_linf": float(np.abs(rho).max()),
                                  "rho_l2": float(np.linalg.norm(rho)), "rho_l4": float(np.sum(rho ** 4) ** 0.25)},
           "waves": []}
    reg = rec["spec"].get("reg_bound")
    if reg:
        out["coefficient_scales"]["regulariser_active_fraction"] = float(np.mean(np.abs(rho) >= 0.99 * reg))
    for I in (0, 1, 2):
        first = 2 * L + (I % 2)
        for ell in range(first, first + 2 * n_waves, 2):
            eta = np.empty(M)
            for k in range(M):
                r = op.rows(I, ell, float(s[k]), k)
                h = kap[k] * (r[0] @ c + 1j * (r[1] @ c))
                eta[k] = abs(1.0 + 1j * h)
            k = int(np.argmax(eta))
            out["waves"].append({"isospin": I, "ell": ell, "max_eta": float(eta[k]), "node": k,
                                 "E_GeV": float(0.14 * np.sqrt(s[k])), "eta": eta.tolist()})
    out["max_eta_all_omitted"] = max(w["max_eta"] for w in out["waves"])
    out["scope"] = "float mixed-pv re-projection of the saved amplitude; diagnostic only"
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("leaf")
    p.add_argument("--ells", type=int, default=3, help="omitted waves per isospin to check")
    p.add_argument("--out")
    a = p.parse_args(argv)
    leaf = Path(a.leaf).resolve()
    rec, c = load_amplitude(leaf)
    rep = omitted_wave_report(rec, c, a.ells)
    rep["source_report"] = str(leaf)
    dest = Path(a.out) if a.out else leaf.parent / "omitted_waves.json"
    dest.write_text(json.dumps(rep, indent=1))
    for w in rep["waves"]:
        print(f"I={w['isospin']} ell={w['ell']:2d} max eta {w['max_eta']:.4f} at node {w['node']} ({w['E_GeV']:.3f} GeV)")
    print(f"max eta over omitted waves: {rep['max_eta_all_omitted']:.4f}; rho_linf {rep['coefficient_scales']['rho_linf']:.4g}; wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
