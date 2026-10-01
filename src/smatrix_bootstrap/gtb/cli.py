"""``python -m smatrix_bootstrap.run gtb <subcommand>`` for the 2403 reproduction."""
from __future__ import annotations

import argparse
import json
from dataclasses import replace

from .spec import GtbSpec


def _solve(a):
    from . import solve
    spec = replace(GtbSpec(), M=a.M, f0_ni=a.ni)
    rep = solve.run_leaf(spec, a.kernels, a.workdir, reg=a.reg, reduce_basis=a.reduce_basis, basis_tol=a.basis_tol,
                         precision=a.precision, nproc=a.nproc, timeout=a.timeout, duality_gap=a.duality_gap,
                         primal_error=a.primal_error, dual_error=a.dual_error, initial_scale=a.initial_scale,
                         keep_pmp=not a.drop_pmp, watson_from=a.watson_from)
    print(json.dumps({k: rep.get(k) for k in ("sdpb_result", "observables", "n_vars", "seconds_total")},
                     indent=1, default=str))
    return 0 if rep.get("sdpb_result", {}).get("terminateReason") == "found primal-dual optimal solution" else 1


def _certify(a):
    from pathlib import Path
    from .dual import certify_leaf
    cert = certify_leaf(a.leaf, a.kernels)
    Path(a.leaf, "certificate.json").write_text(json.dumps(cert, indent=1, default=str) + "\n")
    print(json.dumps({k: cert[k] for k in ("lower", "upper", "l4_feasible")}, indent=1, default=str))
    return 0 if cert["lower"]["certified"] else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="smatrix_bootstrap.run gtb")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("solve", help="one leaf: maximise F1 at F0 = ni/50 * 0.163108")
    s.add_argument("--kernels", required=True, help="kernel export directory (manifest.json + *.txt.gz)")
    s.add_argument("--workdir", required=True)
    s.add_argument("--M", type=int, default=50)
    s.add_argument("--ni", type=int, default=46, help="44/46/48 = blue/green/red")
    s.add_argument("--reg", default="linf-relax", choices=["linf-relax", "linf-inner", "none", "l4-exact"])
    s.add_argument("--reduce-basis", action="store_true")
    s.add_argument("--basis-tol", type=float, default=1e-12)
    s.add_argument("--precision", type=int, default=192)
    s.add_argument("--nproc", type=int, default=8)
    s.add_argument("--timeout", type=int, default=7200)
    s.add_argument("--duality-gap", type=float, default=1e-6)
    s.add_argument("--primal-error", type=float, default=1e-10)
    s.add_argument("--dual-error", type=float, default=1e-10)
    s.add_argument("--initial-scale", type=float, default=1e4)
    s.add_argument("--watson-from", default=None, help="Watson step seeded by this leaf (F0 released)")
    s.add_argument("--drop-pmp", action="store_true", help="delete pmp.json after an optimal termination")
    s.set_defaults(func=_solve)
    c = sub.add_parser("certify", help="V1 (Arb primal) and, for full-space leaves, V2 (dual bound)")
    c.add_argument("--leaf", required=True)
    c.add_argument("--kernels", required=True)
    c.set_defaults(func=_certify)
    a = p.parse_args(argv)
    return a.func(a)
