"""``python -m smatrix_bootstrap.run gtb <subcommand>`` for the 2403 reproduction."""
from __future__ import annotations

import argparse
import json
from dataclasses import replace

from .spec import GtbSpec


def spec_from_args(a) -> GtbSpec:
    """M, L and N_B from the export's manifest; F0 point and objective from the flags; FIELD=VALUE overrides."""
    from dataclasses import fields
    from pathlib import Path
    man = json.loads(Path(a.kernels, "manifest.json").read_text())
    kw = {"M": int(man["M"]), "L": int(man["L"]), "n_b": man["matrices"]["BSpline3"]["shape"][1] // 3,
          "f0_ni": a.ni, "objective": a.objective}
    types = {f.name: type(f.default) for f in fields(GtbSpec)}
    for item in getattr(a, "set", None) or []:
        name, val = item.split("=", 1)
        if name not in types or name in kw:
            raise ValueError(f"--set {name}: not an overridable GtbSpec field")
        t = types[name]
        if t is tuple:
            kw[name] = tuple(val.split(",") if name == "chi_second_norm_waves" else map(float, val.split(",")))
        else:
            kw[name] = val.lower() in ("1", "true", "yes") if t is bool else t(val)
    return replace(GtbSpec(), **kw)


def _solve(a):
    from . import solve
    spec = spec_from_args(a)
    row = None
    if a.lec:
        from .kernels import load
        from .lec import lec_rows
        from .program import Program
        row = lec_rows(Program(replace(spec, objective="F1max"), load(a.kernels), reg="none"))[a.lec]
    rep = solve.run_leaf(spec, a.kernels, a.workdir, reg=a.reg, reduce_basis=a.reduce_basis, basis_tol=a.basis_tol,
                         precision=a.precision, nproc=a.nproc, timeout=a.timeout, duality_gap=a.duality_gap,
                         primal_error=a.primal_error, dual_error=a.dual_error, initial_scale=a.initial_scale,
                         keep_pmp=not a.drop_pmp, watson_from=a.watson_from, objective_row=row,
                         hot_start=a.hot_start)
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


def _analyze(a):
    from pathlib import Path
    from .leafanalysis import analyze
    res = analyze(a.leaf, a.kernels, pole=not a.no_pole)
    Path(a.leaf, "analysis.json").write_text(json.dumps(res, default=float) + "\n")
    print(json.dumps({k: res[k] for k in ("crossings_GeV", "threshold", "radii_fm2", "lbar", "pole")}, default=float))
    return 0


def _sensitivity(a):
    from pathlib import Path
    from .sensitivity import shadow_prices
    res = shadow_prices(a.leaf, a.kernels)
    Path(a.leaf, "sensitivity.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: round(v["elasticity"], 9) for k, v in res.items()}))
    return 0


def _claims(a):
    from pathlib import Path
    from . import claims as K
    load = lambda d: json.loads(Path(d, "analysis.json").read_text())
    init, wat = {}, {}
    for spec in a.point:
        ni, i, w = spec.split(":")
        init[int(ni)], wat[int(ni)] = load(i), load(w)
    res = {}
    c1 = {}
    for spec in a.c1 or []:
        ni, lo, hi = spec.split(":")
        lo_c, hi_c = json.loads(Path(lo).read_text()), json.loads(Path(hi).read_text())
        ok = lo_c["lower"]["certified"] and lo_c["l4_feasible"] == "certified_pass" and hi_c["upper"]["certified"]
        lower, upper = float(lo_c["lower"]["value"]), float(hi_c["upper"]["value"])
        c1[int(ni)] = {"lower": lower, "upper": upper, "paper": K.C1_PAPER[int(ni)],
                       "verdict": K.c1_verdict(lower, upper, K.C1_PAPER[int(ni)]) if ok else "UNCERTIFIED"}
    verdicts = {v["verdict"] for v in c1.values()}
    res["C1"] = {"verdict": "MISSING" if not c1 else ("PASS" if verdicts == {"PASS"} else
                 "FAIL" if "FAIL" in verdicts else "INCONCLUSIVE"), "points": c1}
    res.update(K.phase_claims(wat))
    if 44 in wat and 48 in wat:
        res.update(K.table_claims(wat[44], wat[48]))
    if 48 in wat:
        res["C20"] = K.c20(wat[48])
    res["C21"] = K.c21(init, wat)
    Path(a.out).write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: v["verdict"] for k, v in res.items() if isinstance(v, dict) and "verdict" in v}))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="smatrix_bootstrap.run gtb")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("solve", help="one leaf: maximise F1 at F0 = ni/50 * 0.163108")
    s.add_argument("--kernels", required=True, help="kernel export directory (manifest.json + *.txt.gz)")
    s.add_argument("--workdir", required=True)
    s.add_argument("--M", type=int, default=50, help="ignored: M, L and N_B are read from the export's manifest")
    s.add_argument("--ni", type=int, default=46, help="F0 = ni/50 * 0.163108; 44/46/48 = blue/green/red")
    s.add_argument("--objective", default="F1max", choices=["F1max", "F1min", "custom_max", "custom_min"])
    s.add_argument("--lec", default=None, help="with custom_max/min: a00 a20 a11 b00 b20 b11 a02 a22 r2S r2V r2T")
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
    s.add_argument("--hot-start", default=None, help="start SDPB from this leaf's final checkpoint (same SDP structure)")
    s.add_argument("--set", action="append", metavar="FIELD=VALUE",
                   help="override a GtbSpec field (App. B variations), e.g. eps_chi=1e-3, r=0.02,2,6, reg_bound=500")
    s.set_defaults(func=_solve)
    c = sub.add_parser("certify", help="V1 (Arb primal) and, for full-space leaves, V2 (dual bound)")
    c.add_argument("--leaf", required=True)
    c.add_argument("--kernels", required=True)
    c.set_defaults(func=_certify)
    n = sub.add_parser("analyze", help="phases, threshold parameters, radii, fits, forward amplitudes, rho pole")
    n.add_argument("--leaf", required=True)
    n.add_argument("--kernels", required=True)
    n.add_argument("--no-pole", action="store_true")
    n.set_defaults(func=_analyze)
    e = sub.add_parser("sensitivity", help="shadow prices d opt/d theta of the tolerances from SDPB's multipliers")
    e.add_argument("--leaf", required=True)
    e.add_argument("--kernels", required=True)
    e.set_defaults(func=_sensitivity)
    k = sub.add_parser("claims", help="pre-registered verdicts from analysed points and certificates")
    k.add_argument("--point", action="append", required=True, help="NI:INITIAL_LEAF:WATSON_LEAF")
    k.add_argument("--c1", action="append", help="NI:LOWER_CERTIFICATE_JSON:UPPER_CERTIFICATE_JSON")
    k.add_argument("--out", required=True)
    k.set_defaults(func=_claims)
    a = p.parse_args(argv)
    return a.func(a)
