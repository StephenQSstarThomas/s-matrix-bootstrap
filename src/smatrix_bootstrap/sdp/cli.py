"""Self-contained SDPB CLI. All solves use the same finite 2309 model."""
from __future__ import annotations

import argparse
import os
from dataclasses import asdict
from pathlib import Path

from . import constraints as C
from .sdpb import Settings, generate, run_once, write_json
from .spec import ModelSpec


def main(argv=None):
    p = argparse.ArgumentParser("smatrix_bootstrap.run sdp solve (SDPB only)")
    p.add_argument("--workdir", "--out", dest="workdir", required=True)
    p.add_argument("--M", type=int, default=50)
    p.add_argument("--L", type=int, default=10)
    p.add_argument("--chiral", action="store_true")
    p.add_argument("--chi", choices=["chi-a", "chi-b", "chi-c"], default="chi-b",
                   help="chi-b: one combined 8-dim L2 norm (the authors' released convention); "
                        "chi-a: per-residual box; chi-c: two separate 4-dim L2 norms")
    p.add_argument("--eps-chi", type=float, default=C.EPS_CHI_MAIN)
    p.add_argument("--uv", action="store_true")
    p.add_argument("--uv-parts", nargs="+", choices=["gram", "fesr", "ff"], default=["gram", "fesr", "ff"])
    p.add_argument("--sr", choices=["SR-a", "SR-b", "SR-c", "SR-d"], default="SR-a",
                   help="FESR tolerance packaging: SR-a per-moment box 2e-3 (literal), SR-d per-wave L2 ball 2e-3 "
                        "(authors' code structure with the 2309 value), SR-b/SR-c relative 10%%/20%%")
    p.add_argument("--sr-free", nargs=2, action="append", metavar=("WAVE", "N"), default=[],
                   help="moment-range diagnostic: do not impose the FESR box of this (wave, n); repeatable")
    p.add_argument("--functional", nargs=4, metavar=("KIND", "WAVE", "NODE", "SENSE"),
                   help="replace the support objective by +-one linear functional (kinds ImKH ImS ImF rho SRmom); "
                        "with --sr-free this gives the achievable range of a freed moment")
    p.add_argument("--eps-ff", type=float, default=C.EPS_FF)
    p.add_argument("--eps-sr", type=float, default=C.EPS_SR, help="SR-a box half-width / SR-d ball radius (2309: 2e-3)")
    p.add_argument("--eps-ff-s0", type=float, default=None, help="sensitivity study: eps_FF for the S0 current only")
    p.add_argument("--eps-ff-p1", type=float, default=None, help="sensitivity study: eps_FF for the P1 current only")
    p.add_argument("--mq", choices=["mean", "rms"], default="mean")
    p.add_argument("--ff-factor", choices=["frozen", "node"], default="node")
    p.add_argument("--reg-norm", choices=["none", "linf", "l2", "l4"], default="none",
                   help="M-regularisation of 2103.11484 sec. 3 on the double spectral density: "
                        "linf caps every |rho_{a,ij}| by --reg-bound (no extra variables); "
                        "l2/l4 bound the norm with one/two auxiliary variables per density value")
    p.add_argument("--reg-bound", type=float, default=None,
                   help="Mreg for --reg-norm; fixed by the plateau rule, never by output curves")
    p.add_argument("--cone-scaling", choices=["rownorm", "centrifugal", "none"], default="rownorm")
    p.add_argument("--reduce-basis", action="store_true")
    p.add_argument("--basis-tol", type=float, default=1e-12)
    p.add_argument('--basis-source-report',help='use exactly an accepted source basis; assemble and verify fresh constraints')
    p.add_argument("--operator-dps", type=int, default=40)
    p.add_argument('--scattering-prescription',choices=['sine-cardinal','mixed-pv'],default='mixed-pv',
                   help="mixed-pv is the authors' nodal discretisation (cot kernel, midpoint rule, "
                        "Legendre-Q); sine-cardinal is a declared alternative interpolation")
    p.add_argument("--direction", type=float, nargs=2, default=[1., 0.])
    p.add_argument("--fix-f00", type=float)
    p.add_argument("--points", nargs="+", choices=["tip", "ref", "mid"])
    p.add_argument("--x-tip", type=float)
    p.add_argument("--generate", action="store_true")
    p.add_argument("--start-tol", type=float, default=0., help="0 starts with every native disk")
    p.add_argument("--max-rounds", type=int, default=12)
    p.add_argument("--add-per-round", type=int, default=120)
    p.add_argument("--precision", type=int, default=192)
    p.add_argument("--nproc", type=int, default=8, help="MPI ranks; not a count of machines")
    p.add_argument("--digits", type=int, default=30)
    p.add_argument("--initial-scale", type=float, default=1e4)
    p.add_argument("--max-iterations", type=int, default=500)
    p.add_argument("--duality-gap", type=float, default=1e-6)
    p.add_argument("--primal-error", type=float, default=1e-10)
    p.add_argument("--dual-error", type=float, default=1e-10)
    p.add_argument("--timeout", type=int, default=7200)
    p.add_argument('--checkpoint-interval',type=int,default=300)
    p.add_argument("--skip-mma-audit", action="store_true")
    p.add_argument('--mma-reference',help='reuse a passed same-script Mathematica report and recheck current formulas')
    p.add_argument("--keep", action="store_true", help="all evidence is always kept")
    p.add_argument("--no-crosscheck", action="store_true", help="compatibility; no other solver is called")
    a = p.parse_args(argv)
    if a.M < 2 or a.L < 1 or a.nproc < 1 or a.max_rounds < 1 or a.add_per_round < 1:
        p.error("M>=2, L>=1 and positive process/round counts required")
    if a.reg_norm != "none" and (a.reg_bound is None or not a.reg_bound > 0):
        p.error("--reg-norm needs a positive --reg-bound (Mreg)")
    if a.reg_norm == "none" and a.reg_bound is not None:
        p.error("--reg-bound given without --reg-norm")
    if a.basis_source_report and not a.reduce_basis:p.error('--basis-source-report requires --reduce-basis')
    if a.points and a.fix_f00 is not None:
        p.error("--points and --fix-f00 are distinct selection rules")
    if a.sr_free and not (a.uv and "fesr" in a.uv_parts):
        p.error("--sr-free needs --uv with the fesr part")
    for w, n in a.sr_free:
        if w not in ("S0", "P1") or not n.lstrip("-").isdigit():
            p.error("--sr-free takes WAVE in {S0,P1} and an integer N")
    functional = None
    if a.functional is not None:
        if a.points or a.generate:
            p.error("--functional replaces the objective; it takes no --points and no --generate")
        kind, wave, node, sense = a.functional
        if not node.lstrip("-").isdigit():
            p.error("--functional NODE must be an integer")
        functional = {"kind": kind, "wave": wave, "node": int(node), "sense": sense}
    if a.points and "mid" in a.points and "tip" not in a.points and a.x_tip is None:
        p.error("mid needs this model's accepted tip or an explicitly supplied --x-tip")
    cfg = Settings(**{key: getattr(a, key) for key in asdict(Settings())})
    spec = ModelSpec(M=a.M, L=a.L, chiral=a.chiral, chi_caliber=a.chi,
        eps_chi=a.eps_chi, uv=a.uv, uv_parts=tuple(a.uv_parts), sr_caliber=a.sr,
        eps_ff=a.eps_ff, eps_sr=a.eps_sr, eps_ff_s0=a.eps_ff_s0, eps_ff_p1=a.eps_ff_p1, m_q=C.M_Q if a.mq == "mean" else C.M_Q_RMS,
        ff_frozen_at_s0=a.ff_factor == "frozen", cone_scaling=a.cone_scaling,
        reduce_basis=a.reduce_basis, basis_tol=a.basis_tol, operator_dps=a.operator_dps,
        scattering_prescription=a.scattering_prescription,
        reg_norm=None if a.reg_norm == "none" else a.reg_norm, reg_bound=a.reg_bound,
        sr_free=tuple((w, int(n)) for w, n in a.sr_free))
    if a.basis_source_report:
        from .assembly import load_saved_basis
        load_saved_basis(spec,a.basis_source_report)
    root = Path(a.workdir).resolve()
    if root.exists() and any(root.iterdir()):
        p.error("Use a fresh --workdir; existing evidence is never overwritten")
    root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("SDP_CACHE", str(root / "operator-cache"))
    write_json(root / "selection_rule.json", {"spec": asdict(spec), "points": a.points,
        "direction": a.direction, "fix_f00": a.fix_f00, "x_tip_input": a.x_tip, "functional": functional,
        "basis_source_report":str(Path(a.basis_source_report).resolve()) if a.basis_source_report else None,
        "ref": C.chiral_reference_point()[0], "mid_rule": "(this-set tip + ref)/2",
        "frozen_before_solving": True,
        "source_reading": "2309v3 inputs; chiral norm, FESR packaging and the 2103.11484 "
                          "regulariser are declared conventions recorded in spec",
        "selection_scope": "registered section rule; not proof of nearest boundary point to physical f_pi"})
    if not a.skip_mma_audit:
        from .mma import audit, reuse_audit
        mma = (reuse_audit(a.M,root/'mma',a.mma_reference) if a.mma_reference else audit(a.M, root/'mma'))
        if not mma["passed"]:
            write_json(root / "report.json", {"status": "mma_audit_failed", "mma": mma})
            return 1
        print("Mathematica 2309 audit passed", flush=True)
    jobs = [x for x in ("tip", "ref", "mid") if x in a.points] if a.points else ["support"]
    reports, xtip = {}, a.x_tip
    for label in jobs:
        direction, fixed = a.direction, a.fix_f00
        if label == "tip":
            direction, fixed = (1., 0.), None
        elif label in ("ref", "mid"):
            direction = (0., 1.)
            xr = C.chiral_reference_point()[0]
            fixed = xr if label == "ref" else (xr + xtip)/2
        if a.generate:
            rec = generate(spec, root / label, direction, fixed, cfg,
                           a.start_tol, a.max_rounds, a.add_per_round,basis_source_report=a.basis_source_report)
        else:
            rec, _, _ = run_once(spec, root / label, direction, fixed, cfg,basis_source_report=a.basis_source_report,
                                 functional=functional)
        reports[label] = rec
        write_json(root / "report.json", {"solver": "SDPB", "reports": reports,
                   "accepted": all(r["accepted"] for r in reports.values()) and len(reports)==len(jobs)})
        if not rec["accepted"]:
            return 1
        if label == "tip":
            xtip = rec["verification"]["f00_3"]
    return 0
