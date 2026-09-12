"""Entry point:  ``python -m smatrix_bootstrap.sdp <command> [options]``.

Commands
    selfcheck   run the section 5a mathematical self-checks and cross-checks
    prereg      write preregistration.json (must precede P3)
    solve       one or more boundary directions / representative points
    figures     build the PDF figures and the comparison tables from report.json

A ``solve`` invocation writes exactly one ``report.json`` into ``--out``; every
number in the final report cites one of those files.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from . import constraints as C
from .problem import ModelSpec
from .runner import run_job, sweep_directions

POINTS = ("tip", "ref", "mid")


def _spec(a) -> ModelSpec:
    return ModelSpec(M=a.M, L=a.L, chiral=a.chiral, chi_caliber=a.chi, eps_chi=a.eps_chi,
                     uv=a.uv, sr_caliber=a.sr, eps_ff=a.eps_ff, B=a.B, B_norm=a.B_norm,
                     cone_scaling=a.cone_scaling, sparsify=a.sparsify,
                     reduce_basis=a.reduce_basis, tag=a.tag)


def _jobs(a, x_tip: float | None):
    """Build the (name, direction, extra) list for this invocation."""
    if a.points:
        x_ref, _ = C.chiral_reference_point()
        out = []
        for p in sorted(a.points, key=lambda q: {"tip": 0, "ref": 1, "mid": 2}[q]):
            if p == "tip":
                out.append(("tip", (1.0, 0.0), None))
            elif p == "ref":
                out.append(("ref", (0.0, 1.0), _fix_f00(lambda ctx: x_ref)))
            else:
                # the mid section needs the tip of *this* feasible set; take it
                # from the tip solve in this process, else from --x-tip
                def _mid(ctx, _xt=x_tip, _xr=x_ref):
                    xt = ctx["tip"]["f00_3"] if "tip" in ctx and ctx["tip"].get("f00_3") \
                        is not None else _xt
                    if xt is None:                       # tip solve failed
                        raise _SkipJob("mid needs a tip solve or --x-tip")
                    return 0.5 * (xt + _xr)
                out.append(("mid", (0.0, 1.0), _fix_f00(_mid)))
        return out
    if a.section is not None:
        # the width of the allowed region on the vertical section f00(3) = x:
        # max f11 (direction +y) and min f11 (direction -y)
        x = C.chiral_reference_point()[0] if a.section == "xref" else float(a.section)
        return [("section_hi", (0.0, 1.0), _fix_f00(lambda ctx, _x=x: _x)),
                ("section_lo", (0.0, -1.0), _fix_f00(lambda ctx, _x=x: _x))]
    dirs = sweep_directions(a.ndir, half=a.half)
    lo, hi = (a.slice or "0:%d" % len(dirs)).split(":")
    sel = list(range(len(dirs)))[int(lo):int(hi)]
    return [("dir%03d" % i, dirs[i], None) for i in sel]


class _SkipJob(Exception):
    """Raised when a chained job cannot be set up (e.g. its tip solve failed)."""


def _fix_f00(x):
    """``x`` may be a number or a callable of the accumulated results."""
    def extra(model, ctx):
        v = x(ctx) if callable(x) else x
        return [model.f00 == v]
    return extra


def main(argv=None) -> int:
    p = argparse.ArgumentParser("smatrix_bootstrap.sdp")
    p.add_argument("command", choices=["selfcheck", "prereg", "solve", "figures"])
    p.add_argument("--out", default=None)
    p.add_argument("--M", type=int, default=50)
    p.add_argument("--L", type=int, default=10)
    p.add_argument("--chiral", action="store_true")
    p.add_argument("--chi", default="chi-b", choices=["chi-a", "chi-b", "chi-c"])
    p.add_argument("--eps-chi", type=float, default=C.EPS_CHI_MAIN)
    p.add_argument("--uv", action="store_true")
    p.add_argument("--sr", default="SR-b", choices=["SR-a", "SR-b", "SR-c"])
    p.add_argument("--eps-ff", type=float, default=C.EPS_FF)
    p.add_argument("--B", type=float, default=None)
    p.add_argument("--B-norm", default="l2", choices=["l2", "l4"])
    p.add_argument("--cone-scaling", default="rownorm",
                   choices=["none", "centrifugal", "rownorm"])
    p.add_argument("--sparsify", type=float, default=0.0)
    p.add_argument("--reduce-basis", action="store_true")
    p.add_argument("--ndir", type=int, default=24)
    p.add_argument("--half", action="store_true")
    p.add_argument("--slice", default=None)
    p.add_argument("--points", nargs="*", choices=list(POINTS))
    p.add_argument("--x-tip", type=float, default=None)
    p.add_argument("--section", default=None,
                   help="'xref' or a number: solve max/min f11 on f00(3) = x")
    p.add_argument("--solver", default="CLARABEL")
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--time-limit", type=float, default=7200.0)
    p.add_argument("--tag", default="")
    p.add_argument("--objective", default="plane", choices=["plane", "lambda"])
    p.add_argument("--generate", action="store_true",
                   help="constraint generation over the unitarity disks")
    p.add_argument("--start-tol", type=float, default=1e-6)
    a = p.parse_args(argv)

    if a.command == "selfcheck":
        import subprocess
        import sys
        here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        return subprocess.call([sys.executable, "-m", "pytest", "-q",
                                os.path.join(here, "tests", "sdp")])

    if a.command == "prereg":
        doc = preregistration()
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "preregistration.json"), "w") as fh:
            json.dump(doc, fh, indent=1)
        print(json.dumps(doc, indent=1))
        return 0

    if a.command == "solve":
        if not a.out:
            raise SystemExit("--out is required")
        kw = {"max_iter": a.max_iter, "time_limit": a.time_limit}
        if a.solver == "SCS":
            kw = {"eps": 1e-7, "max_iters": 200000}
        if a.generate:
            kw["start_tol"] = a.start_tol
        run_job(_spec(a), _jobs(a, a.x_tip), a.out, solver=a.solver,
                objective=a.objective, generate=a.generate, **kw)
        print("wrote", os.path.join(a.out, "report.json"))
        return 0

    from .figures import build_all
    build_all(a.out)
    return 0


def preregistration() -> dict:
    """The three under-determined calibers and the acceptance thresholds.

    Frozen before P3; nothing here may be added to, removed or re-tuned
    afterwards (task section 4 and section 7).
    """
    return {
        "written_before": "P3",
        "paper": "arXiv:2309.12402v3",
        "calibers": {
            "chiral_norm": {
                "source": "(3.64) says 'with some norm'",
                "chi-a": "literal L-infinity box: 8 separate |r| <= eps_chi",
                "chi-b": "single 8-dimensional L2 ball (the authors' 2403 code caliber)",
                "chi-c": "two separate 4-dimensional L2 balls",
            },
            "fesr_tolerance": {
                "source": "(3.74) and (2.56) do not close dimensionally",
                "SR-a": "raw absolute: |M_n - T_n| <= 2e-3 per moment",
                "SR-b": "relative 10%: |M_n - T_n| <= 0.10 |T_n|",
                "SR-c": "relative 20%",
            },
            "density_regularisation_B": {
                "source": "absent from the paper; present in the authors' 2403 code as "
                          "norm(rho/Mrho,4) <= 1e2, i.e. ||rho||_4 <= 377500",
                "default": "absent",
                "grid_if_needed": [1e4, 1e5, 377500.0, 1e6],
                "must_report": "whether B is active (||rho||_4 >= 0.9 B)",
            },
        },
        "primary_caliber": {"chiral_norm": "chi-b", "fesr_tolerance": "SR-b"},
        "grid": "3 x 3 calibers on the x_ref section and the three representative points",
        "fixed_inputs": {"M": 50, "L": 10, "s0": C.S0, "alpha_s": 0.4,
                         "m_q_mean": C.M_Q, "m_q_rms": C.M_Q_RMS,
                         "eps_chi_main": C.EPS_CHI_MAIN, "eps_chi_grid": list(C.EPS_CHI_GRID),
                         "eps_sr": C.EPS_SR, "eps_ff": C.EPS_FF,
                         "projection_plane": "(f00(3), f11(3))",
                         "chiral_points": list(C.CHIRAL_POINTS)},
        "representative_points": {
            "tip": "argmax f00(3)",
            "ref": "max f11(3) on the section f00(3) = x_ref = 0.0733214",
            "mid": "max f11(3) on the section f00(3) = (x_tip + x_ref)/2",
        },
        "acceptance": {
            "C1": "f00 extrema in [-2.902, 2.233] +-2%; f11 in [-0.734, 0.079] +-2%",
            "C2": "eps=2e-3 +x end 0.0826 +-5%; x_ref section width 0.00076 +-20%; "
                  "monotone in eps",
            "C3": "subthreshold rms <= 8% of f00(3); S0 zero near 0.425 / 0.305 / none",
            "C4": "delta00(0.9 GeV) in [85,110] deg; delta11(1.2 GeV) <= 25 deg",
            "C5": "x_ref section: upper shrinks 2.3e-4 +-25%, lower rises <= 1.0e-4, "
                  "ratio >= 4; UV +x end 0.0811 +-5%",
            "C6": "90 deg crossing in [795, 845] MeV for all three points, spread "
                  "<= 20 MeV, min eta(P1) >= 0.9",
            "C7": "delta00(1.196 GeV) in [85,110] deg; delta20 in [-40,-15] deg; "
                  "node rms <= 10 deg",
            "C8": "L in {8,10,12}: rho spread <= 20 MeV; M in {45,50,60}: spread "
                  "40-70 MeV with M60 < M45 < M50; five delta00(1.0 GeV) in [85,105]",
        },
        "prediction": "under SR-b/c the ref/mid rho lands at 822 +- 40 MeV, "
                      "min eta(P1) >= 0.9, and the S0 spectrum returns to ~950 MeV",
    }


if __name__ == "__main__":
    raise SystemExit(main())
