"""S3: the caliber-indexed result surface, and the paper's central contrast.

The paper's whole argument is one contrast run twice with the same operators:

    chiral only            ->  P1 has no rho        (Fig. 7)
    chiral + UV sector     ->  P1 acquires a rho    (Fig. 9)

Three of the calibers that enter it are not uniquely fixed by the paper
(TASK_PHYSICS_CALIBERS_ZH.md Q1/Q2/Q3).  Rather than pick one, this driver runs
the contrast on the whole grid and reports where it holds.  If it holds in
*every* cell, the caliber question does not need answering for the paper's
claim; if it holds only in some, the physics side adjudicates among survivors.

Nothing here chooses a caliber, and nothing here is a reproduction verdict.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np


def _one(task):
    """One grid cell, solved by constraint generation.

    Constraint generation is load-bearing, not an optimisation: imposing all
    3LM disks at once is what the previous round found to be unsolvable, while
    starting from the ~100 disks whose operator norm is within start_tol of the
    largest and adding violators converges and *certifies* (see
    runner.solve_generated).  A one-shot solve of the same model returns
    SolverError even for the chiral-only cells that are known to close.
    """
    from smatrix_bootstrap.sdp import constraints as C
    from smatrix_bootstrap.sdp import observables as O
    from smatrix_bootstrap.sdp.problem import ModelSpec
    from smatrix_bootstrap.sdp.runner import solve_generated

    (M, L, chi, sr, B, eps_ff_mult, uv, solver, max_iter, time_limit) = task
    spec = ModelSpec(M=M, L=L, chiral=True, chi_caliber=chi, uv=uv,
                     sr_caliber=sr, eps_ff=C.EPS_FF * eps_ff_mult,
                     B=B, B_norm="l4" if B else "l2", cone_scaling="rownorm")
    rec = {"M": M, "L": L, "chi": chi, "sr": sr if uv else None, "B": B,
           "eps_ff_mult": eps_ff_mult if uv else None,
           "sector": "chiral+UV" if uv else "chiral", "solver": solver}
    t0 = time.time()
    try:
        # Use problem.solver_options' defaults (tol_gap_rel = 1e-8).  Tightening
        # MOSEK's conic relative gap to 1e-10 makes it fail on the very first
        # generation round of problems it otherwise solves and certifies -- a
        # measured, reproducible effect, not a guess.
        res, model, rounds = solve_generated(spec, (1.0, 0.0), solver=solver,
                                             max_iter=max_iter, time_limit=time_limit)
        rec["status"] = res.get("status")
        rec["certified"] = res.get("certified")
        rec["objective"] = res.get("objective")
        rec["generation_rounds"] = len(rounds)
        rec["n_disks_imposed"] = res.get("n_disks_imposed")
        rec["max_violation_at_return"] = res.get("max_violation_at_return")
        rec["stopped_because"] = res.get("stopped_because")
        sol = model.solution()
        if sol.get("c") is not None:
            rec["f00_3"], rec["f11_3"] = res.get("f00_3"), res.get("f11_3")
            rr = O.resonance_report(model.ops, sol["c"])
            rec["resonance"] = {
                w: {k: v[k] for k in ("crossing_90_GeV", "modulus_peak_GeV",
                                      "min_eta_below_1p2GeV", "first_zero_node")}
                for w, v in rr.items()}
            for w, e in (("S0", 0.9), ("P1", 1.2), ("S2", 1.196)):
                full = rr[w]
                rec["resonance"][w][f"delta_at_{e}GeV"] = float(
                    np.interp(e, full["E_GeV"], full["delta_deg"]))
    except Exception as exc:
        rec["status"] = "SolverError"
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:160]}"
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def main(argv=None) -> int:
    p = argparse.ArgumentParser("caliber_grid")
    p.add_argument("--M", type=int, default=20)
    p.add_argument("--L", type=int, default=6)
    p.add_argument("--chi", nargs="+", default=["chi-a", "chi-b", "chi-c"])
    p.add_argument("--sr", nargs="+", default=["SR-a", "SR-b", "SR-c"])
    p.add_argument("--B", nargs="+", type=float, default=[377500.0])
    p.add_argument("--no-B", action="store_true", help="also run with B absent")
    p.add_argument("--eps-ff-mult", nargs="+", type=float, default=[1.0])
    p.add_argument("--solver", default="MOSEK")
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--time-limit", type=float, default=900.0)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    Bs = list(a.B) + ([None] if a.no_B else [])
    tasks = []
    for chi, B in itertools.product(a.chi, Bs):           # link (4): chiral only
        tasks.append((a.M, a.L, chi, a.sr[0], B, 1.0, False, a.solver,
                      a.max_iter, a.time_limit))
    for chi, sr, B, m in itertools.product(a.chi, a.sr, Bs, a.eps_ff_mult):
        tasks.append((a.M, a.L, chi, sr, B, m, True, a.solver,
                      a.max_iter, a.time_limit))

    print(f"{len(tasks)} cells, {a.jobs} workers", flush=True)
    recs = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for rec in ex.map(_one, tasks):
            recs.append(rec)
            r = rec.get("resonance") or {}
            rho = (r.get("P1") or {}).get("crossing_90_GeV")
            print(f"{rec['sector']:10s} chi={rec['chi']} sr={str(rec['sr']):5s} "
                  f"B={'yes' if rec['B'] else 'no ':3s}  {str(rec['status']):20s} "
                  f"cert={str(rec.get('certified')):5s} "
                  f"f00={rec.get('f00_3')} rho={rho}  {rec['seconds']}s", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump({"records": recs}, open(a.out, "w"), indent=2)
    ok = sum(r["status"] in ("optimal", "optimal_inaccurate") for r in recs)
    print(f"\nwrote {a.out}   solved {ok}/{len(recs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
