"""Reproducible, source-led audit calculations; no optimization in this module."""
from pathlib import Path
import hashlib
import json

import numpy as np

from .grid import phi_nodes
from .hilbert import cauchy_offcut_row, cauchy_oncut_real, infinity_constant_row
from .runner import provenance, save_producer


def constant_identity(outdir):
    """A1: compare independent formulas for every node and sine mode."""
    rows = []
    for M in (20, 30, 50, 100):
        p = phi_nodes(M)
        n = np.arange(1, M + 1)
        b = cauchy_offcut_row(M, 0)
        alt = infinity_constant_row(M)
        rows.append({"M": M,
                     "constant_row_max_difference": float(np.max(abs(b - alt))),
                     "stable_tan_max_difference": float(np.max(abs(alt - np.tan(p / 2) / M))),
                     "sine_mode_max_error": float(np.max(abs(
                         alt @ np.sin(np.outer(p, n)) + (-1.) ** n))),
                     "oncut_matrix_max_difference": float(np.max(abs(
                         cauchy_oncut_real(M) - cauchy_oncut_real(M, "infinity")))),
                     "exact_f00_supremum_difference": "0 (operator identity)"})
    doc = {"audit_item": "A1", "provenance": provenance(), "rows": rows,
           "conclusion": "Exactly identical constants, including the Nyquist mode.",
           "optimization_needed": False,
           "scope": "Only a0 is replaced; all offcut/crossed rows remain unchanged.",
           "asymptotic_order": "Identically zero for every M; no nonzero convergence order.",
           "fourth_caliber_needed": False}
    here = Path(__file__).resolve().parent
    doc["audit_code_sha256"] = {name: hashlib.sha256((here / name).read_bytes()).hexdigest()
                                 for name in ("audit.py", "hilbert.py", "grid.py")}
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    doc["producer"] = save_producer(out)
    with (out / "report.json").open("x") as stream:
        json.dump(doc, stream, indent=2)
    return doc


def phase_replay(source_report, outdir):
    """A3: re-read full saved amplitudes, preserving their feasibility status."""
    from .assembly import Operators
    from .observables import phase_shift, interpolate_at
    from .verify import unitarity_report, chiral_report
    source = Path(source_report)
    raw = json.loads(source.read_text())
    hashes = {str(source): hashlib.sha256(source.read_bytes()).hexdigest()}
    rows = []
    for rec in raw["records"]:
        spec = rec["spec"]
        ops = Operators(spec["M"], spec["L"])
        path = source.parent / f"sol_{rec['job']}.npz"
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        c = np.load(path)["c"]
        uni = unitarity_report(ops, c)
        chi = chiral_report(ops, c, spec["chi_caliber"], spec["eps_chi"])
        waves = {w: phase_shift(ops, c, w) for w in ("S0", "S2", "P1")}
        d0 = interpolate_at(waves["S0"], .9)
        d1 = interpolate_at(waves["P1"], 1.2)
        rows.append({"job": rec["job"], "M": spec["M"], "L": spec["L"],
                     "source_certified": rec["result"].get("certified", False),
                     "unitarity": uni, "chiral": chi,
                     "delta00_at_0p9GeV": d0, "delta11_at_1p2GeV": d1,
                     "C4_numerical_window_pass": 85 <= d0 <= 110 and d1 <= 25,
                     "C4_certified": False,  # this replay supplies no support bound
                     "S_reconstruction_max_error": max(float(np.nanmax(abs(
                         ph["eta"] * np.exp(2j * np.deg2rad(ph["delta_deg"])) - ph["S"])))
                         for ph in waves.values()),
                     "waves": {w: {"E_GeV": ph["E_GeV"].tolist(),
                                    "delta_deg": ph["delta_deg"].tolist(),
                                    "eta": ph["eta"].tolist(),
                                    "first_zero_node": ph["first_zero_node"]}
                               for w, ph in waves.items()}})
    doc = {"audit_item": "A3", "provenance": provenance(), "input_sha256": hashes,
           "rows": rows, "phase_convention": waves["S0"]["phase_convention"],
           "scope": "Saved C only; rejected source amplitudes remain diagnostic."}
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    doc["producer"] = save_producer(out)
    with (out / "report.json").open("x") as stream:
        json.dump(doc, stream, indent=2)
    return {"rows": [{k: v for k, v in row.items() if k != "waves"} for row in rows]}


def current_witness(M, outdir, solver="CLARABEL"):
    """A4: artificial strict current-sector witness, not a pion amplitude.

    For S=alpha F/F*, Schur positivity is rho_hat >= 2|F|^2/(1+alpha).
    A slightly stronger quadratic bound constructs a strictly positive Gram.
    No scattering operator, chiral claim or full-problem Slater claim follows.
    """
    import cvxpy as cp
    from flint import acb
    from . import constraints as C, formfactor as FF
    from .grid import s_nodes
    from .hilbert import hilbert_kernel
    from .arbaudit import ArbAudit
    from .problem import solver_options
    alpha, beta = .9, 1.1
    im = cp.Variable((2, M))
    rh = cp.Variable((2, M), nonneg=True)
    K = hilbert_kernel(M)
    idx, bounds, factors = C.ff_asymptotic_bounds(M)
    cons = []
    targets = C.printed_targets()
    for ell, wave in ((0, "S0"), (1, "P1")):
        re = 1 + K @ im[ell]
        cons.append(beta * (cp.square(re) + cp.square(im[ell])) <= rh[ell])
        k2 = FF.kinematic_square_reference(ell, s_nodes(M))
        for n in C.MOMENTS[ell]:
            cons.append((C.moment_row(M, n) * k2) @ rh[ell] / targets[(wave, n)] == 1)
        cons.append(cp.SOC(.9 * bounds[ell] / factors[ell],
                           cp.vstack([re[idx], im[ell, idx]]), axis=0))
    problem = cp.Problem(cp.Minimize(cp.sum(rh)), cons)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "report.json").exists():
        raise FileExistsError(out / "report.json")
    doc = {"audit_item": "A4", "provenance": provenance(), "M": M,
           "scope": "Artificial current-sector witness; S is not from a common scattering C.",
           "alpha": alpha, "beta": beta, "ff_construction_fraction": .9,
           "solver": solver,
           "solver_options": solver_options(solver, max_iter=500, tol_gap_abs=1e-10,
                              tol_gap_rel=1e-10, tol_feas=1e-10, time_limit=3600.0)}
    try:
        problem.solve(solver=solver, **doc["solver_options"])
        doc["solver_status"] = problem.status
        doc["iterations"] = problem.solver_stats.num_iters
        doc["solve_time"] = problem.solver_stats.solve_time
        if im.value is not None:
            imv, rhv = np.asarray(im.value), np.asarray(rh.value)
            F = 1 + imv @ K.T + 1j * imv
            S = np.divide(alpha * F, F.conj(), out=np.zeros_like(F), where=F != 0)
            path = out / "current_witness.npz"
            np.savez_compressed(path, ImF=imv, rho_hat=rhv, S=S)
            doc["solution_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            store = {(ell, i): acb(S[ell, i].real, S[ell, i].imag)
                     for ell in (0, 1) for i in range(M)}
            audit = ArbAudit(M, 10, bits=384)
            doc["arb_audit"] = audit._uv_audit(store, imv, rhv, "SR-b", C.EPS_FF, C.M_Q)
    except Exception as exc:
        doc["error"] = str(exc)
    doc["producer"] = save_producer(out)
    with (out / "report.json").open("x") as stream:
        json.dump(doc, stream, indent=2)
    return {k: v for k, v in doc.items() if k != "arb_audit"}
