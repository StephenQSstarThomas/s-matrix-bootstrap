"""One 2403 leaf: build the program, write the PMP, run SDPB, read back, report.

Every leaf directory keeps pmp.json (deleted on request after acceptance), the SDPB output, solution.npz (free
variables as decimal strings and the lifted raw variables) and report.json with the spec, the settings, the solver
result and the observables.  Verification (V1/V2) is a separate step in certify.py that reads this directory.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from . import kernels, pmpwriter
from .program import Program

SOURCE_FILES = ("program.py", "regulariser.py", "pmpwriter.py", "solve.py", "kernels.py", "spec.py")


def amplitude_basis(prog_rows, tol=1e-12):
    """Right singular vectors of the row-normalised amplitude rows with s > tol * s_max (an inner approximation)."""
    R = np.vstack(prog_rows)
    norms = np.abs(R).max(axis=1)
    R = R[norms > 0] / norms[norms > 0, None]
    _, s, vt = np.linalg.svd(R, full_matrices=False)
    k = int((s > tol * s[0]).sum())
    return vt[:k].T, {"rank": k, "n_rows": int(R.shape[0]), "tol": tol, "s_max": float(s[0]), "s_min_kept": float(s[k - 1])}


def watson_targets_from_leaf(leaf, ks):
    """Targets from any previous leaf (initial or Watson): only its rows and lifted solution are used."""
    from .spec import GtbSpec
    from .watson import targets_from_solution
    prev = json.loads((Path(leaf) / "report.json").read_text())
    sp = prev["spec"]
    pprog = Program(GtbSpec(M=sp["M"], L=sp.get("L", 20), n_b=sp.get("n_b", 40)), ks, reg="none")
    return targets_from_solution(pprog, np.load(Path(leaf) / "solution.npz")["v"])


def _hot_start_seed(seed, prog, reg, reduce_basis, basis, work):
    """Copy a neighbouring leaf's final checkpoint to work/init.ck after checking that the SDP has its structure."""
    import shutil
    s = json.loads((Path(seed) / "report.json").read_text())
    if ((s.get("n_vars"), s.get("reduce_basis"), s.get("regulariser")) != (prog.n_vars, reduce_basis, reg)
            or not (Path(seed) / "sdp.ck").is_dir()
            or basis is not None and not np.array_equal(np.load(Path(seed) / "basis.npy"), basis)):
        raise ValueError(f"hot start from {seed}: the seed's SDP has another structure or no checkpoint")
    shutil.copytree(Path(seed) / "sdp.ck", work / "init.ck", dirs_exist_ok=True)
    return "init.ck"


def _observables(prog, v):
    amp, rho = v[prog.groups["amp"]], v[prog.groups["rho"]]
    S = 1 + 1j * (prog.Rre @ amp + 1j * (prog.Rim @ amp))                # S = 1 + i h at every node and wave
    eta = {w: float(np.abs(S[[prog.row_index(w, j) for j in range(prog.n0)]]).min()) for w in ("S0", "P1", "D0", "S2")}
    return {"F0": float(prog.F0row @ amp), "F1": float(prog.F1row @ amp),
            "rho_linf": float(np.abs(rho).max()), "rho_l4": float((rho ** 4).sum() ** 0.25),
            "rho_l2": float(np.sqrt((rho ** 2).sum())), **{f"min_eta_{w}": e for w, e in eta.items()}}


def run_leaf(spec, kernels_path, workdir, *, reg="linf-relax", reduce_basis=False, basis_tol=1e-12,
             precision=192, nproc=8, timeout=7200, duality_gap=1e-6, primal_error=1e-10, dual_error=1e-10,
             initial_scale=1e4, max_iterations=500, keep_pmp=True, watson_from=None, objective_row=None,
             hot_start=None):
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ks = kernels.load(kernels_path)
    basis, basis_info = None, None
    if reduce_basis:
        probe = Program(spec, ks, reg="none")
        lam = ks.Lam[:, None]
        basis, basis_info = amplitude_basis([probe.Rre / lam, probe.Rim / lam, probe.Rim / lam ** 2, probe.Fchi,
                                             probe.F0row[None], probe.F1row[None]], basis_tol)
    targets = None
    if watson_from is not None:                      # Watson step seeded by a previous leaf (same kernels)
        from dataclasses import replace
        prev = json.loads((Path(watson_from) / "report.json").read_text())
        targets = watson_targets_from_leaf(watson_from, ks)
        np.savez(work / "watson_targets.npz", **{f"{w}_{p}": t[i] for w, t in targets.items()
                                                  for i, p in enumerate(("re", "im"))})
        spec = replace(spec, objective="watson", f0_ni=prev["spec"]["f0_ni"])
    if objective_row is not None:
        np.save(work / "objective_row.npy", np.asarray(objective_row, dtype=float))
    prog = Program(spec, ks, reg=reg, basis=basis, watson_targets=targets, objective_row=objective_row)
    if basis is not None:
        np.save(work / "basis.npy", basis)
    rep = {"spec": asdict(spec), "spec_key": spec.key(), "kernels": str(kernels_path),
           "kernels_manifest": ks.manifest.get("generator"), "regulariser": reg, "reduce_basis": reduce_basis,
           "basis": basis_info, "f0fix": prog.f0fix, "n_vars": prog.n_vars, "pivot": list(prog.pivot) if prog.pivot else None,
           "watson_from": str(watson_from) if watson_from else None,
           "c_dropped": [int(i) for i in prog.c_dropped],
           "settings": {"precision": precision, "nproc": nproc, "timeout": timeout, "duality_gap": duality_gap,
                        "primal_error": primal_error, "dual_error": dual_error, "initial_scale": initial_scale},
           "source_sha256": {f: hashlib.sha256((Path(__file__).parent / f).read_bytes()).hexdigest()
                             for f in SOURCE_FILES}}
    init = _hot_start_seed(hot_start, prog, reg, reduce_basis, basis, work) if hot_start is not None else None
    rep["hot_start"] = str(hot_start) if hot_start is not None else None
    rep["pmp"] = pmpwriter.write(prog.blocks(), prog.objective(), work / "pmp.json", prog.n_vars)
    if init and json.loads((Path(hot_start) / "report.json").read_text())["pmp"]["n_blocks"] != rep["pmp"]["n_blocks"]:
        raise ValueError(f"hot start from {hot_start}: different block count")
    rep["seconds_build"] = time.time() - t0
    res = pmpwriter.solve(work, precision=precision, nproc=nproc, timeout=timeout, duality_gap=duality_gap,
                          primal_error=primal_error, dual_error=dual_error, initial_scale=initial_scale,
                          max_iterations=max_iterations, initial_checkpoint=init)
    rep["sdpb_result"] = {k: v for k, v in res.items() if k not in ("pmp2sdp", "sdpb")}
    rep["pmp2sdp"], rep["sdpb"] = res.get("pmp2sdp"), res.get("sdpb")
    if (work / "out/y.txt").exists():
        y = pmpwriter.read_y(work / "out", prog.n_vars)
        v = prog.lift(np.array(y["float"]))
        np.savez_compressed(work / "solution.npz", y_text=np.array(y["text"]), v=v)
        rep["observables"] = _observables(prog, v)
    rep["seconds_total"] = time.time() - t0
    if not keep_pmp and rep["sdpb_result"].get("terminateReason") == "found primal-dual optimal solution":
        (work / "pmp.json").unlink(missing_ok=True)
    (work / "report.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")
    return rep
