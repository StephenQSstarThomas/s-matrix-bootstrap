"""Degree-0 SDPB PMP writer and solution reader for the 2403 program.

SDPB maximises objective . y subject to normalization . y = 1 and sum_n y_n M_j^(n) >= 0 for every block j.  With
normalization = (1, 0, ..., 0) the variable y0 = 1 carries the constant terms.  Every polynomial has degree 0, i.e.
the constraint is imposed at the collocation data only, exactly as in the finite model.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

_PRE = '{"DampedRational":{"base":"1","constant":"1","poles":[]}'


def num(x, digits: int = 40) -> str:
    """One coefficient as a decimal string; exact zero is "0"."""
    if isinstance(x, str):
        return x
    if hasattr(x, "mid") and hasattr(x, "is_finite"):          # python-flint arb
        if not x.is_finite():
            raise ValueError("nonfinite PMP coefficient")
        return "0" if x == 0 else x.mid().str(digits, radius=False)
    if hasattr(x, "_mpf_"):                                      # mpmath mpf
        from mpmath import nstr
        return "0" if x == 0 else nstr(x, digits, min_fixed=1, max_fixed=0)
    if x == 0:
        return "0"
    return repr(float(x))


def _cell(form, n_vars, digits):
    if form is None:
        return "[" + ",".join(['["0"]'] * n_vars) + "]"
    if len(form) != n_vars:
        raise ValueError(f"form has {len(form)} coefficients, expected {n_vars}")
    if isinstance(form, np.ndarray) and form.dtype == np.float64:
        # repr of a float64 is the shortest string that round-trips; exact zeros print as "0"
        return '[["' + '"],["'.join(["0" if v == 0.0 else repr(v) for v in form.tolist()]) + '"]]'
    return "[" + ",".join('["%s"]' % num(v, digits) for v in form) + "]"


def write(blocks, objective, path, n_vars: int, digits: int = 40) -> dict:
    """Stream the PMP JSON for ``blocks`` (list of program.Block) and ``objective`` (length n_vars)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(objective) != n_vars:
        raise ValueError("objective length differs from n_vars")
    sizes = {}
    with path.open("w") as fh:
        fh.write('{"objective":' + json.dumps([num(v, digits) for v in objective]))
        fh.write(',"normalization":' + json.dumps(["1"] + ["0"] * (n_vars - 1)))
        fh.write(',"PositiveMatrixWithPrefactorArray":[')
        for k, blk in enumerate(blocks):
            fh.write(("," if k else "") + _PRE + ',"polynomials":[')
            for r in range(blk.size):
                fh.write(("," if r else "") + "[" + ",".join(
                    _cell(blk.entry(r, c), n_vars, digits) for c in range(blk.size)) + "]")
            fh.write("]}")
            sizes[blk.size] = sizes.get(blk.size, 0) + 1
        fh.write("]}")
    return {"path": str(path), "bytes": os.path.getsize(path), "n_vars": n_vars,
            "n_blocks": len(blocks), "block_sizes": sizes, "digits": digits}


class _Cfg:
    def __init__(self, nproc, timeout):
        self.nproc, self.timeout = nproc, timeout


def sdpb_args(*, precision, timeout, duality_gap, primal_error, dual_error, initial_scale, max_iterations,
              checkpoint_interval, write_solution, initial_checkpoint=None) -> list:
    """The sdpb command line.  A hot start loads x, X, y, Y from ``initial_checkpoint`` (a neighbouring problem with
    the same block structure) instead of the scaled identity; checkpoints still go to sdp.ck."""
    args = ["sdpb", f"--precision={precision}", "-s", "sdp", "-o", "out", f"--maxIterations={max_iterations}",
            f"--maxRuntime={timeout}", f"--checkpointInterval={checkpoint_interval}",
            f"--dualityGapThreshold={duality_gap}", f"--primalErrorThreshold={primal_error}",
            f"--dualErrorThreshold={dual_error}", f"--writeSolution={write_solution}"]
    if initial_checkpoint:
        return args + ["-i", str(initial_checkpoint), "-c", "sdp.ck"]
    return args + [f"--initialMatrixScalePrimal={initial_scale}", f"--initialMatrixScaleDual={initial_scale}"]


def solve(workdir, precision: int = 192, nproc: int = 8, timeout: int = 7200, duality_gap: float = 1e-6,
          primal_error: float = 1e-10, dual_error: float = 1e-10, initial_scale: float = 1e4,
          max_iterations: int = 500, checkpoint_interval: int = 300, write_solution: str = "x,y,X,Y",
          initial_checkpoint=None) -> dict:
    """pmp2sdp + sdpb in ``workdir`` (which holds pmp.json) through the 2309 launcher; returns out/out.txt.
    ``initial_checkpoint`` is a directory name inside ``workdir``."""
    from ..sdp.sdpb import execute, read_out
    work, cfg = Path(workdir), _Cfg(nproc, timeout)
    pre = execute(work, "pmp2sdp", ["pmp2sdp", f"--precision={precision}", "--input=pmp.json", "--output=sdp"], cfg)
    if pre["returncode"]:
        return {"terminateReason": "pmp2sdp failed", "pmp2sdp": pre}
    args = sdpb_args(precision=precision, timeout=timeout, duality_gap=duality_gap, primal_error=primal_error,
                     dual_error=dual_error, initial_scale=initial_scale, max_iterations=max_iterations,
                     checkpoint_interval=checkpoint_interval, write_solution=write_solution,
                     initial_checkpoint=initial_checkpoint)
    run = execute(work, "sdpb", args, cfg, mpi=nproc > 1)
    out = read_out(work / "out/out.txt")
    out.update({"pmp2sdp": pre, "sdpb": run})
    return out


def read_y(outdir, n_vars: int) -> dict:
    """SDPB's y.txt holds y1..y_{n_vars-1} (y0 is eliminated by the normalisation)."""
    lines = Path(outdir, "y.txt").read_text().split("\n")
    head = lines[0].split()
    vals = [ln.strip() for ln in lines[1:] if ln.strip()]
    if int(head[0]) != n_vars - 1 or len(vals) != n_vars - 1:
        raise ValueError(f"y.txt has {head[0]} rows / {len(vals)} values, expected {n_vars - 1}")
    return {"text": vals, "float": [float(v) for v in vals]}
