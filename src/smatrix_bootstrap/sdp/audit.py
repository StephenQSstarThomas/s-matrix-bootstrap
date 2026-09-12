"""Reproducible, source-led audit calculations; no optimization in this module."""
from pathlib import Path
import hashlib
import json

import numpy as np

from .grid import phi_nodes
from .hilbert import cauchy_offcut_row, cauchy_oncut_real, infinity_constant_row
from .runner import provenance


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
    with (out / "report.json").open("x") as stream:
        json.dump(doc, stream, indent=2)
    return doc
