"""Shadow prices of the 2403 tolerances at an optimal leaf (attribution only; nothing is tuned with them).

Envelope theorem for  opt(theta) = max_y b_theta.[1, y]  s.t.  M_j,theta(y) >= 0  with SDPB's multipliers Z_j >= 0
(the same Z as in dual.py: b_n + sum_j <Z_j, A_jn> = 0):
    d opt / d theta = d_theta b_theta.[1, y*] + sum_j <Z_j, d_theta M_j,theta(y*)> .
The theta-derivative of every block is a central difference of the rebuilt program (the forms are smooth in each
tolerance), at the leaf's own y* and Z.  Elasticities theta d opt/d theta are reported as well.
"""
from __future__ import annotations

import json
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

PARAMETERS = {"eps_chi": None, "eps_sr": 3, "r": 3, "reg_bound": None, "f0_max": None}


def _spec(rep):
    from .spec import GtbSpec
    return GtbSpec(**{f.name: (tuple(rep["spec"][f.name]) if isinstance(rep["spec"][f.name], list)
                               else rep["spec"][f.name]) for f in fields(GtbSpec)})


def _values(prog, yfull):
    """Every block evaluated at [1, y] as a dense symmetric matrix, and the objective value."""
    out = []
    for blk in prog.blocks():
        m = np.zeros((blk.size, blk.size))
        for (r, c), f in blk.forms.items():
            m[r, c] = m[c, r] = float(np.asarray(f, dtype=float) @ yfull)
        out.append(m)
    return out, float(np.asarray(prog.objective(), dtype=float) @ yfull)


def _perturbed(spec, name, k, factor):
    v = getattr(spec, name)
    if k is None:
        return replace(spec, **{name: v * factor})
    v = list(v)
    v[k] *= factor
    return replace(spec, **{name: tuple(v)})


def shadow_prices(leaf, kernels_dir, h: float = 1e-4, which=None) -> dict:
    from .dual import read_dual_matrices
    from .kernels import load
    from .pmpwriter import read_y
    from .program import Program
    leaf = Path(leaf)
    rep = json.loads((leaf / "report.json").read_text())
    spec, ks = _spec(rep), load(kernels_dir)
    basis = np.load(leaf / "basis.npy") if rep.get("reduce_basis") else None
    targets = None
    if (leaf / "watson_targets.npz").exists():
        z = np.load(leaf / "watson_targets.npz")
        targets = {w: (z[f"{w}_re"], z[f"{w}_im"]) for w in {k.rsplit("_", 1)[0] for k in z.files}}
    orow = np.load(leaf / "objective_row.npy") if (leaf / "objective_row.npy").exists() else None
    extra = None
    if rep.get("midpoints") and rep["midpoints"].get("rows"):        # WELLPOSED W3 cones are part of the program
        from .midpoints import load as mload
        extra = mload(rep["midpoints"]["kernels"], rep["midpoints"]["rows"])
    build = lambda sp: Program(sp, ks, reg=rep["regulariser"], basis=basis, watson_targets=targets, objective_row=orow,
                               extra_unitarity=extra)
    base = build(spec)
    yfull = np.r_[1.0, read_y(leaf / "out", base.n_vars)["float"]]
    Z = [np.array(z, dtype=float) for z in read_dual_matrices(leaf / "out", [b.size for b in base.blocks()])]
    out = {}
    for name, n in PARAMETERS.items():
        for k in ([None] if n is None else range(n)):
            label = name if k is None else f"{name}[{k}]"
            theta = getattr(spec, name) if k is None else getattr(spec, name)[k]
            if (which and label not in which) or not np.isfinite(theta):
                continue
            (mp, op), (mm, om) = (_values(build(_perturbed(spec, name, k, 1 + s * h)), yfull) for s in (1, -1))
            d = ((op - om) + sum(float(np.sum(z * (a - b))) for z, a, b in zip(Z, mp, mm))) / (2 * h * theta)
            out[label] = {"theta": theta, "dopt_dtheta": d, "elasticity": theta * d}
    return out
