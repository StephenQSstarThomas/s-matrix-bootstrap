"""Write the paper's cone program as an SDPB polynomial matrix program.

SDPB solves ``maximise b.y`` subject to ``M_j(y) >= 0`` for each block ``j``,
plus one normalisation ``n.y = 1``.  A PMP has no separate constant matrix, so
the normalisation variable ``y0`` (pinned to 1) carries every constant term.

Everything here is **degree-0**: each polynomial is a single coefficient.  That
is deliberate and load-bearing.  A non-zero degree would impose the constraint
for all ``x >= 0``, i.e. unitarity on a *continuum* of energies, which is
strictly stronger than the paper -- the paper imposes it on the M collocation
nodes only (task spec, ironclad rule 3).

Cone translations (each verified against MOSEK on the same ModelSpec by
``scripts/sdp/pmp_crosscheck.py``)::

    |p|^2 <= 2 r            ->  [[2r, p_re, p_im], [p_re, y0, 0], [p_im, 0, y0]] >= 0
    ||v||_2 <= e            ->  arrow  [[e y0, v^T], [v, e y0 I]] >= 0
    a.y <= c                ->  1x1    [c y0 - a.y]
    3x3 Gram (3.68)         ->  the real congruent form of problem._gram_psd

Precision: the main CLI rebuilds source geometry and projects rows in Arb
before emitting 30-digit coefficients. ``operator_dps=17`` is retained only
for historical matrix replay; printing more digits cannot improve those
float64 rows. The SVD coordinate basis remains an explicit numerical choice.
"""
from __future__ import annotations

import json
import os

import numpy as np
from flint import arb, ctx

from . import constraints as C
from . import formfactor as FFM
from .assembly import Operators, load_saved_basis
from .spec import ModelSpec

_PRE = '{"DampedRational":{"base":"1","constant":"1","poles":[]}'

# Secondary linear functionals for the degenerate-face diagnostic.  Each is a
# single node value that is exactly linear in the PMP variables.
FUNCTIONAL_KINDS = {
    "ImKH": "Im(kappa h) = 1 - Re S at one node",
    "ImS": "Im S = Re(kappa h) at one node",
    "ImF": "Im F_ell at one node (UV-sector primal variable)",
    "rho": "rho_hat_ell at one node (UV-sector primal variable)",
    "SRmom": "one FESR moment (3.73) of wave S0/P1; 'node' carries n (S0: 0,1; P1: -1,0)",
    "watson": "linearised unitarity-saturation step (see watson.py): sum over nodes below s0 of Re(t* h) - Im h",
}


def free_moments(spec):
    return {(str(w), int(n)) for w, n in (spec.sr_free or ())}


def check_face(face):
    """``d.(f00, f11) >= value - margin``: the near-optimal slab of a completed support."""
    if face is None:
        return None
    d = np.asarray(face["direction"], dtype=float)
    value, margin = float(face["value"]), float(face["margin"])
    if d.shape != (2,) or not np.all(np.isfinite(d)) or not np.any(d != 0):
        raise ValueError("face direction must be a nonzero finite pair")
    if not (np.isfinite(value) and np.isfinite(margin) and margin > 0):
        raise ValueError("face margin must be a positive finite number")
    return {"direction": d.tolist(), "value": value, "margin": margin}


def check_functional(functional, spec):
    """Validate ``{kind, wave, node, sense}`` against the model; adds ``ell``."""
    if functional is None:
        return None
    kind, wave, node, sense = (functional[k] for k in ("kind", "wave", "node", "sense"))
    if kind not in FUNCTIONAL_KINDS:
        raise ValueError(f"Unknown functional kind {kind!r}; known: {sorted(FUNCTIONAL_KINDS)}")
    if kind == "watson":
        if sense != "max":
            raise ValueError("the Watson step is a maximisation")
        if spec.scattering_prescription != "mixed-pv":
            raise ValueError("node functionals are defined for the nodal (mixed-pv) prescription")
        waves, nodes, targets = functional.get("waves"), functional.get("nodes"), functional.get("targets")
        if not waves or any(w not in ("S0", "P1", "S2") for w in waves):
            raise ValueError("watson waves must be a non-empty subset of S0, P1, S2")
        if not nodes or any((not isinstance(k, (int, np.integer))) or not 0 <= k < spec.M for k in nodes):
            raise ValueError("watson nodes must be node indices of this model")
        if not isinstance(targets, dict) or any(len(targets.get(w, [])) != len(nodes) for w in waves):
            raise ValueError("watson targets must give one (re, im) pair per node for every listed wave")
        for w in waves:
            for pair in targets[w]:
                if len(pair) != 2 or not all(np.isfinite(float(v)) for v in pair):
                    raise ValueError("watson targets must be finite (re, im) pairs")
        weights = functional.get("weights")
        if weights is not None:
            if any(len(weights.get(w, [])) != len(nodes) or any((not np.isfinite(float(v))) or float(v) <= 0 for v in weights[w]) for w in waves):
                raise ValueError("watson weights must be positive finite numbers, one per node and wave")
        out = dict(functional); out.update(kind="watson", sense="max", wave="all", node=0, ell=None,
                                          waves=list(waves), nodes=[int(k) for k in nodes])
        return out
    if wave not in ("S0", "P1"):
        raise ValueError("functional wave must be S0 or P1")
    ell = 0 if wave == "S0" else 1
    if isinstance(node, bool) or not isinstance(node, (int, np.integer)):
        raise ValueError("functional node must be an integer")
    if kind == "SRmom":
        if not (spec.uv and "fesr" in spec.uv_parts):
            raise ValueError("SRmom needs the UV sector with the FESR part")
        if node not in C.MOMENTS[ell]:
            raise ValueError(f"SRmom n must be one of {C.MOMENTS[ell]} for {wave}")
    elif not 0 <= node < spec.M:
        raise ValueError(f"functional node must be an integer in [0, {spec.M})")
    if sense not in ("max", "min"):
        raise ValueError("functional sense must be max or min")
    if kind in ("ImF", "rho") and not spec.uv:
        raise ValueError(f"{kind} is a UV-sector variable; this model has no UV sector")
    if kind in ("ImKH", "ImS") and spec.scattering_prescription != "mixed-pv":
        raise ValueError("node functionals are defined for the nodal (mixed-pv) prescription")
    return {"kind": kind, "wave": wave, "node": int(node), "sense": sense,
            "ell": 0 if wave == "S0" else 1}


def _num(x: float, digits: int = 17) -> str:
    """Coefficient as a JSON string.  Exact zero is printed as "0"."""
    if x == 0.0:
        return "0"
    if isinstance(x, arb):
        if not x.is_finite():
            raise ValueError("Nonfinite PMP coefficient")
        return x.mid().str(digits, radius=False)
    return repr(float(x)) if digits >= 17 else f"{float(x):.{digits}e}"


class Pmp:
    """Assemble and stream one PMP file for a ``ModelSpec`` and a direction."""

    @classmethod
    def from_saved(cls, source_report, direction, fix_f00=None, face=None, functional=None):
        from .precision_pmp import restore
        return restore(cls,source_report,direction,fix_f00,face,functional)

    def __init__(self, spec: ModelSpec, direction=(1.0, 0.0), fix_f00=None,
                 digits: int = 17, basis_source_report=None, face=None, functional=None) -> None:
        if spec.B is not None and spec.B_norm == "l4":
            raise ValueError("l4 density ball is not expressible as a degree-0 PMP "
                             "block of reasonable size; use B_norm='l2' or B=None")
        if spec.reg_norm is not None:
            if spec.reg_norm not in ("linf", "l2", "l4"):
                raise ValueError(f"Unsupported regulariser norm: {spec.reg_norm}")
            if (spec.reg_bound is None or not np.isfinite(spec.reg_bound)
                    or spec.reg_bound <= 0):
                raise ValueError("reg_bound must be a positive finite Mreg")
        elif spec.reg_bound is not None:
            raise ValueError("reg_bound given without reg_norm")
        self.spec, self.direction, self.fix_f00, self.digits = spec, direction, fix_f00, digits
        self.face, self.functional = check_face(face), check_functional(functional, spec)
        saved=load_saved_basis(spec,basis_source_report) if basis_source_report is not None else None
        self.basis_source=None if saved is None else saved[1]
        if spec.scattering_prescription=='sine-cardinal' and spec.operator_dps<=17:
            raise ValueError('Analytic source assembly requires arbitrary-precision coefficients')
        self.ops = ops = Operators(spec.M, spec.L, spec.scattering_prescription, spec.operator_dps)
        ops.set_cone_scaling(spec.cone_scaling, spec.sparsify)
        M = spec.M

        # which unitarity disks are imposed (mirrors problem.Model)
        keep = (np.ones(len(ops.index) * M, dtype=bool) if spec.disk_mask is None
                else np.asarray(spec.disk_mask, dtype=bool).copy())
        if spec.uv and "gram" in spec.uv_parts:
            for ell, I in ((0, 0), (1, 1)):
                a0 = ops.index.index((I, ell))
                keep[a0 * M:(a0 + 1) * M] = False
        self.keep = keep

        # the basis spans only the imposed rows: that is where the cubic saving is
        gram=spec.uv and 'gram' in spec.uv_parts
        V=(ops.adopt_basis(saved[0],mask=keep,include_gram=gram) if saved is not None else
           ops.build_basis(spec.basis_tol,mask=keep,include_gram=gram) if spec.reduce_basis else None)
        self.basis = V
        self.precise = None
        if spec.operator_dps > 17:
            if digits > spec.operator_dps-5:
                raise ValueError("Output digits must leave at least five source guard digits")
            from .precision_pmp import populate
            populate(self)
        else:
            proj = (lambda A: np.asarray(A) @ V) if V is not None else (lambda A: np.asarray(A))
            self.P_re, self.P_im, self.R_im = proj(ops.P_re), proj(ops.P_im), proj(ops.R_im)
            self.f00 = proj(ops.f_proj["f00"][None, :])[0]
            self.f11 = proj(ops.f_proj["f11"][None, :])[0]
            self.chi = proj(ops.chi_rows)
            self.gram = {e: (proj(a[None, :] if a.ndim == 1 else a), proj(b))
                         for e, (a, b) in ops.gram_rows.items()}

        # ---- variable layout:  y0 | a | ImF(2M) | rho_hat(2M)
        self.n_a = self.P_re.shape[1]
        self.i_a = 1
        if spec.uv:
            self.i_ImF = self.i_a + self.n_a
            self.i_rho = self.i_ImF + 2 * M
            self.n_vars = self.i_rho + 2 * M
        else:
            self.i_ImF = self.i_rho = None
            self.n_vars = self.i_a + self.n_a
        # ---- auxiliary variables of the l2/l4 regulariser, appended last so
        # every index above is unchanged: t_i >= r_i^2 (l2, l4) and u_i >= t_i^2 (l4)
        lay = self.ops.lay
        self.n_rho = (lay.r1.stop - lay.r1.start) + (lay.r2.stop - lay.r2.start)
        self.n_aux = {None: 0, "linf": 0, "l2": self.n_rho, "l4": 2 * self.n_rho}[spec.reg_norm]
        self.i_aux = self.n_vars
        self.n_vars += self.n_aux

        self.n_blocks = 0

    # ------------------------------------------------------------------ rows
    def _row(self, **parts) -> np.ndarray:
        """Build a length-``n_vars`` coefficient row from named pieces."""
        v = np.zeros(self.n_vars, dtype=object if self.precise else float)
        if "y0" in parts:
            x = parts["y0"]
            v[0] = arb(str(x)) if self.precise and not isinstance(x,arb) else x
        if "a" in parts:
            v[self.i_a:self.i_a + self.n_a] = parts["a"]
        if "ImF" in parts:
            ell, vec = parts["ImF"]
            off = self.i_ImF + ell * self.spec.M
            v[off:off + self.spec.M] = vec
        if "rho" in parts:
            ell, vec = parts["rho"]
            off = self.i_rho + ell * self.spec.M
            v[off:off + self.spec.M] = vec
        if "aux" in parts:
            index, value = parts["aux"]
            v[self.i_aux + index] = value
        return v

    # ---------------------------------------------------------------- blocks
    def _blocks(self):
        """Yield each PSD block as a list-of-lists of coefficient rows."""
        spec, M = self.spec, self.spec.M
        if self.precise:
            ctx.prec = self.precise.bits
        y0 = self._row(y0=1.0)
        zero = np.zeros(self.n_vars)

        # ---- unitarity disks:  |p|^2 <= 2 r
        for k in np.flatnonzero(self.keep):
            p_re = self._row(a=self.P_re[k])
            p_im = self._row(a=self.P_im[k])
            two_r = self._row(a=2.0 * self.R_im[k])
            yield [[two_r, p_re, p_im], [p_re, y0, zero], [p_im, zero, y0]]

        # ---- chiral (3.64)
        if spec.chiral:
            e = spec.eps_chi
            res = [self._row(a=self.chi[k]) for k in range(self.chi.shape[0])]
            if spec.chi_caliber == "chi-a":
                for r in res:
                    yield [[self._row(y0=e) - r]]
                    yield [[self._row(y0=e) + r]]
            elif spec.chi_caliber == "chi-b":
                yield self._arrow(res, e)
            elif spec.chi_caliber == "chi-c":
                yield self._arrow(res[0::2], e)
                yield self._arrow(res[1::2], e)
            else:
                raise ValueError(spec.chi_caliber)

        # ---- double-spectral-density regulariser (2103.11484 section 3)
        if spec.reg_norm is not None:
            yield from self._regulariser_blocks()

        # ---- density ball (l2 only; see __init__)
        if spec.B is not None:
            lay = self.ops.lay
            if self.basis is not None:
                Vr = np.vstack([self.basis[lay.r1], self.basis[lay.r2]])
                rows = [self._row(a=Vr[i]) for i in range(Vr.shape[0])]
            else:
                idx = list(range(lay.r1.start, lay.r1.stop)) + list(range(lay.r2.start, lay.r2.stop))
                rows = []
                for i in idx:
                    u = np.zeros(self.n_a); u[i] = 1.0
                    rows.append(self._row(a=u))
            yield self._arrow(rows, spec.B)

        # ---- UV sector
        if spec.uv:
            if self.precise:
                K,idx_hi,ffb,ffk,tgt,tol,moments = self.precise.uv_data(spec)
                if spec.sr_caliber in ("SR-a", "SR-d"):      # adjustable eps_SR; the operator module keeps the 2309 value
                    tol = {k: arb(str(spec.eps_sr)) for k in tol}
                rt2 = arb(2).sqrt()
            else:
                K = FFM.hilbert_kernel(M)
                idx_hi, ffb, ffk = C.ff_asymptotic_bounds(M, spec.m_q, spec.eps_ff, spec.ff_frozen_at_s0)
                tgt, tol = C.printed_targets(), C.sr_tolerances(spec.sr_caliber, spec.eps_sr)
                rt2 = np.sqrt(2.0)
            for ell in (0, 1):
                g = np.ones(M) if self.precise else FFM.gram_scale(ell, self.ops.s)
                kin = g if self.precise else FFM.kinematic_factor(ell, self.ops.s)
                re_row, im_row = self.gram[ell]
                re_row = re_row if re_row.ndim == 2 else re_row[None, :]
                for i in range(M):
                    # S = 1 - (im_row.a) + i (re_row.a);  F = (kin/g)(ReF + i ImF)
                    S_re = self._row(y0=1.0, a=-im_row[i])
                    S_im = self._row(a=re_row[i])
                    cF_re = self._row(y0=kin[i] / g[i], ImF=(ell, (kin[i] / g[i]) * K[i]))
                    cF_im = self._row(ImF=(ell, (kin[i] / g[i]) * np.eye(M)[i]))
                    if "gram" in spec.uv_parts:
                        rh = self._row(rho=(ell, np.eye(M)[i]))
                        yield [[y0 + S_re, S_im,      rt2 * cF_re],
                               [S_im,      y0 - S_re, rt2 * cF_im],
                               [rt2 * cF_re, rt2 * cF_im, rh]]
                    if "ff" in spec.uv_parts and i in set(idx_hi):
                        n = list(idx_hi).index(i)
                        b = self._row(y0=ffb[ell] / ffk[ell][n])
                        yield [[b, cF_re, cF_im], [cF_re, b, zero], [cF_im, zero, b]]
                # rho_hat >= 0
                for i in range(M):
                    yield [[self._row(rho=(ell, np.eye(M)[i]))]]
                # FESR (3.73)
                if "fesr" in spec.uv_parts:
                    wave = "S0" if ell == 0 else "P1"
                    residuals = []
                    for n in C.MOMENTS[ell]:
                        if (wave, n) in free_moments(spec):
                            continue                      # moment-range diagnostic: box not imposed
                        mr = moments[ell,n] if self.precise else C.moment_row(M,n)*g**2
                        mom = self._row(rho=(ell, mr))
                        t, d = tgt[(wave, n)], tol[(wave, n)]
                        if spec.sr_caliber == "SR-d":
                            residuals.append(mom - self._row(y0=t))
                        else:
                            yield [[self._row(y0=t + d) - mom]]
                            yield [[mom - self._row(y0=t - d)]]
                    if spec.sr_caliber == "SR-d" and residuals:
                        # per-wave L2 ball of radius eps_SR over the two moments
                        yield self._arrow(residuals, spec.eps_sr)

        # ---- optional section constraint f00(3) = x
        if self.fix_f00 is not None:
            r = self._row(a=self.f00)
            yield [[self._row(y0=self.fix_f00) - r]]
            yield [[r - self._row(y0=self.fix_f00)]]

        # ---- optional near-optimal slab of a completed support (face diagnostic)
        if self.face is not None:
            yield [[self.face_row()]]

    def _regulariser_blocks(self):
        """``|rho_{a,ij}| <= Mreg`` as two 1x1 blocks per double-density value.

        This is the M-regularisation of He-Kruczenski 2103.11484 section 3:
        the nodal Mandelstam parametrisation contains directions that barely
        change the partial waves at the collocation nodes, so unitarity there
        does not bound them; a cap on the primal variables removes them.  The
        rows are exact: with a coordinate basis ``V`` they are the binary64
        rows of ``V`` on the double-density slots, so the blocks are identical
        in the float and Arb assembly paths.

        Each block is divided by ``Mreg`` -- ``1 -/+ rho_i / Mreg >= 0`` -- which
        is the same inequality (a positive scaling of a 1x1 block) but keeps the
        constant entry at 1 like every other block.  With the unscaled form the
        initial primal objective carries ``7550 * initialScale * Mreg`` and SDPB
        diverged (maxComplementarity) at Mreg = 1e6 in the first gate run.
        """
        lay, bound, norm = self.ops.lay, float(self.spec.reg_bound), self.spec.reg_norm
        idx = list(range(lay.r1.start, lay.r1.stop)) + list(range(lay.r2.start, lay.r2.stop))
        one = np.zeros(self.n_vars); one[0] = 1.0
        for k, i in enumerate(idx):
            row = np.zeros(self.n_vars)                # r_i = rho_i / bound
            if self.basis is not None:
                row[self.i_a:self.i_a + self.n_a] = self.basis[i] / bound
            else:
                row[self.i_a + i] = 1.0 / bound
            if norm == "linf":
                upper, lower = -row, row.copy()        # 1 - r_i >= 0, 1 + r_i >= 0
                upper[0] = lower[0] = 1.0
                yield [[upper]]
                yield [[lower]]
                continue
            t = np.zeros(self.n_vars); t[self.i_aux + k] = 1.0
            yield [[t, row], [row, one]]              # t_i >= r_i^2
            if norm == "l4":
                u = np.zeros(self.n_vars); u[self.i_aux + self.n_rho + k] = 1.0
                yield [[u, t], [t, one]]              # u_i >= t_i^2  (=> u_i >= r_i^4)
        if norm in ("l2", "l4"):
            total = one.copy()                        # 1 - sum_i (t_i or u_i) >= 0
            off = self.i_aux + (self.n_rho if norm == "l4" else 0)
            total[off:off + self.n_rho] = -1.0
            yield [[total]]

    def _arrow(self, rows, e: float):
        """``||v||_2 <= e``  ->  [[e y0, v^T], [v, e y0 I]] >= 0."""
        n = len(rows)
        d, zero = self._row(y0=e), np.zeros(self.n_vars)
        blk = [[d] + list(rows)]
        for i in range(n):
            blk.append([rows[i]] + [d if j == i else zero for j in range(n)])
        return blk

    # ------------------------------------------------- face diagnostic rows
    def node_rows(self, wave: str, k: int):
        """``(Re kappa*h, Im kappa*h)`` coefficient rows of ``wave`` at node ``k``.

        The same rows enter the unitarity disk and the Gram block at that node
        (``S = 1 + i kappa h``), so ``Im kappa*h = 1 - Re S`` and
        ``Re kappa*h = Im S`` are exactly linear in the model coordinates.
        """
        I, ell = {"S0": (0, 0), "P1": (1, 1), "S2": (2, 0)}[wave]
        if self.precise:
            ctx.prec = self.precise.bits
            rows = self.precise.rows(I, ell, node=k) * self.precise.kap[k]
            proj = self.precise.projected(rows, self.precise_basis)
            return proj[0], proj[1]
        if wave in ("S0", "P1"):
            re_rows, im_rows = self.gram[ell]
            return re_rows[k], im_rows[k]
        a0 = self.ops.index.index((I, ell)) * self.spec.M + k
        re_row, im_row = self.ops.h_re[a0], self.ops.h_im[a0]
        if self.basis is not None:
            re_row, im_row = re_row @ self.basis, im_row @ self.basis
        return re_row, im_row

    def watson_row(self) -> np.ndarray:
        """``sum_k Re(t_k) Re h_k + (Im t_k - 1) Im h_k`` over the registered waves and nodes (unsigned)."""
        f = self.functional
        from .watson import weight_of
        total = np.zeros(self.n_vars, dtype=object if self.precise else float)
        for wave in f["waves"]:
            for i, ((tr, ti), k) in enumerate(zip(f["targets"][wave], f["nodes"])):
                re_row, im_row = self.node_rows(wave, k)
                wk = weight_of(f, wave, i)
                if self.precise:
                    total = total + self._row(a=arb(str(wk)) * (arb(str(tr)) * re_row + arb(str(ti - 1.0)) * im_row))
                else:
                    total = total + self._row(a=wk * (float(tr) * re_row + float(ti - 1.0) * im_row))
        return total

    def functional_row(self) -> np.ndarray:
        """Unsigned coefficient row of the registered secondary functional."""
        f = self.functional
        if f["kind"] == "watson":
            return self.watson_row()
        if f["kind"] == "SRmom":
            ell, n, M = f["ell"], f["node"], self.spec.M
            if self.precise:
                mr = self.precise.uv_data(self.spec)[6][ell, n]
            else:
                mr = C.moment_row(M, n) * FFM.gram_scale(ell, self.ops.s) ** 2
            return self._row(rho=(ell, mr))
        if f["kind"] in ("ImF", "rho"):
            unit = np.zeros(self.spec.M); unit[f["node"]] = 1.0
            return self._row(**{"rho" if f["kind"] == "rho" else "ImF": (f["ell"], unit)})
        re_row, im_row = self.node_rows(f["wave"], f["node"])
        return self._row(a=im_row if f["kind"] == "ImKH" else re_row)

    def face_row(self) -> np.ndarray:
        """``d.(f00, f11) - (value - margin) >= 0`` as one 1x1 block row."""
        d = ([arb(str(v)) for v in self.face["direction"]] if self.precise
             else np.asarray(self.face["direction"], dtype=float))
        floor = self.face["value"] - self.face["margin"]
        return self._row(a=d[0] * self.f00 + d[1] * self.f11) - self._row(y0=floor)

    # ----------------------------------------------------------------- write
    def objective(self) -> np.ndarray:
        if self.functional is not None:
            row = self.functional_row()
            return row if self.functional["sense"] == "max" else -row
        d = ([arb(str(v)) for v in self.direction] if self.precise else np.asarray(self.direction, dtype=float))
        return self._row(a=d[0] * self.f00 + d[1] * self.f11)

    def write(self, path: str) -> dict:
        """Stream the PMP to ``path``; returns a small summary."""
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        if hasattr(self,'prepared'):
            from .precision_pmp import write_prepared
            return write_prepared(self,path)
        num, nb, sizes = self._num_vec, 0, {}
        with open(path, "w") as fh:
            fh.write('{"objective":')
            fh.write(json.dumps(num(self.objective())))
            norm = np.zeros(self.n_vars); norm[0] = 1.0
            fh.write(',"normalization":' + json.dumps(num(norm)))
            fh.write(',"PositiveMatrixWithPrefactorArray":[')
            for blk in self._blocks():
                if nb:
                    fh.write(",")
                fh.write(_PRE + ',"polynomials":[')
                for r, row in enumerate(blk):
                    fh.write("," if r else "")
                    # polynomials[row][col][variable] is itself a list of
                    # coefficients; degree 0 means exactly one per variable.
                    fh.write("[" + ",".join(
                        "[" + ",".join('["%s"]' % s for s in num(cell)) + "]"
                        for cell in row) + "]")
                fh.write("]}")
                nb += 1
                sizes[len(blk)] = sizes.get(len(blk), 0) + 1
            fh.write("]}")
        self.n_blocks = nb
        return {"path": path, "bytes": os.path.getsize(path), "n_vars": self.n_vars,
                "n_blocks": nb, "block_sizes": sizes,
                "n_unitarity_disks": int(self.keep.sum()),
                "reduce_basis": self.basis is not None,
                "n_a": self.n_a, "n_aux": self.n_aux, "digits": self.digits,
                "operator_precision": getattr(self,"precision_info", {"source_dps":17})}

    def _num_vec(self, v: np.ndarray) -> list[str]:
        return [_num(x, self.digits) for x in v]

    # ------------------------------------------------------- read a solution
    def read_solution(self, y_path: str) -> dict:
        """Map SDPB's ``out/y.txt`` back onto the model variables.

        SDPB eliminates the normalised variable, so ``y.txt`` holds
        ``n_vars - 1`` numbers: our layout minus ``y0``, in order.  Verified on
        the known-answer cone test (variables ``y0, x, z, w``; ``y.txt`` carried
        exactly ``x, z, w``).  The strings are arbitrary precision; we keep the
        full text as well as the float64 cast.
        """
        with open(y_path) as fh:
            head = fh.readline().split()
            vals = [ln.strip() for ln in fh if ln.strip()]
        n = int(head[0])
        if n != self.n_vars - 1 or len(vals) != n:
            raise ValueError(f"y.txt has {head[0]} rows / {len(vals)} values, "
                             f"expected {self.n_vars - 1}")
        y = np.array([float(v) for v in vals])
        a = y[self.i_a - 1:self.i_a - 1 + self.n_a]
        c = (self.basis @ a) if self.basis is not None else a
        out = {"y": y, "y_text": vals, "a": a, "c": c}
        if self.spec.uv:
            M = self.spec.M
            out["ImF"] = y[self.i_ImF - 1:self.i_ImF - 1 + 2 * M].reshape(2, M)
            out["rho_hat"] = y[self.i_rho - 1:self.i_rho - 1 + 2 * M].reshape(2, M)
        return out

    def verify(self, sol: dict) -> dict:
        """Replay all imposed constraint families in the original coordinates.

        This supplies numerical feasibility only. Solver convergence and an
        independent dual certificate are separate and never inferred here.
        """
        if self.precise:
            from .precision_pmp import verify
            return verify(self,sol)
        from .verify import solution_report
        return solution_report(self, sol)
