"""Conic model assembly and solution -- paper section 3, task section 5.

Decision variables
    c          amplitude parameters (3.62), length 1 + 2M + M^2 + M(M+1)/2
    ImF        Im F_ell(s_i), ell = 0,1                      (3.65)
    rho_hat    rescaled current spectral densities, >= 0     (3.69)

Constraints
    unitarity  3 L M rotated second-order cones, |h|^2 <= 2 Im h, h = kappa f
    Gram       2 M realified 3x3 positive semidefinite blocks (3.68)
    FESR       4 linear boxes (3.73)-(3.74)
    FF         2 x (nodes above s0) second-order cones (3.75)
    chiral     (3.64) in one of the pre-registered norms
    B          optional bound on the double spectral densities

Conditioning.  Three reparametrisations are available, all exact:

  * unitarity cone.  For any Lambda > 0, |h|^2 <= 2 Im h is equivalent to
    |h/Lambda|^2 <= 2 (Im h)/Lambda^2 -- multiply through by Lambda^2.  The raw
    rows span 5e-85 .. 34 because of centrifugal suppression at the threshold
    node.  ``cone_scaling`` selects Lambda: "centrifugal" is the analytic
    ((sqrt(s)-2)/(sqrt(s)+2))^(ell/2), "rownorm" the measured row magnitude (the
    one that works; see assembly.set_cone_scaling for why the analytic form
    underflows), "none" leaves the cone alone.
  * Gram congruence.  V_i = diag(1, 1, 1/g_i) with g_i = k_ell(s_i), the (2.33)
    kinematic factor, turns (3.68) into [[1,S,F],[S*,1,F*],[F*,F,rho/k^2]] --
    the (1,3) entry becomes the form factor itself.  Congruence by an invertible
    matrix preserves positive semidefiniteness, so the feasible set is untouched.
  * basis reduction.  Every constraint and the objective is a linear functional
    of c, and they span a subspace of dimension about 1780 out of 3876; c = V a
    is then exact.  Off by default -- it does not converge at M >= 30.

The unitarity disks are normally imposed through the constraint generation of
:func:`smatrix_bootstrap.sdp.runner.solve_generated` rather than all at once:
imposing a subset is a relaxation, and a returned point that satisfies all of
them is optimal for the full problem, which makes the answer certified.
``ModelSpec.disk_mask`` is how the generator passes the current subset in.
"""
from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass

import cvxpy as cp
import numpy as np

from . import constraints as C
from . import formfactor as FFM
from .assembly import PROJECTION_S, Operators
from .grid import M_Q


@dataclass
class ModelSpec:
    M: int = 50
    L: int = 10
    chiral: bool = False
    chi_caliber: str = "chi-b"          # chi-a | chi-b | chi-c
    eps_chi: float = C.EPS_CHI_MAIN
    uv: bool = False                    # Gram + FESR + form-factor asymptotics
    uv_parts: tuple = ("gram", "fesr", "ff")   # for diagnosing infeasibility
    sr_caliber: str = "SR-b"            # SR-a | SR-b | SR-c
    eps_ff: float = C.EPS_FF
    ff_frozen_at_s0: bool = True        # (3.75) factor taken at s0, per the paper text
    m_q: float = M_Q
    B: float | None = None              # bound on (rho1, rho2); None = absent
    B_norm: str = "l2"                  # l2 (stronger, cheap) | l4 (pre-registered)
    cone_scaling: str = "none"          # none | centrifugal | rownorm (all exact)
    sparsify: float = 0.0               # zero entries below this fraction of their row scale
    reduce_basis: bool = False          # exact projection onto span(all rows)
    basis_tol: float = 1e-12
    tag: str = ""
    disk_mask: object = None             # bool array over (wave, node): which disks to impose

    def key(self) -> str:
        return hashlib.sha256(repr(sorted(self.__dict__.items())).encode()).hexdigest()[:16]


class Model:
    """Builds the cvxpy problem once; the objective direction is a Parameter."""

    def __init__(self, spec: ModelSpec) -> None:
        self.spec = spec
        self.ops = Operators(spec.M, spec.L)
        self.ops.set_cone_scaling(spec.cone_scaling, spec.sparsify)
        M = spec.M
        if spec.reduce_basis:
            V = self.ops.build_basis(spec.basis_tol)
            self.a = cp.Variable(V.shape[1], name="a")
            self.c = V @ self.a
            P_re, P_im, R_im = (self.ops.P_re @ V, self.ops.P_im @ V, self.ops.R_im @ V)
            f00r, f11r = self.ops.f_proj["f00"] @ V, self.ops.f_proj["f11"] @ V
            chir = self.ops.chi_rows @ V
            gram = {e: (a @ V, b @ V) for e, (a, b) in self.ops.gram_rows.items()}
            self.basis = V
        else:
            self.a = cp.Variable(self.ops.lay.n, name="c")
            self.c = self.a
            P_re, P_im, R_im = self.ops.P_re, self.ops.P_im, self.ops.R_im
            f00r, f11r = self.ops.f_proj["f00"], self.ops.f_proj["f11"]
            chir = self.ops.chi_rows
            gram = self.ops.gram_rows
            self.basis = None
        cons = []

        # ---- unitarity: ||(p_re, p_im)||^2 <= 2 r,  i.e. ||(2p,2r-1)|| <= 2r+1
        #
        # The 2x2 block of the Gram matrix (3.68) is [[1,S],[S*,1]] >= 0, which
        # is |S| <= 1 -- the same constraint.  Imposing both for S0 and P1 leaves
        # the interior-point method with an exactly redundant pair at every node
        # of those two waves, so when the Gram blocks are on we drop the S0/P1
        # rows from the cone list rather than state the constraint twice.
        keep = (np.ones(len(self.ops.index) * M, dtype=bool)
                if spec.disk_mask is None else np.asarray(spec.disk_mask, dtype=bool).copy())
        if spec.uv and "gram" in spec.uv_parts:
            for ell, I in ((0, 0), (1, 1)):
                a0 = self.ops.index.index((I, ell))
                keep[a0 * M:(a0 + 1) * M] = False
        self.n_unitarity_cones = int(keep.sum())
        p_re = P_re[keep] @ self.a
        p_im = P_im[keep] @ self.a
        r = R_im[keep] @ self.a
        cons.append(cp.SOC(2 * r + 1, cp.vstack([2 * p_re, 2 * p_im, 2 * r - 1]), axis=0))

        # ---- density regularisation (only if requested)
        if spec.B is not None:
            # rho is a linear function of the decision variable in either
            # parametrisation: c = V a in the reduced one, c = a in the full one.
            if self.basis is not None:
                Vr = np.vstack([self.basis[self.ops.lay.r1], self.basis[self.ops.lay.r2]])
                rho_block = Vr @ self.a
            else:
                rho_block = cp.hstack([self.a[self.ops.lay.r1], self.a[self.ops.lay.r2]])
            if spec.B_norm == "l4":
                cons.append(cp.pnorm(rho_block, 4) <= spec.B)
            elif spec.B_norm == "l2":
                # ||rho||_4 <= ||rho||_2, so an ell-2 ball of the same radius is
                # strictly STRONGER than the pre-registered ell-4 bound: if it
                # comes out inactive, so is the ell-4 one, and the optimum is the
                # unregularised one.  It costs one cone instead of the ~7500 that
                # cvxpy's geometric-mean canonicalisation of pnorm(.,4) creates.
                cons.append(cp.norm(rho_block, 2) <= spec.B)
            else:
                raise ValueError(spec.B_norm)
            self.rho_block = rho_block

        # ---- chiral symmetry breaking (3.64)
        if spec.chiral:
            res = chir @ self.a
            e = spec.eps_chi
            if spec.chi_caliber == "chi-a":
                cons += [cp.abs(res) <= e]
            elif spec.chi_caliber == "chi-b":
                cons += [cp.norm(res, 2) <= e]
            elif spec.chi_caliber == "chi-c":
                cons += [cp.norm(res[0::2], 2) <= e, cp.norm(res[1::2], 2) <= e]
            else:
                raise ValueError(spec.chi_caliber)

        # ---- form factor / SVZ sector
        self.ImF = self.rho_hat = None
        if spec.uv:
            K = FFM.hilbert_kernel(M)
            self.ImF = cp.Variable((2, M), name="ImF")
            self.rho_hat = cp.Variable((2, M), name="rho_hat", nonneg=True)
            idx_hi, ffb, ffk = C.ff_asymptotic_bounds(M, spec.m_q, spec.eps_ff,
                                                     spec.ff_frozen_at_s0)
            tgt = C.printed_targets()
            tol = C.sr_tolerances(spec.sr_caliber)
            for ell in (0, 1):
                g = FFM.gram_scale(ell, self.ops.s)   # per-node congruence
                kin = FFM.kinematic_factor(ell, self.ops.s)
                ReF = 1.0 + K @ self.ImF[ell]
                cF_re = cp.multiply(kin / g, ReF)
                cF_im = cp.multiply(kin / g, self.ImF[ell])
                S_re = 1.0 - gram[ell][1] @ self.a
                S_im = gram[ell][0] @ self.a
                rh = self.rho_hat[ell]
                if "gram" in spec.uv_parts:
                    for i in range(M):
                        cons.append(_gram_psd(S_re[i], S_im[i], cF_re[i], cF_im[i], rh[i]))
                # FESR (3.73): moments of the unrescaled rho = g^2 rho_hat
                wave = "S0" if ell == 0 else "P1"
                if "fesr" in spec.uv_parts:
                    for n in C.MOMENTS[ell]:
                        mom = (C.moment_row(M, n) * g ** 2) @ rh
                        cons += [mom - tgt[(wave, n)] <= tol[(wave, n)],
                                 tgt[(wave, n)] - mom <= tol[(wave, n)]]
                # (3.75) on the nodes above s0
                if "ff" in spec.uv_parts:
                    # |k_used * F(s_i)| <= bound; cF_re/cF_im carry F itself
                    # because the Gram congruence divided by k(s_i)
                    for n, i in enumerate(idx_hi):
                        cons.append(cp.SOC(cp.Constant(ffb[ell] / ffk[ell][n]),
                                           cp.hstack([cF_re[i], cF_im[i]])))
        self.constraints = cons
        self.direction = cp.Parameter(2, name="d")
        self.f00 = f00r @ self.a
        self.f11 = f11r @ self.a
        lam_row = self.ops.op.lambda_row()
        self.lam = (lam_row @ V if self.basis is not None else lam_row) @ self.a
        self.extra = []
        self.problem = None

    # ------------------------------------------------------------------
    def finalize(self, extra=(), objective="plane"):
        self.extra = list(extra)
        if objective == "lambda":
            obj = cp.Maximize(self.lam)
        else:
            obj = cp.Maximize(self.direction[0] * self.f00
                              + self.direction[1] * self.f11)
        self.problem = cp.Problem(obj, self.constraints + self.extra)
        return self.problem

    def solve(self, d, solver="CLARABEL", **kw):
        self.direction.value = np.asarray(d, dtype=float)
        t0 = time.time()
        opts = dict(max_iter=500, tol_gap_abs=1e-8, tol_gap_rel=1e-8,
                    tol_feas=1e-8, time_limit=7200.0)
        opts.update(kw)
        if solver == "SCS":
            opts = {"eps": kw.get("eps", 1e-7), "max_iters": kw.get("max_iters", 100000)}
        try:
            self.problem.solve(solver=solver, **opts)
        except cp.error.SolverError as exc:
            return {"status": "SolverError", "message": str(exc),
                    "seconds": time.time() - t0}
        out = {"status": self.problem.status, "seconds": time.time() - t0,
               "objective": _f(self.problem.value),
               "f00_3": _f(self.f00.value), "f11_3": _f(self.f11.value),
               "lambda": _f(self.lam.value),
               "cone_scaling": self.spec.cone_scaling,
               "solver": solver, "direction": list(map(float, d))}
        st = self.problem.solver_stats
        if st is not None:
            out["iterations"] = st.num_iters
            out["solve_time"] = st.solve_time
        if self.a.value is not None:
            out["n_reduced"] = int(self.a.size)
            out["solution_sha256"] = hashlib.sha256(
                np.ascontiguousarray(self.a.value).tobytes()).hexdigest()
        return out

    def solution(self):
        av = self.a.value
        c = None if av is None else (self.basis @ av if self.basis is not None else np.asarray(av))
        d = {"c": c}
        if self.ImF is not None and self.ImF.value is not None:
            d["ImF"] = np.asarray(self.ImF.value)
            d["rho_hat"] = np.asarray(self.rho_hat.value)
        return d


def _f(x):
    return None if x is None else float(x)


def _gram_psd(S_re, S_im, F_re, F_im, rho):
    """Realified (3.68):  [[1,S,cF],[S*,1,cF*],[cF*,cF,rho]] >= 0  as a 6x6."""
    one = cp.Constant(1.0)
    zero = cp.Constant(0.0)
    Re = cp.bmat([[one, S_re, F_re], [S_re, one, F_re], [F_re, F_re, rho]])
    Im = cp.bmat([[zero, S_im, F_im], [-S_im, zero, -F_im], [-F_im, F_im, zero]])
    return cp.bmat([[Re, -Im], [Im, Re]]) >> 0
