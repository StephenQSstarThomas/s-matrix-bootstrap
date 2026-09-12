"""Conic model assembly and solution -- paper section 3, task section 5.

Decision variables
    c          amplitude parameters (3.62), length 1 + 2M + M^2 + M(M+1)/2
    ImF        Im F_ell(s_i), ell = 0,1                      (3.65)
    rho_hat    rescaled current spectral densities, >= 0     (3.69)

Constraints
    unitarity  1500 rotated second-order cones, |h|^2 <= 2 Im h with h = kappa f
    Gram       100 realified 3x3 positive semidefinite blocks (3.68)
    FESR       4 linear boxes (3.73)-(3.74)
    FF         14 second-order cones (3.75)
    chiral     (3.64) in one of the pre-registered norms
    B          optional ell-4 bound on the double spectral densities

Two reversible reparametrisations are applied purely for conditioning, both
exact:

  * unitarity.  For Lambda > 0, |h|^2 <= 2 Im h  <=>  |h/Lambda|^2 <=
    2 (Im h)/Lambda^2 -- multiply through by Lambda^2.  We use
    Lambda_ell(s) = ((sqrt(s)-2)/(sqrt(s)+2))^(ell/2), the centrifugal scale, so
    the rows of the cone no longer span 40 orders of magnitude.
  * Gram.  Congruence by V = diag(1,1,1/g) maps (3.68) to the same matrix with
    cF -> cF/g and rho -> rho/g^2.  Congruence by an invertible matrix preserves
    positive semidefiniteness, so the feasible set is untouched.
"""
from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np

from . import constraints as C
from . import formfactor as FFM
from .grid import M_Q, S0, kappa, lambda_rescale, s_nodes
from .projector import Layout, PartialWaveOperator, ells_for

PROJECTION_S = 3.0            # the plane is (f00(3), f11(3)), section 4.1


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
    m_q: float = M_Q
    B: float | None = None              # bound on (rho1, rho2); None = absent
    B_norm: str = "l2"                  # l2 (stronger, cheap) | l4 (pre-registered)
    cone_scaling: str = "none"          # none | centrifugal | rownorm (all exact)
    sparsify: float = 0.0               # zero entries below this fraction of their row scale
    reduce_basis: bool = False          # exact projection onto span(all rows)
    basis_tol: float = 1e-12
    tag: str = ""

    def key(self) -> str:
        return hashlib.sha256(repr(sorted(self.__dict__.items())).encode()).hexdigest()[:16]


class Operators:
    """Cached partial-wave rows for one (M, L)."""

    _cache: dict = {}

    def __new__(cls, M: int, L: int):
        if (M, L) in cls._cache:
            return cls._cache[(M, L)]
        self = super().__new__(cls)
        self._init(M, L)
        cls._cache[(M, L)] = self
        return self

    def _init(self, M: int, L: int) -> None:
        self.M, self.L = M, L
        self.op = PartialWaveOperator(M)
        self.lay = Layout(M)
        self.s = s_nodes(M)
        self.kap = kappa(self.s)
        self.index = [(I, e) for I in (0, 1, 2) for e in ells_for(I, L)]
        n = self.lay.n
        nw = len(self.index)
        re = np.empty((nw, M, n))
        im = np.empty((nw, M, n))
        lam = np.empty((nw, M))
        for a, (I, e) in enumerate(self.index):
            lam[a] = lambda_rescale(self.s, e)
            for k in range(M):
                r = self.op.rows(I, e, self.s[k], k)
                re[a, k], im[a, k] = r[0], r[1]
        # Fold kappa in.  The cone may then be rescaled by any positive
        # per-row Lambda -- |h|^2 <= 2 Im h  <=>  |h/L|^2 <= 2 (Im h)/L^2, just
        # multiply through by L^2 > 0 -- which is exact for every choice.
        #
        # We nevertheless keep L = 1 by default.  Rescaling does flatten the row
        # norms (raw |h| rows span 5e-85 .. 34 because of centrifugal
        # suppression at the threshold node), but it inflates the *feasible set*
        # by the same factor: the unscaled cone already implies |h| <= 2, i.e. a
        # unit disk centred at i, whereas after division by L ~ 1e-40 the scaled
        # variables range over 1e40 and Clarabel fails with NumericalError at
        # step 0.  Tiny rows are harmless -- they are slack constraints -- while
        # a feasible set spanning 40 decades is not.  ``cone_scaling`` keeps the
        # alternatives available and the row-norm distributions are reported.
        hre = re.reshape(-1, n) * np.repeat(self.kap[None, :], nw, 0).reshape(-1, 1)
        him = im.reshape(-1, n) * np.repeat(self.kap[None, :], nw, 0).reshape(-1, 1)
        nu = np.maximum(np.abs(hre).max(axis=1), np.abs(him).max(axis=1))
        self.lambda_analytic = lam.reshape(-1) ** 2
        self.nu_measured = nu
        ok = self.lambda_analytic > 1e-7
        self.lambda_agreement = nu[ok] / self.lambda_analytic[ok]
        self.h_re, self.h_im = hre, him
        self.row_norm_raw = np.maximum(np.abs(hre), np.abs(him)).max(axis=1)
        self.set_cone_scaling("none", 0.0)
        # unscaled S0 / P1 rows, needed by the Gram blocks
        self.gram_rows = {}
        for ell, I in ((0, 0), (1, 1)):
            a = self.index.index((I, ell))
            self.gram_rows[ell] = (re[a] * self.kap[:, None], im[a] * self.kap[:, None])
        # subthreshold rows
        self.f_proj = {"f00": self.op.rows(0, 0, PROJECTION_S)[0],
                       "f11": self.op.rows(1, 1, PROJECTION_S)[0]}
        self.chi_rows = []
        for sj in C.CHIRAL_POINTS:
            r01, r21 = C.chiral_ratios(sj)
            f00 = self.op.rows(0, 0, sj)[0]
            f11 = self.op.rows(1, 1, sj)[0]
            f20 = self.op.rows(2, 0, sj)[0]
            self.chi_rows.append(f00 - r01 * f11)
            self.chi_rows.append(f20 - r21 * f11)
        self.chi_rows = np.array(self.chi_rows)
        self.basis = None

    def build_basis(self, tol: float = 1e-12) -> np.ndarray:
        """Orthonormal basis V of the span of every row the model can see.

        Each constraint and the objective is a linear functional of c, and they
        span a subspace of dimension k << n.  Writing c = V a is therefore an
        *exact* reparametrisation: directions orthogonal to the span move no
        constraint and no objective, so they can neither help feasibility nor
        improve the optimum.  It removes the 2500-dimensional degeneracy that
        otherwise leaves the interior-point method with a rank-deficient KKT
        system.  (It cannot be used together with the B regulariser, which is a
        norm on c itself -- see ModelSpec.B.)
        """
        rows = [self.P_re, self.P_im,
                self.f_proj["f00"][None, :], self.f_proj["f11"][None, :],
                self.chi_rows]
        R = np.vstack(rows)
        R = R / np.maximum(np.abs(R).max(axis=1, keepdims=True), 1e-300)
        _, sv, Vt = np.linalg.svd(R, full_matrices=False)
        k = int((sv > tol * sv[0]).sum())
        self.basis_singular_values = sv
        self.basis_tol = tol
        self.basis = Vt[:k].T
        return self.basis

    def set_cone_scaling(self, mode: str, drop_tiny: float = 0.0) -> None:
        """Choose the (exact) per-row rescaling of the unitarity cone."""
        if mode == "none":
            lam2 = np.ones_like(self.nu_measured)
        elif mode == "centrifugal":                  # Lambda_ell(s)^2 of (task 5)
            lam2 = np.maximum(self.lambda_analytic, 1e-16)
        elif mode == "rownorm":                      # measured row magnitude
            lam2 = np.where(self.nu_measured > 0.0, self.nu_measured, 1.0)
        else:
            raise ValueError(mode)
        lam1 = np.sqrt(lam2)
        self.cone_scaling = mode
        self.P_re = self.h_re / lam1[:, None]
        self.P_im = self.h_im / lam1[:, None]
        self.R_im = self.h_im / lam2[:, None]
        if drop_tiny > 0.0:
            # Entries this far below their row's scale cannot move a constraint
            # by more than ``drop_tiny`` relative, i.e. far below tol_feas; they
            # only destroy sparsity and conditioning of the KKT system.  The
            # exact operator is kept in ``h_re``/``h_im`` so every solution is
            # re-verified against the unsparsified constraints.
            for Mx in (self.P_re, self.P_im, self.R_im):
                sc = np.abs(Mx).max(axis=1, keepdims=True)
                Mx[np.abs(Mx) < drop_tiny * np.maximum(sc, 1e-300)] = 0.0
        self.sparsify = drop_tiny
        self.row_norm_scaled = np.maximum(np.abs(self.P_re), np.abs(self.R_im)).max(axis=1)
        self.row_scale = lam2


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
        keep = np.ones(len(self.ops.index) * M, dtype=bool)
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
            idx_hi, ffb = C.ff_asymptotic_bounds(M, spec.m_q, spec.eps_ff)
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
                    for i in idx_hi:
                        cons.append(cp.SOC(cp.Constant(ffb[ell] / g[i]),
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
