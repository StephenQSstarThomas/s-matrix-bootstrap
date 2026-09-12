"""Assembly and caching of the partial-wave operator rows for one (M, L).

Split out of :mod:`smatrix_bootstrap.sdp.problem` so that operator assembly and
conic-model construction stay separate responsibilities (and each file stays
inside the repository's size convention).
"""
from __future__ import annotations

import numpy as np

from . import constraints as C
from .grid import kappa, lambda_rescale, s_nodes
from .projector import Layout, PartialWaveOperator, ells_for

PROJECTION_S = 3.0            # the plane is (f00(3), f11(3)), section 4.1


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


