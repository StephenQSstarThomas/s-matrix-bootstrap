"""The rho pole: S_II = 1/S_I (v4 (5.1)), so the second-sheet pole is the zero of S^1_1 on the first sheet in the
lower half plane (PREREGISTRATION amendment 1).  S_I(s) = 1 + i pi sqrt((s-4)/s) f^1_1(s) with f from evaluate.py.
"""
from __future__ import annotations

import numpy as np

M_PI_MEV = 139.57


def find_zero(fun, s0, tol=1e-13, maxit=60):
    """Secant iteration on a complex function from s0 (and s0 (1 + 1e-3))."""
    a, b = complex(s0), complex(s0) * (1 + 1e-3)
    fa, fb = fun(a), fun(b)
    for _ in range(maxit):
        if fb == fa:
            break
        c = b - fb * (b - a) / (fb - fa)
        a, fa, b, fb = b, fb, c, fun(c)
        if abs(b - a) <= tol * abs(b):
            return b
    raise RuntimeError(f"no convergence from {s0}: last {b}, |f| = {abs(fb)}")


def sqrt_s_mev(s):
    """sqrt(s) in MeV with Im <= 0 for a lower-half-plane s."""
    r = np.sqrt(complex(s))
    return (r if r.real >= 0 else -r) * M_PI_MEV


def s_matrix_p1(grid, amp, prescription="interp"):
    """S^1_1(s) on the first sheet as a function of complex s (Im s < 0), from a lifted amplitude vector."""
    from flint import acb
    from .evaluate import f_row, value

    def S(s):
        sa = acb(s.real, s.imag)
        f = value(f_row(grid, 1, 1, ("point", sa), prescription), amp)
        sig = ((sa - 4) / sa).sqrt()
        return complex((1 + acb(0, 1) * acb.pi() * sig * f).mid())
    return S
