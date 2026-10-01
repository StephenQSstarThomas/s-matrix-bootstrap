"""Second-sheet pole = first-sheet zero of S (v4 (5.1)); Phase B Task B4."""
from __future__ import annotations

import numpy as np


def test_zero_of_a_unimodular_s_matrix_is_recovered():
    from smatrix_bootstrap.gtb.pole import find_zero, sqrt_s_mev
    sp = (0.829 / 0.13957 - 0.071j / 0.13957) ** 2          # sqrt(s_p) = 829 - 71 i MeV, m_pi units
    S = lambda s: (s - sp) / (s - np.conj(sp))
    assert abs(abs(S(30.0 + 0j)) - 1) < 1e-14              # unimodular on the real axis
    z = find_zero(S, (0.80 / 0.13957 - 0.04j / 0.13957) ** 2)
    assert abs(z - sp) < 1e-10 * abs(sp)
    w = sqrt_s_mev(z)
    assert abs(w.real - 829) < 1e-6 and abs(w.imag + 71) < 1e-6
