"""Element-wise comparison of two kernel exports in 60-digit decimal arithmetic (gate G1).

For each matrix: the largest relative difference over entries where both are non-zero (``max_rel``), the same
restricted to entries above 1e-30 of the matrix's largest entry (``max_rel_significant``), the number of entries
that are exactly zero in one export only (``zero_mismatch``) and the largest magnitude among those partners
(``max_abs_small``), plus the flat index of the worst relative difference.
"""
from __future__ import annotations

from decimal import Decimal, localcontext

from .kernels import load


def _compare(ea, eb, rel_floor=Decimal("1e-30")):
    with localcontext() as ctx:
        ctx.prec = 60
        da, db = [Decimal(x) for x in ea], [Decimal(x) for x in eb]
        top = max((abs(x) for x in db), default=Decimal(0))
        worst = sig = Decimal(0)
        worst_at, zero_mismatch, small = -1, 0, Decimal(0)
        for k, (a, b) in enumerate(zip(da, db)):
            if a == 0 and b == 0:
                continue
            if a == 0 or b == 0:
                zero_mismatch += 1
                small = max(small, abs(a) + abs(b))
                continue
            r = abs(a - b) / max(abs(a), abs(b))
            if r > worst:
                worst, worst_at = r, k
            if max(abs(a), abs(b)) > rel_floor * top and r > sig:
                sig = r
    return {"n": len(da), "max_rel": float(worst), "max_rel_significant": float(sig), "worst_at": worst_at,
            "zero_mismatch": zero_mismatch, "max_abs_small": float(small), "max_abs": float(top)}


def compare_exports(dir_a, dir_b, names=None) -> dict:
    A, B = load(dir_a), load(dir_b)
    names = names or sorted(set(A.names) & set(B.names))
    out = {}
    for n in names:
        sa, sb = A.manifest["matrices"][n]["shape"], B.manifest["matrices"][n]["shape"]
        out[n] = {"error": f"shape {sa} vs {sb}"} if sa != sb else {"shape": sa, **_compare(A.exact(n), B.exact(n))}
    return out
