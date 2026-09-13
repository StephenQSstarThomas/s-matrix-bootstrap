"""Independent SDPB primal-weight replay against the emitted degree-0 PMP.

For max d.y with A0+Ai.y >= 0, Z>=0 and sum tr(Z Ai)=-d yield an upper
bound sum tr(Z A0). SDPB x stores off-diagonal Z entries twice. Small
stationarity residuals are reported, never discarded on an unbounded domain.
"""
from itertools import combinations
import json
from pathlib import Path

from flint import arb, arb_mat, ctx


def _stream_pmp(path):
    """Read one matrix at a time; M50 JSON does not need to live wholly in RAM."""
    decoder, marker = json.JSONDecoder(), '"PositiveMatrixWithPrefactorArray"'
    with Path(path).open() as stream:
        buf = ""
        while marker not in buf:
            chunk = stream.read(1 << 20)
            if not chunk:
                raise ValueError("PMP matrix array is missing")
            buf += chunk
        head, buf = buf.split(marker, 1)
        yield json.loads(head.rstrip().rstrip(",") + "}")
        buf = buf.split("[", 1)[1]
        while True:
            buf = buf.lstrip(" \r\n\t,")
            if buf.startswith("]"):
                return
            try:
                block, end = decoder.raw_decode(buf)
            except json.JSONDecodeError:
                chunk = stream.read(1 << 20)
                if not chunk:
                    raise ValueError("Truncated PMP matrix array")
                buf += chunk
                continue
            yield block
            buf = buf[end:]


def _psd_status(matrix):
    n = matrix.nrows()
    if n > 12:
        raise ValueError("This audit is scoped to the small unregularized 2309 blocks")
    inconclusive = False
    # Every principal minor, including exact zero: permits semidefinite matrices.
    for size in range(1, n+1):
        for ids in combinations(range(n), size):
            det = arb_mat([[matrix[i, j] for j in ids] for i in ids]).det()
            if det < 0:
                return "certified_fail"
            if not det >= 0:
                inconclusive = True
    return "inconclusive" if inconclusive else "certified_pass"


def replay_dual(workdir, bits=384):
    """Replay Z, its PSD signs and stationarity in ball arithmetic from raw text."""
    ctx.prec = bits
    root = Path(workdir)
    blocks = _stream_pmp(root / "pmp.json")
    header = next(blocks)
    objective = [arb(v) for v in header["objective"]]
    normalization = [arb(v) for v in header["normalization"]]
    if normalization != [arb(1)] + [arb(0)]*(len(objective)-1):
        raise ValueError("Expected y0=1 normalization")
    total = arb_mat(len(objective), 1)
    signs = []
    for j, block in enumerate(blocks):
        poly = block["polynomials"]
        n, rows = len(poly), []
        for col in range(n):
            for row in range(col+1):
                if any(len(p) != 1 for p in poly[row][col]):
                    raise ValueError("Only degree-0 PMP is supported")
                rows.append([arb(p[0]) for p in poly[row][col]])
        vals = (root / f"out/x_{j}.txt").read_text().splitlines()
        shape = list(map(int, vals[0].split()))
        weights = [arb(x.strip()) for x in vals[1:] if x.strip()]
        if shape != [len(rows), 1] or len(weights) != len(rows):
            raise ValueError(f"Invalid x block {j}")
        total += arb_mat(rows).transpose() * arb_mat([[x] for x in weights])
        Z, cursor = arb_mat(n, n), 0
        for col in range(n):
            for row in range(col+1):
                Z[row, col] = Z[col, row] = weights[cursor] / (1 if row == col else 2)
                cursor += 1
        signs.append(_psd_status(Z))
    residual = [total[i, 0] + objective[i] for i in range(1, len(objective))]
    exact = all(r.is_exact() and r == 0 for r in residual)
    upper = total[0, 0] + objective[0]
    psd = all(s == "certified_pass" for s in signs)
    return {"bits": bits, "n_blocks": len(signs),
            "dual_psd_counts": {s: signs.count(s) for s in
                                ("certified_pass", "certified_fail", "inconclusive")},
            "max_stationarity_residual_upper": str(max(r.abs_upper() for r in residual).upper().fmpq()),
            "stationarity_exact": exact, "upper_bound_midpoint": float(upper.mid()),
            "upper_bound_interval": [str(upper.lower().fmpq()), str(upper.upper().fmpq())],
            "emitted_pmp_bound_certified": psd and exact,
            "scope": "emitted rounded PMP only; no transfer through SVD truncation or source rounding",
            "residual_correction_required": not exact}
