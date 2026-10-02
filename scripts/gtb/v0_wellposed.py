"""V0 for P_W kernel exports (WELLPOSED_ZH.md section 5).  Usage: python scripts/gtb/v0_wellposed.py BASE_L20 NEW [OUT_JSON]
- every wave-indexed matrix of NEW (L >= 20) restricted to l < 20 equals BASE string for string;
- intS0/intP1/intD0 of NEW (fesr_norm 'cell') equal BASE (fesr_norm 'node') times (x_n0/x_b)^(n + n_l), to 1e-35;
- every other matrix is identical.
"""
import json
import sys
from decimal import Decimal, getcontext
from pathlib import Path

from smatrix_bootstrap.gtb.kernels import load

getcontext().prec = 60
WAVE = {"Re_sig": "node", "Re_rho": "node", "Im_sig": "node", "Im_rho": "node", "Lam": "node",
        "fchi_sig": "chi", "fchi_rho": "chi", "fF_sig": "wave", "fF_rho": "wave"}
PER = {"node": None, "chi": 4, "wave": 1}
MOMENTS = {"intS0": ((0, 1, 2), 2), "intP1": ((-1, 0, 1), 2), "intD0": ((-2, -1, 0), 3)}


def rows_of(ks, name):
    shape = ks.manifest["matrices"][name]["shape"]
    flat = ks.exact(name)
    c = shape[1] if len(shape) > 1 else 1
    return [flat[i * c:(i + 1) * c] for i in range(shape[0])]


def main(base_dir, new_dir):
    base, new = load(base_dir), load(new_dir)
    M, nb, nn = int(base.manifest["M"]), int(base.manifest["L"]) // 2, int(new.manifest["L"]) // 2
    out = {}
    for name, kind in WAVE.items():
        per = M if kind == "node" else PER[kind]
        rb, rn = rows_of(base, name), rows_of(new, name)
        bad = 0
        for blk in range(3):
            for w in range(nb):
                for j in range(per):
                    bad += rb[(blk * nb + w) * per + j] != rn[(blk * nn + w) * per + j]
        out[name] = {"rows_compared": 3 * nb * per, "mismatches": bad}
    x = [Decimal(t) for t in base.exact("x")]
    n0 = int(base.manifest["n0"])
    from mpmath import cos, mp, pi
    mp.dps = 60
    nu0 = int(base.manifest["nu0"])
    xb = Decimal(str(nu0 + (8 - 2 * nu0) / (1 + cos(n0 * pi / M))))
    for name, (ns, nl) in MOMENTS.items():
        rb, rn = rows_of(base, name), rows_of(new, name)
        worst = Decimal(0)
        for i, n in enumerate(ns):
            f = (x[n0 - 1] / xb) ** (n + nl)
            for a, b in zip(rb[i], rn[i]):
                a, b = Decimal(a), Decimal(b)
                if a != 0:
                    worst = max(worst, abs(b - a * f) / abs(a))
        out[name] = {"max_rel_vs_scaled": float(worst)}
    for name in base.manifest["matrices"]:
        if name in WAVE or name in MOMENTS:
            continue
        out[name] = {"identical": base.exact(name) == new.exact(name)}
    ok = all(v.get("mismatches", 0) == 0 and v.get("identical", True) and v.get("max_rel_vs_scaled", 0) < 1e-35
             for v in out.values())
    return {"base": str(base_dir), "new": str(new_dir), "fesr_norm": new.manifest.get("fesr_norm"), "pass": ok,
            "checks": out}


if __name__ == "__main__":
    res = main(sys.argv[1], sys.argv[2])
    if len(sys.argv) > 3:
        Path(sys.argv[3]).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({"pass": res["pass"], **{k: v for k, v in res["checks"].items() if k.startswith("int")}}))
    sys.exit(0 if res["pass"] else 1)
