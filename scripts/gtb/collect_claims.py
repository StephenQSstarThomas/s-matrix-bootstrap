"""Collect every 2403 claim verdict from certified leaves into one JSON (missing inputs give 'MISSING').

Usage: python scripts/gtb/collect_claims.py RUNS_ROOT OUT_JSON [--variant A|B] [--b-root DIR]
Layout (Phase B, PREREGISTRATION amendments 4-6): RUNS_ROOT/trial_M50/inner{44,46,48}r (C1, initial A),
phaseB/A/watson{ni}i (Watson A), phaseB/C2/F1{max,min}_k{k} and F1{max,min}_tree (C2), phaseB/var/* (C3, C22-C27).
"""
import json
import sys
from pathlib import Path

from smatrix_bootstrap.gtb import appclaims as AC
from smatrix_bootstrap.gtb import claims as K

FIG3 = Path("/playpen1/shiqiu/gtb2403_prep/paper/fig3_shape_dots_extracted.txt")


def _json(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def _v1(leaf):
    """V1-certified F1 of a leaf, or None (missing, uncertified, or not optimal)."""
    c = _json(Path(leaf) / "certificate.json")
    if not c or not c["lower"]["certified"]:
        return None
    return float(c["lower"]["value"])


def collect(root: Path, a_init: dict, a_wat: dict) -> dict:
    out = {}
    c1 = {ni: {"lower": _v1(root / f"trial_M50/inner{ni}r"), "paper": K.C1_PAPER[ni]} for ni in (44, 46, 48)}
    for ni, r in c1.items():
        r["verdict"] = "MISSING" if r["lower"] is None else K.c1_verdict(r["lower"], float("inf"), r["paper"])
    vs = {r["verdict"] for r in c1.values()}
    out["C1"] = {"verdict": "FAIL" if "FAIL" in vs else "PASS" if vs == {"PASS"} else "INCONCLUSIVE", "points": c1}
    c2dir = root / "phaseB/C2"
    ours = {}
    for k in range(2, 51, 4):
        mx = _v1(c2dir / f"F1max_k{k}") if k != 46 else _v1(root / "trial_M50/inner46r")
        ours[k] = {"max": mx, "min": _v1(c2dir / f"F1min_k{k}")}
    tree = {"max": _v1(c2dir / "F1max_tree"), "min": _v1(c2dir / "F1min_tree")}
    out["C2"] = AC.c2(ours, tree, FIG3) if all(v["max"] is not None and v["min"] is not None for v in ours.values()) \
        else {"verdict": "MISSING", "have": {k: v for k, v in ours.items()}}
    var = root / "phaseB/var"
    # exclusion of the tree dot at eps_chi = 1e-3: a V2 bound (relax encoding, contains the l4 problem) on the side
    # that excludes it: max F1 < F1_tree, or min F1 > F1_tree (the F1min leaf's objective is -F1)
    cmax, cmin = _json(var / "C3_max_tree_e1/certificate.json"), _json(var / "C3_min_tree_e1/certificate.json")
    excl = bool(cmax and cmax["upper"]["certified"] and float(cmax["upper"]["value"]) < AC.F1_TREE) or bool(
        cmin and cmin["upper"]["certified"] and -float(cmin["upper"]["value"]) > AC.F1_TREE)
    out["C3"] = AC.c3(tree, {"certified_infeasible": excl})
    if a_wat:
        out.update(K.phase_claims(a_wat))
        out["C21"] = K.c21(a_init, a_wat)
        if 44 in a_wat and 48 in a_wat:
            out.update(K.table_claims(a_wat[44], a_wat[48]))
        if 48 in a_wat:
            out["C20"] = K.c20(a_wat[48])
    an = lambda name: _json(var / f"{name}_watson/analysis.json")
    base = a_wat.get(46)
    if base:
        out["C22"] = {n: (AC.c22(base, an(n)) if an(n) else {"verdict": "MISSING"})
                      for n in ("C22_M60", "C22_L24", "C22_NB46", "C22_Mreg500")}
        r24 = [an("C24_r11_1.5"), base, an("C24_r11_2.5")]
        out["C24"] = AC.c24(r24) if all(r24) else {"verdict": "MISSING"}
        r25 = [an("C25_nsr4"), an("C25_nsr5")]
        out["C25"] = AC.c25(base, r25) if all(r25) else {"verdict": "MISSING"}
    r27 = (an("C27_noP1SR"), an("C27_noD0SR"))
    out["C27"] = AC.c27(*r27) if all(r27) else {"verdict": "MISSING"}
    caps_ok = _v1(var / "C23_minimal_caps") is not None
    inf = [_json(var / f"C23_{n}/infeasibility.json") for n in ("r11_1", "r02_1")]
    certified = all(i and i.get("certified_infeasible") for i in inf)
    out["C23"] = {"verdict": "PASS" if caps_ok and certified else "FAIL" if not caps_ok and _json(
        var / "C23_minimal_caps/certificate.json") else "INCONCLUSIVE", "minimal_caps_feasible": caps_ok,
        "r_equal_1_infeasible": [bool(i and i.get("certified_infeasible")) for i in inf]}
    return out


if __name__ == "__main__":
    root, dest = Path(sys.argv[1]), Path(sys.argv[2])
    variant = sys.argv[sys.argv.index("--variant") + 1] if "--variant" in sys.argv else "A"
    wroot = Path(sys.argv[sys.argv.index("--b-root") + 1]) if "--b-root" in sys.argv else root / "phaseB/A"
    suffix = "i" if variant == "A" else "b"
    a_init = {ni: _json(root / f"trial_M50/inner{ni}r/analysis.json") for ni in (44, 46, 48)}
    a_wat = {ni: d for ni in (44, 46, 48) if (d := _json(wroot / f"watson{ni}{suffix}/analysis.json"))}
    a_init = {ni: d for ni, d in a_init.items() if d and ni in a_wat}
    res = collect(root, a_init, a_wat)
    dest.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: v.get("verdict") if isinstance(v, dict) and "verdict" in v else
                      {kk: vv.get("verdict") for kk, vv in v.items()} if isinstance(v, dict) else v
                      for k, v in res.items()}, indent=1))
