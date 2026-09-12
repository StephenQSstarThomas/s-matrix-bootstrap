#!/usr/bin/env python
"""Collect every report.json, compute the C1-C8 verdicts, and write the report.

Usage: python scripts/sdp/finalize.py --root <run dir> --repo-out <dir in repo>

Large artefacts (solution vectors) stay on the scratch disk; only the JSON
records, the tables, the manifest and the figures are copied into the
repository, so nothing big lands under /home.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil

from smatrix_bootstrap.sdp import claims, report
from smatrix_bootstrap.sdp.figures import build_all


def representative_points(records, chi=None, sr=None, uv=True, eps=None, M=None):
    """{'tip'|'ref'|'mid': observables} for one caliber and one resolution.

    Mixing eps^chi values or resolutions would silently compare curves from
    different problems, so both are pinned; when M is not given, the largest
    resolution that produced points is used.
    """
    cand = []
    for r in records:
        sp = r["spec"]
        if r["job"] not in ("tip", "ref", "mid") or sp["uv"] != uv:
            continue
        if chi and sp["chi_caliber"] != chi:
            continue
        if sr and sp["sr_caliber"] != sr:
            continue
        if eps is not None and abs(sp["eps_chi"] - eps) > 1e-12:
            continue
        if "observables" in r:
            cand.append(r)
    if not cand:
        return {}
    use_M = M if M is not None else max(r["spec"]["M"] for r in cand)
    return {r["job"]: r["observables"] for r in cand if r["spec"]["M"] == use_M}


def by_ml(records):
    out = {}
    for r in records:
        s = r["spec"]
        if r["job"] in ("tip", "ref", "mid") and s["uv"] and "observables" in r:
            out.setdefault((s["M"], s["L"]), {})[r["job"]] = r["observables"]
    return {k: v for k, v in out.items() if len(v) == 3}


def subthreshold_by_eps(records):
    """Chiral-only `ref` point curves keyed by eps (the Fig. 5 colours)."""
    out = {}
    for r in records:
        s = r["spec"]
        if (r["job"] == "ref" and s["chiral"] and not s["uv"]
                and s["chi_caliber"] == "chi-b" and "subthreshold" in r):
            k = round(s["eps_chi"], 12)
            # keep the highest resolution available for each tolerance
            if k not in out or s["M"] >= out[k][0]:
                out[k] = (s["M"], r["subthreshold"])
    return out


def verdicts(records) -> dict:
    from smatrix_bootstrap.sdp.constraints import EPS_CHI_MAIN
    main = representative_points(records, "chi-b", "SR-b", uv=True, eps=EPS_CHI_MAIN)
    chiral_only = representative_points(records, "chi-b", None, uv=False,
                                        eps=EPS_CHI_MAIN)
    ml = by_ml(records)
    return {"C1": claims.c1(records), "C2": claims.c2(records),
            "C3": claims.c3({k: v[1] for k, v in subthreshold_by_eps(records).items()}),
            "C4": claims.c4(chiral_only), "C5": claims.c5(records),
            "C6": claims.c6(main), "C7": claims.c7(main), "C8": claims.c8(ml)}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--repo-out", default=None)
    p.add_argument("--md", default=None)
    a = p.parse_args()

    data = report.collect(a.root)
    v = verdicts(data["records"])
    with open(os.path.join(a.root, "verdicts.json"), "w") as fh:
        json.dump(v, fh, indent=1, default=float)
    try:
        build_all(a.root)
    except Exception as exc:
        print("figure/table build skipped:", exc)
    md = a.md or os.path.join(a.root, "REPORT_SDP_ZH.md")
    report.build(a.root, md, v)
    print("wrote", md)
    for k, d in v.items():
        print(f"  {k}: {d['verdict']:<9s} {d.get('evidence', '')[:110]}")

    if a.repo_out:
        os.makedirs(a.repo_out, exist_ok=True)
        for name in ("manifest.json", "verdicts.json", "tables.json",
                     "preregistration.json"):
            src = os.path.join(a.root, name)
            if os.path.exists(src):
                shutil.copy2(src, a.repo_out)
        for dirpath, _, files in os.walk(a.root):
            for f in files:
                if f == "report.json":
                    rel = os.path.relpath(os.path.join(dirpath, f), a.root)
                    dst = os.path.join(a.repo_out, "reports", rel)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copy2(os.path.join(dirpath, f), dst)
        figs = os.path.join(a.root, "figures")
        if os.path.isdir(figs):
            shutil.copytree(figs, os.path.join(a.repo_out, "figures"),
                            dirs_exist_ok=True)
        print("copied artefacts to", a.repo_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
