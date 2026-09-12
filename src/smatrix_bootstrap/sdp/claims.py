"""Quantitative verdicts for the paper's claims C1-C8 (task section 3).

Each function returns ``{"verdict": pass|fail|not run, "evidence": ...,
"rows": [...]}``.  Thresholds come from the pre-registration and are never
adjusted here; a claim with missing runs reports "not run", never a pass.
"""
from __future__ import annotations

import numpy as np

from . import constraints as C
from .figures import boundary_points, c1_table, load_csv, phase_comparison

PAIRING = {"tip": "red", "mid": "pink", "ref": "light_pink"}


def fig8_reference() -> dict:
    """The Fig. 8 reference numbers, read off the digitised file rather than
    hard-coded: the chiral-only (green) and chiral+UV (cyan) boundaries and the
    four highlighted points."""
    b = load_csv("figure8_boundary.csv")
    p = load_csv("figure8_selected_points.csv")
    out = {}
    for g, key in (("chiral_only_green", "chiral"), ("gauge_cyan", "uv")):
        m = b["group"] == g
        out[key + "_x_end"] = float(b["f00_s3"][m].max())
        out[key + "_f11_min"] = float(b["f11_s3"][m].min())
    for g in ("chiral_reference_black", "red", "pink", "light_pink"):
        m = p["group"] == g
        if m.sum():
            out[g] = (float(p["f00_s3"][m][0]), float(p["f11_s3"][m][0]))
    return out


def _verdict(ok: bool | None) -> str:
    return "not run" if ok is None else ("pass" if ok else "FAIL")


def c1(records) -> dict:
    pure = [r for r in records if not r["spec"]["chiral"] and not r["spec"]["uv"]
            and r["job"].startswith("dir")
            and r["result"].get("f00_3") is not None
            and r.get("verification", {}).get("unitarity", {}).get("feasible", False)]
    if len(pure) < 24:
        return {"verdict": "not run",
                "evidence": f"{len(pure)}/24 verified-feasible directions"}
    Ms = sorted({r["spec"]["M"] for r in pure})
    Ls = sorted({r["spec"]["L"] for r in pure})
    t = c1_table(pure)
    note = ""
    if Ms != [50] or Ls != [10]:
        note = (f"  NOTE: run at M={Ms}, L={Ls}, not the paper's M=50, L=10; the "
                f"digitised Fig. 3 is at M=50 so the comparison carries a "
                f"resolution offset -- see the resolution ladder.")
    return {"verdict": _verdict(t["pass"]), "rows": t["rows"], "M": Ms, "L": Ls,
            "evidence": "; ".join(f"{r['quantity']} {r['ours']:.4f} vs {r['paper']:.4f} "
                                  f"({r['rel_diff']*100:+.2f}%)" for r in t["rows"]) + note}


def _chiral_sets(records):
    out = {}
    for r in records:
        s = r["spec"]
        if s["chiral"] and not s["uv"] and r["result"].get("f00_3") is not None:
            out.setdefault((s["eps_chi"], s["chi_caliber"]), []).append(r)
    return out


def c2(records) -> dict:
    sets = _chiral_sets(records)
    main = sets.get((C.EPS_CHI_MAIN, "chi-b"))
    if not main:
        return {"verdict": "not run", "evidence": "no eps=2e-3 chi-b sweep"}
    pts = boundary_points(main, want="dir")      # sweep directions only
    if len(pts) == 0:
        return {"verdict": "not run",
                "evidence": "no verified-feasible eps=2e-3 chi-b sweep directions"}
    x_end = float(pts[:, 0].max())
    tgt = fig8_reference()["chiral_x_end"]
    Ms = sorted({r["spec"]["M"] for r in main})
    rows = [{"quantity": "+x end at eps=2e-3", "ours": x_end, "paper": tgt,
             "rel_diff": x_end / tgt - 1.0, "pass": abs(x_end / tgt - 1) <= 0.05}]
    ends = {}
    for (eps, cal), rs in sorted(sets.items()):
        if cal != "chi-b":
            continue
        q = boundary_points(rs, want="dir")
        if len(q):
            ends[eps] = float(q[:, 0].max())
    w = section_width(records, chiral=True, uv=False, eps=C.EPS_CHI_MAIN)
    if w:
        rows.append({"quantity": "x_ref section width", "ours": w["width"],
                     "paper": 7.6e-4, "rel_diff": w["width"] / 7.6e-4 - 1.0,
                     "pass": abs(w["width"] / 7.6e-4 - 1) <= 0.20})
    mono = None
    if len(ends) >= 3:
        e = [ends[k] for k in sorted(ends, reverse=True)]
        mono = all(e[i] >= e[i + 1] - 1e-12 for i in range(len(e) - 1))
        rows.append({"quantity": "+x end monotone in eps", "ours": e, "paper": "decreasing",
                     "rel_diff": None, "pass": mono})
    return {"verdict": _verdict(all(r["pass"] for r in rows) if mono is not None else None),
            "rows": rows, "eps_ends": ends, "M": Ms,
            "evidence": f"+x end {x_end:.6f} vs {tgt:.6f} ({100*(x_end/tgt-1):+.2f}%) at M={Ms}; "
                        f"eps ladder { {k: round(v, 6) for k, v in ends.items()} }"}


def section_width(records, chiral: bool, uv: bool, eps=None) -> dict | None:
    """(hi, lo, width) of the x_ref section, from the section_hi/section_lo jobs."""
    hi = lo = None
    for r in records:
        s = r["spec"]
        if s["chiral"] != chiral or s["uv"] != uv:
            continue
        if eps is not None and abs(s["eps_chi"] - eps) > 1e-12:
            continue
        if r["result"].get("f11_3") is None:
            continue
        if r["job"] == "section_hi":
            hi = r["result"]["f11_3"]
        elif r["job"] == "section_lo":
            lo = r["result"]["f11_3"]
    if hi is None or lo is None:
        return None
    return {"hi": hi, "lo": lo, "width": hi - lo}


def c3(by_eps: dict) -> dict:
    """Subthreshold partial waves and the S0 chiral zero (Fig. 5).

    ``by_eps``: {eps: record["subthreshold"]} for the chiral-only representative
    point at that tolerance.  The paper's colours map to
    eps = 2e-3 (green) / 4e-3 (orange) / 6e-3 (blue); the acceptance is a
    pointwise rms below 8% of f00(3), plus the S0 zero at ~0.425 / ~0.305 / none.
    """
    if not by_eps:
        return {"verdict": "not run", "evidence": "no chiral-only subthreshold curves"}
    ref = load_csv("figure5_subthreshold.csv")
    expect_zero = {0.002: 0.425, 0.004: 0.305, 0.006: None}
    rows = []
    for eps, cur in sorted(by_eps.items()):
        s_ours = np.asarray(cur["s"])
        scale = abs(float(np.interp(3.0, s_ours, np.asarray(cur["S0"]))))
        for wave in ("S0", "S2", "P1"):
            m = (ref["epsilon"] == eps) & (ref["wave"] == wave)
            if m.sum() == 0:
                continue
            xs, ys = ref["s"][m], ref["f"][m]
            o = np.argsort(xs)
            pred = np.interp(xs[o], s_ours, np.asarray(cur[wave]))
            d = pred - ys[o]
            rms = float(np.sqrt((d ** 2).mean()))
            rows.append({"eps": eps, "wave": wave, "rms": rms,
                         "rms_over_f00_3": rms / max(scale, 1e-12),
                         "pass": rms <= 0.08 * max(scale, 1e-12)})
        z = _first_zero(s_ours, np.asarray(cur["S0"]))
        tgt = expect_zero.get(eps, "n/a")
        ok = (z is None) if tgt is None else (z is not None and abs(z - tgt) <= 0.15)
        rows.append({"eps": eps, "wave": "S0 zero", "ours": z, "paper": tgt, "pass": ok})
    return {"verdict": _verdict(all(r["pass"] for r in rows)), "rows": rows,
            "evidence": "; ".join(f"eps={r['eps']}: S0 zero {r.get('ours')}"
                                  for r in rows if r["wave"] == "S0 zero")}


def _first_zero(s, f):
    for i in range(len(s) - 1):
        if f[i] * f[i + 1] < 0:
            t = -f[i] / (f[i + 1] - f[i])
            return float(s[i] + t * (s[i + 1] - s[i]))
    return None


def c5(records) -> dict:
    """Upper boundary shrinks much more than the lower one on the x_ref section."""
    x_ref, _ = C.chiral_reference_point()
    chi = [r for r in records if r["spec"]["chiral"] and not r["spec"]["uv"]]
    uv = [r for r in records if r["spec"]["chiral"] and r["spec"]["uv"]]
    if not chi or not uv:
        return {"verdict": "not run", "evidence": "need both chiral and chiral+UV sweeps"}
    ref = fig8_reference()
    pc, pu = boundary_points(chi, want="dir"), boundary_points(uv, want="dir")
    if len(pc) == 0 or len(pu) == 0:
        return {"verdict": "not run",
                "evidence": "need verified-feasible chiral and chiral+UV sweeps"}
    xe_chi = float(pc[:, 0].max())
    xe_uv = float(pu[:, 0].max())
    rows = [{"quantity": "UV +x end", "ours": xe_uv, "paper": ref["uv_x_end"],
             "rel_diff": xe_uv / ref["uv_x_end"] - 1.0,
             "pass": abs(xe_uv / ref["uv_x_end"] - 1) <= 0.05},
            {"quantity": "chiral +x end", "ours": xe_chi, "paper": ref["chiral_x_end"],
             "rel_diff": xe_chi / ref["chiral_x_end"] - 1.0,
             "pass": abs(xe_chi / ref["chiral_x_end"] - 1) <= 0.05},
            {"quantity": "UV/chiral shrink ratio", "ours": xe_uv / xe_chi,
             "paper": ref["uv_x_end"] / ref["chiral_x_end"],
             "rel_diff": (xe_uv / xe_chi) / (ref["uv_x_end"] / ref["chiral_x_end"]) - 1.0,
             "pass": abs((xe_uv / xe_chi) / (ref["uv_x_end"] / ref["chiral_x_end"]) - 1) <= 0.02}]
    return {"verdict": _verdict(all(r["pass"] for r in rows)), "rows": rows,
            "evidence": "chiral +x end %.5f -> UV %.5f (paper %.5f -> %.5f)"
                        % (xe_chi, xe_uv, ref["chiral_x_end"], ref["uv_x_end"])}


def c6(points: dict) -> dict:
    """rho position at the three representative points (Fig. 9)."""
    if not points:
        return {"verdict": "not run", "evidence": "no representative points"}
    rows, cross = [], {}
    for name, obs in points.items():
        p1 = obs.get("P1", {})
        e = p1.get("crossing_90_GeV")
        cross[name] = None if e is None else 1000.0 * e
        rows.append({"point": name, "crossing_MeV": cross[name],
                     "modulus_peak_MeV": None if p1.get("modulus_peak_GeV") is None
                     else 1000 * p1["modulus_peak_GeV"],
                     "min_eta": p1.get("min_eta_below_1p2GeV")})
    have = [v for v in cross.values() if v is not None]
    if len(have) < len(cross):
        return {"verdict": "FAIL", "rows": rows,
                "evidence": "some representative points have no 90-degree crossing"}
    in_band = all(795.0 <= v <= 845.0 for v in have)
    spread = max(have) - min(have)
    eta_ok = all(r["min_eta"] is not None and r["min_eta"] >= 0.9 for r in rows)
    return {"verdict": _verdict(in_band and spread <= 20.0 and eta_ok), "rows": rows,
            "evidence": f"crossings {['%.0f' % v for v in have]} MeV, spread {spread:.0f} MeV, "
                        f"band [795,845] {in_band}, min eta >= 0.9 {eta_ok}"}


def c4(points: dict) -> dict:
    """Chiral only: S0/S2 match experiment, P1 has no rho (Fig. 7)."""
    if not points:
        return {"verdict": "not run", "evidence": "no chiral-only representative point"}
    rows = []
    for name, obs in points.items():
        d00 = _at(obs.get("S0"), 0.9)
        d11 = _at(obs.get("P1"), 1.2)
        rows.append({"point": name, "delta00_0.9GeV": d00, "delta11_1.2GeV": d11,
                     "pass": (d00 is not None and 85 <= d00 <= 110
                              and d11 is not None and d11 <= 25)})
    # The paper's Fig. 7 is the magenta point "closest to the black dot", which
    # is the `ref` point of task section 2 (max f11 on the f00 = x_ref section)
    # and pairs with the digitised light_pink marker.  The others are reported
    # but do not decide the claim.
    decisive = [r for r in rows if r["point"] == "ref"]
    return {"verdict": _verdict(decisive[0]["pass"] if decisive else None),
            "rows": rows, "decisive_point": "ref",
            "evidence": "; ".join("%s: d00(0.9)=%s, d11(1.2)=%s"
                                  % (r["point"],
                                     None if r["delta00_0.9GeV"] is None
                                     else round(r["delta00_0.9GeV"], 1),
                                     None if r["delta11_1.2GeV"] is None
                                     else round(r["delta11_1.2GeV"], 1))
                                  for r in rows)}


def c7(points: dict) -> dict:
    if not points:
        return {"verdict": "not run", "evidence": "no representative points"}
    rows = []
    for name, obs in points.items():
        d00, d20 = _at(obs.get("S0"), 1.196), _at(obs.get("S2"), 1.196)
        cmp00 = phase_comparison(_series(obs["S0"]), "figure10_s0_phases.csv",
                                 PAIRING.get(name)) if obs.get("S0") else {}
        cmp20 = phase_comparison(_series(obs["S2"]), "figure10_s2_phases.csv",
                                 PAIRING.get(name)) if obs.get("S2") else {}
        rows.append({"point": name, "delta00_1.196": d00, "delta20_1.196": d20,
                     "rms00_deg": cmp00.get("rms_deg"), "rms20_deg": cmp20.get("rms_deg"),
                     "pass": (d00 is not None and 85 <= d00 <= 110
                              and d20 is not None and -40 <= d20 <= -15)})
    # C7 is a statement about all three representative points at once
    return {"verdict": _verdict(all(r["pass"] for r in rows) if len(rows) == 3 else None),
            "rows": rows,
            "evidence": "; ".join("%s: d00=%s, d20=%s"
                                  % (r["point"],
                                     None if r["delta00_1.196"] is None
                                     else round(r["delta00_1.196"], 1),
                                     None if r["delta20_1.196"] is None
                                     else round(r["delta20_1.196"], 1))
                                  for r in rows)}


def c8(by_ml: dict) -> dict:
    """M/L stability of the rho position (Fig. 11)."""
    if len(by_ml) < 5:
        return {"verdict": "not run", "evidence": f"{len(by_ml)}/5 (M,L) configurations"}
    rho = {k: [1000 * p["P1"]["crossing_90_GeV"] for p in v.values()
               if p.get("P1", {}).get("crossing_90_GeV")] for k, v in by_ml.items()}
    fixedM = {k: v for k, v in rho.items() if k[0] == 50}
    fixedL = {k: v for k, v in rho.items() if k[1] == 10}
    allL = [x for v in fixedM.values() for x in v]
    allM = [x for v in fixedL.values() for x in v]
    spread_L = max(allL) - min(allL) if allL else None
    spread_M = max(allM) - min(allM) if allM else None
    ok = (spread_L is not None and spread_L <= 20.0
          and spread_M is not None and 40.0 <= spread_M <= 70.0)
    return {"verdict": _verdict(ok), "rho_MeV": rho,
            "evidence": f"L spread {spread_L} MeV (<=20), M spread {spread_M} MeV (40-70)"}


def _at(obs, e_gev):
    if not obs:
        return None
    E = np.asarray(obs["E_GeV"])
    if e_gev < E.min() or e_gev > E.max():
        return None
    return float(np.interp(e_gev, E, np.asarray(obs["delta_deg"])))


def _series(obs):
    return {"E_GeV": np.asarray(obs["E_GeV"]), "delta_deg": np.asarray(obs["delta_deg"])}
