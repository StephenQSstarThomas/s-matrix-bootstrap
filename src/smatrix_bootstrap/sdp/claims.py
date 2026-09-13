"""Quantitative verdicts for the paper's claims C1-C8 (task section 3).

Each function returns ``{"verdict": pass|fail|not run, "evidence": ...,
"rows": [...]}``.  Thresholds come from the pre-registration and are never
adjusted here; a claim with missing runs reports "not run", never a pass.
"""
from __future__ import annotations

import numpy as np

from . import constraints as C
from .figures import accepted_support, boundary_points, c1_table, load_csv, phase_comparison
from .assembly import recorded_basis_key

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


def c1(records, ladder: dict | None = None) -> dict:
    pure = [r for r in records if not r["spec"]["chiral"] and not r["spec"]["uv"]
            and r["job"].startswith("dir")
            and r["result"].get("f00_3") is not None
            and r.get("verification", {}).get("unitarity", {}).get("feasible", False)]
    if len(pure) < 24:
        ev = f"{len(pure)}/24 verified-feasible directions"
        if ladder:
            ok = [r for r in ladder.get("rows", []) if r.get("objective") is not None]
            if ok:
                ev += ("; the +x tip alone, from the resolution ladder: "
                       + ", ".join("M=%d %.4f (%+.1f%%)"
                                   % (r["M"], r["objective"], 100 * r["rel_to_paper"])
                                   for r in ok))
        return {"verdict": "not run", "evidence": ev}
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
        if s["chiral"] and not s["uv"] and accepted_support(r) and r["result"].get("f00_3") is not None:
            out.setdefault((s["eps_chi"], s["chi_caliber"],recorded_basis_key(r)), []).append(r)
    return out


def c2(records, eps_ladder: dict | None = None) -> dict:
    """All registered epsilon tips and the main x_ref width in one finite model."""
    records = [r for r in records if accepted_support(r) and r['spec']['chiral']
               and not r['spec']['uv']]
    contracts = {(recorded_basis_key(r),tuple(sorted((k,repr(v)) for k,v in r['spec'].items()
                 if k not in ('eps_chi','tag','disk_mask')))) for r in records}
    if len(contracts)!=1:
        return {'verdict':'not run','evidence':'one common accepted chiral source contract required'}
    ends = {}
    for r in records:
        d = np.asarray(r['result'].get('direction'),dtype=float)
        if r.get('fix_f00') is None and d.shape==(2,) and d[0]>0 and d[1]==0:
            ends[r['spec']['eps_chi']] = r['result']['f00_3']
    w = section_width(records,True,False,C.EPS_CHI_MAIN)
    if not set(C.EPS_CHI_GRID).issubset(ends) or w is None:
        return {'verdict':'not run','evidence':'six accepted epsilon tips and x_ref upper/lower supports required',
                'eps_ends':ends}
    x_end,tgt = ends[C.EPS_CHI_MAIN],fig8_reference()['chiral_x_end']
    e = [ends[k] for k in sorted(C.EPS_CHI_GRID,reverse=True)]
    rows = [{'quantity':'+x end at eps=2e-3','ours':x_end,'paper':tgt,
             'pass':abs(x_end/tgt-1)<=.05},
            {'quantity':'x_ref section width','ours':w['width'],'paper':7.6e-4,
             'pass':abs(w['width']/7.6e-4-1)<=.20},
            {'quantity':'+x end monotone in eps','ours':e,
             'pass':all(e[i]>=e[i+1]-1e-12 for i in range(len(e)-1))}]
    return {'verdict':_verdict(all(r['pass'] for r in rows)),'rows':rows,'eps_ends':ends,
            'evidence':'same-contract SDPB source supports; external ladder summaries alone are insufficient'}


def section_width(records, chiral: bool, uv: bool, eps=None) -> dict | None:
    """(hi, lo, width) of the x_ref section, from the section_hi/section_lo jobs."""
    hi = lo = None
    contracts = set()
    for r in records:
        s = r["spec"]
        if s["chiral"] != chiral or s["uv"] != uv:
            continue
        if eps is not None and abs(s["eps_chi"] - eps) > 1e-12:
            continue
        if r["result"].get("f11_3") is None:
            continue
        if r['job'] not in ('section_hi','section_lo') or not accepted_support(r):
            continue
        xr = C.chiral_reference_point()[0]
        direction = [0.,1.] if r['job']=='section_hi' else [0.,-1.]
        if r.get('fix_f00') is None or abs(r['fix_f00']-xr)>1e-12 or r['result'].get('direction')!=direction:
            continue
        contracts.add((recorded_basis_key(r),tuple(sorted((k,repr(v)) for k,v in s.items() if k not in ('tag','disk_mask')))))
        if r["job"] == "section_hi":
            hi = r["result"]["f11_3"]
        elif r["job"] == "section_lo":
            lo = r["result"]["f11_3"]
    if hi is None or lo is None or len(contracts)!=1:
        return None
    return {"hi": hi, "lo": lo, "width": hi - lo}


def c3(by_eps: dict) -> dict:
    """Subthreshold partial waves and the S0 chiral zero (Fig. 5).

    ``by_eps``: {eps: record["subthreshold"]} for the chiral-only representative
    point at that tolerance.  The paper's colours map to
    eps = 2e-3 (green) / 4e-3 (orange) / 6e-3 (blue); the acceptance is a
    RMS against the same-epsilon paper curve below 8% of f00(3), not a
    linearity budget, plus the S0 zero at ~0.425 / ~0.305 / none.
    """
    if set(by_eps) != {0.002,0.004,0.006}:
        return {"verdict": "not run", "evidence": "no chiral-only subthreshold curves"}
    if any(not accepted_support(r) or not r['spec']['chiral'] or r['spec']['uv'] for r in by_eps.values()):
        return {'verdict':'not run','evidence':'accepted chiral-only source amplitudes required'}
    identities = {(tuple(sorted((k,repr(v)) for k,v in r['spec'].items()
                   if k not in ('eps_chi','tag','disk_mask'))),tuple(r.get('direction',r.get('result',{}).get('direction',[]))),
                   r.get('fix_f00'),recorded_basis_key(r)) for r in by_eps.values()}
    if len(identities)!=1 or any(r['spec']['eps_chi']!=eps for eps,r in by_eps.items()):
        return {'verdict':'not run','evidence':'same source model and representative objective across epsilon required'}
    ref = load_csv("figure5_subthreshold.csv")
    expect_zero = {0.002: 0.425, 0.004: 0.305, 0.006: None}
    rows = []
    for eps, record in sorted(by_eps.items()):
        cur = record.get('subthreshold',record)
        from .observables import subthreshold_coverage
        missing = subthreshold_coverage(cur,ref['s'][ref['epsilon'] == eps])
        if missing:
            return {'verdict':'not run','evidence':'C3 coverage: '+missing,'epsilon':eps}
        exact = record.get('verification',{}).get('f00_3',cur.get('f00_3'))
        if exact is None:
            return {'verdict':'not run','evidence':'C3 needs directly evaluated f00(3), not display interpolation'}
        s_ours = np.asarray(cur["s"])
        scale = abs(float(exact))
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
    records = [r for r in records if accepted_support(r) and r['spec']['chiral']]
    contracts = {(recorded_basis_key(r),tuple(sorted((k,repr(v)) for k,v in r['spec'].items()
                 if k not in ('uv','tag','disk_mask','uv_parts','sr_caliber','eps_ff','m_q','ff_frozen_at_s0')))) for r in records}
    chi = [r for r in records if not r['spec']['uv']]
    uv = [r for r in records if r['spec']['uv']]
    if any(set(r['spec']['uv_parts'])!={'gram','fesr','ff'} for r in uv):
        return {'verdict':'not run','evidence':'C5 requires full UV constraints'}
    uv_contracts = {tuple(sorted((k,repr(v)) for k,v in r['spec'].items()
                    if k not in ('tag','disk_mask'))) for r in uv}
    if not chi or not uv or len(contracts)!=1 or len(uv_contracts)!=1:
        return {"verdict": "not run", "evidence": "need both chiral and chiral+UV sweeps"}
    ref = fig8_reference()
    pc, pu = (boundary_points([r for r in sector if r.get('fix_f00') is None and r['result'].get('direction')==[1.,0.]])
              for sector in (chi,uv))
    if len(pc) == 0 or len(pu) == 0:
        return {"verdict": "not run",
                "evidence": "need verified-feasible chiral and chiral+UV sweeps"}
    xe_chi = float(pc[:, 0].max())
    xe_uv = float(pu[:, 0].max())
    rows = [{"quantity": "UV +x end", "ours": xe_uv, "paper": ref["uv_x_end"],
             "rel_diff": xe_uv / ref["uv_x_end"] - 1.0,
             "pass": abs(xe_uv / ref["uv_x_end"] - 1) <= 0.05}]
    sc = section_width(chi,True,False,C.EPS_CHI_MAIN)
    su = section_width(uv,True,True,C.EPS_CHI_MAIN)
    if sc and su:
        upper,lower = sc['hi']-su['hi'],su['lo']-sc['lo']
        rows.append({'quantity':'xref asymmetric shrink','upper':upper,'lower':lower,
                     'pass':abs(upper/2.3e-4-1)<=.25 and 0<=lower<=1e-4 and upper>=4*lower})
    return {"verdict": _verdict(all(r["pass"] for r in rows) if sc and su else None), "rows": rows,
            "evidence": "chiral +x end %.5f -> UV %.5f (paper %.5f -> %.5f)"
                        % (xe_chi, xe_uv, ref["chiral_x_end"], ref["uv_x_end"])}


def c6(points: dict) -> dict:
    """rho position at the three representative points (Fig. 9)."""
    if set(points) != set(PAIRING):
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
    # The registered section ref is a declared proxy, not proof of the paper's
    # nearest boundary point. Other representatives do not decide this gate.
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
                              and d20 is not None and -40 <= d20 <= -15
                              and cmp00.get("rms_deg") is not None and cmp00["rms_deg"] <= 10
                              and cmp20.get("rms_deg") is not None and cmp20["rms_deg"] <= 10)})
    # C7 is a statement about all three representative points at once
    return {"verdict": _verdict(all(r["pass"] for r in rows) if set(points)==set(PAIRING) else None),
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
    required = {(50,8),(50,10),(50,12),(45,10),(60,10)}
    if not required.issubset(by_ml):
        return {"verdict": "not run", "evidence": f"{len(by_ml)}/5 (M,L) configurations"}
    by_ml = {k:by_ml[k] for k in required}
    names = set(by_ml[(50,10)])
    if not names or not names.issubset(PAIRING) or any(set(v)!=names for v in by_ml.values()):
        return {'verdict':'not run','evidence':'same named representatives required across all five configurations'}
    if any(p.get('P1',{}).get('crossing_90_GeV') is None for v in by_ml.values() for p in v.values()):
        return {'verdict':'FAIL','evidence':'a required configuration has no rho crossing'}
    rho = {k: [1000 * v[name]['P1']['crossing_90_GeV'] for name in sorted(names)] for k,v in by_ml.items()}
    fixedM = {k: v for k, v in rho.items() if k[0] == 50}
    fixedL = {k: v for k, v in rho.items() if k[1] == 10}
    spread_L = [max(v[i] for v in fixedM.values())-min(v[i] for v in fixedM.values()) for i in range(len(names))]
    spread_M = [max(v[i] for v in fixedL.values())-min(v[i] for v in fixedL.values()) for i in range(len(names))]
    ok = all(v<=20. for v in spread_L) and all(40.<=v<=70. for v in spread_M)
    ok = ok and all(a<b<c for a,b,c in zip(rho[(60,10)],rho[(45,10)],rho[(50,10)]))
    s0 = [_at(p.get('S0'),1.) for v in by_ml.values() for p in v.values()]
    if any(d is None for d in s0):
        return {'verdict':'not run','evidence':'all required S0 values at 1 GeV are needed'}
    ok = ok and all(d is not None and 85<=d<=105 for d in s0)
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
