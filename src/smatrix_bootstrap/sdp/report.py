"""Assemble REPORT_SDP_ZH.md and manifest.json from the report.json artefacts.

Every number in the report is traced back to the ``report.json`` that produced
it, so the manifest is the recomputation index the task asks for.  Claims that
were not run are printed as "not run", never as a pass.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np

from . import constraints as C
from .figures import ENERGY_AXIS_NOTE, c1_table, load_reports

CLAIMS = {
    "C1": "Pure-unitarity region shape (Fig. 3)",
    "C2": "Chiral constraints collapse the region onto f11 = -f00/15 (Fig. 4)",
    "C3": "Subthreshold partial waves near-linear, S0 chiral zero moves (Fig. 5)",
    "C4": "Chiral only: S0/S2 agree with experiment, P1 has no rho (Fig. 7)",
    "C5": "FESR+FF shrink the upper boundary, not the lower (Fig. 8)",
    "C6": "rho appears in P1 at all three representative points (Fig. 9)",
    "C7": "S0/S2 coincide at low energy and spread at high energy (Fig. 10)",
    "C8": "M/L stability (Fig. 11, Appendix A)",
}

NOT_DONE = [
    "continuous unitarity was not verified: unitarity is imposed only on the "
    "M x L x 3 collocation points of (3.61), never between them",
    "no pole was proven: the rho is read off a 90-degree crossing of the "
    "interpolated phase shift, which is not an analytic continuation to the "
    "second sheet",
    "no M -> infinity convergence was established: five (M, L) pairs are a "
    "sensitivity check, not a limit",
    "the double spectral densities carry no Mandelstam support restriction, "
    "because the paper states none",
    "the FESR cutoff bias of section 5a.8 is reported, not corrected",
]


def collect(root: str) -> dict:
    recs = load_reports(root)
    by_stage: dict[str, list] = {}
    for r in recs:
        by_stage.setdefault(r["_file"].split(os.sep)[0], []).append(r)
    return {"records": recs, "by_stage": by_stage}


def solver_stats(recs) -> dict:
    ok = [r for r in recs if r["result"].get("status") in ("optimal", "optimal_inaccurate")]
    it = [r["result"].get("iterations") for r in ok if r["result"].get("iterations")]
    se = [r["result"].get("seconds") for r in ok if r["result"].get("seconds")]
    status: dict[str, int] = {}
    for r in recs:
        status[r["result"].get("status", "?")] = status.get(r["result"].get("status", "?"), 0) + 1
    return {"n_solves": len(recs), "status_counts": status,
            "median_iterations": float(np.median(it)) if it else None,
            "median_seconds": float(np.median(se)) if se else None,
            "total_seconds": float(sum(se)) if se else 0.0}


def b_activity(recs) -> dict:
    rows = []
    for r in recs:
        v = r.get("verification", {})
        if "B" in v:
            rows.append({"job": r["job"], "file": r["_file"], "B": v["B"],
                         "B_norm": v["B_norm"], "rho_l2": v["rho_l2"],
                         "rho_l4": v["rho_l4"], "active": v["B_active"]})
    return {"n": len(rows), "any_active": any(r["active"] for r in rows), "rows": rows}


def unitarity_worst(recs) -> dict:
    worst, where = -1.0, None
    for r in recs:
        u = r.get("verification", {}).get("unitarity")
        if not u:
            continue
        if u["max_eta_minus_1"] > worst:
            worst, where = u["max_eta_minus_1"], {"job": r["job"], "file": r["_file"],
                                                  "wave": u["max_eta_wave"]}
    return {"max_eta_minus_1": worst, "at": where}


def manifest(root: str, data: dict) -> dict:
    out = {"generated": datetime.now(timezone.utc).isoformat(), "root": root,
           "entries": []}
    for r in data["records"]:
        out["entries"].append({
            "job": r["job"], "report_json": r["_file"],
            "M": r["spec"]["M"], "L": r["spec"]["L"],
            "chiral": r["spec"]["chiral"], "uv": r["spec"]["uv"],
            "chi_caliber": r["spec"]["chi_caliber"], "sr_caliber": r["spec"]["sr_caliber"],
            "eps_chi": r["spec"]["eps_chi"], "B": r["spec"]["B"],
            "status": r["result"].get("status"),
            "objective": r["result"].get("objective"),
            "f00_3": r["result"].get("f00_3"), "f11_3": r["result"].get("f11_3"),
            "solution_sha256": r["result"].get("solution_sha256")})
    return out


def render(root: str, verdicts: dict | None = None) -> str:
    data = collect(root)
    recs = data["records"]
    verdicts = verdicts or {}
    pure = [r for r in recs if not r["spec"]["chiral"] and not r["spec"]["uv"]]
    L = []
    A = L.append
    A("# REPORT_SDP_ZH — He–Kruczenski 2309.12402v3, SDP 直解路线\n")
    A(f"生成时间 (UTC): {datetime.now(timezone.utc).isoformat()}\n")
    A(f"结果根目录: `{root}`\n")
    A("\n## 1. 核心 claim 判定表\n")
    A("| # | claim | 判定 | 依据 |")
    A("|---|---|---|---|")
    for k, v in CLAIMS.items():
        d = verdicts.get(k, {})
        A(f"| {k} | {v} | {d.get('verdict', '未运行')} | {d.get('evidence', '—')} |")
    A("\n## 2. 求解器统计\n")
    st = solver_stats(recs)
    A(f"- 求解次数: {st['n_solves']}")
    A(f"- 状态分布: {st['status_counts']}")
    A(f"- 迭代中位数: {st['median_iterations']}, 单次秒数中位数: {st['median_seconds']}")
    A(f"- 机器时间合计: {st['total_seconds']:.0f} s")
    A("\n## 3. 密度正则化 B 的活跃性\n")
    b = b_activity(recs)
    A(f"- 带 B 的求解: {b['n']}；**任一活跃: {b['any_active']}**")
    if b["rows"]:
        A("\n| job | B | 范数 | \\|\\|rho\\|\\|_2 | \\|\\|rho\\|\\|_4 | 活跃 |")
        A("|---|---|---|---|---|---|")
        for r in b["rows"][:20]:
            A(f"| {r['job']} | {r['B']:.3e} | {r['B_norm']} | {r['rho_l2']:.3e} "
              f"| {r['rho_l4']:.3e} | {r['active']} |")
    A("\n## 4. 幺正性事后复验（未经任何重缩放的原式）\n")
    u = unitarity_worst(recs)
    A(f"- 全部解中最大 `max eta - 1` = {u['max_eta_minus_1']:.3e}，位置: {u['at']}")
    A("\n## 5. FESR 目标值独立重算 (5a.9)\n")
    au = C.fesr_target_audit()
    A(f"- m_q 算术平均 = {au['m_q_mean_MeV']:.3f} MeV，均方根 = {au['m_q_rms_MeV']:.3f} MeV")
    A("\n| 波 | n | 论文 (2.56) raw | 重算 (m_q 平均) | 比 | 重算 (m_q 均方根) | 比 |")
    A("|---|---|---|---|---|---|---|")
    for r in au["rows"]:
        A(f"| {r['wave']} | {r['n']} | {r['printed_raw']:.6e} | {r['recomputed_mq_mean']:.6e} "
          f"| {r['ratio_mean']:.4f} | {r['recomputed_mq_rms']:.6e} | {r['ratio_rms']:.4f} |")
    A("\n论文打印值对应 m_q 取均方根；取算术平均则 S0 两个矩差 9.3%。主线用打印值。\n")
    if pure:
        A("\n## 6. C1: 纯幺正区域与数字化 Fig.3 对照\n")
        try:
            t = c1_table(pure)
            A("| 量 | 本工作 | 论文 (数字化) | 相对差 | 通过 (±2%) |")
            A("|---|---|---|---|---|")
            for r in t["rows"]:
                A(f"| {r['quantity']} | {r['ours']:.6f} | {r['paper']:.6f} "
                  f"| {r['rel_diff']*100:+.2f}% | {r['pass']} |")
        except Exception as exc:
            A(f"(C1 表未能生成: {exc})")
    A(f"\n## 7. 能量轴口径\n\n- {ENERGY_AXIS_NOTE}\n")
    A("\n## 8. 本任务未做的事\n")
    for s in NOT_DONE:
        A(f"- {s}")
    A("")
    return "\n".join(L)


def build(root: str, out_md: str, verdicts: dict | None = None) -> None:
    data = collect(root)
    with open(os.path.join(root, "manifest.json"), "w") as fh:
        json.dump(manifest(root, data), fh, indent=1, default=float)
    with open(out_md, "w") as fh:
        fh.write(render(root, verdicts))
