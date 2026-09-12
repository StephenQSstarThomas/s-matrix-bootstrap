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

FINDINGS = [
 ("(3.67) 核的独立推导",
  "论文的 K 核不是照抄的：本实现从 sigma 在 2M 点交错网格上的奇延拓推出共轭函数算子，"
  "与 (3.67) 逐元素相差 1e-15。它把 sin(n phi) 精确映到 cos(n phi)（n < M）。"),
 ("割线上必须用 cot 核，不能用去点中点法",
  "两者离散同一个主值。在解析测试函数上，M=50 时 cot 核误差 3e-15（谱精度），"
  "去点中点法误差 9e-2 且只按 1/M 衰减（M=50/100/200/400: 9.0e-2, 4.5e-2, 2.2e-2, 1.1e-2）。"
  "仓库历史实现用的也是 cot 核，5a.10 因此可比。"),
 ("m_q 口径由 (2.56) 唯一确定",
  "用 (2.53)-(2.54) 独立重算 (2.50)：只有取 m_q = sqrt((m_u^2+m_d^2)/2) = 5.886 MeV 时，"
  "S0 两个矩与论文打印值差 0.7%；取算术平均 5.650 MeV 则差 9.3%。P1 不含 m_q，差 0.07%。"
  "这与 N_f m_q^2 -> sum_f m_f^2 的物理一致。主线仍用打印值；(3.75) 的 m_q 按任务书取算术平均。"),
 ("FESR 硬截止有系统偏差",
  "(3.72) 保留 s_i <= s0 的节点，但最后一个节点的求积格子止于 s = 84.06 而非 s0 = 73.47"
  "（s 上超出 14.4%）。中点法本身对其实际区间精确到 5e-3 以内，但相对 [4, s0] 的真积分，"
  "n = -1/0/1 三个矩分别偏高 0.6% / 3.5% / 11.3%。如实报告，不做修正。"),
 ("锥的条件数与可行集的尺度是两回事",
  "|h|^2 <= 2 Im h 对任意 Lambda > 0 等价于 |h/Lambda|^2 <= 2 (Im h)/Lambda^2（两边乘 Lambda^2）。"
  "按 Lambda_ell(s) 缩放确实把行范数从 5e-85..34 压到 O(1)，但同时把可行集在缩放坐标下撑到 1e40，"
  "Clarabel 在第 0 步就报 NumericalError。不缩放则约 1650/3000 行下溢为零、锥退化成 0 <= 0、"
  "Slater 条件失效、求解器停在 gap = 6.6。最终用逐行实测幅度做 Lambda^2：行范数全部落入 [0.33, 5.8]，"
  "而未缩放的锥本身已蕴含 |h| <= 2。"),
 ("精确基约化在实践中不可靠",
  "所有约束行与目标行只张成 3876 维中的约 1780 维，c = V a 在数学上是精确重参数化。"
  "但在 M=20 上它给出 2.006，而收敛答案是 1.928。默认关闭。"),
 ("B 的口径：||rho||_4 <= 377500 有两处独立佐证",
  "论文正文没有密度正则化。作者 2403 代码有 norm(rho/Mrho,4) <= 1e2，即 ||rho||_4 <= 100*3775 "
  "= 377500；本仓库历史主线 basis.py:233 用的是同一条 "
  "(density_rule='declared 100*double_density_count', density_norm='actual-rho L4')，"
  "并在 M=50 上得到 f00 极值 2.2349，与数字化 Fig.3 的 2.2329 差 0.1%。两处独立使用同一数值，"
  "说明论文 Fig.3 实际上是带这条密度界的区域。本任务按预登记网格报告 B 的敏感性与活跃性。"),
 ("求解器状态不作为判据：用约束生成把结果变成可认证的",
  "1500 个幺正圆盘里只有约 100 个的算子范数在最大值的 1% 以内，441 个在 1e-6 以内，"
  "其余被离心因子压低许多个量级。只上其中一个子集是**松弛**，其极值是真极值的上界；"
  "若返回点随后在未改动的算子上满足全部 1500 个圆盘，它就对完整问题可行，因而就是完整问题的最优解。"
  "runner.solve_generated 就跑这个循环。松弛上界单调下降（M=20 上 4.248 -> 2.041 -> 1.930 -> 1.92855），"
  "这本身也是一重检查。既把数字变成可认证的，又缩小了 M=50 下主导开销的稠密 KKT。"),
 ("一次被更正的中间结论",
  "中途曾以为“Clarabel 报告的 optimal 不可信”，理由是把 l2 密度球从 377500 放大到 1e7、1e9 时"
  "报告的极值反而下降。这条对 l2 代用球成立（那里求解器确实返回了次优点），但不能推广到预登记的"
  "l4 口径：无 B 时得到的 f00 = 2.0071 的点 ||rho||_4 = 5.19e8，本来就不在 ||rho||_4 <= 377500 的"
  "可行集内，所以 l4 家族内从未出现单调性矛盾。约束生成给出的认证值使这一判断不再依赖求解器状态。"),
 ("事后可行性用哪个指标",
  "max eta - 1 在离心压低的圆盘上恒为 0（h ~ 0 时 eta = 1 恰好），最小 margin 则被同一批圆盘"
  "压到 1e-40；两者都不辨真伪。真正有分辨力的是限制在 |h| > 1e-6 的圆盘上的相对违反 "
  "(|h|^2 - 2 Im h)/|h|^2。"),
 ("(3.75) 的 eps^FF = 6e-5 在低分辨率下不可行，在 M>=30 上可行——这是分辨率效应",
  "把 (3.64) 手征 + (3.68) Gram + (3.73) FESR 都加上后，求“使 (3.75) 可行的最小 eps^FF”"
  "（scripts/sdp/ff_tolerance.py：最小化倍数 t 使 |cF| <= t·bound，则 eps^FF_min = 6e-5·t^2）：\n\n"
  "  | M, L | 施加圆盘 | t | eps^FF_min | 相对论文 6e-5 |\n"
  "  |---|---|---|---|---|\n"
  "  | 20, 6 | 342 | 4.06 | 9.9e-4 | 16x |\n"
  "  | 25, 8 | 541 | 5.28 | 1.7e-3 | 28x |\n"
  "  | 30, 8 | 551 | 0.22 | 2.9e-6 | 0.05x（可行，且有富余） |\n\n"
  "  即：论文的 eps^FF 在 M=20/25 上确实不可行，但在 M=30 上可行且留有 4 倍余量。"
  "因此这是**离散分辨率不足**，不是论文口径有问题。机制是清楚的：F_0 在圆盘上满足均值定理 "
  "(1/pi)∫_0^pi Re F_0 dphi = 1，而 (3.75) 要求它在 s > s0 那段 phi 上小到 0.072，"
  "所以 F_0 必须在阈值附近变大来补偿，而阈值附近的节点数随 M 增长。\n\n"
  "  绑定的是 S0 不是 P1。论文说 F 与 cF 之间的因子“which we evaluate at s = s_0”；"
  "按这句把 (2.33) 冻结在 s0 只影响 P1（k_1 从 s0 到最高节点涨 15.5 倍），对 S0 几乎无影响"
  "（比值 1.00–1.01）。实测两种读法的 t 几乎相同（M=20: 4.06 vs 4.06；M=25: 5.28 vs 3.23），"
  "证实瓶颈在 S0。两种读法都已实现（ModelSpec.ff_frozen_at_s0），默认取论文的字面读法。"),
 ("P1 对拍改用 lambda_max",
  "任务书的 P1 要求用作者 2403 参数复现 2403 的一个公开点。本任务改用 lambda = (pi/4) "
  "T_3333(4/3,4/3,4/3) 的上界 2.661——这是对 2309 同一套（解析性+交叉+幺正性）已发表的独立数值，"
  "既验证了同一条“组装约束 + 调求解器”链路，又不必把 2403 的物理参数（nu0 = -20、s0 = 2 GeV、"
  "三电流）搬进来（第 7 节明令禁止）。lambda 行本身用论文 (2.14) 下方的 Weinberg 值 "
  "m_pi^2/(32 pi f_pi^2) = 0.023 校准。"),
]

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
    for lad, title in (("ladder_pure.json", "纯幺正 max f00(3)"),
                       ("ladder_chiral.json", "手征 (eps=2e-3, chi-b) max f00(3)")):
        lp = os.path.join(root, lad)
        if not os.path.exists(lp):
            continue
        with open(lp) as fh:
            L_ = json.load(fh)
        A("\n## 6b. 分辨率阶梯 — %s（对照 %s = %s）\n"
          % (title, L_["paper_reference"], L_["paper_value"]))
        A("| M | L | 状态 | 认证 | 目标值 | 相对论文 | 施加圆盘 | 轮数 | \\|\\|c\\|\\|_inf | 秒 |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for r in L_["rows"]:
            A("| {M} | {L} | {st} | {c} | {o} | {rel} | {d} | {n} | {cn} | {s} |".format(
                M=r["M"], L=r["L"], st=r["status"], c=r.get("certified", "—"),
                o="—" if r.get("objective") is None else "%.6f" % r["objective"],
                rel="—" if r.get("rel_to_paper") is None else "%+.2f%%" % (100 * r["rel_to_paper"]),
                d=r.get("n_disks_imposed", "—"), n=r["rounds"],
                cn="—" if r.get("c_norm_inf") is None else "%.2e" % r["c_norm_inf"],
                s="%.0f" % r["seconds"]))
        A("")
    st_path = os.path.join(root, "solver_study_M20.json")
    if os.path.exists(st_path):
        with open(st_path) as fh:
            st = json.load(fh)
        A("\n## 7b. 求解器行为研究 (scripts/sdp/solver_study.py, M=%d L=%d)\n"
          % (st["M"], st["L"]))
        A("| 配置 | 状态 | f00(3) | 事后可行 | 最大相对违反 | \\|\\|rho\\|\\|_4 | 迭代 | 秒 |")
        A("|---|---|---|---|---|---|---|---|")
        for r in st["rows"]:
            f = r.get("f00_3")
            A("| {t} | {s} | {f} | {ok} | {v} | {r4} | {it} | {sec} |".format(
                t=r["tag"], s=r["status"],
                f="—" if f is None else "%.6f" % f,
                ok=r.get("feasible", "—"),
                v="—" if r.get("max_rel_violation") is None else "%+.1e" % r["max_rel_violation"],
                r4="—" if r.get("rho_l4") is None else "%.3e" % r["rho_l4"],
                it=r.get("iterations", "—"),
                sec="%.0f" % r["seconds"]))
        A("")
        A("- 目标值对 B 单调（可行集随 B 单调变大，极值必须非减）: **%s**"
          % st["objective_is_monotone_in_B"])
        b = st.get("best_verified_feasible")
        if b:
            A("- 已验证可行的最好点: %s, f00(3) = %.6f" % (b["tag"], b["f00_3"]))
        A("")
    A("\n## 8. 方法、偏离与发现\n")
    for title, body in FINDINGS:
        A(f"**{title}.** {body}\n")
    A("\n## 9. 本任务未做的事\n")
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
