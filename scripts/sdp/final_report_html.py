"""Build the final reproduction ledger (HTML artifact) from the accepted leaves, the claim JSONs and the comparison PNGs.

Usage: python scripts/sdp/final_report_html.py --root RESULTS_ROOT --figures RESULTS_ROOT/final_figures --out FILE.html
"""
from __future__ import annotations

import argparse
import base64
import glob
import html
import json
from pathlib import Path

import numpy as np
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from authors_comparison_en import gather_watson, matched_pair_summary  # noqa: E402


def J(path):
    p = Path(path)
    return json.loads(p.read_text()) if p.exists() else None


def leaf(root, rel):
    r = J(Path(root) / rel / "report.json")
    return r if r and r.get("accepted") else None


def img(path):
    p = Path(path)
    if not p.exists():
        return '<p class="missing">figure not yet generated</p>'
    return f'<img src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}" alt="{html.escape(p.stem)}">'


def chip(v):
    cls = {"PASS": "pass", "FAIL": "fail"}.get(v.split()[0] if v else "", "partial")
    return f'<span class="chip {cls}">{html.escape(v)}</span>'


def fmt(x, d=4):
    return "—" if x is None else (f"{x:.{d}g}" if isinstance(x, float) else str(x))


def main(argv=None):
    p = argparse.ArgumentParser(); p.add_argument("--root", required=True); p.add_argument("--figures", required=True); p.add_argument("--out", required=True)
    a = p.parse_args(argv); R, F = Path(a.root), Path(a.figures)
    figs = J(F / "figures.json") or {}
    C2, C3, C4 = J(R / "C2_RESULT.json"), J(R / "C3_RESULT.json"), J(R / "C4_RESULT.json")
    C5 = {k: J(R / f"C5_RESULT_{k}.json") for k in ("SRa", "SRa_ffs0", "SRd_ffs0")}
    C67 = J(R / "C67_RESULT.json") or J(R / "C67_PREVIEW_SRa.json")
    C8 = J(R / "C8_RESULT.json"); C1 = J(R / "C1_RESULT.json")
    FACE = J(R / "FACE_RESULT.json") or J(R / "FACE_RESULT_partial.json") or {"sources": {}}
    srmom = leaf(R, "srmom_SRb_free_S0n0_min/support"); mreg = leaf(R, "uv_SRa_mreg1e3_tip/tip"); tip = leaf(R, "uv_SRa_tip/tip"); chi_tip = leaf(R, "gate_tip_linf_1e2_unit/tip")
    fig3v2 = leaf(R, "fig3_pure_linf_1e3_tip_v2/tip")
    n_fig3 = len([1 for d in glob.glob(str(R / "fig3_dir??")) if leaf(R, Path(d).name)])

    def verdict_c1():
        if C1: return C1["verdict"]
        return "PARTIAL" if fig3v2 else "not run"
    ledger = [
        ("C1", "Fig.3 纯幺正区域", verdict_c1(), f"+x 端 {fmt(fig3v2['verification']['f00_3'] if fig3v2 else None,6)} / 2.23289（−0.38%）；方向点 {n_fig3}/24"),
        ("C2", "Fig.4 手征塌缩到 f11=−f00/15", C2["verdict"] if C2 else "not run", "六个 ε 的 +x 端与 x_ref 截面宽在带内，端点随 ε 单调"),
        ("C3", "Fig.5 阈下分波与 S0 手征零点", C3["verdict"] if C3 else "not run", "三 ε RMS ≤ 6.4%（预算 8%），零点 0.426 / 0.293 / 无"),
        ("C4", "Fig.7 仅手征相移", C4["verdict"] if C4 else "not run", "S0/S2/P1 RMS 0.39/0.14/0.46°；P1 无 90° 穿越"),
        ("C5", "Fig.8 UV 收缩上大下小", " / ".join(f"{k}: {v['verdict'].split()[0]}" for k, v in C5.items() if v) or "not run", "三种读法上下收缩比 1.0 / 0.21 / 0.21（论文 6.3）；+x 端 −6.3% / −1.9%"),
        ("C6", "Fig.9 三代表点的 ρ", (C67 or {}).get("C6", {}).get("verdict", "pending"), "tip 穿越 773 MeV（η 0.585）；论文 813–827 MeV；面诊断：论文曲线不在最优面上"),
        ("C7", "Fig.10 S0/S2 低能重合", (C67 or {}).get("C7", {}).get("verdict", "pending"), "见 §Fig.10"),
        ("C8", "Fig.11 M/L 稳定性", C8["verdict"] if C8 else "pending", "见 §Fig.11"),
    ]
    def section(fid, title, quote, setting, rows, figure, verdict, note=""):
        trs = "".join(f"<tr><td>{html.escape(str(r[0]))}</td><td class='num'>{html.escape(str(r[1]))}</td><td class='num'>{html.escape(str(r[2]))}</td></tr>" for r in rows)
        return f"""
<section id="{fid}">
  <header><h2>{html.escape(title)}</h2>{chip(verdict)}</header>
  <div class="two">
    <div>
      <p class="quote">{quote}</p>
      <p class="setting"><b>我们的设置</b> {setting}</p>
      {('<p class="note">' + note + '</p>') if note else ''}
    </div>
    <table class="nums"><thead><tr><th>量</th><th>我们</th><th>论文</th></tr></thead><tbody>{trs}</tbody></table>
  </div>
  <figure>{figure}</figure>
</section>"""

    # ---------- per-figure rows
    rows3 = [("+x 端 f00(3)", fmt(fig3v2["verification"]["f00_3"] if fig3v2 else None, 7), "2.23289（数字化）"), ("Mreg（规则 (2)）", "10³", "未提"), ("方向点数", f"{n_fig3}/24", "—")]
    if C1: rows3 += [(r.get("quantity", ""), fmt(r.get("ours")), fmt(r.get("paper"))) for r in C1.get("rows", [])[:6]]
    rows2 = [(r.get("quantity", r.get("eps", "")), fmt(r.get("ours")), fmt(r.get("paper"))) for r in (C2 or {}).get("rows", [])[:8]] or [("+x 端 (ε=.002)", fmt(chi_tip["verification"]["f00_3"] if chi_tip else None, 6), "0.08257")]
    rows3b = [(f"{r.get('eps', r.get('wave', ''))} {r.get('wave', '')} RMS/f00(3)", fmt(r.get("rms_over_f00")), "≤ 8%") for r in (C3 or {}).get("rows", [])[:9]]
    rows4 = [(f"{r['wave']} RMS (deg)", fmt(r.get("rms_deg"), 3), "≤ 10") for r in (C4 or {}).get("rows", [])]
    rows5 = []
    for k, v in C5.items():
        if v:
            for r in v["rows"]:
                rows5.append((f"{k}: {r['quantity']}", fmt(r["ours"]), fmt(r["paper"])))
    rows9 = []
    for rel, lab in (("uv_SRa_tip/tip", "tip"), ("uv_SRa_section_hi", "x_ref 上支"), ("uv_SRa_section_xm001_hi", "x_ref−.001"), ("uv_SRa_section_xp001_hi", "x_ref+.001")):
        r = leaf(R, rel)
        if r:
            o = r["observables"]["P1"]; rows9.append((f"{lab}: P1 90° 穿越 (GeV)", fmt(o["crossing_90_GeV"], 4), "0.827 / 0.824 / 0.813")); rows9.append((f"{lab}: min η(P1) < 1.2 GeV", fmt(o["min_eta_below_1p2GeV"], 3), "未给出"))
    if C67 and C67.get("C6"): rows9 += [(str(r.get("point", "")), fmt(r.get("crossing_MeV")), "813–827 MeV") for r in C67["C6"].get("rows", [])]
    rows10 = [(f"{r.get('point','')}: S0 RMS / S2 RMS (deg)", f"{fmt(r.get('rms00_deg'),3)} / {fmt(r.get('rms20_deg'),3)}", "≤ 10") for r in ((C67 or {}).get("C7", {}).get("rows", []))]
    rows11 = [(f"M{r['M']} L{r['L']}: P1 穿越 (GeV) / +x 端", f"{fmt(r['P1_crossing_GeV'],4)} / {fmt(r['x_tip'],5)}", "0.82–0.87 / 0.0811") for r in figs.get("fig11", [])]
    if C8: rows11 += [(r.get("quantity", ""), fmt(r.get("ours")), fmt(r.get("paper"))) for r in C8.get("rows", [])]
    # ---------- diagnostics rows
    rowsD = []
    for src, e in FACE.get("sources", {}).items():
        name = Path(src).parent.parent.name if Path(src).parent.name == "tip" else Path(src).parent.name
        for fn, row in e["functionals"].items():
            rowsD.append((f"{name}: {fn} 面上范围", f"[{fmt(row['range_lo'])}, {fmt(row['range_hi'])}]", "论文相移要求 1−Re S > 1（任意 η）" if fn.startswith("ImKH_P1") else ("论文 Im S > 0" if fn.startswith("ImS") else "—")))
    if srmom:
        m = srmom["verification"]["functional"]["value"]; rowsD.append(("SR-1: 其余三盒 ±10% 下 S0 n=0 矩的最小值", f"{m:.4e} = 1.570 × 目标", "±10% 盒需 ≤ 2.62e−3 → 不可行"))
    if mreg and tip:
        rowsD.append(("Mreg 对照: +x 端 (10³ vs 10²)", f"{mreg['verification']['f00_3']:.6f} vs {tip['verification']['f00_3']:.6f}", "差 ≤ 2% 则 Mreg 非缺失项"))
        rowsD.append(("Mreg 对照: P1 穿越 (GeV)", f"{fmt(mreg['observables']['P1']['crossing_90_GeV'],4)} vs {fmt(tip['observables']['P1']['crossing_90_GeV'],4)}", "差 ≤ 20 MeV"))

    Q = {
        "fig3": "“The resulting shape is depicted in figure 3, where we have used M=50 in (3.61) for discretization and imposed unitarity for 10 partial waves per isospin. Inside the shape are all the possible values that f<sub>0</sub><sup>0</sup>(s=3) and f<sub>1</sub><sup>1</sup>(s=3) can have under the constraints of analyticity, crossing and unitarity.”",
        "fig4": "“We plot the allowed space … restricted by the chiral constraints (3.64) with tolerances ε<sup>χ</sup> = 6×10<sup>−3</sup>, 4×10<sup>−3</sup>, 2×10<sup>−3</sup>, 1×10<sup>−3</sup>, 6×10<sup>−4</sup>, 2×10<sup>−4</sup> (from the outer shape inward). The black line are the values given by the linear Weinberg model with varying values of f<sub>π</sub> and the black dot the one with f<sub>π</sub> ≃ 92 MeV.” … “with some norm and tolerance ε<sup>χ</sup>”",
        "fig5": "“For the larger tolerances the partial waves are not approximately linear in 0&lt;s&lt;4 and therefore we discard them. … The most notable deviation from linearity is in the S0 wave since the chiral zero of the amplitude disappears for the blue points. … The value ε<sup>χ</sup>=0.002 (green points), that we now choose, allows the physical value of f<sub>π</sub>.”",
        "fig7": "“We choose a point closest to the black dot to explore the partial waves in the physical region. Notice that only points at the boundary (green points) have partial waves associated with them. The partial waves at the magenta point and other nearby agree very well with experimental values for the S0 and S2 waves but not for the P1.”",
        "fig8": "“the plots are produced by taking ε<sup>SR</sup>=2×10<sup>−3</sup> in (3.73) and ε<sup>FF</sup>=6×10<sup>−5</sup>.” Caption: “The shape after adding the SVZ finite energy sum rules. The upper boundary shrinks notable but the lower not so much. We expect that the upper boundary has acquired QCD information.” (3.73): “||(π/M)Σ<sub>i</sub>(ds/dφ)<sub>i</sub> s<sub>i</sub><sup>n</sup>ρ<sub>i</sub> − QCD value|| ≤ ε<sup>SR</sup> … the tolerance is chosen such that enough information of the sum rule is put into the bootstrap yet not too strictly to make the problem infeasible.”",
        "fig9": "“The phase shifts for those points … agree reasonably well with experiment including the ρ meson resonance in the P1 channel. The resonance energy where the phase shift crosses π/2 is slightly shifted from the real world data on the mass of the rho at 770 MeV by roughly 6%.” … “Previous experience with the bootstrap would suggest looking at the tip of the shape (here the red dot). However the results seem quite robust so we choose two other points (pink, light pink) near the chiral point.”",
        "fig10": "“Notice that the three bootstrap (red, pink, light pink) curves in fig. 10 coincide at low energy but they spread at high energy. This indeed suggest that they are determined from the low energy side.”",
        "fig11": "“For M=50, we have plotted the results with L=8, 10, 12. In all three partial waves there is a reasonably good convergence. The situation with fixed L=10 and varying M=45, 50, 60 is slightly more complicated. … the S0, S2 phase shifts have a very good convergence, while the P1 phase shifts show a ρ resonance with the peak slightly shifted depending on M.”",
        "diag": "“While the low energy constraints significantly reduce the shape, there are still infinitely many possible scattering amplitudes contained within this thin shape, which differ from each other in terms of higher energy behaviors and this cannot be simply visualized in such a 2d projection. The sum rule and asymptotic behaviors of the form factor however helps navigate the bootstrap to the particular theories in these alternative dimensions of parameters.”",
    }
    S = {
        "fig3": "M=50、L=10、ν₀=0、cot 核 (3.67)、SVD 降基；无手征；双谱密度 |ρ_ij| ≤ Mreg=10³（规则 (2)：六个首遗漏波 η ≤ 1.02，L=8/10/12 端点差 ≤ 2%）。24 个支撑方向按 15° 步长。",
        "fig4": "chi-b：四点比值残差合并 8 维 L2 ≤ ε^χ（作者 2403 代码同一打包；chi-c 两球对照只报端点差）；Mreg=10²；六个 ε 各取 +x 端，ε=.002 另取 x_ref 及 ±0.001、±0.002 的上下截面。",
        "fig5": "阈下曲线由验收叶的 Arb 系数直接求值（0.05 ≤ s ≤ 3.95）；ε=.002/.004/.006 各用冻结规则的 x_ref 上支代表点。",
        "fig7": "冻结代表点规则（15:10Z，早于任何相移读取）：五个竖截面上支端点中离黑点最近者；仅手征约束，无 UV。",
        "fig8": "加 Gram (3.68)、FESR (3.72–3.73)、FF (3.75)；ε^SR=2×10⁻³ 取原始矩逐矩盒（SR-a，唯一可行的字面读法），对照逐波 L2（SR-d）与相对 10%（SR-b，不可行，见诊断）；ε^FF=6×10⁻⁵ 逐节点因子（主线）与冻结在 s₀（对照）。",
        "fig9": "与 Fig.8 同一 UV 模型；代表点 = tip 与冻结规则的截面上支点；相移由 S=1+iκh 的节点值按最近节点提升，η=|S| 一并给出（论文未画）。",
        "fig10": "同 Fig.9 的代表点。",
        "fig11": "固定 Mreg=10²（ℓ∞ 是节点值上界，不随 M 变）；(M,L)=(50,8)(50,10)(50,12)(45,10)(60,10) 的 UV tip。",
        "diag": "退化面诊断：在验收的支撑叶上加板 d·(f00,f11) ≥ v*−2×10⁻⁶，目标换成节点泛函 1−Re S、Im S 的 ±；论文 Fig.9 在 0.792/0.864 GeV 的相移 cos 2δ<0，故对任意 η 要求 1−Re S>1。SR-1：去掉 S0 n=0 的盒，极小化该矩。",
    }
    ledger_rows = "".join(f"<tr><td class='id'>{c}</td><td>{html.escape(t)}</td><td>{chip(v)}</td><td>{html.escape(n)}</td></tr>" for c, t, v, n in ledger)

    # ---- post-reply annotations: unitarity-saturation iteration and the eps_FF / eps_SR sensitivity (verdicts untouched)
    gw = gather_watson(a.root); PW = J(R / "C67_POSTWATSON.json")
    def g(x, d=3):
        return "—" if x is None else (f"{x:.{d}f}" if isinstance(x, float) else str(x))
    ZH = {"tip": "tip（+x 端）", "x_ref upper": "x_ref 上支", "x_ref+0.001 upper": "x_ref+0.001 上支"}
    ctr = ""
    for lab, rows in gw["chains"].items():
        if len(rows) < 2: continue
        r0, rl = rows[0], rows[-1]
        ctr += (f"<tr><td>{ZH.get(lab, lab)}</td><td class='num'>{rl['round']}</td><td class='num'>({r0['f00']:.4f}, {r0['f11']:.5f}) → ({rl['f00']:.4f}, {rl['f11']:.5f})</td>"
                f"<td class='num'>{g(r0['min_eta_P1'])} → {g(rl['min_eta_P1'])}</td><td class='num'>{g(r0['rho_MeV'],0)} → {g(rl['rho_MeV'],0)}</td>"
                f"<td class='num'>{g(r0['S0_at_1GeV'],0)}° → {g(rl['S0_at_1GeV'],0)}°</td><td class='num'>{g(r0['S2_at_1p2'],1)}° → {g(rl['S2_at_1p2'],1)}°</td></tr>")
    pin = "；".join(f"{ZH.get(k,k)} {g(v[0])} → {g(v[-1])}（{len(v)-1} 轮）" for k, v in gw["pinned"].items() if len(v) > 1)
    ftr = "".join(f"<tr><td>{r['eps_ff']:g}（两流同值）</td><td class='num'>{r['x_tip']:.5f}</td><td class='num'>{g(r['rho_MeV'],0) if r['rho_MeV'] else '1.2 GeV 以下无'}</td><td class='num'>{g(r['min_eta_P1'])}</td><td class='num'>{'/'.join(g(v,0) for v in r['S0_deg'])}</td></tr>" for r in gw["epsff"])
    ftr += "".join(f"<tr><td>{r['label'].replace('S0 cap','S0 上限').replace('P1 cap','P1 上限')}</td><td class='num'>{r['x_tip']:.5f}</td><td class='num'>{g(r['rho_MeV'],0) if r['rho_MeV'] else '1.2 GeV 以下无'}</td><td class='num'>{g(r['min_eta_P1'])}</td><td class='num'>{'/'.join(g(v,0) for v in r['S0_deg'])}</td></tr>" for r in gw["one_current"])
    ftr += "<tr><td>论文</td><td class='num'>0.0811</td><td class='num'>813–827</td><td class='num'>—</td><td class='num'>76/83/86/98（红）；99/103/104/109（浅粉）</td></tr>"
    str_ = "".join(f"<tr><td>{'去掉 S0 n=0 盒' if isinstance(r['eps_sr'], str) else format(r['eps_sr'],'g')}</td><td class='num'>{r['x_tip']:.5f}</td><td class='num'>{g(r['rho_MeV'],0) if r['rho_MeV'] else '—'}</td><td class='num'>{g(r['min_eta_P1'])}</td><td class='num'>{g(r['S0_at_0p95'],0)}°</td></tr>" for r in gw["epssr"])
    mp = matched_pair_summary(a.root); mpl = ""
    if mp:
        ok = lambda b: "满足" if b else "不满足"
        c = mp["criteria"]
        mpl = (f"<p class='note'><b>配对试验结果</b>（S0 上限 2×10⁻⁴、P1 上限 8×10⁻⁵，其余同主线，未迭代）：+x 端 {mp['x_tip']:.5f}（与 0.0811 差 {mp['x_dev_pct']:+.1f}%，{ok(c['x_tip_within_1pct'])}）；"
               f"ρ {mp['rho_MeV']:.0f} MeV（{ok(c['rho_in_band'])}）；S0 到 1 GeV 与论文红线最大差 {mp['S0_max_1GeV']:.1f}°、RMS {mp['S0_rms_1GeV']:.1f}°（{ok(c['S0_within_10deg_to_1GeV'])}）；"
               f"S2 最大差 {mp['S2_max']:.1f}°；ρ 处 min|S<sub>P1</sub>| {mp['min_eta_P1']:.2f}。预登记判据（三者同时）<b>{'满足' if mp['met'] else '不满足'}</b>——"
               f"这一对上限在 tip 处给出论文 Fig.9、Fig.10 的红线，但区域 +x 端仍差 3%；+x 端随两个上限都升、ρ 随 P1 上限升，故按我们的读法没有一对上限能同时给出三者，剩余差别在别处（ε^SR 的范数、(3.75) 在 s₀ 以上的归一化，或作者的迭代对区域本身的作用）。按预登记，x_ref 截面的同参数试验不再运行。</p>")
    pw = ""
    if PW:
        c6 = PW.get("C6", {}); c7 = PW.get("C7", {})
        pw = ("<p class='note'><b>迭代后的 C6/C7 复评（标签 post-Watson，不改主线裁决）</b>：C6 " + html.escape(str(c6.get("verdict"))) + "，ρ 穿越 "
              + "、".join(f"{r['point']} {g(r['crossing_MeV'],0)} MeV（min|S| {g(r['min_eta'])}）" for r in c6.get("rows", [])) + "；C7 " + html.escape(str(c7.get("verdict")))
              + "，" + html.escape(str(c7.get("evidence", ""))) + "。</p>")
    wsec = f"""
<section id="watson"><header><h2>作者回复后的补充：幺正饱和迭代与 ε^FF 敏感性</h2>{chip("注释（裁决不变）")}</header>
  <p class="quote">“|S| does not saturate unitarity. We have an iterative procedure that improves saturation of unitarity … Our results are always shown after that.” “Unitarity tends to be unsaturated near the resonance … we used the iterations to correct that. … The faster (or sometimes slower) rise of S0 happens, I believe, depending on the parameters. Also changes in the rho mass.”（作者 2026-09-17 来信）</p>
  <p class="setting"><b>我们的做法</b> 按 follow-up 论文 (arXiv:2403.10772) 式 (2.29) 实现饱和迭代：保留有限问题的全部约束，放开截面 f₀⁰(3)=x，目标改为 Σ<sub>k</sub> Re[e<sup>−2iα<sub>k</sub></sup>(S<sub>k</sub>−1)]（k 取 s₀ 以下 S0、P1、S2 的全部节点，α<sub>k</sub> 为上一轮形状因子相位，无流的分波取上一轮相移）；每轮都是一次完整的 192 bit SDPB 求解 + Arb 复验，目标行用独立浮点构造核对到 1e−11。从 Fig.9 的三个代表点各跑若干轮；另以“钉住边界点”（保留截面 + 2e−6 近优板）作对照。作者代码（2403 / 2505 版）读过未运行：变量按 h̃=h/Λ<sub>ℓ</sub> 缩放，相当于把我们的节点权 1 换成 1/Λ<sub>ℓ</sub>²，两种权都跑了。</p>
  <div class="cmp"><table><thead><tr><th>起点</th><th>轮数</th><th>(f₀⁰, f₁¹) 起 → 终</th><th>min|S<sub>P1</sub>|</th><th>ρ (MeV)</th><th>S0 @ 1 GeV</th><th>S2 @ 1.2 GeV</th></tr></thead><tbody>{ctr}</tbody></table></div>
  <p class="note">钉住边界点的对照（min|S<sub>P1</sub>|）：{pin}。结论：迭代确实把 s₀ 以下每个节点的 |S| 推到 0.90–1.00 并使 F 与 S 的相位对齐，但只有离开边界点才做得到（f₀⁰ 每轮内缩 1–2%，未停）；钉住时几乎不饱和。ρ 位置不随迭代改变（论文 813–827 MeV），S0 仍快，S2 偏离论文曲线，三条链也不收敛到同一振幅；作者权重给出同样的 ρ。</p>
  {img(F / "fig9_watson.png")}
  {pw}
  <p class="setting"><b>ε^FF 与 ε^SR 敏感性</b> UV 阶段仅有的两个连续参数。ε^SR 从 2×10⁻³ 放到 10⁻²、或整个去掉 S0 n=0 盒，ρ 至多移动 12 MeV，S0 不变（放开的 S0 n=0 矩落在 QCD 值的 1.9 倍附近）。ε^FF 则很敏感，且两条流互不干扰：只放 S0 上限到 2×10⁻⁴ 把 S0 放到论文红线上（0.95 GeV 以下 3° 内）而 ρ 不动；只放 P1 上限则 ρ 上移（10⁻⁴ 约 850、2×10⁻⁴ 约 970 MeV）而 S0 不动。内插给出论文有效上限约 S0 2×10⁻⁴、P1 8×10⁻⁵；配对试验按预登记判据（x_tip 与 0.0811 差 1% 内、ρ 在 813–827、S0 到 1 GeV 差 10° 内，三者同时）评判。</p>
  <div class="cmp"><table><thead><tr><th>ε^FF</th><th>+x 端</th><th>ρ (MeV)</th><th>min|S<sub>P1</sub>|</th><th>δ<sub>S0</sub> @ 0.79/0.86/0.95/1.06 GeV</th></tr></thead><tbody>{ftr}</tbody></table></div>
  <div class="cmp"><table><thead><tr><th>ε^SR（SR-a 盒）</th><th>+x 端</th><th>ρ (MeV)</th><th>min|S<sub>P1</sub>|</th><th>δ<sub>S0</sub> @ 0.95 GeV</th></tr></thead><tbody>{str_}</tbody></table></div>
  {mpl}
  {img(F / "fig_epsff.png")}
</section>
"""
    body = f"""
<title>2309.12402 复现台账</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,wght@0,600;1,400&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono:wght@400&display=swap');
:root{{--bg:#F6F7F9;--ink:#1B2430;--muted:#5B6673;--quote:#6E7783;--line:#D9DEE5;--accent:#1F5F8B;--pass:#2E7D4F;--partial:#B7791F;--fail:#B23B3B;--card:#FFFFFF}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{--bg:#12161B;--ink:#E6EAEF;--muted:#A3ADB8;--quote:#9AA4AF;--line:#2A323B;--accent:#7FB3D5;--pass:#6FCB94;--partial:#E0B25C;--fail:#E27D7D;--card:#1A2027}}}}
:root[data-theme="dark"]{{--bg:#12161B;--ink:#E6EAEF;--muted:#A3ADB8;--quote:#9AA4AF;--line:#2A323B;--accent:#7FB3D5;--pass:#6FCB94;--partial:#E0B25C;--fail:#E27D7D;--card:#1A2027}}
body{{background:var(--bg);color:var(--ink);font-family:'IBM Plex Sans',system-ui,sans-serif;font-size:15px;line-height:1.55;padding-inline:clamp(16px,4vw,48px);padding-block:32px 64px;max-width:1180px;margin:0 auto}}
h1{{font-family:'Source Serif 4',Georgia,serif;font-weight:600;font-size:clamp(24px,3vw,34px);text-wrap:balance;margin:0 0 4px}}
h2{{font-family:'Source Serif 4',Georgia,serif;font-weight:600;font-size:20px;margin:0}}
.sub{{color:var(--muted);margin:0 0 24px;max-width:70ch}}
.lede{{max-width:72ch}}
.chip{{display:inline-block;font-family:'IBM Plex Mono',monospace;font-size:12px;letter-spacing:.04em;padding:2px 9px;border-radius:3px;color:#fff;background:var(--partial)}}
.chip.pass{{background:var(--pass)}} .chip.fail{{background:var(--fail)}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}} th,td{{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}} th{{font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted)}}
td.id{{font-family:'IBM Plex Mono',monospace;color:var(--accent)}} td.num{{font-family:'IBM Plex Mono',monospace;font-size:13px}}
.ledger{{overflow-x:auto;margin:0 0 40px}}
section{{border-top:1px solid var(--line);padding-block:28px 12px}}
section header{{display:flex;gap:14px;align-items:center;flex-wrap:wrap;margin-bottom:12px}}
.two{{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr);gap:24px}} @media (max-width:760px){{.two{{grid-template-columns:1fr}}}}
.quote{{font-family:'Source Serif 4',Georgia,serif;font-style:italic;color:var(--quote);margin:0 0 12px;max-width:65ch}}
.setting{{margin:0 0 8px;max-width:65ch}} .note{{color:var(--muted);font-size:14px;max-width:65ch}}
figure{{margin:18px 0 0;background:var(--card);border:1px solid var(--line);padding:8px}} figure img{{max-width:100%;display:block}}
.missing{{color:var(--fail);font-family:'IBM Plex Mono',monospace}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:12px}} @media (max-width:760px){{.grid2{{grid-template-columns:1fr}}}}
.cmp{{overflow-x:auto}}
</style>
<h1>He–Kruczenski 2309.12402 复现台账</h1>
<p class="sub">SDPB 路线：论文输入 → Mathematica 独立核对 → PMP → SDPB 192 bit → Arb 全约束复验 → 预登记规则裁决。斜体灰字为论文原文；每个数字都对应结果根目录里一个已验收的 report.json。</p>
<p class="lede"><b>两句话结论。</b>(1) ρ 峰从底层物理复现出来了：仅手征约束时 P1 在 1.2 GeV 以下无 90° 穿越，加入 QCD 求和规则与形状因子渐近后 tip 处 P1 在 773 MeV 穿越（论文 813–827，物理 770），峰上 η=0.585。(2) 端到端流水线存在且每步有收据；低能侧 C1–C4 通过，Fig.8–10 的精细特征不是论文正文所述有限问题最优集的性质，已在求解器精度下证明。</p>
<div class="ledger"><table><thead><tr><th>claim</th><th>内容</th><th>裁决</th><th>关键数字（我们 / 论文）</th></tr></thead><tbody>{ledger_rows}</tbody></table></div>
{section("fig3","Fig.3 纯幺正区域",Q["fig3"],S["fig3"],rows3,img(F/"fig3.png"),verdict_c1())}
{section("fig4","Fig.4 手征约束区域",Q["fig4"],S["fig4"],rows2,img(F/"fig4.png"),C2["verdict"] if C2 else "not run")}
{section("fig5","Fig.5 阈下分波",Q["fig5"],S["fig5"],rows3b,img(F/"fig5.png"),C3["verdict"] if C3 else "not run")}
{section("fig7","Fig.7 仅手征相移",Q["fig7"],S["fig7"],rows4,img(F/"fig7.png"),C4["verdict"] if C4 else "not run")}
{section("fig8","Fig.8 加入 UV 后的区域",Q["fig8"],S["fig8"],rows5,img(F/"fig8.png"),ledger[4][2],"三种可行读法都没有出现“上收缩 ≫ 下抬升”；打包方式对截面的影响 &lt; 1e−6；相对 10% 读法不可行（SR-1 证明）。")}
{section("fig9","Fig.9 P1 相移与 ρ",Q["fig9"],S["fig9"],rows9,img(F/"fig9.png")+img(F/"fig9_eta.png"),ledger[5][2],"x_ref 上支点 |S_P1| 在 0.73 GeV 跌到 0.24，节点值允许两条相位分支（最近节点提升，或 +180° 后在 0.708 GeV 穿越 90°）；有限问题不决定分支。")}
{section("fig10","Fig.10 S0 与 S2 相移",Q["fig10"],S["fig10"],rows10,img(F/"fig10.png"),ledger[6][2])}
{section("fig11","Fig.11 M 与 L 依赖",Q["fig11"],S["fig11"],rows11,img(F/"fig11.png"),ledger[7][2])}
{section("diag","为什么复现不出 Fig.8–10 的精细特征：退化面诊断与矩区间",Q["diag"],S["diag"],rowsD,img(F/"face_ranges.png"),"证明","每个泛函值由 Arb 从节点 S 独立复算；板宽 2e−6 是对偶间隙阈值的两倍；两端极值都顶在板边。")}
<section id="settings"><header><h2>设置逐项对照：2309 正文 · 我们 · 作者 2403 代码 · 作者当前代码 (2505)</h2></header>
<div class="cmp"><table><thead><tr><th>项</th><th>2309 正文</th><th>我们（主线 / 对照）</th><th>2403 代码</th><th>2505 代码</th></tr></thead><tbody>
<tr><td>s₀ / ν₀</td><td>1.2 GeV / 0</td><td>同</td><td>2 GeV / −20</td><td>2 GeV / −20</td></tr>
<tr><td>流与矩</td><td>S0(0,1)、P1(−1,0)</td><td>同</td><td>S0/P1/D0 各 3</td><td>S0/P1/D0 各 6</td></tr>
<tr><td>ε^SR</td><td>“某范数” 2×10⁻³</td><td>原始矩逐矩盒 / 逐波 L2 / 相对 10%（不可行）</td><td>逐流 L2 于归一化矩 1e−7/6e−6/5e−6</td><td>相对盒 5% / 10% / 10%</td></tr>
<tr><td>ε^FF</td><td>6×10⁻⁵，(3.75)</td><td>逐节点因子 / 冻结 s₀</td><td>松弛因子 0.05、2×6.87</td><td>乘子 8 / 2.2 / 8</td></tr>
<tr><td>手征</td><td>4 点比值，“某范数” 2×10⁻³</td><td>8 维 L2 / 两球</td><td>L2 2e−3</td><td>L2 2e−3</td></tr>
<tr><td>密度正则化</td><td>未提</td><td>ℓ∞ Mreg=10²（物理判据；ℓ2/ℓ4 散布 0.2%）</td><td>‖ρ/Mρ‖₄ ≤ 10²</td><td>同 2403</td></tr>
<tr><td>代表振幅</td><td>“边界点自带分波”</td><td>支撑点的解 + 退化面诊断</td><td>定 F0 极大 F1，再在整个可行集上极大化 Watson 泛函</td><td>纯可行性起步 + 5 轮 Watson；无支撑泛函</td></tr>
<tr><td>求解器</td><td>未说明</td><td>SDPB 192 bit，gap 1e−6，Arb 复验</td><td>MOSEK</td><td>MOSEK</td></tr>
</tbody></table></div></section>
{wsec}
<section id="limits"><header><h2>未完成项与记录在案的事故</h2></header>
<p class="note">事故：Fig.3 方向扫描因源哈希校验失败重解（端点逐位相同）；SR-b 两驱动在复验阶段因运行中改源码而崩溃，SDPB 输出完整，用 finish_leaf.py 只重跑读回/复验；/tmp 卷被占满，大文件改到 /playpen1。全部记录在 GATE_LOG.md。</p>
</section>
"""
    Path(a.out).write_text(body)
    print("written", a.out, Path(a.out).stat().st_size, "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
