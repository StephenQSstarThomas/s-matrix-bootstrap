# 2309.12402 复现台账（供审核）

分支 `sdpb-2309-regularised`，定稿 2026-09-14 15:45Z。本目录是审核用的自足参考包：

| 路径 | 内容 |
|---|---|
| `index.html` | 最终台账页（自足 HTML，图已内嵌）：两句话结论、八条 claim 台账、每张论文图的原文（斜体灰字）/我们的设置/数字/原图 vs 我们并排图、退化面与 SR-1 证明、2309 正文 · 我们 · 2403 代码 · 2505 代码的设置对照表。在线副本：https://claude.ai/code/artifact/db9c21c9-6e7e-44c8-9d59-370cbdca79e4 |
| `figures/fig3.png … fig11.png, fig9_eta.png, face_ranges.png` | 并排对照图（左/上为论文 PDF 渲染，右/下为我们的验收叶 + 数字化论文曲线）；`figures.json` 为图中数字 |
| `receipts/C1_RESULT.json … C8_RESULT.json` | 八条 claim 的预登记规则裁决（C5 三种读法、C6/C7 合于 `C67_RESULT.json`） |
| `receipts/FACE_RESULT.json` | 退化面诊断 14 个叶子的范围与 F1–F6 裁决 |
| `receipts/UV_REPRESENTATIVE_SELECTION.json`, `IR_REPRESENTATIVE_SELECTION.json` | 冻结代表点规则的选点收据（先冻结再读相移） |
| `receipts/MREG_STAR_M50.json`, `MREG_STAR_PURE.json` | 正则化尺度的物理判据选择 |
| `receipts/GATE_PREREGISTRATION.md`, `GATE_LOG.md`, `GATE_RESULT.md` | 预登记规则与带时间戳的全过程日志（含事故与双会话说明） |
| `receipts/queue.log` | 全部求解任务的启动记录 |

每个 `report.json` 叶子（PMP、SDPB 输出、Arb 复验）在结果根目录 `/playpen1/shiqiu/sdpb_regularised_20260913`（仓库内符号链接 `results/runs/sdpb_regularised_20260913`），未纳入 git（数 GB）。

生成脚本：`scripts/sdp/final_figures.py`（图）、`scripts/sdp/final_report_html.py`（台账页）、`scripts/sdp/c1_eval.py … c8_eval.py`、`scripts/sdp/face_eval.py`（裁决）。叙述性报告：仓库根 `REPORT_SDPB_2309_ZH.md`；计划与预登记：`PLAN_SDPB_2309_ZH.md`；代码分级：`TRIAGE_ZH.md`。

未完成项（不改任何裁决）：C1 的 Mreg=10⁴ 敏感性注释（纯幺正 tip 重解与 −x/+y/−y 三个支撑）在结果根目录中运行，完成后补入 `index.html` 与 `GATE_LOG.md`。
