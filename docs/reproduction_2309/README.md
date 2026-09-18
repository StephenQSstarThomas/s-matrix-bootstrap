# 2309.12402 复现台账（供审核）

分支 `sdpb-2309-regularised`，定稿 2026-09-14 15:45Z；2026-09-18 按作者回复增补幺正饱和迭代与 ε 敏感性（注释，裁决不变）。本目录是审核用的自足参考包：

| 路径 | 内容 |
|---|---|
| `REPRODUCTION_COMPARISON_EN.pdf` | **给作者的英文对照文件（LaTeX 编译）**：方法与工具链、五个论文未定项及读法、八张图的逐图对照（论文原文灰色斜体 + 我们的数字 + 并排图）、退化面与 SR-1 诊断、§5 幺正饱和迭代（作者来信灰色斜体 + 三条链的表与图）、§6 ε^FF 敏感性、向作者请求代码。同内容的 `.tex`（可重编译）、`.md`、自足 `.html` 并列 |
| `gtb2309_sdpb_reproduction_code.zip` | **给作者的代码包**：`src/`（算子推导与 PMP 写出）、`scripts/`（求解驱动、评估、图与文档）、`tests/`、数字化参考曲线、英文推导笔记、启动脚本示例；含英文 `OVERVIEW.md`（流水线与运行方法）与 `DERIVATIONS.md`（论文方程 → 代码映射）。不含任何运行输出 |
| `index.html` | 最终台账页（自足 HTML，图已内嵌）：两句话结论、八条 claim 台账、每张论文图的原文（斜体灰字）/我们的设置/数字/原图 vs 我们并排图、退化面与 SR-1 证明、2309 正文 · 我们 · 2403 代码 · 2505 代码的设置对照表、作者回复后的饱和迭代与 ε 敏感性一节。在线副本：https://claude.ai/code/artifact/db9c21c9-6e7e-44c8-9d59-370cbdca79e4 |
| `figures/fig3.png … fig11.png, fig9_eta.png, face_ranges.png, fig9_watson.png, fig_epsff.png` | 并排对照图（左/上为论文 PDF 渲染，右/下为我们的验收叶 + 数字化论文曲线）；`figures.json` 为图中数字 |
| `receipts/C1_RESULT.json … C8_RESULT.json` | 八条 claim 的预登记规则裁决（C5 三种读法、C6/C7 合于 `C67_RESULT.json`） |
| `receipts/FACE_RESULT.json` | 退化面诊断 14 个叶子的范围与 F1–F6 裁决 |
| `receipts/UV_REPRESENTATIVE_SELECTION.json`, `IR_REPRESENTATIVE_SELECTION.json` | 冻结代表点规则的选点收据（先冻结再读相移） |
| `receipts/MREG_STAR_M50.json`, `MREG_STAR_PURE.json` | 正则化尺度的物理判据选择 |
| `receipts/GATE_PREREGISTRATION.md`, `GATE_LOG.md`, `GATE_RESULT.md` | 预登记规则与带时间戳的全过程日志（含事故与双会话说明） |
| `receipts/queue.log` | 全部求解任务的启动记录 |

每个 `report.json` 叶子（PMP、SDPB 输出、Arb 复验）在结果根目录 `/playpen1/shiqiu/sdpb_regularised_20260913`（仓库内符号链接 `results/runs/sdpb_regularised_20260913`），未纳入 git（数 GB）。

生成脚本：`scripts/sdp/final_figures.py`（图）、`scripts/sdp/final_report_html.py`（台账页）、`scripts/sdp/c1_eval.py … c8_eval.py`、`scripts/sdp/face_eval.py`（裁决）。叙述性报告：仓库根 `REPORT_SDPB_2309_ZH.md`；计划与预登记：`PLAN_SDPB_2309_ZH.md`；代码分级：`TRIAGE_ZH.md`。

已补齐：C1 的 Mreg=10⁴ 敏感性注释（GATE_LOG 09-17 18:05Z）。迭代与敏感性的全部叶子在结果根目录 `watson_*`、`sens_*`；`REPORT_SDPB_2309_ZH.md` §8 为中文叙述。s₀ 敏感性等剩余事项列在仓库根 `TODO_ZH.md`。
