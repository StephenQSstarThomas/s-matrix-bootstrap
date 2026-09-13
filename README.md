# He–Kruczenski 2309.12402 复现（SDPB 路线）

目标：用 SDPB 严谨复现 [Bootstrapping gauge theories, arXiv:2309.12402v3](references/2309.12402v3.pdf) 的 Fig.3–11，即八条 claim C1–C8。分支 `sdpb-2309-regularised`。

**入口文档（只有这两份是现行的）**

- [`PLAN_SDPB_2309_ZH.md`](PLAN_SDPB_2309_ZH.md)：主线合同、预登记裁决规则、精度预算、七个阶段的流水线与验收、门槛结果。
- [`TRIAGE_ZH.md`](TRIAGE_ZH.md)：代码、结果目录、参考资料、旧文档的四级分类，以及五条已被证伪命题的正确版本。

其余顶层文档（STATUS、SCIENCE、HANDOFF、REPORT_SDP、REVISED_CLAIMS、REPRODUCTION_GUIDE、AGENTS）是历史记录，顶部横幅说明其失效之处；与 PLAN 冲突时以 PLAN 为准。

**主线合同一句话**：论文的有限问题（M=50、L=10、mixed-pv 节点离散、chi-b 合并 8 维手征范数、ε^χ=0.002、两个流、四个 FESR 矩、(3.75) 形状因子界）加上作者方法论文 2103.11484 §3 要求的双谱密度正则化 |ρ_ij| ≤ Mreg，Mreg 由物理判据（未施加分波幺正、L=8/10/12 稳定）定为 10²；SDPB 192 bit，gap 10⁻⁶，全部原约束用 Arb 复验。

**当前状态（2026-09-13）**：门槛已通过。ε=0.002 的 +x 端 0.08276 对论文 0.08257（+0.22%），L=8/10/12 端点跨 0.18%，M=30/50 差 0.1%，ℓ∞ 与 ℓ2 正则化差 0.20%。正在进行阶段 2（Fig.3/4 区域）与阶段 0b 的 ℓ4 对照。结果根目录 `results/runs/sdpb_regularised_20260913`（符号链接到 `/playpen1`），其中 `GATE_RESULT.md`、`GATE_LOG.md` 是门槛的判定与全过程记录。

**运行**

```bash
# 安装（可编辑）
python -m pip install -e '.[test]'
# 一次求解：M50/L10、手征 ε=0.002、正则化 Mreg=1e2、+x 端
PYTHONPATH=src python -m smatrix_bootstrap.run sdp solve --workdir RUNS/tip \
  --M 50 --L 10 --chiral --chi chi-b --eps-chi 0.002 \
  --scattering-prescription mixed-pv --reduce-basis --reg-norm linf --reg-bound 1e2 \
  --operator-dps 40 --digits 30 --precision 192 --nproc 8 --duality-gap 1e-6 \
  --points tip --mma-reference <已通过的 M=50 Mathematica 核对 report.json>
# 复用同一 PMP 只换目标（截面、方向扫描）
PYTHONPATH=src python -m smatrix_bootstrap.run sdp support --source-report RUNS/tip/tip/report.json \
  --direction 0 1 --fix-f00 0.07332139057293464 --out RUNS/section_hi
# 未施加分波的幺正性诊断
PYTHONPATH=src python scripts/sdp/omitted_waves.py RUNS/tip/tip/report.json
# 测试
python -m pytest -q tests/sdp        # SDPB 路线，约 240 项
python -m pytest -q legacy/tests     # 退役 Newton 包自己的 107 项
```

依赖：SDPB 3.1.0（Docker 镜像 `bootstrapcollaboration/sdpb:3.1.0`，启动器 `scripts/sdpb/sdpb.sh`）、Wolfram Engine 15（仅独立公式核对，`scripts/mma/`）、python-flint、mpmath。不使用 MOSEK/CVX。

**代码布局**：`src/smatrix_bootstrap/sdp/` 是全部现行代码；`src/smatrix_bootstrap/run.py` 只分发 `sdp`。退役的 Newton 包在 `legacy/smatrix_bootstrap_newton/`（见 `legacy/README.md`），只作历史结果的复验与 `sdp/crosscheck.py` 的独立对照。`references/` 存论文源码、数字化图数据、作者 2403/2505 代码快照（提交 801684d）与 Córdoba 2511.11513 快照。
