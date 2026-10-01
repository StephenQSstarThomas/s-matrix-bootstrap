# 2309 验证层审查与 2403 的分级验证定义（2026-10-01）

本文件回答两个问题：2309 复现的"验证"到底证明了什么；2403 复现用什么标准称一个结论"已证明"。审查只读 2309 分支的代码与结果根 `/playpen1/shiqiu/sdpb_regularised_20260913`，不改任何 2309 裁决。

## 1. 2309 的验证实际证明了什么

| 层 | 2309 的实现 | 证据 | 结论 |
|---|---|---|---|
| 公式 | `scripts/mma/audit_2309.wls`：在 Mathematica 中独立重算公式（核、运动学因子、FESR 目标） | 每个叶子的 `mma/` 子目录与 "Mathematica 2309 audit passed" | 公式一致性检查，不是对解的证书 |
| 原问题可行性（节点处） | `sdp/arbaudit.py`：在 Arb 中重建全部节点约束；Legendre Q 经 Arb 的 ₂F₁（带误差半径）；对每个松弛给出 certified_pass / certified_fail / inconclusive | 模块文档字符串："Each UV slack has a certified pass/fail/inconclusive verdict … These are finite-node checks, not continuum or optimization certificates."；`report.json` 的 `verification.constraint_checks` 全为 true，`unitarity.max_absolute_violation_all ≈ −1e−88` | **严格**：被报告的振幅确实满足有限节点约束。这使"区域至少包含这个点"成为证明 |
| 最优性 | 无。`sdp/dual.py` 只做弱对偶重放，要求等式残差恰为零，从未成立 | 所有叶子 `certified: false`、`convergence.independent_dual_certified: false` | 区域边界、支撑值、极大极小都只是 SDPB 声称的最优值（对偶间隙 1e−6），没有严格上界 |
| 降基 | SVD 把变量限制在约束行张成的子空间（`reduce_basis: True`，1780–2185 个变量，全空间 3877） | `full_M50_L10_linf_1e2_tip` 是唯一的全空间叶子：0.0827775 对降基 0.0827577（+0.02%） | 降基是内逼近：只能使极大偏小、极小偏大。作为"可行点存在"成立；作为上界或排除不成立 |
| 排除类结论 | SR‑1（相对 10% 读法不可行）、退化面 F1′/F3′/F6 | `srmom_SRb_free_S0n0_min` 与全部 `face_*` 叶子：`reduce_basis: True`，`n_vars 1980`，`certified: false` | 两个缺口叠加（降基 + 无对偶证书），因此是强数值证据，**不是证明**。2309 报告中"求解器精度下证明"的措辞过强 |
| 相移读取 | 节点值，最近节点提升，90° 穿越在相邻节点间线性插值 | `sdp/observables.py` `phase_from_S`、`crossing_energy` | 读取约定；节点间距在 ρ 附近约 60 MeV，读数精度受此限制 |

注：2309 的数值证据本身没有被推翻。需要修正的是措辞：凡是"上界""排除""不可行"的结论，在补做全空间求解与对偶证书之前，只能称为数值证据。这一修正对 2309 台账的影响另行提交（需用户同意）。

## 2. 2403 的分级验证（2309 回溯共用）

- **V0 公式**：三方一致（我们的 Mathematica 核、Python/Arb 独立重建、作者 notebook 无头运行的参照导出），见计划 Task A6。
- **V1 原问题可行性（严格）**：在 Arb 中用独立重建的核复算被报告振幅的全部约束：幺正 3×3、复 3×3 Gram 的全部 7 个主子式、手征两个范数、三个 FESR 球、33 个 FF 模长、B 样条系数非负、ℓ4 正则化（精确的 ℓ4 范数）。结论只针对有限节点问题。
- **V2 对偶界（严格上界）**：对 max b·y s.t. C_j + Σ_n y_n A_jn ⪰ 0，取 SDPB 的 X_j，在 Arb 中把等式残差修复为零（在 A_jn 张成的方向上解 Gram 方程组，浮点 LU 起步、Arb 迭代精化），用 Arb 的 LDLᵀ 给出 λ_min(X̃_j) 的下界，然后
  max ≤ Σ_j ⟨C_j, X̃_j⟩ + Σ_j max(0, −λ_min(X̃_j)) · T̄_j，
  其中 T̄_j 是 Tr M_j(y) 在可行集上的先验上界，由推导文件 §12 的显式界给出。证书必须针对**全空间**问题：降基问题的对偶界不约束全问题。
- **V3 连续性抽查**：节点之间的幺正与 Gram，用推导文件规定的同一插值；只报告，不称证书。

台账中每个数字标出所达到的级别：可行点 = V1；上界、排除、不可行 = V2（或全空间求解 + V2）；V0 是前提。

## 3. 回溯 2309 的清单（Phase B）

全空间求解 + V2：手征 tip（`gate_tip_linf_1e2_unit`）、UV tip（`uv_SRa_tip`）、SR‑1（`srmom_SRb_free_S0n0_min`）、退化面的排除叶子（`face_SRa_tip_ImKH_P1_38_max`、`face_SRa_tip_ImKH_P1_39_max`、`face_SRa_xrefhi_ImKH_P1_38_min`、`face_SRa_tip_ImS_P1_38_max`、`face_SRa_ffs0_tip_ImKH_P1_38_max`、`face_SRa_ffs0_tip_ImKH_P1_39_max`）。逐条重新标注"已证明"或"数值证据"。
