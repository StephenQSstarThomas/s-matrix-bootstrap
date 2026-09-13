# 仓库整理分类（TRIAGE，2026-09-13）

分支 `sdpb-2309-regularised`。本文件回答"哪些是必要的推导、底座、核心函数；哪些有问题或不必存在；哪些有误导性"。分类标准：

- **CORE**：新主线直接依赖或复用的代码、数据、推导。
- **EVIDENCE**：保留作历史证据与出处，不进主线，不删。
- **MISLEADING**：包含已被证伪的结论或非论文的选择，保留但必须带说明；本文件即说明。
- **DISPOSABLE**：冒烟、中断、重复或可再生的大体积产物；阶段 1 只加索引不删，删除另行决定。

物理与工具链事实见 `PLAN_SDPB_2309_ZH.md`。本文件依据三份只读审计（SDPB/PMP 模块、结果目录与文档、旧包与脚本）以及作者代码对照。

## 1. 已证伪的五条命题（各文档中出现处见 §5）

| 编号 | 被证伪的说法 | 正确版本与依据 |
|---|---|---|
| a | "2309 没有任何正则化，密度界 B 是从 2403 搬来的，应删除" | 2309 的数值方法是作者自己的 2103.11484，该文 §3 要求对双谱密度范数加界并"always consider regularized primal problems"；2403/2505 代码与 Córdoba 2511 都有此界。删除它使有限问题病态 |
| b | "去掉 B 后 M=50 首次认证、区域/相移失败是截断或精度问题" | 无 B 的全部 SDPB 解 ‖ρ‖₄≈10¹⁶–10¹⁷，手征 +x 端 0.2–0.8 对论文 0.0826；同 ε 下 L10/L12 端点差 3.9 倍。失败是缺正则化的直接后果 |
| c | "M=50 在所有配置下不收敛；双精度求解器给的只是下界" | 本仓库自己的 Clarabel 阶梯（mixed-PV、chi-b、B 在场）M=20–40 给 0.0842/0.0825/0.0826/0.0825/0.0825，与论文 0.08257 一致。注意：该阶梯用的是 Clarabel（cvxpy），不是 MOSEK；此前口头称"MOSEK 阶梯"有误 |
| d | "mixed-PV 不是单一解析函数，需换成 sine-cardinal" | mixed-PV 就是作者的离散化：cot 核与作者 Mathematica 生成的矩阵逐元相同，中点权重、Legendre-Q、ρ₂ 打包一致。sine-cardinal 是 Córdoba 式有限正弦基，非作者做法；实测两者在同合同下端点 0.4975 对 0.5042，无关根因 |
| e | "chi-c 两个分离球、SR-a 逐矩绝对盒是 2309 字面合同；合并 L2 是输出辅助" | 作者 2403/2505 与 Córdoba 都用一个合并 8 维 L2；FESR 用归一化矩逐波 L2；2309 自身 Fig.5 残差指纹（1.000007/1.00002/1.00015）也指向合并 L2 |

派生的误导项：ε^χ=0.00021 的"fπ 校准"实验（改论文输入 10 倍去补偿区域爆炸）；L13 截断诊断；"严格排除"证书（Newton 集合内 10⁻⁶ 量级薄层的命题）；障碍中心代表点（内点，论文明说只有边界点有分波）；320 bit 与 10⁶⁹ 条件数作为"必需精度"。

## 2. 代码分类

### 2.1 `src/smatrix_bootstrap/sdp/`（SDPB 路线，37 文件）

| 文件 | 分类 | 说明 |
|---|---|---|
| `grid.py` `hilbert.py` `legendreq.py` `projector.py` `precision.py` `formfactor.py` `constraints.py` `assembly.py` `spec.py` | CORE | 论文 (2.7)(2.9)(3.58–3.75) 的独立推导与装配；cot 核由奇延拓推出并与 (3.67) 逐元一致；Q_ℓ 对 mpmath/Arb/Mathematica；FESR 目标由 (2.50) 重算 |
| `pmp.py` `precision_pmp.py` `sdpb.py` `checkpoint.py` `mma.py` `cli.py` `__main__.py` | CORE | PMP 写出、Arb 装配、SDPB 执行、检查点、Mathematica 核对、入口。正则化块已加入 `pmp.py`；`cli.py` 默认改为 chi-b、mixed-pv |
| `verify.py` `arbaudit.py` `observables.py` `evaluation.py` `accuracy.py` `crosscheck.py` `amplitude_export.py` | CORE | 全原约束复验（float 与 Arb）、相移/η 读取、独立角积分、系数导出。`crosscheck.py` 四处函数内引用旧包核作对照，旧包搬迁后需重指向 |
| `figures.py` `claims.py` `delivery.py` `ir_selection.py` `uv_selection.py` | CORE，需修订 | 图与 C1–C8 判定保留；C4/C7 门槛按冻结的源图比较规则修订；`uv_selection` 的 near 规则与 `ir_selection` 的竖截面规则是代理，主线改为二维最近黑点规则 |
| `sine.py` | EVIDENCE | 声明的替代插值族，非作者离散化；退出主线默认，保留为可选对照 |
| `problem.py` `runner.py` | EVIDENCE | cvxpy/Clarabel/MOSEK 路径；是 M≤40 阶梯（事实 c）的产生者，保留以复现该证据，不在 SDPB 主线 |
| `dual.py` | EVIDENCE | 只对已写出的 PMP 做弱对偶重放，残差必须精确为零才认证，实际从不成立 |
| `report.py` | DISPOSABLE | 已退役的报告生成器桩 |

### 2.2 `src/smatrix_bootstrap/*.py`（Newton 旧包，23 模块，6040 行）

| 文件 | 分类 | 说明 |
|---|---|---|
| `kernels.py`（PV 行、外部 Legendre-Q、区域几何）、`operators.py` 207–312 行（网格与角核）、`model.py` 29–172 行（幺正、运动学、Gram、FESR、FF 算子） | CORE-可复用 | 独立测试过（`tests/test_kernels.py` 14–106、168–194、209–259；`test_mainline.py` 42–89），且是 `sdp/crosscheck.py` 的对照目标 |
| `certificates.py` `imaginary.py`（`support_outer`）、`sampling.py`（`analytic_joint_audit`） | EVIDENCE | Newton 报告里 `joint_audit`/`outer` 块的产生者，用于重验历史证据 |
| `linear.py` `merit.py` `scattering.py` `gauge.py` `ir.py` `probes.py` `quotient.py` `conic.py` | MISLEADING | 自研 barrier-Newton、障碍中心代表点、LP 对偶、带 B 球的 conic 嵌入 |
| `endpoints.py` `analytic.py` `basis.py` `quadrature.py` `spectra.py` | MISLEADING | T0=0、五尾、FF 端点零、analytic-cardinal 处方及其求积；`analytic.py` 的 GMOR 输入审计与局部反例是 EVIDENCE 片段 |
| `analysis.py` `io.py` `run.py` `__init__.py` | MISLEADING | Newton 报告契约、竖截面/凸包边缘选点、17 个 Newton 子命令。只有 `run.py` 前 7 行的 `sdp` 分发与一个可导入的 `__init__` 必须保留 |
| `tests/test_kernels.py` `test_linear.py` `test_mainline.py` | EVIDENCE | 107 项，全部通过，不引用 sdp；随旧包迁入 `legacy/` |

阶段 1 动作：旧包整体 `git mv` 到 `legacy/smatrix_bootstrap_newton/`，保留薄 `src/smatrix_bootstrap/__init__.py` 与 `run.py` 的分发；`sdp/crosscheck.py` 四处引用改为可选（缺失则跳过对照）；`pyproject` 去掉 highspy/scs/clarabel 基础依赖；console script 指向 `smatrix_bootstrap.sdp.__main__:main`。

### 2.3 `scripts/`

| 文件 | 分类 | 说明 |
|---|---|---|
| `scripts/sdpb/sdpb.sh` `scripts/mma/wolfram.sh` `scripts/mma/audit_2309.wls` `scripts/mma/audit_sine_family.wls` | CORE | 启动器与独立公式核对 |
| `scripts/sdp/pmp_verify_blocks.py` | CORE，需更新 | 唯一的 PMP 块与模型逐块对拍；补正则化块与 chi-b 后并入测试 |
| `scripts/sdp/ladder.py` `solver_study.py` `caliber_grid.py` `ff_tolerance.py` `watson_probe.py` | EVIDENCE | 事实 c 的产生者与求解器研究；都用 cvxpy 与 B=377500 |
| `scripts/sdp/drive.py` `pmp_run.py` `finalize.py` | DISPOSABLE | 已不能运行、重复或退役桩 |

## 3. 结果目录分类

### 3.1 CORE（新主线直接依赖）

| 目录 | 内容 |
|---|---|
| `results/runs/sdp_reproduction_20260912` | Clarabel 阶梯与 Fig.4 六 ε（事实 c）；`preregistration.json` 关于 B 的措辞已过时 |
| `results/runs/stage_B_reference_20260906` `stage_C_reference_20260908` | Fig.3/4/7 数字化数据的早期版本 |
| `results/runs/stage_B_pv_M50_L10_prepare_20260907` | PV M50/L10 `amplitude.npz`，SDP 清单引用 240 次 |
| `results/runs/stage_B1_source_method_audit_20260907` | ν₀=0 配点核映射审计 |
| `results/runs/mainline_alignment_20260909`（`D1_hard_M50`、`F_paper/*_current_hard`、χ 指纹、归一化转移） | 电流/FESR/FF 算子与来源审计；其 Newton 结果已被取代 |
| `results/runs/stage_B_delivery_20260907`、`stage_F_mainline_20260909/*_prepare` | 区域索引与五组 (M,L) 算子 |
| `results/evidence/72_midpoint_PV_kernel_construction.md`、`MAJOR_SOURCE_NORM_AUDIT_20260912/` | 作者离散化的重建与 FESR 打包的来源审计 |
| `results/runs/sdpb_regularised_20260913`（符号链接到 `/playpen1`） | 新主线结果根 |

### 3.2 MISLEADING（保留并标注）

`results/runs/sdpb_mainline_20260912`（49 GB，无 B/chi-c/SR-a/sine/ε=0.00021 全线；其中 `foundations_normalization_review`、`optimization_gap_control`、`author_implementation_search`、检查点与导出工具的记录仍有用）、`engineering_20260912`（"去掉 B 才认证"的标题结论；其装配精度、rownorm、PMP 块身份、SDPB 调参 JSON 有效）、`major_claims_20260912`、`physical_consistency_20260911`、`sequential_reproduction_20260910`、`foundations_review_20260910`、`followup_research_20260910`、`expert_response_20260910`、`stage_E_mainline_20260908`、`stage_E_closure_20260909`、`results/evidence/SAME_MODEL_IR_UV_EXCLUSION_20260911.md`、`MAJOR_CLAIM_CERTIFICATES_20260912/`、`results/reports/expert_review_20260910`（数字来自障碍中心）。

### 3.3 EVIDENCE

`sdp_audit_20260912`、`toolchain_decision_20260912`、`audit_20260912`（关于 B 的判断除外）、`claims_alignment_20260909`、`stage_C/D_mainline_20260908`、各 `stage_B3_*`/`stage_B_pv_*`/`stage_B_regulator_protocol`（B 阶梯）、`results/` 根目录 9 月 5–6 日的有限基证书链（`verification_*`、`provenance_*`、`gauge_joint_*` 等）、`results/evidence/` 的推导章节与快照、`references/` 的背景论文。

### 3.4 DISPOSABLE（约 5.9 GB 实占，暂不删）

`pv_remaining_{quotient,rank}_20260906`（2.45 GB）、`stage_B_dense_*`（1.7 GB）、`physical_consistency_20260911/{interlaced_*,fine_prepare,window_prepare}`（1.7 GB）、`sequential_reproduction_20260910/F_*_prepare`（0.6 GB）、`results/analytic_source_storage`（5 GB，sine 族源行）、约 200 个 `stage_B_*` sine 冒烟/中断目录、`results/` 根的可再生 conic/npz/Julia 导出（约 1 GB，均有哈希）。

## 4. `references/`

CORE：`2309.12402v3.*`、`figure*_*.csv` 与 metadata、`manifest.json`、`upstream-gauge-theory-bootstrap/`（作者 2403/2505 代码，事实 a/e 的出处，目前被 `.gitignore` 排除，阶段 1 改为记录其提交哈希并保留副本）、`cordoba-discrete-gtb-2511.11513/` 与 `2511.11513v1.pdf`。EVIDENCE：`1708.06765`、`2203.02421`、两篇 SDP 数值参考。

## 5. 顶层文档

八份文档在 09-12 都被加上同一段横幅："论文没有的密度正则化 B"、"SDPB-only、MOSEK/CVX 不在链路"、"HANDOFF 是唯一入口"。横幅在事实 a–e 下全部失效；讽刺的是各文件的旧正文（Newton 合同 `--density-limit 377500`、`combined-l2` 默认）在这两点上是对的。

| 文件 | 分类 | 阶段 1 处置 |
|---|---|---|
| `PLAN_SDPB_2309_ZH.md` | CORE | 唯一主线参照 |
| `TRIAGE_ZH.md` | CORE | 本文件 |
| `README.md` | 重写 | 一页入口：目标、分支、主线合同、计划与整理文档、运行方式 |
| `HANDOFF_PHYSICS_AUDIT_ZH.md` | MISLEADING | 顶部加"已被 PLAN 取代"说明；保留 Docker/SDPB/MPI、检查点、导出等工程段落的引用 |
| `REPORT_SDP_ZH.md` | 正文 CORE，横幅 MISLEADING | 横幅改为指向本文件 §1 |
| `REVISED_CLAIMS_ZH.md` | EVIDENCE | 加勘误：B、chi-a 建议、"MOSEK 未安装" |
| `STATUS.md` `SCIENCE.md` `REPRODUCTION_GUIDE_ZH.md` `AGENTS.md` | EVIDENCE | 横幅改为指向 PLAN 与 TRIAGE；`SCIENCE.md` §135–284 的参数表与算子推导继续作为参考 |

## 6. 阶段 1 的执行清单

1. 顶层横幅统一替换；README 重写；HANDOFF 降级。
2. 旧包与旧测试迁入 `legacy/`；`crosscheck.py` 引用改可选；`pyproject` 依赖收缩；全部 sdp 测试通过。
3. `sine.py` 与 `problem.py`/`runner.py` 保留但从 CLI 默认与文档主线移出。
4. `references/upstream-gauge-theory-bootstrap` 记录提交哈希并取消 ignore，`engineering_20260912`、`toolchain_decision_20260912` 纳入 git。
5. 不删除任何结果目录；DISPOSABLE 清单交课题组决定。
