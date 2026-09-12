# Root已审阅并应用十模块整理

2026-09-10：已切换到repo，旧15模块和三测试完整保存在[证据](../../../evidence/core_before_ten_modules_20260910.json.gz)。实际79项测试通过；新代码重验原M50联合点得到完全相同支持区间，并重新生成全部1500原生分波后通过幺正性检查。[应用记录](APPLIED.json)、[实际验证](verification.json)。以下保留隔离阶段的原审查记录，其“未切换”描述适用于当时。

# 十模块暂存整理：逐文件限额与验证已通过，待 root 审阅

最终源码位于 `/tmp/smatrix_core10_review/src/smatrix_bootstrap`，恰好十个真实模块。全部文件满足 350 行 / 24 KiB；三测试文件满足各 350 行。**没有切换或修改仓库当前 src/tests，没有追加科学优化。** 原生产源码/测试逐字节未变。

[完整 diff](CORE10_DIFF.patch) 对照原 15 模块与三测试；[本轮限额收束 diff](LIMIT_CLOSURE_DIFF.patch) 对照上次已审阅的十模块草案。函数迁移清单见 `INTERFACE_MAP.json`。所有计算都在源码中，没有 exec/eval、代理模块或数据中的可执行实现。

## 最终尺寸

| 文件 | 行数 | 字节 | 通过 |
|---|---:|---:|---|
| `__init__.py` | 348 | 22883 | 是 |
| `analysis.py` | 319 | 24412 | 是 |
| `imaginary.py` | 314 | 23839 | 是 |
| `io.py` | 290 | 24313 | 是 |
| `kernels.py` | 350 | 24412 | 是 |
| `linear.py` | 328 | 24147 | 是 |
| `model.py` | 350 | 24549 | 是 |
| `operators.py` | 346 | 24505 | 是 |
| `quotient.py` | 330 | 24236 | 是 |
| `run.py` | 342 | 24198 | 是 |

合计 **3317 行、241494 字节**，较原 4609 行、327199 字节减少 1292 行、85705 字节。机器记录为 `FILE_LIMITS.json`。测试文件为 229 / 254 / 219 行。

## 收束范围与职责

- 保留完整 PV A–F、全部幅度/当前变量、原式 audit、B normal_dual、D Clarabel 联合初始化、E 完整 Newton / 具名 Watsonian、F 完整 Phase I / 分辨率运输及所有 Fig.3–11 主图、eta 图和原图数值比较。
- PV 已直接继承 object，保留 canonical row ID、全部原生物理行、阈下行、常数坐标和完整 C_flat；没有删掉旧 sine 父类而遗失其公共方法。
- I/O 基础与共享报告构造位于包入口；物理输入/曲线比较位于 model，主图与预声明选点位于 analysis/io。PSD、原式审计与坐标变换归入相应数值模块。kernels 中没有绘图代码。
- 重复 imports 合并；部分已由报告字段表达的说明文字缩短，原说明保存在 `FUNCTION_CONTRACTS.md` 和原快照。没有新增多语句挤行，没有删数值字段以达限额。
- `center_report`、`phase_one_report`、联合/区域/选点/F 报告构造已抽取，数据仍使用原问题的 M、全部维数与原验收状态。新增一个回归检查覆盖报告 helpers 的 4076 / 4063 维记录，避免函数抽取后引用外层 M。

主要接口：JointProblem/PhaseOneProblem 在 operators，完整 Newton 在 linear；joint_hull_candidate 与 normal_dual 在 quotient；joint_audit 在 kernels，support_outer 与 Watsonian 目标在 imaginary；current_operators 与 compare/resolution 比较在 model。统一 CLI 子命令保留，完整迁移以 `INTERFACE_MAP.json` 为准。

## 数学验证：已经重做，旧整函数宣称撤回

最后一次从暂存 cwd 运行，并明确确认包来自暂存：

```bash
cd /tmp/smatrix_core10_review
PYTHONPATH=/tmp/collocation_arb:/tmp/smatrix_core10_review/src \
SMATRIX_WORKSPACE=/home/shiqiu/s-matrix-bootstrap OPENBLAS_NUM_THREADS=1 \
/home/shiqiu/miniconda3/bin/python -m pytest -q --disable-warnings
```

**79 passed in 10.55s**：原 78 个移植后的独立检查通过，新增报告构造维数回归检查通过。见 `pytest.log`。测试包括独立角积分、原生跳跃、完整实际密度包装/坐标改变、手征两种分组、Gram/FESR/FF、原式支持、Phase I Farkas、Watsonian 恒等式与分辨率运输。

`MATH_BODY_COMPARISON.json` 是重新生成的 18 项 AST 检查，全部为 true。比对忽略 imports/docstrings，并仅还原已列明的文字缩短；变量、常数、运算、数值控制流继续严格比较：

- B 内层 Newton 循环与径向中心化；E 内层 Newton 循环。
- JointProblem 全部方法；barrier_terms 完整函数。
- PhaseOneProblem 的初始化、values/state/delta/advance/slope/raw/original_duals。
- joint_audit / current_operators 的完整计算，以及三种完整坐标变换。

**没有继续声称整个 center_joint_support AST 不变。** 它的外围报告构造已抽取；Phase I 的零对偶构造也使用共享 helper。Newton 的 QR/CG、线搜索、μ 停止/下降、φ/gradient/Hessian、全部约束与原式证书均未改变。短文字变化与外围改动在 JSON 中明确列出。`STATIC_GLOBALS_CHECK.json` 没有未解析的全局变量。

## F 缺数据报告：四组真实输入已集成

使用 root 提供的四个已完成原生求值，没有构造小模型代替它们：

- `F_paper/M50_L8_evaluate`
- `E3_paper_tip_evaluate`
- `F_paper/M45_L10_evaluate`
- `F_paper/M60_L10_evaluate`

最终暂存代码产物位于 `integration_partial_F_final/`：

- `resolution_partial.json`
- `fig11_partial.pdf/png`
- `fig11_eta_partial.pdf/png`
- `integration_checks.json`

结果为 `status=incomplete_resolution_comparison`、`complete=false`、`missing_configurations=[[50,12]]`。四组共 **537** 项已保存原生主波检查。固定 L=10 的三 M 子比较完整；固定 M=50 的 L 子比较不完整，图标题分别标识。没有 L12 曲线、虚构见证或全 Fig.11 完成宣称。

五组全集仍是完成条件。每个提供点仍须完整 +x 原式支持，gap 至多 1e-4；参数给出更松 gap 也不会放宽这个比较门槛。实际输入的重复配置、未知配置及人为置为 2e-4 的支持间隙均被拒绝，记录见 integration_checks。数据/图片只写入暂存目录。

## 退役代码和额外图片

此前已批准退役：finite-sine 振幅、五尾、FG/多项式矩锥、旧全空间 conic、CG checkpoint 导入、embed、regulators。保持 separate/combined、hard/clipped、printed/eq250 与 raw/normalized 的具名审计。原证明输入与失败没有转移为 PV 验收。

本轮依 root 转达的用户要求，另退役三类非论文主图的绘制：

| 不再生成的图片 | 保留的数据 |
|---|---|
| compare 的 `*_differences.pdf/png` | `comparison.json`、`summary.csv`、原图对应与差值统计 |
| `ff_guided_phases.pdf/png` | phases/evaluation 中 FF-guided lift、nearest phase、原始 S/F 与来源 |
| `watson_diagnostics.pdf/png` | 全部 Watson/两 pion 谱分数/投影运动数值与 lineage |

Fig.3–11 主图、eta 图、原图曲线及数值比较仍保留。旧运行已有图片不删除。最终 `RETIRED_OUTPUTS.json` 列明范围。

## 旧 PV 数据与 CLI 重放

现有 coefficients/joint/current/selection/report、barrier_state 和 phase_I_state 仍可读。已只读检查 M50 ref 的 3876 / 4076 / 4062 数据，以及 M45 Phase I 缓存；完整向量前缀相同，缓存字段保留。见 `COMPATIBILITY_CHECK.json`。`--help` 已通过；root 计划的正式 M50 audit/原生求值与正式 F partial 尚待其审阅后执行。

旧报告不能把 parameters 全字典机械翻成 CLI。重放需移除以下退休参数：

- `embed` / `regulators`；`--prescription finite-sine-cardinal`；`--unitarity-scope strengthened`；`--infinity zero`。
- `--density-cone`；旧 `--coefficient-key raw_candidate/working_candidate/incumbent`。使用已经声明坐标的完整 PV C_flat JSON。
- pure/chiral 的旧全空间 `--solver clarabel`。Clarabel 仍保留 D 的 gauge `--joint-feasibility`。
- 非原生 PV 采样（sampling_factor≠1、非空 unitarity_energies/additional_constraints、scan_nodes≠0）不能静默重放为同一有限问题。

`--native-backend` 仍用于 D；完整数值坐标接口保留；angular_order=24、sampling_factor=1、scan_nodes=0 等无效旧默认字段继续接收。separate/clipped 必须明确指定，不能继承为 combined/hard 的主线证书。计算时继续明确指定当前 PV preparation/current-preparation，保留路径身份检查。

## 证据保留与切换边界

`BASELINE_PROOF_SOURCE.json.gz` 保存原 15 模块和三测试的完整非执行源码快照；`RETIRED_SOURCE.json`、`RETIRED_TESTS.json` 及 `RETIRED_OUTPUTS.json` 列明退休/移植范围。`ACTIVE_SOURCE_UNCHANGED.json` 全部为 true。

本次整理验收完成，**生产切换留给 root review**。没有通过新优化、参数变化、放宽物理约束或省略未完成的 L12 来完成源整理。
