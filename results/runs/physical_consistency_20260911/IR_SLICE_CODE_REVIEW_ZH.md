# IR-only 解析切片通路的独立代码审阅

日期：2026-09-11。范围：`ir.py`、`scattering.py` 的 barrier／selection／capture，`certificates.py::analytic_scattering_audit`，CLI 与 UV 图的 IR baseline 身份；同时只读核对相关坐标、渐近行及采样函数。本次未改生产源码、未跑优化、未启动子代理。主线在审阅期间自行修正了以下已反馈问题。

## 发现及修复状态

| 级别 | 实际问题与触发 | 本次最后核对状态 |
|---|---|---|
| P1 | `coherent_boundary` 在赋值前读取局部 `prescription`；IR 使用 `--asymptotic-zeros` 即触发 `UnboundLocalError` | **已修复**：检查已移到 `support_data/source_contract` 之后 |
| P1 | 零 T₀、默认 `dual_seconds>0`、barrier 尚未达到 gap 时调用 `normal_dual(infinity_zero=True)`；后者明确只支持 sampled/free，抛出的 `ValueError` 未被捕获，使已取得的 IR 候选运行变成 inconclusive | **已修复**：零 T₀ 切片明确跳过不适用的 LP fallback，保留其已变换的 barrier covectors |
| P2 | IR selector 只重验当前 C 的可行性，却沿用 report 的 center 标志、targets、lower/upper；同目录另一个可行 C 可被误标为原受支持中心 | **已修复**：已保存 `center_converged.npz` 的完整 C，selector 对 C／candidate／实际中心逐位核对，并重算该 C 的目标、fixed-x 与支持距离 |
| P2 | IR support 未将读取的 preparation/report.json 与 amplitude.npz 加入 `inputs`；supervisor 无法检测3123行算子变化，单靠路径不能鉴定约束身份 | **已修复**：运行开始保存两个输入 hash；零 T₀ selector 要求并复核这些 hash |
| P2 | 下游评价原先可失去中心身份：evaluate 直接附加 selection.json，baseline 仅核对模型字段；选点后替换 C，可能让另一个同模型可行振幅沿用旧中心身份 | **已修复并复核（所有新 IR selection）**：selection 保存完整系数文件 SHA256；evaluate 比较并记录该 hash；`_phase_data` 核对 evaluation／selection／当前 C 三者 hash |

最后一项没有表示本轮已有振幅被替换。最新代码已闭合所指出的具体错配：新 selector 在完成 C／candidate／中心逐位核对后写入 `coefficient_sha256`；evaluate 开始记录当前 C 的 hash，比较 selection 的绑定，并将 C 与 selection 加入 `report.inputs`，同时把 C hash 写入评价；IR baseline 经 `_phase_data` 再比较 evaluation／selection／当前 C。因此选点后或求值后替换 C 都会被拒绝。旧 selection 缺少 hash 时仍保留历史兼容路径，不能追认为具有这项新中心绑定。本次列出的五项缺陷，在各自说明的范围内均已修复；未扩展一般安全框架。

## 数学与数值通路核对

**零 T₀ 换元未发现错误。** 原 C_flat 先除非对角 Q 的打包因子，再以原生 S0/S2 的虚部替代单谱自由坐标，随后乘对应 κ；它们在保留幺正盘中均位于[0,2]。删除的仅是 T₀ 坐标，所有密度方向仍保留。`P.raw` 使用完整逆变换和打包因子，并精确返回 T₀=0。新增采样不改变原生前缀，因此锚点索引仍正确。

用既有 M3/L2 准备算子做小规模变换核对，得到 C 往返最大误差 `2.17×10⁻¹⁸`、H 前向值差 `8.35×10⁻¹⁸`；T₀ 精确为0。该检查仅检验数值变换，没有重新优化或产生新的物理见证。

**五条渐近抛物约束与 joint 共用定义。** 原变量中为 `2 I·C−(R·C)²≥0`；数值缩放采用 `x=(R·C)/√w`、`y=(I·C)/w`，不改变可行集。其 barrier 梯度、Hessian 特征及方向导数一致。对偶恢复使用 `u/√w`、`v/w`，并保留 `v<0` 时抛物域支持 `−u²/(2v)`；`audit_ir` 把相应原 C_flat 残差和支持常数交给原支持核验，没有漏掉缩放。

既有严格首谐波点的小检查中，原／变换变量的对偶配对差 `6.79×10⁻¹⁹`，支持常数差 `7.68×10⁻²⁴`。这些是浮点一致性检查；严谨支持仍由 Arb 对原保存 H 及精确渐近行重新核验。

主线新增的吸收坐标梯度误差宽度

\[
2\sum_{j<\mathrm{free}}|r_j|+2B\|r_{\rho}\|_{4/3}
\]

对自由坐标盒[0,2]和半径 B 的实际密度 L4 球是有效的保守宽度，修正了只看 Newton decrement 而不控制支持残差的缺口。它是求解器停机控制，不能替代最终支持审计。

**可行恢复与中心标记目前一致。** `recover_ir` 先验证返回 C，已经通过时逐位保留；否则只沿已认证种子到候选的凸线段寻找恢复点，并对实际舍入后的完整 C 再做审计。任何恢复都撤销 `representative_center_converged`。线段搜索的浮点可行性只负责候选搜索，不被当成最后证书。

## 采样、审计与物理范围

- `prepared_sampling` 保留原生1500盘前缀并读取完整附加能量／波标签。barrier、χ行偏移与 `audit_ir` 使用总行数 n，未发现把3123盘截回1500盘的索引错误。
- `analytic_scattering_audit` 对全部 n 行使用系数包络；原生点用精确 midpoint 能量，额外点用所声明 binary64 能量。T₀、五条渐近条件、两条 χ 球及密度界均独立检查；不读取 FF、Gram 或 FESR 输入。
- IR baseline 的 M/L、准备路径、解析处方、密度界、χ、T₀ 和渐近条件与 joint 一致时，比较的是同一散射约束切片。`current.json` 不存在的检查防止直接误用 joint 文件夹；新 IR selection 的完整 C 身份另由上述中心逐位核对及下游 hash 链保留。
- 保存 float64 H 的支持界与含解析行包络的 primal 结论仍分开，这是正确的范围区分。五条高能条件是必要领先系数约束，3123盘是有限采样；两者都不自动认证阈值邻域、中间全部能量或无限自旋。

审阅采用 `PYTHONPATH=/tmp/collocation_arb:src /home/shiqiu/miniconda3/bin/python` 运行上述小矩阵算术。没有把主线106项测试或小 CLI 结果冒称为本审阅独立执行的测试，也未据此宣布新的 IR 科学代表已经交付。
