# 指定原输入的 PV/midpoint 历史清单

2026-09-11。机器清单：[PV_SAME_INPUTS_INVENTORY.json](PV_SAME_INPUTS_INVENTORY.json)。只读取本轮已定位的 STATUS／SCIENCE 历史链及其指定文件；没有追加目录扫描、读取压缩旧源码、恢复 legacy 模块、改动原文件或运行科学优化。

**可确认：这组精确设置已有完整联合可行见证及一个 ref 支持，但尚未定位到同设置完整三代表的独立选择／绝对相位交付链。** 清单中的“通过”来自历史 producer 的原式报告；本次核对文件身份、参数、完整数组与保存状态，没有重新运行这些证书。

## 1. 匹配的固定设置及算子

所有下列匹配项采用：mixed `pv-midpoint`、M50/L10、1500个散射盘、3876维 C、两个 separate-L2 χ 球（ε=.002）、printed 四矩、每矩 raw绝对容差.002、hard-midpoint、实际双密度 L4 上限377500、`infinity=free`，没有新增端点／高能必要条件或额外采样。原 FF 归一化与高能平方界保持：F₀(0)=F₁(0)=1，向量 cap=3×10⁻⁵，标量 cap≈1.954438775510204×10⁻⁷。

仓库根为 `/home/shiqiu/s-matrix-bootstrap`。JSON 为每个文件同时保存仓库相对路径、绝对路径、是否存在、字节数和 SHA256。

| 用途 | 仓库相对目录 | 确切入口 |
|---|---|---|
| H准备 | `results/runs/stage_B_pv_M50_L10_prepare_20260907` | [report.json](/home/shiqiu/s-matrix-bootstrap/results/runs/stage_B_pv_M50_L10_prepare_20260907/report.json)、[amplitude.npz](/home/shiqiu/s-matrix-bootstrap/results/runs/stage_B_pv_M50_L10_prepare_20260907/amplitude.npz) |
| 电流准备 | `results/runs/mainline_alignment_20260909/D1_hard_M50` | [report.json](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D1_hard_M50/report.json)、[current_data.npz](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D1_hard_M50/current_data.npz) |
| 旧设置及几何规则 | `results/runs/mainline_alignment_20260909` | [HARD_MAINLINE.json](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/HARD_MAINLINE.json) |

电流准备还包括同目录 `gram_linear.npz`、`moment_linear.npz`、`ff_cap_linear.npz`，全部已记录哈希。H为3011×3876：实部1500行、虚部1500行及11辅助行。其原报告明确 `analytic_function_family_defined=False`、`integral_enclosures=False`、原生PV cot核及节点跳跃、原生点以外未定义物理求值。这是本仓库声明的混合配点实现，不是后来 `analytic-cardinal` 1500盘基线，也不冒称恢复了作者全部未公开离散细节。

旧 producer 早于新增 `asymptotic_zeros/ff_endpoint_order` 等接口；JSON保留这些原字段的缺失状态，不伪造为当时显式写出的False。无新增条件的判断依据为原H的 `infinity_equality_imposed=False`、`asymptotic_necessary_conditions=False`、sampled/free合同与保存的原电流模型。

## 2. 完整可复用数据点与支持范围

三个 producer 均位于 `results/runs/mainline_alignment_20260909/`。

| Producer | 保存的 (f₀⁰(3),f₁¹(3)) | 已有范围 |
|---|---|---|
| [D3_hard_hull](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D3_hard_hull/report.json) | (.060683815494175064, −.0038092113023219737) | 完整联合可行见证；目标[0,0]，保存幅度凸包上的Clarabel电流SDP；不是全空间支持最优代表 |
| [E_tip_path_03](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_tip_path_03/report.json) | (.0971070800504441, −.006266874434231021) | 联合可行；+x支持区间[.09710708005043221,.10032306820041371]，gap=.0032159881499815024，原报告 `support_optimality_certified=False` |
| [E_ref_path_03](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_ref_path_03/report.json) | (.07298305467748878, −.004291624605377557) | 联合可行；方向[35.430585138789816,500]的支持区间[.4400200297491743,.4408397133258436]，gap=.0008196835766692834，通过其原设gap=.001 |

确切点文件：

- D3：[C](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D3_hard_hull/coefficients.json)、[current](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D3_hard_hull/current.json)、[joint](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/D3_hard_hull/joint.npz)。
- tip：[C](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_tip_path_03/coefficients.json)、[current](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_tip_path_03/current.json)、[joint](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_tip_path_03/joint.npz)。
- ref：[C](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_ref_path_03/coefficients.json)、[current](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_ref_path_03/current.json)、[joint](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E_ref_path_03/joint.npz)。

三项均已逐位核对：C长3876、joint长4076；current中的两组50维ImF与两组50维R分别等于joint对应段。其所有文件哈希见机器清单。

**旧ref是输出辅助坐标。** 原请求 `fixed_x=.07298305468121799` 来自HARD_MAINLINE中Fig.8两个near-black标记相对黑点的比例，再以140/92参考缩放；没有使用相移，但仍使用了论文输出。下一轮独立选择不会继承该坐标。tip的角色规则是完整+x方向；D3不承担几何代表角色。

## 3. 实际检查点 μ 与历史中心区别

| Producer | 原报告 start_mu | 实际 barrier_state.npz μ | 最后历史 `central_converged=True` μ |
|---|---:|---:|---:|
| D3_hard_hull | .01（通用参数） | 不存在该缓存；不能据通用参数推定实际Newton μ | 不适用：Clarabel可行性见证 |
| E_tip_path_03 | 1.5625×10⁻⁵ | **4.8828125×10⁻⁷** | 9.765625×10⁻⁷ |
| E_ref_path_03 | 3.90625×10⁻⁶ | **1.220703125×10⁻⁷** | 2.44140625×10⁻⁷ |

两份Newton缓存均在各producer目录的 `barrier_state.npz`，保存z形状为4062，包含mu、reference_point与约束状态；确切绝对路径、哈希、字段及μ的float十六进制表示均已记录。最后一条history都标为未收敛，不能把缓存μ当成已经完成的中心μ。

两份 `center_path.npz` 各有5个历史中心。此次逐位比较发现，**各producer最终交付的joint以及C+ImF前缀均不与其5个历史中心中的任何一个逐位相同**。这不撤销其原式可行／支持范围，但禁止直接把交付C挂到这些历史中心标签上。缓存存在也不等于当前求解器可以免检恢复；恢复前仍须核对数值坐标和同一算子身份。

ref的 `center_observables.json` 最后一对历史中心（μ=4.8828125×10⁻⁷→2.44140625×10⁻⁷）报告P1最大模π相位变化.15678308190195234°、新中心原生最小η=.2326442055079171。它是中心间诊断，**不是最终交付C的完整绝对相位曲线或新中心认证**。

## 4. 完整相位交付的设置差别

- STATUS链接的 [E3_paper_phases](/home/shiqiu/s-matrix-bootstrap/results/runs/mainline_alignment_20260909/E3_paper_phases/phases.json)具有完整三点相位，但采用 **combined-L2**；其ref/mid仍使用Fig.8标记比例，历史范数选择也参考过Fig.5输出指纹。
- 更早 [E3_phases_final](/home/shiqiu/s-matrix-bootstrap/results/runs/stage_E_mainline_20260908/E3_phases_final/phases.json)采用 **clipped-phi** 电流准备；signature未显式保存chiral_norm，伴随文档说明两个χ球。它不满足本次指定的hard-midpoint。
- 已定位的hard/separate tip/ref目录没有独立selection.json或evaluation.json；这条已命名历史链尚不能被记成所问设置的完整三代表相位交付。此处没有声称已穷尽整个仓库的所有未知目录。

## 5. 下一主线决定与保留范围

按根线程已定次序：先完成当前固定Watson goal1，再保留上述原H和current准备不变，补这个mixed-PV相同输入下的独立几何三代表链；暂不进入Watson goal2。旧ref的论文标记坐标不继承。

这些完整C/current/joint可以作为随后重新核验的数值数据；不得把旧支持最优性移交新的代表、范数或截止。尤其不能导入analytic-cardinal的T₀/FF端点等式、领先高能约束或连续认证前提。所有旧文件、失败检查点和历史标签均保持原样。
