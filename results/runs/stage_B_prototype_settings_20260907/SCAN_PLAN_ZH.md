# Sampled 有限原型的区域执行计划

[scan_plan.json](scan_plan.json)包含七个区域各四条初始单位法向，共28个任务；没有启动任何优化。[旧加强计划](../stage_B23_region_scan_plan_20260907/scan_plan.json)保留原样。机器计划中的 `execution.source_snapshot` 留待 root 填入，同时替换环境中的冻结源码路径。

| 共同问题 | 本计划固定值 |
|---|---|
| 准备算子 | `stage_A_M50_L10_20260906`，原50个节点、每isospin 10个分波，共1500个散射盘 |
| 振幅变量 | 完整3876个原始 C_flat 系数，含自由 T0；不删双密度方向 |
| 约束范围 | `unitarity_scope=sampled`、`infinity=free`，无追加物理行、FG约束或五尾约束 |
| 密度界 | actual rho1与rho2上三角的普通 L4，B=377500；这是明确的条件选择 |
| 手征设置 | 两个独立四维L2球；ε=.006/.004/.002/.001/.0006/.0002；pure省去手征球 |
| 数值起步 | barrier，start_mu=1e−5；solver600秒、dual60秒、总cap900秒；BLAS/OMP/MKL均8线程 |
| 对照后端 | 同题 Clarabel/faer、subtracted、SOC；不得借此改变物理约束或移用加强模型的对偶上界 |

有限sine-cardinal插值及两个独立L2手征范数仍是条件定义，不宣称恢复了原文没有明确公布的全部处方。已有准备矩阵直接复用；计划的832-bit/order40参数不会重建或替换缓存算子。

手征warm使用 `stage_B_FG_repaired_M50_v2_20260907/coefficients.json`；pure使用 `stage_B_recover_pure_M50_raw_20260907/coefficients.json`。它们只是已保存并经既往审核的起点，尤其在较小ε下必须重新检查；不导入旧加强模型的工作集或支持证书。

执行前置条件是**同题有限B1的原始坐标支持上下界通过，且源码及条件条款清楚**。先完成ε=.002，其余手征容差及pure继续。`run regions`按共同scope/网格验收支持记录，输出保守内外界距离、下一建议方向和预算，并检查六容差嵌套。原有误差目标保留：度量空间距离≤.01、单方向度量gap≤.0025，同时原始gap≤1e−4；四方向本身不是完成条件。

ROI、离节点、高spin及高能尾的物理诊断和已有违例单独保存并标明对应振幅，不能把其结果伪装成当前施加的有限约束。连续幺正证明不作为这份原型计划的执行门槛；`completion_claimed`仍为false，也不因区域汇总而宣称后续Fig.5–11已经复现。
