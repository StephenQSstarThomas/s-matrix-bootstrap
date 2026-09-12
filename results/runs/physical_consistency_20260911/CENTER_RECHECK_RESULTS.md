# 八个冻结联合中心的新检查点标准重验

2026-09-11。独立验证结果；输入来自 [CENTER_RECHECK_MANIFEST.json](CENTER_RECHECK_MANIFEST.json)，全部明细、命令、来源与输出哈希见 [CENTER_RECHECK_RESULTS.json](CENTER_RECHECK_RESULTS.json)。

**八项全部通过：零 Newton 迭代，有 `joint_stationary_point` 且 `point_updated_after_check=False`；中心、支持和原式 primal 同时通过。** 每项新旧完整 C、全部 ImF、有矩界 R 逐值相同。自由高 R 也全部逐值相同，本轮没有任何自由谱提升变化。完整 joint 点的逐值身份亦为 8/8 相同。

每次仅一个 `boundary --mode gauge --solver centered --solver-seconds 0 --require-center` CLI，`SMATRIX_BLAS_THREADS=1`，worker 确认为 BLAS 1 线程；八项依次执行，每项监督上限 180 秒。显式使用检查点实际 `--start-mu`，保留原物理标志、原 direction/fixed_x/gap/χ障碍权重。fine ref/mid 始终使用原 `[0,500]` 目标及原固定 x，没有用最终对偶支持平面方向替代。未查看相移、调参、修改源码或旧点，也没有续跑任何检查。

| 原 producer／重验报告 | 检查点 μ | Newton 迭代 | 新检查 decrement | 线性相对残差 | 本次支持 gap | CLI 用时（秒） | 冻结点新标准 |
|---|---:|---:|---:|---:|---:|---:|---|
| [fine_tip_support_02](center_recheck_fine_tip_support_02/report.json) | 6.103515625e-09 | 0 | 4.13076e-17 | 1.97125e-11 | 5.96403561e-05 | 119.36 | 通过 |
| [fine_UV_ref_support_02](center_recheck_fine_UV_ref_support_02/report.json) | 1.953125e-07 | 0 | 1.32843e-17 | 6.31389e-08 | 0.0017952968 | 85.83 | 通过 |
| [fine_UV_mid_support_02](center_recheck_fine_UV_mid_support_02/report.json) | 1.953125e-07 | 0 | 2.29225e-17 | 2.07568e-11 | 0.00182350887 | 82.59 | 通过 |
| [F_M50_L8_support_05](center_recheck_F_M50_L8_support_05/report.json) | 1.220703125e-08 | 0 | 9.21694e-18 | 9.40619e-08 | 6.86247176e-05 | 24.11 | 通过 |
| [E1_tip_center](center_recheck_E1_tip_center/report.json) | 1.220703125e-08 | 0 | 2.04893e-17 | 5.16878e-08 | 7.43718037e-05 | 29.74 | 通过 |
| [F_M50_L12_support_02](center_recheck_F_M50_L12_support_02/report.json) | 1.220703125e-08 | 0 | 3.31488e-17 | 9.46688e-08 | 7.91944322e-05 | 36.21 | 通过 |
| [F_M45_L10_support_01](center_recheck_F_M45_L10_support_01/report.json) | 1.220703125e-08 | 0 | 1.32993e-17 | 5.35702e-08 | 6.38755154e-05 | 20.09 | 通过 |
| [F_M60_L10_support_02](center_recheck_F_M60_L10_support_02/report.json) | 1.220703125e-08 | 0 | 3.09642e-17 | 6.7876e-08 | 9.62321849e-05 | 141.27 | 通过 |

八次均在各自 180 秒上限内完成，没有超时。原 manifest 的 133 个文件哈希在重验后全部保持，所有新 run 已登记输入的当前哈希也匹配；各 run 的 `source_changes`、`input_changes` 均为空，八次保存的核心源码哈希一致。F 五配置沿原比较文件的明确中心链绑定，M50/L10 使用原 `E1_tip_center`，没有被新 fine tip 替代。

固定谱的划分来自对应 `current_data.npz` 的高能节点索引，并逐项与实际 `moment_linear.npz` 的零／非零谱列核对。有矩界 R 数量为 M50 的 86、M45 的 76、M60 的 102；自由 R 数量分别为 14、14、18。JSON 单列全部自由谱索引、变化列表及最大差值，八项变化列表均为空。

这次通过是对保存点按新规则实际完成的核查，不是把旧“步前 decrement 标记步后点”追认为同一规则。结论只针对这些冻结点在实际 μ、原目标及约化障碍下的数值中心资格，不证明原完整全 Gram 障碍中心唯一，也不识别唯一 QCD 振幅。boundary 中的原式 finite-H primal/support 已重验；本次没有另外运行解析 source-audit，既有同 C 的解析证据保留其原范围。未作 continuum 认证或新的物理曲线结论。

原 manifest、原生产结果、PROGRESS、STATUS 等均未改动；新结果仅归档在上述独立检查目录及本汇总。
