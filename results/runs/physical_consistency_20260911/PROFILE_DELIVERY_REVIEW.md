# Profile 交付接口小修改审阅

2026-09-11。对照 [修改前快照清单](profile_delivery_before/manifest.json)及其 `__init__.py`、`spectra.py`、`test_kernels.py`，只读检查当前三文件。未修改源码、测试、生产点或物理参数；未运行优化、求值或新相移计算。

**本次发现的两项输出标签问题已在当前代码关闭。** 没有发现阻断当前正常单个 Watson 原生交付路径的其余问题。此结论是接口审阅，不是新的中心、支持、解析或物理结果验收。

| 审阅项 | 结果 |
|---|---|
| 共用 `validate_evaluation_identity` | 三类消费者 `_phase_data`、`direct_profiles`、`plot_profiles` 均在使用已保存振幅结果之前调用；正常 Watson 的 C/current/joint 三类 hash 齐全。 |
| P1 强度与当前谱图 | 强度为同一保存 S 的 `abs(S-1)^2/4`；新增谱图复用 evaluation 中的原生 current diagnostics，不重选点或重算目标。单个带原生诊断的 Watson profile 可以调用现有谱图绘制器；不声称极点、谱饱和或连续认证。 |
| 原生网格被写成新增能量检查 | 已改为 `declared_samples_checked`；scope 只声明所给采样能量，不再声称存在离节点新增能量。 |
| 联合 Watson 被标为 Fig.7／不使用相位选择 | `plot_profiles` 按 `current_preparation` 区分原 IR Fig.7 与保存的 QCD 联合振幅；选择来源标志读取 selection 的相位使用标志或 Watson iteration。 |
| 几何 `gauge-phases` | 原有 exactly-three、同模型、未混合系数、tip/mid/ref 各一及共同网格检查仍保留；单个 Watson 的新增通路没有削弱该入口。 |
| 腐败回归 | 原 Watson 接口测试分别改变 C/current/joint 文件 bytes，要求三个消费者均拒绝；共覆盖九个入口组合。原数学检查仅作等价压行，没有删除公式或断言。 |

## 身份与读取范围

这里只核对已有 producer／evaluation 所登记的 hash，不追认 legacy 缺失 hash 的输入已获得完整身份认证。helper 保留此前 `_phase_data` 的可选字段行为：未登记的 C/current/joint hash 不会凭空补成历史认证。本次正常 Watson lineage 具有三类 hash，没有发现该正常路径缺项。

hash 校验和后续 JSON 读取不是原子文件系统快照。既有 CLI 的开始／结束输入核对能够拒绝持续改写的输入；本次没有输入改动或混读证据，也没有扩展原子文件系统架构。当前科学交付继续使用冻结 producer，验收以最终 run report 的来源／输入状态为准。

[前一次原生交付重放](profile_delivery_native_replay/report.json)已完成，`exit_code=0`、无超时、`source_changes=[]`、`input_changes=[]`；它保留修补前的旧标签，不能当作新标签版本的重放。主线正在复验完整 106 项测试，并将用新目录完成最终 single-native replay。本审阅没有独立重跑该套测试，也没有将这些待完成检查写成已通过。

## 已核对的文件

修改前的三文件 SHA-256 均匹配快照清单；当前三文件通过语法解析。当前行数为 337／270／350，两个核心文件分别为 23059／22007 bytes，符合现有文件限额。

| 当前文件 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/__init__.py` | `eab9d89f9c5aec186369fb48cef62f20db87a20a16512f8c8a6e7349a2a4473a` |
| `src/smatrix_bootstrap/spectra.py` | `10b7bc464cb216b00674c95fa4d1218543a7f8737a0707b20ea7e2b55d53de95` |
| `tests/test_kernels.py` | `244f62c280f75853ccdc9ec09be632a37f2ea3c6c243ced58f9049c9787b6615` |

原始审阅副本现以只读 .py.txt 数据保存；字节SHA未变，路径映射见 `CURRENT_REVIEW_SNAPSHOT_INDEX.json`。它们不是新的可执行模块。
