# 固定 IR 振幅与独立 UV 供体的短审阅

2026-09-11。只读检查 `gauge.py` 的 dual 分支、`certificates.py::joint_audit` 固定 C 分支、`run.py::calculation_inputs` 及 [SAME_MODEL_IR_UV_EXCLUSION_20260911.md](../../evidence/SAME_MODEL_IR_UV_EXCLUSION_20260911.md)。没有改源码、运行优化或重新积分。

**未发现供体入口或固定 C Farkas 计算的代数／实现错误。排除结论须明确包含当前模型的 FF 二阶端点条件。**

1. **入口正确固定了全部 C。** `dual --fixed-amplitude` 从独立 `interior_coefficients` 供体目录读取 joint 候选及电流对偶，随即以指定 IR JSON 覆盖完整 p 维散射段，目标设为[0,0]。CLI 也禁止固定振幅与非零目标、fixed-x、ray 混用。IR 目录不再需要 joint.npz。已只读核对真实输出：其3876个散射系数与原 IR JSON 逐位相同；原 IR 目录仍没有 joint.npz 或 current.json。

2. **任意供体候选不影响排除证书的量词。** `joint_audit` 对供体 Gram 对偶重新作 PSD 检查／修复，并以当前算子重算 FF／FESR 支持及所有电流残差界。固定 C 分支将完整振幅残差直接代入该 C，清零散射和 χ 的待消元协向量，再求零目标上界。候选 v、R 可供 primal 诊断，并不被假定为该 IR 点的可行电流。负上界若成立，即与任何假设可行扩展必须满足的 `0≤U_fixed(C)` 矛盾。供体原来的非零支持方向不会成为新目标；`support_optimality_certified=False` 也避免把该排除输出当成新支持代表。

3. **自由高能谱没有被偷加上界。** 当前审计将相应 Gram 对偶的末行／列置零以保留 PSD，并要求自由高能谱残差非正；受矩约束的低能谱与 ImF 残差使用原约束导出的有限界。FF 端点依赖分量通过仿射等式消去。负数的大小不是物理距离，整体缩放对偶就能改变它。

4. **必须列明 FF 端点前提。** 实际报告 `current_model.ff_endpoint_order=2`，因此本次允许的电流扩展额外满足两个 FF 在 `z=−1` 的值与一阶导数均为零。固定 C 上界计算确实使用了这些等式及其常数项。证明已有“端点仿射消元”说明，但开头的 UV 条件清单应明确列出这四个等式。结论是该完整 IR C 没有当前加强模型的 UV 扩展，不能凭此证书排除去掉这些端点条件的原2309有限模型，也不排除全部 IR 振幅或所有无ρ振幅。

5. **真实记录与来源范围一致。** [旧入口失败](fixed_IR_dual_before_fix/report.json)确为要求 IR/joint.npz 导致；[原式复核](interlaced_fixed_IR_dual/report.json)给 `−0.1850389032871912`，[解析复核](interlaced_fixed_IR_analytic/report.json)给 `−0.18503890326560488`，两者均为固定 C 排除。第二份重新使用解析行包络、精确运动学／K 和打印输入，不是搬用第一份 H 的结论。已有同模型其它联合可行点与这种单个 C 的排除不矛盾。

`calculation_inputs` 在工作进程消费文件前捕获 H、全部电流准备文件、输入 C 和供体 report/joint 的 hash，结束时以该初始快照合并，供 supervisor 检查变更；无需在 IR 目录制造电流文件。解析复核报告的10项已记录输入 hash 全部与当前文件匹配。原式 `interlaced_fixed_IR_dual` 运行早于这项补丁，其历史 `inputs` 仍为空；新机制不能追认为该旧运行已经记录了读取前快照。该历史局限不应从报告中抹除。

审阅不以106项测试或模块限额替代数学证书，也不据此扩大为全模型不可行、共振必然性或完整 Fig.3–11 交付。
