# Watson 目标重放与谱系接口的窄代码审阅

2026-09-11。范围：`sampling.replay_support_data/analytic_joint_audit`、`io.watson_lineage`、`gauge` 调用、`calculation_inputs`、evaluate／`_phase_data` 的身份检查及相关回归。审阅期间没有修改生产源码或测试，没有运行优化；只读复核后写入本报告。

**本次接口审阅当时发现的窄代码问题均已修复。** 此历史结论只针对目标、点、对偶和来源的传递，不代表正式 Watson 科学结果已经完成，也不覆盖下文补记的求解器修补。

| 原问题 | 最终核对结果 |
|---|---|
| 解析审计仅保留 projection，Watson 被写成零目标 projection，并丢目标／functional／谱系 | 从源 `objective_kind` 读取固定目标，不受 CLI 默认 projection 误导；原目标传给 `joint_audit` 和 `joint_result`，objective.npz 被原样复制，成功重放保留相应谱系和诊断 |
| objective.npz 初始 witness 可能被混同为已优化对偶 | 支持重验读取源 joint.npz 的最新对偶，包含 Gram／FF／FESR／渐近协向量；不以 objective.npz 的初始 witness 替代它 |
| partial 源的最新对偶因未达到支持 gap 而被清零 | `same and clean` 控制读取最新对偶；`retained` 单独控制已完成支持的继承标记。partial 回归同时验证最新 kR、目标和 metadata 保留 |
| 本次解析 primal 未通过时仍丢 Watson 身份 | 目标读取／复制已前移；所有正常完成的解析审计均记录 `objective_kind/functional/source_coefficients`。失败不进入 selection 生成分支 |
| lineage 继承旧点的中心标记及系数哈希 | 新记录的中心／支持／可行状态来自本次 result；C、current、joint 哈希取当前输出。解析中心保留还要求 `center_acceptance` 通过，即完整 C、ImF、低能谱不变，仅允许自由高能谱增加 |
| 同一目标的审计、续算可能误计为新 Watson 步骤 | 保持原 `functional.source_coefficients`，以该固定旧点的迭代号生成当前号；重复重放同一目标不会再次增加。新步骤才从新源点构造目标；解析重放没有调用目标生成器重算 Q |
| 输入快照遗漏 objective／selection／current／update、直接旧点及实际 profile 消费链 | `io.calculation_inputs` 覆盖参数输入、profile_runs、baseline、相关 current_data，以及 objective 指向的直接旧点；只跟随直接引用，不递归扫描历史 |
| helper 搬迁后漏掉非标准名称的实际源文件 | 参数来源与 objective 引用现均先加入精确 source Path，再展开相邻文件；`center_last_coefficients.json` 等名称有对应回归 |

## 已读取的证据

- [负谱控制的解析审计](watson_negative_spectrum_audit/report.json)为 `analytic_joint_primal_unresolved`、`analytic_primal_feasible=False`、`objective_kind=watson`；其 functional 与源控制相同，objective.npz 字节 hash 相同，且输出没有 selection.json。它证明失败分支保留任务身份，同时没有生成合格代表。
- [goal_00 初始目标控制](fine_watson_goal_00/report.json)是**零预算的一步控制，返回原可行点且不授予新中心资格**。`solver_seconds=0`，但运行时间在步后检查，`report.solver.iterations=1`；返回C仍为原可行incumbent。`joint_feasible=True`、`objective_kind=watson`，但 `solver.representative_center_converged=False`、`support_optimality_certified=False`。此前“零迭代”的措辞已由[预算控制更正](WATSON_BUDGET_CONTROL_CORRECTION.json)纠正，可行性没有被提升为已完成中心或支持证明。
- 已读新增回归：源 objective 中的初始 kR 为−99、joint 中最新 kR 为7，重放取7；partial 同样保留7及固定目标；新谱系不继承旧中心／旧哈希；重复调用迭代号不增加；coefficients 和 current-preparation 均为空的真实 profile 参数形态仍捕获直接输入；非标准源文件也进入快照。
- 三个测试文件均为350行。旧“大于所需下界的自由高能谱必须保持”断言已并入既有完整联合审计检查，原测试体保存在 [WATSON_TEST_CONSOLIDATION.json](WATSON_TEST_CONSOLIDATION.json)。独立目标恒等式和既有数学检查没有被元数据测试取代。
- 主线报告106项全测再次通过；本审阅没有独立重跑该测试集。正 partial／projection 重放及正式 Watson 计算不纳入本报告的完成声明。

## 验收范围

固定旧 Q 所定义的线性目标支持 gap、该目标下 Newton 中心收敛、Watsonian 新旧 S/F/Q 的收敛是不同结论。此次接口保留前两类已有证据及后续诊断所需的引用；不能由支持 gap 推出相位误差、近弹性或固定点收敛。

解析重放仍分别报告解析 primal 与保存 float64 H 的目标支持。高能自由谱提升不改变完整 C、FF 或低能矩预算；任何其它主动变量改变都不能继承原中心资格。lineage 保留初始几何角色，但不会把后续释放投影坐标的 Watson 点重新称为原几何边界代表。

`__init__` 的 lineage、run 的输入采集及 io 的 Fig.11 绘图入口保留薄 wrapper；实现迁移服务于文件限额，没有新增计算入口。本次未扩展一般安全框架或物理主线推理。


## 后续更正：预算控制与中心标记

正式首轮 [fine_watson_01_support_01](fine_watson_01_support_01/report.json)因linear residual约1.08257e−5停止，状态inconclusive，不是不可行证书。主线正在修补：允许经检查确有下降的不精确方向，保持原中心门槛，并保存刚通过检查的点。本次更正时新代码仍在106项测试中，尚未重启；本审阅没有将前次接口测试通过移用于新修补代码。

旧联合实现以步前的小Newton decrement标记步后点。以上谱系及完整C身份核验不等于交付点按新中心规则重验；旧中心不能追称已经通过新标准。其独立primal、support和完整C身份证据仍有效，原记录不改写。
