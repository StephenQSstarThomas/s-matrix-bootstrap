# Fine模型中心身份与来源审计

结论：fine_tip_support_02→fine_tip_analytic的完整C、全部ImF及有矩界R逐位相同，仅有允许的自由高能R向上提升；fine IR的保存中心、native/candidate、选择和求值使用同一完整C，恢复未改幅度。七个run都没有记录源/输入漂移，全部已登记输入哈希及源码快照均复核匹配。

这里“已收敛中心”只指实际使用的**约化障碍问题的数值中心**。自由高能R的有限可行提升不保证是原4076变量全Gram障碍的唯一中心，也不唯一识别QCD振幅。本审计没有重算Newton/Hessian或证明中心唯一性，没有读新UV相移、运行优化或改动生产文件。

## Tip链

完整C为3876项，原联合点为4076项；实际检查点z为4062项、Gram缓存86块，14个自由高能R不在该数值中心变量中。另有五个端点坐标按原规则固定。

逐项核对 raw joint_path_candidate、最后center_path点、joint.npz、coefficients.json及解析核验输出。raw路径到交付点的活动C/ImF/有界R变化数为0；最后中心μ与保存检查点、末条history及Newton日志一致。解析报告明确引用fine_tip_support_02作为中心与支持来源。

| support→analytic的部分 | 改变项数 | 最大绝对变化 |
|---|---:|---:|
| C | 0 | 0 |
| ImF | 0 | 0 |
| bounded_R | 0 | 0 |
| free_high_R | 12 | 2.57177516089e-09 |

12项R变化都位于FF高能节点且对应矩算子列全零，方向均为增加；不涉及低能/有界R。全部100Gram仍由原式与解析primal审计检查，允许的谱提升不能转写为“完整4076变量点逐位相同”。

tip的原式支持区间为[0.09317173379309758,0.093231374149149601]，gap=5.9640356052e-05；解析阶段重放相同上下界，并通过3582盘、100Gram、2χ、4FESR、14FF、L4、5端点等式和5高能条件。有限H支持与解析primal是两份不同范围的证据，不是全能量证书。

## IR链

center_converged、center_best、center_last所存完整C与native_candidate、candidate、support coefficients、selected系数均逐位相同；C0精确为0，recovery.repaired=False。center_converged记录μ=1.220703125e-07、iteration=29、decrement=1.816513427e-15，has_newton_gradient=True。selection引用该中心文件并保存其SHA和完整C的SHA。

原式支持区间[0.5161863611099164,0.517255595032318]，gap=0.001069233922402；法向y分量为500，纵向边界间隙为2.138467844803e-06，除ε=.002后的相对预算量为0.001069233922402。后两个量不能与未归一化支持gap任意互换；此例数值关系由保存的方向和ε决定。

fine_IR_ref_analytic不生成另一份C文件：它通过输入路径及SHA对fine_IR_ref_02的C独立核验。fine_IR_selected直接复制support C，因此严格的文件关系是“同一support C分别接受解析审计和选择”，不是解析目录产生了新振幅。解析primal通过3582盘、2χ、L4、1个T0等式和5个高能条件；IR没有FF/电流，不能要求它与UV一样有5个端点等式或100Gram。

selection命令参数中的endpoint/tail默认False没有改变模型：实际选择从support记录继承infinity=zero、tail_conditions_applied=True及fine_prepare附加行，并按这些声明重验。不要把消费者CLI默认值当成实际模型签名。

## 求值身份与覆盖范围

两份evaluation.json都引用fine_IR_selected/coefficients.json，相同coefficient_sha256与selection副本均匹配；没有形成另一份幅度。

| 求值 | 能量数×波数 | 物理幺正样本 | 状态与范围 |
|---|---:|---:|---|
| fine_IR_evaluate | 44×3 | 132 | sampled_passed；仅这些给定物理样本 |
| fine_IR_subthreshold | 18×3 | 0 | unresolved指没有适用的物理幺正样本；它求值s=.1–3.75的阈下幅度 |

阈下文件的unresolved不代表上述3582行解析primal失败，也不能反过来把阈下求值称为额外的物理幺正通过。两者都没有检查全部能量或所有遗漏自旋。本次只读取身份、标签和覆盖字段，没有检查新UV相位。

## 源码、输入与附证

七份report均exit_code=0、timeout=False、source_changes=[]、input_changes=[]。每份source.json中21模块的文本SHA与其声明、report.core_sources一致，七份run也是同一producer版本。44项已登记输入引用的SHA全部匹配；其余中心、candidate、selection、evaluation和证明产物的哈希另存于JSON。没有把“目前源码与旧版本不同”混称为当时运行漂移。

无漂移结论限定为已记录并复核的来源。IR选择inputs未单独枚举它读取的所有H/candidate文件；其中心SHA、C SHA及source-report链提供额外绑定，本次也记录了相关当前产物哈希，但不能追认成当时已监控所有未登记读取。当前完整C身份与所用证明的适用性没有发现失配。

补充的[FINE_OPERATOR_PREFIX_REPLAY](FINE_OPERATOR_PREFIX_REPLAY.json)与[FINE_CURRENT_OPERATOR_REPLAY](FINE_CURRENT_OPERATOR_REPLAY.json)已有一致性检查全通过，本次仅检查其记录及输入哈希，没有重新做矩阵比较。[FINE_TIP_SUPPORT_NESTING](FINE_TIP_SUPPORT_NESTING.json)只给两个嵌套的保存float64-H有限程序的+x最优值差上界6.68429115469×10⁻⁵；它既不是相移误差界，也不是精确解析源优化支持的转移。

全部细节、12项谱提升和来源哈希：[FINE_CENTER_IDENTITY_AUDIT.json](FINE_CENTER_IDENTITY_AUDIT.json)。
