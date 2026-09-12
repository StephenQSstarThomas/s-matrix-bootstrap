# 领先高能幺正条件实现只读审阅

结论：五项解析行、正尺度变换、障碍梯度/Hessian、原尺度对偶回传及支持常数的符号，与 `ANALYTIC_TAIL_UNITARITY_20260911.md` 一致。审阅发现一个直接阻断原式核验的局部变量覆盖错误，root已在保留真实失败记录后修复。修后没有再发现影响当前 `tail_current → tail_transfer_02 → tail_phase_I_01` 正常硬IR/UV Phase I链的尾部约束丢失。仍有两类旁路/后续交付问题：固定C及Watson续算遗漏保存的尾部对偶，跨分辨率比较的缩减模型签名不区分是否施加这五式。

本报告不判定正在运行的Phase I成功，也不增加尾部物理假设。没有运行优化、改src/tests、执行旧快照或读取相移；仅核对给定推导和实现，并执行短小的数值/接口检查。

## 1. 已修复：对偶字典被度量权重列表覆盖

原 `certificates.joint_audit`用d保存输入对偶字典，之后在 `affine_norm_support`循环中再次给d赋权重列表。新增 `asymptotic_audit(z[:p],M,d)`因此把list传给预期字典的函数，调用 `.get`时报错。打印FESR输入下U>0，确实会执行此覆盖。

真实CLI失败保留在 [tail_transfer/report.json](tail_transfer/report.json)：status为inconclusive，error为`'list' object has no attribute 'get'`。root还在既有M3真实解析算子测试中复现，再把局部列表改名为 `metric_weights`。该变化修复变量生存期，不改变约束、尺度或证书公式。

我检查的后续源码已经包含修复，因此没有执行旧producer复现。独立M3合法原生自由列布局的 `joint_audit`随后完整返回五个尾部裕量；该零振幅/任意电流fixture的joint primal为false，这是预期的接口测试，不能当作物理见证。root报告的完整103项测试由root执行，本审阅不把它计为自己的完整测试运行。

## 2. 未修复的重放遗漏：固定C与Watson续算未载入尾部对偶

`certificates.analytic_fiber_audit`加载joint.npz时仍只读取旧七组键，缺少 `asymptotic_kR` / `asymptotic_kI`。`gauge.gauge_support`中Watson恢复objective及随后重载joint witness，也沿用旧键集合，丢弃已保存的尾部乘子。

这两处**不会关闭尾部primal约束**：`joint_audit`仍依据current元数据计算并检查全部五式；缺失乘子按照合法零对偶处理，所有上界仍重新计算。确定影响是不能原样重放同一个对偶证书，可能损失已取得的支持精度，或失去依赖尾部乘子的排除证据。不能仅据“缺键”声称重算上界无效。当前Phase I不通过这两条加载路径。

应将两组可选尾部键随已有键一起载入、保存和继续；如旧证书无该键，应明确记为零乘子。正常 `gauge dual`、同M `transfer`、`sampling.analytic_joint_audit`的保留支持分支都已加载两组尾部键；不同M转移使用新的零对偶属于明确的数值初始化，不是原证书重放。

## 3. 未修复的模型签名遗漏：跨M/L比较可混淆端点与端点加五式

`certificates.model_signature`没有显式包含 `current_model.asymptotic_unitarity`或FF端点阶数，仍写出遗留 `tail_conditions_applied=False`。正常同模型核验目前依靠不同current preparation路径区分新旧变体，所以当前主求解和同路径重放没有串错。

但 `model.resolution_compare`的缩减签名主动移除了current preparation路径，同时也没有加入这两个物理约束标志。对真实 `endpoint_witness_analytic/report.json`和 `tail_transfer_02/report.json`提取的元数据分别为false/true，完整签名不同，差别只有current preparation路径；现有缩减签名却完全相同。这证明比较身份检查遗漏了新条件。

没有执行完整比较，也不能把这两个同M配置直接混入Fig.11：另有重复配置检查。可达的后续风险是跨不同M/L时，一个配置仅有端点、另一个有五项领先条件，仍被称为“物理输入相同”。应从**current_model的布尔字段**将 `asymptotic_unitarity`及 `ff_endpoint_order`放入完整签名和缩减签名；不要读取source-audit同名的诊断dict来代替物理开关，也不要把五项领先条件混称为已完成全尾部认证。

## 行系数、尺度和对偶核对

`asymptotic_rows`的五个未缩放裕量依次为

\[
3\alpha_1+2\alpha_2,\quad\alpha_2,\quad
(4r+q)-\frac{\pi}{5}(\alpha_1+4\alpha_2)^2,\quad
(r-q)-\frac{\pi}{3}(\alpha_1-\alpha_2)^2,\quad
(r+q)-\frac{\pi}{5}(\alpha_1+\alpha_2)^2.
\]

因此与推导/`asymptotic_margins`最后三项的关系分别为乘5、3、5；符号等价，数值裕量的口径不可混淆。Q的非对角C=2Q约定已正确计入：`2*I`对Q上三角C的系数统一为a_i a_j。前两行R=0，正确表示半直线约束；其障碍生成u=0、v<0，支持常数为0。

`configure_asymptotic_constraints`先将原C行乘ds，再乘完整振幅逆换元 `P.inverse @ P.endpoint_inverse[:p,:p]`，对应同一原始系数；之后应用正尺度w。w使用变换前的实际双密度行选取正数，仅影响数值尺度。若x'=x/√w、y'=y/w，则δ'=δ/w，障碍产生的原尺度乘子为

\[
u=\frac{2\mu x'}{\sqrt w\,\delta'}=\frac{2\mu x}{\delta},\qquad
v=-\frac{2\mu}{w\,\delta'}=-\frac{2\mu}{\delta}<0.
\]

代码的 `asymptotic_duals`实现了这两个除法。独立检查取w=[1,3,17,100,10^8]，最大|Δu|=6.938893903907228×10^−18、Δv=0。对μ=.03125，每个抛物域都满足h−ux−vy=μ，最大浮点偏差约10^−17；其中h=−u²/(2v)。这同时核对了尺度回传、支持常数和残差符号。

`asymptotic_audit`返回残差 `−uR−vI` 与支持h；`joint_audit`把前者加到目标站立性残差、后者加到支持常数，再交给原式支持核验，符号正确。v=0仅u=0允许有限支持，其它情况明确报错。缺少两组乘子时使用(0,0)，合法但不等于保留原证书。

## 硬约束、Phase I、保存和缓存

- 两种problem的state都通过 `augment_asymptotic`加入五个 `−log(2y−x²)`。增量梯度为(2xR−2I)/δ，Hessian分解为该梯度外积加2RRᵀ/δ，新增features/rhs一致。两种slope都包含相同的导数；advance重新按新z计算尾部裕量，固定x的乘子计算也通过完整slope计入尾部。
- Phase I的tau只放松原有UV部分；新增尾部裕量没有加tau，direction的tau列补零。硬IR种子的回退循环用 `asymptotic_margins(candidate,M)>0`筛选，且原解析小λ种子属于给定推导保证的严格内点。有限散射、χ及L4检查仍保留。
- `PhaseOneProblem.duals`给candidate primal检查返回零对偶；这不会漏查五式。真正的原式Farkas调用使用 `original_duals`，从super.duals取得尾部两组乘子，保存在 `phase_I_duals.npz`。普通支持点的joint.npz、joint_path_candidate.npz通过`**duals`完整保留所有键。
- checkpoint只保存z及原有传播值，没有另存尾部state；恢复时 `state`重新由z计算尾部值，所以不需要继承旧的尾部缓存数组。正常路径先核对current元数据与CLI旗标，且旧/新current目录不同；旧缓存不能正常跨路径进入。可额外记录约束版本和尾部开关以加强错误诊断，但未发现当前正常恢复会省略五式。
- `validate_endpoint_model`同时核对端点阶数、函数族、infinity及尾部布尔开关；CLI要求 `--asymptotic-unitarity`必须伴随 `--asymptotic-zeros`。真实tail_current记录阶数2、两开关true、infinity=zero、source_physical_rows=1593，与当前主路径一致。
- `joint_audit`无条件按已启用元数据检查五个解析裕量；`sampling.analytic_joint_audit`把给定独立式加入解析primal检查，并记录数量5与各项裕量。`joint_result`保存完整current_model，正常严格成功分支沿用该模型。一般有限支持和严格原式核验都没有因尾部对偶缺省而放宽primal。
- 候选segment恢复器仍使用较宽的原障碍来寻找恢复比例，但接受前会再次执行包含尾部五式的原式audit；失败则回退已验seed，并由center_acceptance取消改变活动坐标的中心资格。此策略可能保守，不构成尾部约束的静默放宽。

## 范围

上述实现与给定的五项领先必要条件一致。只有在T0=0且五式均严格时，所给推导才支持每个固定保留分波的最终高能正裕量；当前代码保留“未给有限crossover界、未给无限自旋统一界”的范围。≥0、有限采样可行或Phase I进展都不能升级为全能量证明。

## 输入指纹

审阅时间：2026-09-11T06:44:35.607164+00:00。代码指纹对应root已修复metric_weights、当前主计算重启后的版本；原失败producer另列。

| 输入 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/endpoints.py` | `6b4dacd438e750ca61f285744cc0a152d16e3090581e22c77962f23b867005db` |
| `src/smatrix_bootstrap/operators.py` | `5850eb2add0c8a8638b7b92bf4868b4d0a9015e74be3737a0b9f4006ae1a1c0c` |
| `src/smatrix_bootstrap/linear.py` | `ee754804d8e175dfc1940108c2c599d936a70c38f3ee04f0033b7008b1a92e75` |
| `src/smatrix_bootstrap/gauge.py` | `fdc94ef2e59e4bb24ca86b4ce00ff965733b6961ad1c0cf9add96a091a5d5ff7` |
| `src/smatrix_bootstrap/certificates.py` | `7471d68841501b3af52d7cc34cd9c5971bf2c4dd5d7c5c9399ba644c54ecca66` |
| `src/smatrix_bootstrap/sampling.py` | `311d7536c2f75b772afc1597951902ea608e15fefd90483f2337445a370ab188` |
| `src/smatrix_bootstrap/__init__.py` | `9e65f9de80eaceb49fe41feb6d29c8335cafa24f7fa588f46d1b63cd66ff6cbb` |
| `src/smatrix_bootstrap/run.py` | `5c54a004cd77920e946f57904e55e3003593ee6af047240eb2f212561413753a` |
| `src/smatrix_bootstrap/quotient.py` | `3ec58dc689af80cf12b564c74e40cd91af2189cda1789686694cc9691ef8c16d` |
| `src/smatrix_bootstrap/model.py` | `749d522b66404da47357cff2ab93f923f3d8d1024cbdfb210f5c84a74023b8a1` |
| `src/smatrix_bootstrap/imaginary.py` | `fbf4ac702058f5c14733b572c7f3091db9cb65209ca2638968c3578f00a9b488` |
| `results/evidence/ANALYTIC_TAIL_UNITARITY_20260911.md` | `8725bf42554962bd6507249868c0887f02c14df75b3f1f853df033233663ce62` |
| `results/runs/physical_consistency_20260911/tail_current/report.json` | `4bec45338c47bf01d58bd2417df58dffb6af0ab4cac9c63f6fc4faecead76178` |
| `results/runs/physical_consistency_20260911/tail_current/current_data.npz` | `ff8d5809c68515899780eb4793c4d2bb42dbf95d62fdb805d3b60044d071eac9` |
| `results/runs/physical_consistency_20260911/tail_transfer/report.json` | `d7e63f6623c25cfa413237bd7163ffebe004daa9bd5b1aae384f737f9b169f16` |
| `results/runs/physical_consistency_20260911/tail_transfer/source.json` | `70a6ea46566618942c611260a527dfb8479902f95be2e2c3d72b2ae1faaa11a7` |
| `results/runs/physical_consistency_20260911/tail_transfer_02/report.json` | `025064668c4981b3e76c1228197fdc7fa3706baffb5617978e237444960dec72` |
| `results/runs/physical_consistency_20260911/endpoint_witness_analytic/report.json` | `2c75794d6085f113b90c852ebe12e04e3c15b426494f5ee4e249d3fe9005a6cf` |
| `results/runs/physical_consistency_20260911/window_prepare/report.json` | `8884fb8c0a17fc9b5173d33ab96b7b3b8726d2f8850518c5fde8fe4557d944aa` |
| `results/runs/physical_consistency_20260911/window_prepare/amplitude.npz` | `2ab1720aa2c8127dcebef42a14d395c57414e91e6660e26f00684eeb2565c528` |
