# 端点实现有界只读审阅

审阅结论：当前显式 `endpoint_current` → `resolution --joint-feasibility` → 正常联合/解析核验主路径，未发现会丢失端点条件或使已实施的仿射换元、对偶常数失效的缺陷。当前输入明确为双零点阶数2、`asymptotic_zeros=true`、`infinity=zero`，新增取样总数1593。此结论是代码与小规模代数核对，不代表正在运行的Phase I已经成功，也不代表存在全能量物理解。

有两项固定振幅严格重放问题，应在启用该旁路前修复；另有一项可构造的阶数身份错误。第1项在完整控制流中最终被维度检查拒绝，不能称为已证误接受。审阅期间未修改src/tests、未启动优化、未查看相移，未中断root的主运行。原算法只作为非执行源码文本对照。

## 1. 中优先级：增行H的固定振幅解析重放抽行错误，最终被布局检查拒绝

位置：`certificates.analytic_fiber_audit`约229–233行；`joint_audit`约31、42、56、60、63行。

该入口读取 `amplitude.npz['energies']` 并重复 `3L` 次，仍得到旧的1500个原生取样。`joint_audit`立即用 `n=len(kap)` 抽取 `H[n+ids]` 为原生电流Gram所需的Im f行。然而当前H的真实虚部起点是1593。虽然 `analytic_preparation` 分支随后利用 `prepared_sampling` 重建正确kap，错误的rows及其误差包络已经形成。

对真实 `window_prepare/amplitude.npz` 的只读检查得到：H形状为3197×3876，旧n=1500，正确n=1593；100个原生电流Gram的虚部行全部不同，最大系数差为8.36316931860088。这是抽行错误，不是舍入误差。

完整控制流中的最终行为是拒绝，必须与抽行错误区分。`joint_audit`末端的 `if fixed_amplitude`会执行 `kr=ki=np.zeros(n)`，把保存的扩展对偶重新覆盖成入口旧n=1500的长度；此时analytic分支已把kap重建为1593，`support_outer`的布局检查随后抛出 `ValueError: Finite complete native scattering layout required`。按此真实尾部调用执行的小检查已复现异常。此前仅用保存的1593长度kR与新kap单独执行布局检查，不能代表该可达控制流；该初始推断已撤回。

因此，确定影响是**当前增行固定C严格重放无法完成**，不是已证存在错误负上界或误接受。当前端点Phase I及普通 `sampling.analytic_joint_audit`使用正确n，不走此入口。仍应在抽行之前统一布局，避免让错误计算持续到末端才失败。

修正应先在 `analytic_fiber_audit` 用 `prepared_sampling` 构造所有kap；并在 `joint_audit`入口核对 `H.shape==(2*len(kap)+11,p)` 与current的物理行数量，再抽行。

## 2. 中优先级：固定振幅解析入口漏检端点声明

位置：`basis.audit_cardinal`约241–246行、`certificates.analytic_fiber_audit`约217–235行。

`source-audit`根据旧report的 `fixed_amplitude` 自动转入 `analytic_fiber_audit`；该入口没有调用 `validate_endpoint_model`。因此当输入current元数据是阶数2、用户省略 `--asymptotic-zeros` 时，CLI参数仍记录 `infinity=free`，而 `joint_audit`内部依据current元数据默默按端点子集给界。这个界即使数学上对端点子集合法，也不能作为原free模型的界解释。该入口结果还没有完整复制 `current_model`/端点阶数来消除歧义。

此问题属于真实入口的声明缺口：原生网格端点输入可到达该结果；当前增行输入会先因第1项的最终布局保护失败，修正布局后仍需关闭这个身份缺口。当前主Phase I入口已经调用validator，不受影响。应在固定振幅解析入口建立resolved prescription、调用同一个validator，并显式记录current模型及端点阶数；同时核对旧report的模型身份，不能仅核对路径和χ/B/epsilon。

## 3. 中优先级：非零阶数布尔比较不足以落实双零点承诺

位置：`endpoints.validate_endpoint_model`约42–47行。

当前检查只比较 `bool(ff_endpoint_order)` 与 `infinity=='zero'`，并在缺少resolved prescription时默认使用analytic。因此可构造 `asymptotic_zeros=True`、`infinity=zero`、current阶数1仍返回True的输入；本审阅用SimpleNamespace和阶数1元数据复现了这一行为。

当前真实 `endpoint_current` 是新CLI生成的阶数2，所以这不是当前运行中的模型切换，也不是已发生的生产错误。安全间隙应严格核对阶数等于 `2 if args.asymptotic_zeros else 0`，并从已解析身份或显式 `args.prescription` 获取函数族，不默认认作analytic。模型签名最好附带阶数/完成版本；目前 `infinity`与新的current preparation路径已能区分现有正常旧/新生产模型。

## 仿射问题与证书检查

`configure_endpoint_coordinates`的五个固定坐标对应T0与两通道各两个FF条件。T0枢轴选自自由振幅列，并避开x及第二投影坐标；另外四个枢轴是FF列。矩阵列变换同时作用于A/R/I/J/F/W和cost，C重新取A的χ行；逆变换记录在 `endpoint_inverse`。密度坐标本身不变，因此原L4障碍及其梯度/Hessian不需要另作密度换元。先构建可逆线性坐标，再把五个坐标锁定为(0,1,0,1,0)，表达相同的仿射切片。

`linear`从Newton方向删除这五列；固定x时额外删除x列，射线模式明确拒绝。Phase I保留单独tau变量，并在端点投影后重新选择能容纳当前Gram/FF/矩的tau。`PhaseOneProblem.raw`先去掉tau，再调用统一raw完成，维度一致。

数值Newton使用浮点K与长双精度端点矩阵，是精确模型的舍入数值表示；`completed_raw_point`把C0写为精确0，并以保存的前M−2个独立FF浮点数重新定义两个依赖变量。严格 `joint_audit` 与 `sampling.analytic_joint_audit`都会重新调用相同的 `complete_imf`，不会把末尾两个显示float作为独立精确参数。二者使用精确K。原生Newton与严格完成之间的舍入差由候选原式核验处理；数值中心日志不是该精确模型的独立数学中心证书。

端点对偶等式为 r·v=λ·t+(r_f−λE_f)·v_f，t=(1,0)，故 `endpoint_residual`返回的常数λ0应加到c0。对M=3、5、8，以384-bit Arb执行了双约束残差、一般r的完整恒等式及依赖显示值不变性检查，全部通过。特别令r=−3b+5d得到常数−3且剩余残差包络含0，确认 `c0 += constant` 的符号。`support_outer(... infinity_zero=True)`省去T0残差并要求C0精确零，符合同一切片。

## 严格陈述、输出身份与缓存

- `__init__.joint_result`从args保留 `infinity`，完整复制current元数据；`resolution_report`复用该结果。`sampling.analytic_joint_audit`成功分支也通过同一结果函数保存；其严格结果写出5个端点等式与依赖变量显示值说明。未发现当前正常链将zero改回free。
- 严格联合、解析联合核验都使用精确FF完成。`spectra.direct_profiles`的无限远零值依据已保存阶数定义；一般谱图使用float作显示，与本次授权一致。`analytic.audit_reconstruction`中仍有直接使用全部显示ImF的无限远公式，但该历史入口明确拒绝非PV坐标，正常analytic-cardinal端点点不能进入；不把这个不可达的旧分支列为当前端点错误。
- 当前CLI先将显式端点旗标解析为 `infinity=zero`，再选择续算μ。续算检查两套preparation路径、infinity、χ/B、section和目标；`center_joint_support`用保存reference_point重建同一坐标变换后读取z与传播值。当前current路径固定且阶数2，旧free缓存不能正常跨到新current路径。建议在checkpoint另存端点阶数/完成版本/枢轴以便更早拒绝身份错误，但本审阅未发现当前正常续算路径发生静默丢失。
- ROI/代表签名保留 `infinity`与current preparation路径；旧free与新zero不能通过正常身份比较。显式阶数纳入签名仍是合理加固，尤其若未来提供单零点/双零点两个变体。
- `joint_hull_candidate`仍是历史有限凸包/固定C初始化器，没有把五个端点等式送入其conic求解；正常端点主入口使用完整resolution Phase I。该旁路最终仍经过端点原式audit，不能据一个候选或conic状态声称端点联合可行。若未来面向用户开放端点fixed-amplitude求解，应明确补齐或拒绝此未实现的初始化约束。

以上没有重跑全套102项测试；已有测试数由root提供，本审阅只计入自己实际执行的小规模Arb等式检查、真实输入布局检查和源码路径审阅。

## 审阅输入指纹

审阅时间：2026-09-11T06:14:34.315004+00:00。原 `analytic_fiber_audit` 与 window_support_01快照的AST相同：True；抽行遗留不由本次FF代数变换新引入。

| 输入 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/endpoints.py` | `3c289a25507ef9c7087c8be904705cdb110dfe5d0d61be9aaffdcb7ccae6cf5a` |
| `src/smatrix_bootstrap/operators.py` | `2786c843c58d98619f90c346c463c676a10df806ab4b9bf65e1c231f8648e9d1` |
| `src/smatrix_bootstrap/linear.py` | `b7268ea3055de6f5a0296cce780cf78d9b11a97f0574edd85256232b5b6020a8` |
| `src/smatrix_bootstrap/run.py` | `2d5bad079703c6b041378ab734e91b5c1e3964371716a65cc356e7a23cb89175` |
| `src/smatrix_bootstrap/gauge.py` | `f4494bef9128d2e17753b67d65fb5b4dfbc171d0ef08a0f087e642307eb16017` |
| `src/smatrix_bootstrap/certificates.py` | `879d23b7b772909ede6b997926ae1e7b0dac61d9f0b1a61d878bb41e9c20c6fd` |
| `src/smatrix_bootstrap/sampling.py` | `5fd7b5b7e406cc6df1bbbdbb962773ef86a947d3857cd226fc5827e873314ba2` |
| `src/smatrix_bootstrap/__init__.py` | `9e65f9de80eaceb49fe41feb6d29c8335cafa24f7fa588f46d1b63cd66ff6cbb` |
| `src/smatrix_bootstrap/quotient.py` | `3ec58dc689af80cf12b564c74e40cd91af2189cda1789686694cc9691ef8c16d` |
| `src/smatrix_bootstrap/basis.py` | `5eda7c08a5cea679958b47980c9972ada51358f6352698d8db10db91bd5abec3` |
| `src/smatrix_bootstrap/imaginary.py` | `fbf4ac702058f5c14733b572c7f3091db9cb65209ca2638968c3578f00a9b488` |
| `src/smatrix_bootstrap/spectra.py` | `07489253c120fe925489a5c4560980a30128012ab61972c04aa1ee3adfd528a4` |
| `src/smatrix_bootstrap/analytic.py` | `f7b6ffb9bd6b520f00dec8653ff0864da6fe46eaa0ed5814bef7bd3c5c22a727` |
| `results/runs/physical_consistency_20260911/ENDPOINT_DESIGN_ZH.md` | `16226a73719bd1dd19509ed92e77782afe4d44d803eafb82e97273b3e6eb3879` |
| `results/runs/physical_consistency_20260911/PLAN_ZH.md` | `5436ed5ae8d89ffef6afea39768058ba519fdf109f9de0819fa948c7f184cf5c` |
| `results/runs/physical_consistency_20260911/endpoint_current/report.json` | `a5370e9439190aed21e5e45141deac42c73d5a88f7799c1a2a0a401a2cf7e61d` |
| `results/runs/physical_consistency_20260911/endpoint_current/current_data.npz` | `ff8d5809c68515899780eb4793c4d2bb42dbf95d62fdb805d3b60044d071eac9` |
| `results/runs/physical_consistency_20260911/window_prepare/report.json` | `8884fb8c0a17fc9b5173d33ab96b7b3b8726d2f8850518c5fde8fe4557d944aa` |
| `results/runs/physical_consistency_20260911/window_prepare/amplitude.npz` | `2ab1720aa2c8127dcebef42a14d395c57414e91e6660e26f00684eeb2565c528` |
| `results/runs/physical_consistency_20260911/window_support_01/source.json` | `b38c06c10c7e6ad0e724590e9aa6cfe1568b5a193a9cc97bf128667de71adceb` |
