# 全谱锥模型与SCS接口只读复核

结论：未发现本版变量保留、暖启动变换、SCS行排列/返回逆排列或结果接受逻辑的阻断错误。conic现保留全部2M个谱变量和全部2M个Gram，前一版删除自由谱后闭边界range条件丢失的问题不再适用于该完整锥模型。Newton默认变量集合保持不变。

本复核没有运行新求解、测试或修改src/tests。主线报告的单个SCS控制和pack检查通过仅作为已有验证记录；106项完整测试在本次请求时仍在运行，本报告不将其写成我重新完成的结果。

## 1. 全变量保留及默认兼容

`JointProblem`新增参数放在签名末尾且默认False。默认 `retained_spectra=bounded`，active和keep与原Newton相同。conic显式传True后，retained_spectra为全部2M个True，active包含连续的p+4M原变量，keep包含全部2M个索引。J/g实际循环和PSD cone数量都以这个keep生成，并非只改元数据。

bounded仍保存真实的矩支持掩码，没有为自由R添加新的上界、目标惩罚或矩约束。端点切片只消去既有五个固定坐标；全部剩余C/ImF/R方向保留。L4辅助w与正列缩放沿用前次已核对的同一提升。实际M50的全部100个Gram现在直接进入原始锥问题，所以闭边界的range相容条件由完整PSD约束自身施加。

最终audit仍可按原规则增加无上界高R以修复数值舍入，但不会免检任何Gram。浮点组装/SCS容差仍不等于严格原式或解析源误差认证。

## 2. 暖启动是建议点，不是额外约束

`conic_warm_point`按正确逆序恢复数值坐标：从完整raw active值开始，振幅除ds，解P.inverse取得absorptive坐标，再解endpoint_inverse，最后把声明的固定端点坐标设回目标值。这正是P.raw所用逆映射的反向；显示FF依赖值的微小数值差或端点投影只影响建议点，不改变Q的可行域。

密度w先取rho²，再在其二范数小于1时增加统一正余量。若原density已在L4内域，三角不等式保证这个增量仍保留全局w球内域，并令局部SOC严格；若warm本身不可行，代码没有把它认证为可行解。最终Q固定offset、scale、columns的还原顺序正确，辅助w不会被误当作原C或谱。

SCS调用前显式将x转成连续float64，再构造原排列的s=b−Ax，按SCS行顺序传入，并给y=0。这使暖数据满足线性关系；s未必在K中，而SCS仍需实际求解。x没有被固定，没有缩小可行域。Clarabel分支不使用warm，warm_start_used=False与实际一致。

前个控制中“仅x暖启动＋longdouble”和本版“显式float64＋完整x/y/s”两项同时变化；现有一次通过只能支持组合设置工作，不能归因之前残差8单独来自哪一项，也不能预告真实M50会收敛。

## 3. SCS排列及dual/slack还原

`scs_layout`首先放8行FESR非负锥，随后按照blocks顺序放各SOC，最后放全部PSD块。scattering宽4、density_squares/FF/tail宽3，χ和全局density_norm各自作为一个完整SOC；零长度FF块被跳过。由compile产生的blocks本身是全行不重叠分区，这个构造覆盖全部行；既有测试另验证它是全行置换。

原Gram排列为(11,√2·12,22,√2·13,√2·23,33)，SCS下三角按列排列为(11,√2·21,√2·31,22,√2·32,33)。对实对称3×3矩阵，置换[0,1,3,2,4,5]正确，不需要再次乘√2。

定义order使新行i来自旧行order[i]。输入A[order]、b[order]、initial_slack[order]一致。返回时 `dual[order]=sol['y']`、`slack[order]=sol['s']`是正确的全局逆排列；不能改成直接再索引order。x没有行维度重排，按原列映射恢复即可。保存的conic_candidate中的dual/slack因此属于原conic_A/b行序；它们没有被冒充为原式支持证书，原式审计使用独立零对偶。

## 4. current、reference与状态

派发前仍核对H/current路径、分辨率、UV配置、endpoint/tail旗标和完整C/NPZ一致性。读取旧Newton缓存时只从同operator路径的reference_point构造新坐标/congruence，没有把旧z/tau/部分Gram缓存当作全谱解。

有一个需要明确的**条件性兼容限制**：保留全谱后，JointProblem正对角检查也覆盖原来被省略的自由R。旧Newton reference可能在这些位置为0，却仍满足当时的构造要求；若优先读取这种旧reference，新构造会在编译前明确报错。未检查出当前指定生产输入存在这个问题，也不将其列为当前阻断。它是参考congruence的构造条件，不能解释成物理不可行；完整正对角的reference或等价正congruence可用于后续处理，无需改变物理集合。

SCS与Clarabel都只是针对同一个已保存Q的算法。native status、残差、版本和暖启动标志被保存；非有限x明确失败，有限x必须经P.raw和完整joint_audit。solved/inaccurate/infeasible文本均没有替代原式接受。

method维持phase_I=True、support_search=False、representative_center_converged=False、infeasibility_certified=False。resolution_report因此强制零目标方向，不产生support最优性或代表selection；初始化器再次显式关闭support_optimality。新的全锥缓存身份仍不伪装为Newton中心。原式通过后还需主线既定的解析源核验。

本结论限定本版代码和已给定的锥接口约定，不包含下一次真实求解的性能或物理成功判断。

## 复核输入指纹

时间：2026-09-11T09:23:27.807968+00:00。

| 文件 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/operators.py` | `74cb905a1e9e881e3f3e5f7843d44e7e096e8c57e672f8f7abdbae336b78fc4e` |
| `src/smatrix_bootstrap/conic.py` | `90e83ce342f67dbd35b845204dd7f1f7927bcc50bf920835dc16284902f64944` |
| `src/smatrix_bootstrap/endpoints.py` | `6b4dacd438e750ca61f285744cc0a152d16e3090581e22c77962f23b867005db` |
| `src/smatrix_bootstrap/gauge.py` | `4fbfb5b647cafd1e3feea523a89c31d4051c469a833eba0a20d0b4b819f97829` |
| `src/smatrix_bootstrap/run.py` | `b8c26f63dd637a36c1b97fdda612a19f1a3bca445501c8057e4dc4fd4e3f83a1` |
| `src/smatrix_bootstrap/certificates.py` | `726151e3225eda13f82260d589b4a503c1ebd3afdbdddd05831e440a7efc66a4` |
| `src/smatrix_bootstrap/__init__.py` | `9e65f9de80eaceb49fe41feb6d29c8335cafa24f7fa588f46d1b63cd66ff6cbb` |
| `tests/test_mainline.py` | `ca86f1daa9ddd9d38f158ec63c66fb1cb1038b06f4b444d9b8b79934a8054050` |
