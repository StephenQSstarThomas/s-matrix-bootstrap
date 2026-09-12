# Free-T0 有界 IR 坐标：限定只读审阅

2026-09-11。对照 [修改前快照](PV_bounded_coordinates_before/manifest.json)核对 factory、scattering caller、支持／宽度及 center raw 恢复。未修改源码、当前 worker 或原输入，未运行优化、求值或新数值检查。

**未发现 κ 重复缩放、排列、dual/width 或原 C 恢复的阻断问题。** 107 项测试通过是主线提供的结果；本审阅只阅读相关测试并核对源码和语法。

## 坐标和 κ

令 h=κf。factory 的前 2M 个坐标依次是 M 个原生 S0 的 Im h 和 M 个原生 S2 的 Im h；第 2M+1 个坐标是原生 S0 的 Re h，锚点索引 `(M//2)*3*L`。由散射盘，前者在 [0,2]，后者在 [−1,1]。

absorptive_change 后先把 T0 所在坐标替换为所选 Re h，再排列为 `[两组 Im h, Re h, density]`。虚部坐标乘 κ，同时 A 列和逆变换除以同一 κ。实部坐标在替换时已经包含 κ，不再乘第二次。因而新 A 仍输出 f：原生虚部轴是 e_j/κ_j，实部锚是 e_(2M)/κ_anchor；scattering 的 raw=False 分支再乘 κ，才得到 h。没有将 κ²f 误作盘变量。

T0 没有删除或置零，通过完整逆变换恢复；密度坐标及实际密度打包不变。原生 S0 实部的 T0 系数 2.5 与正 κ 检查保证实部 pivot 非零。小 H 测试覆盖所有列的 H 映射、原生轴、非零 T0 往返和密度行保持；浮点恒等式按其声明精度检查，最终仍需原 H 审计。

## 支持、宽度及 caller

- free=2M+1 且 raw=False 只在容器具有 `bounded_constant=True` 时放行；正常 caller 的该标志由本 factory 提供。旧 zero-T0 factory 仍为 free=2M。
- 盒支持 `2*sum(max(r[:2M],0))+abs(r[2M])` 正确处理 [0,2] 与 [−1,1] 的不同中心；不能把最后一项写成前两组相同的正部。
- 差值宽度 `2*sum(abs(r[:free]))+2*B*norm(rho_residual,4/3)` 对所有 free 坐标使用长度为 2 的区间，正确。density 未被变换，原 L4 对偶范数仍适用。
- 盘与 χ 的对偶仍按约束行编号；无需再把这些行乘子乘一次坐标逆变换。最后 `audit_ir` 使用原 H 和完整 raw C 重新审计，内部盒支持估计不替代该接受步骤。
- P 作为现有 `asymptotic` 参数容器传入，但没有 `asymptotic_R`；augment/slope 会保持原结果／返回零，不会给 PV 暗加高能尾约束。
- 固定 x 继续以变换后的 A[-2] 定义同一线性截面；后续消元、方向恢复和原坐标审计的接口一致。

须记录一个既有数值联动：raw=False 使用 `mu_factor=.5`，原 raw=True 使用 `.1`。新路径因此同时采用减半续算；它不改变物理集合，但不能把未来全部数值改善只归因于坐标换元这一项。

## 中心与失败快照

capture 与最终交付都使用同一个 P.raw。保存的 `coefficients` 是完整原 C，`numerical_coordinates='absorptive-free-T0'` 明确区别于零 T0 坐标。`center_identity` 优先使用保存的完整 coefficients，不会把新 z 按旧 unsubtracted 规则解释；缺少可识别 raw 记录时不会授予中心资格。

此前失败 NPZ 字段备注已关闭：现在具有 initialization_only=True、center_converged=False、M、prescription 和 numerical_coordinates；相邻 JSON 也明确仅作初始化。失败文件不进入既有中心选择名单。

## 文件身份

三个旧文件均匹配快照清单。当前三文件通过语法解析，两个核心文件为 234／306 行、18307／23039 bytes，测试文件 350 行。

| 当前文件 | SHA-256 |
|---|---|
| ir.py | `bd6fb3d54b523cbf266da65f677e9364c1393aaae0f29f114bd315e60c6e5c17` |
| scattering.py | `29d6d8355d1f84632656df7d7de9968e7bb5458748c3d3646805861214fa54eb` |
| test_linear.py | `b3eabdae425084bb155b4d366117e4771737b2fbd843012432684dd63b518890` |
