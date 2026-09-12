# Arb 局部下降核验：限定只读审阅

2026-09-12。审阅 [LOCAL_DESCENT_PROOF_ZH.md](LOCAL_DESCENT_PROOF_ZH.md)、merit.py、linear.py 的新增调用及相关测试；对照 [修改前清单](merit_before/manifest.json)。未运行优化、求值或状态重建，未修改源码或生产输入。

**在声明的普通、无额外渐近障碍的缓存局部模型范围内，未发现阻断性错误。** 密度区间求根、Gram 数值约定和严格下降判定与推导相符；该结论不将数值中心提升为原 H 的精确驻点或物理唯一性结论。

## 密度及区间 Newton

`_binary` 使用 longdouble 的整数比，保留输入的实际二进制数值；后续除以 B、求根与求和由 Arb 包络。对 r=1−Σρ⁴>0，F(0)=−r，且 F′∈[1+N/2,1+N]，故初始区间 [r/(1+N),r/(1+N/2)] 包含真根。

每轮先在固定 midpoint 包住 F，再在当前区间包住 F′，与 interval-Newton 像取交集。相关量的区间依赖只使包络更宽，不破坏包含性；16 轮提前结束仍保留根包络。非正／无法判定的密度严格性不被接受。

梯度 4ρv/δ、对角 Hessian 2/d+8ρ²/(δ+2d²)，以及 rank²/den 项均与隐式微分／Schur 消元一致。密度部分极小化的第三导数中，v″ 项由 `(H w)_v=0` 消去；不存在漏加这一项的问题。

## Gram、其它锥与尺度

Gram 解包采用与 `imaginary.symmetric` 相同的 `sqrt(longdouble(2))`，再将该**数值常数**转成精确二进制有理数，避免混用理想 √2 与保存的 svec 约定。Arb 检查三个顺序主子式为正，随后以 T=G⁻¹dG 计算 −tr(T) 和 tr(T²)。不需要假设数值 svec 在精确实数算术中仍严格等距。

散射、χ、FF 的 first/curved 项对应局部正二次裕量，moment 项对应 1−r²；χ 权重只乘一次，并要求 w≥1。cost·dz/μ 的 μ 只除一次。总曲率须被包为非负，t>0 且区间半径 t√h 严格小于 1；充分下降式的下端点严格为正才产生 `certified_descent=True`。

## Caller 与证据范围

- 新 Arb 路径是最后一层 fallback：只在普通 joint、原线性残差精度合格、decrement<1e−6，且现有下降读数未获接受时调用。Phase I 被 caller 排除；helper 另拒绝 Phase-I moment_slack 和额外 asymptotic 项。
- 正的原函数差及原有浮点 self-concordant fallback 均保留，所以不能将所有接受步都写成“已获 Arb 局部下降证书”。实际方法由日志区分。
- 证明对象使用缓存盘／χ／FF 的正裕量和二进制局部方向，定义局部一致的平移二次锥；Gram 和 density 也按该明确局部数据解释。它不包住缓存到原 H、坐标回传或最终浮点存储的全部误差。
- `certified_descent` 不设置 `converged`。普通中心仍须通过 1e−14 数值门槛与 fresh QR 检查；完整 C 的原式 primal/dual 审计仍独立。
- 定量区间证据应引用日志保存的 `proof.decrease_lower` 等区间记录。caller 的 `barrier_merit_decrease` 是该数值向 float 的便利表示，不应单独当作向外舍入的严格下界。

## 现有测试与来源

已阅读独立 log1p 真值与反方向拒绝检查，以及完整联合模型有限变分对密度、Gram、FF、moment、χ 的检查。14 项针对检查通过由主线报告；本审阅未独立重跑测试，也不提前声明正在执行的全套检查结果。

当前 Python 文件均通过语法解析；merit.py 为 92 行／5595 bytes，linear.py 为 229 行／18122 bytes，测试文件 350 行。

| 审阅对象 | SHA-256 |
|---|---|
| merit.py | `dfccd05c8852c07b114e7fed7a1aa95245c2bb85da797367383bc24c2b3f4cb3` |
| linear.py | `5881abe58f1328fab3589f75320e25b0ccca1add5acee3b55809d7ae38d66a96` |
| test_mainline.py | `9bb4c09fc2f3c3e48c44f7ecc7e4805924b49df0dc22790682632c5a618e8159` |
| LOCAL_DESCENT_PROOF_ZH.md | `cb4f61bbec57c94ae296d5810d19896bd04f9f3ed84281d1ab89e20024309e06` |
