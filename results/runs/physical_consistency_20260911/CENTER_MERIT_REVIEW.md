# 近中心下降值保护：限定只读审阅

2026-09-11。仅审阅[完整推导](../../evidence/LOCAL_BARRIER_DESCENT_BOUND_20260911.md)、`linear.py` 的 helper/caller 及 [修改前快照清单](center_merit_before/manifest.json)对应测试。未改源码、状态或物理输入，未运行新优化、求值或状态重建。

**推导及普通 joint 调用范围一致，未发现阻断当前 Watson polish 的问题。** 发现的 Phase I 调用范围问题保留为明确待办；当前 `polish_03` 继续使用其已保存的源码，不中断或改写该 run。

## 推导核对

- 各障碍能写成仿射正定矩阵的 −logdet，或其和；散射／FF／手征球的箭形矩阵写法正确。静态合同变换只添常数。w≥1 的正权与各项相加保持常数 2 的自协调不等式；仿射换元、等式切片和线性目标不破坏该结论。
- 密度部分确实是现有 SOC 提升障碍的严格部分极小化。令 γ(t)=(q+th,v*(q+th))、p=γ′，由 f_v=0 得 (Hp)_v=0；于是 b''=pᵀHp，求导产生的 2γ''ᵀHp 为零，b'''=D³f[p,p,p]。v* 的二阶导数没有被遗漏。该结论针对精确极小化的数学障碍，不把有限次浮点求根误差自动包含在内。
- 沿线二阶导上界积分给出 ta−ω*(tb)，再用 ω*(x)≤x²/[2(1−x)] 得到实现所用保守界。条件 t>0、a>0、b²≥0、tb<1 必不可少，helper 均检查；b=0 时公式自然返回 ta。

## Caller 核对

Φ=barrier−cost·z/μ。caller 传入 `a=model_decrement`、`b²=model_curvature`、`t=step`；a 已含 μ，Hessian 中线性目标无贡献，没有重复缩放。future 来自既有严格可行线搜索；普通 joint 的中心分支在此前已经结束，fallback 不会授予中心资格。

fallback 仅在方向通过原 1e−6 线性残差门槛、decrement<1e−6 且原始函数差落入消去尺度估计时触发。普通 joint 在局部界≤0、非法半径或非有限界时仍失败；超出消去尺度的明显负下降不会被覆盖。普通中心 1e−14、fresh-QR 重验、支持 gap 1e−4 和物理约束均未更改。

**待办：本轮结束后在 fallback 条件加 `not phase_one`，保留 Phase I 原有行为。** 当前 helper 也可被 Phase I 调用，但 Phase I 的辅助 centered 分支不在此处之前结束，后续 `not centered` 条件可能跳过非正局部界的失败门，例如零方向导致 step=0／局部界 −inf。因此“非正局部界必失败”目前只对本次普通 joint 路径成立。此问题不影响正在运行的 Watson polish，未据此改动运行中源码。

## 数值证据的范围

64εmach 乘函数值量级是选择何时避开相减的估计，不是总浮点误差上界；它不覆盖梯度／曲率计算、缓存误差或密度 δ 最多 12 次 Newton 更新的全部求根误差。浮点求出的 `model_decrease_lower` 是数值保护，不能称为 Arb 区间下降证书。最后原式与解析接受仍各自独立执行。

已阅读大常数平移的一维 −log 测试：高精度 Arb 的真实下降独立于 helper，并严格大于其正值；大函数值直接相减确已舍入为零。这验证了该例子的稳定求值用途，不构成当前高维点的区间证明。本审阅未独立重跑 106 项测试，不更新主线测试状态。

## 源码与证明身份

修改前两个文件的哈希匹配快照清单。审阅时 `linear.py` 227 行、测试文件 350 行，均通过语法解析。

| 对象 | SHA-256 |
|---|---|
| 修改前 linear.py | `d6138aa0b6770bbb55e1e7e5c83c5425ddd6228512fe8791c9f13cdcca7a003e` |
| 本次 linear.py | `dd45e60ae4bd3138a041f50a287b3a954ff9ee7da8b7c6cf5e744329466e9142` |
| 本次 test_mainline.py | `031af4cbc56ef682174ef19dbacdeda2b31f448ffd910365b394639020a29005` |
| LOCAL_BARRIER_DESCENT_BOUND_20260911.md | `8047e2fe5788090e8504588c66c718a234f85c96b9330c691c0b8d8d986c76b1` |

主线关闭记录：polish_03运行结束后已在fallback条件加入 `not phase_one`，保留Phase I原流程；106项全套再次通过，见 `center_merit_scoped_pytest.log`。普通Watson路径不受此范围修补影响。
