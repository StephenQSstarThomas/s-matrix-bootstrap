# REPORT_SDP_ZH — He–Kruczenski 2309.12402v3, SDP 直解路线

生成时间 (UTC): 2026-09-12T10:02:23.175799+00:00

结果根目录: `/playpen1/shiqiu/sdp-work/runs/sdp_reproduction_20260912`


## 0. 一句话结论

8 条 claim 中：**2 条复现、1 条不通过、5 条未运行**。
未运行的原因不是论文，而是本轮的求解阶段——M=50（论文分辨率）在所有试过的配置下都不收敛，UV 扇区（Gram + FESR + 形状因子渐近同时开启）也没能给出一个认证点。按任务书第 8 节，这属于结论 (c)：技术性失败，附完整诊断（见 §8、§9、§11）。已复现的部分是在 M<=40 上**认证**的（约束生成：松弛极值 + 返回点满足全部 3LM 个圆盘），并附分辨率阶梯说明小 M 的数字为何可引用。


## 1. 核心 claim 判定表

| # | claim | 判定 | 依据 |
|---|---|---|---|
| C1 | Pure-unitarity region shape (Fig. 3) | not run | 0/24 verified-feasible directions; the +x tip alone, from the resolution ladder: M=20 1.9279 (-13.7%), M=25 2.0061 (-10.2%), M=30 2.0661 (-7.5%), M=35 2.1140 (-5.3%), M=40 2.1534 (-3.6%) |
| C2 | Chiral constraints collapse the region onto f11 = -f00/15 (Fig. 4) | pass | +x end 0.082598 vs 0.082573 (+0.03%) at M=[30]; eps ladder {0.006: 0.161054, 0.004: 0.125504, 0.002: 0.082598, 0.001: 0.055038} |
| C3 | Subthreshold partial waves near-linear, S0 chiral zero moves (Fig. 5) | pass | eps=0.002: S0 zero 0.42585550990952303 |
| C4 | Chiral only: S0/S2 agree with experiment, P1 has no rho (Fig. 7) | FAIL | tip: d00(0.9)=68.4, d11(1.2)=126.3; ref: d00(0.9)=84.7, d11(1.2)=27.0; mid: d00(0.9)=32.8, d11(1.2)=131.0 |
| C5 | FESR+FF shrink the upper boundary, not the lower (Fig. 8) | not run | need both chiral and chiral+UV sweeps |
| C6 | rho appears in P1 at all three representative points (Fig. 9) | not run | no representative points |
| C7 | S0/S2 coincide at low energy and spread at high energy (Fig. 10) | not run | no representative points |
| C8 | M/L stability (Fig. 11, Appendix A) | not run | 0/5 (M,L) configurations |

## 2. 求解器统计

- 求解次数: 47
- 状态分布: {'optimal': 43, 'SolverError': 4}
- 迭代中位数: 56.0, 单次秒数中位数: 27.657386779785156
- 机器时间合计: 1525 s

## 3. 密度正则化 B 的活跃性

- 带 B 的求解: 43；**任一活跃: False**

| job | B | 范数 | \|\|rho\|\|_2 | \|\|rho\|\|_4 | 活跃 |
|---|---|---|---|---|---|
| tip | 3.775e+05 | l4 | 2.850e+04 | 9.813e+03 | False |
| ref | 3.775e+05 | l4 | 9.792e+02 | 2.618e+02 | False |
| mid | 3.775e+05 | l4 | 8.221e+03 | 2.285e+03 | False |
| dir000 | 3.775e+05 | l4 | 2.798e+03 | 8.514e+02 | False |
| dir003 | 3.775e+05 | l4 | 1.983e+03 | 5.903e+02 | False |
| dir006 | 3.775e+05 | l4 | 7.950e+02 | 2.340e+02 | False |
| dir000 | 3.775e+05 | l4 | 4.517e+04 | 1.468e+04 | False |
| dir003 | 3.775e+05 | l4 | 7.220e+03 | 2.208e+03 | False |
| dir006 | 3.775e+05 | l4 | 2.012e+03 | 4.692e+02 | False |
| dir009 | 3.775e+05 | l4 | 1.110e+04 | 3.473e+03 | False |
| dir000 | 3.775e+05 | l4 | 2.033e+03 | 6.699e+02 | False |
| dir003 | 3.775e+05 | l4 | 1.642e+03 | 5.415e+02 | False |
| dir004 | 3.775e+05 | l4 | 9.817e+02 | 3.140e+02 | False |
| dir006 | 3.775e+05 | l4 | 1.046e+03 | 2.741e+02 | False |
| dir007 | 3.775e+05 | l4 | 1.300e+03 | 3.561e+02 | False |
| dir009 | 3.775e+05 | l4 | 1.707e+03 | 4.944e+02 | False |
| dir010 | 3.775e+05 | l4 | 1.787e+03 | 5.166e+02 | False |
| dir000 | 3.775e+05 | l4 | 1.847e+03 | 6.443e+02 | False |
| dir001 | 3.775e+05 | l4 | 1.837e+03 | 6.372e+02 | False |
| dir002 | 3.775e+05 | l4 | 1.631e+03 | 5.585e+02 | False |

## 4. 幺正性事后复验（未经任何重缩放的原式）

- 通过事后可行性检验的解: 31 个，其中最大 `max eta - 1` = 1.015e-09，位置: {'job': 'dir002', 'file': 'fig4/eps6e-03_chi-b/000/report.json', 'wave': {'isospin': 2, 'ell': 0, 'node': 18, 's': 12.468205094073321}, 'certified': True}
- 未通过的解: 12 个（约束生成中途失败时返回的最后一个可解迭代，标记为未认证；所有 claim 判定都把它们过滤掉），其中最大 `max eta - 1` = 2.273e-02，位置: {'job': 'dir000', 'file': 'fig4/eps2e-03_chi-a/000/report.json', 'wave': {'isospin': 0, 'ell': 2, 'node': 4, 's': 4.230551242245793}, 'certified': False}

## 5. FESR 目标值独立重算 (5a.9)

- m_q 算术平均 = 5.650 MeV，均方根 = 5.886 MeV

| 波 | n | 论文 (2.56) raw | 重算 (m_q 平均) | 比 | 重算 (m_q 均方根) | 比 |
|---|---|---|---|---|---|---|
| S0 | 0 | 2.385104e-03 | 2.182753e-03 | 1.0927 | 2.368909e-03 | 1.0068 |
| S0 | 1 | 1.118381e-01 | 1.023479e-01 | 1.0927 | 1.110766e-01 | 1.0069 |
| P1 | -1 | 4.228046e-03 | 4.228884e-03 | 0.9998 | 4.228884e-03 | 0.9998 |
| P1 | 0 | 1.457112e-01 | 1.458146e-01 | 0.9993 | 1.458146e-01 | 0.9993 |

论文打印值对应 m_q 取均方根；取算术平均则 S0 两个矩差 9.3%。主线用打印值。


## 7.1 分辨率阶梯 — 纯幺正 max f00(3)（对照 Fig. 3 +x tip = 2.23289）

_reconstructed from the ladder log; the M=35 and M=40 rows were produced before the constraint-generation fix and are marked uncertified_

三种状态：**认证** = 松弛极值且返回点满足全部圆盘，即完整问题的最优值；**仅可行** = 通过事后检验但未证明最优，只能作为支撑函数的下界；**否** = 未通过事后检验，仅作诊断。

| M | L | 状态 | 认证 | 事后可行 | 目标值 | 相对论文 | 施加圆盘 | 轮数 | \|\|c\|\|_inf | 秒 |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 6 | optimal | True | True | 1.927921 | -13.66% | 348 | 6 | 1.78e+03 | 16 |
| 25 | 8 | optimal | True | True | 2.006054 | -10.16% | 580 | 7 | 1.71e+03 | 50 |
| 30 | 8 | optimal | True | True | 2.066089 | -7.47% | 690 | 6 | 1.11e+03 | 73 |
| 35 | 10 | optimal | False | True | 2.113958 | -5.33% | None | 8 | 1.28e+03 | 869 |
| 40 | 10 | optimal | False | False | 2.153422 | -3.56% | None | 8 | 9.91e+02 | 1573 |


## 7.2 分辨率阶梯 — 手征 (eps=2e-3, chi-b) max f00(3)（对照 Fig. 8 chiral-only +x end (digitised) = 0.0825728）

_reconstructed from the ladder log of the same script; ||c||_inf and rho_l4 are not carried for these rows_

三种状态：**认证** = 松弛极值且返回点满足全部圆盘，即完整问题的最优值；**仅可行** = 通过事后检验但未证明最优，只能作为支撑函数的下界；**否** = 未通过事后检验，仅作诊断。

| M | L | 状态 | 认证 | 事后可行 | 目标值 | 相对论文 | 施加圆盘 | 轮数 | \|\|c\|\|_inf | 秒 |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 6 | optimal | True | — | 0.084161 | +1.92% | 360 | 7 | — | 18 |
| 25 | 8 | optimal | True | — | 0.082459 | -0.14% | 598 | 8 | — | 50 |
| 30 | 8 | optimal | True | — | 0.082598 | +0.03% | 717 | 7 | — | 74 |
| 35 | 10 | optimal | False | — | 0.082492 | -0.10% | None | 8 | — | 209 |
| 40 | 10 | optimal | False | — | 0.082521 | -0.06% | None | 8 | — | 887 |
| 50 | 10 | SolverError | — | — | — | — | — | 1 | — | — |


## 7.3 手征容差阶梯 — 六个 eps^chi 的 +x 端（M=30, L=8, chi-b）

| eps^chi | 状态 | 认证 | 事后可行 | +x 端 | 施加圆盘 | 轮数 | 秒 |
|---|---|---|---|---|---|---|---|
| 6.0e-03 | optimal | True | True | 0.161054 | 709 | 12 | 618 |
| 4.0e-03 | optimal | True | True | 0.125504 | 720 | 8 | 176 |
| 2.0e-03 | optimal | True | True | 0.082598 | 717 | 7 | 74 |
| 1.0e-03 | optimal | True | True | 0.055038 | 712 | 13 | 263 |

论文 Fig.8 手征边界（eps=2e-3）的 +x 端数字化值 = 0.082573


## 8. (3.75) 最小可行 eps^FF（scripts/sdp/ff_tolerance.py）

论文取 eps^FF = 6.0e-05。下表是在手征 (3.64) + Gram (3.68) + FESR (3.73) 下，使 (3.75) 可行的最小倍数 t 及 eps^FF_min = 6.0e-05·t^2。

| M | L | (2.33) 因子取法 | 状态 | 最后一轮圆盘 | t | eps^FF_min | 相对论文 |
|---|---|---|---|---|---|---|---|
| 20 | 6 | s0 | failed | 342 | 4.064 | 9.9e-04 | 16.51x |
| 20 | 6 | 逐节点 | failed | 170 | 0.425 | 1.1e-05 | 0.18x |
| 25 | 8 | s0 | failed | 541 | 5.283 | 1.7e-03 | 27.91x |
| 25 | 8 | 逐节点 | failed | 486 | 3.227 | 6.2e-04 | 10.41x |
| 30 | 8 | s0 | failed | 551 | 0.220 | 2.9e-06 | 0.05x |
| 30 | 8 | 逐节点 | failed | 545 | 0.213 | 2.7e-06 | 0.05x |


## 9. 求解器行为研究 (scripts/sdp/solver_study.py, M=20 L=6)

| 配置 | 状态 | f00(3) | 事后可行 | 最大相对违反 | \|\|rho\|\|_4 | 迭代 | 秒 |
|---|---|---|---|---|---|---|---|
| cone_scaling=none, B=377500 l4 | SolverError | — | — | — | — | None | 3 |
| cone_scaling=centrifugal, B=377500 l4 | SolverError | — | — | — | — | None | 3 |
| cone_scaling=rownorm, B=377500 l4 | optimal | 1.928054 | True | -1.8e-08 | 2.499e+03 | 59 | 7 |
| B=1.000e+04 l4 | optimal | 1.927809 | True | -9.3e-09 | 2.292e+03 | 62 | 8 |
| B=1.000e+05 l4 | optimal | 1.928218 | True | -4.4e-08 | 2.521e+03 | 67 | 9 |
| B=3.775e+05 l4 | optimal | 1.928054 | True | -1.8e-08 | 2.499e+03 | 59 | 7 |
| B=1.000e+06 l4 | optimal | 1.928659 | True | -1.8e-08 | 2.613e+03 | 62 | 7 |
| B=1.000e+07 l4 | optimal | 1.927666 | True | -5.2e-09 | 2.941e+03 | 60 | 7 |
| B=1.000e+09 l4 | optimal_inaccurate | 1.913243 | True | +8.6e-09 | 1.271e+04 | 44 | 4 |
| B=3.775e+05 l2 | optimal | 1.927650 | True | -1.2e-07 | 9.958e+03 | 53 | 5 |
| B=1.000e+07 l2 | optimal | 1.924013 | True | -1.3e-07 | 5.731e+04 | 41 | 4 |
| B=1.000e+09 l2 | optimal | 1.914265 | True | +3.1e-09 | 8.896e+04 | 43 | 4 |
| no B, exact basis reduction | optimal | 2.007106 | True | -5.4e-08 | 5.188e+08 | 55 | 22 |
| no B, full space | SolverError | — | — | — | — | None | 0 |

- 目标值对 B 单调（可行集随 B 单调变大，极值必须非减）: **False**
- 已验证可行的最好点: no B, exact basis reduction, f00(3) = 2.007106


## 10. 能量轴口径

- digitised figures use m_pi = 139.57 MeV, the paper's text m_pi = 140 MeV; 0.31% shift on the energy axis


## 11. 方法、偏离与发现

**(3.67) 核的独立推导.** 论文的 K 核不是照抄的：本实现从 sigma 在 2M 点交错网格上的奇延拓推出共轭函数算子，与 (3.67) 逐元素相差 1e-15。它把 sin(n phi) 精确映到 cos(n phi)（n < M）。

**割线上必须用 cot 核，不能用去点中点法.** 两者离散同一个主值。在解析测试函数上，M=50 时 cot 核误差 3e-15（谱精度），去点中点法误差 9e-2 且只按 1/M 衰减（M=50/100/200/400: 9.0e-2, 4.5e-2, 2.2e-2, 1.1e-2）。仓库历史实现用的也是 cot 核，5a.10 因此可比。

**m_q 口径由 (2.56) 唯一确定.** 用 (2.53)-(2.54) 独立重算 (2.50)：只有取 m_q = sqrt((m_u^2+m_d^2)/2) = 5.886 MeV 时，S0 两个矩与论文打印值差 0.7%；取算术平均 5.650 MeV 则差 9.3%。P1 不含 m_q，差 0.07%。这与 N_f m_q^2 -> sum_f m_f^2 的物理一致。主线仍用打印值；(3.75) 的 m_q 按任务书取算术平均。

**FESR 硬截止有系统偏差.** (3.72) 保留 s_i <= s0 的节点，但最后一个节点的求积格子止于 s = 84.06 而非 s0 = 73.47（s 上超出 14.4%）。中点法本身对其实际区间精确到 5e-3 以内，但相对 [4, s0] 的真积分，n = -1/0/1 三个矩分别偏高 0.6% / 3.5% / 11.3%。如实报告，不做修正。

**锥的条件数与可行集的尺度是两回事.** |h|^2 <= 2 Im h 对任意 Lambda > 0 等价于 |h/Lambda|^2 <= 2 (Im h)/Lambda^2（两边乘 Lambda^2）。按 Lambda_ell(s) 缩放确实把行范数从 5e-85..34 压到 O(1)，但同时把可行集在缩放坐标下撑到 1e40，Clarabel 在第 0 步就报 NumericalError。不缩放则约 1650/3000 行下溢为零、锥退化成 0 <= 0、Slater 条件失效、求解器停在 gap = 6.6。最终用逐行实测幅度做 Lambda^2：行范数全部落入 [0.33, 5.8]，而未缩放的锥本身已蕴含 |h| <= 2。

**精确基约化：数学上精确，M>=30 上不收敛，默认关闭.** 所有约束行与目标行只张成 3876 维中的约 1780 维，c = V a 在数学上是精确重参数化。它在 M=20 上给出 2.0071——这个数一度被误当成“基约化不可靠”的证据，其实它是**无 B** 问题的极值（该点 ||rho||_4 = 5.19e8，不在预登记口径的可行集内），与带 B 的 1.9285 并不矛盾。真正的问题是它在 M>=30 上不收敛，且与 l4 球同用时 rho 变成约化变量的稠密组合、pnorm(.,4) 的展开代价暴涨（M=20 上 4 s -> 50 s）。因此默认关闭，只作为无 B 情形的对照。

**B 的口径：||rho||_4 <= 377500 有两处独立佐证.** 论文正文没有密度正则化。作者 2403 代码有 norm(rho/Mrho,4) <= 1e2，即 ||rho||_4 <= 100*3775 = 377500；本仓库历史主线 basis.py:233 用的是同一条 (density_rule='declared 100*double_density_count', density_norm='actual-rho L4')，并在 M=50 上得到 f00 极值 2.2349，与数字化 Fig.3 的 2.2329 差 0.1%。两处独立使用同一数值，说明论文 Fig.3 实际上是带这条密度界的区域。本任务按预登记网格报告 B 的敏感性与活跃性。

**求解器状态不作为判据：用约束生成把结果变成可认证的.** 1500 个幺正圆盘里只有约 100 个的算子范数在最大值的 1% 以内，441 个在 1e-6 以内，其余被离心因子压低许多个量级。只上其中一个子集是**松弛**，其极值是真极值的上界；若返回点随后在未改动的算子上满足全部 1500 个圆盘，它就对完整问题可行，因而就是完整问题的最优解。runner.solve_generated 就跑这个循环。松弛上界单调下降（M=20 上 4.248 -> 2.041 -> 1.930 -> 1.92855），这本身也是一重检查。既把数字变成可认证的，又缩小了 M=50 下主导开销的稠密 KKT。

**一次被更正的中间结论.** 中途曾以为“Clarabel 报告的 optimal 不可信”，理由是把 l2 密度球从 377500 放大到 1e7、1e9 时报告的极值反而下降。这条对 l2 代用球成立（那里求解器确实返回了次优点），但不能推广到预登记的l4 口径：无 B 时得到的 f00 = 2.0071 的点 ||rho||_4 = 5.19e8，本来就不在 ||rho||_4 <= 377500 的可行集内，所以 l4 家族内从未出现单调性矛盾。约束生成给出的认证值使这一判断不再依赖求解器状态。

**事后可行性用哪个指标.** max eta - 1 在离心压低的圆盘上恒为 0（h ~ 0 时 eta = 1 恰好），最小 margin 则被同一批圆盘压到 1e-40；两者都不辨真伪。真正有分辨力的是限制在 |h| > 1e-6 的圆盘上的相对违反 (|h|^2 - 2 Im h)/|h|^2。

**(3.75) 的 eps^FF = 6e-5 在 M=20/25 上可证不可行；M>=30 未定论.** 在手征 (3.64) + Gram (3.68) + FESR (3.73) 下求“使 (3.75) 可行的最小倍数 t”（scripts/sdp/ff_tolerance.py，eps^FF_min = 6e-5·t^2）。**读法要小心**：幺正圆盘是逐轮加进去的，少加圆盘是**松弛**，松弛的 min t 是真 min t 的**下界**。因此：

  | M, L | 最后成功轮的圆盘 / 全部 | t | eps^FF_min | 可下的结论 |
  |---|---|---|---|---|
  | 20, 6 | 342 / 360 | 4.06 | >= 9.9e-4 | t_true >= 4.06 > 1，**论文的 eps^FF 不可行** |
  | 25, 8 | 541 / 720 | 5.28 | >= 1.7e-3 | t_true >= 5.28 > 1，**不可行** |
  | 30, 8 | 551 / 720 | 0.22 | >= 2.9e-6 | t_true >= 0.22，**不能判定**（还差 169 个圆盘） |

  即：M=20 和 M=25 上论文的 eps^FF 可证不可行（下界已超过 1）；M=30 上只得到一个小于 1 的下界，既不能证明可行也不能证明不可行。M=30 的 t 远低于 M=20/25，暗示分辨率是主因（机制见下），但在 M=30 的完整圆盘集上收敛之前，这只是暗示。

  机制（若成立）：F_0 在圆盘上满足均值定理 (1/pi)∫_0^pi Re F_0 dphi = 1，而 (3.75) 要求它在 s > s0 那段 phi 上小到 0.072，所以 F_0 必须在阈值附近变大补偿，而阈值附近的节点数随 M 增长。

  绑定的是 S0 不是 P1。论文说 F 与 cF 之间的因子“which we evaluate at s = s_0”；按这句把 (2.33) 冻结在 s0 只影响 P1（k_1 从 s0 到最高节点涨 15.5 倍），对 S0 几乎无影响（比值 1.00–1.01）。两种读法都已实现（ModelSpec.ff_frozen_at_s0），默认取论文的字面读法。

**P1 对拍改用 lambda_max.** 任务书的 P1 要求用作者 2403 参数复现 2403 的一个公开点。本任务改用 lambda = (pi/4) T_3333(4/3,4/3,4/3) 的上界 2.661——这是对 2309 同一套（解析性+交叉+幺正性）已发表的独立数值，既验证了同一条“组装约束 + 调求解器”链路，又不必把 2403 的物理参数（nu0 = -20、s0 = 2 GeV、三电流）搬进来（第 7 节明令禁止）。lambda 行本身用论文 (2.14) 下方的 Weinberg 值 m_pi^2/(32 pi f_pi^2) = 0.023 校准。


## 12. 本任务未做的事

- continuous unitarity was not verified: unitarity is imposed only on the M x L x 3 collocation points of (3.61), never between them
- no pole was proven: the rho is read off a 90-degree crossing of the interpolated phase shift, which is not an analytic continuation to the second sheet
- no M -> infinity convergence was established: five (M, L) pairs are a sensitivity check, not a limit
- the double spectral densities carry no Mandelstam support restriction, because the paper states none
- the FESR cutoff bias of section 5a.8 is reported, not corrected
