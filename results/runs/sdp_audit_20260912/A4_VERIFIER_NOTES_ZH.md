# A4：原式 UV 的独立 Arb 复验器修补

日期：2026-09-12。范围是 `sdp/arbaudit.py` 的 UV 复验以及正、负和无法判定控制。
本项不恢复旧 Newton 求解路径、不扫描物理参数，也不宣称完成 Fig.3–11。

## 1. 原故障与修补边界

旧 `_uv_audit` 将返回三个量的 `ff_asymptotic_bounds` 解包为两个量。
M50、`S=0, ImF=0, rho_hat=3` 的最小调用原样得到：

```text
ValueError: too many values to unpack (expected 2)
```

除此之外，旧复验器还复用了浮点运动学因子、目标和矩行，未传递冻结 FF 因子的模型选项；
FESR 使用 `abs_lower <= tolerance`，不能证明整球满足不等式；Gram 只选出一个
区间可比较的最小值，没有汇总全部主子式的判定。

现在 UV 复验直接用 Arb 重建节点、权重、K、运动学平方因子、打印目标及容差。
模型向量的 binary64 数值按其精确二进制有理数读取；若输入本身为 Arb，保留其不确定性。
`audit(..., ff_frozen_at_s0=True)` 向 `_uv_audit(..., frozen_at_s0=True)` 传递设置。
默认不变。另有显式 `eps_ff`、`m_q` 的真实调用控制，避免复验默认量覆盖模型输入。

## 2. 从原文到复验公式

一手来源为 `references/2309.12402v3-source/prd_submission_2.tex`：
`curF` 及 `h19/h20`（约701行）、`srnumbers`、`srnum/h43`（约1022行）、
`FFasym` 和冻结因子的括注（约1039–1050行）；对应论文 (2.33)、(2.56)、
(3.66–3.75)。不调用旧主线算子产生这些 UV 行。

采用 (m_\pi=1\)，(s_0=(1200/140)^2=3600/49\)。
节点是 (s_i=4/\cos^2[(2i+1)\pi/(4M)]\)，
代码 `w_i` 是 ((ds/d\phi)_i/M\)，因此原始 FESR 权重是 **\(\pi w_i\)**。
只保留严格可判定的 `s_i <= s0` 节点；若当前精度不能区分截止位置则拒绝继续。
M50时，低能43节点，高能7节点。

由 `curF` 平方独立得到

\[
k_0^2(s)=\frac{3\sqrt{1-4/s}}{256\pi^5},\qquad
k_1^2(s)=\frac{(s-4)\sqrt{1-4/s}}{384\pi^5}.
\]

保持 `ReF = 1 + K ImF` 和 `rho = k(s)^2 rho_hat`。
每个原始 Hermitian Gram 块的七个主子式分别为

\[
1,\ 1,\ k^2\widehat\rho,\ 1-|S|^2,\
k^2(\widehat\rho-|F|^2),\ k^2(\widehat\rho-|F|^2),\
k^2\{\widehat\rho(1-|S|^2)-2|F|^2+2\operatorname{Re}[S(F^*)^2]\}.
\]

所有主子式非负等价于该有限 Hermitian 块半正定。
这里利用正的 (k^2\) 因子表达**原始未缩放**主子式，避免额外的平方根运算；
没有用求解器矩阵的最小特征值替代原式检查。

打印的四个 raw FESR 目标用精确十进制常数重建为

| 波、矩阶数 | raw目标 |
|---|---|
| S0，0 | (4.4187\times10^{-7}s_0^2\) |
| S0，1 | (2.82014\times10^{-7}s_0^3\) |
| P1，−1 | (5.75484\times10^{-5}s_0\) |
| P1，0 | (2.69948\times10^{-5}s_0^2\) |

对矩 (m=\pi\sum_{s_i\le s_0}w_i s_i^n\rho_i\)，检查
`slack = tolerance - abs(m-target)`。SR-a采用raw绝对`.002`，SR-b/SR-c分别采用
目标绝对值的`.1`/`.2`；它们是三个具名合同，不能互相改写身份。
这里遵循打印输出方程，并未根据曲线或算出的矩回选质量。

FF分别检查平方余量

\[
2m_q^2\epsilon^{FF}-k_0(s_*)^2|F_0(s_i)|^2,\qquad
\epsilon^{FF}/2-k_1(s_*)^2|F_1(s_i)|^2.
\]

冻结设置取 (s_*=s_0\)，显式关闭时取 (s_*=s_i\)；Gram一直使用每节点因子。
这保留现有 `ModelSpec` 的冻结解释，不宣称作者原始程序已恢复。
默认 (epsilon^{FF}=6\times10^{-5}\)，

\[
m_q=\frac{4+7.3}{2\times140}=\frac{113}{2800}.
\]

默认质量由这个精确算术平均值重建；非默认质量和容差按传入数值的十进制表示读取，
该约定也写入 JSON。不把标量打印矩暗中改成RMS质量的重算矩。

## 3. 严格判定和序列化

对每个余量球，只有整球非负才是 `certified_pass`；整球严格负才是
`certified_fail`；跨零则是 `inconclusive`。每一类的`*_ok`都要求全部行严格通过。
没有“尚未证明违反，所以通过”的退路。

报告逐项保存 `ball` 可读字符串，以及 `lower/upper` **精确有理数端点字符串**。
字符串中的微小负数不会被 binary64 下溢成零，端点也不作十进制向内舍入。
旧的浮点极值仅用于显示；`form_factor_max_excess`现在明确标注为平方电流超额。
严格结论来自逐项余量及计数，而不是这些显示数值。

## 4. 控制和实际见证证据

新增16个UV控制，保留原文件6个既有交叉检查。独立验证包含：

- `S=0,F=1,rho_hat=3`：700/700 Gram主子式通过，4/4 FESR和14/14名义FF界违反。
- `rho_hat=2±1e-100`：512位下，100个行列式的符号分别严格正、负；
  含`2`的输入球产生无法判定。
- `±1e-400`以及跨零球：小于binary64范围的符号被保留；精确零合法通过。
- 四个目标和三类容差与独立有理数计算比较；四个矩与150位mpmath按原始`curF`
  因子及Jacobian的独立计算比较。FESR容差边缘`±1e-100`及跨边缘球单独控制。
- 冻结/每节点FF设置、`audit`传参、非默认质量/容差均有真实调用控制。

验证命令：

```bash
env PYTHONPATH=/tmp/collocation_arb:src OPENBLAS_NUM_THREADS=4 \
  /home/shiqiu/miniconda3/bin/python -m pytest -q tests/sdp/test_sdp_crosscheck.py
```

结果：**22 passed in 0.92s**。`git diff --check`通过；核心文件306行/15164字节，
该测试文件212行/10472字节，均在单文件限额内。

主代理另由统一CLI产生独立人工电流见证，证据为
[A4_current_witness/report.json](A4_current_witness/report.json)和
[current_witness.npz](A4_current_witness/current_witness.npz)。它采用M50、SR-b、
冻结FF、默认算术平均质量；Clarabel19次迭代仅用于构造测试输入。
384位原式重验计数为700/700主子式、4/4 FESR、14/14 FF，全部通过。
最小严格余量下端的显示值分别为约

\[
1.4069698530\times10^{-12},\quad
2.3851042066\times10^{-4},\quad
3.7134335903\times10^{-8}.
\]

人工见证的S并非共同散射系数C产生，因此它只验证电流子系统和复验器；
不能称为完整UV散射可行点、支持证书或论文共振结果。
散射聚合摘要、手征浮点行和完整原式primal/dual证明仍属于后续A6范围。
