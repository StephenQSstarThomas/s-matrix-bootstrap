# Watson 固定目标平台：纯代数与现有日志审阅

2026-09-11。对象是 `fine_watson_01_support_03` 最后 20 条 Newton 记录及其生产实现。没有重建状态、计算新方向、调用优化或求值，没有修改源码或接受门槛。下列建议尚未执行；后续诊断应放在既有 CLI 内，由主线在安全间隙处理。

目前没有发现可直接解释平台的线搜索符号错误、密度 Schur Hessian 漏项或严格约化域中的无界自由直线。代码支持的是一个有界、无精确零曲率方向的严格约化问题；这不排除很长、病态的数值通道。尚未确定当前平台的数值根因。

## 已保存的证据

[生产报告](fine_watson_01_support_03/report.json)和[日志](fine_watson_01_support_03/progress.jsonl)显示，最后 20 条为 iteration 31–50，μ=2.44140625e−7。decrement 范围 2.0212885274–2.9420682924，步长范围 1.7632455030–1.9980816997；最后一条 decrement=2.0230434331、step=1.9980816997。线性相对残差为 8.60647e−8–1.33532e−7。全部使用 fresh QR，没有命中新增的不精确方向分支。

日志 decrement 与 model_decrement 的最大相对差为 5.38104e−11，与 model_curvature 的最大相对差为 2.15383e−10。这排除了这些记录上明显的方向回传错位，但不是对原始函数值的独立有限差分核验。

20 步报告的 `barrier_merit_decrease` 总计 40.6988887687。该值使用 `state/future` 中的缓存线性量与递推 slack；不能直接改称“从同一 z 重新乘 H/current 后的原始 Φ 下降”。model_decrement 和 Hessian 特征也共享这些状态。密度和高能尾项则在 `state` 中从 z 重算，不能把所有族都归为缓存。

该 run 在前一 μ=4.8828125e−7 已保存一个合格中心；当前 μ 的最后 history 明确 `central_converged=False`。报告交付的已中心点与最后 `barrier_state.npz` 不是同一阶段资格，不能用报告顶层中心标签接受最后检查点。

最后检查点 SHA-256：`067565ee39c1f83d6c0127433e6d7cd49ddb2f236f695bb69669be8d815b49c9`。

## 线搜索与导数约定

`linear.center_joint_support` 最小化 Φ(z)=barrier(z)−cost·z/μ。特征定义满足 Hessian=FᵀF，`F.T @ rhs` 是障碍负梯度；所以 `g=cost/μ+F.T@rhs` 是 −∇Φ，解 FᵀF·dz=g 的符号正确。

`P.slope(state,dv,dz)−cost·dz/μ` 是沿线 Φ 的导数。扩张时仅把严格可行且导数为负的点接受为 lo；二分再以不可行或非负导数更新 hi，返回下降侧 lo。`old.value−new.value+step*cost·dz/μ` 与同一状态约定下的 −ΔΦ 相符。没有发现这里的符号不一致。

散射与 χ slack 的 `advance` 是对应仿射量的精确二次更新公式；FF、moment 和 Gram 使用仿射更新。实数代数中它与重新求值相同，浮点递推漂移仍是未在当前 μ 独立检查的数值问题。

## 密度 Schur 项没有明显的漏项

令 r 为代码中除以 B 的真实密度坐标，aᵢ=rᵢ²。密度障碍由提升变量 v 消元：

\[
\min_v\left[-\sum_i\log(v_i-r_i^2)-\log(1-v^Tv)\right].
\]

记 δ=1−vᵀv、dᵢ=vᵢ−aᵢ。一阶条件给 dᵢ=δ/(2vᵢ)，及 vᵢ=(aᵢ+√(aᵢ²+2δ))/2；这正是 `quotient.barrier_terms` 所解的标量方程。消元后梯度是 2rᵢ/dᵢ=4rᵢvᵢ/δ。

消元 Hessian 为 D+uuᵀ，其中

\[
D_i=\frac{2}{d_i}+\frac{8a_i}{\delta+2d_i^2},\qquad
u_i=\frac{-4r_iv_i}{\delta+2d_i^2}
\left(1+\sum_j\frac{\delta}{\delta+2d_j^2}\right)^{-1/2}.
\]

这与代码的 diagonal、rank 和两类特征行一致，原 z 坐标再分别除以 B。正的 rank-one 项是消元后的正确项；其向量整体负号不改变 Hessian。密度对角项在严格域内为正，不能产生精确的自由密度方向。

数值实现最多做 12 次 δ Newton 更新，没有独立记录最终标量方程残差。未见证据说明此处在当前检查点失效；记录该残差可作为未来 CLI 重算诊断的一部分，不能据此直接归因平台。

## 端点、自由列与严格约化域

1. `absorptive_change` 及两个投影坐标变换只替换振幅前 `free=1+2M` 行；密度坐标不变。端点的 T0 pivot 仍在此前 free 部分，另外四个 pivot 在 FF 部分，因此密度 Hessian 无需再乘一个被遗漏的密度变换。
2. `configure_endpoint_coordinates` 对 A/R/I/J/F/W 及 cost 同时应用同一可逆列变换，并重置 C 行视图。Newton 删除五个固定端点列；`free-sum(pivot<free)` 对应删除后的前段位置。`radial_coordinates` 只在其 eligible free 列回传增量，密度不变。现有日志的原特征方向曲率与解空间曲率一致到上述量级，未看到回传错误的直接证据。
3. 有矩界 R 的保留掩码由实际正矩权重产生；每个保留 R 在 Gram 对角线上非负，并出现在至少一条有有限上界的正权 FESR 中，故逐项有界。高能自由 R 被删除时，对应 Gram 也被删除，这是声明过的约化障碍；其原散射 disk 和高能 FF cap 仍存在。
4. 原生 S0/S2 的虚部恒等式覆盖两个单谱块。密度有界且 T0=0 后，这些原生 disk 给单谱有限界。每个 ImF 节点或者有保留 Gram 的直接虚部项，或者有高能 FF cap 的直接虚部项；有限原生节点上 kᵢ 非零，所以 ImF 也有有限界。因而剩余 C、ImF、保留 R 都有界，可逆变换不会生成新的无界直线。

更直接的 Hessian 零空间检查是：若一个允许增量使全部特征行为零，密度正对角块先强迫 dρ=0；T0 端点和原生 S0/S2 虚部二次特征再强迫所有振幅增量为零；低能 Gram 与高能 FF 二次特征强迫全部 dImF=0；最后 Gram 的谱对角变化强迫保留的 dR=0。因此实数代数中的约化 Hessian 没有非零零向量。该结论依赖现有节点／权重覆盖关系，不给出浮点最小奇异值的下界，也不排除近冗余方向。

## 最小后续鉴别建议

先保留主线正在运行的 step≤1 对照。主线已回报前四步 decrement 约 2.02→2.47→2.83→2.95，单步 merit 仍约 1.97；这不支持简单“α≈2 往复”的解释，不能把 2 附近的平台直接判为线搜索错误。

若仍需鉴别，等运行结束后在既有 boundary 加最小显式缓存重算记录：同一检查点、同一 reference、同一 μ，比较缓存与从 z 全重算得到的各族值／正 slack，并记录密度 δ 的标量残差。不要先重算 QR 或新方向，不改变约束、中心门槛或自动移动点。旧 [WATSON_CACHE_STATE_AUDIT.json](WATSON_CACHE_STATE_AUDIT.json)的较大 μ 结论不能转给本检查点。

若当前 μ 缓存也一致，再优先在已有 Newton 步内记录实际 dz 的各族 Dikin 曲率份额及各族障碍变化；这可使用已经计算的特征和沿线状态，不需新的 solver 或参数扫描。它能区分哪一族提供约 2 的障碍下降、哪一族决定长步边界，而单个总 decrement 无法判断。

## 源码身份

本次读取的六个文件与 support_03 报告保存的 producer 哈希逐项相同：

| 文件 | SHA-256 |
|---|---|
| linear.py | `24fcc257524fe6216315989826f900399a6d055638dfbd05a2203b6402d68484` |
| operators.py | `74cb905a1e9e881e3f3e5f7843d44e7e096e8c57e672f8f7abdbae336b78fc4e` |
| quotient.py | `32b669bec9e0b9af1c21bf464f2806411cdca87fffa6ff0593d7bb77dd3428e9` |
| endpoints.py | `6b4dacd438e750ca61f285744cc0a152d16e3090581e22c77962f23b867005db` |
| __init__.py | `453f182db710c6bf96dc18ffa46ce0ae8d6b8a4ce3e8e7983d0d2cb50f860eec` |
| imaginary.py | `9ad338b253d2941310401d245a24cba4ab546926e5c9c630f80c7bd647c1828c` |
