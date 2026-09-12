# 旧 clipped 问题：第二轮 tip 固定 Watson 目标的中心路径进度

2026-09-09。仅完成已经启动的 [W2_tip_path_01](W2_tip_path_01/report.json) 这一批；未追加第二批、其他角色或新的 Watson 目标，未修改任何源码或测试。

**最终原式区间为[0.2667777111904631, 0.27379687443024087]，gap=0.007019163239777748，大于要求的 .001。** 保存的4076变量通过全部原联合约束，但第二轮固定目标的支持尚未闭合，不能计作完成的第二轮 Watsonian 迭代。

该批属于旧的 `clipped-phi` FESR 问题。任务期间主线转向 `hard-midpoint`，因此依指示让已经运行的这一批正常收尾并停止续算；其可行性和上界不转移到新的 hard 主线。

## 冻结目标与初值身份

- 目标来自 [W2_tip_01/objective.npz](W2_tip_01/objective.npz)，按 `--resume-objective` 读取，没有重新计算 Q。
- 该第二轮目标原先由 [literal_W_tip_path_03](../stage_E_closure_20260909/literal_W_tip_path_03/report.json) 的完整第一轮 tip 解生成。目标的三通道、43节点、129项及 Lambda^-2 权重保持不变。
- 本批开始和结束的目标系数 SHA-256 相同：

```text
7bb74b8017651b6a0da232fcc9a8bbfa80568cdd351915d1289068f358e73377
```

- 数值 warm start 使用 [literal_W_tip_path_01](../stage_E_closure_20260909/literal_W_tip_path_01/report.json) 的 `barrier_state.npz`。缓存实际 mu=1.5625e-5，保存4062个工作变量及完整4076变量参考点。该目录较早的 mu=3.125e-5 中心已收敛，而缓存所在的最后 mu=1.5625e-5 当时尚未完全中心化；本批从这个真实缓存继续。
- 两者的散射/current准备、PV/M50/L10、B377500、两个 epsilon_chi=.002 的 separate-l2 球、原 FESR/FF 设置均相同。使用 BLAS2、χ障碍权重750、原精确可行扩展线搜索和现有QR刷新；没有约束或求解器变体。

## 实际中心与原式界

本批完成52个 Newton 步、36次QR、累计129次CG迭代。求解器计时309.796秒，监督器总计312.563秒；在既定300秒求解预算后的当前步和审计完成时正常退出，没有达到420秒监督期限。

| mu | 原式中心目标下界 | 同一对偶的原式上界 | 该中心的gap | 完全中心化 |
|---|---:|---:|---:|---|
| 1.5625e-5 | .2562014336301192 | .3183739701849491 | .062172536554829905 | 是 |
| 7.8125e-6 | .26130774195548206 | .29322487280315945 | .03191713084767739 | 是 |
| 3.90625e-6 | .2639832456845837 | .2803418410700904 | .016358595385506713 | 是 |
| 1.953125e-6 | .2654059875293796 | .27379687443024087 | .00839088690086126 | 是 |
| 9.765625e-7 | .2661555133123448 | 3775.371016993647 | 3775.104861480335 | 否，预算结束 |

mu 仅在完整中心后减半。最后未中心化点的宽对偶没有覆盖较早的有效上界；最终保留 mu=1.953125e-6 的上界。

旧可行点的目标高于这批所得中心，因此最终 [joint.npz](W2_tip_path_01/joint.npz) 的完整 point 与 [W2_tip_01/joint.npz](W2_tip_01/joint.npz) 逐分量相同，仍为4076变量。其坐标约为(.09249593704790876, −.006047397630659428)，原式目标下界 .2667777111904631。将这个较好的旧原始点与本批较好的对偶合并复核，得到本文开头的最终区间。

[joint_audit.json](W2_tip_path_01/joint_audit.json) 中 `primal_feasible`、`joint_primal_feasible`、`amplitude_primal_feasible` 均为 true；原1500散射盘、两χ球、actual-rho L4、100个current Gram、四项FESR和14个FF cap均按原式检查。`support_optimality_certified=false`。

## 已保存的可观测量记录

[watson_update.json](W2_tip_path_01/watson_update.json) 给出本批最终保存点相对于冻结来源的三波和电流变化。三波最大复 S 变化均为0，是因为保留了同一个更好的旧原始点，**不是第二轮优化已经收敛或到达 Watson 固定点的证据**。该文件没有修改相位分支。

[center_observables.json](W2_tip_path_01/center_observables.json) 另保存本批已经中心化并通过原式审计的相邻中心比较；未中心化的最后一级不计入该文件的中心变化序列。

| 已完成中心之间的mu变化 | P1最大复S变化 | P1最大相位变化（模pi，度） | 新中心P1最低eta |
|---|---:|---:|---:|
| 1.5625e-5 → 7.8125e-6 | .08877465358548908 | 2.5490085589114813 | .983939313714774 |
| 7.8125e-6 → 3.90625e-6 | .04060448314197537 | 1.164557746269437 | .9973957700329156 |
| 3.90625e-6 → 1.953125e-6 | .0189875966536088 | .5440756304282213 | .9992920744542974 |

这些变化量记录同一个冻结目标的数值中心路径；不是相移误差包络、极点证明或新的相移选点依据。没有用它们重新排名任何代表振幅。

## 可重放参数与结论范围

本批使用如下唯一入口；记录用于重放，不表示授权继续旧 clipped 问题：

```bash
OPENBLAS_NUM_THREADS=2 PYTHONPATH=/tmp/collocation_arb:src \
/home/shiqiu/miniconda3/bin/python -m smatrix_bootstrap.run boundary \
  --mode gauge --solver centered --objective watson --resume-objective \
  --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
  --current-preparation results/runs/stage_D_mainline_20260908/D1_current_M50 \
  --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --density-limit 377500 --chiral-tolerance .002 --chiral-norm separate-l2 \
  --chiral-barrier-weight 750 --fesr-cutoff clipped-phi \
  --moment-source printed --sr-error raw-absolute --mq-rule arithmetic-mean \
  --coefficients results/runs/mainline_alignment_20260909/W2_tip_01/coefficients.json \
  --interior-coefficients results/runs/stage_E_closure_20260909/literal_W_tip_path_01/coefficients.json \
  --direction 0 0 --start-mu 1.5625e-5 --gap .001 \
  --solver-seconds 300 --seconds 420 \
  --output results/runs/mainline_alignment_20260909/W2_tip_path_01
```

结果说明：使用合适mu的完整联合中心缓存后，同一第二轮目标恢复了稳定中心化，并获得显著更紧的原式上界。它没有证明第二轮固定目标已完成，更没有完成新的hard主线。原 `source.json`、最终对偶、完整原始向量及 mu=9.765625e-7 的未完成中心缓存均保留；当前没有本任务运行中的worker，也不追加下一批。
