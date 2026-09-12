# B3 combined 自适应补算：首个有界任务无几何进展

2026-09-09。仅完成 [adaptive_001](adaptive_001/report.json)，随后依“只有原式界/几何有进展才继续”的条件暂停。未启动其余5个预算内候选、epsilon=.004/.006任务或任何求解器变体，未修改源码。

## 本次固定问题

- epsilon=.002，combined-L2，PV/M50/L10、sampled/free、B377500。
- 法向由 [regions_initial](regions_initial/regions.json) 的最大剩余几何缺口给出：(-6.126745046352974,-265.30534472194347)。
- 完整初值为 [retarget_071/coefficients.json](retarget_071/coefficients.json)，请求支持gap=.0024999999999999996。
- 沿用原成功B的 `start_mu=.01`、barrier/unsubtracted，BLAS2；solver-seconds180、dual-seconds60、seconds300。完整CLI命令保存在 [manifest.json](manifest.json) 的 `adaptive_records`。

## 实际结果

| 项目 | 结果 |
|---|---:|
| 监督器总时间 | 263.5027秒 |
| barrier时间/步数 | 186.1055秒 / 19步 |
| 最后mu | .01 |
| 最后Newton减量 | 331.0725906957 |
| 最后CG次数/线性相对残差 | 64 / .0219778280 |
| 最终原式下界 | 1.1155708333314698 |
| 最终原式上界 | 684299.6135204639 |
| 原式点可行性 | 通过 |
| 请求支持精度 | 未通过 |

barrier未完成第一个mu中心；其数值工作点目标约.04086072，低于已有初值约1.11557084，因此恢复步骤保留了更好的旧可行点，并作原式检查。最终 `coefficients.json` 不能被误认为最新工作中心；工作中心材料另保存在 `center_first.npz`、`center_best.npz` 和 `center_first_coefficients.json`。

随后固定法向的normal-dual LP在约61.39秒达到时间限制，没有得到解，转用了保守轴向回退。该LP状态不是原combined物理集合不可行的证明；已有163份新可行点及本次原式可行点仍成立。

## 重新汇总与停止理由

[regions_001](regions_001/regions.json) 合并38条旧pure、163条通过的combined重验及本次1条有效但很宽的支持，共202条。三条旧χ舍入拒绝没有重跑，仍在原目录及manifest中保留；本次汇总输入仅采用已合格记录。

epsilon=.002的最大几何距离上界在补算前后完全相同：

```text
before = 0.0872544560134244
after  = 0.0872544560134244
```

`geometry_ready=false`，下一建议法向也完全不变。这一次既未产生更好的原始点，也没有收紧当前区域所需的外界。故不能把“运行成功退出且点可行”当作B3推进，亦不应在相同输入上连续重复剩余任务。当前保持暂停状态，等待主流程对这一明确的中心求解未完成情况作下一次有界决定。
