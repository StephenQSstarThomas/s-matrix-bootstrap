# B2最小纯散射见证：范围与结果

**已有一个新解析模型下的纯散射见证，足以在+x方向支持加入手征条件后的强收缩。** 这是一条只验算的顺序链，未执行优化或查看相移。

输入为 `results/runs/stage_B_prototype_pure_M50_blas8_20260907` 的全部3876项C和完整1500+1500+8项候选对偶。模型由本轮 `A2_M50_L10` 的 analytic-cardinal H决定，M50/L10、实际双密度L4界B=377500、原生幺正盘和自由T0保持；pure模式没有χ或电流约束。旧数据与新模型共享未减除nodal C坐标，不作DST、π或整体2转换，历史JSON标签保持不变。输出C经逐项比较与输入完全相同。

| 新核验 | 结果 | 证据 |
|---|---|---|
| 保存有限H上的纯散射primal／候选对偶重验 | primal通过；x=f00(3)的下界 1.606567808854946；重算上界 1248.019462175265 | [B2_pure_replay](B2_pure_replay/report.json) |
| 解析源行区间primal | 1500个散射盘、0条χ和1项L4全部通过；最小裕量 [1.487334199652e-85 +/- 6.91e-98] | [B2_pure_analytic](B2_pure_analytic/source_audit.json) |

解析目标记录为 x∈[1.6065678089 +/- 7.74e-11]，双密度L4范数约 105418.482491，低于B。旧 `values`、`margins`、可行性标志和上下界没有被转移；候选权重仅用于对当前H重新计算残差及支持界。

对照同一保存H、M/L及B下已经核验的[手征+x上界](B3_IR_tip_replay/report.json) 0.1059280880042031，这个纯散射点的x约1.60657，明显超出手征集允许的范围。因而在该具名有限程序中，纯散射集与手征集在+x方向有严格且显著的差别。这一方向性见证支持Fig.3→4的收缩主线；它没有给出完整pure区域、区域面积变化或纯散射最优端点。

重算的pure上界约1248，支持间隙很宽，报告明确 `support_optimality_certified=false`。这里只使用合法primal给出的下界，不将其写成已收敛的pure极值。解析primal验证的是原生节点；手征全局支持上界仍属于保存的有限H语句，两者的数值层次保持区分。没有额外证明连续全能量／全自旋幺正性。

两步均经 `python -m smatrix_bootstrap.run`，分别使用 `dual --mode pure --direction 1 0 --dual-seconds 0` 和不传current的 `source-audit --mode pure`；完整参数、输入哈希及代码快照在各自report中。均以nice15、`SMATRIX_BLAS_THREADS=1`运行，每步180秒封顶，实际约 6.13／7.73秒。两份报告均无source/input变更，未修改src或tests。
