# Phase-I初值修补审阅

限定审阅未发现本次L12新初始化的阻断问题。helper复制完整float64点，仅改最后2M个R中的非正项；C、ImF和正R保持，active/keep与所有R优化列不变。k²只是数值初始尺度，不是新增下界。

τ<0时，正定的G+τI及FF/矩辅助不等式逐项推出原式约束，最终仍以原式audit为接受门。有效v2缓存的reference_point/z覆盖备用初始化，helper不重置缓存；故续算日志中的replaced_indices不能解释成缓存被裁正。

实际失败_01没有任何Phase-I缓存。修补后的_02取得τ=−0.01319656824100681，原始完整primal通过。它是新的联合见证，尚须独立+x支持中心才用于相位。

审阅：center_rounding_audit；验证：verification_PV_phase_I_seed.json；失败、源码和负/零谱RED均保留。
