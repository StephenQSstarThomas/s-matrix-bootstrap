# E 原三点选择与投影精度复核

2026-09-09。只读核对原图标记、原 E2 冻结记录、实际生产报告及后续精度报告；未修改代码、未启动求解、未读取压缩历史归档。

## 结论

原三点均有完整联合可行性证据，选择先于相移；但原 mid 的自身投影问题没有收敛，其几何近边界保证来自借用 ref 的支持平面。现有最紧相关界只证明该 mid 在同 x 处距真正上边界不超过3.3244e-5，不能将这个几何数值转成 P1 相移误差。

后续 `precision_comparison.json` 检查的是投影等式释放后的冻结 Watson ref 目标，没有检查原 ref/mid 的投影极值精度。它不能排除原 mid 的未完成中心、原代表振幅或由其 FF 冻结的新 Watson 目标对三点比较造成的影响。

最小且最有信息的补算是：**保持原 mid 的 x、目标、全部物理输入不变，从其真实投影 checkpoint 继续，把原问题支持 gap 收紧至 .001，再比较同一原生节点上的完整振幅与 P1。** 不重新选择物理参数，也不按改善相移的程度选择新代表点。

## 原论文选点与本项目规则的差别

直接来源为 [Fig.8 标记读数](../../../references/figure8_selected_points.csv)及 [v3 正文](../../../references/2309.12402v3-source/prd_submission_2.tex)第1148–1149行。正文选择红色 tip，另外两个粉色点靠近黑色手征点；没有给出最近距离度量、目标法向或精确数值选择算法。

| 原图标记 | x=f00(3) | y=f11(3) | x−原黑点x |
|---|---:|---:|---:|
| 黑色手征参考 | .071313036955 | −.004754210915 | 0 |
| 浅粉 | .070983968456 | −.004420014278 | −.000329068499 |
| 粉色 | .075329811085 | −.004759010786 | .004016774130 |
| 红色 tip | .081124968389 | −.005252986472 | .009811931434 |

本项目在查看新相移前采用了具名几何规则：红点为 +x 支持；ref 固定在由140/92输入计算的 xref=.07332139057293466；mid 固定在 (xref+x_tip)/2=.08628472452183995。冻结记录见 [PRESELECTED_POSITIONS](../stage_E_mainline_20260908/E2_select/PRESELECTED_POSITIONS.json)、[E2 selection](../stage_E_mainline_20260908/E2_representatives/selection.json)，规则来源见 [SOURCE_E](../stage_E_mainline_20260908/SOURCE_E_ZH.md)。

| 本项目角色 | x | y | 规则 |
|---|---:|---:|---|
| tip | .099248058471 | −.006404544597 | 完整 +x 支持 |
| ref | .073321390571 | −.004319635045 | 物理 xref 截面的上侧 |
| mid | .086284724520 | −.005265816457 | 物理 xref 至本项目 tip 的中点截面上侧 |

两个区别必须保留：

1. 140/92输入计算的 xref 比原图黑点横坐标大 .002008353618；本项目没有通过重标物理输入消除这个差别。
2. 原粉点位于“黑点到 tip”横向距离的约40.94%，浅粉位于约−3.35%；本项目取0%和50%。本项目两非tip点的横向间距为 .012963333949，原图为 .004345842630，约相差三倍。这与本项目更长的允许区域有关。

这些是横坐标与规则的比较，不是假定某个未公布的二维距离度量。原文没有唯一选点算法，因此本项目规则不是已证明的代码错误；但它也不是原论文两个 near-black 点的数值身份复现。三曲线接近程度的差别尚混合了选择位置与振幅求解精度的影响，不能直接解读成对原论文稳健性的反证。

## 实际生产链与精度

| 角色 | 振幅实际生产来源 | 证书来源与已证范围 |
|---|---|---|
| tip | [E1_tip_016](../stage_E_mainline_20260908/E1_tip_016/report.json) | 原式 +x 区间[.09924805847073383,.09925920163675225]，gap=1.1143166e-5；其中心已收敛 |
| ref | [E1_ref_upper_warm_002](../stage_E_mainline_20260908/E1_ref_upper_warm_002/report.json) | 固定 xref、目标500y；最后 mu=6.25e-7 中心已收敛；全局法向(35.473409171440345,500)，支持gap=.003384815105620498 |
| mid | [E1_mid_upper_warm_002](../stage_E_mainline_20260908/E1_mid_upper_warm_002/report.json) | 固定 xmid、目标500y；30步、约431.48秒后仍在mu=6.25e-7，`central_converged=false`；自身上界约3.1946e7，不是有用的最优性证书 |

mid 的选点记录指向 [E1_mid_upper_bound](../stage_E_mainline_20260908/E1_mid_upper_bound/report.json)。这一运行的 `command=dual`、`optimization_performed=false`：它保留 mid 的完整原振幅，并借用 ref 的完整有效对偶。该借界合法，但不等于 mid 自身中心已经收敛。

本次另直接比较了三份 `E2_representatives/<role>/joint.npz` 的完整 point 数组：tip/ref/mid 与上述实际生产者的4076个分量分别逐项完全相同。因此，来源问题不在冻结时混合或替换了振幅，而在选取位置及原 mid 的求解精度。

二者使用同一个全局平面

```text
35.473409171440345*x + 500*y <= 0.44452698153157977.
```

对于任意通过原式审计的点 p=(x,y)，若 n_y>0，则其自身横坐标处的真正上边界 y_max(x) 满足

```text
0 <= y_max(x)-y <= (U-n·p)/n_y.
```

因此：

| 角色 | 该点全局支持下界 | 全局支持gap | 同x向上距离上限 |
|---|---:|---:|---:|
| ref | .44114216642595927 | .003384815105620498 | 6.769630211241e-6 |
| mid | .4279051096872424 | .016621871844337388 | 3.324374368867e-5 |

这些是距离的**上限**，没有证明实际距离等于上限。mid 的界约为 ref 的4.91倍。ref/mid 的实际 x 与预定 x 分别只差约−1.93e-12/−2.33e-12，因此主要问题不是横坐标舍入。

`E1_mid_upper_bound` 使用 `--gap .02`，所以记录中的 `support_optimality_certified=true` 与 .01662 的间隙并不矛盾；该布尔值表示通过当次声明的容差，不能当作已达到 ref 的 .005 标准，更不能当作相移精度保证。原 [E_RESULT](../stage_E_mainline_20260908/E_RESULT_ZH.md)已明确保留 mid 未完全中心化及几何界不等于相移精度的限制。

最新 [E1_regions_certified](../stage_E_closure_20260909/E1_regions_certified/regions.json)仍采用这两个原始上侧记录；本轮新增的下支界没有替代原 mid 的投影求解。

## 对 P1 和后续 Watson 结果能排除什么

原 [E3 审计](../stage_E_mainline_20260908/E3_AUDIT_ZH.md)记录：tip/ref 的P1首次90°线性读数约812.36/699.96MeV；mid在既定分支中没有上穿，最低eta约.174。该审计也已声明，本项目三点间距更大，尚未分离位置、未定作者设置、数值极值余量和离散化的贡献。

线性投影目标的支持gap不控制完整4076变量之间的距离，也没有给出 S 节点值的变化界或极值振幅唯一性。相位在较小的 |S| 附近还会更敏感。因此，现有几何界**不能排除**未完成的原 mid 投影求解对 P1 有显著影响；也**没有证明**其误差足以解释全部差异。

后续的 [precision_comparison.json](../stage_E_closure_20260909/precision_comparison.json)与 [E3_ref_precision/report.json](../stage_E_closure_20260909/E3_ref_precision/report.json)明确设置为：

```text
objective=watson, resume_objective=true, fixed_x=None, direction=[0,0].
```

它完成的是冻结旧 Q 后的 ref 线性 Watson 目标：gap=.0006468069643142949，P1最大变化约.08271°，90°读数约变化.029MeV。其结论限于这一次后续目标、这一中心的精度。它没有重算原 ref 的500y截面，也没有重算原 mid 的500y截面；原 mid 若因投影精度改变 FF/S，由其生成的后续冻结 Watson 目标也可能不同。固定旧 Watson 目标内的精度检查不能检验这个上游影响。

## 最小补算：只继续原 MID 投影问题

优先 mid 的理由是：它有原三点中最弱的上侧几何界，且自身中心明确未完成。这是基于既有求解证据确定的信息优先级，不是按相移改善程度挑点。先收紧同一个 mid，能避免同时改变位置规则与求解精度。

真实起点必须是 `E1_mid_upper_warm_002/coefficients.json` 及同目录 `barrier_state.npz`，**不是**仅借界的 `E1_mid_upper_bound`。只读检查确认该缓存含完整 `reference_point`、工作 `z` 和余量，mu=6.25e-7，工作 x=.08628472452183995。

建议一次有界续算命令如下；本文未执行：

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=/tmp/collocation_arb:src \
/home/shiqiu/miniconda3/bin/python -m smatrix_bootstrap.run boundary \
  --mode gauge --solver centered --objective projection \
  --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
  --current-preparation results/runs/stage_D_mainline_20260908/D1_current_M50 \
  --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --density-limit 377500 --chiral-tolerance .002 --chiral-norm separate-l2 \
  --moment-source printed --sr-error raw-absolute \
  --fesr-cutoff clipped-phi --mq-rule arithmetic-mean \
  --coefficients results/runs/stage_E_mainline_20260908/E1_mid_upper_warm_002/coefficients.json \
  --fixed-x .08628472452183995 --direction 0 500 \
  --start-mu 6.25e-7 --gap .001 --solver-seconds 360 --seconds 480 \
  --output results/runs/claims_alignment_20260909/mid_projection_precision_001
```

新目标容差对应同x向上距离不超过2e-6，约比当前 mid 的保证收紧16.6倍。除求解容差外，原投影问题、位置和物理输入均不变；不调用Watsonian、不加ray、不重新作受限凸包优化。

验收先看原式联合可行性和支持区间，再按原生节点、同一相位分支比较旧新完整 S、eta、P1相移与90°括区。保留旧 E2 三向量；把新结果标为原投影精度复核，不自动替换主图，也不挑更像实验的解。若该有界批次未闭合，则继续将原投影精度列为未决，不能拿后续 Watson 的小gap替代它。
