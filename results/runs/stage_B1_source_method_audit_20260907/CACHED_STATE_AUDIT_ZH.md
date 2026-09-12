# 同算子梯度与缓存状态审计

对象：ε=.0002 strip+，PV M50/L10、sampled/free、B377500、两组四维L2。只读取现有点/双权重并做算子乘积；未重复求解Newton、改核心或扫描参数。精确数值及方法范围保存在 [CACHED_STATE_DIAGNOSTIC.json](CACHED_STATE_DIAGNOSTIC.json)。

**新增匹配记录已定量排除首个坏中心的缓存漂移及算子回传作为主因：它们的支持尺度误差仅约10⁻¹⁰，而实际stationarity缺陷为3.587。** 小Newton decrement没有控制原坐标对偶恢复的精度。下文保留先前best-dual/final-primal分解及其限制，第6节给出匹配复核；第7节已用保存的dz完成直接因果对照，支持修复near-center分支过早break。

## 1. 实际间隙在哪里

使用 [precision结果](../stage_B3_pv_M50_eps0002_strip_plus_precision_20260907/report.json) 已保存的合法双权重及最终可行点，按原float64 H视为精确dyadic的规则，以Arb384重新计算残差、自由σ/T₀消元及配对（约4.8秒）：

| 分量 | 数值 |
|---|---:|
| 原H支持gap | 1.32340755581229 |
| density support B‖rρ‖₄⁄₃ | 1.4204921341994422 |
| 实际配对 rρ·ρ | .09708700344487866 |
| density Hölder gap | **1.3234051307545635** |
| 全部physical disk互补gap | 2.39698292119083×10⁻⁶ |
| 两组χ互补gap | 2.80748075744608×10⁻⁸ |

因此几乎全部间隙来自 B‖rρ‖₄⁄₃−rρ·ρ；Hölder配对率仅约.06835。普通LD对照还给出rρ与ρ³的L2余弦约.06590。[center结果](../stage_B3_pv_M50_eps0002_strip_plus_center_20260907/report.json)同样有density gap≈1.31504266，disk+χ仅约2.41×10⁻⁶。

ρ₂上三角非对角坐标的处理与主计算一致：H作用于C_flat，ρ₂的C_flat非对角为actualρ₂的2倍；密度对偶范数采用actualρ坐标，故相应残差乘2。没有把该配对混成完整对称矩阵的重复元素范数。

## 2. 不能把记录中的两个μ混为一个点

这两个run保存的最优dual来自 **μ=10⁻⁹**；candidate/native_candidate里的系数却是最终 **μ≈10⁻¹¹** 的点。没有保存最优dual对应的z和缓存。

日志还纠正一个停止口径：precision在μ=10⁻⁹时确实触发central_converged，decrement≈8.08×10⁻⁹；最终μ≈10⁻¹¹则以decrement≈1.194随300秒预算收尾，没有center。center版本对应最后decrement≈.0211，也未在10⁻¹¹中心化。因而不能描述为“已在10⁻¹¹中心化却有最佳gap1.3”。

μ=10⁻⁶已足以复现问题：precision的decrement≈7.19×10⁻¹²，而该μ内估gap≈3.495。不需要为了诊断重新跑完整μ序列。

## 3. 缓存偏移假说的数学边界

若缓存可写成x=Rz+δx、y=Iz+δy、χ=Cz+δχ，在当前点将δ视为常数，则其矩阵导数仍是R/I/C。对任意该缓存点，若梯度真正按同一算子平衡，应满足

\[
c-\kappa H_R^T k_R-\kappa H_I^T k_I-C^T y
=\frac{\mu}{B}\nabla_r\Phi_\rho(r),\qquad r=\rho/B.
\]

右侧在自由σ/T₀坐标为零；接近L4球边界时密度方向趋向ρ³。因此，单纯常数偏移并不自动产生与ρ³严重失配的密度stationarity残差。偏移可以通过dual加权进入原式互补误差，但这和密度残差恢复缺口须分别测量。

实际代码的radial步骤会同时缩放z及缓存，且重建slack；浮点误差也会逐步新增，因此不能未经测量把全部历史误差简化成一个固定δ。现有不匹配点对也不允许将第1节分解直接倒写成“最佳μ当时的缓存误差”。

## 4. Fᵀb与原算子直接梯度对照

没有历史缓存时，可以由保存的合法dual在μ=10⁻⁹反演对应的圆盘／χ barrier点，构造一个明确标记的代表性缓存状态。对圆盘，令w=‖k‖、D=√(μ²+w²)+μ，则G=(0,1)+k/D，物理slack=2μ/D。χ球同理，χ=ε²y/[√(μ²+ε²‖y‖²)+μ]。这是数学控制，不声称恢复了作者或旧run未保存的实际缓存。

在这个相同状态上调用literal `barrier_terms`、同CSR `extended_product`，比较c+μFᵀb与原A的直接梯度。密度项在两边使用同一z以相消。约3.95秒得到：

- 最大分量差4.44×10⁻¹⁶；
- 将该差通过相同自由σ/T₀消元后，B‖差ρ‖₄⁄₃≈8.52×10⁻¹⁰；
- 反演再导出的dual最大相对偏差约8.15×10⁻¹⁸。

因此，没有发现这一代表性状态中的Fᵀb分解符号、actualρ packing或κ回传代数错误足以解释1.3缺口。最后float64导出也不是本次主因：已有LD内估2.10806910057与独立Arb2.10806911620仅差约1.56×10⁻⁸。

## 5. 本次采用的短匹配重放证据

建议原raw/unsubtracted同题、从原strip+系数、起点μ=10⁻⁶，60秒以内只抓首个已center点。记录同一时刻、转float之前的NPZ：mu、z、t/c/inverseκ、state.x/y/slack/chi/chiral_slack/density_gradient、未自由修正kr0/ki0/y、修正后kr/ki/rρ；保留decrement/curvature/linear_residual/iteration，最好再保留已经算好却在central分支被丢弃的dz。无需保存F或额外求解矩阵。

直接诊断只需：

1. 以NPZ里精确表示的LD z重新计算Arb原H@z，测量缓存δ、slack内部一致性及其dual加权偏移；初始化Cflat JSON只作warm入口，不能替代该精确点。
2. 同z/state比较Fᵀb、原H直接梯度与μ∇density/B。
3. 分开量化未修正free残差、自由修正对rρ的贡献以及修正后的ρKKT缺陷。
4. 在这个同一点分解全部互补gap；若有dz，只评估已有full-step候选而不重新解Newton，以判断过早丢弃近中心步骤是否影响dual平衡。

**现在不应盲目宣布refresh、升导出dtype或再减小μ就是修复。** 匹配记录将决定最小修正究竟是原H仿射刷新、同算子的梯度/自由消元精度，还是当前中心停止条件与支持恢复精度之间的不匹配。所有最终接受仍由原H独立内外界判定。

## 6. 匹配first_center的直接证据

[同题短重放的center_first.npz](../stage_B3_pv_M50_eps0002_matched_20260907/center_first.npz)已保存同μ、同z的longdouble状态。μ=10⁻⁶，iteration=1，decrement=7.19060787×10⁻¹²，curvature=7.19060769×10⁻¹²，记录的preconditioned linear residual≈2.01858×10⁻⁶。

该审计NPZ的旧`y`字段被state.y覆盖；按照原代码唯一确定的公式2μχ/chiral_slack重建dual_y，未改写原NPZ。z按其LD二进制尾数／指数精确转换到Arb，不用初始化JSON里的float点替代。

| 同点诊断 | Arb384结果 |
|---|---:|
| 原算子μ-scaled梯度最大分量 | .000304919645649609 |
| 未修正free残差最大分量 | .000239488255431882 |
| 原梯度与保存μgradient最大差 | 3.98×10⁻¹⁶ |
| 两者差经free消元的B‖·‖₄⁄₃ | 7.31×10⁻¹⁰ |
| **真实stationarity经free消元的B‖·‖₄⁄₃** | **3.587250490260957** |
| 保存μgradient同一投影cost | 3.587250490433087 |
| 缓存偏移被raw dual加权后的绝对上界 | 1.45003×10⁻¹⁰ |
| 原H最小physical裕量 | 3.90684×10⁻⁸⁴ >0 |
| 原H两组χ平方裕量 | 2.29044×10⁻¹³；2.85999×10⁻¹³ |

LD逐项链检查还给出：μ∇density/B的B范数≈.10131421，实际修正后rho_residual的B范数≈3.59228526；二者之差的cost≈3.58725049。自由修正前密度stationarity缺陷cost≈377.80085，free修正贡献cost≈376.95691。gd缓存与同z重算逐项一致。

所以，原H梯度忠实反映了保存的梯度，缓存误差远不足解释该缺口；大量误差确实存在于待平衡梯度里。这个坏中心只发生过1次已接受迭代，也不需要“长期缓存累计漂移”才能出现。

最小后续辨因不是换μ或backend：取得代码已经计算却丢弃的dz，在同点检查μFᵀF dz−μgradient的原坐标投影误差，并只评估既有dz的α=1候选。若full step显著修复KKT，则中心分支过早break；若线性Newton残差本身仍大，则当前CG/QR验收量不足。任何修复都必须同时控制与请求support gap相联系的原算子残差，不能只凭Hessian度量里的decrement认定对偶已经平衡。该判断不更改物理约束或把未通过的支持问题改判完成。

## 7. 已计算full Newton步的因果对照：最小修复有了证据

[step重放的center_first.npz](../stage_B3_pv_M50_eps0002_step_20260907/center_first.npz)额外保存`newton_step`，与先前首中心同μ、同iteration、同decrement；dual_y字段已独立保存。本审计只评估这条**已计算**的方向，没有重做QR、CG或任何优化。

将α=1代入原有二次slack更新，并保持仿射增量缓存，不额外做radial_center：

| 量 | 丢弃这一步时 | 采用该full step后 |
|---|---:|---:|
| stationarity投影cost B‖P(μ∇Φ−c)‖₄⁄₃ | 3.58725049 | **5.36158×10⁻⁵** |
| 原梯度最大分量（μ-scaled） | .000304920 | 4.50281×10⁻⁹ |
| LD支持gap | 3.49505719 | .00408676963 |

该步的线性Newton残差μFᵀF dz−μgradient，经同一free投影后的cost≈**.000118063**，低于本次请求gap=.0025。因此这一步的inexact solve已足够；主要问题是near-zero分支直接break，丢弃了对对偶平衡仍然重要的更新。虽然局部Hessian度量中decrement仅7.19×10⁻¹²，数值坐标中的max|dz|仍为.0161768。

随后按**生产相同的float64导出方式**，调用既有原H `support_outer` 做Arb384验证，得到：

- primal_feasible=True；
- lower=.7811192482688103，upper=.7852060189863017；
- gap=.00408677071749；
- density support=.101314411910，已与期望μ∇density的B范数≈.10131421一致；
- 原H最小physical裕量为正，约3.90684×10⁻⁸⁴。

μ仍是10⁻⁶，剩余约.00409属于该μ的barrier量，**尚未满足本次.0025请求**；应继续既定μ序列，而不是把这一步当作完整B1或区域完成。

另一项对照尤其重要：同一个full-step点若简单fresh LD重算Rx、Iz、Cz，stationarity投影cost反而约3516.97。故不能用低精度仿射refresh替换原来的增量缓存；其小的值误差在窄slack下会重新扰动dual。原H的独立接受检查继续保留。

**最小实现修复方向**是：近中心时执行已经计算的可行α=1 Newton步，使用原二次缓存更新，再检查／退出；不要提前丢弃，也不要未经检查再做radial缩放或fresh LD刷新。若记录更新后的checkpoint，步前mu_gradient/dz不能伪装成步后梯度。现有μ序列、物理集合、L4规范及最终原H支持验收均不改变。

该结果只证明这个明确故障的修复依据，不证明所有后续方向或六个区域已经收敛；其它运行仍需实际内外界验收。

## 8. 已实现并实测的最小修复

root在raw分支将e=μgradient按同一native S0/S2与实S0锚点消去free残差，得到e·z=Σδk·G+rρ·ρ。精确算术下，任意两可行点之间该线性误差的变化不超过

\[
2\left(\sum_i\sqrt{\delta k_{R,i}^2+\delta k_{I,i}^2}+B\|r_\rho\|_{4/3}\right).
\]

这里单位圆盘的中心项在两点相减时消失；忽略χ／尾等额外约束只会使宽度更保守。实现将此数值误差宽度与请求gap/4比较：若未满足，执行已计算的strict可行full Newton步，保留二次缓存，跳过radial和fresh matvec，下一迭代重算梯度；满足后才接受中心。它是数值停止预算，最终证书仍由Arb原H裁决。

[同题实际M50修复结果](../stage_B3_pv_M50_eps0002_central_fixed_20260907/report.json)为lower=.7842651766942954、upper=.7847169067405557，gap≈4.51730×10⁻⁴<.0025；μ只推进到10⁻⁷。旧run更高的已验primal仍可用于同题内包，不用较差新primal替换它。

这项修复没有更改物理约束、μ序列、密度范数、变量数或新增求解后端。后续B3继续使用已验证冻结源 `/tmp/smatrix_B_central_fix_rekhgjje`；个别支持成功不代替六个区域的整体几何验收。
