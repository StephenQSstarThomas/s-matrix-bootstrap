# 本轮 Fig.11：当前计算进度

五组统一采用 PAPER_MAINLINE 的 combined-L2 χ=.002、hard-midpoint、原打印 raw 四矩容差.002、原 FF 输入；B(M)=100[M²+M(M+1)/2]。代表点统一为全幅度和全部电流方向的 +f00(3) 支持，无 Watsonian 和相位选点。已有 M50/L10 基准为 [paper_tip_path_02](../paper_tip_path_02/report.json)。

## M50/L8 已通过支持验收

复用原 L8 H，新建 [hard 电流算子](M50_L8_current_hard/report.json)，重新审计完整转移点后推进支持，没有继承旧可行性或旧证书。旧同 L8 保存点的 combined χ 超限，不能直接替代新模型内点。

初始小 μ 中心化缓慢，按明确决定从 [原式可行点](M50_L8_tip_01/report.json) 把数值 μ 提高16倍至1.5625e−5并重新中心；物理集合不变。[recenter_02](M50_L8_recenter_02/report.json) 连续完成五个真实中心，恢复半 μ 路径。再由实际 cache μ4.8828125e−7续算至 [M50_L8_tip_03](M50_L8_tip_03/report.json)，完整原式支持区间为 [.11780316557327108,.11789381781898335]，gap=9.06522457123e−5≤1e−4。

保存目标 x=.11780316557316134、y=−.007544025388345892；[selection.json](M50_L8_tip_03/selection.json) 在任何新 F 相移前冻结。全部原式联合约束通过，所有幅度方向保留。上界来自已完成 μ3.0517578125e−8 的中心，最后 μ1.52587890625e−8 只改善可行下界后即达到总gap；不误称最后中心已完成。

## M50/L12 尚无可复用真内点

[new hard current](M50_L12_current_hard/report.json) 已准备。按预定先后只检查 paper_tip_path_01/center_path.npz 四个已保存的 M50/L10 真中心，不根据相移选点。其 μ7.8125e−6、3.90625e−6、1.953125e−6、9.765625e−7 对应40、48、57、79个新增 L12 高波盘负裕量；最差−1.3300、−1.4252、−1.5528、−1.6431，均为明显原式违例，不是接近舍入误差的边界。

[候选诊断](L12_center_candidates.json) 保留逐点结果；最早 μ7.8125e−6 的完整点另经 [CLI resolution原式审计](M50_L12_center_transfer/report.json)，因散射盘不通过而未取得联合可行性。新增违例例如s≈649.79、I0、ell22，因此不能把 L10 真中心直接称为新 L12 内点。

随后按明确决定对旧同L12完整振幅只作固定7/8缩放，不缩放F(0)或旧current。[scaled_audit](M50_L12_scaled_audit/report.json) 通过新的combined χ=.00185632343497443及全部1800盘，ρnorm328914.35<B，未触发额外缩放。以这一完整振幅建立[单列池](M50_L12_initializer_pool/regions.json)，仅求hard电流的[joint_seed](M50_L12_joint_seed/report.json)包含188个求解变量；后端PrimalInfeasible，最终原式joint_audit未通过，幅度仍通过但返回current的FESR/FF/部分Gram未通过。此为单列初始化失败，不是完整L12联合问题不可行的证明；未将受限池用作最终支持。

L12 尚未启动完整支持，本代理没有启动M45/M60。本轮没有生成或查看任何新 F 相移；本代理没有运行中的生产进程。早先 M50_L8_tip_02 按暂停指令72秒终止、标记inconclusive，并未用于后来重新中心的两个批次。

## L12 完整电流 Phase I：第一批

按主线决定执行 [M50_L12_phase_I_01](M50_L12_phase_I_01/report.json)：完整参考电流来自最早L10中心的新H/current转移，已验7/8同L12振幅作为interior-coefficients；既有规则加入少量analytic seed形成严格散射内点。全部振幅方向开放，原散射/combined χ/L4始终为硬约束，仅current Grams/FF/FESR由辅助τ暂时放松。

24步、313秒后，实际phase_I_state.npz μ=.05676183183295568、τ=14.290899476629951；首个中心尚未完成，末Newton递减量7.7289。原始散射审计通过，但原始联合审计未通过；infeasibility_certified=false，不能宣称完整模型不可行，也不能把此点当作D联合见证。只执行所授权的一批并停止，尚未启动L12支持；没有新增缩放、父池或参数扫描。

同题续算 [phase_I_02](M50_L12_phase_I_02/report.json) 再完成29步、324秒；实际μ仍.05676183183295568，τ由14.29089948升至14.33840466，首个辅助中心仍未完成。递减量从7.73降至约2.03后长期停留2附近，末段步长约2；原始amplitude通过，joint未通过，仍无不可行性证明。因未出现授权条件要求的μ中心或τ向原始可行方向进展，未自行启动下一批。当前保留phase_I_02的完整实际状态；不更换初值或方法、不改物理约束。

按后续明确决定，[phase_I_03](M50_L12_phase_I_03/report.json) 仍使用phase_I_02同一完整cache，仅将辅助障碍μ降低至.005676183183295568（原值/10），作为推进可行性的候选策略；不宣称上一μ已中心，不调整物理约束，不用于最终支持。26步、313秒后实际τ=1.2822217613103417，较14.3384明显降低，但中心仍未完成，后段递减量约2；原始amplitude通过、joint未通过、infeasibility_certified=false。本批后停止，没有进一步减μ或追加批次。

## L12 预定四级短辅助序列

继续同一个完整cache，按明确指定的四个辅助μ，每级solver90秒/总180秒，逐次CLI，未追求上一μ的精确中心。每份report记录其辅助候选生成用途、请求/实际μ、真实τ及center列表；没有改变任何物理约束、容差或选点输入。结果见 [L12_auxiliary_mu_sequence.json](L12_auxiliary_mu_sequence.json)。

| 运行 | 实际μ | τ | 原amplitude | 原joint | 已完成中心 |
|---|---:|---:|---|---|---|
| M50_L12_phase_I_04 | 0.00056761831833 | 0.328211493408 | True | False | [] |
| M50_L12_phase_I_05 | 5.6761831833e-05 | 0.111108539625 | True | False | [] |
| M50_L12_phase_I_06 | 5.6761831833e-06 | 0.0696672104618 | False | False | [] |
| M50_L12_phase_I_07 | 5.6761831833e-07 | 0.0584160272722 | False | False | [] |

四级τ持续下降，未触发连续两级不改善的提前停止；四级上限用完后停止，没有继续减μ或新批次。没有一份通过原式joint，亦没有不可行性证明，故仍没有L12新D见证。06的原式amplitude失败来源是C_flat回代combinedχ=.002000000000126287，比.002超约1.26e−13；其缓存χ范数低于.002、全部盘最小裕量正、ρnorm小于B。审计保留失败，未用缓存严格性代替原式或放宽容差。

## 旧同L12早期完整点的有限复用检查

只检查旧stage_F的M50_L12_phase_I_11（该组唯一已joint-feasible Phase I/seed保存点）、tip_01、tip_02，未扫描其它分辨率或历史归档。三者新combined χ分别为.0021076382、.0020843962、.0021215125，均>.002；hard四矩也均不全通过，最大绝对raw误差分别.0119404482、.0053659747、.0044733931，均>.002。故这三个旧同L早期点均不能直接绕过新Phase I；没有进入CLI重复审计或新增缩放/父池/优化。完整数值与实际保存μ见 [L12_early_same_resolution_screen.json](L12_early_same_resolution_screen.json)。这个结论只排除列明的三个点，不排除完整模型可行性。

## 当前Phase I振幅的单次current fiber检查

只对phase_I_07的原C_flat作固定(1−1e−8)向内收缩，不缩放current/F(0)，不改ε。[fiber_audit](M50_L12_phase_I_fiber_audit/report.json) 原式通过全部1800盘、χ=.001999999979983、ρnorm358268.05<B，未触发额外缩放。以这一完整振幅建单列pool，执行一次[fiber_current](M50_L12_phase_I_fiber_current/report.json)，solver限30秒/总90秒，188变量；后端PrimalInfeasible，返回current的原式joint审计未通过，幅度继续通过。

实际未通过项见 [L12_phase_I_fiber_result.json](L12_phase_I_fiber_result.json)：四矩、FF及部分Gram失败。该报告只说明此次固定振幅的电流候选没有通过，不能据后端状态证明完整联合问题不可行，更不是最终支持。此次限额内工作已完成并停止，没有再找点/缩放/扫描或看相移。

## phase_I_07同μ继续推进的一批

[phase_I_08](M50_L12_phase_I_08/report.json)从07原完整cache继续，起始/最终μ均5.676183183295568e−7；没有使用fiber内缩点、手动跳降、换起点或改变物理输入。35步、315秒后τ=.05803372692755798（07末.05841602727216649），末递减量10.1392，仍无已完成中心，因此τ不是辅助最优值或不可行性下界。

原式joint未通过；原C_flat的χ近边界回代仍未通过：384-bit Arb对保存H和C_flat的独立点积给χ范数.0020000000000205030215，超ε约2.05e−14；longdouble粗算在这个量级不足以确定符号。全部盘裕量保持正，ρnorm<B。原zero-objective对偶上界=9.415050765205531e12，为很松的正上界，infeasibility_certified=false，不能证明原问题为空。完整缓存、原joint审计和原zero-objective对偶记录均保留在该run。指定单批已结束，没有自动续批或查看相移。
