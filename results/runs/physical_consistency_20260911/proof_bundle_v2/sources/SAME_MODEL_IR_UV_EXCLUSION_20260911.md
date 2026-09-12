# 同模型IR→UV的严格分离：复用已有对偶，不另作优化

## 对象和结果

IR点是在物理参考xref=.07332139057293466处按最大化f11(3)预选的完整3876系数振幅。其数值中心、原式支持和解析primal均通过，且交付C与保存中心逐位相同。它与UV程序共用同一个3123盘准备、T0=0、五个领先高能必要条件、两条separate-L2手征球(.002)和实际双密度L4(B=377500)。只有UV加入归一化FF、电流Gram、打印FESR/raw绝对.002及原高能FF界；本份证书还明确采用F0、F1在z=−1处至少二阶零的端点完成，即两者的值与一阶导数均为零。这些FF端点等式不是原2309有限样本界自动包含的约束，不能删除后直接搬用本负界。

[IR支持](../runs/physical_consistency_20260911/interlaced_IR_ref_02/report.json)、[IR解析核验](../runs/physical_consistency_20260911/interlaced_IR_ref_analytic/source_audit.json)、[相位前选择](../runs/physical_consistency_20260911/interlaced_IR_selected/selection.json)。这里尚不借助该IR点的相移是否有ρ来选点或作证明。

复用已完成的[UV参考截面支持](../runs/physical_consistency_20260911/interlaced_UV_ref_support_01/report.json)的电流对偶，对这个IR点作零目标Farkas复核，得到

\[
U_{\rm fixed,H}=-0.1850389032871912<0.
\]

进一步计入原生解析行包络、精确运动学及FF/FESR输入，得到

\[
\boxed{U_{\rm fixed,analytic}=-0.18503890326560488<0.}
\]

后一个数是向外舍入的严格上界，完整包络和来源见[解析Farkas记录](../runs/physical_consistency_20260911/interlaced_fixed_IR_analytic/report.json)。两次复核均没有新优化。这是**该完整IR振幅没有满足所声明UV条件的电流扩展**的证书，不是完整联合域不可行——同一域已有其他完整可行见证。

## 为什么可以复用UV支持的对偶

固定C后，电流Gram具有形式

\[
G_q(C,v,R)=G_{q,0}+G_{q,C}(C)+G_{q,v}(v)+G_{q,R}(R)\succeq0.
\]

对任意已审计的半正定对偶Y_q，Tr(Y_qG_q)≥0。把这些不等式与FESR区间的支持函数、FF圆盘的支持函数相加，并对尚未精确消去的电流变量残差使用独立上界，可得

\[
0\le c_0+a\cdot C+r_v\cdot v+r_R\cdot R
\le c_0+a\cdot C+h_{\mathcal B_v}(r_v)+h_{\mathcal B_R}(r_R)
=:U_{\rm fixed}(C).
\]

这些残差上界来自正矩权重、Gram所给R≥|mathcal F|²、高能FF界和FF归一化；自由高能谱没有虚构上界，其残差符号必须合法。FF端点依赖变量按原仿射关系精确消去。必要时对保存的浮点Y作向外PSD修复并重新计算整个界，不能仅信任求解器状态。

原UV线性支持证书使用的正是这些电流协向量，并再加入散射、手征和L4的支持函数以消去自由C。固定C后不必重新寻找协向量；把其余C变量消元改为直接代入，就得到一个有效的零目标上界。等价地，原UV支持平面若严格排除该IR点，其差距可以分解为上述固定C上界与已满足的散射/手征/密度约束的非负支持余量，因此已有电流对偶已经提供一个可用的分离方向。

最终验收仍依靠完整的重新计算。本例得到U_fixed<0；若真的存在电流扩展，则必须同时满足0≤U_fixed<0，矛盾。对偶的负数大小可以通过整体缩放改变，不应解释为物理距离或ρ质量误差。

## 原式与解析复核的区别

第一份证书使用保存的float64 H。第二份在同一个完整解析函数族中包住固定C的原生分波值，并用精确K、运动学、打印十进制矩和hard-midpoint权重复核电流证书；没有把旧H证书直接搬成解析结论。

作为复核时的任意电流候选，程序可读取UV点的v、R；这不声称它们对IR点可行。零目标上界对所有允许v、R成立，结论不依赖这个任意候选。IR选择目录没有加入joint.npz或current.json，因而后续IR图仍明确只有散射/手征输入。

[接口修复前的失败](../runs/physical_consistency_20260911/fixed_IR_dual_before_fix/report.json)保留：旧dual入口错误要求IR目录也有joint.npz。新入口用单独的UV对偶供体并固定所给IR的C，正确执行已有数学证书；[原式重验](../runs/physical_consistency_20260911/interlaced_fixed_IR_dual/report.json)和解析重验分别保存。随后补丁加入工作进程读取前的物理算子输入哈希。新解析复核有完整算子哈希；较早的原式供体复核仍保留其inputs为空的历史记录，不能追认成当时已经记录了读取前哈希。源码快照及新解析证书的完整输入保持可复核。

## 不能由此推出的结论

这只排除一个预定IR振幅；不排除整个IR集合，不证明所有无ρ振幅都不可能，也不证明UV允许域只有一个ρ极点。ρ是否在所选UV点形成正确而稳定的结构，仍由同一完整振幅的相移、强度、FF、电流谱及分辨率检验回答。


## 删除新增FF端点条件后的独立复核

为避免把上述加强模型证书搬到较弱模型，又单独复用了保留的原基线E2_ref_verified对偶，固定同一个新IR完整C，在原A2_M50_L10/D1_M50_L10准备上重新计算。原1500散射行及十一辅助行与加强准备的原生部分逐位相同，连解析误差包络也相同，见[直接重验](../runs/physical_consistency_20260911/ORIGINAL_NATIVE_PREFIX_REPLAY.json)；没有跨PV／analytic处方转移证据。

这次电流模型不含新增的FF端点零条件，也不含新增高能散射约束；F(0)=1、原电流Gram、打印/raw/hard的四矩及原高能FF样本界保留。原式上界为−.1484359517640356，进一步[计入解析误差](../runs/physical_consistency_20260911/original_UV_fixed_new_IR_analytic/report.json)后仍有

\[
\boxed{U_{\rm original\ UV,analytic}=-0.1484359517403425<0.}
\]

因此原有限UV条件已足以排除这个新IR振幅，不依赖新增FF二阶端点条件。这里是另一份完整重算的证书，不是把较强条件下的负上界直接沿用；两份对偶不同，其负数大小也不用于度量条件强弱。它仍只排除该C，不证明所有无ρ振幅都被排除。
