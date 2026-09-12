# 原生峰值证书：限定只读审阅

2026-09-11。只读核对 `native_peak_bounds`、`native_rho_certificate`、gauge_phases 集成及现有 peak 测试。没有运行优化或新求值，没有读取任何新 PV 相位结果，没有修改源码或当前 worker。

**原生 PV、order=0、精确 float 节点别名的正常路径中，区间数学未发现阻断问题。** 两项通用入口范围缺口已向主线确认，留待当前生产 worker 结束后修补；不在运行中修改代码。

## 数学与数据核对

- `native_peak_bounds` 比较完整 Arb 区间的邻点差，只有左右两差都严格为正才授予离散峰；没有用 midpoint 比较代替证明。现有宽区间回归正确地区分可证明峰和无法判定情形。
- `decode_real_ball(None)` 在本仓库的序列化协议中表示精确零；`encode_real_ball` 对非有限值直接拒绝，因此这里并非把缺失／无限区间解读为零。
- P1 位于已检查的 `[S0,S2,P1]` 波序列第三列。由 S=1+iκf，强度公式 π²(1−4/q)(Re f²+Im f²)/4 正确；阈值 q=4 给精确零。平方和的区间依赖可能变宽，不会因此漏掉真实值而误授峰。
- FF 在工作精度中重算 `pv_matrix(M)`，使用同一冻结 current.json 的向量 ImF，计算 `(1+K·ImF)²+ImF²`，没有改用另一份 ImF 或浮点 current 矩阵。重新采用精确 K 不自动转移原浮点矩阵的可行性，代码已声明该范围。
- 若中心样本严格高于两个端点样本，则任何在该闭区间连续、且匹配这些区间值的完成，都在内部取得最大值。该条件性结论正确；它不证明存在物理连续散射完成、不证明峰唯一或极点。`continuous_scattering_completion_proved=False`、`pole_proved=False` 和不证明峰缺失的标记适当。输出浮点能量与 bracket 用于展示，不是极点质量区间或精确峰位置。
- gauge_phases 的三角色、共同模型／网格检查仍保留；IR baseline 无 current 时只生成强度峰项。此证书本身不替代代表中心、支持 gap 或原式 primal 的独立验收。

## 已确认、待运行结束后修补的范围问题

1. **节点匹配需与 provider 完全一致。** 当前 `allclose(rtol=1e−13)` 比 PV 与 analytic provider 的“完全相同 float alias 才映射到原生节点”规则更宽。合法的 analytic 近节点 off-node 求值可能被此入口吸入，再混用原生 q 与 FF 节点。待改为 provider 使用的完全相同 float 别名匹配，拒绝近似邻点；不更改已有节点或数据。
2. **本工具限制为 PV/order0。** 当前入口未拒绝 analytic 端点分支，直接读取显示用 dependent ImF，未做该分支定义所要求的 exact completion。主线已决定不扩展此工具：非 PV 在 gauge_phases 返回 `not_applicable` 说明，保留其既有图交付；PV 非零 endpoint order 则拒绝。当前正常 PV/order0 不受该遗漏影响。

以上是证书入口身份与适用范围问题，不是新增物理条件。两项仍是待办，不能把本报告写成修补已落地；当前 worker 的输入、源码与结果继续冻结。

## 验证与来源

已阅读扩展的区间 peak 测试；主线报告测试数仍为 107。本审阅未独立重跑测试。当前三个文件通过语法解析，行数分别为 296／311／350；两个核心文件 22405／23799 bytes。

| 审阅时文件 | SHA-256 |
|---|---|
| certificates.py | `d42788ee0c51af6a6aa4f5c9ef8406efffe629ecc460c78cfecb6d05a1b35ac0` |
| analysis.py | `ef81b00243e8c181d621218cd3239933789d8ae088e4b48b66a18620a7faa563` |
| test_mainline.py | `f9a232cfa6306f6eab8de9a167608f6e54b5ebed6c3044f76401739d959f194b` |

主线处理：上述两项已在当前worker结束后修正。新入口严格比较原生float别名，near-node不会替换能量；非PV返回not_applicable，PV非零FF端点模型拒绝。已有peak测试加入近节点误认回归。PV_peak_scope_pytest.log：107 passed。随后添加同一CLI中的已证支持平面排除，避免重复固定C优化。
