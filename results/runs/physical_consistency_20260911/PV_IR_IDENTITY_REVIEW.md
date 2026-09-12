# Free-T0 IR 恢复与中心身份：限定审阅

2026-09-11。对照 [修改前快照](PV_IR_identity_before/manifest.json)与当前 `ir.py`，只读核对既有恢复和选择调用。未修改源码、worker 输入或旧 C，未运行新计算；主线报告的 106 项测试通过未由本审阅重复运行。

**当前 `require_center=True` 路径未发现阻断问题。** 原式通过的 candidate 由 `recover_ir` 直接原样返回；不再无条件乘 1−10⁻⁸。因此这里不会改变完整 C 或其已有固定截面坐标，也不再因这次人为缩小而失去与保存中心的身份对应。

`zero=True` 才要求 candidate 和 seed 的 T0 均为零。caller 显式传入 `zero=args.infinity=='zero'`；free-T0 PV 不被误施加该等式。零 T0 输入和解析坐标的上游检查仍在，原式 `audit_ir` 也使用相同 zero 标志。新增非零 T0 小 H 回归直接验证可行 C 被保留，不把该测试说成生产代表已通过。

当 candidate 未通过原式核验时，helper 仍只沿 seed–candidate 的完整系数段恢复，每个最终候选再调用 certify。任何此类修复均返回 `repaired=True`，caller 清除中心标志；不能把修复点称为原中心。helper 本身不是 fixed_x 投影器：段上的精确截面保持以两端在同一截面为前提，浮点舍入不由它消除。正常通过分支没有新增坐标变动；选择阶段还要求交付 C 逐值匹配一个保存的收敛中心及 candidate。

## 建议保留旧 non-center 支持行为

建议在当前 worker 结束后的安全间隙，把新增 free-T0 路径限定到 `require_center=True`，同时保留既有 `zero=True` 路径：

```text
if zero or require_center:
    recover_ir(..., zero=zero)
else:
    recover_candidate_segment(...)  # 既有 sampled/free 支持恢复
```

新 `recover_ir` 对已可行 candidate 立即返回，不比较 seed 的目标值；旧 helper 的 finish 会保留目标更好的已认证 incumbent，并比较径向和系数段恢复结果。实际已有 `test_zero_native_candidate_keeps_the_verified_nontrivial_incumbent` 覆盖“零 candidate 可行，但正目标 incumbent 更好”的情形。该测试仍直接调用旧 helper，因此不能证明新统一 caller 保留了原选择行为。

这是 non-center 支持的候选质量／兼容性问题，不是当前 require-center 路径的可行性漏洞。后者应优先保持实际保存中心的完整 C，不能仅为改善目标值换成未认证为同一中心的 incumbent。恢复旧 non-center 分支只保留其历史算法行为，不意味着其径向收缩突然具有固定截面不变的保证。

此范围调整尚未在本审阅中实施；不在运行中改代码或转移当前结果的来源身份。

## 身份记录

两个修改前文件均匹配快照清单；当前文件通过语法解析。

| 对象 | SHA-256 |
|---|---|
| 修改前 ir.py | `9fe83590027ee9bd611c954c1c7c9b7f333c469ab7cd242cd39ba97160fd419b` |
| 当前 ir.py（209 行） | `89beab1b11ad83be9f43db7ec664362cf756107e32b19f354bfe997789499423` |
| 当前 test_linear.py（350 行） | `4e93641ced14b3774e9c6a75bf430ffb2d92a71798ff30207adcbc8a9f4afb28` |
