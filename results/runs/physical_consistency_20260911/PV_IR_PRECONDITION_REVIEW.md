# PV IR 散射预条件与失败快照：限定审阅

2026-09-11。仅对照 [修改前快照](PV_IR_precondition_before/manifest.json)检查 scattering、IR 恢复分支及对应参数化测试。未修改代码、运行优化或重建科学状态。

**未发现方向回传、曲率或 gradient-width 的正确性阻断项。** 残差门槛保持 1e−6；不能通过 fresh QR 的方向仍抛错，新增预条件没有授予失败状态中心资格。

## 坐标与宽度

- 固定 x 先消元。以 X=A[-2]、v=−X[1:]/X[0] 及 E=[v;I] 表示允许增量，F 与 cost 分别成为 F·E 和 Eᵀc；这保持 X·dz=0。
- `radial_coordinates` 再同步变换 F 和 cost，PCG 的方程、作用量与残差均使用这一当前矩阵。先 `restore`，再补回 dz0=v·dz，顺序正确。新坐标中的 F·d、cost·d 分别等于原 full_features·dz、原 c·dz，因此方向曲率和 decrement 没有遗漏回传因子。
- `full_features` 保留换元前矩阵；radial helper 构造新矩阵，不就地覆盖它。near 分支用原 c 和 full_featuresᵀb 重算完整负梯度，避免将预条件坐标梯度误用于原支持宽度。
- 固定 x 时减去 g0·X/X[0]，是在该截面上等价的协向量代表，并非声称欧氏正交投影；乘 μ 后再交给原坐标 width 计算，尺度正确。密度位置使用 free−int(section)，与删除第一列后的排列一致。线搜索实际使用恢复后的 dz。

新增 section=False/True 参数化保留原测试并检查组合换元的完整 F 作用、cost 点积与密度增量不变。现有数学检查仅等价压行。本审阅只阅读测试与解析语法，不独立宣称正在运行的整套测试已通过。

## IR 恢复分支

已核对 `zero or args.require_center` 分支：零 T0 或需中心身份时使用 recover_ir；普通 free-T0 non-center 支持恢复使用 legacy helper。此前关于较好 incumbent 的兼容性建议已实施。非零 T0 不被错误加入零常数等式，当前可行中心 C 的原样保留行为继续有效。

## 失败快照字段

fresh 残差失败时，在任何候选方向回传或状态更新之前保存当前 z、μ、缓存量及线性残差等诊断，然后抛错。raw(z) 使用已有反变换和密度打包，保存的是当前完整 C，不是失败方向移动后的新点。

`failed_newton_coefficients.json` 明确 `initialization_only=True`、`center_converged=False`。`center_identity` 只查询既有四个中心文件名，不读取新增 failed 文件；故当前调用不会把失败点接受为中心。

非阻断字段备注：`failed_newton_state.npz` 尚未显式带相同资格标志与坐标说明，解释它仍需相邻 JSON／worker。若后续要求该 NPZ 独立可读，可补齐这些字段；该建议不要求修改运行中源码，也不影响上述拒绝行为。

## 文件身份

修改前三文件均匹配快照清单；当前三文件通过语法解析。两个核心文件分别为 305／211 行、22745／16679 bytes，测试文件为 350 行。

| 当前文件 | SHA-256 |
|---|---|
| scattering.py | `657a60a013a6fe636c5308d345e68f6ad58c680e94f82bb0d39fdee0ad2781ef` |
| ir.py | `e4bb92ec3a8a962820790314ab3472e6566288b1fcd7273515963dc9e2a54537` |
| test_mainline.py | `6f98eadb5676137d00e39a9bf54b6e56703f4c9345b147655dcc0354b370ba34` |
