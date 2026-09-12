# 最后一组小补丁只读复核

结论：本轮χ行索引与固定fiber对偶长度修正正确；单位步长分支保持同一物理集合与原式接受判据，原默认行为不变。发现一个非当前_05的CLI漏口和一个报告描述遗漏，均不影响当前实际Phase I路径。未运行计算、优化或测试，未修改src/tests；完整106测试通过是root提供的验证结果。

- **χ行与n作用域正确。** `quotient.py:136`的 `H[-11:-3]`恰取最后11行中的8条χ行，不受新增散射行数影响，排除末尾的常数诊断与两目标行。固定fiber在同一 `if frozen`路径的第125行先定义 `n=(len(H)-11)//2`，第159行才传给 `zero_joint_duals(M,n,high)`；中间没有重新赋值n，非frozen路径也不读取这个局部n，故没有未定义变量或旧3ML长度残留。

- **toy回归保持原排除逻辑。** 两个现有测试都扩展extra=False/True；新增Re/Im散射行均减半，能量8→16后，其S值是原S与1的凸组合，系数为sqrt(3/8)，因此不会凭新增行制造另一个散射不可行原因。固定fiber测试还显式要求 `amplitude_primal_feasible`，并检查保存的kR/kI长度为实际n，使“只排除这个固定振幅的电流扩展”仍是检验目标。χ混合测试分别验证两父点违χ、共同混合满足χ，覆盖了尾部切片位置。

- **单位步长控制流保留可行与下降检查。** `linear.py:114–117`先检查a=1的trial严格可行且Φ′(1)<0，才把lo设为1并因新旗标停止扩张。对该凸障碍沿线函数，这保证[0,1]仍下降；若a=1不可行或斜率非负，就走原[0,1]二分，保留严格可行且导数非正的左端点。新旗标只限制搜索区间，不修改任何约束、端点、尾项、μ或最终原式核验。Phase I仍是它原有的τ放松内域，不是另外放宽原物理接受条件。store_true默认False，原扩张循环保持不变。

- **中低优先级CLI漏口。** `run.py:249`允许 `resolution --joint-feasibility --unit-newton-step`，但没有排除同时传入 `--profile-runs`。`operators.resolution_run:326`优先进入resolution_compare，这时单位步长旗标会被非Newton路径静默忽略。静态控制流可达；本次没有执行该命令。建议允许resolution的条件再加 `not args.profile_runs`，或统一拒绝profile_runs与joint_feasibility组合。其它prepare、evaluate、dual、普通chiral/pure及joint hull路径会被现检查拒绝。当前_05没有profile_runs，不受此漏口影响。

- **低优先级报告描述。** 启用单位步长的joint boundary仍经 `__init__.center_report`写出 `step_policy='Exact feasible line minimum with expansion'`。应按旗标改为[0,1]限制内的线最小化描述。完整parameters已经保存unit_newton_step，故来源未丢失；当前_05使用Phase I报告，不写这个字段。

本复核不预判单位步长改善居中速度，也不把测试通过当作物理成功。

## 输入指纹

复核时间：2026-09-11T08:44:22.803201+00:00。

| 文件 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/quotient.py` | `32b669bec9e0b9af1c21bf464f2806411cdca87fffa6ff0593d7bb77dd3428e9` |
| `src/smatrix_bootstrap/linear.py` | `3879f8c3d3faf8cd016f63649a5f14da197b4a59b1740d0d3a2642107ab33105` |
| `src/smatrix_bootstrap/run.py` | `75fde99da8579fcb1852f4d451372f2546db4b3ffd0646265f497e0874e0fb1d` |
| `src/smatrix_bootstrap/operators.py` | `5850eb2add0c8a8638b7b92bf4868b4d0a9015e74be3737a0b9f4006ae1a1c0c` |
| `src/smatrix_bootstrap/__init__.py` | `9e65f9de80eaceb49fe41feb6d29c8335cafa24f7fa588f46d1b63cd66ff6cbb` |
| `tests/test_mainline.py` | `45aab039e48800e7ad76cc1dd6e356a2d6b7d3ae6bf0615de3f93971f45ad929` |
