# B 阶段交付索引与重放（条件有限几何完成）

本索引只整理已有计算的路径、参数和结果，刷新时没有启动优化或区域聚合。最终索引见 [support_index.json](support_index.json)，便于逐条浏览的表见 [support_index.csv](support_index.csv)。所有索引路径均相对仓库根目录；204 份完整 C_flat 系数留在原支持目录，没有复制进交付包。

## 当前交付范围

| 部分 | 已保存汇总 | 当前结论 |
|---|---|---|
| B2 pure | [aggregate_031](../stage_B2_PV_B377500_region_20260907/aggregate_031/report.json) | 38 个有效支持记录、36 个不同法向；内外包距离上界 **0.009675773346601383 ≤ 0.01**，条件有限区域完成。 |
| B3 六个 ε | [final](../stage_B3_PV_B377500_region_20260907/final/report.json) | 166 个有效支持记录；六个窗口均达到既定距离上界 ≤0.01，条件有限几何完成。 |

B 阶段的固定合并交付为 [七区域报告](regions/report.json)、[区域数据](regions/regions.json)、[PDF](regions/regions.pdf)、[PNG](regions/regions.png) 和 [比较数据](regions/comparisons.json)。合并报告含全部 204 份有效支持，`status=geometry_ready`，B2/B3 的 geometry flags 均为 true，六个 ε 未发现嵌套矛盾。所有求解 worker 已结束。这里的完成仅指逐项声明的条件有限几何；报告中的 `completion_claimed=false`、`physical_unitarity_accepted=false` 保留，不解释为原 v3 数值身份或连续幺正性已经恢复。

| B3 ε | 支持记录 | metric 距离上界 |
|---|---:|---:|
| 0.006 | 29 | 0.009812577971420494 |
| 0.004 | 27 | 0.009729681547872668 |
| 0.002 | 28 | 0.009849593699275495 |
| 0.001 | 27 | 0.009566841580099160 |
| 0.0006 | 25 | 0.009587225269165866 |
| 0.0002 | 30 | 0.009073726190648930 |

B2 固定快照是 `aggregate_031`，B3 固定快照是 `final`；本次重放以这些路径及合并报告为准。[B2 manifest](../stage_B2_PV_B377500_region_20260907/manifest.json) 和 [B3 manifest](../stage_B3_PV_B377500_region_20260907/manifest.json) 保留调度历史；其中旧的 `latest_summary`、进度状态不能覆盖已经落盘的最终报告。

B2 独立交付包括 [结果说明](../stage_B2_PV_B377500_region_20260907/B2_RESULT_ZH.md)、[pure PDF](../stage_B2_PV_B377500_region_20260907/pure_region.pdf)、[PNG](../stage_B2_PV_B377500_region_20260907/pure_region.png) 和 [原图比较](../stage_B2_PV_B377500_region_20260907/paper_comparison.json)。B2、B3 单独汇总的整体 `status=incomplete` 表示未包含全部七域；应分别读取其中的 `B2_geometry_ready=true`、`B3_geometry_ready=true`。合并七域报告已为 `geometry_ready`。

## 重放的同一物理问题

固定使用 PV-midpoint / Legendre-Q 源处方、M=50、L=10、原始 50 个能量节点的 1500 个幺正盘、3876 个原始 C_flat 变量、普通组合密度 l4 范数 B=377500、`unitarity_scope=sampled`、`infinity=free`。共同 preparation 为 `results/runs/stage_B_pv_M50_L10_prepare_20260907`，无 additional constraints、解析全局 FG 或五个尾条件。B3 采用两个分别施加的四维 L2 手征约束，ε 为 0.006、0.004、0.002、0.001、0.0006、0.0002。

这是固定且逐项声明的条件模型；这里的两个四维 L2 手征约束不等同于后续作者代码的 combined-eight 范数。有效支持区间来自保存的 float64 有限算子 H 的可行下界与 Arb 上界；没有把求解器状态或自报 gap 当作证书。历史粗区间可以参与内外包，索引保留 `pointwise_tolerance_met=false` 和真实 gap；区域完成判据仍为保守整体距离 ≤0.01。B3 的几何窗口为 x≥0，距离用既定手征 metric；保存的方向可能不是原坐标中的欧氏单位向量，重放时不得自行归一化。

有限采样通过不等于连续幺正性证明，也不等于恢复了原论文全部未公开处方。`physical_unitarity_accepted=false` 和原图比较的限定继续保留。原图 marker 比较是独立的图像读数比较，不用于选择 B 或求解参数。

## 仅重建区域与图

在仓库根目录执行下面示例会调用唯一计算入口 `regions`，读取最终合并报告中的同一支持列表，不启动 bootstrap 优化。使用当前已包含精确叉积角度排序修复的 `src`；该修复只影响几何后处理，优化器冻结源没有改变。输出目录必须是尚不存在的新目录。

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/tmp/collocation_arb:src python - <<'PY'
import json, subprocess, sys
from pathlib import Path
root = Path.cwd()
selected = root / 'results/runs/stage_B_delivery_20260907/regions/report.json'
report = json.loads(selected.read_text())
runs = report['parameters']['support_runs']
subprocess.run([sys.executable, '-m', 'smatrix_bootstrap.run', 'regions',
    '--support-runs', *dict.fromkeys(runs), '--nodes', '50',
    '--reference-directory', 'results/runs/stage_B_reference_20260906',
    '--output', 'results/runs/stage_B_delivery_replay_regions_NEW',
    '--seconds', '60'], cwd=root, check=True)
PY
```

`regions` 的物理模型签名来自各支持报告，而非其自身 argparse 中未使用的 `mode`、`infinity` 等默认参数。查看汇总的 `model_signature` 和每个区域，不应据这些无关默认值把 sampled/free 结果读成 strengthened/zero。原聚合输入列表包含所有当时送入的记录，拒收原因留在 `rejected`；索引只列实际进入区域的已验证记录。

## 按报告精确选择 boundary / dual 的输入

每个索引条目给出 `report`、输出 `coefficients`、`candidate_npz`、`source_snapshot` 及原始 `input_coefficients`。**完整重放参数在该支持的 `report.parameters`**；`command` 是位置参数，其余非空字段对应同名 CLI 选项，将下划线换成连字符。列表逐项传递，true 布尔值只传旗标，false / null / 空列表省略。必须把 `output` 改为新的目录。

`boundary` 重放会再次优化，应使用原参数中的 warm 输入，而不是自动将本次输出系数当作原输入。`dual` 重放会再次计算原式支持界，并可能运行声明的 dual LP；必须保留原 `parameters.coefficients` 的目录及其配套 `candidate.npz`。不能只拿一个输出系数文件替代原 dual 输入，更不能跨源处方、B、采样网格、手征 ε 或 strengthened/sample scope 移用上界。

以下是**需要重新计算时才执行**的完整参数重放模板。它从该运行已经保存的 `source.json/source_snapshot_utf8` 恢复十个源文件到临时包；不依赖易失的历史 `/tmp` 路径，不改当前核心，也不复制任何大系数。这里恢复的是被选中现有 B 运行的源快照，不能据此恢复退休的另一处方执行树。依赖库仍须由当前环境提供；源相同不保证不同机器上的浮点迭代逐字相同。

```bash
python - <<'PY'
import json, os, subprocess, sys, tempfile
from pathlib import Path
root = Path('/home/shiqiu/s-matrix-bootstrap')
selected = root / 'results/runs/stage_B2_PV_B377500_region_20260907/support_037_adaptive/report.json'
report = json.loads(selected.read_text())
params = dict(report['parameters'])
params['output'] = str(root / 'results/runs/stage_B_delivery_replay_support_NEW')
argv = [sys.executable, '-m', 'smatrix_bootstrap.run', params.pop('command')]
for key, value in params.items():
    if key == 'worker' or value is None or value is False or value == []:
        continue
    argv.append('--' + key.replace('_', '-'))
    if value is not True:
        argv.extend(map(str, value if isinstance(value, list) else [value]))
sources = json.loads((selected.parent / 'source.json').read_text())['source_snapshot_utf8']
with tempfile.TemporaryDirectory(prefix='smatrix_B_replay_') as frozen:
    package = Path(frozen) / 'smatrix_bootstrap'
    package.mkdir()
    for name, content in sources.items():
        (package / name).write_text(content)
    env = dict(os.environ, SMATRIX_WORKSPACE=str(root), PYTHONDONTWRITEBYTECODE='1',
               PYTHONPATH='/tmp/collocation_arb:' + frozen,
               OPENBLAS_NUM_THREADS='8', OMP_NUM_THREADS='8', MKL_NUM_THREADS='8')
    subprocess.run(argv, cwd=root, env=env, check=True)
PY
```

把 `selected` 换成索引中的 boundary 或 dual 报告即可按该次输入重放。`regions` 不写 `source.json`，应使用上方专门的聚合模板。source 快照与报告共同保留当时的数值算法、预算及参数；支持系数索引仅作导航，不替代原报告针对保存的 float64 有限算子的可行性、Arb enclosure 与原始来源记录。
