# E重放：同一物理问题与已保存向量

全部计算仍经 `python -m smatrix_bootstrap.run`，不新增实验脚本。依赖已保存的B散射准备、D电流准备和原E2三向量；这些位置在[replay.json](replay.json)中明确列出。该文件包含完整argv，原生43节点不用手抄。

本工作区环境：

```sh
export PYTHONPATH=/tmp/collocation_arb:src
export SMATRIX_BLAS_THREADS=2
```

`python`应使用装有NumPy/SciPy/python-flint的环境；本次使用 `/home/shiqiu/miniconda3/bin/python`。输出目录须未存在，可在argv中改为新的results/runs路径。

直接重画最终图，无需优化：

```python
import json, subprocess
jobs = json.load(open("results/runs/stage_E_closure_20260909/replay.json"))
subprocess.run(jobs["plots"], check=True)
```

`jobs["regions"]`重算7份有效支持的Fig.8和Arb截面；`jobs["evaluate"]`重算三个冻结完整向量的原生主波。它们复用已经完成的科学结果。

重新求每个Watsonian固定目标，使用`jobs["watson_start"]`。目标始终来自原E2对应角色，数值初值为完整analytic center；新目标默认mu=.001。三个问题独立，默认逐个运行，不以输出相移挑点。

如果时限到而原式gap尚未达到.001，用下面的同目标续算模板；替换PREVIOUS_RUN/NEXT_RUN。省略`--interior-coefficients`和`--start-mu`时，自动复用本目标的checkpoint及mu。显式start-mu始终优先。

```sh
python -m smatrix_bootstrap.run boundary --mode gauge --solver centered --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 --current-preparation results/runs/stage_D_mainline_20260908/D1_current_M50 --prescription pv-midpoint --unitarity-scope sampled --infinity free --density-limit 377500 --chiral-tolerance .002 --chiral-norm separate-l2 --fesr-cutoff clipped-phi --moment-source printed --sr-error raw-absolute --mq-rule arithmetic-mean --objective watson --resume-objective --coefficients results/runs/PREVIOUS_RUN/coefficients.json --gap .001 --solver-seconds 600 --seconds 720 --output results/runs/NEXT_RUN
```

`--resume-objective`不会重新计算目标FF相位。省略它表示新的一次Watsonian目标更新；未完成的固定目标会被阻止直接进入下一次更新。本次只完成一次作者2403映射，不自动迭代至任意人为相位阈值。

验收看report.json中的all_original_constraints_checked、outer.primal_feasible和support_optimality_certified；求解器状态标签不能替代原式界。完整向量在joint.npz/coefficients.json/current.json，固定目标在objective.npz，续算状态在barrier_state.npz。原始三解与精度对照分别保留，不按相移改写冻结结果。

原生产命令、数值预算和源码还保存在各run的report.json/source.json中；失败路径见[运行清单](ATTEMPTS.json)。重放当前程序复现同一有限问题和协议，新的时限分段可能返回不同的近最优向量；原保存向量用于确定的图形与物理比较。没有启动F的M/L扫描。
