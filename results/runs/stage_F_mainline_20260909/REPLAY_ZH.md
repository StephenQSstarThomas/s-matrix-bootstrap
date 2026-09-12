# F重放：单一入口与完整数据

主入口始终是`python -m smatrix_bootstrap.run`。标准环境安装见仓库README/pyproject；本机生产运行用`/home/shiqiu/miniconda3/bin/python`，`PYTHONPATH=/tmp/collocation_arb:src`及`SMATRIX_BLAS_THREADS=2`。

## 重画已验收的Fig.11

输出目录须取新名字；此命令读取五组完整审核和新原生相移，不启动优化。

```sh
python -m smatrix_bootstrap.run resolution --profile-runs results/runs/stage_F_mainline_20260909/M50_L8_evaluate results/runs/stage_F_mainline_20260909/M50_L10_evaluate results/runs/stage_F_mainline_20260909/M50_L12_evaluate results/runs/stage_F_mainline_20260909/M45_L10_evaluate results/runs/stage_F_mainline_20260909/M60_L10_evaluate --gap .0001 --output results/runs/NEW_F_FIG11
```

[replay.json](replay.json)列出每组完整系数、电流、联合向量、原生求值、实际支持生产者及其全部参数。[fig11/resolution.json](fig11/resolution.json)包含完整相移/η、支持区间、原式审核入口和原图90°读数。

## 重新计算一组完整支持

下面从已经通过全部原约束的M60联合起点重新计算；它不是缩小后的振幅子空间，全部5551个幅度系数和240个电流变量对应的约束都保留。

```sh
SMATRIX_BLAS_THREADS=2 PYTHONPATH=/tmp/collocation_arb:src /home/shiqiu/miniconda3/bin/python -m smatrix_bootstrap.run boundary \
  --mode gauge --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --nodes 60 --waves 10 --density-limit 543000 --chiral-tolerance .002 --chiral-norm separate-l2 \
  --preparation results/runs/stage_F_mainline_20260909/M60_L10_prepare \
  --current-preparation results/runs/stage_F_mainline_20260909/M60_L10_current \
  --fesr-cutoff clipped-phi --moment-source printed --sr-error raw-absolute --mq-rule arithmetic-mean \
  --coefficients results/runs/stage_F_mainline_20260909/M60_L10_phase_I_04/coefficients.json \
  --solver centered --direction 1 0 --start-mu .00001 --gap .0001 \
  --solver-seconds 600 --seconds 720 --output results/runs/NEW_F_M60_SUPPORT
```

内层预算是一次续算段，不是验收本身。继续时把`--coefficients`指向上段输出，并读取其`barrier_state.npz`中实际的`mu`传给`--start-mu`；正常支持沿半μ路径推进。只有`joint_feasible`和原式`support_optimality_certified`同时通过才选为Fig.11代表。不要用一个未完成段的最后临时对偶替换已经核验的较好支持界。

## 新M/L的必要顺序

1. `prepare`建立该M/L全部原生散射行；`prepare-current`绑定同一个H，固定上述UV选项。M50基准复用B/D原有准备，其余四组准备已在本目录。
2. `resolution --coefficients OLD/coefficients.json`作数值运输并重新审核。M变化时不继承可行性；L改变也检查所有新增盘。
3. 需要起点时用`resolution --joint-feasibility --solver centered`。散射/χ/L4始终保持，只对电流子系统用单τ Phase I。可用`--interior-coefficients`提供同分辨率的完整B振幅作数值起点；不按相移选择。新问题初始μ按τ尺度自动定标，续算读取`phase_I_state.npz`的μ。该缓存和Phase-I对偶不作为普通+x支持缓存/证书。
4. 全联合起点通过后，用`boundary --mode gauge --solver centered --direction 1 0`求完整支持。再用`resolution`原式复核生成选点，`evaluate --primary-waves --bits 768`在该M原生节点加阈值求S0/S2/P1，最后由上面的统一命令比较五组。

对新计算固定B(M)=100[M²+M(M+1)/2]，不改变物理容差、追加Watsonian或按相移挑点。精确的已交付向量是数据证据；重新优化可能给出支持容差内不同的完整振幅，不承诺逐位身份或相移误差界。

失败的投影池及旧初始化代码保存在各次`source.json`的`source_snapshot_utf8`中，实际输入和失败未删除；活跃代码只保留最终完整Phase I。无需恢复旧分支即可加载、原式审核或重画已经保存的结果。
