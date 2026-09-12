# E阶段重放

全部计算使用 `python -m smatrix_bootstrap.run`。复用已交付的B散射算子与D电流算子；物理参数与D保持一致。每次输出须使用新的 `results/runs` 子目录。

```sh
export PYTHONPATH=/tmp/collocation_arb:src
export SMATRIX_BLAS_THREADS=2

e_boundary() {
  python -m smatrix_bootstrap.run boundary --mode gauge --solver centered \
    --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
    --current-preparation results/runs/stage_D_mainline_20260908/D1_current_M50 \
    --prescription pv-midpoint --unitarity-scope sampled --infinity free \
    --density-limit 377500 --chiral-tolerance .002 \
    --fesr-cutoff clipped-phi --moment-source printed \
    --sr-error raw-absolute --mq-rule arithmetic-mean "$@"
}

# 严格截面初值：两套完整联合解同时组合散射和电流变量
e_boundary --joint-feasibility --direction 0 0 --fixed-x .07332139057293466 \
  --coefficients results/runs/stage_E_mainline_20260908/E1_ref_ray_balanced_002/coefficients.json \
  --interior-coefficients results/runs/stage_E_mainline_20260908/E1_ref_ray_interior/coefficients.json \
  --seconds 90 --output results/runs/NEW_E_REF_INITIAL

# 随后的支持优化保留全部振幅方向
e_boundary --direction 0 500 --fixed-x .07332139057293466 \
  --coefficients results/runs/NEW_E_REF_INITIAL/coefficients.json \
  --start-mu 1e-5 --gap .005 --solver-seconds 360 --seconds 480 \
  --output results/runs/NEW_E_REF_UPPER
```

若该批未收敛，从新目录的 `coefficients.json` 续算同一截面，将 `--start-mu` 设为同目录 `barrier_state.npz` 保存的mu。保持物理输入与算子相同；每一步都有缓存，超时目录不能作为已验收结果。目标反向时优先使用同截面的零目标全变量障碍中心：去掉`--joint-feasibility`，用`--direction 0 0 --fixed-x X --solver centered`求中心；再从其缓存优化正/负y。零目标不计支持最优性。验收读取 `joint_feasible`、`outer` 的原式上下界及 `support_optimality_certified`，不按求解器名称或减量单独验收。

中间上支预定横坐标为 `.08628472452183996`。初值由已验收的参考上支与冻结的 `E1_tip_016` 完整联合解构成，仍使用上述 `--joint-feasibility --interior-coefficients` 接口；随后以 `[0,500]` 优化全部变量。参考下支固定原xref，目标为 `[0,-500]`。两者的目标和位置在查看新相移前已确定。

相移使用与C相同的43个M50原生物理节点及阈值，完整能量列表保存在 `E3_tip_evaluate/report.json:parameters.energies`。调用 `evaluate --nodes 50 --waves 10 --bits 384 --primary-waves --prescription pv-midpoint --unitarity-scope sampled --infinity free --coefficients ... --energies ...`，不引入物理离节点续接。`gauge-regions` 汇总实际联合支持与B的ε=.002基线；`select-gauge --coefficients E2_select/tip/coefficients.json` 在新输出目录重用冻结tip并固定另两点；`gauge-phases` 汇总同样三套完整解的δ、η与Fig.9–10。

精确的交付目录、各段计算参数和结果见[E报告](E_RESULT_ZH.md)及[机器摘要](E_RESULT.json)。每次运行已保留实际生产源码；初值凸组合不作为最终代表振幅，参考图和实验数据不进入优化目标。

## 从已交付数据重放图与三点相移

以下只读取已经冻结的解，输出放入新目录；每个计算仍通过同一入口。

```sh
python -m smatrix_bootstrap.run gauge-regions \
  --region-summary results/runs/stage_B_delivery_20260907/regions/regions.json \
  --support-runs \
  results/runs/stage_E_mainline_20260908/E1_tip_016 \
  results/runs/stage_E_mainline_20260908/E1_ref_upper_warm_002 \
  results/runs/stage_E_mainline_20260908/E1_mid_upper_bound \
  results/runs/stage_E_mainline_20260908/E1_ref_ray_balanced_004 \
  results/runs/stage_E_mainline_20260908/E1_ref_lower_centered_002 \
  results/runs/stage_E_mainline_20260908/E1_ref_ray_slack_001 \
  --output results/runs/NEW_E_REGIONS

python -m smatrix_bootstrap.run select-gauge \
  --region-summary results/runs/NEW_E_REGIONS/regions.json \
  --coefficients results/runs/stage_E_mainline_20260908/E2_representatives/tip/coefficients.json \
  --output results/runs/NEW_E_SELECTION

python - <<'PY_REPLAY'
import json, subprocess, sys
from pathlib import Path
base=Path('results/runs/stage_E_mainline_20260908')
energies=json.loads((base/'E3_tip_evaluate/report.json').read_text())['parameters']['energies']
for role in ('tip','mid','ref'):
    subprocess.run([sys.executable,'-m','smatrix_bootstrap.run','evaluate',
        '--nodes','50','--waves','10','--bits','384','--primary-waves',
        '--prescription','pv-midpoint','--unitarity-scope','sampled','--infinity','free',
        '--coefficients',str(Path('results/runs/NEW_E_SELECTION')/role/'coefficients.json'),
        '--energies',*map(str,energies),'--seconds','120',
        '--output','results/runs/NEW_E_PHASE_'+role],check=True)
PY_REPLAY

python -m smatrix_bootstrap.run gauge-phases \
  --profile-runs results/runs/NEW_E_PHASE_tip results/runs/NEW_E_PHASE_mid results/runs/NEW_E_PHASE_ref \
  --phase-reference results/runs/stage_C_reference_20260908/phase_reference.csv \
  --output results/runs/NEW_E_PHASES
```

`dual --mode gauge --coefficients CANDIDATE --interior-coefficients DONOR`连同同一散射/电流准备与物理参数，可重用有效联合对偶给候选验界。中间点采用此方式得到同x的y误差上限3.3244e−5；不把未收敛内部中心冒称精确极值。当前三套相移与图形差异都是交付的一部分，重放不会把它们修饰成论文曲线。
