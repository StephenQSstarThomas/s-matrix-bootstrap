# D阶段重放

在repo根目录使用当前Python环境，必要时设置`PYTHONPATH=/tmp/collocation_arb:src`。每次输出必须为新的results/runs目录。物理设置固定如下；不要按运行结果改动。

```sh
python -m smatrix_bootstrap.run prepare-current \
  --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
  --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --fesr-cutoff clipped-phi --moment-source printed \
  --sr-error raw-absolute --mq-rule arithmetic-mean --bits 384 \
  --output results/runs/NEW_D_CURRENT

python -m smatrix_bootstrap.run boundary --mode gauge --joint-feasibility \
  --solver clarabel --native-backend faer \
  --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
  --current-preparation results/runs/NEW_D_CURRENT \
  --region-summary results/runs/stage_B_delivery_20260907/regions/regions.json \
  --coefficients results/runs/stage_C_mainline_20260908/C2_fig6/coefficients.json \
  --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --density-limit 377500 --chiral-tolerance .002 --direction 0 0 \
  --fesr-cutoff clipped-phi --moment-source printed \
  --sr-error raw-absolute --mq-rule arithmetic-mean \
  --solver-seconds 60 --seconds 120 --output results/runs/NEW_D_JOINT

python -m smatrix_bootstrap.run dual --mode gauge \
  --preparation results/runs/stage_B_pv_M50_L10_prepare_20260907 \
  --current-preparation results/runs/NEW_D_CURRENT \
  --coefficients results/runs/NEW_D_JOINT/coefficients.json \
  --prescription pv-midpoint --unitarity-scope sampled --infinity free \
  --density-limit 377500 --chiral-tolerance .002 --direction 0 0 \
  --fesr-cutoff clipped-phi --moment-source printed \
  --sr-error raw-absolute --mq-rule arithmetic-mean --bits 768 \
  --seconds 60 --output results/runs/NEW_D_AUDIT
```

读取`joint_feasible`及`outer.primal_feasible`验收。`sampled_feasible`只反映振幅部分；原生solver状态不能代替原式检查。完整变量为`joint.npz:point`，振幅为`coefficients.json`，电流为`current.json`，逐约束结果为`joint_audit.json`。零方向不会被报告为支持最优性认证。

`--joint-feasibility`使用B全部同设置可行幅度凸包，显式重新施加χ球，再联立电流。它只构造可行见证。去掉此选项并给非零方向才进入全变量联合支持计算；本次D未交付该支持最优性结果。

新核检查用`evaluate --nodes 50 --waves 10 --bits 768`、同一系数、四个手征能量与原50个物理节点；确切54个能量和所有参数保存在`D3_fresh_native_768/report.json:parameters`。不使用未定义的PV物理离节点延拓。
