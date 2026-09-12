# C阶段重放

所有科学计算仍用 `python -m smatrix_bootstrap.run`。从仓库根目录运行，本工作区Python为 `/home/shiqiu/miniconda3/bin/python`，环境可设 `SMATRIX_WORKSPACE=/home/shiqiu/s-matrix-bootstrap PYTHONPATH=/tmp/collocation_arb:src`。每次使用新的results/runs输出目录。

1. `select --chiral-tolerance .006`（另取.004/.002）读取既有B区域，在固定物理xref选择上边界内包插值，并检查原完整有限约束；保存coefficients.json与selection.json。不是新的支持优化。
2. `evaluate --coefficients <选择目录>/coefficients.json --primary-waves --energies ...` 输出三波。C1采用.01、.05、.10、…、3.95、3.99共81点；只做阈下求值。
3. `profiles --profile-kind subthreshold --profile-runs <三个C1求值目录>` 生成Fig.5和偏离／零点括区记录。
4. `profiles --profile-kind selection --profile-runs <绿色C1求值目录>` 生成Fig.6并原样复制完整系数。C2_fig6/coefficients.json是本次唯一主代表。
5. C3的evaluate仍使用该C2文件；能量列表为4加M50 preparation/amplitude.npz中满足.140√s≤1.2的原生节点。不能改成任意物理离节点。
6. `profiles --profile-kind phase --profile-runs <C3求值目录> --phase-reference results/runs/stage_C_reference_20260908/phase_reference.csv` 生成Fig.7、η和描述性对照误差。

本轮各命令均显式使用 `--unitarity-scope sampled --infinity free`；select/evaluate上限120秒，profiles上限60秒，实际单次约1–4秒。完整实际参数、输出路径和运行时间见[C_RESULT.json](C_RESULT.json)的runs列表；它保留所有绘图版本。权威图入口：C1_fig5、C2_fig6_display、C3_fig7_final；权威代表系数始终是C2_fig6/coefficients.json，显示版副本与之完全相同。

不要为图形差异重新筛选相移、改变规范或改用sine插值。Fig.7参考CSV是PDF读取数据，拟合线只在原显示范围内参与对照；没有外推或重新拟合。物理核与相移归一化见SCIENCE，B阶段的具体失败仍保留。
