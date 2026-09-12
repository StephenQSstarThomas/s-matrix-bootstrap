# PV四组F算子与完整初值复用清单

本次仅核对已命名的F_paper报告及其直接引用的准备目录、文件SHA与shape。没有读取相移／曲线，没有优化、生成算子或恢复旧producer。全量绝对／相对路径、44个文件SHA-256及维度见[PV_F_REUSE_INVENTORY.json](PV_F_REUSE_INVENTORY.json)。

四组H/current均已存在，声明设置匹配本轮目标：pv-midpoint、angular_order=24、infinity=free、不添加端点／高能条件或额外采样；printed四矩各raw绝对.002、hard-midpoint、原质量／耦合／FF界。B(M)=100[M²+M(M+1)/2]。H/current可原样供新的separate-L2 .002核验使用；这不表示作者全部离散约定已恢复。

所有下列完整联合初值均来自**combined-L2 .002**，不是当前separate-L2的既有完成记录。只复用数值，不继承旧支持、中心、最优性或代表点身份。

| M/L | B(M) | H rows shape | 完整C／joint维数 | 完整初值目录（F_paper下） | 历史联合状态 |
|---|---:|---|---|---|---|
| 50/8 | 377500 | 2411×3876 | 3876／4076 | M50_L8_tip_03 | joint_feasible，support=true |
| 50/12 | 377500 | 3611×3876 | 3876／4076 | M50_L12_phase_I_08 | joint_candidate，joint_feasible=false |
| 45/10 | 306000 | 2711×3151 | 3151／3331 | M45_L10_selected | joint_feasible，support=true |
| 60/10 | 543000 | 3611×5551 | 5551／5791 | M60_L10_selected | joint_feasible，support=true |

H目录前缀为 `results/runs/stage_F_mainline_20260909/`；四组依次为 `M50_L8_prepare`、`M50_L12_prepare`、`M45_L10_prepare`、`M60_L10_prepare`，矩阵文件均为 `amplitude.npz`。

电流与初值目录前缀为 `results/runs/mainline_alignment_20260909/F_paper/`。每组必须配套使用其 `{Mxx_Lyy}_current_hard`，内含 `current_data.npz`、`gram_linear.npz`、`moment_linear.npz`、`ff_cap_linear.npz`。四组Gram矩阵shape依次为600×4076、600×4076、540×3331、720×5791；moment为4行；FF矩阵依次为42、42、42、54行。

上述每份初值目录均同时保存 `coefficients.json`、`current.json`、`joint.npz`；ImF与rho shape均为2×M，joint保留完整联合变量。M45/M60的selected报告是同分辨率原算子重验，分别直接读取tip_03／tip_05；通用CLI参数中的默认nodes=50不能代替实际M和矩阵shape。

L12保存的Phase-I tau=.05803372692755798，joint_feasible=false，support=false，无不可行证明。其完整变量只是待恢复可行性的候选，不能称为联合见证，也不能把旧combined问题的失败推广为当前separate问题不可行。

M50 L8/L12的三个稀疏电流矩阵SHA相同，但current_data.npz的SHA不同。保持各自H/current目录配对，不因M相同而替换。所有支持和中心资格均需在本轮问题及实际完整点上重新确认；本清单不宣称新的可行性或F完成。
