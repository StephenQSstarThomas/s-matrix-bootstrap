# He–Kruczenski pion / gauge bootstrap reproduction

最新执行合同：[combined-l2 + hard-midpoint](../../../results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json)。χ分组由原Fig.5三色残差指纹及作者后续代码共同支持；B/C正在重新验算，D新联合见证已通过，E/F新主线未完成。CLI默认χ范数已改为combined；历史separate复放须显式指定。后文旧阶段交付作为历史对照。

当前主线更新：D–F按原文统一π/M权重使用`hard-midpoint`重新对齐；[新D见证](../../../results/runs/mainline_alignment_20260909/D_RESULT_ZH.md)已通过，E支持正在计算。旧clipped结果保留，原论文核心claims仍未全部闭合；详见[STATUS](../../../STATUS.md)。

更新：2026-09-09。本项目复现[2309.12402v3](../../../references/2309.12402v3.pdf)的Fig.3–11。A–F的有限计算链已交付；[F报告](../../../results/runs/stage_F_mainline_20260909/F_RESULT_ZH.md)给出五组完整UV支持及[Fig.11](../../../results/runs/stage_F_mainline_20260909/fig11/fig11.pdf)。全部原约束、669项新原生主波检查和120项测试通过。五组均保留P1的90°上穿；三点稳健性、L依赖及中高能S0仍有核心claim缺口，不能标作全部通过。最新[逐阶段审计](../../../results/runs/claims_alignment_20260909/CLAIMS_ZH.md)和[逐曲线对照](../../../results/runs/claims_alignment_20260909/CURVE_COMPARISON_ZH.md)给出具体差异。先读[STATUS](../../../STATUS.md)、[SCIENCE](../../../SCIENCE.md)和[中文路线](../../../REPRODUCTION_GUIDE_ZH.md)。

| 文件 | 当前责任 |
|---|---|
| `__init__.py` | 公开类型、稀疏乘积、锥缩放与原生后端适配 |
| `model.py` | 幺正性、UV／FESR、Weinberg/相移、幅度求值与障碍函数 |
| `kernels.py` | sine／PV核、解析Legendre-Q角投影、原坐标支持外界 |
| `operators.py` | sine／PV散射provider、能量行、算子准备及固定算子的regulator比较 |
| `imaginary.py` | 高spin必要条件、原生散射锥与联合原式审计 |
| `linear.py` | 散射障碍求解、共同幅度／电流SDP、Newton线性代数 |
| `quotient.py` | 精确换元、Gram／moment锥、散射支持法向及数值换元 |
| `gauge.py` | 电流算子、全联合障碍函数及作者Watsonian目标 |
| `joint.py` | 完整联合初始化、可行Newton路径与原式支持界 |
| `selection.py` | C/E几何选点、Watsonian来源及原图曲线对应 |
| `figures.py` | C/E区域、相移、非弹性度及Fig.5–10数据／绘图 |
| `resolution.py` | F分辨率运输、完整电流Phase I、Fig.11与原式比较 |
| `analysis.py` | 支持计算编排、真实系数恢复、独立验算 |
| `io.py` | 处方／系数运输、有效支持汇总与区域 |
| `run.py` | 唯一计算入口与运行时间限制 |

当前15个核心文件按科学职责拆分；本轮用户允许在原十文件基础上增加少量必要模块。每个≤350行、24KiB；最多三个测试文件，每个≤350行。历史rank/PV执行已退休，源码与测试在[归档](../../../results/evidence/pre_global_sdp_core_20260907.tar.gz)，证明输入和失败保留。

```sh
python -m pip install -e '.[test]'
python -m smatrix_bootstrap.run --help
python -m pytest -q -p no:cacheprovider

# 准备作者collocation数值原型；保留自由常数
python -m smatrix_bootstrap.run prepare --prescription pv-midpoint \
  --nodes 50 --waves 10 --bits 384 --processes 1 \
  --unitarity-scope sampled --infinity free --output results/runs/NEW_PREPARATION

# 固定本轮条件模型；B取自已声明的后续作者代码设置，不按原图拟合
python -m smatrix_bootstrap.run boundary --preparation results/runs/NEW_PREPARATION \
  --unitarity-scope sampled --infinity free --solver barrier \
  --mode chiral --chiral-tolerance .002 --direction 1 0 \
  --density-limit 377500 --solver-seconds 600 --seconds 720 \
  --output results/runs/NEW_BOUNDARY
```

本工作区可用`PYTHONPATH=/tmp/collocation_arb:src`。活跃命令为`status/prepare/prepare-current/boundary/evaluate/dual/embed/regions/regulators/select/profiles/gauge-regions/select-gauge/gauge-phases/resolution/compare`。`--solver`选择barrier或clarabel；联合E/F支持使用centered（gauge模式的barrier为同一入口），原生线性后端可选auto/qdldl/faer；每次运行写入新的`results/runs`目录。`evaluate`检查保存的振幅，`dual`独立核验支持界，`regions`汇总同一设置下的有效支持。旧工作集派发、`--constraint-generation/--working-state`与独立`recover`入口已退休；源与测试保存在[D前快照](../../../results/evidence/before_stage_D_joint_20260908.tar.gz)，历史证明与失败保留。B求解内部的原式系数恢复仍在。

M50系数JSON保存3876个系数，并声明`unsubtracted sine-cardinal C_flat`或`unsubtracted PV-midpoint C_flat`。`--native-coordinates subtracted`仅改变求解器内部坐标；只有`--infinity zero`分支才用一条等式施加T0_unsub=0，free分支保留该自由度。双密度不变，输出恢复原坐标。`--subtracted`则用于明确声明的减除系数文件，两者不要混用。


原生FG条件采用标准加权SOS的精确全区间PSD表示；M50是六个25阶Gram块。独立Arb审计验证moment Toeplitz及νD/B原坐标残差。ρ上限、κ、手征／UV定义见SCIENCE，数值换元不改变它们；高spin必要条件仍不等于完整连续幺正性。

实际支持结果、regulator平台与原图差异统一见[STATUS](../../../STATUS.md)。原生solver状态、有效有限支持界和物理审核状态分别报告；准备出矩阵、通过实现测试或产生图文件均不计作论文交付完成。

`prepare-current`已完成真实M50的100个电流Gram、四矩和14个FF锥准备。`boundary --mode gauge --current-preparation ...`联立同一幅度与两种电流；`dual --mode gauge`从保存的完整联合解独立复核。D2固定clipped-phi、打印矩及raw绝对容差，见[D1/D2结果](../../../results/runs/stage_D_mainline_20260908/D1_D2_RESULT_ZH.md)。B1最新来源审计及原型／加强分支的区别见[STATUS](../../../STATUS.md)。

`--prescription`区分`finite-sine-cardinal`与`pv-midpoint`；`auto`读取保存的实际处方。PV只允许`sampled/free`，明确拒绝未定义的物理离节点求值和sine解析embedding。`--unitarity-scope strengthened`保留sine的加强路线。两种处方不会混用秩、可行性或支持证书。

`regulators --support-runs ...`对同一有限支持问题的不同密度上限汇总Arb支持区间、预设两decade窗口及可导出的图；它不选择最贴近论文的上限，也不把单方向平台推广到全部区域。

B计划已收敛为固定有限程序、必要区域查询和逐图差异报告，见[当前计划](../../../results/runs/stage_B_plan_audit_20260907/PLAN_ZH.md)。未改善M50的CVXOPT分支及其依赖已从活跃代码移除；历史运行及源码快照保留。

B最终交付见[报告](../../../results/runs/stage_B_delivery_20260907/B_RESULT_ZH.md)、[七域图](../../../results/runs/stage_B_delivery_20260907/regions/regions.pdf)和[系数／重放](../../../results/runs/stage_B_delivery_20260907/REPLAY_ZH.md)：204份有效支持，七域几何≤.01；原v3设置身份与连续幺正性未由此建立。

本轮[Fig.4与C总报告](../../../results/runs/stage_C_mainline_20260908/C_RESULT_ZH.md)、[Fig.5](../../../results/runs/stage_C_mainline_20260908/C1_fig5/profiles.pdf)、[Fig.6](../../../results/runs/stage_C_mainline_20260908/C2_fig6_display/profiles.pdf)、[Fig.7](../../../results/runs/stage_C_mainline_20260908/C3_fig7_final/profiles.pdf)及[单入口重放](../../../results/runs/stage_C_mainline_20260908/REPLAY_ZH.md)已交付。C2原样复用C1绿色完整振幅，C3未根据相移重新选点；0次新增bootstrap优化、计算顺序执行、98项测试通过。E/F计算与差异审计均已交付；尚未复现的形态与稳健性claims继续具名保留。

D阶段[报告](../../../results/runs/stage_D_mainline_20260908/D_RESULT_ZH.md)、[完整4076变量](../../../results/runs/stage_D_mainline_20260908/D3_hull_009/joint.npz)及[重放](../../../results/runs/stage_D_mainline_20260908/REPLAY_ZH.md)已交付。`--joint-feasibility --direction 0 0`从保存的同设置B幅度凸包构造联合见证，并按原1500盘、χ、L4、100 Gram、四矩及FF重新验收；它不替代E1全空间非零方向支持优化。`joint_feasible`才是完整联合可行性标志。全部101项测试及结构检查通过。

E阶段使用`boundary --mode gauge --solver centered`求解全部4076个原始变量对应的联合问题。`--fixed-x X --direction 0 500`求上截面，`--direction 0 -500`求下截面；截面乘子给出完整问题的支持法向。`--ray --direction 1 0`沿保存初值的投影射线求边界。同一截面、物理设置及坐标布局可由上一运行的`coefficients.json`及同目录`barrier_state.npz`续算。数值路径只改变坐标、障碍权重与Newton步，最终仍检查全部原约束。未奏效的列分解／受限对偶分支已退休，保留在[源快照](../../../results/evidence/before_E_final_solver_cleanup_20260909.tar.gz)。当前模块职责见上表。

本轮使用`SMATRIX_BLAS_THREADS=2`，实际线程数写入report。联合支持保存`barrier_state.npz`，F的初始化保存独立`phase_I_state.npz`。内层预算结束后可保留已核验的joint_feasible点，支持验收仍只看原式gap；外层超时或执行错误才记为inconclusive。

E阶段见[结果与差异](../../../results/runs/stage_E_mainline_20260908/E_RESULT_ZH.md)、[Fig.8](../../../results/runs/stage_E_mainline_20260908/E2_display_final/fig8_selected.pdf)、[Fig.9](../../../results/runs/stage_E_mainline_20260908/E3_phases_final/fig9.pdf)、[Fig.10](../../../results/runs/stage_E_mainline_20260908/E3_phases_final/fig10.pdf)及[η](../../../results/runs/stage_E_mainline_20260908/E3_phases_final/inelasticities.pdf)。`centered --fixed-x X --direction 0 0`（不加`--joint-feasibility`）求完整障碍中心供后续支持使用；零目标不算支持最优性。`dual --mode gauge --interior-coefficients DONOR`可用同算子的有效联合对偶给另一完整候选验界。近边界误差、相移精度和π分支不混为一谈。

E的[重放说明](../../../results/runs/stage_E_closure_20260909/REPLAY_ZH.md)与[完整argv](../../../results/runs/stage_E_closure_20260909/replay.json)可直接重画区域、重新求值或重求固定Watsonian目标。`--objective watson`建立新的旧FF/S目标，默认mu=.001；`--resume-objective`继续同一目标并自动读取checkpoint mu，显式`--start-mu`优先。新旧投影位置、η、谱份额和FF相位诊断一起保留；不通过更换相位分支或物理容差制造rho。

F的[单入口重放](../../../results/runs/stage_F_mainline_20260909/REPLAY_ZH.md)与[逐图图谱](../../../results/runs/stage_F_mainline_20260909/FIGURE_ATLAS_ZH.md)提供完整入口。`resolution --joint-feasibility`使用全空间电流Phase I，散射/χ/L4始终保持；`resolution --profile-runs ...`汇总五组已验收结果。M变化的运输仅是数值初值，必须在新问题上验回。

`compare --comparison-manifest PATH`从保存的profiles及参考CSV产生明确来源的差值、RMS和模180°诊断；不会改变幅度、输入或相位分支。见[本轮manifest](../../../results/runs/claims_alignment_20260909/comparison_manifest.json)。
