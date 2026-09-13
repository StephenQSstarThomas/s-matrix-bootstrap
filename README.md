# He–Kruczenski pion / gauge bootstrap reproduction

> ## ⚠ 读之前：这个仓库有两条平行主线
>
> **(1) Newton 主线**（`results/runs/mainline_alignment_20260909/`、`major_claims_20260912/`）——
> 本文件描述的就是它。它已被 `REVISED_CLAIMS_ZH.md` 的只读审计判定为
> **"解的是一个自建的近似模型，不是论文的问题"**（论文没有的密度正则化 B、比字面更紧的 L2 手征球、
> 以障碍中心而非极值点作代表、自研 barrier–Newton 求解器）。它的 P0–P5 计划已不再执行。
>
> **(2) SDP 直解路线**（`src/smatrix_bootstrap/sdp/`）——当前在做的。求解器已改为 **SDPB-only**
> （任意精度），MOSEK/MATLAB/CVX 都不在链路上。
>
> **当前进度、运行方法、精度判断与待补缺口的唯一入口是
> [`HANDOFF_PHYSICS_AUDIT_ZH.md`](HANDOFF_PHYSICS_AUDIT_ZH.md)。**
> 本文件保留为历史记录，其中与下列事实冲突的陈述以 handoff 为准：
>
> - 本文件的「当前统一设置 combined-l2 + hard-midpoint / B=377500」属于 Newton 主线，**不是** SDP 路线的设置。
> - SDP 路线采用 2309 的**字面设定：无密度正则化 B**。实测：删掉从 2403.10772 搬来的 `‖ρ‖₄ ≤ 377500`
>   之后，M=50 纯幺正才第一次拿到认证点（带着它会在 981 个圆盘处求解失败并返回高 60% 的值）。
>

复现[2309.12402v3](references/2309.12402v3.pdf)的Fig.3–11：有限散射约束 → 手征对照 → 电流与QCD求和规则 → ρ相移 → 分辨率。**SDP 直解路线（当前在做的）先读 [HANDOFF](HANDOFF_PHYSICS_AUDIT_ZH.md)。** Newton 主线的历史记录见 [STATUS](STATUS.md)、[SCIENCE](SCIENCE.md)、[运行路线](REPRODUCTION_GUIDE_ZH.md)。

当前统一设置为[combined-l2 + hard-midpoint](results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json)，物理输入与选点先于新相移固定。C三步及D联合见证已完成；E三点支持与相移已生成，三条均有P1上穿，但三点稳健性和部分弹性/形态差异尚未闭合。B区域精度及F缺项按STATUS如实保留；E1已证明两侧收缩，但未复现强非对称。旧separate/clipped的完成状态不转移。

| 当前数据 | 入口 |
|---|---|
| Fig.5–7 | [C结果](results/runs/mainline_alignment_20260909/C_RESULT_ZH.md) |
| 电流/QCD共同见证 | [D结果](results/runs/mainline_alignment_20260909/D_RESULT_ZH.md) |
| Fig.9–10 | [E三点与差异](results/runs/mainline_alignment_20260909/E3_PAPER_RESULT_ZH.md) |
| Fig.11四组已完成及L12缺项 | [F部分交付](results/runs/mainline_alignment_20260909/F_paper/F_RESULT_ZH.md) |
| 求解器为何这样选择 | [原文、后续实现和实测依据](results/runs/mainline_alignment_20260909/SOLVER_DECISION_ZH.md) |

唯一计算入口：

```bash
export PYTHONPATH=/tmp/collocation_arb:src
export SMATRIX_BLAS_THREADS=2
/home/shiqiu/miniconda3/bin/python -m smatrix_bootstrap.run --help
/home/shiqiu/miniconda3/bin/python -m pytest -q -p no:cacheprovider
```

生产数据写入新的`results/runs`目录，具体argv保存在各report.json。PV支持显式使用`--prescription pv-midpoint --unitarity-scope sampled --infinity free`；χ默认combined-l2，历史两球复放必须显式指定separate-l2。原式可行性、支持误差、物理曲线及分辨率分别判定。

| 文件 | 当前职责 |
|---|---|
| `__init__.py` | 公共类型、数值运输和结果序列化 |
| `model.py` | 物理归一化、电流算子、观测量及原图数值比较 |
| `kernels.py` | 独立PV振幅、解析核与联合原式审计 |
| `operators.py` | 网格/算子、联合状态及分辨率运输 |
| `imaginary.py` | 散射支持编排、原式外界与PSD工具 |
| `linear.py` | B/E完整Newton与CG求解 |
| `quotient.py` | 密度障碍、FESR目标、换元、B对偶与D初始化 |
| `analysis.py` | C/E选点、区域/相移与Phase I启动 |
| `io.py` | 区域记录及Fig.3–11/η图 |
| `run.py` | 单一CLI及联合计算派发 |

十个活跃模块均≤350行/24KiB；三个测试文件保留79项独立数学/实现检查。原始Newton数学循环、联合约束与数据格式保留，sine/FG/五尾、旧全空间conic/CG/embed等执行分支已退休；[整理审查](results/runs/mainline_alignment_20260909/core10_review/REVIEW_ZH.md)和[原源码/测试证据](results/evidence/core_before_ten_modules_20260910.json.gz)完整保留来源。额外差值/FF-lift/Watson调试图停止自动生成，对应数值与诊断仍在。

[整理前README](results/runs/mainline_alignment_20260909/README_BEFORE_CONSOLIDATION.md)保存历轮交付与旧接口说明；历史源码和证明输入留在原位置。完整计算或测试通过不等于原论文全部物理主张通过，有限节点幺正性也不等于连续认证。
