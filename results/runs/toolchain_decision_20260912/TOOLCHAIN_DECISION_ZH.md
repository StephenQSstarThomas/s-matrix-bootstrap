# 工具链判定：MATLAB 不必需；Mathematica 作高精度核生成；SDPB 作高精度基准臂

2026-09-12，分支 `sdp-reproduction`。本文件只判定**软件平台**，不判定物理结论是否复现。
§3 是在定位到 M=50 不收敛的真实病因之后对 SDPB 的改判，以 §3 为准。
所有断言都附本目录内的探针产物或可复跑的命令。

## 结论速览

| 平台 | 本机状态 | 本任务是否必需 | 判定 |
|---|---|---|---|
| MOSEK 11.2.4 | 许可有效至 2027-09-12，`HOSTID=DEMO`（非节点锁）；Python 接口已装并实跑 | 必需（主求解器） | **已具备** |
| Python / CVXPY 1.9.2 | 已装，`installed_solvers()` 含 MOSEK | 必需（建模层） | **已具备** |
| Mathematica 15.0.0 | 容器 `wolframresearch/wolframengine:15.0.0` + `/home/shiqiu/Licensing/mathpass`，`$LicenseType = Professional`，无到期日 | 非必需，但是 §5a.10 交叉核对的**唯一**可行途径 | **已具备，应当启用** |
| MATLAB (≥R2019b) + CVX 2.2.2 | 未安装、无许可 | **不必需** | **判定放弃** |
| SDPB 3.1.0 | 镜像已在本地，端到端 1024 bit 验证通过 | 不是通量臂，但**是**高精度基准臂 | **基准臂**（§3 改判） |
| Mathematica 自带 SDP | 可用但无 `WorkingPrecision` | 不可用作求解器 | **只作核生成**（§4） |

## 1. MATLAB 为什么不必需

作者公开仓库的分工是 `Mathematica → MATLAB → Mathematica`。逐个看这三段实际承担什么：

- **第一段（Mathematica）产出物理内容。** `src/mathematica/gtb_qcd_01_generate_matrices.nb`
  的导出单元（用本目录 `mma_probe.wls` 的同法解压 notebook 可读到）把全部核写成
  **JSON**，不是 `.mat`：`partial_waves.json`（B1/B2/A1/A2/Bhat1/Bhat2）、`chisb.json`
  （fsigmachi/frhochi）、`form_factors.json`（Lambda3/KFF3）、`pion_coupling.json`
  （Asymsigma/Asymrho）、`spectral_density.json`（atorhoS0/P1/D0）、`sum_rules.json`
  （intS0/intP1/intD0/qcdfesr）。JSON 由 Python 直接读，无需任何 MATLAB。
- **第二段（MATLAB）不产出物理内容。** `papers/arxiv-2403.10772/GTB_numerics.m` 与
  `src/matlab/optimize_core.m` 是 MATLAB 源码（`%` 注释、`importdata`、`cvx_begin sdp`）。
  它做三件事：`importdata` 读上一段的核、用 CVX 声明锥约束、`cvx_solver mosek` 交给
  MOSEK。CVX 是建模 DSL，数值全部由 MOSEK 完成。**MATLAB 在这条链上是传输层，不是数学。**
- **第三段（Mathematica）只画图。** 本任务的图由 Matplotlib 原生输出。

判定的实质内容是：CVX 里那四个非平凡构造在 CVXPY 里逐个存在，且 MOSEK 解出的值与
解析值一致（`cone_probe.json`，四项全 `optimal`）：

| `optimize_core.m` 的 CVX 写法 | CVXPY 等价写法 | 探针结果 |
|---|---|---|
| `variable Bmat(3,3,3*n0) hermitian` + `>= 0` | `cp.Variable((3,3), hermitian=True)`, `X >> 0` | optimal |
| `norms([Reht,Imht],2,2) <= sqrt(2*Imhh)`（幺正盘，旋转二阶锥） | `cp.norm(z[:2],2) <= cp.sqrt(2*z[2])` | optimal，2.0 对 1.99999999972 |
| `norm(rho/Mrho,4) <= 1e2`（密度正则化） | `cp.pnorm(r,4) <= 1` | optimal，5^0.75 相对误差 1e-9 |
| CVX 内部二阶锥 | `cp.SOC(t,u)` | optimal，精确 5.0 |

`src/smatrix_bootstrap/sdp/problem.py` 已经在用 `cp.SOC` / `cp.pnorm(...,4)` / PSD 这三样
搭同一个锥程序，`A4_current_witness_mosek/report.json` 记录 MOSEK `optimal`、25 次迭代。

所以装 MATLAB 只能换来"用作者原生的那一层传输代码"，换不来任何新的数学或新的数值能力。
放弃 MATLAB 的代价是：报告里不能写"在作者原生环境内跑通"，只能写"在与作者相同的
求解器（MOSEK）上、用等价的建模层（CVXPY 替 CVX）跑通"。这个差别要在报告里明写。

## 2. Mathematica 已经可用，且应当从"可选"提升为"必做"

任务书 §5a.10 写的是"若可用 Mathematica 容器"。现在它可用：

```
MachineID  6503-44442-37175      <- 与 /home/shiqiu/Licensing/mathpass 首字段匹配
Version    15.0.0 for Linux x86 (64-bit) (May 26, 2026)
LicenseType Professional，LicenseExpirationDate {Infinity,0,0,0,0,0}（无到期）
符号积分、LegendreQ、40 位任意精度、读取作者两个 notebook（82 / 71 个输入单元）全部通过
```

完整探针见 `mma_probe.json`，脚本 `mma_probe.wls`，运行器 `scripts/mma/wolfram.sh`。
两个关键实现细节（否则会静默失败）：

- `--network none`：既离线，又让容器的 `$MachineID` 稳定，这是容器绑定的 mathpass 每次
  都能激活的原因。
- `-e HOME=/tmp`：以 `--user $(id -u)` 运行时若 HOME 不可写，`wolframscript` 仍会求值并
  写出文件，但**stdout 全部丢失**，看起来像空跑。

这件事改变了 §5a.10 的性质。此前"和作者的核对比"没有可行途径（`.mat` 文件不在仓库里，
只有生成它们的 notebook），所以交叉核对只能退化为"和本仓库旧算子对比"——而旧算子正是
被审计判定为在自建模型内打转的那一支，用它交叉核对等于自证。现在可以取到作者侧的核，
§5a.10 应当从可选提升为**硬门**：新实现的 H 与作者 notebook 导出的 JSON 核逐行比较，
差异 ≤1e-8 或逐项解释，否则不进入 §6。

注意两条边界，不能越：作者 notebook 是 2403 的参数（nu0=-20、s0=2 GeV、三电流），按
任务书 §7 只能借"约束如何组装、如何缩放"的模板，物理参数一律不搬；且 §1b 规定
Mathematica 只能是**核对对象**，不是实现来源——`src/smatrix_bootstrap/sdp/` 的算子仍必须
从论文 §2 公式独立推导。

## 3. SDPB：从"不上主线"改判为"基准臂"，理由是找到了真正的病因

初判基于成本，写的是不上主线。随后的审计把 M=50 不收敛的病因定位到**装配精度**而非求解器：
Λ₁₉(s₁)² = 1.03e−80，从 O(1) 中间量得到它需要约 95 位有效数字，双精度只有 16 位；
不缩放则约 1650/3000 行退化成 `0 ≤ 0`，按 Λ 缩放则可行集撑到 1e40。两种失败
`REPORT_SDP_ZH.md` §11 都实测到了。SDPB 全程任意精度，在这个动态范围下**不需要任何缩放
技巧**，因此它是双精度臂的独立裁判。成本判断不变（小时级/次 vs 分钟级/次），所以定位是
**基准臂**而非通量臂。详见 `TASK_SDP_PIPELINE_ZH.md` §1、§4。

已验证（全部本机实跑）：

- 镜像 `bootstrapcollaboration/sdpb:3.1.0`（SDPB 3.1.0，Boost 1_82 / Elemental 0.88-dev /
  FLINT 3.2.0-dev / GMP 6.3.0），2.07 GB，已在本地。
- 官方 1d 端到端例子：`pmp2sdp --precision=1024` + `sdpb --precision=1024` →
  `found primal-dual optimal solution`，dualityGap 7.8e−31，primalError 2.0e−87。
- **锥编码配方**（`sdpb_cone_encoding_pmp.json`）：把旋转二阶锥 x²+z² ≤ 2w 写成 3×3 PSD、
  常数项由归一化变量 y₀ 承载、线性不等式写成 1×1 块。已知答案 x = 2，SDPB 512 bit 返回
  2.0000…0007，gap 3.7e−31。
- MPI 路径 `-n 4` 正常（残留 `.ck` 与不同进程数不兼容，换并发前先删）。
- `scripts/sdpb/sdpb.sh` 运行器已签入并验证。

三个坑：SDPB 容器**不能**加 `--network none`（MPI 初始化永久挂起）；`docker pull` 慢但可用，
必要时走 `/playpen1/shiqiu/sdp-work/tools/pull_image.py`；PMP 必须用 degree-0 常数多项式，
非零次数会把幺正性施加到连续 s 上，那是比论文更强的问题。

## 4. MMA 自带的 SDP 不能替代 SDPB

`mma_solver_probe.json`：`SemidefiniteOptimization` 与 `ConicOptimization` 都**不接受
`WorkingPrecision`**（实测返回 `$Failed`），`Options[]` 里也只有 MaxIterations / Method /
PerformanceGoal / Tolerance——双精度封顶，买不到精度。`ConicOptimization` 能表达
NormCone + SemidefiniteCone 的混合（实测可解），但在随机生成的中等规模锥程序上直接
`$Failed`（n=50..1000 四档全失败）。

所以 MMA 在本项目的角色是**高精度核生成与交叉校准**，不是求解器。这恰好也是作者的用法：
`gtb_qcd_01_generate_matrices.nb` 里 `B1 = b/ΛI; Bhat1 = b/ΛI^2;` 把 Λ 的除法在高精度里做完，
导出 O(1) 的 JSON 核，下游只看到良态矩阵。

## 5. MATLAB 判定不变：不必需

第 1、2 节的论证不受影响——CVX 是建模层，MOSEK 是求解器，四个构造在 CVXPY 下逐个等价且
实跑一致。本轮新增的 SDPB 臂同样不经过 MATLAB。

## 6. 本轮任务书

`TASK_SDP_PIPELINE_ZH.md`。严格判定版 `TASK_SDP_REPRODUCTION_ZH.md` 及其已冻结的
`preregistration.json` 本轮不改。
