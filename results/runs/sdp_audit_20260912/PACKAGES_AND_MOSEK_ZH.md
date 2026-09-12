# MOSEK 配置与端到端复现的依赖核对

2026-09-12，工作分支 `sdp-reproduction`。

**当前 Python SDP 主线的安装依赖已齐，MOSEK 许可证和实际求解均已验证。原样执行作者公开的 Mathematica → MATLAB/CVX/MOSEK 工具链，还缺两个可用软件平台：Mathematica/Wolfram 内核与 MATLAB。** 包齐全不等于论文复现已经通过；共同 C 的 M50 IR/UV 主线及其严格验收仍在进行。

## 已完成的配置

- 使用用户提供的 `/home/shiqiu/mosek.lic`，配置到 `/home/shiqiu/mosek/mosek.lic`；目录权限 700，目标文件权限 600。许可证内容不进入源码、日志或结果包。
- 在 `/home/shiqiu/miniconda3/bin/python` 对应环境安装 `Mosek==11.2.4`、`python-flint==0.9.0`。不再需要 `/tmp/collocation_arb` 或临时 MOSEK 包目录。
- 本仓库在该环境中完成 editable install；`pyproject.toml` 的 `sdp` extra 声明 CVXPY、MOSEK、mpmath。可以直接运行下列入口，不必额外设置 PYTHONPATH：

```bash
OPENBLAS_NUM_THREADS=4 /home/shiqiu/miniconda3/bin/python \
  -m smatrix_bootstrap.run sdp selfcheck
```

- 修复 `--solver MOSEK` 将 Clarabel 关键字传入 MOSEK 的错误。现在通过 CVXPY 官方 `mosek_params` 接口传递线程数、最大迭代数、时间限制、原/对偶可行容差和相对间隙容差。实际选项保存到结果。优化算法由 MOSEK 提供。
- 入口导入审计没有加载 `linear/scattering/merit/ir/gauge/quotient/conic` 等历史自研求解模块。独立旧算子交叉测试不因此被当作生产求解路径。

## 已完成的实际验证

1. CVXPY 调用 MOSEK 求解一个 3×3 PSD 安装测试：`optimal`，目标约 .9999999999895，迹约束误差约 7.8e−12。这是软件安装测试，不是物理结果。
2. 新增真实 MOSEK 接口回归测试，先复现原有 `Invalid keyword-argument(s)`，修复后通过。
3. `python -m smatrix_bootstrap.run sdp selfcheck`：**155 passed，39.77s**，包含数学核对、MOSEK 接口及 Arb 符号控制。
4. [M50 电流见证](A4_current_witness_mosek/report.json)：MOSEK `optimal`，25 次迭代；384-bit Arb 严格通过 **700 主子式 + 4 FESR + 14 FF**。该见证明确仅属人工电流子系统、SR-b、冻结 FF 因子；不是共同散射 C，也不证明论文联合问题的 Slater 条件。
5. [依赖闭包](PYTHON_DEPENDENCIES.json)：检查 11 个直接依赖和共32个安装包，**无缺失、无声明版本冲突**。
6. [入口及 PDF 检查](ENTRY_AND_PDF_CHECK.json)：PDF 实际生成并重新打开成功；SDP 入口没有加载上述旧求解模块。

主线直接包版本：NumPy 2.4.4、SciPy 1.15.3、CVXPY 1.9.2、MOSEK 11.2.4、Clarabel 0.11.1、SCS 3.2.11、python-flint 0.9.0、mpmath 1.3.0、Matplotlib 3.10.9、pytest 9.0.3、HiGHS/highspy 1.15.1。此为本次实际环境，不声称恢复了作者未公布的软件版本。

## 作者原生工具链还需要什么

作者[公开仓库说明](../../../references/upstream-gauge-theory-bootstrap/README.md)明确要求 Mathematica≥12、MATLAB≥R2019b、CVX和MOSEK。2403 快照的 `.m` 第51–52行明确使用 `cvx_begin sdp` 和 `cvx_solver mosek`；`.nb` 负责产生矩阵。2309 本身未发布程序，不能把后续参数设置当成原文参数。

| 组件 | 当前状态 | 剩余操作 |
|---|---|---|
| Mathematica / 可运行原 notebook 的 Wolfram 内核 | PATH及常用安装目录未发现 | 需要安装位置和有效授权，再实际验证 notebook 求值 |
| MATLAB（≥R2019b） | PATH及常用安装目录未发现 | 需要安装位置和有效授权 |
| CVX 2.2.2 | 官方完整包已下载并解包，记录 SHA256 | MATLAB 可用后运行 `cvx_setup`；当前不能称为已启用 |
| MOSEK 11.2.4 MATLAB Toolbox | 官方 Linux 包已下载，SHA256与厂商一致，MEX文件齐全 | MATLAB 可用后运行 `mosekdiag` 并验证 CVX 调用 |

CVX 路径：`/playpen1/shiqiu/sdp-work/tools/cvx-2.2.2/cvx`。

MOSEK MATLAB 接口：`/playpen1/shiqiu/sdp-work/tools/mosek-11.2.4/mosek/11.2/toolbox/r2019b`。

收到两个软件平台的可用安装/授权后，可直接执行下列正式安装验证；这段命令尚未在 MATLAB 内运行：

```matlab
addpath('/playpen1/shiqiu/sdp-work/tools/mosek-11.2.4/mosek/11.2/toolbox/r2019b');
mosekdiag;
cd('/playpen1/shiqiu/sdp-work/tools/cvx-2.2.2/cvx');
cvx_setup;
cvx_solver mosek;
```

CVX 当前官方开源版不再要求另行申请 CVX Professional 许可证；MOSEK仍使用已提供的有效许可证。依据：[CVX官方发布说明](https://cvxr.com/cvx/download/)、[MOSEK MATLAB接口安装文档](https://docs.mosek.com/latest/toolbox/install-interface.html)。Octave不是受支持的替代方案，不用它伪称已跑通作者原环境。

图 PDF 由 Matplotlib 原生输出，当前交付合同的 Markdown 报告和图 PDF 不要求额外安装 TeX/Pandoc。后续若扩展成排版后的完整报告PDF，可另作交付工具选择，不将其设为物理计算前置条件。
