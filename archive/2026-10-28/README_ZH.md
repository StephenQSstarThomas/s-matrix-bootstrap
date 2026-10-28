# 暂停归档：代码、公式与精简结果

归档日期：2026-10-28。用户决定暂停项目，保留代码、必要的中间公式与结果，删除大体积主线运行产物。该目录是研究资料归档，不是求解器检查点备份；没有重新开展科学计算，也没有把未完成的科学结论改成完成。

## GitHub 核查

清理前远端只有 `main`（`f467fbf95f85acd7389658504cd6126332cdc966`）与 `sdpb-2309-regularised`（`27f03c462be834efe553166bcc731b0500c06e5a`）。2309 分支包含 SDPB 求解和结果分析脚本、Mathematica 公式核对脚本、测试、报告以及内含 `DERIVATIONS.md` 和运行示例的代码包。

本地 `sdpb-2403-gtb` 的 `7517945` 比该 2309 分支多 77 个提交，原先没有上传。新的 `archive/paused-2026-10-28` 分支以该提交为基础，补齐当前代码、2403 推导、学习笔记和本目录的精简结果。

作者代码继续通过 `.gitmodules` 固定在 `hyfysics/gauge-theory-bootstrap` 的 `801684d9ece3de2918a20145178f569459081098`，不改为复制进归档包。归档时已核对该提交仍是作者仓库的 `main` 和 `arxiv-2505.19332` 标签；`arxiv-2403.10772` 标签指向 `c8409331d5df911f592e53e5a43628418a39b213`。

## 代码和公式入口

| 内容 | 位置 |
|---|---|
| 两条计算入口 | `python -m smatrix_bootstrap.run sdp …`、`python -m smatrix_bootstrap.run gtb …` |
| 2309 离散算子、约束、求解和验证 | `src/smatrix_bootstrap/sdp/`、`scripts/sdp/` |
| 2309 方程到程序映射、英文推导与运行示例 | `docs/reproduction_2309/gtb2309_sdpb_reproduction_code.zip` 内的 `DERIVATIONS.md`、`derivations/`、`OVERVIEW.md`、`examples/` |
| 2403 数学推导、先验界及证书定义 | `docs/gtb2403/DERIVATION_ZH.md`、`VERIFICATION_AUDIT_ZH.md`、`WELLPOSED_ZH.md` |
| 2403 独立核和锥规划 | `src/smatrix_bootstrap/gtb/` |
| Mathematica 核与作者 notebook 核对 | `mma/gtb2403/GTBKernels.wl`、`build_kernels.wls`、`authors_reference.wls` |
| 复现学习笔记 | `docs/gtb_learning_notes/GTB_REPRODUCTION_NOTES_ZH.tex` 及 PDF |
| 科学结果及局限 | `REPORT_SDPB_2309_ZH.md`、`REPORT_SDPB_2403_ZH.md`、两条主线的计划和预登记 |

代码与推导已保留，并不意味着论文全部结论已经复现。2309 的部分最优性和排除结论只有数值证据，范围见 2403 验证层审查；2403 Phase B 仍有未完成项，以原始报告、证书和台账为准。

## 精简结果包

这里采用两层归档。GitHub 保存完整代码、推导和关键结果；本机 `/playpen-hdd2/shiqiu/s-matrix-bootstrap-evidence-2026-10-28` 另保存约 880 MiB 的完整小证据包，共 57,652 个文件，包含更多历史试算记录和精确降基矩阵。后者没有上传。GitHub 中的 `local_full_retention_manifest.json.gz`、`local_full_archive_checksums.json` 和 `local_full_verified.json` 记录该本地备份的清单与校验结果。本地 Git 分支 `archive/full-local-evidence-2026-10-28` 也保留了完整小证据快照，仅为本机备份。

`retention_manifest.json` 逐文件记录本次 GitHub 精简包的原路径、字节数和 SHA-256。`archive_checksums.json` 给出每个压缩包或分卷的 SHA-256；大于 80 MiB 的包分卷保存，以避免 GitHub 单文件大小限制。包内同内容的大文件可用 tar 硬链接去重，原始字节不改写。

| 包名 | 原目录 | 保留内容 |
|---|---|---|
| `local_results` | `/home/shiqiu/s-matrix-bootstrap/results` | Git 尚未收录的推导、源码片段和小型结果报告；原已跟踪结果仍留在 Git 中，其余小型候选解在本地完整证据包中 |
| `sdpb2309_regularised` | `/playpen1/shiqiu/sdpb_regularised_20260913` | 2309 关键叶子报告、裁决、最终图、日志、启动脚本、全部最终解和求解器坐标／摘要 |
| `sdpb2309_mainline` | `/playpen1/shiqiu/sdpb_mainline_20260912` | 更早主线的根目录总结、推导和源码片段 |
| `gtb2403` | `/playpen1/shiqiu/gtb2403_runs` | 当前报告、证书、图表、论文数字化数据、最终解、核生成清单和 M=10 小型核 |

2403 的三个 C1 内界点 `trial_M50/inner{44,46,48}r`、三个 `phaseB/A/watson{44,46,48}i` 的精确降基矩阵保存在本地完整证据包中，不在 GitHub 精简包中。GitHub 包保留 `trial_M50/relax_full`、`wellposed/W40_46_relax`、`wellposed/W60_46_relax` 和小型验证用例的 V2 矩阵 `X_matrix_*.txt`。W40/W60 和 M=50 大型核需根据保存的清单重新生成。只有 GitHub 包时，降基叶子的证书是历史收据，不能声称不恢复原降基矩阵就能原样完整重验；最终解和分析结果仍可查阅。

没有保留巨大的 PMP 输入、SDPB 二进制分块、MPI 检查点、热启动目录、大多数大矩阵和重复输入打包。恢复资料包不会恢复这些文件，也不能继续原来的热启动任务。需从代码重新生成输入后冷启动求解。

`/playpen1/shiqiu/gtb2403_prep` 的约 65 MiB 论文及公式准备材料另保留在本机；关键 Fig.3 提取数据另附于本目录 `paper_data/`。如在另一台机器重建，请取得对应论文版本及其源文件，并使用已保留的提取脚本。许可证、个人 Python 环境和 Docker 数据目录不属于本次删除范围，也未上传到 GitHub。

## 校验及解包

先按 `archive_checksums.json` 校验文件 SHA-256。未分卷包可直接解开：

```bash
mkdir -p /absolute/path/restore/gtb2403
tar -xzf archive/2026-10-28/gtb2403.tar.gz -C /absolute/path/restore/gtb2403
```

如果对应包实际为分卷，按文件名顺序合并后解开：

```bash
cat archive/2026-10-28/gtb2403.tar.gz.part* | tar -xz -C /absolute/path/restore/gtb2403
```

包内路径相对于上表的原目录。结果里的绝对路径保留原样，以免改变来源和哈希；运行分析时显式传入解包后的 `--leaf`、`--kernels` 等路径。依赖其他叶子绝对路径的流程需恢复对应原目录，或在新机器建立兼容软链接。旧自动队列和旧检查点续算不属于恢复目标。

## 从代码重新计算

```bash
git clone --recurse-submodules --branch archive/paused-2026-10-28 \
  https://github.com/StephenQSstarThomas/s-matrix-bootstrap.git
cd s-matrix-bootstrap
python -m pip install -e '.[test]'
```

数值依赖及归档时的版本见 `environment_versions.json`。SDPB 使用 `bootstrapcollaboration/sdpb:3.1.0`；Mathematica 核使用 `wolframresearch/wolframengine:15.0.0`，需要有效的本机许可证，脚本支持 `MMA_LICENSE_DIR`。代码中的部分历史分析脚本仍含本机绝对路径，需按实际解包位置设置。

2403 核从 `mma/gtb2403/build_kernels.wls` 生成。示例仅说明参数顺序；精确重建某次运行时，以对应 `kernels/*/manifest.json` 中的 M、L、精度、基数和约定为准：

```bash
mkdir -p /absolute/path/new-runs/kernels
scripts/mma/wolfram.sh mma/gtb2403/build_kernels.wls \
  rw:/absolute/path/new-runs/kernels:/out -- 10 ours_M10 80 20 40 3 node
python -m smatrix_bootstrap.run gtb solve \
  --kernels /absolute/path/new-runs/kernels/ours_M10 \
  --workdir /absolute/path/new-runs/fresh_leaf --ni 46 \
  --objective F1max --reg linf-relax --precision 128 --nproc 1
```

M=10 是工具链检查，不代替论文规模结果。2309 的完整示例见原 README 和英文代码包；2403 的实际参数、三方核核对、Phase B 和精度约定见现有计划及预登记。

归档前的轻量检查能导入 Python 模块及加载 M=10 核；Docker 新容器启动发生过超时，原目录也出现相同现象。本次未修复或重测完整容器求解链，因此不声称当前机器可立即完成全部重算。
