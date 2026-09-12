# 历史IR代表的中心身份复核

结论：蓝、橙、绿三份最终IR代表都不是其已保存数值中心的同一完整C。`center_best.npz`按原坐标恢复后逐项等于`native_candidate.npz`；随后全部3876个C都被执行一次float64的(1−10⁻⁸)收缩。`candidate.npz`、`coefficients.json`和最终selection一致保存这个收缩后的点，而run与selection仍标记representative_center_converged=True。现有中心日志不能自动转移到这个改变后的C。

本审计没有修改任何历史输入、源码、选择、状态或图，也没有重新优化、计算Hessian或检查相移。历史来源保持原位；补充的机器可读结果见[IR_CENTER_IDENTITY_REVIEW.json](IR_CENTER_IDENTITY_REVIEW.json)。

## 实际身份链

最终来源由SEQUENCE及Fig.5保存的profile元数据共同定位：blue→C1_blue，orange→C1_orange，green→C1_green_balanced。三者都是M50/L10、subtracted数值坐标、actual-rho打包、同一A2_M50_L10算子，固定xref=0.07332139057293466。ε分别为.006/.004/.002。

使用现存转换函数重放center_best.z的subtracted→unsubtracted转换、float64舍入和ρ2非对角存储因子；所涉及五个函数的AST与各run保存producer均完全一致。没有执行旧源码快照。每份center_best的迭代数与run末次迭代一致、μ与末条history一致、has_newton_gradient=True；转换结果与native_candidate完全逐项相等。另行核对center_first与其初始coefficients文件也相等，且与best不同，未把初始点误认作最终点。

| 颜色 | center_best μ | 迭代 | 保存decrement | 改变的C数 | 最大绝对C差 | 相对L2变化 |
|---|---:|---:|---:|---:|---:|---:|
| blue | 1e-06 | 100 | 5.5240765e-15 | 3876/3876 | 0.003197618178 | 1.000000005e-08 |
| orange | 1e-06 | 99 | 2.0777012e-19 | 3876/3876 | 0.003182942106 | 1.000000005e-08 |
| green | 1e-07 | 111 | 1.3944738e-14 | 3876/3876 | 0.003340187424 | 1.000000005e-08 |

三份均满足精确的程序结果身份 `candidate == np.asarray((1-1e-8)*native_candidate, float64)`。这里指浮点乘法的保存结果，不把它写成无舍入的实数等式。变化的绝对值约.003来自原C系数较大；相对变化虽小，仍不是同一个幅度，也不是中心证明。

恢复元数据三者相同：candidate_fraction=1、radial_scale=1、incumbent_retained=False，final_contraction=.99999999、repaired=True。结合原producer控制流，这对应原candidate先通过audit，再进入finish的无条件收缩；不是失败后的segment或seed fallback。

## 固定截面也发生小幅移动

| 颜色 | native x−xref | delivered x−xref | delivered x−native x |
|---|---:|---:|---:|
| blue | -1.849381759e-12 | -7.315964778e-10 | -7.297470961e-10 |
| orange | 6.068784364e-12 | -7.277621561e-10 | -7.338309405e-10 |
| green | -2.894906537e-14 | -7.238797062e-10 | -7.238507571e-10 |

以上坐标来自当时保存的完整H·C值。收缩及逐系数舍入将交付x移离声明截面约7.2–7.3×10⁻¹⁰；不能借“量小”把恢复前固定截面的中心记录当作交付点的中心核验。该表也没有把浮点坐标近似当成严格等式证书。

## 标签与仍成立的证据范围

`recover_candidate_segment.finish`即使原候选通过也执行收缩；其历史coherent_boundary在返回barrier的solver记录后，只追加primal_recovery，并没有访问或清除representative_center_converged。最终selection依赖这个旧标签又写入True。这里需要收窄的是**中心身份声明**：保存的微小decrement属于恢复前数值坐标中的中心记录，当前没有在收缩后完整C上重新居中或验证同一Φ的中心条件。不能声称已证明它不是任何参数下的中心，但也不能称它已获本记录的中心证据。

三份历史finite primal、finite support和analytic primal都报告通过；解析核验所记录的coefficients输入哈希与现在的交付文件一致。因此这些证据仍对应收缩后的真实振幅，不因中心标签过强而自动作废。本次读取并认证其输入身份，没有重新计算那些证明。几何选取和后续图仍可说明使用了这三个固定、可行且有有限支持界的振幅；“经过数值中心后作微小收缩的可行支持点”比“同一个已验收敛中心”更符合保存证据。

最终C2_blue/orange/green_selection的C逐项等于各candidate；Fig.5的三项profile正是读取这些selection系数。selection中的coefficients_mixed=False可以表示selection阶段没有再混合各来源，不能据此否认此前生产恢复中有零幅度权重。三份完整交付C的身份没有在support→selection→Fig.5链中再次变化。

新实现若任何恢复改变了完整C，应取消中心资格，除非对改变后的点另有同一问题的中心核验。当前审计不改写旧文件，也不以无条件收缩的旧点代替新零T0切片的代表。

## 输入指纹

转换源码指纹记录的是本次重放时的版本；具体数学函数另与旧producer作AST比对。其余均为保留的历史输入。

| 输入 | SHA-256 |
|---|---|
| `src/smatrix_bootstrap/quotient.py` | `32b669bec9e0b9af1c21bf464f2806411cdca87fffa6ff0593d7bb77dd3428e9` |
| `src/smatrix_bootstrap/kernels.py` | `790ebcdba291c80f17c4909a31ddf77ad6373748c2af41c7ce1a0c28b309172a` |
| `src/smatrix_bootstrap/operators.py` | `74cb905a1e9e881e3f3e5f7843d44e7e096e8c57e672f8f7abdbae336b78fc4e` |
| `results/runs/sequential_reproduction_20260910/C1_blue/report.json` | `5f27a666c54aacdfe2f2d6ce3bdb78a7eca30a74b5622cfe2cae80330b588adb` |
| `results/runs/sequential_reproduction_20260910/C1_blue/source.json` | `ab8a05d9391719f347655cee9619b51b07181cd6a1a3d67db6f405257d6a1e9f` |
| `results/runs/sequential_reproduction_20260910/C1_blue/native_candidate.npz` | `550b642fd3232c88affa6901614be8746921761be4bbab334529d95d131e2c26` |
| `results/runs/sequential_reproduction_20260910/C1_blue/candidate.npz` | `94eb3be2454debbda1f39d9e30500e0e41e4f209f470d533824c7788c822b26e` |
| `results/runs/sequential_reproduction_20260910/C1_blue/center_best.npz` | `e86833a7c5d98bf9f96b148bb04f5fc69b453ce164eb71259f3fc84acc4c79c9` |
| `results/runs/sequential_reproduction_20260910/C1_blue/center_first.npz` | `9d7557bd086fedecd960da6b54299725eec0433c426030fe16ff329c87e8a4a7` |
| `results/runs/sequential_reproduction_20260910/C1_blue/center_first_coefficients.json` | `d89e6aa340391de7a4ed2976c8b6b15858b2cd21507771c1e5a728cf704c70f8` |
| `results/runs/sequential_reproduction_20260910/C1_blue/coefficients.json` | `a0abe00d61d0d325a5e54ce5328065039e5f7a5c175d3e789d9da6b2c306e887` |
| `results/runs/sequential_reproduction_20260910/C1_blue/progress.jsonl` | `e0465ca052f31f175a2d56c4514a6fb5ad1748266cd3b5dff2b9ac05c1d6c446` |
| `results/runs/sequential_reproduction_20260910/C2_blue_selection/coefficients.json` | `a0abe00d61d0d325a5e54ce5328065039e5f7a5c175d3e789d9da6b2c306e887` |
| `results/runs/sequential_reproduction_20260910/C2_blue_selection/selection.json` | `27aa497364a1cdd05086e99c08b92235e9fa878abc160ca8d3d17e6e5e5b0a12` |
| `results/runs/sequential_reproduction_20260910/C1_blue_analytic/report.json` | `91ddee8220190c11d0c2061f65b069623454a873948de95c58170ecda25265f9` |
| `results/runs/sequential_reproduction_20260910/C1_orange/report.json` | `4f067d4a20214187cc193482805b4ad8c4a669fddf7d226840b4a8bb00a2a2ab` |
| `results/runs/sequential_reproduction_20260910/C1_orange/source.json` | `f7759800974f55ed86645be51a88c60abd6a0a124201fcb4e96bf5da65fb2252` |
| `results/runs/sequential_reproduction_20260910/C1_orange/native_candidate.npz` | `0ce657195658aada1e2f6b810a1d3c2ea08256f5cf00f12a72ffb5b0dc7c7c87` |
| `results/runs/sequential_reproduction_20260910/C1_orange/candidate.npz` | `4ec38812d345e965fda7fec52bdf11dc884088e09bdbc69c7f70f3dcdf006026` |
| `results/runs/sequential_reproduction_20260910/C1_orange/center_best.npz` | `61be92a39e91cf02f71c45a72dfedf37cd2bb2b2af8859d08bec8b0837c7b712` |
| `results/runs/sequential_reproduction_20260910/C1_orange/center_first.npz` | `de098bd87d8b1c25974997d8fd2e97a75d0b52d5b40eea575ab9ceefc7d0e155` |
| `results/runs/sequential_reproduction_20260910/C1_orange/center_first_coefficients.json` | `6e5057ec318808c87c6bdd13adcd0e1325e3564edb6bc5aa5d62113fe0820465` |
| `results/runs/sequential_reproduction_20260910/C1_orange/coefficients.json` | `ded38d4430c2c97f36508a20abdc348207e21dfe91f233c86362a78b8ed91c02` |
| `results/runs/sequential_reproduction_20260910/C1_orange/progress.jsonl` | `2adce5d1fc1e2439d9fa8b3da85b5d764130d6c3dbdd5b9c45da87b7595c514e` |
| `results/runs/sequential_reproduction_20260910/C2_orange_selection/coefficients.json` | `ded38d4430c2c97f36508a20abdc348207e21dfe91f233c86362a78b8ed91c02` |
| `results/runs/sequential_reproduction_20260910/C2_orange_selection/selection.json` | `44922238ab6d8da8705e170fd50d9aef2e5d9e67acc11eeef75099331e9fdecd` |
| `results/runs/sequential_reproduction_20260910/C1_orange_analytic/report.json` | `c9a31b4e744db01795aacba5964ddb42adb5778f62b1b6ce3deae60bf4467f85` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/report.json` | `38d8ecc07a6dab00192a52b1eba706cf0ec4f20a7225157db69b228c19385b4e` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/source.json` | `a2d87ec2eaea9e40ea713e40132abf3709c83ac5391e1b52daf65303dacdff50` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/native_candidate.npz` | `41d786389e346104e7fefd8d2075646fef4350cfdcf2229284d56cc9ba9c179f` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/candidate.npz` | `8158b68a8dee049625abd53fb3bf8b0685b5b8da41d79634a760e2fec1823134` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/center_best.npz` | `5976056ba480a5cef379889cbcc487590471f844703f35005aaf2995cfc9ccb4` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/center_first.npz` | `6470c8b945de4fb76b207a4a6703f06b535d00cc6e96b906be5cbf39861779e5` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/center_first_coefficients.json` | `f1d246d44719db17749616ad39d1f00dacf0e0e182b9eb2bfeb59d5c7ab4c8b2` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/coefficients.json` | `d2de43c36f82ef74ff17b506d60111e6ecef36e0172bfa12e2c25a6759245cc2` |
| `results/runs/sequential_reproduction_20260910/C1_green_balanced/progress.jsonl` | `814f84ac697f7c8c14fa1991a2166dd908f4df308e6fc6dcef6354ac3fe2159e` |
| `results/runs/sequential_reproduction_20260910/C2_green_selection/coefficients.json` | `d2de43c36f82ef74ff17b506d60111e6ecef36e0172bfa12e2c25a6759245cc2` |
| `results/runs/sequential_reproduction_20260910/C2_green_selection/selection.json` | `ca3ec7b8b1c07fde9981a78bfb14124e334b62b0716ee37b5ca9a07daf739692` |
| `results/runs/sequential_reproduction_20260910/C1_green_analytic/report.json` | `b460247c3435a2204aa3b34de302d4ce7d3f8203e4b678da7fac2fc51eda8023` |
| `results/runs/sequential_reproduction_20260910/C1_fig5/profiles.json` | `ae996983f07f66d3ee67b8f615135b9594c467ab58d1b37fb579c8207a72dfa9` |
| `results/runs/sequential_reproduction_20260910/SEQUENCE.json` | `01e1b9739e8f045db1adc9f54a55cdd3d67906f8e72b343876e4df3f0552d7ba` |
