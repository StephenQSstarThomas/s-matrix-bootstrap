# 专家审阅补充包 · 2026-09-10

本包内嵌于 He–Kruczenski 2309.12402v3 中文专家审阅报告。
主 PDF 为30页；全部数值是仓库已有结果的快照，本次未追加 bootstrap 求解。

## 阅读与使用

1. 先读主 PDF 的第3页状态表、第4–11页方法、第12–26页原图对比与缺口。
2. 当前模型以 results/runs/mainline_alignment_20260909/PAPER_MAINLINE.json 为准。
   PV/midpoint、sampled/free T0、combined八维L2、hard-midpoint FESR、四个raw独立误差。
3. src/smatrix_bootstrap 含十个活跃核心模块；tests 含三个测试文件。
   Python项目安装及统一入口参见 README.md 与 REPRODUCTION_GUIDE_ZH.md。
4. 本包保留主要汇总、曲线数据、代表解及失败记录；JSON中的原始绝对路径忠实保留。
   在其他机器使用时应按本地解包根目录解析相应仓库相对路径。
5. 本包不是完整 results/runs 的镜像；大型H/current矩阵、所有中间日志和历史证明producer仍在仓库原路径。
   PDF第30页给出主要索引；不凭本包缺少某个历史输入推断该记录不存在。
6. report_source/report.tex 与 assets/ 可用 Tectonic 重建正文；中文字体为 Droid Sans Fallback，西文为 DejaVu。
   该重建不自动恢复PDF附件。本包已经提供重排后的矢量图，不修改任何原数值。
7. 原图来自 references/2309.12402v3-source/；原作者完整振幅系数未恢复。
   figure*.csv 为原图读数，插值只是比较／显示，不是优化目标。

## 状态边界

A–F全部claims尚未闭合：B3整体区域精度未完成；E的强非对称、三点相移稳健性及部分S0/弹性性质存在差异；F仅4/5配置完成。
(50,12)末候选没有通过联合原式，也没有不可行性证书。F固定L10的三M齐全只表示数据齐全。
历史separate/clipped/Watsonian计算与当前模型分开解读，不能拼成当前全通过。
本包源代码是项目实现，不是未公开的2309原作者程序。有限幺正盘与支持界不是连续域认证。

## 文件索引

清单见 FILES.txt。报告只重排已有数据；内／外区域以保存顶点的凸包边界作图，不按存储顺序连线。
