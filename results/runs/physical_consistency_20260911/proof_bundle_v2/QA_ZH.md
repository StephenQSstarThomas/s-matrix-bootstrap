# 证明与适用范围伴册 v2：构建与 QA

- 输出：[PROOFS_ZH.pdf](PROOFS_ZH.pdf)，33页，包含原 manifest 的12篇证明及指定新增3篇。
- 标题／副标题明确是证明与适用范围伴册，不是完整物理复现成功报告；前言注明 Watson 正式结果尚未完成。
- 原12篇来源及另外2篇新增来源均未改动。依明确补充指示，仅将 Watson 文中第3节标题与秩二例句的两处“Gram 正定”校正为“Gram 半正定”；模型、方程、结论及其它文本不变。
- 当前15篇原文快照在 `sources/`。校正前原文、SHA、已构建版本及变更记录在 [terminology_history/correction.json](terminology_history/correction.json)。所有当前来源与快照均匹配 manifest；原12篇与旧 manifest 的SHA全部相同。
- 旧28页 PDF保持不变。v2沿用 A4、10pt、18mm页边距、Droid Sans Fallback中文字体及DejaVu西文字体。
- 其余处理仅涉及新增封面／范围说明、换页、局部链接路径和中文斜体字体；原12篇沿用旧版数学／路径排版。首次链接重复换算已修正，初次TeX与日志保留。
- Pandoc与Tectonic退出码均为0，最终执行额外2次TeX重排。没有缺字、overfull、undefined-control或字体缺形警告。25个局部链接目标均存在，15个章节书签页码全部对应正文。
- 视觉检查第1、2、29、30、31、32、33页，覆盖目录、新Watson公式与半正定术语、实际SOC/density Schur约化障碍、峰存在性表格。未见裁切、重叠或缺字；最终页渲染与已检查PNG一致。
- 构建使用nice=19、I/O idle class及单线程环境变量，逐步执行。未运行科学优化或修改生产源码。

来源SHA、工具与字体SHA、工具版本、命令及日志SHA见 [manifest.json](manifest.json)；机器QA见 [qa.json](qa.json)；最终日志见 [compile5.log](compile5.log) 和 [PROOFS_ZH.log](PROOFS_ZH.log)。

保留的非致命TeX提示：`Package amsmath Warning: Foreign command \over;`；`Package amsmath Warning: Foreign command \atopwithdelims;`；`LaTeX Warning: Label(s) may have changed. Rerun to get cross-references right.`。已强制重排并独立核对章节书签／正文页码，不隐藏这些提示。

最终PDF SHA256：`6280627f6191a7ff33cdaa93d097bb0dd8356fc1737b39cbbe145918bf51e503`。
旧28页PDF SHA256：`b458ca0092a5960d9d0295dda4fc884f600e5fbff817a568ffe2be578f849b09`。
