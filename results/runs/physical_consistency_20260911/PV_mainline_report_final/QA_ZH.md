# PV主线正式报告交付检查

成品：[REPORT_ZH.pdf](REPORT_ZH.pdf)，18页。首页明确“有限链已完成；论文定量核心仍有major差异”，并区分有限IR→UV信号与F五组完成、完整论文定量成功、描述性90°读数、严格离散峰和极点。

科学正文冻结SHA为`4ba7d903c94155799bc4e29f208aee1d5548ed08cd3cb2fc9a77b7ff92c487f8`。原文件与只读快照仍逐字节一致；正文经清单所列链接重定位及数学记号等价转写后可完整重建，未删减推导。18项A1–F3台账全部保留。八份指定PDF图件原字节保持；Fig.8使用PV_E1_regions_complete/regions.pdf。

全部18页已视觉检查，覆盖正文、公式、长表、结论和八份图。无缺字或overfull溢出警告；保留1处非溢出的underfull行间距提示，已目视确认无截断。27个PDF书签、23个相对文件链接及正文15个本地链接均有效；另有1个外部作者README链接。

构建仅通过内联Pandoc/TeX命令执行，单线程环境、nice19；未创建或恢复活动build.py。命令见logs/commands_1.json、logs/commands_2.json和logs/BUILD_COMMANDS.txt。首次PDF元数据宏加载顺序问题已仅在排版头文件修正，原失败日志保留，最终编译为logs/compile_2.log。

技术来源读取自[PV_FINAL_VERIFICATION.json](sources/PV_FINAL_VERIFICATION.json)：107测试对应的21个当前模块SHA一致，模块行/字节限制与3个测试文件行数限制通过；27个CLI报告完成，142个已声明输入哈希保持，29份IR上游报告及C另已核对。该校验明确记录旧regions汇总接口不产source.json，原记录保留null；未补造producer，也不把该接口范围当作物理major问题。排版端未重跑测试或物理计算。

PDF SHA-256：`1a0361352682665879605ea42bacf2d83e68407c94a201bab90eea6a6f9e98c4`。

完整交付仅指本文及已声明有限链的证据齐备；科学结论以冻结正文的最终判定和边界为准。
