# 主线综述工作稿：构建与 QA

- [REVIEW_ZH.pdf](REVIEW_ZH.pdf)，6页。标题、页眉及正文均明确是工作稿，Watson固定目标仍在计算中；没有新增最终结果或成功声明。
- 原始Markdown快照保存在 [sources/MAINLINE_REVIEW_ZH.md](sources/MAINLINE_REVIEW_ZH.md)。源文件科学文本未修改；图片 [images/fig9.png](images/fig9.png) 是既有原图的逐字节副本。
- 使用已验证的Pandoc/Tectonic配方：A4、10pt、18mm页边距、Droid Sans Fallback中文及DejaVu字体；单线程环境变量、nice=19与I/O idle class顺序构建。
- 仅作编译层布局：将原标题用于两行PDF标题并加工作稿标记、标题层级平移、重定位相对链接、固定既有图的位置与等比尺寸。表格保持10pt原布局，没有为压缩宽度缩小科学文字。
- Pandoc和Tectonic退出码均为0，Tectonic额外重排2次。初稿新增页眉破折号的缺字已改为受支持的中点；初稿日志及TeX保留。
- 最终日志无缺字、overfull、underfull、undefined-control或字体缺形警告。保留一条非致命labels-rerun提示；所有目录书签已逐项核对到正文页码。
- 全部6页已视觉检查；CJK、Gram与相位公式、两张数值表、第7节完整claim表及第4页fig9/图例/坐标轴均可读，无裁切或重叠。
- 所有本地链接目标存在；外部URL原样保留，未扩展联网核查。

完整来源／图片SHA、工具与字体版本、构建命令及日志SHA见 [manifest.json](manifest.json)；机器QA见 [qa.json](qa.json)。最终日志为 [compile2.log](compile2.log) 与 [REVIEW_ZH.log](REVIEW_ZH.log)。该draft目录保留，主线最终结果应另建final目录。

源Markdown SHA256：`c890e9949c5624809c328b418bb00e20331959540dd8fc113ea721e23404212c`。
PDF SHA256：`29152d8731b9ddfda394c7ab1ac9dc9dc29ee339ae5fc59567d28e49b2100a98`。
