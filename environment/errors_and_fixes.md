# 运行问题与处理 / Execution issues

此表记录实际遇到的问题；未通过的尝试不记为 PASS。成功运行后的完整日志在 logs/；process_metrics.csv 保留了04的失败退出记录及之后成功记录。

| 现象 | 原因与处理 | 状态 |
|---|---|---|
| 官方 git clone 连接重置 | 改用同一官方仓库 codeload ZIP，固定源码提交；没有更换方法或数据 | PASS：源码安装/包加载 |
| raw counts 的全部262spots没有对应annot | 官方annot仅覆盖260个过滤后spots；pos与counts排序也不同。最终保留现有ID，过滤后按名字匹配；禁止按位置重命名 | PASS：真实输入与回归检查 |
| correlationPlot 离散值传入连续轴 | reshape2把纯数字Topic名转数值，但文字层用字符；只在显示矩阵副本添加Topic前缀，不改相关矩阵 | PASS：两相关图均导出 |
| fgsea worker找不到fgseaMultilevelCpp | 下游默认并行的包库/worker加载冲突；明确使用SerialParam与nproc=1，不改算法 | PASS：层代理GSEA及05串行运行 |
| Reactome官方映射下载403 | current和官方版本下载受访问限制；保留BLOCKED。官方Zenodo链接提供整包约4.9GB，本轮没有为单一小示例下载整库 | BLOCKED：不伪造或用人类库替代 |
| fgsea有9条无法估计p/NES | 正负排序统计分布不平衡，保留NA与警告，不改成0或降阈值造阳性 | PARTIAL：8935/8944条有p值 |
| Windows继承LC_*=C.UTF-8启动警告 | 当前Windows locale不接受该名称；任务已实际跑通。没有修改系统区域设置 | 警告保留 |
| 图中文字乘号/连接符显示异常 | 最终绘图标签使用ASCII * 与 -，仅重新绘图；论文LDA缓存复用，不重算已通过模型 | PASS：最终PDF文字审查 |
| 自动PDF字级审查将Cairo Tf=1当作1pt | Cairo同时通过文本变换矩阵放大；用读取实际变换后字级的PyMuPDF复核，27图最小8pt | PASS：实际字级；原简化审查工具局限已记录 |

安装只向本项目独立库新增8个包，未升级R或改写其他项目依赖。调试中曾出现检查条件过严、路径/可执行文件定位与工具读取目录错误，均在运行/文件审查时修正；不把这些尝试当作模型或数据验证。
