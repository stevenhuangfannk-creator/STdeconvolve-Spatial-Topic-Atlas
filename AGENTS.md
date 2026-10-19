# 协作规范 / Research collaboration

本仓库是独立的 STdeconvolve 公开方法学示例。只能操作本仓库；不得读取或改写其他科研项目的正式结果。本阶段禁止使用种植体周围炎私人数据。

- 先读 README、来源登记和实际代码。明确输入、输出、参数与验收条件，再作最小改动。
- R 脚本按顶部说明 → 集中可调参数 → 依赖 → 数据检查 → 分析/绘图 → 预览与导出排列；中文注释，标明输入矩阵方向和单位。
- 00_Start.R 统一参数；独立包库 environment/R-library；不自动安装包，不升级全局 R 或包库。Windows 16 GB 默认单核。
- 保留原始数据及已确认结果；仅按官方方法处理真实数据。不得为图形相似修改注释、筛除细胞或伪造结果。
- Topic 是潜在表达成分，不能仅凭少量 marker 称为新细胞类型。theta 是模型比例，不是绝对细胞数。空间邻近不证明通讯。
- 重要 GitHub 文档中文为主，保留英文标题和术语；标识符、命令、路径与机器字段不翻译。
- 修复数据对齐或模型问题后运行 tests/validate_outputs.R；检查真实图件，更新复现登记、运行日志和报告。
- PASS 只表示真实执行并验证通过；代码已写但未运行用 NOT_RUN；依赖/数据受阻用 BLOCKED；部分复现用 PARTIAL。
- 升降级 R、大规模依赖调整、更换官方数据、改用其他算法、预计超过 16 GB 或修改其他仓库前，必须先报告并等待用户决定。
- 新的真实科研分析由用户在 VS Code 中执行和审查；本仓库已授权的官方小数据复现除外。

English summary: keep edits surgical; use real public observations and official functions; distinguish implemented, executed and validated states; do not touch private research projects.
