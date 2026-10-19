# 依赖版本 / Dependency versions

实际执行：Windows，R 4.6.1（UCRT），Bioconductor 3.23，串行。官方源码固定提交 `13b1953263a5ce34e2471c4defebedbda44c7d08`。

| 包 | 本次版本 |
|---|---|
| STdeconvolve | 1.3.2 |
| topicmodels | 0.2-17 |
| fgsea | 1.38.0 |
| BiocParallel | 1.46.0 |
| Matrix | 1.7-5 |
| mgcv | 1.9-4 |
| ggplot2 | 4.0.3 |
| scatterpie | 0.2.6 |
| slam | 0.1-56 |
| clue | 0.3-68 |
| patchwork | 1.3.2 |
| ragg | 1.5.2 |
| svglite | 2.2.2 |
| AnnotationDbi | 1.74.0 |
| org.Mm.eg.db | 3.23.0 |
| GO.db | 3.23.1 |
| digest | 0.6.39 |

完整已安装包快照见 installed_packages.csv，当前载入会话见 sessionInfo.txt。安装器固定STdeconvolve源码，其他依赖按可得仓库解析；本快照不是自动恢复的renv锁文件。其他机器/未来版本需按快照审查，不能宣称跨所有版本逐数值重演已验证。

独立库新增：NLP、tm、fastmatch、topicmodels、fgsea、STdeconvolve、org.Mm.eg.db、GO.db。未修改其他版本的R或用户原包库。
