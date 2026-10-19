# 图件复现登记表 / Figure reproduction registry

来源审查日期：2026-10-19。正文共 Fig.1–5；Supplementary Information 共 S1–S22（60 页）。下表已覆盖全部主要图和补图。状态必须根据本仓库实际输出更新，不能将“已读原图/写好脚本”记为 PASS。

- **PASS**：指定范围实际运行且输出/数值检查通过；只有明确声明的范围通过。
- **PARTIAL**：局部复现、有参数/数据/排版差异，不声称完整复现。
- **NOT_RUN**：尚未运行或未取得可审核的运行证据。
- **BLOCKED**：已确认受到数据、依赖或环境限制；不能仅因不在本轮范围就标 BLOCKED。

路径均相对仓库根目录；图件格式默认 PNG/PDF/SVG。`—` 表示未实现/未生成，不能视为已有文件。项目 Figure01–07 是教学模块编号，与原论文 Fig.1–5 不是同一编号系统。

## A. 原论文主图 / Main figures

| 原图 | 原始链接 | 名称及面板范围 | 数据 | 方法/官方函数 | 可复现范围 | 本轮状态 | R 脚本 | 输出路径 | 与原图差异 | 未复现原因/后续 |
|---|---|---|---|---|---|---|---|---|---|---|
| Fig.1a–c | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/1) | 方法输入、特征选择、LDA、θ/β | 示意流程 | LDA 方法示意 | 可绘制有署名的教学流程图 | PARTIAL（教学流程） | README 流程图 | README.md | 已绘制Mermaid教学流程，不是原图重绘 | 不属于量化数据复现 |
| Fig.2a–h | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/2) | MERFISH 聚合模拟、真值相关、不同参考 RMSE | MPOA MERFISH 12 sections，模拟 ST | `fitLDA/getBetaTheta/getCorrMtx`；参考方法；RMSE | 可另下载原始数据及基准代码续作 | NOT_RUN | — | — | 本轮 mOB 不能替代真值模拟 | 未执行大规模多方法 benchmark |
| Fig.3a–l | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/3) | 年龄扰动、混合解卷积与聚类/BayesSpace | aged/young mammary scRNA 模拟；MERFISH cortex | LDA、聚类、BayesSpace、RMSE | 需专门模拟/比较模块 | NOT_RUN | — | — | 不能用 mOB 空间图冒充 | 本轮无上述模拟与 BayesSpace |
| Fig.4a–d | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/4) | mOB 混合饼图、X7、Sox11/Nrep、Allen ISH | mOB replicate 8；原论文 255 OD genes、K=12 | `preprocess/fitLDA/getBetaTheta/vizAllTopics/vizTopic/vizGeneCounts` | 相同公开切片、260×255、K2:20与固定12已真实运行 | PARTIAL | `01b_Paper_mOB.R` | `figures/official_reproduction/paper_mOB/` | 当前min perplexity K17；K12 raw rare=0；Topic/颜色未对齐2022模型 | Fig4d ISH未复现；当前Topic7尚无原模型谱匹配 |
| Fig.4e–g | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/4) | Visium、DBiT-seq、Slide-seq 空间组成 | mouse brain、E11 embryo、cerebellum | LDA + `vizAllTopics` | 可按原始数据各自复现 | NOT_RUN | — | — | mOB 不代表这些平台的真实运行结果 | 优先完成小数据核心；未下载这些数据 |
| Fig.5a–h | [原图](https://www.nature.com/articles/s41467-022-30033-z/figures/5) | 乳腺癌免疫 Topic、GO GSEA、肿瘤模拟 | breast cancer ST + immune-excluded/infiltrated 模拟 | LDA、空间可视化、表达 logFC、GO GSEA | 需乳腺癌数据和独立模拟 | NOT_RUN | — | — | 本仓库 mOB 富集不等于癌症 Fig.5 | 未执行癌症或免疫浸润模拟 |

## B. 全部补图 / Supplementary figures

共用来源 **SI**：[Supplementary Information](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-022-30033-z/MediaObjects/41467_2022_30033_MOESM1_ESM.pdf)，PDF 页码从第一页算起。官方下载若重定向，可通过 [Europe PMC 附件包](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9055051/supplementaryFiles) 获取同一文件。下表函数仅列对应能力；原文未给专用包函数的模拟、统计和作图，不虚构官方封装。

| 原图/页码 | 原始来源 | 图件名称 | 数据 | 方法/官方函数 | 可复现性 | 状态 | 脚本 | 输出 | 差异 | 未复现原因/后续 |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 / 28 | SI | 12 个 bregmas 的 MPOA 真值及解卷积 | MERFISH 聚合 ST | LDA + 空间饼图 | 可按官方数据续作 | NOT_RUN | — | — | 无本轮对应结果 | 本轮以 mOB 为核心 |
| S2 / 29–30 | SI | 模拟 MPOA 的表达/组成与真值匹配 | 同 S1 | `getCorrMtx`、rank GSEA | 可续作 | NOT_RUN | — | — | 无真值精度结果 | 需重建真值数据 |
| S3 / 31 | SI | 更换/缺失单细胞参考后的预测 | MERFISH 与参考子集 | SPOTlight/RCTD 等比较 | 需多方法环境 | NOT_RUN | — | — | 无参考基准输出 | 本轮不运行参考 benchmark |
| S4 / 32 | SI | aged/young 模拟表达恢复 | mammary scRNA 模拟 | LDA、表达相关 | 可续作 | NOT_RUN | — | — | 无年龄扰动结果 | 需原始 scRNA 与模拟 |
| S5 / 33 | SI | cortex MERFISH 单细胞聚类 | mouse coronal MERFISH | 单细胞聚类与空间定位 | 可续作 | NOT_RUN | — | — | 非 mOB | 未执行该上游聚类 |
| S6 / 34 | SI | mOB 聚类、K 诊断、marker 与层比较 | mOB #8 | LDA、t-SNE、相关、logFC | 可部分复现；须论文预处理 | PARTIAL | `01b_Paper_mOB.R`; `04_Topic_Annotation.R` | `figures/official_reproduction/paper_mOB/`; `results/models/paper_mOB/K_diagnostics.csv` | 已运行论文K网格/marker图型；未重建t-SNE与全原面板 | 不声明S6完整复现 |
| S7 / 35–36 | SI | mOB 生物学重复稳定性 | replicates 2、5、8、12 | β Pearson、Topic 配对 | 可按官方教程续作 | PARTIAL | `06_Multi_Sample_Analysis.R` | `figures/extended_analysis/06*.png`; `results/multi_sample/` | 真实Rep8/Rep2，共同232genes固定K8独立拟合；不是原四样本分析 | Rep5/12及论文全参数对照未运行 |
| S8 / 37 | SI | mOB 不同解卷积方法相关 | mOB + MOB scRNA | θ Pearson、参考型比较 | 需匹配参考与多方法 | NOT_RUN | — | — | 未运行方法对照 | 不以 layer proxy 替代 scRNA |
| S9 / 38 | SI | 参考缺失 OEC 的影响 | mOB + 缺失 OEC reference | 参考型解卷积比较 | 可续作 | NOT_RUN | — | — | 未重建缺失参考实验 | 需参考环境与原数据 |
| S10 / 39 | SI | cortex 参考导致 VLMC 错分 | mOB + cortex scRNA | 参考型解卷积与 markers | 可续作 | NOT_RUN | — | — | 未运行错配参考实验 | 同上 |
| S11 / 40–41 | SI | 分辨率变化时的真值恢复 | MERFISH 多分辨率聚合 | LDA、RMSE、比例相关 | 可续作 | NOT_RUN | — | — | 无分辨率 benchmark | 需多尺度模拟 |
| S12 / 42–44 | SI | Visium 单 Topic 与 markers | Visium mouse brain | `vizTopic/vizGeneCounts` | 可按 Visium 教程续作 | NOT_RUN | — | — | mOB 不替代 Visium | 未执行该数据 |
| S13 / 45 | SI | DBiT-seq embryo 的 Topics | DBiT-seq E11 lower tail | LDA、空间图、marker 表达 | 可续作 | NOT_RUN | — | — | 无该数据结果 | 未执行该数据 |
| S14 / 46 | SI | Purkinje/Bergmann glia 与 RCTD | Slide-seq cerebellum | STdeconvolve/RCTD 比较 | 需 Slide-seq 与参考 | NOT_RUN | — | — | 无细胞类型比较 | 未执行该数据 |
| S15 / 47 | SI | 运行时间与内存 | Visium 子采样 | microbenchmark、profmem | 可另作基准 | NOT_RUN | `00_Start.R`（仅运行日志） | `environment/run_metrics.csv`（非论文基准） | 本机耗时/R 堆内存不同于论文系统/总分配 | 运行日志不称为 S15 复现 |
| S16 / 48–49 | SI | 乳腺癌 ST 的解卷积 | breast cancer ST | LDA、相关、空间图 | 可续作 | NOT_RUN | — | — | 无癌症结果 | 未执行癌症数据 |
| S17 / 50 | SI | 癌症 X3/X13 表达与定位 | 同 S16 | markers、空间图、logFC | 可续作 | NOT_RUN | — | — | 无上述 Topic | 同上 |
| S18 / 51 | SI | 癌症 X15 免疫表达 | 同 S16 | marker、表达、GO | 可续作 | NOT_RUN | — | — | mOB 富集不能替代 | 同上 |
| S19 / 52 | SI | spot 数量对精度影响 | 模拟 ST 子采样 | LDA、真值 RMSE | 可另作小规模模拟 | NOT_RUN | — | — | 无 spot 数量曲线 | 未执行该基准 |
| S20 / 53–54 | SI | 均匀混合比例的失败示例 | luminal/macrophage 模拟 | LDA、表达/比例相关 | 可按官方失败教程续作 | NOT_RUN | — | — | 无模拟失败输出 | 文档解释不能当运行结果 |
| S21 / 55 | SI | 难以分离细微 neuronal subtypes | MERFISH 76 types/subtypes | K=76 与真值相关 | 可续作但需资源估计 | NOT_RUN | — | — | 未执行高 K | 当前小数据/内存策略不含此项 |
| S22 / 56 | SI | 入模 gene 数对精度影响 | 模拟 ST 的 OD gene 子集 | 随机子集、LDA、RMSE | 可另作小规模模拟 | NOT_RUN | — | — | 无基因数量敏感性曲线 | 不将本次232基因拟合称为该基准 |

补充 Data 1 是癌症 Topic 15 的显著 GO terms；不是 mOB GO 结果，也没有作为本项目富集输出使用。

## C. 官方教程与本仓库教学图 / Tutorial figures

下表状态已结合实际脚本退出状态、数值合同、图件三格式和视觉审查更新；PASS只适用于声明范围。

| 本项目编号/原始位置 | 原图链接 | 名称 | 数据 | 方法/官方函数 | 可复现性 | 状态 | 脚本 | 目标输出 stem | 差异 | 未完成项 |
|---|---|---|---|---|---|---|---|---|---|---|
| Figure01 / Getting Started K 图 | [教程](https://jef.works/STdeconvolve/getting_started.html) | K、perplexity、rare Topics、α | mOB quick corpus | `fitLDA/perplexityPlot` | 可直接运行 | PASS（指定计算/图型） | `02_Topic_Number_Selection.R` | `figures/official_reproduction/Figure01_Official_K_Selection` | 统一字体/尺寸；数值以本轮表为准 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure01 / 本项目扩展 | 同上为依据 | 不同种子匹配稳定性 | 同一 mOB、固定 K | β Pearson + Hungarian one-to-one | 小规模扩展 | PASS（指定计算/图型） | `02_Topic_Number_Selection.R` | `figures/extended_analysis/Figure01_Seed_Stability` | 四seed配对均值βr .758–1；非生物学重复 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure02 / Getting Started 空间图 | [教程](https://jef.works/STdeconvolve/getting_started.html) | 全 Topic 饼图 | mOB | `vizAllTopics` | 可直接运行 | PASS（指定计算/图型） | `03_Spatial_Visualization.R` | `figures/official_reproduction/Figure02_All_Topics` | 模型版本/配色/尺寸可有差异 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure02 / Getting Started 单 Topic | 同上 | 各 Topic 空间比例 | mOB | `vizTopic` | 可直接运行 | PASS（指定计算/图型） | `03_Spatial_Visualization.R` | `figures/official_reproduction/Figure02_Topic_Maps` | 本项目采用统一 panel 网格 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure03 / 本项目扩展 | 同上为依据 | Topic × gene 热图 | mOB β | `topGenes`、ggplot2 heatmap | 可直接运行 | PASS（指定计算/图型） | `03_Spatial_Visualization.R` | `figures/extended_analysis/Figure03_Topic_Gene_Heatmap` | 图为本项目设计，非差异表达检验 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure03 / 官方 marker 展示方法 | 同上 | 测得 marker counts 空间图 | mOB counts | `vizGeneCounts` | 可直接运行 | PASS（指定计算/图型） | `03_Spatial_Visualization.R` | `figures/official_reproduction/Figure03_Marker_Count_Maps` | 本次 top β genes 不等于论文 X7 marker | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure04 / Annotation 教程 | [教程](https://jef.works/STdeconvolve/celltype_annotation.html) | β/θ 与组织层 proxy 相关 | mOB + layer labels | `getCorrMtx/correlationPlot` | 可直接运行 | PASS（指定计算/图型） | `04_Topic_Annotation.R` | `figures/official_reproduction/Figure04_*` | proxy 来自同一 ST，不是独立 scRNA 真值 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure04 / Annotation GSEA | 同上 | layer markers 富集与候选注释 | mOB β + layer proxy sets | `annotateCellTypesGSEA` | 可运行，依赖版本字段需核验 | PARTIAL（官方方法+显式参数调整） | `04_Topic_Annotation.R` | `results/annotation/` | gseaParam=.1来自固定vignette；minSize15/串行为项目调整；同切片代理 | 独立scRNA参考注释NOT_RUN |
| Figure05 / 本项目扩展 | [原论文 Fig.5e 为 GSEA 先例](https://www.nature.com/articles/s41467-022-30033-z/figures/5) | GO/KEGG/Reactome/GSEA | mOB β + 正式基因集 | 下游工具，不虚构官方 API | 需数据库和ID映射 | PARTIAL | `05_Functional_Enrichment.R` | `figures/extended_analysis/05A_ORA_dotplot` 等 | GO/KEGG已跑；GO ORA6项、KEGG ORA0项；GSEA9条NA | Reactome官方403，BLOCKED；保留全部未显著/NA结果 |
| Figure06 / Combining 教程 | [教程](https://jef.works/STdeconvolve/combining_datasets.html) | 真正跨样本 Topic 对应 | 真实独立 mOB replicates | `getCorrMtx/lsatPairs` | 需其他真实公开样本 | PASS（指定计算/图型） | `06_Multi_Sample_Analysis.R` | `figures/extended_analysis/06A_cross_sample_similarity` 等 | 真实260/279spots，232共同genes；配对r .472–.993；n=2仅描述 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |
| Figure07 / 本项目扩展 | [官方标准方法范围](https://jef.works/STdeconvolve/) | Topic 空间邻域/共定位 | mOB θ + 坐标 | 空间权重、置换统计；非官方内置 | 可设计小规模检验 | PASS（指定计算/图型） | `07_Spatial_Neighborhood.R` | `figures/extended_analysis/07A_spatial_knn_graph` 等 | 对称6NN、999次整条组成联合标签置换；双侧p与明确BH；不控制层/梯度 | 已核对真实表、图件三格式和矩阵合同；科学边界见解读手册 |

## D. 如何更新登记状态

1. 运行编号脚本，保存 stdout/stderr、参数、版本与 elapsed time。
2. 检查图文件三格式均能打开、字体/图例无截断；核对输入规模、K、θ/β 行和、spot ID 顺序。
3. 将本表 C 对应行更新为 PASS 或 PARTIAL，并写明实际路径/差异；若只有局部结果，不更新 A/B 为完整 PASS。
4. 将未运行项目保留 NOT_RUN；数据/依赖确实失败后才记录 BLOCKED 与原始错误。

本表不借用原论文的图作为本仓库实际结果。引用原图用于学习时应注明原作者及 CC BY 4.0/第三方图像的各自授权；本仓库生成图依本次真实模型解释。

27幅实际输出的逐文件索引见 [figure_manifest.csv](../figures/figure_manifest.csv)，逐图十项解读见 [04_图件解读指南.md](04_图件解读指南.md)。本表全部状态已按本次真实运行更新。
