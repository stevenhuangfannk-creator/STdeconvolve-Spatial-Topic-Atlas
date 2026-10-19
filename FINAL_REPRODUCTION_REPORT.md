# 最终复现报告 / Final reproduction report

审查与执行日期：2026-10-19（Asia/Singapore）。总体状态 **PARTIAL**：官方mOB核心、论文参数支线和可运行扩展已真实执行，但不声称整篇论文或所有数据库完全复现。未使用PI私人数据，未更改其他科研项目。

## 1. 官方资料审查

**PASS（审查范围）**：核对原始论文DOI 10.1038/s41467-022-30033-z、Methods/Results/Discussion、60页SI、S1–S22全部图注、官方Getting Started/Annotation/Combining/Failure教程与函数源码。论文信息、各方法对照及精确链接见docs/01_论文精读.md。主图Fig1–5、补图S1–S22与教程图均登记在docs/Figure_Reproduction_Registry.md。

## 2. 软件与环境

**PASS（本机环境）**：已有Windows R4.6.1（UCRT）、Bioconductor3.23。STdeconvolve1.3.2固定官方package提交13b1953263a5ce34e2471c4defebedbda44c7d08。向项目独立R-library新增8包，没有升级R，没有改写原有Seurat/Bioconductor库。依赖版本、会话和完整安装包快照在environment/。未来安装器的其他依赖不完全锁版，跨机器逐数值重演尚未验证，不能把版本快照称为已验证renv恢复。

## 3. 数据取得与合同

- **PASS**：官方包内mOB原始15928genes×262spots；坐标262，组织层annot260。pos顺序不同，按真实Spot ID对齐。标签仅为五个组织层，不是细胞组成真值。
- **PASS**：快速教程清理后14828×260，OD232×260；论文参数支线260×255，与论文规模一致。
- **PASS**：官方Rep2 TSV下载成功，跨样本分析保留279spots，与Rep8的260spots共享232genes。没有把同切片拆分冒充生物重复。
- **PASS**：GO与KEGG小鼠来源/版本、映射、背景与SHA256已记录。原始数据/模型实际保留在本地，Git排除大数据、二进制模型和包库；data/data_manifest.csv保留8个输入/关键对象校验值。数据库和第三方数据许可保留来源条款。

## 4. 真实运行情况

| 模块 | 状态 | 实际证据与范围 |
|---|---|---|
| 00环境/参数 | PASS | 依赖加载、参数快照、单worker/串行下游 |
| 01官方快速流程 | PASS | 全K2–9 LDA，seed0，选K8；θ/β、模型与对齐索引已保存 |
| 02模型诊断/种子 | PASS（探索性扩展） | seed0/1/2/3，β Hungarian匹配；跨两类解平均βr约.758、θr约.790，不称完全稳定 |
| 03空间/基因图 | PASS | 官方pie、全部8Topic raw图、7个不同top1 gene图及β热图 |
| 04层代理注释 | PARTIAL（官方方法+参数调整） | 相关与官方代理GSEA已实际运行；gseaParam=.1来自固定vignette，minSize15/串行为项目显式设置；独立scRNA注释NOT_RUN |
| 01b论文mOB | PARTIAL（复现图型） | 原文预处理260×255；K2–20已跑，固定K12；实际min perplexity K17、K12 raw rare=0；Fig4a–c图型，不含Allen ISH与原Topic谱匹配 |
| 05GO/KEGG ORA | PASS（扩展） | 219/232可靠映射；GO背景207/1046sets，KEGG背景115/72sets；GO ORA6项q<.05均Topic4，KEGG ORA0项 |
| 05ranked GSEA | PARTIAL | 8944条中9条p/NES为NA，警告完整保留；没有将NA改为显著/不显著的确定结论 |
| 05Reactome | BLOCKED | 官方current/版本下载403；没有使用人类基因库替代，没有下载约4.9GB官方整包来完成一个小模块 |
| 06真实切片比较 | PASS（描述性扩展） | 两独立K8模型，共同232genes；配对r .472–.993；不做组间推断 |
| 07空间统计 | PASS（探索性扩展） | 对称6NN、999联合标签置换，Moran/邻域双侧p及明确BH；Moran I .089–.665；不控制层/梯度 |
| 输出回归检查 | PASS | ID匹配、θ/β行和、论文255/K12、UNASSIGNED边界、Moran直接公式、ORA范围、配对唯一性及27三格式完整性 |

## 5. 每幅图的状态与位置

计算成功与原图精确复现分开判定。下表PASS仅为本项目指定计算/图型；原论文图型维持PARTIAL。逐图十项解读在docs/04_图件解读指南.md，图件/source-data映射见figures/figure_manifest.csv。

| 图件stem | 运行/输出 | 原图关系 | PNG路径（另有PDF/SVG） |
|---|---|---|---|
| Figure01_Official_K_Selection | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure01_Official_K_Selection.png` |
| Figure01_Raw_K_Diagnostics | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/Figure01_Raw_K_Diagnostics.png` |
| Figure01_Topic_Proportions | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/Figure01_Topic_Proportions.png` |
| Figure01_Seed_Stability | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/Figure01_Seed_Stability.png` |
| Figure02_All_Topics | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure02_All_Topics.png` |
| Figure02_Topic_Maps | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure02_Topic_Maps.png` |
| Figure03_Topic_Gene_Heatmap | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/Figure03_Topic_Gene_Heatmap.png` |
| Figure03_Marker_Count_Maps | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure03_Marker_Count_Maps.png` |
| Figure04_Layer_Beta_Correlation | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure04_Layer_Beta_Correlation.png` |
| Figure04_Layer_Theta_Correlation | 实际执行/三格式 | PASS指定范围 | `figures/official_reproduction/Figure04_Layer_Theta_Correlation.png` |
| 05A_ORA_dotplot | 实际执行/三格式 | PARTIAL模块 | `figures/extended_analysis/05A_ORA_dotplot.png` |
| 05B_ORA_barplot | 实际执行/三格式 | PARTIAL模块 | `figures/extended_analysis/05B_ORA_barplot.png` |
| 05C_enrichment_map | 实际执行/三格式 | PARTIAL模块 | `figures/extended_analysis/05C_enrichment_map.png` |
| 05D_GSEA_pathway_heatmap | 实际执行/三格式 | PARTIAL模块 | `figures/extended_analysis/05D_GSEA_pathway_heatmap.png` |
| 06A_cross_sample_similarity | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/06A_cross_sample_similarity.png` |
| 06B_descriptive_composition | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/06B_descriptive_composition.png` |
| 06C_matched_spatial_topics | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/06C_matched_spatial_topics.png` |
| 07A_spatial_knn_graph | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/07A_spatial_knn_graph.png` |
| 07B_Moran_I | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/07B_Moran_I.png` |
| 07C_Topic_neighbor_association | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/07C_Topic_neighbor_association.png` |
| 07D_descriptive_distance_decay | 实际执行/三格式 | PASS指定范围 | `figures/extended_analysis/07D_descriptive_distance_decay.png` |
| Paper_mOB_K_Selection | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_mOB_K_Selection.png` |
| Paper_Figure4a_K12_Composition | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_Figure4a_K12_Composition.png` |
| Paper_Figure4b_Current_Topic7 | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_Figure4b_Current_Topic7.png` |
| Paper_mOB_All_Topic_Maps | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_mOB_All_Topic_Maps.png` |
| Paper_Figure4c_Sox11_Nrep_Counts | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_Figure4c_Sox11_Nrep_Counts.png` |
| Paper_mOB_Marker_Beta_Profiles | 实际执行/三格式 | PARTIAL图型 | `figures/official_reproduction/paper_mOB/Paper_mOB_Marker_Beta_Profiles.png` |

共27幅图、81个图件文件。PNG约300dpi，SVG XML可解析，PDF实际最小字级8pt；R后端独立完成所有图件，视觉审查与自动工具局限见environment/figure_QA.md。修复了Unicode设备标签显示问题，仅重新绘图、复用论文模型。

## 6. 耗时与内存

首次01约99.6秒R模块时间、02约63.4秒、03约9.0秒、04约2.3秒；论文01b首次约1176.9秒（19.6分钟）；05约26.4秒、06约52.1秒、07约2.6秒。绘图修复后的缓存重跑另外记录，不替代首次训练时长。

R堆高水位：01约371.6MB、论文分支449.6MB、05约573.9MB。Windows直接x64 R进程采样峰值：论文首次约591.6MB、05约910.7MB、06约723.7MB、07约431.5MB。最早01/02父进程约7MB仅为启动包装器，不能用作模型内存；其真实R堆指标另列。未监控整机和所有子进程总峰值，不宣称这些数字就是总内存上界。process_metrics.csv与run_metrics.csv保留原始指标/范围；默认单核，没有大型并行。

## 7. 错误及修复

environment/errors_and_fixes.md记录数据对齐、ggplot数字Topic名尺度冲突、fgsea默认worker加载冲突、下载403、NA和字体审查问题。04的失败退出记录仍在进程指标中，最终成功结果不抹去失败尝试。没有升级/降级R或替换算法。

## 8. 文档与可复用性

README包含方法背景、原论文、LDA术语与数学/流程图、可复制安装/运行命令、逐步参数合同、实际图件及八项解读、联合应用、真实案例和FAQ。docs/01–07及扩展方法说明齐全；图件手册逐一覆盖27幅并提供10项阅读内容。00–07与01b为可编辑R脚本，按输入/参数/检查/分析/导出组织。新正式数据由研究者在VS Code执行与审查，PI章全部是未来设计，没有编造PI结果。

## 9. GitHub状态

**PASS（发布）**：[公开仓库](https://github.com/stevenhuangfannk-creator/STdeconvolve-Spatial-Topic-Atlas)已创建，main首次推送成功。基准提交`3856a3a`的199个远程blob与本地索引完全一致；README的8张嵌入图均返回HTTP200。Markdown与图件相对路径本地检查通过，全部81个图件文件包含在远程tree中。发布检查记录见environment/remote_validation.csv；最终报告/检查记录随后随文档提交推送，最终HEAD以GitHub页面及git远程为准。原始数据、本机包库、大模型及私人路径没有上传。

## 10. 尚未完成与后续建议

1. Reactome可靠小鼠映射受403限制：待官方小文件恢复或提供带来源/版本的合法文件，再独立重跑05；不以其他物种库代替。
2. GSEA9条NA需先检查基因排序分布与小背景的适用性，不能只增加排列次数或降低阈值来追求阳性。
3. 原论文MERFISH真值、参考方法比较、癌症/其他平台、完整四切片和Allen ISH尚未运行，登记NOT_RUN；随机种子稳定性不能代替它们。
4. 当前Topic7在论文分支同时具有较高Sox11/Nrep β，是线索，但未完成与2022原模型谱的对应审查，不直接借用论文X7身份。
5. 真实PI/Healthy应先确定生物重复、共同特征、平台/深度与批次、Topic可比性及组织结构受限零模型，再做局部可验证分析。机械表达、空间邻近、Topic身份都不等于已证实刚度/通讯/新亚群。

本次没有把只生成脚本、语法检查、下载资料或官方原图当作真实分析PASS。
