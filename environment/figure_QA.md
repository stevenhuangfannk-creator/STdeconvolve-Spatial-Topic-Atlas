# 图件 QA / Figure quality record

图件合同：每图的结论、证据、变量、重复单位与边界见 docs/04_图件解读指南.md。所有图为真实定量图（quantitative grid），没有AI生成数据或替代量化面板。核心为一个切片，seed为计算重复；两真实切片只作描述。

## 实际检查

- 27个PNG已由R逐个读取并用于视觉审查；PNG签名、三格式齐全检查通过，文件分辨率约300dpi。
- 27个SVG均通过XML解析；27个PDF可打开并提取可编辑文字。实际变换后的最小字级8pt，无已提取文字越过页面边界。
- 全部面板已检查标签、legend、坐标、背景与数值表对应；最终修复了ASCII以外乘号/连接符在本机设备上的异常。01b只复用原模型重绘，首次19.6分钟拟合日志另存。
- 图为教学尺寸，非某期刊最终单/双栏排版；缩小时应再次检查。官方组成图保留原教程的定性rainbow颜色，细组织层边较浅，此项属于官方风格保留。连续量采用顺序/发散色标，扩展矩阵不使用连续rainbow。

## 自动工具的人工复核

source_preflight.json保留了nature-figure静态预检原输出，不将其中FAIL隐藏或改写成自动PASS。FONT-FAMILY与导出项在单脚本出现FAIL是因为统一Arial字体和ragg/cairo_pdf/svglite导出由00_Start.R的cfg与save_figure提供，静态工具没有跨文件解析；已逐文件核对真实三格式、实际字体与字级。COLOR-MAP提示的rainbow限于用户要求接近官方的定性Topic/组织层配色，不作为连续色图使用。LOG-GUARD等提示由输入非负计数、log1p或显式伪计数和实际输出检查闭合；没有修改数据以消除警告。

技能中的audit_pdf_text.py对Cairo仅读取Tf=1，没有应用Tm文本变换，会错误报告1pt。实际PDF使用文本矩阵缩放，已用PyMuPDF读取最终变换后span大小复核，见pdf_text_audit.csv；最小8pt。该复核只读取文档对象，没有用Python重新渲染或修改任何R图。

## 每幅图记录

除种子箱图明确描述Topic间分布外，所有图不绘生物学重复误差条：它们展示单模型、单切片或两个未分组切片；不虚构置信区间。05GSEA与07统计在图注/手册明确校正范围和零假设。种子箱体不是供体不确定性，06匹配星号不是显著性。

| 图件 | 实际最小字级pt | 可提取文字span数 | 文字越界 | 最终检查 |
|---|---:|---:|---|---|
| figures/extended_analysis/05A_ORA_dotplot.pdf | 8.8 | 40 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/05B_ORA_barplot.pdf | 8.8 | 24 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/05C_enrichment_map.pdf | 8.54 | 22 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/05D_GSEA_pathway_heatmap.pdf | 8.8 | 38 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/06A_cross_sample_similarity.pdf | 8.54 | 90 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/06B_descriptive_composition.pdf | 8.8 | 19 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/06C_matched_spatial_topics.pdf | 8.8 | 70 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/07A_spatial_knn_graph.pdf | 8.8 | 18 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/07B_Moran_I.pdf | 8.8 | 18 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/07C_Topic_neighbor_association.pdf | 8.8 | 82 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/07D_descriptive_distance_decay.pdf | 8.8 | 40 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/Figure01_Raw_K_Diagnostics.pdf | 8.8 | 32 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/Figure01_Seed_Stability.pdf | 8.8 | 15 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/Figure01_Topic_Proportions.pdf | 8.8 | 17 | 无检测到 | 已核对面板/图例与边界 |
| figures/extended_analysis/Figure03_Topic_Gene_Heatmap.pdf | 8.0 | 50 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure01_Official_K_Selection.pdf | 8.8 | 17 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure02_All_Topics.pdf | 8.8 | 17 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure02_Topic_Maps.pdf | 8.0 | 82 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure03_Marker_Count_Maps.pdf | 8.0 | 63 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure04_Layer_Beta_Correlation.pdf | 11.0 | 60 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/Figure04_Layer_Theta_Correlation.pdf | 11.0 | 60 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_Figure4a_K12_Composition.pdf | 8.8 | 21 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_Figure4b_Current_Topic7.pdf | 8.8 | 21 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_Figure4c_Sox11_Nrep_Counts.pdf | 8.8 | 19 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_mOB_All_Topic_Maps.pdf | 8.0 | 121 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_mOB_K_Selection.pdf | 8.8 | 53 | 无检测到 | 已核对面板/图例与边界 |
| figures/official_reproduction/paper_mOB/Paper_mOB_Marker_Beta_Profiles.pdf | 8.8 | 38 | 无检测到 | 已核对面板/图例与边界 |
