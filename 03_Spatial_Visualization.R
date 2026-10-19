# 目的：使用官方空间绘图函数展示 Topic 组成、单 Topic 分布及代表性基因。
# 输入：01 的 core.rds；beta 是潜在主题的 gene 概率 * 1000，不是测量表达量。
# 输出：Figure 02/03 系列 PNG（300 dpi）/PDF/SVG、Topic 色表与 top markers。
# 科学边界：Topic 是潜在表达成分；高 beta gene 不是经差异检验确认的细胞类型 marker。
# 运行：仓库根目录 source("03_Spatial_Visualization.R")。

# ---- 集中可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
pie_radius <- 0.4                # 官方 mOB 教程单位：原始空间坐标。
spot_size <- 3                  # ggplot2 点大小。
map_columns <- 3                # 多面板每行图数。
n_top_markers <- 10              # 输出每个 Topic 的 top beta genes 数量。
n_heatmap_markers <- 5           # 热图用每个 Topic top 5 的并集。
gene_winsorize <- 0              # 官方 gene count 图不截尾；0 表示保留实际 counts。

# ---- 固定步骤：读取与统一颜色 / Input and colours ----
stopifnot(file.exists("data/processed/core.rds"))
core <- readRDS("data/processed/core.rds")
stopifnot(identical(rownames(core$theta_raw), rownames(core$pos)),
          identical(rownames(core$theta_raw), names(core$annot)))
topic_ids <- colnames(core$theta_raw)
topic_cols <- stats::setNames(grDevices::rainbow(length(topic_ids)), topic_ids)
layer_cols <- stats::setNames(grDevices::rainbow(nlevels(core$annot)), levels(core$annot))
utils::write.csv(data.frame(topic = names(topic_cols), colour = unname(topic_cols)),
                 "results/theta_beta/topic_palette.csv", row.names = FALSE)

# ---- 官方空间图 / Official spatial visualisation ----
display_ids <- colnames(core$theta)
p_pie <- STdeconvolve::vizAllTopics(
  core$theta, core$pos, topicCols = unname(topic_cols[display_ids]),
  groups = core$annot, group_cols = layer_cols, r = pie_radius,
  plotTitle = "mOB: Topic composition (display threshold 5%)") +
  ggplot2::labs(caption = "Official tutorial: low fractions are set to zero and each spot is renormalised. Borders indicate annotated tissue layers.")
save_figure(p_pie, "figures/official_reproduction/Figure02_All_Topics", 10, 5.5)
# 单 Topic 图使用 perc.filt=0 的原始 theta，确保展示所有 K 个 Topic。
topic_plots <- lapply(topic_ids, function(id) {
  STdeconvolve::vizTopic(core$theta_raw, core$pos, topic = id,
                        size = spot_size, stroke = 0.15, low = "white",
                        high = topic_cols[id], plotTitle = paste("Topic", id)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 9),
                    legend.text = ggplot2::element_text(size = 8),
                    plot.title = ggplot2::element_text(size = 11))
})
p_topics <- patchwork::wrap_plots(topic_plots, ncol = map_columns) +
  patchwork::plot_annotation(title = "mOB spatial Topic fractions",
                             caption = "Official vizTopic function; unfiltered theta on the 0-1 proportion scale.")
save_figure(p_topics, "figures/official_reproduction/Figure02_Topic_Maps", 12, 3.4 * ceiling(length(topic_ids) / map_columns))

# ---- 代表性基因 / Top beta genes (project extension) ----
top_genes <- STdeconvolve::topGenes(core$beta, n = min(n_top_markers, ncol(core$beta)))
marker_table <- do.call(rbind, lapply(names(top_genes), function(id) {
  data.frame(topic = id, rank = seq_along(top_genes[[id]]),
             gene = names(top_genes[[id]]), beta_scaled = as.numeric(top_genes[[id]]),
             beta_probability = as.numeric(top_genes[[id]]) / core$beta_scale)
}))
utils::write.csv(marker_table, "results/annotation/top_beta_genes.csv", row.names = FALSE)
heat_genes <- unique(unlist(lapply(top_genes, function(x) head(names(x), n_heatmap_markers))))
heat <- as.data.frame(as.table(core$beta[, heat_genes, drop = FALSE]), stringsAsFactors = FALSE)
names(heat) <- c("topic", "gene", "beta_scaled")
heat$topic <- factor(heat$topic, levels = topic_ids)
heat$gene <- factor(heat$gene, levels = rev(heat_genes))
p_heat <- ggplot2::ggplot(heat, ggplot2::aes(topic, gene, fill = log1p(beta_scaled))) +
  ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(name = "log(1 + beta * 1000)") +
  ggplot2::labs(x = "Topic", y = NULL, title = "Topic gene profiles: top beta genes",
                caption = "Project extension: union of top 5 genes per Topic; no differential-expression test or gene-wise rescaling.") +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), panel.grid = ggplot2::element_blank())
save_figure(p_heat, "figures/extended_analysis/Figure03_Topic_Gene_Heatmap", 8, max(5, length(heat_genes) * 0.18))

# 每 Topic top 1 的基因，用官方 vizGeneCounts 展示真实测得的原始 count 分布。
# 它与 beta 概率是不同对象；测序深度可能影响 count 图，不能解释为 Topic 特异性检验。
map_genes <- unique(vapply(top_genes, function(x) names(x)[1], character(1)))
gene_df <- data.frame(core$pos, t(as.matrix(core$counts[map_genes, rownames(core$pos), drop = FALSE])),
                      check.names = FALSE)
gene_plots <- lapply(map_genes, function(gene) {
  STdeconvolve::vizGeneCounts(gene_df, gene = gene, size = spot_size,
                             stroke = 0.1, winsorize = gene_winsorize,
                             plotTitle = gene) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8),
                    plot.title = ggplot2::element_text(size = 11))
})
p_genes <- patchwork::wrap_plots(gene_plots, ncol = map_columns) +
  patchwork::plot_annotation(title = "Measured counts of top beta genes",
                             caption = "Official vizGeneCounts function; each panel has its own raw-count scale. Library depth is not removed.")
save_figure(p_genes, "figures/official_reproduction/Figure03_Marker_Count_Maps", 12, 3.4 * ceiling(length(map_genes) / map_columns))
record_step("03_Spatial_Visualization", start_time,
            paste(length(topic_ids), "Topic maps;", length(map_genes), "unique top-gene maps"))
