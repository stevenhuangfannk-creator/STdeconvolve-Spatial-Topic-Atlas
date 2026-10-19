# 目的：使用论文 mOB 预处理与固定 K=12，部分复现论文 Figure 4a-c 的图型。
# 输入：官方 STdeconvolve::mOB；不读取或覆盖快速教程 core.rds。
# 输出：data/processed/paper_mOB.rds；results/models/paper_mOB；figures/official_reproduction/paper_mOB。
# 科学边界：PARTIAL 图型复现；当前 Topic 编号未对齐论文，不能宣称逐像素复现或独立层验证。
# Figure 4d 是 Allen Brain Atlas ISH 原图，本脚本不生成替代图，记录 NOT_RUN。
# 来源：论文 Supplementary Methods（mOB Rep8）、官方 combining_datasets.html 的 mOB 步骤。
# 运行：仓库根目录 source("01b_Paper_mOB.R")；可独立于 01-07 执行。

# ---- 集中可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
paper_ks <- 2:20                 # 论文搜索范围；Combining 网页示例为 2:15，此处依论文 2:20。
paper_k <- 12L                   # 论文固定选择：低 perplexity 范围且 mean<5% 的 rare Topic 数为 0。
paper_seed <- cfg$seed            # 当前复现固定 seed 0，不保证与 2022 年原始 Topic 编号相同。
marker_genes <- c("Sox11", "Nrep") # 论文 Figure 4c 的真实 marker，不用其他基因代替。
pie_radius <- 0.4
map_columns <- 3
highlight_topic <- "7"           # 仅当前模型 Topic 7；未与论文 X7 进行跨模型匹配。
paper_settings <- list(Ks = paper_ks, K = paper_k, seed = paper_seed,
                        min.reads = 100, min.lib.size = 100, min.detected = 1,
                        removeAbove = NA, removeBelow = NA, nTopOD = NA,
                        od.genes.alpha = 0.05, gam.k = 5, ncores = 1L)
model_dir <- "results/models/paper_mOB"
table_dir <- "results/theta_beta/paper_mOB"
figure_dir <- "figures/official_reproduction/paper_mOB"
for (d in c(model_dir, table_dir, figure_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
saveRDS(paper_settings, file.path(model_dir, "parameters.rds"))
capture.output(dput(paper_settings), file = file.path(model_dir, "parameters.txt"))

# ---- 官方 preprocess：参数名已与函数定义核对 / Official preprocessing ----
data("mOB", package = "STdeconvolve")
stopifnot(!is.null(names(mOB$annot)), all(colnames(mOB$counts) %in% rownames(mOB$pos)))
# 官方 preprocess 即使 plot=FALSE 也会设置 par；用临时空设备避免产生空白 Rplots.pdf。
grDevices::pdf(file = NULL)
paper_input <- STdeconvolve::preprocess(
  dat = t(as.matrix(mOB$counts)), extractPos = FALSE,
  selected.genes = NA, nTopGenes = NA, genes.to.remove = NA,
  removeAbove = NA, removeBelow = NA, min.reads = 100, min.lib.size = 100,
  min.detected = 1, ODgenes = TRUE, nTopOD = NA, od.genes.alpha = 0.05,
  gam.k = 5, verbose = TRUE, plot = FALSE)
grDevices::dev.off()
# preprocess$corpus 为 spot * gene，与快速教程保存的 gene * spot corpus 方向不同。
paper_corpus <- paper_input$corpus
stopifnot(all(rowSums(paper_corpus) > 0), all(rownames(paper_corpus) %in% names(mOB$annot)))
paper_pos <- as.data.frame(mOB$pos[rownames(paper_corpus), c("x", "y"), drop = FALSE])
paper_annot <- mOB$annot[rownames(paper_corpus)]
dimension_check <- data.frame(actual_spots = nrow(paper_corpus), actual_genes = ncol(paper_corpus),
                               published_spots = 260, published_genes = 255,
                               matches_published = nrow(paper_corpus) == 260 && ncol(paper_corpus) == 255)
utils::write.csv(dimension_check, file.path(model_dir, "preprocessing_dimensions.csv"), row.names = FALSE)
if (!dimension_check$matches_published) warning("当前预处理维度与论文 260 * 255 不同；保留实际输出，不为匹配图件修改数据。")

# ---- 论文 K 网格；保留已完成的同参数拟合 / Model fitting ----
model_path <- file.path(model_dir, "ldas.rds")
if (file.exists(model_path)) {
  paper_ldas <- readRDS(model_path)
  stopifnot(identical(attr(paper_ldas, "paper_settings"), paper_settings))
  message("复用已完成、参数一致的论文 mOB K 网格模型。")
} else {
  paper_ldas <- STdeconvolve::fitLDA(paper_corpus, Ks = paper_ks, seed = paper_seed,
                                    ncores = 1L, perc.rare.thresh = cfg$rare,
                                    plot = FALSE, verbose = TRUE)
  attr(paper_ldas, "paper_settings") <- paper_settings
  saveRDS(paper_ldas, model_path)
}
raw_by_k <- lapply(paper_ldas$models, function(m) {
  STdeconvolve::getBetaTheta(m, perc.filt = 0, verbose = FALSE)$theta
})
paper_diagnostics <- data.frame(K = paper_ks,
                                 in_sample_perplexity = as.numeric(paper_ldas$perplexities),
                                 alpha = vapply(paper_ldas$models, function(m) methods::slot(m, "alpha"), numeric(1)),
                                 rare_topics_raw = vapply(raw_by_k, function(t) sum(colMeans(t) < cfg$rare), integer(1)),
                                 fixed_paper_K12 = paper_ks == paper_k,
                                 current_minimum_perplexity = paper_ks %in% paper_ldas$minOptK)
utils::write.csv(paper_diagnostics, file.path(model_dir, "K_diagnostics.csv"), row.names = FALSE)
utils::write.csv(data.frame(paper_fixed_K = paper_k,
                            current_min_perplexity_K = paste(paper_ldas$minOptK, collapse = ";"),
                            current_rare_topics_at_K12 = paper_diagnostics$rare_topics_raw[paper_diagnostics$K == paper_k],
                            decision = "K12固定对应论文选择；不替换为本次最小perplexity。低perplexity范围无官方数值阈值，勿伪造自动选择规则。"),
                 file.path(model_dir, "K_selection_record.csv"), row.names = FALSE)
p_k <- STdeconvolve::perplexityPlot(paper_ldas, perc.rare.thresh = cfg$rare) +
  ggplot2::labs(caption = "Paper-style preprocessing; K=2-20; fixed paper K=12. Perplexity is in-sample; official rare counts use display filtering.")
save_figure(p_k, file.path(figure_dir, "Paper_mOB_K_Selection"), 10, 5)

# ---- 固定 K=12 的结果，不改快速教程对象 / Fixed paper K ----
paper_model <- STdeconvolve::optimalModel(paper_ldas, opt = paper_k)
paper_raw <- STdeconvolve::getBetaTheta(paper_model, perc.filt = 0, betaScale = cfg$beta_scale)
paper_display <- STdeconvolve::getBetaTheta(paper_model, perc.filt = cfg$perc_filt, betaScale = cfg$beta_scale)
stopifnot(identical(rownames(paper_raw$theta), rownames(paper_pos)))
saveRDS(list(corpus = paper_corpus, pos = paper_pos, annot = paper_annot,
              theta_raw = paper_raw$theta, theta = paper_display$theta, beta = paper_raw$beta,
              model = paper_model, settings = paper_settings, status = "PARTIAL_FIGURE_REPRODUCTION"),
        "data/processed/paper_mOB.rds")
utils::write.csv(data.frame(spot_id = rownames(paper_raw$theta), paper_raw$theta, check.names = FALSE),
                 file.path(table_dir, "theta_raw_K12.csv"), row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(paper_display$theta), paper_display$theta, check.names = FALSE),
                 file.path(table_dir, "theta_display_K12.csv"), row.names = FALSE)
utils::write.csv(data.frame(topic = rownames(paper_raw$beta), paper_raw$beta, check.names = FALSE),
                 file.path(table_dir, "beta_scaled_K12.csv"), row.names = FALSE)
paper_cols <- stats::setNames(grDevices::rainbow(paper_k), colnames(paper_raw$theta))
paper_layer_cols <- stats::setNames(grDevices::rainbow(nlevels(paper_annot)), levels(paper_annot))
p_pie <- STdeconvolve::vizAllTopics(paper_display$theta, paper_pos,
                                    topicCols = unname(paper_cols[colnames(paper_display$theta)]),
                                    groups = paper_annot, group_cols = paper_layer_cols, r = pie_radius,
                                    plotTitle = "Paper mOB: current K=12 composition") +
  ggplot2::labs(caption = "Partial reproduction of Fig. 4a graph type; Topic labels and colours are not aligned to the 2022 fit.")
save_figure(p_pie, file.path(figure_dir, "Paper_Figure4a_K12_Composition"), 10, 5.5)
if (highlight_topic %in% colnames(paper_display$theta)) {
  highlight_cols <- rep("grey90", ncol(paper_display$theta))
  highlight_cols[colnames(paper_display$theta) == highlight_topic] <- "black"
  p_highlight <- STdeconvolve::vizAllTopics(paper_display$theta, paper_pos,
                                            topicCols = highlight_cols, groups = paper_annot,
                                            group_cols = paper_layer_cols, r = pie_radius,
                                            plotTitle = paste("Current Topic", highlight_topic, "(unmatched to published X7)")) +
    ggplot2::labs(caption = "Fig. 4b-style display only: current Topic 7 is not established as the paper's X7.")
  save_figure(p_highlight, file.path(figure_dir, "Paper_Figure4b_Current_Topic7"), 10, 5.5)
}
maps <- lapply(colnames(paper_raw$theta), function(id) {
  STdeconvolve::vizTopic(paper_raw$theta, paper_pos, topic = id, size = 3,
                        stroke = 0.15, high = paper_cols[id], plotTitle = paste("Topic", id)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
})
save_figure(patchwork::wrap_plots(maps, ncol = map_columns) +
              patchwork::plot_annotation(title = "All 12 current Topics: raw theta"),
            file.path(figure_dir, "Paper_mOB_All_Topic_Maps"), 12, 13)

# ---- Figure 4c 的真实基因；缺失时记录而不替换 / Measured marker counts ----
marker_status <- data.frame(gene = marker_genes, present_in_counts = marker_genes %in% rownames(mOB$counts),
                             present_in_beta = marker_genes %in% colnames(paper_raw$beta))
utils::write.csv(marker_status, file.path(table_dir, "marker_availability.csv"), row.names = FALSE)
available_genes <- marker_genes[marker_status$present_in_counts]
if (length(available_genes)) {
  gene_df <- data.frame(paper_pos, t(as.matrix(mOB$counts[available_genes, rownames(paper_pos), drop = FALSE])), check.names = FALSE)
  gene_maps <- lapply(available_genes, function(g) {
    STdeconvolve::vizGeneCounts(gene_df, gene = g, size = 3, stroke = 0.1,
                               winsorize = 0, plotTitle = g)
  })
  save_figure(patchwork::wrap_plots(gene_maps, ncol = length(available_genes)) +
                patchwork::plot_annotation(title = "Published mOB markers: measured raw counts",
                                           caption = "Fig. 4c graph types; each gene has its own count scale. Fig. 4d ISH images are not reproduced."),
              file.path(figure_dir, "Paper_Figure4c_Sox11_Nrep_Counts"), 10, 5)
}
profile_genes <- marker_genes[marker_status$present_in_beta]
if (length(profile_genes)) {
  profile <- as.data.frame(as.table(paper_raw$beta[, profile_genes, drop = FALSE]), stringsAsFactors = FALSE)
  names(profile) <- c("topic", "gene", "beta_scaled")
  utils::write.csv(profile, file.path(table_dir, "published_marker_beta_profiles.csv"), row.names = FALSE)
  p_profile <- ggplot2::ggplot(profile, ggplot2::aes(topic, beta_scaled, fill = topic)) +
    ggplot2::geom_col() + ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::scale_fill_manual(values = paper_cols, guide = "none") +
    ggplot2::labs(x = "Current Topic (unmatched to paper)", y = "Beta probability * 1000",
                  title = "Sox11 and Nrep across current K=12 Topic profiles")
  save_figure(p_profile, file.path(figure_dir, "Paper_mOB_Marker_Beta_Profiles"), 10, 4.5)
}
utils::write.csv(data.frame(panel = c("Fig4a", "Fig4b", "Fig4c", "Fig4d"),
                            status = c("PARTIAL", "PARTIAL", "PARTIAL", "NOT_RUN"),
                            boundary = c("相同数据和K12的图型复现；颜色与Topic顺序未对齐", "当前Topic7未证明等同于论文X7", "真实Sox11/Nrep counts；图型复现，非逐像素校准", "Allen Brain Atlas ISH外部图像未复现；不以表达图替代")),
                 file.path(model_dir, "figure_reproduction_status.csv"), row.names = FALSE)
record_step("01b_Paper_mOB", start_time,
            paste("PARTIAL Fig4a-c graph types; actual", nrow(paper_corpus), "spots", ncol(paper_corpus), "genes; fixed K12; Fig4d NOT_RUN"))
