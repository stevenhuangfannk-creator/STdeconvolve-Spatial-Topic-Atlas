# 目的：在两个真实公开切片 Rep8（内置 mOB）与 Rep2 间匹配独立拟合的 Topic。
# 输入：01 core.rds；官方 Combining Datasets 提供的 Rep2 count TSV。
# 输出：共同基因、两个独立模型、Pearson 相似性、LSAP 对应、切片平均 theta 与图件。
# 科学边界：这是固定 K/共同基因的小规模扩展，非完整官方四切片实验；两个切片不支持组间推断。
# 运行：核心复现通过后在根目录 source("06_Multi_Sample_Analysis.R")。
# 依赖：00_Start.R 的 STdeconvolve/topicmodels/clue；digest 用于原始文件校验。

# ---- 可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
rep2_url <- "https://www.spatialresearch.org/wp-content/uploads/2016/07/Rep2_MOB_count_matrix-1.tsv"
rep2_file <- "data/raw/Rep2_MOB_count_matrix-1.tsv"
allow_download <- TRUE             # 数据不可得时停止并记录 BLOCKED，不用同一切片拆分替代。
refit_models <- FALSE              # 默认复用本脚本保存且参数相同的模型；改参数须 TRUE。
min_lib_size <- 100                # Rep2 spot 总 counts >100；原始 counts 下的 QC。
min_reads <- 1                     # Rep2 gene 总 counts >1。
min_detected <- 1                  # Rep2 gene 在 >1 个 spot 中检出。
seed_rep8 <- cfg$seed              # 固定独立模型种子。
seed_rep2 <- cfg$seed + 1L

# ---- 输入与来源 / Inputs ----
if (!requireNamespace("digest", quietly = TRUE)) stop("06 缺少 digest；不会自动安装。")
core <- readRDS("data/processed/core.rds")
K <- core$selected_k              # 共用核心选出的 K；不是对 Rep2 重新选择的最优 K。
out <- "results/multi_sample"
fig <- "figures/extended_analysis"
download_problem <- tryCatch({
  if (!file.exists(rep2_file)) {
    if (!allow_download) stop("下载关闭且没有 Rep2 原始文件。")
    utils::download.file(rep2_url, rep2_file, mode = "wb", quiet = TRUE)
  }
  if (file.info(rep2_file)$size < 1000) stop("Rep2 文件过小，不是可信 count TSV。")
  NULL
}, error = function(e) conditionMessage(e))
if (!is.null(download_problem)) {
  utils::write.csv(data.frame(status = "BLOCKED", source = rep2_url, reason = download_problem),
                   file.path(out, "status.csv"), row.names = FALSE)
  stop("官方 Rep2 不可获取：", download_problem, "。未更换数据。")
}
utils::write.csv(data.frame(sample = c("Rep8", "Rep2"),
  source = c("STdeconvolve::mOB; official Combining Datasets identifies it as Rep8", rep2_url),
  retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  sha256 = c(digest::digest(core$counts, algo = "sha256"), digest::digest(rep2_file, file = TRUE, algo = "sha256")),
  hash_scope = c("serialized cleaned Rep8 counts", "original Rep2 TSV")),
  file.path(out, "data_provenance.csv"), row.names = FALSE)
rep2_spot_gene <- as.matrix(utils::read.delim(rep2_file, header = TRUE, row.names = 1,
                                              check.names = FALSE, stringsAsFactors = FALSE))
stopifnot(is.numeric(rep2_spot_gene), all(is.finite(rep2_spot_gene)),
          all(rep2_spot_gene >= 0), all(rep2_spot_gene == floor(rep2_spot_gene)),
          !anyDuplicated(rownames(rep2_spot_gene)), !anyDuplicated(colnames(rep2_spot_gene)))
rep2_counts <- STdeconvolve::cleanCounts(t(rep2_spot_gene), min.lib.size = min_lib_size,
  min.reads = min_reads, min.detected = min_detected, plot = FALSE, verbose = TRUE)
# 仅沿用 Rep8 核心的已选基因，再与 Rep2 实测基因取交集；不做同名补零。
# 这是有 Rep8 特征选择偏倚的示范。正式多样本设计应预先确定共同选基因规则。
shared_genes <- intersect(colnames(core$beta), rownames(rep2_counts))
rep8 <- t(as.matrix(core$counts[shared_genes, rownames(core$pos), drop = FALSE]))
rep2 <- t(as.matrix(rep2_counts[shared_genes, , drop = FALSE]))
shared_genes <- shared_genes[colSums(rep8) > 0 & colSums(rep2) > 0]
rep8 <- rep8[, shared_genes, drop = FALSE]
rep2 <- rep2[, shared_genes, drop = FALSE]
# 共同特征空间内总 counts=0 的 spot 无法拟合；明确保存排除清单。
empty8 <- rownames(rep8)[rowSums(rep8) == 0]
empty2 <- rownames(rep2)[rowSums(rep2) == 0]
excluded <- rbind(data.frame(sample = rep("Rep8", length(empty8)), spot = empty8),
                  data.frame(sample = rep("Rep2", length(empty2)), spot = empty2))
utils::write.csv(excluded, file.path(out, "zero_shared_gene_spots.csv"), row.names = FALSE)
rep8 <- rep8[rowSums(rep8) > 0, , drop = FALSE]
rep2 <- rep2[rowSums(rep2) > 0, , drop = FALSE]
stopifnot(length(shared_genes) > K, nrow(rep8) > K, nrow(rep2) > K,
          identical(colnames(rep8), colnames(rep2)))
utils::write.csv(data.frame(gene = shared_genes), file.path(out, "shared_genes.csv"), row.names = FALSE)
utils::write.csv(data.frame(sample = c("Rep8", "Rep2"), spots = c(nrow(rep8), nrow(rep2)),
  shared_genes = length(shared_genes), K = K, seed = c(seed_rep8, seed_rep2)),
  file.path(out, "input_summary.csv"), row.names = FALSE)
params <- list(K = K, seed = c(Rep8 = seed_rep8, Rep2 = seed_rep2),
  min_lib_size = min_lib_size, min_reads = min_reads, min_detected = min_detected,
  genes = shared_genes, counts_sha256 = c(digest::digest(rep8, algo = "sha256"), digest::digest(rep2, algo = "sha256")))

# ---- 相同特征空间上的独立模型 / Independent fitting ----
models <- results <- list()
for (sample_name in c("Rep8", "Rep2")) {
  x <- if (sample_name == "Rep8") rep8 else rep2
  model_file <- file.path(out, paste0(sample_name, "_independent_model.rds"))
  if (file.exists(model_file) && !refit_models) {
    cached <- readRDS(model_file)
    if (!identical(cached$params, params)) stop("06 缓存参数已改变，请确认后设 refit_models=TRUE。")
    models[[sample_name]] <- cached$model
  } else {
    # 官方 fitLDA 的 Kneed 步骤面向多 K；固定单 K 使用其内部完全相同的 LDA 控制参数。
    models[[sample_name]] <- topicmodels::LDA(slam::as.simple_triplet_matrix(x), k = K,
      control = list(seed = unname(params$seed[sample_name]), verbose = 0, keep = 1, estimate.alpha = TRUE))
    saveRDS(list(model = models[[sample_name]], params = params), model_file)
  }
  results[[sample_name]] <- STdeconvolve::getBetaTheta(models[[sample_name]], perc.filt = 0, betaScale = 1)
}
correlation <- STdeconvolve::getCorrMtx(as.matrix(results$Rep8$beta),
                                       as.matrix(results$Rep2$beta), type = "b", thresh = NULL)
stopifnot(all(is.finite(correlation)))
matched <- STdeconvolve::lsatPairs(correlation)
pairs <- data.frame(Rep8_topic = rownames(correlation)[matched$rowix],
                    Rep2_topic = colnames(correlation)[matched$colsix],
                    pearson_r = correlation[cbind(matched$rowix, matched$colsix)])
utils::write.csv(correlation, file.path(out, "topic_beta_pearson.csv"))
utils::write.csv(pairs, file.path(out, "LSAP_topic_pairs.csv"), row.names = FALSE)
composition <- rbind(data.frame(sample = "Rep8", aligned_topic = pairs$Rep8_topic,
  original_topic = pairs$Rep8_topic, mean_theta = colMeans(results$Rep8$theta)[pairs$Rep8_topic]),
  data.frame(sample = "Rep2", aligned_topic = pairs$Rep8_topic,
  original_topic = pairs$Rep2_topic, mean_theta = colMeans(results$Rep2$theta)[pairs$Rep2_topic]))
utils::write.csv(composition, file.path(out, "descriptive_mean_composition.csv"), row.names = FALSE)
saveRDS(list(results = results, correlation = correlation, pairs = pairs,
             composition = composition, params = params), file.path(out, "multi_sample.rds"))

# ---- 相似性、组成、真实空间位置 / Figures ----
cor_df <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
names(cor_df) <- c("Rep8_topic", "Rep2_topic", "pearson_r")
cor_df$matched <- paste(cor_df$Rep8_topic, cor_df$Rep2_topic) %in% paste(pairs$Rep8_topic, pairs$Rep2_topic)
heat <- ggplot2::ggplot(cor_df, ggplot2::aes(Rep2_topic, Rep8_topic, fill = pearson_r)) +
  ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%s", pearson_r, ifelse(matched, "*", ""))), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B35858", limits = c(-1, 1), name = "Pearson r") +
  ggplot2::coord_equal() + ggplot2::labs(title = "Independent Topic profiles: Rep8 vs Rep2",
    subtitle = "* LSAP one-to-one assignment; assignment does not establish cell identity", x = "Rep2 Topic", y = "Rep8 Topic")
save_figure(heat, file.path(fig, "06A_cross_sample_similarity"), width = 8, height = 7)
comp <- ggplot2::ggplot(composition, ggplot2::aes(aligned_topic, mean_theta, fill = sample)) +
  ggplot2::geom_col(position = "dodge", width = 0.75) +
  ggplot2::scale_fill_manual(values = c(Rep8 = "#547B99", Rep2 = "#B98968")) +
  ggplot2::labs(title = "Mean fitted Topic fraction across two slices", subtitle = "Descriptive comparison: one slice per sample, no group-level inference",
    x = "Aligned Topic (Rep8 label)", y = "Mean raw theta")
save_figure(comp, file.path(fig, "06B_descriptive_composition"), width = 9, height = 5)
rep2_parts <- strsplit(rownames(rep2), "x", fixed = TRUE)
if (all(lengths(rep2_parts) == 2L)) {
  pos2 <- t(vapply(rep2_parts, as.numeric, numeric(2)))
  if (all(is.finite(pos2))) {
    spatial <- list()
    for (j in seq_len(nrow(pairs))) {
      spatial[[paste0("Rep8", j)]] <- data.frame(sample = "Rep8", topic = pairs$Rep8_topic[j],
        x = core$pos[rownames(rep8), "x"], y = core$pos[rownames(rep8), "y"], theta = results$Rep8$theta[, pairs$Rep8_topic[j]])
      spatial[[paste0("Rep2", j)]] <- data.frame(sample = "Rep2", topic = pairs$Rep8_topic[j],
        x = pos2[, 1], y = pos2[, 2], theta = results$Rep2$theta[, pairs$Rep2_topic[j]])
    }
    spatial <- do.call(rbind, spatial)
    sp <- ggplot2::ggplot(spatial, ggplot2::aes(x, y, colour = theta)) + ggplot2::geom_point(size = 1.1) +
      ggplot2::facet_grid(sample ~ topic) + ggplot2::coord_equal() +
      ggplot2::scale_colour_gradient(low = "grey93", high = "#2166AC", limits = c(0, 1)) +
      ggplot2::labs(title = "Matched Topics in each slice's original coordinates", subtitle = "Slices are not spatially registered; coordinate units are array units", x = "Array x", y = "Array y")
    save_figure(sp, file.path(fig, "06C_matched_spatial_topics"), width = max(10, nrow(pairs) * 2), height = 7)
  }
}
utils::write.csv(data.frame(status = "PASS", scope = "Two real public slices; fixed K; shared selected genes; descriptive only"),
                   file.path(out, "status.csv"), row.names = FALSE)
record_step("06_Multi_Sample_Analysis", start_time, "Rep8/Rep2 independent fixed-K fit; shared features; no replicate-level test")
