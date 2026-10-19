# 目的：保存官方 K 诊断图，并扩展未阈值化 Topic 占比和重复随机种子稳定性。
# 输入：01 生成的 core.rds、ldas.rds；不重新执行预处理或 K 网格搜索。
# 输出：results/models/K_*.csv、seed_stability_*.csv/RDS、Figure 01 系列。
# 科学边界：perplexity 在拟合语料上计算，不是 held-out 泛化性能；稳定性不证明生物真实性。
# 运行：仓库根目录 source("02_Topic_Number_Selection.R")。

# ---- 集中可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
repeat_seeds <- cfg$repeats       # 固定选择出的 K，仅重复种子 1、2、3；基准为 seed 0。
rare_fraction <- cfg$rare        # 原始 theta 中跨 spot 平均比例 < 0.05 的 Topic 定义为 rare。

# ---- 固定步骤：输入与官方诊断 / Official diagnostics ----
stopifnot(file.exists("data/processed/core.rds"), file.exists("results/models/ldas.rds"))
core <- readRDS("data/processed/core.rds")
ldas <- readRDS("results/models/ldas.rds")
ks <- vapply(ldas$models, function(m) methods::slot(m, "k"), numeric(1))
raw_means <- lapply(ldas$models, function(m) {
  colMeans(STdeconvolve::getBetaTheta(m, perc.filt = 0, verbose = FALSE)$theta)
})
k_table <- data.frame(K = ks, in_sample_perplexity = as.numeric(ldas$perplexities),
                      alpha = vapply(ldas$models, function(m) methods::slot(m, "alpha"), numeric(1)),
                      rare_topics_raw_theta = vapply(raw_means, function(x) sum(x < rare_fraction), integer(1)),
                      selected_by_minimum = ks == core$selected_k)
utils::write.csv(k_table, "results/models/K_diagnostics.csv", row.names = FALSE)
means_table <- do.call(rbind, lapply(seq_along(ks), function(i) {
  data.frame(K = ks[i], topic = names(raw_means[[i]]), mean_fraction = unname(raw_means[[i]]))
}))
utils::write.csv(means_table, "results/models/K_topic_mean_fractions_raw.csv", row.names = FALSE)
# 官方 perplexityPlot 内部使用默认阈值 theta；保留原函数图，并另给原始 theta 的诊断。
p_official <- STdeconvolve::perplexityPlot(ldas, perc.rare.thresh = rare_fraction)
p_official <- p_official + ggplot2::labs(caption = "Official function; perplexity is in-sample. Rare-topic counts use its default filtered theta.")
save_figure(p_official, "figures/official_reproduction/Figure01_Official_K_Selection", 9, 5)
p_perplexity <- ggplot2::ggplot(k_table, ggplot2::aes(K, in_sample_perplexity)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::geom_vline(xintercept = core$selected_k, linetype = "dashed", colour = "#B2182B") +
  ggplot2::scale_x_continuous(breaks = ks) +
  ggplot2::labs(x = "Number of topics (K)", y = "In-sample perplexity",
                title = "Model fit across K", subtitle = "Dashed line: minimum perplexity in this tested grid")
p_rare <- ggplot2::ggplot(k_table, ggplot2::aes(K, rare_topics_raw_theta)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = ks) +
  ggplot2::labs(x = "Number of topics (K)", y = "Topics with mean raw fraction < 0.05", title = "Rare topics before display filtering")
save_figure(p_perplexity + p_rare + patchwork::plot_layout(ncol = 2),
            "figures/extended_analysis/Figure01_Raw_K_Diagnostics", 11, 4.5)
p_means <- ggplot2::ggplot(means_table, ggplot2::aes(factor(K), mean_fraction, group = interaction(K, topic))) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.10, height = 0, seed = cfg$seed), size = 2) +
  ggplot2::geom_hline(yintercept = rare_fraction, linetype = "dashed", colour = "#B2182B") +
  ggplot2::labs(x = "K", y = "Mean raw theta across spots", title = "All fitted Topic proportions", caption = "Each point is a Topic; raw theta has no 5% display threshold.")
save_figure(p_means, "figures/extended_analysis/Figure01_Topic_Proportions", 7, 4.5)

# ---- 固定 K 的种子稳定性 / Project extension ----
# fitLDA 的 Kneed 诊断面向多 K；单 K 重复用其内部相同 topicmodels::LDA 调用与控制参数。
stopifnot(!anyDuplicated(c(core$seed, repeat_seeds)))
seed_models <- list()
seed_models[[as.character(core$seed)]] <- core$optLDA
for (s in repeat_seeds) {
  message("固定 K=", core$selected_k, "，重复 seed=", s)
  seed_models[[as.character(s)]] <- topicmodels::LDA(
    slam::as.simple_triplet_matrix(t(as.matrix(core$corpus))), k = core$selected_k,
    control = list(seed = s, verbose = 0, keep = 1, estimate.alpha = TRUE))
}
seed_profiles <- lapply(seed_models, function(m) {
  STdeconvolve::getBetaTheta(m, perc.filt = 0, betaScale = 1, verbose = FALSE)
})
seed_pairs <- utils::combn(names(seed_profiles), 2, simplify = FALSE)
matched_rows <- lapply(seed_pairs, function(s) {
  a <- seed_profiles[[s[1]]]
  b <- seed_profiles[[s[2]]]
  corr <- STdeconvolve::getCorrMtx(a$beta, b$beta, type = "b", verbose = FALSE)
  stopifnot(all(is.finite(corr)))
  pairs <- STdeconvolve::lsatPairs(corr)  # 官方 Hungarian 一对一匹配，避免直接按 Topic 编号比较。
  ix <- pairs$rowix
  jx <- pairs$colsix
  data.frame(seed_a = as.integer(s[1]), seed_b = as.integer(s[2]),
             topic_a = rownames(corr)[ix], topic_b = colnames(corr)[jx],
             beta_pearson = corr[cbind(ix, jx)],
             theta_pearson = vapply(seq_along(ix), function(i) {
               stats::cor(a$theta[, ix[i]], b$theta[, jx[i]])
             }, numeric(1)))
})
matches <- do.call(rbind, matched_rows)
utils::write.csv(matches, "results/models/seed_stability_topic_matches.csv", row.names = FALSE)
seed_summary <- aggregate(cbind(beta_pearson, theta_pearson) ~ seed_a + seed_b, matches, mean)
utils::write.csv(seed_summary, "results/models/seed_stability_pair_means.csv", row.names = FALSE)
saveRDS(list(models = seed_models, profiles = seed_profiles, matches = matches,
              selected_k = core$selected_k, matching = "Pearson beta + Hungarian one-to-one"),
        "results/models/seed_stability.rds")
matches$seed_pair <- paste(matches$seed_a, matches$seed_b, sep = " vs ")
p_stability <- ggplot2::ggplot(matches, ggplot2::aes(seed_pair, beta_pearson)) +
  ggplot2::geom_boxplot(outlier.shape = NA, fill = "#D9EAF3", width = 0.55) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.10, height = 0, seed = cfg$seed), size = 1.7) +
  ggplot2::coord_cartesian(ylim = c(-1, 1)) +
  ggplot2::labs(x = "Random-seed pair", y = "Matched beta Pearson correlation",
                title = paste0("Seed stability at fixed K = ", core$selected_k),
                caption = "Project extension: every point is one Hungarian-matched Topic; seed stability is not external validation.")
save_figure(p_stability, "figures/extended_analysis/Figure01_Seed_Stability", 9, 5)
record_step("02_Topic_Number_Selection", start_time,
            paste("fixed K", core$selected_k, "seed repeats", paste(repeat_seeds, collapse = ",")))
