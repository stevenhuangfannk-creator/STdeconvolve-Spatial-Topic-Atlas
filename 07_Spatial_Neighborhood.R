# 目的：以原始 theta 研究 Topic 空间自相关、邻域关联及距离衰减；本模块为项目扩展。
# 输入：01 core.rds 的 theta_raw 与 pos；不使用显示阈值截断后的 theta。
# 输出：对称 kNN 图、Moran's I、跨 Topic 邻域统计、999 次随机标签置换及诊断图。
# 科学边界：探索性无条件 random-label 零假设；不控制组织层/空间梯度，不能解释为细胞通讯。
# 运行：在仓库根目录 source("07_Spatial_Neighborhood.R")；只处理小规模官方示例。
# 依赖：00_Start.R 的 ggplot2；矩阵和置换计算使用 base R，串行运行。

# ---- 可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
k_neighbors <- 6L                 # 每 spot 最近邻数；对称化取任一方向存在的边。
n_permutations <- 999L            # 联合打乱全部 Topic 组成向量的 spot 标签；最小 p=0.001。
q_cutoff <- 0.05                  # Moran 跨所有 Topic BH；邻域跨所有唯一 Topic 对（含对角）BH。
distance_bins <- 8L               # 距离衰减图分箱数；等宽、仅作描述，无显著性声明。
max_example_spots <- 5000L        # 本示例使用 dense 距离矩阵；更大数据需另建稀疏图并评估内存。

# ---- 输入与空间图 / Inputs and graph ----
core <- readRDS("data/processed/core.rds")
theta <- as.matrix(core$theta_raw)
pos <- as.matrix(core$pos[rownames(theta), c("x", "y"), drop = FALSE])
stopifnot(identical(rownames(theta), rownames(pos)), all(is.finite(theta)),
          all(is.finite(pos)), max(abs(rowSums(theta) - 1)) < 1e-8,
          nrow(theta) > k_neighbors, nrow(theta) <= max_example_spots,
          k_neighbors >= 1, n_permutations >= 99, !anyDuplicated(as.data.frame(pos)))
n <- nrow(theta)
out <- "results/neighborhood"
fig <- "figures/extended_analysis"
distance <- as.matrix(stats::dist(pos))
diag(distance) <- Inf
adjacency <- matrix(0, n, n, dimnames = list(rownames(theta), rownames(theta)))
for (i in seq_len(n)) {
  # 等距离按原始 spot 顺序打破并列；保留固定 k，规则明确且可复现。
  neighbors <- order(distance[i, ], seq_len(n))[seq_len(k_neighbors)]
  adjacency[i, neighbors] <- 1
}
adjacency <- 1 * ((adjacency + t(adjacency)) > 0)
diag(adjacency) <- 0
S0 <- sum(adjacency)                # 对称二元权重总和；每条无向边计算两个方向。
stopifnot(S0 > 0, all(rowSums(adjacency) > 0))
edge_index <- which(upper.tri(adjacency) & adjacency == 1, arr.ind = TRUE)
edges <- data.frame(spot_from = rownames(theta)[edge_index[, 1]],
                    spot_to = rownames(theta)[edge_index[, 2]],
                    array_distance = distance[edge_index])
utils::write.csv(edges, file.path(out, "symmetric_knn_edges.csv"), row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(theta), x = pos[, 1], y = pos[, 2],
                            degree = rowSums(adjacency)), file.path(out, "graph_nodes.csv"), row.names = FALSE)
params <- list(k = k_neighbors, graph = "symmetric union kNN; binary weights; zero diagonal",
  units = "original array coordinates; not micrometers", permutations = n_permutations, seed = cfg$seed,
  null = "joint whole-composition random labeling across all observed spots; no spatial/stratum restriction",
  moran_BH_family = "all Topics", pair_BH_family = "all unique Topic pairs including same-Topic pairs")
saveRDS(params, file.path(out, "parameters.rds"))

# ---- 观测统计与联合随机标签零分布 / Statistics and null ----
topic_sd <- apply(theta, 2, stats::sd)
if (any(topic_sd == 0)) stop("存在常量 Topic，Moran/标准化邻域统计不可定义；不静默剔除。")
z <- sweep(sweep(theta, 2, colMeans(theta), "-"), 2, topic_sd, "/")
denominator <- colSums(z^2)
moran <- n / S0 * colSums(z * (adjacency %*% z)) / denominator
expected_moran <- -1 / (n - 1)
# H_ab = 1/S0 * sum_ij w_ij z_ia z_jb；正值表示相邻位置相对均值共同升高。
# 对称二元图使 H_ab=H_ba；不是直接的细胞接触频率或配体受体分数。
neighbor_association <- crossprod(z, adjacency %*% z) / S0
expected_pair <- -crossprod(z) / (n * (n - 1))
pair_index <- which(upper.tri(neighbor_association, diag = TRUE), arr.ind = TRUE)
observed_pair <- neighbor_association[pair_index]
null_pair_mean <- expected_pair[pair_index]
perm_moran <- matrix(NA_real_, n_permutations, ncol(theta), dimnames = list(NULL, colnames(theta)))
perm_pair <- matrix(NA_real_, n_permutations, nrow(pair_index))
set.seed(cfg$seed)
for (r in seq_len(n_permutations)) {
  # 同一个 permutation 用于所有 Topic，保持同一 spot 内的组成闭合与 Topic 间相关结构。
  zp <- z[sample.int(n), , drop = FALSE]
  weighted <- adjacency %*% zp
  perm_moran[r, ] <- n / S0 * colSums(zp * weighted) / denominator
  perm_pair[r, ] <- (crossprod(zp, weighted) / S0)[pair_index]
}
moran_p <- vapply(seq_along(moran), function(i)
  (1 + sum(abs(perm_moran[, i] - expected_moran) >= abs(moran[i] - expected_moran))) /
    (n_permutations + 1), numeric(1))
pair_p <- vapply(seq_along(observed_pair), function(i)
  (1 + sum(abs(perm_pair[, i] - null_pair_mean[i]) >= abs(observed_pair[i] - null_pair_mean[i]))) /
    (n_permutations + 1), numeric(1))
moran_result <- data.frame(topic = colnames(theta), Moran_I = unname(moran),
  null_expectation = expected_moran, permutation_mean = colMeans(perm_moran),
  pval_two_sided = moran_p, padj_BH = p.adjust(moran_p, "BH"), n_spots = n)
pair_result <- data.frame(topic_a = colnames(theta)[pair_index[, 1]],
  topic_b = colnames(theta)[pair_index[, 2]], H = observed_pair,
  null_expectation = null_pair_mean, permutation_mean = colMeans(perm_pair),
  permutation_sd = apply(perm_pair, 2, stats::sd), pval_two_sided = pair_p,
  padj_BH = p.adjust(pair_p, "BH"))
utils::write.csv(moran_result, file.path(out, "Moran_I_permutation.csv"), row.names = FALSE)
utils::write.csv(pair_result, file.path(out, "Topic_neighbor_permutation.csv"), row.names = FALSE)
saveRDS(list(moran = moran_result, pairs = pair_result, adjacency = adjacency,
  moran_null = perm_moran, pair_null = perm_pair, params = params), file.path(out, "spatial_statistics.rds"))

# ---- 距离衰减：分箱后的同 Topic 标准化乘积 / Descriptive distance analysis ----
all_pairs <- which(upper.tri(distance), arr.ind = TRUE)
pair_distance <- distance[all_pairs]
breaks <- seq(0, max(pair_distance) * (1 + 1e-10), length.out = distance_bins + 1L)
distance_group <- cut(pair_distance, breaks = breaks, include.lowest = TRUE)
distance_summary <- do.call(rbind, lapply(seq_len(ncol(z)), function(j) {
  product <- z[all_pairs[, 1], j] * z[all_pairs[, 2], j]
  do.call(rbind, lapply(levels(distance_group), function(level) {
    ix <- which(distance_group == level)
    data.frame(topic = colnames(z)[j], distance_bin = level, n_spot_pairs = length(ix),
               mean_distance = if (length(ix)) mean(pair_distance[ix]) else NA_real_,
               mean_standardized_product = if (length(ix)) mean(product[ix]) else NA_real_)
  }))
}))
utils::write.csv(distance_summary, file.path(out, "descriptive_distance_decay.csv"), row.names = FALSE)

# ---- 白底图件 / Figures ----
node <- data.frame(x = pos[, 1], y = pos[, 2], degree = rowSums(adjacency))
edge_plot <- data.frame(x = pos[edge_index[, 1], 1], y = pos[edge_index[, 1], 2],
                        xend = pos[edge_index[, 2], 1], yend = pos[edge_index[, 2], 2])
graph <- ggplot2::ggplot() + ggplot2::geom_segment(data = edge_plot,
  ggplot2::aes(x, y, xend = xend, yend = yend), colour = "grey80", linewidth = 0.25) +
  ggplot2::geom_point(data = node, ggplot2::aes(x, y, colour = degree), size = 1.5) +
  ggplot2::scale_colour_gradient(low = "#A6BDD3", high = "#2166AC") + ggplot2::coord_equal() +
  ggplot2::labs(title = paste0("Symmetric ", k_neighbors, "-nearest-neighbor graph"), subtitle = "Exploratory graph; original array coordinates", x = "Array x", y = "Array y")
save_figure(graph, file.path(fig, "07A_spatial_knn_graph"), width = 7, height = 6)
moran_plot <- ggplot2::ggplot(moran_result, ggplot2::aes(topic, Moran_I)) +
  ggplot2::geom_hline(yintercept = expected_moran, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(colour = padj_BH < q_cutoff), size = 3) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#2166AC"), name = paste("BH q <", q_cutoff)) +
  ggplot2::labs(title = "Topic spatial autocorrelation", subtitle = paste0("Moran's I; ", n_permutations, " joint random-label permutations; two-sided p; BH over all Topics"),
    x = "Topic", y = "Moran's I")
save_figure(moran_plot, file.path(fig, "07B_Moran_I"), width = 8, height = 5)
pair_plot <- rbind(pair_result, transform(pair_result[pair_result$topic_a != pair_result$topic_b, ],
                                         topic_a = topic_b, topic_b = topic_a))
pair_plot$star <- ifelse(pair_plot$padj_BH < q_cutoff, "*", "")
heat <- ggplot2::ggplot(pair_plot, ggplot2::aes(topic_a, topic_b, fill = H)) +
  ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = star), size = 5) +
  ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B35858", midpoint = 0) +
  ggplot2::coord_equal() + ggplot2::labs(title = "Cross-Topic spatial-neighbor association",
    subtitle = paste0("* BH q < ", q_cutoff, " across unique pairs incl. diagonal; two-sided joint random-label test"),
    x = "Topic at spot i", y = "Topic at neighboring spot j")
save_figure(heat, file.path(fig, "07C_Topic_neighbor_association"), width = 8, height = 7)
decay <- ggplot2::ggplot(distance_summary, ggplot2::aes(mean_distance, mean_standardized_product)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
  ggplot2::geom_line(colour = "#547B99") + ggplot2::geom_point(ggplot2::aes(size = n_spot_pairs), colour = "#547B99") +
  ggplot2::facet_wrap(~topic, ncol = 3) +
  ggplot2::labs(title = "Descriptive Topic association across spatial distances",
    subtitle = "Spot pairs are dependent; curves carry no independent-pair significance test",
    x = "Mean array distance (not micrometers)", y = "Mean standardized same-Topic product", size = "Spot pairs")
save_figure(decay, file.path(fig, "07D_descriptive_distance_decay"), width = 10, height = 7)
record_step("07_Spatial_Neighborhood", start_time,
  sprintf("%d spots; k=%d symmetric union; %d permutations; unconditioned random-label null; exploratory", n, k_neighbors, n_permutations))
