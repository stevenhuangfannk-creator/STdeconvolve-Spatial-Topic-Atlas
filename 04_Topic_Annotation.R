# 目的：复现官方组织层代理相关性/GSEA，输出可审查的 Topic 候选注释。
# 输入：01 的 counts、beta、theta_raw、mOB 组织层 annot；无独立单细胞参考。
# 输出：results/annotation 的 markers、相关性、GSEA、候选注释表及 Figure 04。
# 科学边界：组织层不是细胞类型；同一切片生成代理参考，属于循环参照示例，不能当独立验证。
# 运行：仓库根目录 source("04_Topic_Annotation.R")。
# 依据：固定源码 vignette 的 Annotation Strategy 1 / 2；安全索引、minSize=15 和串行为项目调整。

# ---- 集中可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
n_top_genes <- 10                # 候选审查表的高 beta genes 数。
proxy_min_counts <- 10           # 官方教程：组织层聚合 counts > 10。
proxy_min_log2fc <- 1            # 官方教程：相对其他层平均聚合 counts 的 log2FC > 1。
proxy_top_n <- 200               # 官方教程：最多取 log2FC 排序前 200 个基因。
gsea_qval <- 0.05                # fgsea 校正 p 值阈值，非独立参考置信度。
gsea_exponent <- 0.1             # 固定版本 vignettes/vignette.Rmd 的 gseaParam，网页示例可能不同。
gsea_min_size <- 15              # 项目显式调整，非 fgsea 默认；避免很小的层代理基因集。

# ---- 固定步骤：组织层代理参考 / Official tissue-layer proxy ----
stopifnot(file.exists("data/processed/core.rds"))
core <- readRDS("data/processed/core.rds")
stopifnot(identical(names(core$annot), rownames(core$theta_raw)))
annot <- droplevels(core$annot)
mobProxyTheta <- stats::model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
colnames(mobProxyTheta) <- levels(annot)
mobProxyGexp <- core$counts[, rownames(mobProxyTheta), drop = FALSE] %*% mobProxyTheta
saveRDS(list(theta = mobProxyTheta, gene_counts = mobProxyGexp,
              source = "Same-slice mOB annotated tissue-layer aggregate; not scRNA-seq"),
        "results/annotation/official_layer_proxy.rds")

# 官方 Pearson 相似性；一个层可以关联多个 Topic，不强制将 K 个 Topic 唯一映射到 5 层。
beta_cor <- STdeconvolve::getCorrMtx(as.matrix(core$beta), t(as.matrix(mobProxyGexp)), type = "b")
theta_cor <- STdeconvolve::getCorrMtx(core$theta_raw, mobProxyTheta, type = "t")
utils::write.csv(data.frame(topic = rownames(beta_cor), beta_cor, check.names = FALSE),
                 "results/annotation/topic_layer_beta_correlations.csv", row.names = FALSE)
utils::write.csv(data.frame(topic = rownames(theta_cor), theta_cor, check.names = FALSE),
                 "results/annotation/topic_layer_theta_correlations.csv", row.names = FALSE)
# ggplot2 兼容：仅给显示矩阵加文字前缀，避免官方 melt 将纯数字 Topic 名转成连续轴，
# 而 annotation 的字符坐标落在同一连续轴时报错；原始相关系数和保存表不变。
beta_display <- beta_cor
theta_display <- theta_cor
rownames(beta_display) <- paste("Topic", rownames(beta_display))
rownames(theta_display) <- paste("Topic", rownames(theta_display))
p_beta <- STdeconvolve::correlationPlot(beta_display, colLabs = "Topic", rowLabs = "Tissue-layer proxy",
                                        title = "Transcriptional Pearson correlation", annotation = TRUE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
p_theta <- STdeconvolve::correlationPlot(theta_display, colLabs = "Topic", rowLabs = "Annotated tissue layer",
                                         title = "Spatial-fraction Pearson correlation", annotation = TRUE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
save_figure(p_beta, "figures/official_reproduction/Figure04_Layer_Beta_Correlation", 9, 5)
save_figure(p_theta, "figures/official_reproduction/Figure04_Layer_Theta_Correlation", 9, 5)

# ---- 官方 GSEA 层代理注释 / Official layer-marker demonstration ----
# 这里的基因集从真实官方 counts 和已给出的组织层标签生成，不手工创造细胞类型 marker。
gexp <- t(as.matrix(mobProxyGexp))
mobProxyLayerMarkers <- list()
for (i in seq_len(nrow(gexp))) {
  high_genes <- colnames(gexp)[gexp[i, ] > proxy_min_counts]
  log2fc <- sort(log2(gexp[i, high_genes] /
                       colMeans(gexp[-i, high_genes, drop = FALSE])), decreasing = TRUE)
  log2fc <- head(log2fc, proxy_top_n)
  mobProxyLayerMarkers[[rownames(gexp)[i]]] <- names(log2fc)[log2fc > proxy_min_log2fc]
}
saveRDS(mobProxyLayerMarkers, "results/annotation/official_proxy_layer_markers.rds")
marker_coverage <- data.frame(layer = names(mobProxyLayerMarkers),
                              markers_from_full_counts = lengths(mobProxyLayerMarkers),
                              markers_in_beta_corpus = vapply(mobProxyLayerMarkers, function(g) sum(g %in% colnames(core$beta)), integer(1)))
utils::write.csv(marker_coverage, "results/annotation/proxy_marker_coverage.csv", row.names = FALSE)
set.seed(cfg$seed)
layer_gsea <- STdeconvolve::annotateCellTypesGSEA(
  beta = core$beta, gset = mobProxyLayerMarkers, qval = gsea_qval,
  gseaParam = gsea_exponent, minSize = gsea_min_size,
  BPPARAM = BiocParallel::SerialParam(), nproc = 1)
saveRDS(layer_gsea, "results/annotation/official_layer_proxy_gsea.rds")
gsea_rows <- lapply(names(layer_gsea$results), function(id) {
  x <- as.data.frame(layer_gsea$results[[id]])
  if (!nrow(x)) return(NULL)
  if ("leadingEdge" %in% names(x)) x$leadingEdge <- vapply(x$leadingEdge, paste, collapse = ";", character(1))
  data.frame(topic = id, x, check.names = FALSE)
})
gsea_table <- do.call(rbind, gsea_rows)
if (is.null(gsea_table)) gsea_table <- data.frame(topic = character(), pathway = character(), padj = numeric(), NES = numeric())
utils::write.csv(gsea_table, "results/annotation/official_layer_proxy_gsea_significant.csv", row.names = FALSE)

# ---- 候选注释审查表 / Reviewable candidate table ----
top <- STdeconvolve::topGenes(core$beta, n = min(n_top_genes, ncol(core$beta)))
top_layer_beta <- max.col(beta_cor, ties.method = "first")
top_layer_theta <- max.col(theta_cor, ties.method = "first")
candidate <- data.frame(topic = rownames(core$beta),
                        high_beta_genes = vapply(top, function(x) paste(names(x), collapse = ";"), character(1)),
                        best_beta_layer_proxy = colnames(beta_cor)[top_layer_beta],
                        beta_pearson = beta_cor[cbind(seq_len(nrow(beta_cor)), top_layer_beta)],
                        best_spatial_layer = colnames(theta_cor)[top_layer_theta],
                        theta_pearson = theta_cor[cbind(seq_len(nrow(theta_cor)), top_layer_theta)],
                        gsea_positive_layers = vapply(rownames(core$beta), function(id) {
                          x <- as.data.frame(layer_gsea$results[[id]])
                          if (!nrow(x)) return("")
                          paste(x$pathway[x$ES > 0], collapse = ";")
                        }, character(1)),
                        proposed_cell_type = "UNASSIGNED",
                        annotation_status = "CANDIDATE_LAYER_ASSOCIATION_ONLY",
                        review_note = "同切片组织层代理；不能证明独立细胞类型，需要独立参考和实验审查。")
utils::write.csv(candidate, "results/annotation/topic_annotation_review.csv", row.names = FALSE)
layer_means <- t(vapply(levels(annot), function(layer) {
  colMeans(core$theta_raw[annot == layer, , drop = FALSE])
}, numeric(ncol(core$theta_raw))))
colnames(layer_means) <- colnames(core$theta_raw)
utils::write.csv(data.frame(layer = rownames(layer_means), layer_means, check.names = FALSE),
                 "results/annotation/raw_theta_means_by_layer.csv", row.names = FALSE)
status <- data.frame(analysis = c("official_tissue_layer_correlation", "official_tissue_layer_proxy_GSEA", "independent_scRNA_cell_type_annotation"),
                      status = c("PASS", "PASS", "NOT_RUN"),
                      evidence_scope = c("同一切片真实组织层标签关联；非独立验证", "按官方教程从同切片聚合counts导出代理基因集；非独立验证", "未提供独立单细胞表达谱和可靠细胞类型标记；不生成细胞类型结论"))
utils::write.csv(status, "results/annotation/annotation_status.csv", row.names = FALSE)
record_step("04_Topic_Annotation", start_time,
            "Official same-slice layer proxy correlation/GSEA executed; independent cell-type annotation NOT_RUN")
