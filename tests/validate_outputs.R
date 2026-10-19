# 目的：验证真实输出的矩阵合同、ID 对齐、空间统计与导出完整性，不重跑模型。
# 回归点：mOB 的 annot 只覆盖260而非原始262 spots；pos顺序不同，必须按ID匹配。
# 运行：Rscript --vanilla tests/validate_outputs.R（仓库根目录）。
source("00_Start.R")
for (script in list.files(pattern = "^[0-9].*\\.R$")) parse(script)
core <- readRDS("data/processed/core.rds")
data("mOB", package = "STdeconvolve")
stopifnot(ncol(mOB$counts) == 262, nrow(mOB$pos) == 262, length(mOB$annot) == 260,
          !identical(colnames(mOB$counts), rownames(mOB$pos)),
          identical(names(core$annot), rownames(core$pos)),
          identical(rownames(core$theta_raw), rownames(core$pos)),
          identical(colnames(core$theta_raw), rownames(core$beta)),
          identical(unname(core$annot), unname(mOB$annot[rownames(core$pos)])),
          all(as.matrix(core$pos) == as.matrix(mOB$pos[rownames(core$theta_raw), ])),
          all(is.finite(core$theta_raw)), all(core$theta_raw >= 0),
          all(core$beta >= 0), max(abs(rowSums(core$theta_raw) - 1)) < 1e-8,
          max(abs(rowSums(core$theta) - 1)) < 1e-8,
          max(abs(rowSums(core$beta) - cfg$beta_scale)) < 1e-6)
paper <- readRDS("data/processed/paper_mOB.rds")
stopifnot(nrow(paper$corpus) == 260, ncol(paper$corpus) == 255,
          ncol(paper$theta_raw) == 12, max(abs(rowSums(paper$theta_raw) - 1)) < 1e-8)
ann <- read.csv("results/annotation/topic_annotation_review.csv")
stopifnot(all(ann$proposed_cell_type == "UNASSIGNED"))
sp <- readRDS("results/neighborhood/spatial_statistics.rds")
stopifnot(isSymmetric(sp$adjacency), all(diag(sp$adjacency) == 0),
          all(rowSums(sp$adjacency) > 0), nrow(sp$moran_null) == 999,
          all(sp$moran$pval_two_sided >= 1/1000),
          all(sp$moran$padj_BH >= sp$moran$pval_two_sided - 1e-12))
# 直接公式与保存的 Moran I 比较；避免归一化/方向错误。
theta <- core$theta_raw
centered <- sweep(theta, 2, colMeans(theta), "-")
manual_moran <- nrow(theta) / sum(sp$adjacency) *
  colSums(centered * (sp$adjacency %*% centered)) / colSums(centered^2)
stopifnot(max(abs(manual_moran - sp$moran$Moran_I)) < 1e-10)
ora <- read.csv("results/enrichment/ORA_all_results.csv")
stopifnot(all(ora$overlap <= ora$marker_n), all(ora$set_size <= ora$background_n),
          all(ora$padj >= ora$pval - 1e-12), all(ora$padj <= 1))
ms <- read.csv("results/multi_sample/LSAP_topic_pairs.csv")
stopifnot(!anyDuplicated(ms$Rep8_topic), !anyDuplicated(ms$Rep2_topic),
          all(is.finite(ms$pearson_r)))
pngs <- list.files("figures", pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
stopifnot(length(pngs) > 0)
for (p in pngs) {
  stem <- sub("\\.png$", "", p)
  paths <- paste0(stem, c(".png", ".pdf", ".svg"))
  stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 1000),
            identical(readBin(p, "raw", n = 8), as.raw(c(137,80,78,71,13,10,26,10))))
}
writeLines(paste("PASS: ID alignment, theta/beta, paper255/K12, annotation boundaries,",
                 "Moran formula, ORA, one-to-one matching and", length(pngs), "complete figure triples."),
           "environment/validation_results.txt")
cat(readLines("environment/validation_results.txt"), "\n")
