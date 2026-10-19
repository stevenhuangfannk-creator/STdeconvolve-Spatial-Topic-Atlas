# 目的：按官方 quick tutorial 用公开 mOB 数据复现预处理、K=2:9 拟合及 theta/beta 提取。
# 输入：STdeconvolve::mOB（原始整数 counts：gene × spot；pos；组织层 annot）。
# 输出：data/raw、data/processed/core.rds、results/models、results/theta_beta。
# 科学边界：这是官方快速教程复现，不是论文 Figure 4 的 K=12 全参数复现。
# 运行：在仓库根目录的 VS Code R 会话 source("01_Official_Reproduction.R")。
# 依赖：00_Start.R 检查 STdeconvolve；实际版本见 environment/sessionInfo.txt。

# ---- 集中可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
min_lib_size <- 100              # spot 总 counts > 100；官方 cleanCounts 使用严格大于。
min_reads <- 1                   # gene 跨 spot 总 counts > 1。
min_detected <- 1                # gene 在 > 1 个 spot 中检出。
remove_above <- 1.0              # 去掉在 >= 100% spots 检出的 gene。
remove_below <- 0.05             # 去掉在 <= 5% spots 检出的 gene。
od_alpha <- 0.05                 # 过度离散基因的 BH 校正显著性阈值。
n_top_od <- 1000                 # 官方默认，实际少于此值时保留全部过度离散基因。

# ---- 官方数据读取与检查 / Input checks ----
data("mOB", package = "STdeconvolve")
stopifnot(all(c("counts", "pos", "annot") %in% names(mOB)))
raw_counts <- mOB$counts
stopifnot(!is.null(rownames(raw_counts)), !is.null(colnames(raw_counts)),
          !anyDuplicated(rownames(raw_counts)), !anyDuplicated(colnames(raw_counts)),
          all(is.finite(raw_counts)), all(raw_counts >= 0),
          all(raw_counts == floor(raw_counts)),
          all(c("x", "y") %in% colnames(mOB$pos)),
          all(colnames(raw_counts) %in% rownames(mOB$pos)),
          !is.null(names(mOB$annot)), !anyDuplicated(names(mOB$annot)),
          all(names(mOB$annot) %in% colnames(raw_counts)))
# 官方原始 counts/pos 有 262 spots，而组织层 annot 只标注其中 260 个。
# 坐标顺序也可能与 counts 不同；所有对齐都使用已有 spot ID，不按位置重命名。
raw_annot <- mOB$annot
saveRDS(mOB, "data/raw/mOB_official.rds")
utils::write.csv(data.frame(annotation_source_index = seq_along(raw_annot),
                            count_source_index = match(names(raw_annot), colnames(raw_counts)),
                            spot_id = names(raw_annot),
                            tissue_layer = as.character(raw_annot)),
                 "data/processed/annotation_source_index.csv", row.names = FALSE)

# ---- 官方预处理 / Official preprocessing ----
counts <- STdeconvolve::cleanCounts(raw_counts, min.lib.size = min_lib_size,
                                    min.reads = min_reads, min.detected = min_detected,
                                    verbose = TRUE, plot = FALSE)
corpus <- STdeconvolve::restrictCorpus(counts, removeAbove = remove_above,
                                       removeBelow = remove_below, alpha = od_alpha,
                                       nTopOD = n_top_od, plot = FALSE, verbose = TRUE)
stopifnot(nrow(corpus) >= 2, ncol(corpus) >= 2,
          all(Matrix::colSums(corpus) > 0), all(colnames(corpus) %in% names(raw_annot)))
pos <- as.data.frame(mOB$pos[colnames(corpus), c("x", "y"), drop = FALSE])
annot <- raw_annot[colnames(corpus)]
stopifnot(identical(rownames(pos), colnames(corpus)),
          identical(names(annot), colnames(corpus)), all(is.finite(as.matrix(pos))))
qc <- data.frame(stage = c("official_raw", "cleanCounts", "restrictCorpus"),
                 genes = c(nrow(raw_counts), nrow(counts), nrow(corpus)),
                 spots = c(ncol(raw_counts), ncol(counts), ncol(corpus)))
utils::write.csv(qc, "data/processed/preprocessing_summary.csv", row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(pos), pos, tissue_layer = as.character(annot)),
                 "data/processed/spot_metadata.csv", row.names = FALSE)
saveRDS(list(counts = counts, corpus = corpus, pos = pos, annot = annot),
        "data/processed/official_inputs.rds")

# ---- 官方 LDA / Official model fitting ----
# fitLDA 接收 spot × gene 的整数计数；串行、固定随机种子。
# 上游 fitLDA 的 fitCorpus 与 testCorpus 相同；perplexity 是拟合语料内诊断。
ldas <- STdeconvolve::fitLDA(t(as.matrix(corpus)), Ks = cfg$Ks, seed = cfg$seed,
                            perc.rare.thresh = cfg$rare, ncores = cfg$ncores,
                            plot = FALSE, verbose = TRUE)
saveRDS(ldas, "results/models/ldas.rds")
optLDA <- STdeconvolve::optimalModel(models = ldas, opt = "min")
selected_k <- methods::slot(optLDA, "k")
raw_result <- STdeconvolve::getBetaTheta(optLDA, perc.filt = 0,
                                         betaScale = cfg$beta_scale)
display_result <- STdeconvolve::getBetaTheta(optLDA, perc.filt = cfg$perc_filt,
                                             betaScale = cfg$beta_scale)
theta_raw <- raw_result$theta
theta <- display_result$theta
beta <- raw_result$beta
# perc.filt=0.05 会将低比例置零、重新归一化并可能丢弃所有 spot 都为零的 Topic。
# theta_raw 供统计，theta 仅用于与官方教程一致的显示；beta 为概率 × 1000。
stopifnot(identical(rownames(theta_raw), rownames(pos)),
          identical(colnames(theta_raw), rownames(beta)),
          max(abs(rowSums(theta_raw) - 1)) < 1e-8,
          max(abs(rowSums(beta) - cfg$beta_scale)) < 1e-6)
core <- list(counts = counts, corpus = corpus, pos = pos, annot = annot,
             theta = theta, theta_raw = theta_raw, beta = beta, optLDA = optLDA,
             selected_k = selected_k, seed = cfg$seed,
             beta_scale = cfg$beta_scale, display_threshold = cfg$perc_filt)
saveRDS(core, "data/processed/core.rds")
utils::write.csv(data.frame(spot_id = rownames(theta_raw), theta_raw, check.names = FALSE),
                 "results/theta_beta/theta_raw.csv", row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(theta), theta, check.names = FALSE),
                 "results/theta_beta/theta_display_filtered.csv", row.names = FALSE)
utils::write.csv(data.frame(topic = rownames(beta), beta, check.names = FALSE),
                 "results/theta_beta/beta_scaled.csv", row.names = FALSE)
utils::write.csv(data.frame(topic = rownames(beta), beta / cfg$beta_scale, check.names = FALSE),
                 "results/theta_beta/beta_probability.csv", row.names = FALSE)
record_step("01_Official_Reproduction", start_time,
            sprintf("official mOB; %d genes x %d spots; selected K=%d; in-sample perplexity", nrow(corpus), ncol(corpus), selected_k))
message("官方快速教程核心完成；选择 K = ", selected_k,
        "。Topic 尚未经过细胞类型参考验证。")
