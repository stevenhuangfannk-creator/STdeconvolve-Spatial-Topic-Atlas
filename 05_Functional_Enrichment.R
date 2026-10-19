# 目的：解释潜在 Topic 的 GO BP、KEGG、Reactome 功能；本模块为项目扩展。
# 输入：01 的 core.rds；小鼠 SYMBOL；官方数据库或有来源记录的 mouse Entrez GMT。
# 输出：逐 Topic 的 ORA/fgsea 全表、映射与背景、数据库来源、4 类图件。
# 科学边界：所选过度离散基因背景很小，富集是条件性、探索性解释，不证明通路激活。
# 运行：核心复现通过后，在仓库根目录 source("05_Functional_Enrichment.R")。
# 依赖：AnnotationDbi、org.Mm.eg.db、GO.db、fgsea、digest；不自动安装。

# ---- 可调参数 / Parameters ----
source("00_Start.R")
start_time <- Sys.time()
ora_top_n <- 40L                   # 每 Topic 最多使用前 40 个高表达、正向特异基因。
marker_beta_min <- 5               # beta 概率 ×1000 的阈值，实际换算为 core$beta_scale。
rank_pseudocount <- 1e-8           # 概率单位；log2(beta/其余 Topic 均值) 的伪计数。
set_min_size <- 5L                 # 与实际背景交集后的通路最小基因数；小示例不使用默认15。
set_max_fraction <- 0.90           # 剔除覆盖背景 >90% 的过宽通路。
q_cutoff <- 0.05                   # BH 阈值；每个 Topic、每个库、每种检验分别校正。
plot_top_n <- 8L                   # 诊断图每库最多展示通路数，可能包含不显著项。
map_max_nodes <- 24L               # Enrichment Map 仅取显著 ORA 中 q 最小的通路。
map_jaccard <- 0.25                # 通路背景基因集的 Jaccard 连边阈值。
allow_download <- TRUE             # 只下载数据库文件到本项目；失败记录 BLOCKED。
reactome_mapping <- "data/raw/Reactome_mouse_NCBI.tsv" # 可预先放置官方六列 mouse 映射。
reactome_gmt <- ""                # 可选：mouse Entrez GMT；必须同时填写来源 URL/版本。
reactome_gmt_source <- ""         # GMT 物种应为 Mus musculus，禁止直接用人类 SYMBOL。
reactome_gmt_version <- ""

# ---- 输入与依赖检查 / Input checks ----
deps <- c("AnnotationDbi", "org.Mm.eg.db", "GO.db", "fgsea", "digest")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("05 缺少依赖：", paste(missing, collapse = ", "))
core <- readRDS("data/processed/core.rds")
beta_prob <- core$beta / core$beta_scale
stopifnot(nrow(beta_prob) > 1, !anyDuplicated(colnames(beta_prob)),
          all(is.finite(beta_prob)), all(beta_prob >= 0))
out <- "results/enrichment"
fig <- "figures/extended_analysis"
params <- list(ora_top_n = ora_top_n, marker_beta_min = marker_beta_min,
               rank_pseudocount = rank_pseudocount, min_size = set_min_size,
               max_fraction = set_max_fraction, q = q_cutoff, seed = cfg$seed,
               correction = "BH within topic/library/test", species = "Mus musculus")
saveRDS(params, file.path(out, "parameters.rds"))

# ---- SYMBOL -> Entrez 映射 / Identifier mapping ----
org <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
valid_symbols <- intersect(colnames(beta_prob), AnnotationDbi::keys(org, keytype = "SYMBOL"))
mapping <- unique(AnnotationDbi::select(org, keys = valid_symbols, keytype = "SYMBOL",
                                      columns = "ENTREZID"))
mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
# 多对一/一对多全部记录并排除，避免同一基因在超几何检验中被重复计数。
mapping$unique_symbol <- !duplicated(mapping$SYMBOL) & !duplicated(mapping$SYMBOL, fromLast = TRUE)
mapping$unique_entrez <- !duplicated(mapping$ENTREZID) & !duplicated(mapping$ENTREZID, fromLast = TRUE)
mapping$used <- mapping$unique_symbol & mapping$unique_entrez
utils::write.csv(mapping, file.path(out, "symbol_entrez_mapping.csv"), row.names = FALSE)
used <- mapping[mapping$used, , drop = FALSE]
if (nrow(used) < 2 * set_min_size) stop("可靠映射的背景太小，请检查物种/ID，不能替换成全基因组背景。")
b <- beta_prob[, used$SYMBOL, drop = FALSE]
colnames(b) <- used$ENTREZID
universe <- colnames(b)
provenance <- data.frame(library = character(), status = character(), source = character(),
                         version = character(), retrieved_at = character(), sha256 = character(), note = character())

# 用途：缓存官方小文件；输入 URL 与本项目目标路径，输出实际文件；不替换失败数据。
fetch_database <- function(url, dest) {
  if (!file.exists(dest)) {
    if (!allow_download) stop("下载关闭，缺少：", dest)
    utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  }
  if (!file.exists(dest) || file.info(dest)$size < 10) stop("数据库文件为空：", dest)
  dest
}

# ---- 三个独立数据库 / Gene-set sources ----
libraries <- list()
go_map <- AnnotationDbi::select(org, keys = universe, keytype = "ENTREZID",
                               columns = c("GOALL", "ONTOLOGYALL"))
go_map <- go_map[!is.na(go_map$GOALL) & go_map$ONTOLOGYALL == "BP", , drop = FALSE]
go_names <- AnnotationDbi::select(GO.db::GO.db, keys = unique(go_map$GOALL),
                                 keytype = "GOID", columns = "TERM")
libraries$GO_BP <- list(sets = split(go_map$ENTREZID, go_map$GOALL),
                        names = stats::setNames(go_names$TERM, go_names$GOID))
provenance <- rbind(provenance, data.frame(library = "GO_BP", status = "READY",
  source = "Bioconductor org.Mm.eg.db GOALL; GO.db (BP including ancestors)",
  version = paste(packageVersion("org.Mm.eg.db"), packageVersion("GO.db"), sep = ";"),
  retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  sha256 = digest::digest(list(go_map, go_names), algo = "sha256"), note = "Serialized mapping SHA256"))

for (lib in c("KEGG", "Reactome")) {
  result <- tryCatch({
    if (lib == "KEGG") {
      urls <- c("https://rest.kegg.jp/link/pathway/mmu", "https://rest.kegg.jp/list/pathway/mmu")
      files <- c("data/raw/KEGG_mmu_link.tsv", "data/raw/KEGG_mmu_pathway_names.tsv")
      for (i in seq_along(files)) fetch_database(urls[i], files[i])
      links <- utils::read.delim(files[1], header = FALSE, stringsAsFactors = FALSE)
      labels <- utils::read.delim(files[2], header = FALSE, stringsAsFactors = FALSE)
      sets <- split(sub("^mmu:", "", links[[1]]), sub("^path:", "", links[[2]]))
      term_names <- stats::setNames(sub(" - Mus musculus.*$", "", labels[[2]]),
                                    sub("^path:", "", labels[[1]]))
      source <- paste(urls, collapse = "; ")
      version <- "KEGG REST current; cached content hashes below"
    } else if (nzchar(reactome_gmt)) {
      if (!nzchar(reactome_gmt_source) || !nzchar(reactome_gmt_version))
        stop("使用 GMT 必须填写来源、版本并确认小鼠 Entrez ID。")
      lines <- strsplit(readLines(reactome_gmt, warn = FALSE), "\t", fixed = TRUE)
      if (!all(lengths(lines) >= 3)) stop("GMT 至少三列：名称、描述、基因…")
      sets <- stats::setNames(lapply(lines, function(x) x[-c(1, 2)]), vapply(lines, `[`, "", 1))
      term_names <- stats::setNames(names(sets), names(sets))
      source <- reactome_gmt_source
      version <- reactome_gmt_version
      files <- reactome_gmt
    } else {
      source <- "https://download.reactome.org/current/NCBI2Reactome_All_Levels.txt"
      version <- "Reactome current; Mus musculus, all hierarchy levels"
      if (!file.exists(reactome_mapping)) {
        full <- fetch_database(source, "data/raw/NCBI2Reactome_All_Levels.txt")
        # 逐块读取全物种文件，只保留小鼠；不一次载入大表，也不删除原始文件。
        input <- file(full, open = "rt")
        output <- file(reactome_mapping, open = "wt")
        repeat {
          chunk <- readLines(input, n = 50000L, warn = FALSE)
          if (!length(chunk)) break
          writeLines(grep("\tMus musculus$", chunk, value = TRUE), output)
        }
        close(input)
        close(output)
      }
      react <- utils::read.delim(reactome_mapping, header = FALSE, quote = "",
                                 stringsAsFactors = FALSE)
      if (ncol(react) != 6 || !all(react[[6]] == "Mus musculus"))
        stop("Reactome 必须是官方六列 NCBI 映射且第6列全为 Mus musculus。")
      sets <- split(as.character(react[[1]]), react[[2]])
      term_names <- stats::setNames(react[[4]], react[[2]])
      files <- reactome_mapping
    }
    if (!length(sets) || !any(universe %in% unlist(sets))) stop("数据库与小鼠 Entrez 背景没有交集。")
    libraries[[lib]] <- list(sets = sets, names = term_names)
    provenance <<- rbind(provenance, data.frame(library = lib, status = "READY", source = source,
      version = version, retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      sha256 = paste(vapply(files, digest::digest, "", file = TRUE, algo = "sha256"), collapse = ";"),
      note = "Cached file SHA256; an existing file retains its recorded source"))
    TRUE
  }, error = function(e) {
    provenance <<- rbind(provenance, data.frame(library = lib, status = "BLOCKED",
      source = if (lib == "KEGG") "https://rest.kegg.jp" else "https://reactome.org/download-data",
      version = NA_character_, retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      sha256 = NA_character_, note = conditionMessage(e)))
    message(lib, " BLOCKED: ", conditionMessage(e))
    FALSE
  })
}
utils::write.csv(provenance, file.path(out, "database_provenance.csv"), row.names = FALSE)

# ---- ORA 与有方向的排序富集 / Enrichment ----
ora_list <- gsea_list <- list()
coverage <- list()
for (lib in names(libraries)) {
  original_sets <- lapply(libraries[[lib]]$sets, unique)
  background <- intersect(universe, unique(unlist(original_sets)))
  sets <- lapply(original_sets, intersect, y = background)
  sets <- sets[lengths(sets) >= set_min_size & lengths(sets) <= floor(length(background) * set_max_fraction)]
  coverage[[lib]] <- data.frame(library = lib, beta_genes = ncol(beta_prob), mapped_genes = ncol(b),
                                annotated_background = length(background), tested_sets = length(sets))
  utils::write.csv(data.frame(entrez_id = background), file.path(out, paste0(lib, "_background.csv")), row.names = FALSE)
  saveRDS(sets, file.path(out, paste0(lib, "_tested_gene_sets.rds")))
  if (!length(sets)) next
  for (i in seq_len(nrow(b))) {
    other_mean <- colMeans(b[-i, , drop = FALSE])
    ranking <- log2((b[i, ] + rank_pseudocount) / (other_mean + rank_pseudocount))
    ranking <- sort(ranking[background], decreasing = TRUE)
    # ORA markers 按 beta 排序，且特异性 >0；所有被测通路均参与 BH，包括零 overlap。
    candidates <- names(sort(b[i, background], decreasing = TRUE))
    markers <- head(candidates[b[i, candidates] > marker_beta_min / 1000 & ranking[candidates] > 0], ora_top_n)
    utils::write.csv(data.frame(entrez_id = names(ranking), log2_topic_vs_other = unname(ranking),
                               ora_marker = names(ranking) %in% markers),
                     file.path(out, paste0(lib, "_Topic", rownames(b)[i], "_ranking.csv")), row.names = FALSE)
    overlap <- vapply(sets, function(x) length(intersect(x, markers)), integer(1))
    p <- stats::phyper(overlap - 1L, lengths(sets), length(background) - lengths(sets),
                       length(markers), lower.tail = FALSE)
    ora <- data.frame(library = lib, topic = rownames(b)[i], pathway = names(sets),
                       description = unname(libraries[[lib]]$names[names(sets)]),
                       overlap = overlap, set_size = lengths(sets), marker_n = length(markers),
                       background_n = length(background), pval = p, padj = p.adjust(p, "BH"),
                       overlap_entrez = vapply(sets, function(x) paste(intersect(x, markers), collapse = ";"), ""))
    ora$description[is.na(ora$description)] <- ora$pathway[is.na(ora$description)]
    ora$gene_ratio <- if (length(markers)) ora$overlap / length(markers) else 0
    ora_list[[paste(lib, i)]] <- ora
    set.seed(cfg$seed + i)
    gs <- as.data.frame(fgsea::fgseaMultilevel(pathways = sets, stats = ranking,
      minSize = set_min_size, maxSize = floor(length(background) * set_max_fraction),
      eps = 0, nproc = 1, BPPARAM = BiocParallel::SerialParam()))
    if (nrow(gs)) {
      gs$padj <- p.adjust(gs$pval, "BH")
      gs$leadingEdge <- vapply(gs$leadingEdge, paste, "", collapse = ";")
      gs$library <- lib
      gs$topic <- rownames(b)[i]
      gs$description <- unname(libraries[[lib]]$names[gs$pathway])
      gsea_list[[paste(lib, i)]] <- gs
    }
  }
}
coverage <- do.call(rbind, coverage)
utils::write.csv(coverage, file.path(out, "mapping_coverage.csv"), row.names = FALSE)
ora <- do.call(rbind, ora_list)
gsea <- do.call(rbind, gsea_list)
if (is.null(ora)) stop("没有可测试的通路，请检查 mapping_coverage.csv。")
utils::write.csv(ora, file.path(out, "ORA_all_results.csv"), row.names = FALSE)
if (!is.null(gsea)) utils::write.csv(gsea, file.path(out, "GSEA_all_results.csv"), row.names = FALSE)
saveRDS(list(ORA = ora, GSEA = gsea, coverage = coverage, provenance = provenance),
        file.path(out, "functional_enrichment.rds"))

# ---- 诊断图 / Diagnostic figures ----
# 同一库内取最小 q 的通路并展示所有 Topic；不会只保留“漂亮”的正结果。
selected <- unique(unlist(lapply(split(ora, ora$library), function(x) {
  x <- x[order(x$padj, x$pval, -x$overlap), ]
  head(unique(paste(x$library, x$pathway, sep = "::")), plot_top_n)
})))
show <- ora[paste(ora$library, ora$pathway, sep = "::") %in% selected, ]
show$label <- paste0(show$library, ": ", vapply(show$description, function(x) paste(strwrap(x, 45), collapse = "\n"), ""))
show$significant <- show$padj < q_cutoff
dot <- ggplot2::ggplot(show, ggplot2::aes(topic, label)) +
  ggplot2::geom_point(ggplot2::aes(size = overlap, colour = -log10(pmax(padj, 1e-300)), shape = significant)) +
  ggplot2::scale_colour_gradient(low = "grey80", high = "#2166AC", name = "-log10(BH q)") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = paste("BH q <", q_cutoff)) +
  ggplot2::labs(title = "ORA diagnostic: selected terms include non-significant results", x = "Topic", y = NULL, size = "Overlap genes")
save_figure(dot, file.path(fig, "05A_ORA_dotplot"), width = 11, height = max(6, length(unique(show$label)) * 0.3))
bar_data <- show[!duplicated(paste(show$library, show$pathway)), ]
bar_data <- do.call(rbind, lapply(split(show, interaction(show$library, show$pathway)), function(x) x[which.min(x$padj), ]))
bar <- ggplot2::ggplot(bar_data, ggplot2::aes(-log10(pmax(padj, 1e-300)), reorder(label, -padj))) +
  ggplot2::geom_col(fill = "#547B99") + ggplot2::geom_vline(xintercept = -log10(q_cutoff), linetype = 2) +
    ggplot2::labs(title = "ORA diagnostic: best Topic per term (selection shown)", subtitle = paste0("Dashed line: BH q = ", q_cutoff, "; full results in CSV"), x = "-log10(BH q)", y = NULL)
save_figure(bar, file.path(fig, "05B_ORA_barplot"), width = 10, height = max(6, nrow(bar_data) * 0.3))

sig <- ora[ora$padj < q_cutoff & ora$overlap > 0, ]
sig <- sig[order(sig$padj), ]
sig$key <- paste(sig$library, sig$pathway, sep = "::")
nodes <- head(sig[!duplicated(sig$key), ], map_max_nodes)
if (nrow(nodes)) {
  angle <- seq(0, 2 * pi, length.out = nrow(nodes) + 1)[seq_len(nrow(nodes))]
  nodes$x <- cos(angle); nodes$y <- sin(angle)
  edges <- data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric(), jaccard = numeric())
  if (nrow(nodes) > 1) for (i in seq_len(nrow(nodes) - 1L)) for (j in seq.int(i + 1L, nrow(nodes))) {
    if (nodes$library[i] != nodes$library[j]) next
    a <- libraries[[nodes$library[i]]]$sets[[nodes$pathway[i]]]
    z <- libraries[[nodes$library[j]]]$sets[[nodes$pathway[j]]]
    a <- intersect(a, universe); z <- intersect(z, universe)
    jac <- length(intersect(a, z)) / length(union(a, z))
    if (jac >= map_jaccard) edges <- rbind(edges, data.frame(x = nodes$x[i], y = nodes$y[i],
      xend = nodes$x[j], yend = nodes$y[j], jaccard = jac))
  }
  map <- ggplot2::ggplot() + ggplot2::geom_segment(data = edges,
    ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = jaccard), colour = "grey75") +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, size = overlap, colour = library)) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x * 1.13, y * 1.13, label = pathway), size = 3) +
    ggplot2::coord_equal(clip = "off", xlim = c(-1.4, 1.4), ylim = c(-1.3, 1.3)) +
    ggplot2::theme_void(base_family = cfg$font) +
    ggplot2::labs(title = "Significant ORA term map", subtitle = "Edges: within-library gene-set Jaccard; labels identify pathways in CSV")
} else {
  map <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = paste0("No ORA terms meet BH q < ", q_cutoff, "\nSmall selected-gene universe; see all-result diagnostics"), size = 4) +
    ggplot2::theme_void(base_family = cfg$font) + ggplot2::labs(title = "Enrichment Map: no significant nodes")
}
save_figure(map, file.path(fig, "05C_enrichment_map"), width = 10, height = 8)
if (!is.null(gsea) && nrow(gsea)) {
  gsea$key <- paste(gsea$library, gsea$pathway, sep = "::")
  choose <- unique(unlist(lapply(split(gsea, gsea$library), function(x) head(unique(x$key[order(x$padj, -abs(x$NES))]), plot_top_n))))
  gh <- gsea[gsea$key %in% choose, ]
  gh$label <- paste(gh$library, gh$pathway, sep = ": ")
  gh$star <- ifelse(!is.na(gh$padj) & gh$padj < q_cutoff, "*", "")
  heat <- ggplot2::ggplot(gh, ggplot2::aes(topic, label, fill = NES)) + ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 4) +
    ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B35858", midpoint = 0) +
    ggplot2::labs(title = "Ranked pathway enrichment (fgsea)", subtitle = paste0("* BH q < ", q_cutoff, " within Topic/library; unmarked cells are not significant"), x = "Topic", y = NULL)
  save_figure(heat, file.path(fig, "05D_GSEA_pathway_heatmap"), width = 9, height = max(6, length(unique(gh$label)) * 0.28))
}
record_step("05_Functional_Enrichment", start_time,
            paste("Exploratory; libraries:", paste(provenance$library, provenance$status, collapse = "; ")))
