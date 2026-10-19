# 目的：在项目独立包库安装缺少的依赖；不升级已有包，不改变全局 R 环境。
# 运行：在仓库根目录用选定的 Rscript --vanilla environment/install_dependencies.R
# 已验证环境：Windows / R 4.6.1；其他版本须先检查二进制与 Bioconductor 配对。
# 此脚本只在用户主动调用时安装。主分析脚本不调用它。
stopifnot(file.exists("00_Start.R"), getRversion() >= "4.1")
lib <- "environment/R-library"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 120)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib)
}
repos <- BiocManager::repositories()
deps <- c("topicmodels", "BiocParallel", "Matrix", "mgcv", "ggplot2", "scatterpie",
          "viridis", "slam", "clue", "reshape2", "fgsea", "patchwork", "svglite", "ragg",
          "AnnotationDbi", "digest")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  kind <- if (.Platform$OS.type == "windows") "binary" else "source"
  available <- available.packages(repos = repos, type = kind)
  if (!all(missing %in% rownames(available))) {
    stop("依赖包不可用；请停下审查环境，不自动升级 R 或改用其他算法。")
  }
  install.packages(missing, repos = repos, lib = lib, type = kind)
}
stopifnot(all(vapply(deps, requireNamespace, logical(1), quietly = TRUE)))
# GO 注释数据库是纯数据包；Bioconductor 不提供 Windows 二进制时直接安装源码数据。
# 不涉及编译，也不升级已有 AnnotationDbi；缺少兼容版本则停止审查。
annotation_pkgs <- c("org.Mm.eg.db", "GO.db")
missing_annotation <- annotation_pkgs[!vapply(annotation_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_annotation)) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("GO 扩展需要 AnnotationDbi；请先审查依赖后安装到独立包库。")
  }
  install.packages(missing_annotation, repos = unname(repos["BioCann"]), lib = lib, type = "source")
}
# 固定官方源码版本；纯 R 包本体可从源码安装，不需要改写官方函数。
sha <- "13b1953263a5ce34e2471c4defebedbda44c7d08"
archive <- tempfile(fileext = ".zip")
download.file(paste0("https://codeload.github.com/JEFworks-Lab/STdeconvolve/zip/", sha),
              archive, mode = "wb")
unpack <- tempfile()
dir.create(unpack)
unzip(archive, exdir = unpack)
source_dir <- list.dirs(unpack, recursive = FALSE, full.names = TRUE)
install.packages(source_dir, repos = NULL, type = "source", lib = lib)
stopifnot(requireNamespace("STdeconvolve", quietly = TRUE))
cat("安装完成：", as.character(packageVersion("STdeconvolve")), "\n")
