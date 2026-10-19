# 目的：统一项目参数、独立包库、依赖检查与图件导出。
# 输入：从本仓库根目录启动 R / VS Code R Interactive。
# 输出：项目内目录、参数快照、真实 sessionInfo；不自动安装任何包。
# 科学边界：本项目只使用官方公开数据；Topic 编号不代表已确认的细胞类型。

# ---- 可调参数 / Parameters ----
cfg <- list(
  root = normalizePath(Sys.getenv("STDECONVOLVE_PROJECT", unset = getwd()), winslash = "/"),
  seed = 0L, Ks = 2:9, ncores = 1L, rare = 0.05,
  perc_filt = 0.05, beta_scale = 1000, repeats = c(1L, 2L, 3L),
  dpi = 300, font = "Arial",
  upstream_sha = "13b1953263a5ce34e2471c4defebedbda44c7d08"
)

# ---- 固定步骤 / Fixed steps ----
stopifnot(file.exists(file.path(cfg$root, "00_Start.R")))
setwd(cfg$root)
Sys.setenv(TZ = "Asia/Singapore", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
dirs <- c("data/raw", "data/processed", "results/models", "results/theta_beta",
          "results/annotation", "results/enrichment", "results/multi_sample",
          "results/neighborhood", "figures/official_reproduction",
          "figures/extended_analysis", "environment/R-library", "environment/logs")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(cfg$root, "environment/R-library"), .libPaths()))
required <- c("STdeconvolve", "topicmodels", "BiocParallel", "ggplot2",
              "patchwork", "svglite", "ragg", "clue")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("缺少依赖：", paste(missing, collapse = ", "),
                          "。请先阅读 environment/install_dependencies.R，不会自动安装。")
set.seed(cfg$seed)
BiocParallel::register(BiocParallel::SerialParam())
if (!interactive()) options(device = function(...) grDevices::pdf(file = NULL))
ggplot2::theme_set(ggplot2::theme_classic(base_size = 11, base_family = cfg$font) +
                    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA)))

# 用途：同一图对象导出三个格式；输入为 ggplot / patchwork / grid grob，stem 不含扩展名。
save_figure <- function(plot, stem, width = 8, height = 6) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height,
                 units = "in", dpi = cfg$dpi, bg = "white", device = ragg::agg_png)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                 units = "in", bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width, height = height,
                 units = "in", bg = "white", device = svglite::svglite)
  if (interactive()) print(plot)
}

# 用途：记录脚本耗时与 R 堆内存高水位；不是操作系统 RSS 或整机峰值。
record_step <- function(step, start_time, notes = "") {
  mem <- gc()
  row <- data.frame(step = step, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                    elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
                    r_heap_max_mb = sum(mem[, 6]), notes = notes)
  log <- "environment/run_metrics.csv"
  utils::write.table(row, log, sep = ",", row.names = FALSE,
                     col.names = !file.exists(log), append = file.exists(log), qmethod = "double")
}

# 参数快照不含机器私人路径；重复 source 不会清除已加载分析对象。
saveRDS(cfg[names(cfg) != "root"], "environment/parameters.rds")
capture.output(dput(cfg[names(cfg) != "root"]), file = "environment/parameters.txt")
capture.output(sessionInfo(), file = "environment/sessionInfo.txt")
ip <- installed.packages()
utils::write.csv(ip[, c("Package", "Version", "Built"), drop = FALSE],
                 "environment/installed_packages.csv", row.names = FALSE)
message("项目就绪：", R.version.string, "; STdeconvolve ", packageVersion("STdeconvolve"),
        "; ncores=", cfg$ncores)
