#!/usr/bin/env Rscript

# Generate transcription-factor bar plots across comparisons and clusters.
#
# Inputs
#   --input-dir: directory containing `TF_scores_long-format_*.csv` files (required).
#   --output-dir: directory for PDF outputs (required).
#   --metric: statistic name to plot (default: `mlm`).
#   --tfs: optional comma-separated list of TFs; defaults to intersection across comparisons.
#   --clusters: cluster range/list (default: `0:29`).
#   --order: comma-separated comparison order (default: `PvD,PvM,DvM`).
#   --width/--height: PDF size (default: 12 x 6).
#   --label-size: font size for p-value stars (default: 3).
#   --bar-outline: bar outline thickness (default: 0.25).
#
# Dependencies: optparse, ggplot2, dplyr, tidyr, readr, stringr, forcats

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

metric_aliases <- function(x) {
  switch(tolower(x),
         mlm = c("mlm"),
         wmean_norm = c("wmean_norm", "norm_wmean"),
         norm_wmean = c("wmean_norm", "norm_wmean"),
         wmean_corr = c("wmean_corr", "corr_wmean"),
         corr_wmean = c("wmean_corr", "corr_wmean"),
         wmean = c("wmean"),
         tolower(x))
}

choose_col <- function(nms, candidates) {
  ix <- which(tolower(nms) %in% tolower(candidates))
  if (length(ix)) nms[ix[1]] else NA_character_
}

p_to_stars <- function(p) {
  ifelse(is.na(p) | p > 5e-2, "", "*")
}

comparison_from_path <- function(path) {
  nm <- basename(path)
  m <- regmatches(nm, regexec("^TF_scores_long-format_([^-_]+).*\\.csv$", nm))[[1]]
  if (length(m) >= 2) m[2] else tools::file_path_sans_ext(nm)
}

read_one <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  if (!nrow(df)) return(NULL)
  nms <- names(df)
  tf_col <- choose_col(nms, c("source", "tf"))
  stat_col <- choose_col(nms, c("statistic", "method", "metric"))
  cond_col <- choose_col(nms, c("condition"))
  score_col <- choose_col(nms, c("wmean_norm", "norm_wmean", "wmean", "wmean_corr", "corr_wmean", "score", "activity", "estimate", "value"))
  pval_col <- choose_col(nms, c("p_value", "pvalue", "p_val"))
  cluster_col <- choose_col(nms, c("cluster", "cluster_id"))
  if (any(is.na(c(tf_col, stat_col, score_col, pval_col, cluster_col)))) stop("Missing columns in ", basename(path))
  tibble(
    TF = df[[tf_col]],
    statistic = as.character(df[[stat_col]]),
    score = suppressWarnings(as.numeric(df[[score_col]])),
    p_value = suppressWarnings(as.numeric(df[[pval_col]])),
    cluster = suppressWarnings(as.integer(df[[cluster_col]])),
    comparison = comparison_from_path(path)
  )
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Input directory"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--metric", type = "character", default = "mlm"),
    make_option("--tfs", type = "character", default = NULL),
    make_option("--clusters", type = "character", default = "0:29"),
    make_option("--order", type = "character", default = "PvD,PvM,DvM"),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 6),
    make_option("--label-size", type = "double", default = 3),
    make_option("--bar-outline", type = "double", default = 0.25)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_tf_activity_bars.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) stop("--input-dir and --output-dir are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(args$`input-dir`, pattern = "^TF_scores_long-format_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No input CSVs found")

  raw <- bind_rows(lapply(files, read_one))
  if (!nrow(raw)) stop("No data read from input files")

  metric_labels <- metric_aliases(args$metric)
  dat_metric <- raw %>% mutate(stat_lower = tolower(statistic)) %>% filter(stat_lower %in% metric_labels)

  dat_agg <- dat_metric %>%
    group_by(TF, comparison, cluster) %>%
    summarise(value = mean(score, na.rm = TRUE), pmin = min(p_value, na.rm = TRUE), .groups = "drop")

  comparisons <- intersect(trimws(strsplit(args$order, ",")[[1]]), unique(dat_metric$comparison))
  cluster_levels <- if (grepl(":", args$clusters)) {
    bounds <- as.integer(strsplit(args$clusters, ":")[[1]])
    as.integer(seq(bounds[1], bounds[2]))
  } else {
    as.integer(trimws(strsplit(args$clusters, ",")[[1]]))
  }

  tfs <- if (!is.null(args$tfs) && nzchar(args$tfs)) {
    trimws(strsplit(args$tfs, ",")[[1]])
  } else {
    shared_upper <- dat_metric %>% distinct(TF, comparison) %>% mutate(TF_upper = toupper(TF)) %>%
      group_by(TF_upper) %>% summarise(n_comp = n_distinct(comparison), .groups = "drop") %>%
      filter(n_comp == length(comparisons)) %>% pull(TF_upper)
    dat_metric %>% mutate(TF_upper = toupper(TF)) %>% filter(TF_upper %in% shared_upper) %>%
      arrange(TF_upper) %>% distinct(TF_upper, TF) %>% pull(TF)
  }
  if (!length(tfs)) stop("No TFs satisfy criteria for plotting")

  tf_comp <- tidyr::expand_grid(TF = tfs, comparison = comparisons)
  dat_complete <- tf_comp %>% tidyr::expand_grid(cluster = cluster_levels) %>%
    left_join(dat_agg, by = c("TF", "comparison", "cluster")) %>%
    mutate(comparison = factor(comparison, levels = comparisons),
           cluster = factor(cluster, levels = cluster_levels),
           value = ifelse(is.na(value), 0, as.numeric(value)),
           stars = vapply(pmin, p_to_stars, character(1)),
           dir = case_when(value > 0 ~ "up", value < 0 ~ "down", TRUE ~ "zero"))

  pal <- c(up = "firebrick4", down = "dodgerblue4", zero = "grey70")
  for (tf in tfs) {
    d <- dat_complete %>% filter(TF == tf)
    max_abs <- max(abs(d$value), na.rm = TRUE)
    y_lim <- if (is.finite(max_abs) && max_abs > 0) c(-max_abs, max_abs) * 1.05 else c(-1, 1)

    p <- ggplot(d, aes(x = cluster, y = value, fill = dir)) +
      geom_col(width = 0.75, color = "black", linewidth = args$`bar-outline`) +
      scale_fill_manual(values = pal, guide = "none") +
      geom_text(aes(label = stars, vjust = ifelse(value >= 0, -0.2, 1.2)), size = args$`label-size`, na.rm = TRUE) +
      facet_wrap(~ comparison, nrow = 1) +
      coord_cartesian(ylim = y_lim) +
      labs(title = paste0(tf, " - ", args$metric, " by cluster"), x = "Cluster", y = args$metric) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"), panel.grid.major.x = element_line(), panel.grid.minor.y = element_blank())

    ggsave(file.path(args$`output-dir`, paste0("TFbar_", tf, "_", args$metric, ".pdf")),
           p, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
