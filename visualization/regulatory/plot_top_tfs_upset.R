#!/usr/bin/env Rscript

# Summarise transcription factor activity across clusters; produce dot+bar plots and optional UpSet diagrams.
#
# Inputs
#   --input-dir: directory containing `TF_scores_long-format_*_(wmean|mlm).csv` (required).
#   --output-dir: directory for output PDFs (required).
#   --file-pattern: regex to match input files (default: `^TF_scores_long-format_.*_(wmean|mlm)\.csv$`).
#   --alpha: adjusted p-value threshold (default: 0.05).
#   --top-n: comma-separated list (e.g. `10,20`).
#   --clusters: cluster range/list for dot plot (default: `0:29`).
#   --make-upset: whether to emit UpSet plots (default: TRUE).
#   --width/--height: PDF size (default: 12 x 7).
#
# Dependencies: optparse, ggplot2, dplyr, tidyr, stringr, patchwork

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

comparison_from_path <- function(path) {
  nm <- basename(path)
  m <- regmatches(nm, regexec("^TF_scores_long-format_([^-_]+)_(wmean|mlm)\\.csv$", nm))[[1]]
  if (length(m) >= 2) m[2] else tools::file_path_sans_ext(nm)
}

pick_col <- function(nms, candidates) {
  ix <- which(tolower(nms) %in% tolower(candidates))
  if (length(ix)) nms[ix[1]] else NA_character_
}

metric_order <- c("mlm", "wmean", "wmean_norm", "wmean_corr")

read_tf_file <- function(path, alpha, fallback_p) {
  df <- read.csv(path, check.names = FALSE)
  if (!nrow(df)) stop("Empty TF file: ", basename(path))

  nms <- names(df)
  tf_col <- pick_col(nms, c("source", "tf"))
  cluster_col <- pick_col(nms, c("cluster", "cluster_id"))
  method_col <- pick_col(nms, c("statistic", "method", "metric"))
  score_col <- pick_col(nms, c("wmean", "wmean_norm", "norm_wmean", "wmean_corr", "corr_wmean", "score"))
  padj_col <- pick_col(nms, c("p.adjust", "padj", "qvalue", "fdr"))
  pval_col <- pick_col(nms, c("p_value", "pvalue", "p_val", "p.value"))
  rep_col <- pick_col(nms, c("replicate", "rep"))

  dat <- df %>% mutate(
    TF = .data[[tf_col]],
    cluster = as.integer(as.character(.data[[cluster_col]])),
    method = tolower(if (!is.na(method_col)) as.character(.data[[method_col]]) else NA_character_),
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_,
    padj = if (!is.na(padj_col)) suppressWarnings(as.numeric(.data[[padj_col]])) else NA_real_,
    pval = if (!is.na(pval_col)) suppressWarnings(as.numeric(.data[[pval_col]])) else NA_real_
  )

  if (!all(is.na(dat$padj))) {
    dat <- dat %>% mutate(sig = !is.na(padj) & padj <= alpha)
    sig_method <- paste0("padj ≤ ", alpha)
  } else if (!all(is.na(dat$pval))) {
    dat <- dat %>% mutate(sig = !is.na(pval) & pval <= fallback_p)
    sig_method <- paste0("p ≤ ", fallback_p)
  } else if (!all(is.na(dat$score))) {
    dat <- dat %>% mutate(sig = !is.na(score) & score != 0)
    sig_method <- "score ≠ 0 (fallback)"
  } else {
    stop("No significance metric in ", basename(path))
  }

  chosen_method <- intersect(metric_order, unique(na.omit(dat$method)))[1]
  dir_df <- dat %>% filter(is.na(chosen_method) | method == chosen_method) %>%
    group_by(TF, cluster) %>% summarise(signed_sum = sum(score, na.rm = TRUE), .groups = "drop")

  agg <- dat %>% group_by(TF, cluster) %>% summarise(
    sig = any(sig, na.rm = TRUE),
    agg_abs_stat = sum(abs(score), na.rm = TRUE),
    .groups = "drop"
  ) %>%
    left_join(dir_df, by = c("TF", "cluster")) %>%
    mutate(signed_sum = ifelse(is.na(signed_sum), 0, signed_sum), comparison = comparison_from_path(path))
  attr(agg, "sig_method") <- sig_method
  attr(agg, "direction_method") <- chosen_method
  agg
}

make_bar_dot <- function(dat, top_n, clusters_all) {
  term_counts <- dat %>% group_by(TF) %>% summarise(
    n_clusters = sum(sig, na.rm = TRUE),
    sum_abs_stat = sum(agg_abs_stat, na.rm = TRUE),
    n_pos = sum(sig & signed_sum > 0),
    n_neg = sum(sig & signed_sum < 0), .groups = "drop") %>%
    mutate(dir_class = case_when(n_pos > 0 & n_neg == 0 ~ "activated",
                                 n_neg > 0 & n_pos == 0 ~ "deactivated",
                                 TRUE ~ "mixed")) %>%
    filter(n_clusters > 0) %>% arrange(desc(n_clusters), desc(sum_abs_stat)) %>% slice_head(n = top_n)
  if (!nrow(term_counts)) return(list(plot = NULL, counts = term_counts))
  dat_top <- dat %>% filter(TF %in% term_counts$TF) %>%
    left_join(term_counts, by = "TF") %>% mutate(TF = reorder(TF, n_clusters),
                                                  cluster = factor(as.character(cluster), levels = clusters_all),
                                                  cluster_dir = case_when(sig & signed_sum > 0 ~ "up",
                                                                          sig & signed_sum < 0 ~ "down",
                                                                          TRUE ~ "none"))
  p_dot <- ggplot(dat_top, aes(x = cluster, y = TF)) +
    geom_point(aes(alpha = sig, colour = cluster_dir), size = 2.2) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
    scale_colour_manual(values = c(up = "firebrick4", down = "dodgerblue4", none = "grey70"), guide = "none") +
    labs(title = "Clusters per TF", x = "cluster", y = NULL) +
    theme_bw(base_size = 11)
  p_bar <- ggplot(term_counts, aes(x = reorder(TF, n_clusters), y = n_clusters, fill = dir_class)) +
    geom_col() + geom_text(aes(label = n_clusters), hjust = -0.1, size = 2.8) +
    scale_fill_manual(values = c(activated = "firebrick4", deactivated = "dodgerblue4", mixed = "grey50"), guide = "none") +
    coord_flip(clip = "off") +
    labs(title = "# clusters significant", x = NULL, y = "# clusters") +
    theme_bw(base_size = 11) + theme(plot.margin = margin(10, 20, 10, 10))
  list(plot = (p_dot | p_bar) + plot_layout(widths = c(1.25, 1)), counts = term_counts)
}

maybe_upset <- function(dat, top_terms, outfile, clusters_fixed) {
  if (!length(top_terms) || !nzchar(outfile) || !requireNamespace("UpSetR", quietly = TRUE)) return()
  pdf_open <- FALSE
  tryCatch({
    dat_top <- dat %>% filter(TF %in% top_terms)
    cl_all <- clusters_fixed
    memb <- dat_top %>% filter(sig, !is.na(cluster)) %>% mutate(cluster = as.integer(cluster)) %>%
      distinct(TF, cluster) %>% mutate(val = 1L) %>% pivot_wider(names_from = cluster, values_from = val, values_fill = 0L)
    if (!nrow(memb)) return()
    needed_cols <- as.character(cl_all)
    for (mc in setdiff(needed_cols, names(memb))) memb[[mc]] <- 0L
    memb <- memb[, c(setdiff(names(memb), needed_cols), needed_cols), drop = FALSE]
    df <- as.data.frame(memb[, needed_cols, drop = FALSE])
    if (!nrow(df) || sum(df, na.rm = TRUE) == 0) return()
    rownames(df) <- make.names(memb[[1]], unique = TRUE)
    colnames(df) <- paste0("cl_", colnames(df))
    pdf(outfile, width = 12, height = 6)
    pdf_open <- TRUE
    UpSetR::upset(df, sets = colnames(df), keep.order = TRUE,
                  empty.intersections = "on", mb.ratio = c(0.6, 0.4), sets.bar.color = "grey30",
                  mainbar.y.label = "TFs per cluster", sets.x.label = "TF count per cluster")
    dev.off()
    pdf_open <- FALSE
  }, error = function(e) {
    if (isTRUE(pdf_open)) try(dev.off(), silent = TRUE)
    if (file.exists(outfile)) unlink(outfile)
    warning("Skipping UpSet output for ", basename(outfile), ": ", conditionMessage(e))
  })
}

main_script <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Input directory"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--file-pattern", type = "character", default = "^TF_scores_long-format_.*_(wmean|mlm)\\.csv$"),
    make_option("--alpha", type = "double", default = 0.05),
    make_option("--top-n", type = "character", default = "10,20"),
    make_option("--clusters", type = "character", default = "0:29"),
    make_option("--make-upset", action = "store_true", default = TRUE),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 7),
    make_option("--upset-height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_top_tfs_upset.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) stop("--input-dir and --output-dir are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(args$`input-dir`, pattern = args$`file-pattern`, full.names = TRUE)
  if (!length(files)) stop("No TF files matched pattern")

  clusters_to_show <- if (grepl(":", args$clusters)) {
    bounds <- as.integer(strsplit(args$clusters, ":")[[1]])
    seq(bounds[1], bounds[2])
  } else {
    as.integer(trimws(strsplit(args$clusters, ",")[[1]]))
  }
  clusters_fixed <- clusters_to_show
  top_ns <- as.integer(trimws(strsplit(args$`top-n`, ",")[[1]]))

  for (path in files) {
    dat <- read_tf_file(path, args$alpha, fallback_p = args$alpha)  # fallback uses same alpha
    for (tn in top_ns) {
      res <- make_bar_dot(dat, tn, clusters_to_show)
      if (is.null(res$plot)) next
      pdf_file <- file.path(args$`output-dir`, sprintf("TF_top%02d_%s.pdf", tn, comparison_from_path(path)))
      n_terms <- nrow(res$counts)
      pdf_height <- max(5, min(14, n_terms * 0.45 + 2))
      pdf(pdf_file, width = args$width, height = pdf_height)
      print(res$plot)
      dev.off()
      if (isTRUE(args$`make-upset`)) {
        upset_file <- file.path(args$`output-dir`, sprintf("TF_top%02d_%s_UpSet.pdf", tn, comparison_from_path(path)))
        maybe_upset(dat, res$counts$TF, upset_file, clusters_fixed)
      }
    }
  }
}

if (sys.nframe() == 0L) {
  main_script()
}
