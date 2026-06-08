#!/usr/bin/env Rscript

# Plot cluster distributions along binned spatial axes with optional LOESS smoothing.
#
# Inputs
#   --input-dir: directory containing combined and per-axis binned CSVs plus cluster colour table (required).
#   --output-dir: destination for PDF plots (required).
#   --identity-groups: semicolon-separated definitions (e.g. `NN=0:7;GABA=8:19;Glut=20:29`).
#   --hex-map: CSV with columns `cluster` and `custom_hex` (default: `cluster-hex-colors.csv` in input).
#   --line-size: line width for raw line plots (default: 0.8).
#   --smooth-line-size: line width for LOESS plots (default: 0.8).
#   --loess-span: span parameter for LOESS smoothing (default: 0.18).
#
# Dependencies: optparse, readr, dplyr, tidyr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

parse_groups <- function(spec) {
  pieces <- strsplit(spec, ";")[[1]]
  groups <- lapply(pieces, function(piece) {
    kv <- strsplit(piece, "=")[[1]]
    name <- kv[1]
    range <- kv[2]
    if (grepl(":", range)) {
      bounds <- as.integer(strsplit(range, ":")[[1]])
      clusters <- seq(bounds[1], bounds[2])
    } else {
      clusters <- as.integer(strsplit(range, ",")[[1]])
    }
    setNames(list(as.character(clusters)), name)
  })
  Reduce(function(x, y) c(x, y), groups, init = list())
}

build_plots <- function(df, clusters, hex_colors, axis, identity_name, output_dir, line_size, smooth_line_size, loess_span) {
  df <- df %>% filter(cluster %in% clusters) %>% mutate(cluster = factor(cluster, levels = clusters))
  if (!nrow(df)) return()
  axis_min <- min(df$bin_lower, na.rm = TRUE)
  axis_max <- max(df$bin_lower, na.rm = TRUE)
  axis_breaks <- seq(axis_min, axis_max, by = 0.5)

  base_plot <- ggplot(df, aes(x = bin_lower, y = count, group = cluster, colour = cluster)) +
    scale_color_manual(values = hex_colors[clusters]) +
    scale_x_continuous(breaks = axis_breaks, labels = function(x) sprintf("%.2f", x)) +
    labs(title = paste(identity_name, axis, "distribution"), x = paste0(axis, " bin (mm)"), y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))

  p <- base_plot + geom_line(size = line_size)
  ggsave(filename = file.path(output_dir, paste0("lineplot_", identity_name, "_", axis, ".pdf")),
         plot = p, width = 6, height = 3)

  p_smooth <- base_plot +
    geom_smooth(method = "loess", span = loess_span, se = FALSE, size = smooth_line_size) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste(identity_name, axis, "smoothed distribution"), x = paste0(axis, " (mm)"))
  ggsave(filename = file.path(output_dir, paste0("smoothed_lineplot_", identity_name, "_", axis, ".pdf")),
         plot = p_smooth, width = 6, height = 3)
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory containing binned CSVs"),
    make_option("--output-dir", type = "character", help = "Directory for line plots"),
    make_option("--identity-groups", type = "character", default = "NN=0:7;GABA=8:19;Glut=20:29",
                help = "Cluster group definitions [default %default]"),
    make_option("--hex-map", type = "character", default = NULL,
                help = "CSV with cluster hex colours"),
    make_option("--line-size", type = "double", default = 0.8,
                help = "Line width for raw plots [default %default]"),
    make_option("--smooth-line-size", type = "double", default = 0.8,
                help = "Line width for LOESS plots [default %default]"),
    make_option("--loess-span", type = "double", default = 0.18,
                help = "Span parameter for LOESS [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step02_coordinate_line_plots_loess.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) {
    stop("--input-dir and --output-dir are required")
  }
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  hex_candidates <- c(
    args$`hex-map`,
    file.path(args$`input-dir`, "cluster-hex-colors.csv"),
    file.path(args$`input-dir`, "cluster_hex_colors.csv"),
    "input/live_data/metadata/cluster_hex_colors.csv",
    "input/live_data/metadata/cluster-hex-colors.csv"
  )
  hex_candidates <- hex_candidates[!is.null(hex_candidates) & nzchar(hex_candidates)]
  hex_path <- hex_candidates[file.exists(hex_candidates)][1]
  if (is.na(hex_path)) hex_path <- hex_candidates[1]
  if (!file.exists(hex_path)) stop("Hex colour CSV not found: ", hex_path)
  hex_df <- read_csv(hex_path, show_col_types = FALSE)
  names(hex_df) <- sub("^\ufeff", "", names(hex_df))
  if (!all(c("cluster", "custom_hex") %in% names(hex_df))) {
    stop("Hex colour CSV must contain 'cluster' and 'custom_hex' columns")
  }
  hex_colors <- setNames(hex_df$custom_hex, as.character(hex_df$cluster))

  identity_groups <- parse_groups(args$`identity-groups`)

  for (axis in c("x", "y", "z")) {
    combined_path <- file.path(args$`input-dir`, paste0("combined_binned_", axis, ".csv"))
    if (!file.exists(combined_path)) {
      warning("Missing file ", combined_path, "; skipping axis")
      next
    }
    combined <- read_csv(combined_path, show_col_types = FALSE)
    long_df <- combined %>%
      pivot_longer(-bin, names_to = "cluster", values_to = "count") %>%
      mutate(bin = factor(bin, levels = unique(bin)),
             bin_lower = as.numeric(sub("\\[(.+),.*", "\\1", bin)))

    for (group_name in names(identity_groups)) {
      clusters <- identity_groups[[group_name]]
      build_plots(long_df, clusters, hex_colors, axis, group_name, args$`output-dir`,
                  args$`line-size`, args$`smooth-line-size`, args$`loess-span`)
    }
  }
}

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

if (sys.nframe() == 0L) {
  main()
}
