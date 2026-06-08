#!/usr/bin/env Rscript

# Create a DEG dot plot for a manually specified gene list.
#
# Inputs
#   --deg-csv: differential expression CSV with columns `gene`, `cluster`, `pct.1`, `pct.2`, `avg_log2FC`, `p_val_adj` (required).
#   --genes: comma-separated gene list (required).
#   --output-pdf: destination PDF (required).
#   --max-dot-size: maximum dot size (default: 2.5).
#   --width/--height: PDF dimensions (default: 10 x 3).
#   --p-threshold: adjusted p-value for outlining (default: 0.001).
#   --logfc-cap: cap for log2 fold-change colour scale (default: 3.5).
#   --size-mode: `pct_diff` (abs(pct.1 - pct.2)) or `pct1` (default: `pct_diff`).
#
# Dependencies: optparse, ggplot2, dplyr, readr

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "DEG CSV"),
    make_option("--genes", type = "character", help = "Comma-separated genes"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--max-dot-size", type = "double", default = 2.5),
    make_option("--width", type = "double", default = 10),
    make_option("--height", type = "double", default = 3),
    make_option("--p-threshold", type = "double", default = 0.001),
    make_option("--logfc-cap", type = "double", default = 3.5),
    make_option("--size-mode", type = "character", default = "pct_diff")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_deg_dot_manual_selection.R --deg-csv FILE --genes G1,G2 --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$genes, args$`output-pdf`)))) stop("All inputs required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  genes <- trimws(strsplit(args$genes, ",")[[1]])
  data <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  size_col <- if (args$`size-mode` == "pct1") "pct.1" else "pct_diff"
  if (!all(c("gene", "cluster", "pct.1", "pct.2", "avg_log2FC", "p_val_adj") %in% names(data))) {
    stop("Input CSV missing required columns")
  }
  data <- data %>% mutate(pct_diff = abs(pct.1 - pct.2))

  data_filtered <- data %>%
    filter(gene %in% genes) %>%
    mutate(cluster = as.numeric(as.character(cluster)),
           dot_outline = ifelse(p_val_adj < args$`p-threshold`, "black", "white"),
           avg_log2FC = pmax(pmin(avg_log2FC, args$`logfc-cap`), -args$`logfc-cap`)) %>%
    arrange(factor(gene, levels = genes), cluster)

  p <- ggplot(data_filtered, aes(x = factor(cluster), y = factor(gene, levels = genes))) +
    geom_hline(aes(yintercept = as.numeric(factor(gene, levels = genes))),
               colour = "gray70", linewidth = 0.3) +
    geom_point(aes(size = .data[[size_col]], fill = avg_log2FC, colour = dot_outline), shape = 21) +
    scale_size_continuous(range = c(1, args$`max-dot-size`), limits = c(0, 1)) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-args$`logfc-cap`, args$`logfc-cap`)) +
    scale_color_identity() +
    labs(x = "Cluster", y = "Gene", size = size_col, fill = "Avg Log2FC") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          legend.position = "right", legend.box = "vertical",
          panel.grid = element_blank())

  ggsave(args$`output-pdf`, p, width = args$width, height = args$height)
  message("Dot plot saved to ", args$`output-pdf`)
}

if (sys.nframe() == 0L) {
  main()
}
