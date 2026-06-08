#!/usr/bin/env Rscript

# Generate dot plots from a FindAllMarkers CSV using manually supplied gene lists.
#
# Inputs
#   --deg-csv: FindAllMarkers CSV with columns `gene`, `cluster`, `pct.1`, `avg_log2FC`, `p_val_adj` (required).
#   --genes: comma-separated list of genes (required).
#   --output-pdf: destination PDF (required).
#   --max-dot-size: maximum dot size (default: 2.5).
#   --width/--height: PDF size (default: 10 x 5).
#   --p-threshold: adjusted p-value threshold for outline (default: 0.05).
#   --logfc-cap: cap for log2 fold-change colour scale (default: 3.5).
#   --palette-low/mid/high: colours for gradient (defaults: dodgerblue, white, firebrick).
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
    make_option("--deg-csv", type = "character", help = "FindAllMarkers CSV"),
    make_option("--genes", type = "character", help = "Comma-separated gene list"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--max-dot-size", type = "double", default = 2.5),
    make_option("--width", type = "double", default = 10),
    make_option("--height", type = "double", default = 5),
    make_option("--p-threshold", type = "double", default = 0.05),
    make_option("--logfc-cap", type = "double", default = 3.5),
    make_option("--palette-low", type = "character", default = "dodgerblue"),
    make_option("--palette-mid", type = "character", default = "white"),
    make_option("--palette-high", type = "character", default = "firebrick")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_findallmarkers_dot_manual.R --deg-csv FILE --genes G1,G2 --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$genes, args$`output-pdf`)))) {
    stop("--deg-csv, --genes, and --output-pdf are required")
  }
  if (!file.exists(args$`deg-csv`)) stop("FindAllMarkers CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  genes <- trimws(strsplit(args$genes, ",")[[1]])
  data <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  stopifnot(all(c("gene", "cluster", "pct.1", "avg_log2FC", "p_val_adj") %in% names(data)))

  data_filtered <- data %>%
    filter(gene %in% genes) %>%
    mutate(cluster = as.numeric(as.character(cluster)),
           dot_outline = ifelse(p_val_adj < args$`p-threshold` & avg_log2FC > 0, "black", "white")) %>%
    mutate(avg_log2FC = pmax(pmin(avg_log2FC, args$`logfc-cap`), -args$`logfc-cap`)) %>%
    arrange(factor(gene, levels = genes), cluster)

  p <- ggplot(data_filtered, aes(x = factor(cluster), y = factor(gene, levels = genes))) +
    geom_hline(aes(yintercept = as.numeric(factor(gene, levels = genes))), colour = "gray70", linewidth = 0.3) +
    geom_point(aes(size = pct.1, fill = avg_log2FC, colour = dot_outline), shape = 21) +
    scale_size_continuous(range = c(1, args$`max-dot-size`), limits = c(0, 1)) +
    scale_fill_gradient2(low = args$`palette-low`, mid = args$`palette-mid`, high = args$`palette-high`,
                         midpoint = 0, limits = c(-args$`logfc-cap`, args$`logfc-cap`)) +
    scale_color_identity() +
    labs(x = "Cluster", y = "Gene", size = "Pct.1", fill = "Avg Log2FC") +
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
