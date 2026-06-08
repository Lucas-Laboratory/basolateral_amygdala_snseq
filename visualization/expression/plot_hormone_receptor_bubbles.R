#!/usr/bin/env Rscript

# Generate hormone receptor bubble plots from a FindAllMarkers table.
#
# Inputs
#   --deg-csv: differential expression CSV with columns `gene`, `cluster`, `pct.1`, `pct.2`, `avg_log2FC`, `p_val_adj` (required).
#   --genes: comma-separated list of genes to plot (required).
#   --output-dir: directory for output PDFs (required).
#   --width/--height: PDF size (default: 7 x 2).
#   --max-dot-size: maximum bubble size (default: 5).
#   --logfc-cap: cap for log2FC colour scale (default: 3.5).
#   --p-threshold: adjusted p-value threshold for outlines (default: 0.05).
#   --colors: comma-separated gradient colours low,mid,high (default: dodgerblue,white,firebrick).
#
# Dependencies: optparse, ggplot2, dplyr, readr, scales

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

geom_mean <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(0)
  exp(mean(log(x)))
}

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "DEG CSV"),
    make_option("--genes", type = "character", help = "Comma-separated genes"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--width", type = "double", default = 7),
    make_option("--height", type = "double", default = 2),
    make_option("--max-dot-size", type = "double", default = 5),
    make_option("--logfc-cap", type = "double", default = 3.5),
    make_option("--p-threshold", type = "double", default = 0.05),
    make_option("--colors", type = "character", default = "dodgerblue,white,firebrick")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_hormone_receptor_bubbles.R --deg-csv FILE --genes G1,G2 --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$genes, args$`output-dir`)))) stop("All inputs required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  genes <- trimws(strsplit(args$genes, ",")[[1]])
  data <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  required_cols <- c("gene", "cluster", "pct.1", "pct.2", "avg_log2FC", "p_val_adj")
  if (length(setdiff(required_cols, names(data)))) stop("CSV missing required columns")

  colors <- trimws(strsplit(args$colors, ",")[[1]])
  if (length(colors) != 3) stop("--colors must provide three values")

  data <- data %>% mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-324, p_val_adj),
                          pct_diff = abs(pct.1 - pct.2))

  for (current_gene in genes) {
    df_gene <- data %>% filter(.data$gene == current_gene)
    if (!nrow(df_gene)) {
      message("Skipping ", current_gene, ": no rows found")
      next
    }
    df_gene <- df_gene %>%
      arrange(desc(pct.1), cluster) %>%
      mutate(cluster = factor(as.character(cluster), levels = unique(as.character(cluster))),
             bubble_size = scales::rescale(-log10(p_val_adj), to = c(1, args$`max-dot-size`)),
             outline_color = ifelse(p_val_adj < args$`p-threshold`, "black", "white"),
             avg_log2FC = pmax(pmin(avg_log2FC, args$`logfc-cap`), -args$`logfc-cap`))

    y_max <- max(df_gene$pct.1, na.rm = TRUE) * 1.1
    p <- ggplot(df_gene, aes(x = cluster, y = pct.1, size = bubble_size)) +
      geom_point(aes(fill = avg_log2FC, colour = outline_color), shape = 21, stroke = 0.5) +
      scale_fill_gradientn(colours = colors, limits = c(-args$`logfc-cap`, args$`logfc-cap`), oob = scales::squish) +
      scale_color_identity() +
      scale_size(range = c(1, args$`max-dot-size`)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "gray80") +
      coord_cartesian(ylim = c(-0.1, y_max)) +
      theme_minimal(base_size = 10) +
      labs(title = paste("Bubble Plot for", current_gene), x = "Cluster", y = "pct.1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(args$`output-dir`, paste0("BubblePlot_", current_gene, ".pdf")),
           p, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
