#!/usr/bin/env Rscript

# Generate volcano plots for selected clusters from a DEG table using manual label overrides.
#
# Inputs
#   --deg-csv: DEG CSV (cluster-wise) with columns `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2`, `pct_diff`, `cluster`, `gene` (required).
#   --output-pdf: destination PDF (required).
#   --clusters: comma-separated cluster IDs to plot (default: all).
#   --min-pct-diff: minimum |pct_diff| (default: 0.1).
#   --min-p-adj: adjusted p-value cutoff (default: 0.01).
#   --min-logfc: minimum |log2FC| (default: 1).
#   --label-genes: optional comma-separated list of genes to label.
#   --labels-per-direction: number of auto-selected labels per direction (default: 20).
#   --patterns-to-remove: optional regex patterns to drop genes.
#   --colors: comma-separated ns, low_pos, low_neg, high_pos, high_neg (default: gray90,pink2,lightblue2,firebrick,dodgerblue3).
#   --width/--height: PDF size (default: 4 x 4).
#
# Dependencies: optparse, ggplot2, dplyr, readr, ggrepel

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

select_labels <- function(df, min_p_adj, min_logfc, min_pct_diff, label_genes, labels_per_direction) {
  if (length(label_genes)) return(df %>% filter(gene %in% label_genes))
  significant <- df %>% filter(p_val_adj <= min_p_adj, abs(avg_log2FC) >= min_logfc, abs_pct_diff >= min_pct_diff)
  pos <- significant %>% filter(avg_log2FC > 0) %>% arrange(desc(vector_length)) %>% head(labels_per_direction)
  neg <- significant %>% filter(avg_log2FC < 0) %>% arrange(desc(vector_length)) %>% head(labels_per_direction)
  bind_rows(pos, neg)
}

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "DEG CSV"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--min-pct-diff", type = "double", default = 0.1),
    make_option("--min-p-adj", type = "double", default = 0.01),
    make_option("--min-logfc", type = "double", default = 1),
    make_option("--label-genes", type = "character", default = NULL),
    make_option("--labels-per-direction", type = "integer", default = 20),
    make_option("--patterns-to-remove", type = "character", default = NULL),
    make_option("--colors", type = "character", default = "gray90,pink2,lightblue2,firebrick,dodgerblue3"),
    make_option("--width", type = "double", default = 4),
    make_option("--height", type = "double", default = 4)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_deg_volcano_manual_selection.R --deg-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$`output-pdf`)))) stop("--deg-csv and --output-pdf are required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  df <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  required_cols <- c("avg_log2FC", "p_val_adj", "pct.1", "pct.2", "pct_diff", "cluster", "gene")
  if (length(setdiff(required_cols, names(df)))) stop("CSV missing required columns")

  patterns_default <- c("^mt-", "^Rpl", "^Mrpl", "^Rps", "^Nduf", "\\d+Rik$", "^Gm\\d+",
                        "os$", "os\\d+$", "^[0-9]+$", "^[0-9]{4,}$", "^[A-Z][0-9]{4,}$",
                        "^[A-Z]{2}[0-9]{4,}(\\\\.[0-9]+)?$", "^D[0-9]{1,2}E.*$", "^[0-9]{5}$")
  if (!is.null(args$`patterns-to-remove`) && nzchar(args$`patterns-to-remove`)) {
    patterns_default <- unique(c(patterns_default, trimws(strsplit(args$`patterns-to-remove`, ",")[[1]])))
  }
  df <- df %>% filter(!grepl(paste(patterns_default, collapse = "|"), gene)) %>% filter(pct.1 > 0)
  df <- df %>% mutate(p_val_adj = ifelse(p_val_adj == 0, 1.5e-303, p_val_adj),
                      neg_log_padj = -log10(p_val_adj),
                      abs_pct_diff = abs(pct_diff),
                      vector_length = sqrt(avg_log2FC^2 + neg_log_padj^2))

  clusters <- if (!is.null(args$clusters) && nzchar(args$clusters)) trimws(strsplit(args$clusters, ",")[[1]]) else sort(unique(df$cluster))
  colors <- trimws(strsplit(args$colors, ",")[[1]])
  if (length(colors) != 5) stop("--colors must supply five values")
  names(colors) <- c("ns", "low_pos", "low_neg", "high_pos", "high_neg")
  label_genes <- if (!is.null(args$`label-genes`) && nzchar(args$`label-genes`)) trimws(strsplit(args$`label-genes`, ",")[[1]]) else character()

  pdf(args$`output-pdf`, width = args$width, height = args$height)
  for (cl in clusters) {
    cluster_data <- df %>% filter(cluster == cl)
    if (!nrow(cluster_data)) next
    cluster_data <- cluster_data %>% mutate(color = case_when(
      (p_val_adj > args$`min-p-adj` | abs(avg_log2FC) < args$`min-logfc`) ~ colors["ns"],
      (abs_pct_diff < args$`min-pct-diff` & avg_log2FC > 0) ~ colors["low_pos"],
      (abs_pct_diff < args$`min-pct-diff` & avg_log2FC < 0) ~ colors["low_neg"],
      (abs_pct_diff >= args$`min-pct-diff` & avg_log2FC > 0) ~ colors["high_pos"],
      (abs_pct_diff >= args$`min-pct-diff` & avg_log2FC < 0) ~ colors["high_neg"]
    )) %>% mutate(order = case_when(color == colors["ns"] ~ 1,
                                    color %in% colors[c("low_pos", "low_neg")] ~ 2,
                                    TRUE ~ 3)) %>% arrange(order)

    labeled <- select_labels(cluster_data, args$`min-p-adj`, args$`min-logfc`, args$`min-pct-diff`, label_genes, args$`labels-per-direction`)
    x_lim <- c(min(cluster_data$avg_log2FC, na.rm = TRUE) * 1.1, max(cluster_data$avg_log2FC, na.rm = TRUE) * 1.1)
    y_lim <- c(0, max(cluster_data$neg_log_padj, na.rm = TRUE) * 1.1)

    p <- ggplot(cluster_data, aes(x = avg_log2FC, y = neg_log_padj, colour = color)) +
      geom_point(size = 1) +
      geom_hline(yintercept = -log10(args$`min-p-adj`), linetype = "dashed", colour = "gray20", linewidth = 0.3) +
      geom_vline(xintercept = c(-args$`min-logfc`, args$`min-logfc`), linetype = "dashed", colour = "gray20", linewidth = 0.3) +
      scale_color_identity() +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      theme_minimal(base_size = 10) +
      labs(title = paste("Cluster", cl), x = "log2 Fold Change", y = "-log10 adjusted p-value")

    if (nrow(labeled)) {
      p <- p + geom_label_repel(data = labeled, aes(label = gene), size = 2.5, colour = "gray20",
                                fill = "white", segment.colour = "gray50", segment.size = 0.2,
                                max.overlaps = 50, box.padding = 0.15)
    }
    print(p)
  }
  dev.off()
  message("Saved manual DEG volcano plots to ", args$`output-pdf`)
}

if (sys.nframe() == 0L) {
  main()
}
