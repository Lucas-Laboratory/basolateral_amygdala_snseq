#!/usr/bin/env Rscript

# Generate volcano plots for all clusters in a FindAllMarkers CSV.
#
# Inputs
#   --deg-csv: FindAllMarkers CSV with columns `avg_log2FC`, `p_val_adj`, `pct.1`, `pct_diff`, `cluster`, `gene` (required).
#   --output-dir: output directory (required).
#   --min-pct-diff: minimum |pct_diff| for high-significance colouring (default: 0.5).
#   --min-p-adj: adjusted p-value threshold (default: 0.01).
#   --min-logfc: minimum |log2FC| threshold (default: 2.5).
#   --clusters: optional comma-separated list of clusters; otherwise iterate all.
#   --patterns-to-remove: optional comma-separated regex patterns for gene removal.
#   --label-genes: optional comma-separated list of gene symbols to label.
#   --labels-per-direction: number of positive/negative labels when auto-selected (default: 20).
#   --colors: comma-separated values for ns, low_pos, low_neg, high_pos, high_neg (default: gray90,pink2,lightblue2,firebrick,dodgerblue3).
#   --plot-width/--plot-height: PDF size (default: 4x4).
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
    make_option("--deg-csv", type = "character", help = "FindAllMarkers CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--min-pct-diff", type = "double", default = 0.5),
    make_option("--min-p-adj", type = "double", default = 0.01),
    make_option("--min-logfc", type = "double", default = 2.5),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--patterns-to-remove", type = "character", default = NULL),
    make_option("--label-genes", type = "character", default = NULL),
    make_option("--labels-per-direction", type = "integer", default = 20),
    make_option("--colors", type = "character", default = "gray90,pink2,lightblue2,firebrick,dodgerblue3"),
    make_option("--plot-width", type = "double", default = 4),
    make_option("--plot-height", type = "double", default = 4)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_findallmarkers_volcano.R --deg-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$`output-dir`)))) stop("--deg-csv and --output-dir are required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  df <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  names(df) <- sub("^\ufeff", "", names(df))
  if (!"pct_diff" %in% names(df) && all(c("pct.1", "pct.2") %in% names(df))) {
    df <- df %>% mutate(pct_diff = pct.1 - pct.2)
  }
  required_cols <- c("avg_log2FC", "p_val_adj", "pct.1", "pct_diff", "cluster", "gene")
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

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
  if (length(colors) != 5) stop("--colors must provide five values")
  names(colors) <- c("ns", "low_pos", "low_neg", "high_pos", "high_neg")
  label_genes <- if (!is.null(args$`label-genes`) && nzchar(args$`label-genes`)) trimws(strsplit(args$`label-genes`, ",")[[1]]) else character()

  pdf_path <- file.path(args$`output-dir`, sprintf("plot_volcano_FindAll_mpd-%.3f_mp-%.3f_mfc-%.2f.pdf",
                                                   args$`min-pct-diff`, args$`min-p-adj`, args$`min-logfc`))
  pdf(pdf_path, width = args$`plot-width`, height = args$`plot-height`)

  for (cl in clusters) {
    cluster_data <- df %>% filter(cluster == cl)
    if (!nrow(cluster_data)) next
    cluster_data <- cluster_data %>%
      mutate(color = case_when(
        (p_val_adj > args$`min-p-adj` | abs(avg_log2FC) < args$`min-logfc`) ~ colors["ns"],
        (abs_pct_diff < args$`min-pct-diff` & avg_log2FC > 0) ~ colors["low_pos"],
        (abs_pct_diff < args$`min-pct-diff` & avg_log2FC < 0) ~ colors["low_neg"],
        (abs_pct_diff >= args$`min-pct-diff` & avg_log2FC > 0) ~ colors["high_pos"],
        (abs_pct_diff >= args$`min-pct-diff` & avg_log2FC < 0) ~ colors["high_neg"]
      )) %>%
      mutate(order = case_when(
        color == colors["ns"] ~ 1,
        color %in% colors[c("low_pos", "low_neg")] ~ 2,
        TRUE ~ 3
      )) %>% arrange(order)

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
  message("Saved FindAll volcano plots to ", pdf_path)
}

if (sys.nframe() == 0L) {
  main()
}
