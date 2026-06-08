#!/usr/bin/env Rscript

# Generate volcano plots for differential expression tables.
#
# Inputs
#   --deg-csv: differential expression CSV with columns `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2`, `pct_diff`, `cluster`, `gene` (required).
#   --output-dir: directory for PDF output (required).
#   --comparison-name: label used in plot filenames (default derived from CSV name).
#   --min-pct-diff: minimum |pct_diff| to classify as high-percentage (default: 0.1).
#   --min-p-adj: adjusted p-value threshold (default: 0.01).
#   --min-logfc: minimum absolute log2 fold-change (default: 0.5).
#   --clusters: optional comma-separated list of cluster IDs; otherwise iterate all.
#   --patterns-to-remove: optional comma-separated regex patterns to exclude genes.
#   --label-genes: optional comma-separated list of gene symbols to label (overrides auto-labelling).
#   --labels-per-direction: number of positive/negative genes to label when auto-selecting (default: 20).
#   --colors: comma-separated colors for ns, low_pct_pos, low_pct_neg, high_pct_pos, high_pct_neg (default: gray90,pink2,lightblue2,firebrick,dodgerblue3).
#   --plot-width/--plot-height: PDF dimensions (default: 4x4).
#
# Dependencies: optparse, ggplot2, dplyr, readr, ggrepel

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

select_labels <- function(df, min_p_adj, min_fc, min_pct_diff, label_genes, labels_per_direction) {
  if (length(label_genes)) {
    return(df %>% filter(gene %in% label_genes))
  }
  significant <- df %>%
    filter(p_val_adj <= min_p_adj,
           abs(avg_log2FC) >= min_fc,
           abs_pct_diff >= min_pct_diff)
  pos <- significant %>% filter(avg_log2FC > 0) %>% arrange(desc(vector_length)) %>% head(labels_per_direction)
  neg <- significant %>% filter(avg_log2FC < 0) %>% arrange(desc(vector_length)) %>% head(labels_per_direction)
  bind_rows(pos, neg)
}

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "Differential expression CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--comparison-name", type = "character", default = NULL),
    make_option("--min-pct-diff", type = "double", default = 0.1),
    make_option("--min-p-adj", type = "double", default = 0.01),
    make_option("--min-logfc", type = "double", default = 0.5),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--patterns-to-remove", type = "character", default = NULL),
    make_option("--label-genes", type = "character", default = NULL),
    make_option("--labels-per-direction", type = "integer", default = 20),
    make_option("--colors", type = "character", default = "gray90,pink2,lightblue2,firebrick,dodgerblue3"),
    make_option("--plot-width", type = "double", default = 4),
    make_option("--plot-height", type = "double", default = 4)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "deg_volcano_plot.R --deg-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-csv`, args$`output-dir`)))) {
    stop("--deg-csv and --output-dir are required")
  }
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  df <- read_csv(args$`deg-csv`, show_col_types = FALSE)
  required_columns <- c("avg_log2FC", "p_val_adj", "pct.1", "pct.2", "pct_diff", "cluster", "gene")
  missing <- setdiff(required_columns, names(df))
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

  clusters <- if (!is.null(args$clusters) && nzchar(args$clusters)) {
    trimws(strsplit(args$clusters, ",")[[1]])
  } else {
    sort(unique(df$cluster))
  }

  colors <- trimws(strsplit(args$colors, ",")[[1]])
  if (length(colors) != 5) stop("--colors must provide five comma-separated values")
  names(colors) <- c("ns", "low_pos", "low_neg", "high_pos", "high_neg")

  label_genes <- if (!is.null(args$`label-genes`) && nzchar(args$`label-genes`)) trimws(strsplit(args$`label-genes`, ",")[[1]]) else character()
  comparison <- args$`comparison-name` %||% tools::file_path_sans_ext(basename(args$`deg-csv`))
  pdf_path <- file.path(args$`output-dir`, sprintf("plot_volcano_%s_mpd-%.3f_mp-%.3f_mfc-%.2f.pdf",
                                                   comparison, args$`min-pct-diff`, args$`min-p-adj`, args$`min-logfc`))
  pdf(pdf_path, width = args$`plot-width`, height = args$`plot-height`)

  for (cluster_id in clusters) {
    cluster_data <- df %>% filter(cluster == cluster_id)
    if (nrow(cluster_data) == 0) next
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

    labeled <- select_labels(cluster_data, args$`min-p-adj`, args$`min-logfc`,
                             args$`min-pct-diff`, label_genes, args$`labels-per-direction`)
    x_lim <- c(min(cluster_data$avg_log2FC, na.rm = TRUE) * 1.1,
               max(cluster_data$avg_log2FC, na.rm = TRUE) * 1.1)
    y_lim <- c(0, max(cluster_data$neg_log_padj, na.rm = TRUE) * 1.1)

    p <- ggplot(cluster_data, aes(x = avg_log2FC, y = neg_log_padj)) +
      geom_point(aes(color = color), size = 1, alpha = 0.8) +
      scale_color_identity() +
      geom_vline(xintercept = c(-args$`min-logfc`, args$`min-logfc`), colour = "gray20", linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = -log10(args$`min-p-adj`), colour = "gray20", linetype = "dashed", linewidth = 0.3) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      theme_classic(base_size = 10) +
      labs(title = paste("Cluster", cluster_id), x = "log2 Fold Change", y = "-log10 adjusted p-value")

    if (nrow(labeled)) {
      p <- p + geom_text_repel(data = labeled,
                               aes(label = gene), size = 2.5,
                               colour = "gray20", segment.colour = "gray50", segment.size = 0.2,
                               max.overlaps = Inf)
    }
    print(p)
  }

  dev.off()
  message("Saved volcano plots to ", pdf_path)
}

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

if (sys.nframe() == 0L) {
  main()
}
