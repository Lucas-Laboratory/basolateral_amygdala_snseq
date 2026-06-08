#!/usr/bin/env Rscript

# Generate cluster-by-gene heatmaps ordered by a STRING-derived ranking.
#
# Inputs
#   --deg-dir: directory containing DEG CSV files with columns `cluster`, `gene`, `avg_log2FC` (required).
#   --string-csv: CSV from previous step with ordered gene symbols (required).
#   --output-dir: directory for heatmap PDFs (required).
#   --palette-length: number of colours in the diverging palette (default: 100).
#   --value-range: absolute value range for heatmap colour breaks (default: 2.5).
#
# Output
#   One PDF per DEG input, sorted by the gene order supplied in --string-csv.
#
# Dependencies: optparse, readr, dplyr, tidyr, pheatmap

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})

build_heatmap <- function(deg_path, ordered_genes, output_path, palette_length, value_range) {
  deg <- read_csv(deg_path, show_col_types = FALSE)
  required_cols <- c("cluster", "gene", "avg_log2FC")
  if (!all(required_cols %in% names(deg))) {
    stop("File ", basename(deg_path), " missing columns: ", paste(setdiff(required_cols, names(deg)), collapse = ", "))
  }

  matrix_data <- deg %>%
    filter(gene %in% ordered_genes) %>%
    group_by(gene, cluster) %>%
    summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop") %>%
    mutate(gene = factor(gene, levels = ordered_genes)) %>%
    pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0) %>%
    arrange(gene) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  if (!nrow(matrix_data)) {
    warning("No overlapping genes for ", basename(deg_path), "; skipping heatmap")
    return(invisible(NULL))
  }

  col_order <- sort(as.numeric(colnames(matrix_data)))
  if (!any(is.na(col_order))) {
    matrix_data <- matrix_data[, as.character(col_order), drop = FALSE]
  }

  palette <- colorRampPalette(c("dodgerblue", "white", "firebrick"))(palette_length)
  breaks <- seq(-value_range, value_range, length.out = palette_length + 1)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  pdf(output_path, width = 5, height = 12)
  on.exit(dev.off(), add = TRUE)
  pheatmap(matrix_data,
           color = palette,
           breaks = breaks,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize = 4)
}

main <- function() {
  option_list <- list(
    make_option("--deg-dir", type = "character", help = "Directory of cluster DEG CSVs"),
    make_option("--string-csv", type = "character", help = "Ordered gene list CSV"),
    make_option("--output-dir", type = "character", help = "Directory for heatmap PDFs"),
    make_option("--palette-length", type = "integer", default = 100,
                help = "Number of colours in the palette [default %default]"),
    make_option("--value-range", type = "double", default = 2.5,
                help = "Absolute value range for colour scale [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step06_create_ranked_heatmap.R --deg-dir DIR --string-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`deg-dir`, args$`string-csv`, args$`output-dir`)))) {
    stop("--deg-dir, --string-csv, and --output-dir are required")
  }
  if (!dir.exists(args$`deg-dir`)) stop("DEG directory not found: ", args$`deg-dir`)
  if (!file.exists(args$`string-csv`)) stop("STRING order CSV not found: ", args$`string-csv`)

  ordered_genes <- read_csv(args$`string-csv`, show_col_types = FALSE)$gene
  ordered_genes <- unique(ordered_genes)
  if (!length(ordered_genes)) stop("Ordered gene list is empty")

  deg_files <- list.files(args$`deg-dir`, pattern = "\\.csv$", full.names = TRUE)
  if (!length(deg_files)) stop("No DEG CSV files found in ", args$`deg-dir`)

  for (file_path in deg_files) {
    output_pdf <- file.path(args$`output-dir`, paste0(tools::file_path_sans_ext(basename(file_path)), "_heatmap.pdf"))
    tryCatch(
      build_heatmap(file_path, ordered_genes, output_pdf, args$`palette-length`, args$`value-range`),
      error = function(err) warning("Failed to build heatmap for ", basename(file_path), ": ", err$message)
    )
  }
}

if (sys.nframe() == 0L) {
  main()
}
