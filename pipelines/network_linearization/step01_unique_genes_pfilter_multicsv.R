#!/usr/bin/env Rscript

# Collect unique genes across DEG tables using significance and fold-change thresholds.
#
# Inputs
#   --input-dir: directory containing DEG CSV files (required).
#   --output-csv: destination CSV listing unique genes (required).
#   --pattern: regex to match filenames (default: `\\.csv$`).
#   --pval-cutoff: adjusted p-value threshold (default: 0.01).
#   --logfc-cutoff: absolute log2 fold-change threshold (default: 1).
#   --gene-column: column holding gene symbols (default: `gene`).
#   --pval-column: column with adjusted p-values (default: `p_val_adj`).
#   --logfc-column: column with log2 fold-change values (default: `avg_log2FC`).
#
# Output
#   CSV with one column `gene` listing unique genes that meet criteria in any file.
#
# Dependencies: optparse, readr, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory containing DEG CSVs"),
    make_option("--output-csv", type = "character", help = "Destination CSV"),
    make_option("--pattern", type = "character", default = "\\.csv$",
                help = "Regular expression to select files [default %default]"),
    make_option("--pval-cutoff", type = "double", default = 0.01,
                help = "Adjusted p-value cutoff [default %default]"),
    make_option("--logfc-cutoff", type = "double", default = 1,
                help = "Absolute log2 fold-change cutoff [default %default]"),
    make_option("--gene-column", type = "character", default = "gene"),
    make_option("--pval-column", type = "character", default = "p_val_adj"),
    make_option("--logfc-column", type = "character", default = "avg_log2FC")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step01_unique_genes_pfilter_multicsv.R --input-dir DIR --output-csv FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-csv`)))) {
    stop("--input-dir and --output-csv are required")
  }
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)

  files <- list.files(args$`input-dir`, pattern = args$pattern, full.names = TRUE)
  if (!length(files)) stop("No files matched pattern '", args$pattern, "' in ", args$`input-dir`)

  genes <- unique(unlist(lapply(files, function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE)
    required <- c(args$`gene-column`, args$`pval-column`, args$`logfc-column`)
    if (!all(required %in% names(df))) {
      stop("File ", basename(path), " missing columns: ", paste(setdiff(required, names(df)), collapse = ", "))
    }
    df %>%
      filter(!is.na(.data[[args$`pval-column`]]), !is.na(.data[[args$`logfc-column`]])) %>%
      filter(.data[[args$`pval-column`]] <= args$`pval-cutoff`,
             abs(.data[[args$`logfc-column`]]) >= args$`logfc-cutoff`) %>%
      pull(.data[[args$`gene-column`]])
  })))

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble(gene = sort(genes)), args$`output-csv`)
  message("Identified ", length(genes), " unique genes")
}

if (sys.nframe() == 0L) {
  main()
}
