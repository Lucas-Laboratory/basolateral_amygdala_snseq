#!/usr/bin/env Rscript

# Summarise barcode counts per cluster (and optionally per sample) from a tidy CSV table.
#
# Inputs
#   --input-csv: CSV containing barcode-level metadata (required).
#   --barcode-column: column holding cell barcodes (default: `Barcode`).
#   --cluster-column: column holding cluster identifiers (default: `Cluster`).
#   --sample-delimiter: character used to split sample identifiers from barcodes (default: `_`).
#   --output-csv: destination for the summary table (default: `<input_basename>_cluster_counts.csv`).
#
# Output
#   CSV with columns `cluster`, `sample`, `barcode_count`.
#
# Dependencies: optparse, readr, dplyr, stringr, rlang

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(rlang)
})

summarise_barcodes <- function(input_csv, barcode_col, cluster_col, delimiter, output_csv) {
  if (!file.exists(input_csv)) {
    stop("Input CSV not found: ", input_csv)
  }

  message("Reading ", input_csv)
  data <- readr::read_csv(input_csv, show_col_types = FALSE)

  if (!all(c(barcode_col, cluster_col) %in% names(data))) {
    missing <- setdiff(c(barcode_col, cluster_col), names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  barcode_sym <- sym(barcode_col)
  cluster_sym <- sym(cluster_col)

  data <- data %>%
    mutate(
      .sample = str_split_fixed(.data[[barcode_col]], delimiter, n = 2)[, 1]
    )

  per_sample <- data %>%
    count(!!cluster_sym, .sample, name = "barcode_count") %>%
    rename(cluster = !!cluster_sym, sample = .sample)

  totals <- data %>%
    count(!!cluster_sym, name = "barcode_count") %>%
    mutate(sample = "Total") %>%
    rename(cluster = !!cluster_sym)

  summary_df <- bind_rows(per_sample, totals) %>%
    arrange(cluster, sample)

  message("Writing summary to ", output_csv)
  readr::write_csv(summary_df, output_csv)
  invisible(summary_df)
}

main <- function() {
  option_list <- list(
    make_option("--input-csv", type = "character", help = "Path to barcode metadata CSV"),
    make_option("--barcode-column", type = "character", default = "Barcode",
                help = "Column containing barcode identifiers [default %default]"),
    make_option("--cluster-column", type = "character", default = "Cluster",
                help = "Column containing cluster identifiers [default %default]"),
    make_option("--sample-delimiter", type = "character", default = "_",
                help = "Delimiter splitting sample from barcode [default %default]"),
    make_option("--output-csv", type = "character", default = NULL,
                help = "Path to write summary CSV (default derived from input filename)")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "count_barcodes_per_cluster.R --input-csv FILE [options]")
  args <- parse_args(parser)

  input_csv <- args$`input-csv`
  barcode_col <- args$`barcode-column`
  cluster_col <- args$`cluster-column`
  delimiter <- args$`sample-delimiter`
  output_csv <- args$`output-csv`

  if (is.null(input_csv) || !nzchar(input_csv)) stop("--input-csv is required")
  if (is.null(output_csv) || !nzchar(output_csv)) {
    output_csv <- file.path(dirname(input_csv), paste0(tools::file_path_sans_ext(basename(input_csv)), "_cluster_counts.csv"))
  }

  summarise_barcodes(input_csv, barcode_col, cluster_col, delimiter, output_csv)
}

if (sys.nframe() == 0L) {
  main()
}
