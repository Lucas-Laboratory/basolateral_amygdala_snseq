#!/usr/bin/env Rscript

# Summarise barcode counts per cluster across a predefined set of sample prefixes.
#
# Inputs
#   --input-csv: CSV containing at minimum `Barcode` and `Cluster` columns (required).
#   --output-csv: destination CSV for the summary table (required).
#   --prefixes: comma-separated list of barcode prefixes (default: inferred from all barcodes).
#   --delimiter: delimiter separating prefix from barcode suffix (default: `_`).
#
# Output
#   CSV reporting total counts per cluster plus columns `<prefix>_barcodes_count`.
#
# Dependencies: optparse, readr, dplyr, tidyr, stringr

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

build_summary <- function(data, prefixes, delimiter) {
  if (!all(c("Barcode", "Cluster") %in% names(data))) {
    stop("Input CSV must include 'Barcode' and 'Cluster' columns")
  }

  data <- data %>%
    mutate(
      Cluster = as.character(.data$Cluster),
      prefix = str_split_fixed(.data$Barcode, delimiter, n = 2)[, 1]
    )

  if (is.null(prefixes)) {
    prefixes <- sort(unique(data$prefix))
  }

  counts <- data %>%
    count(Cluster, prefix) %>%
    complete(Cluster, prefix = prefixes, fill = list(n = 0))

  totals <- data %>%
    count(Cluster, name = "total_barcodes_count")

  wide_counts <- counts %>%
    mutate(prefix = paste0(prefix, "_barcodes_count")) %>%
    pivot_wider(names_from = prefix, values_from = n)

  summary <- totals %>%
    full_join(wide_counts, by = "Cluster") %>%
    arrange(Cluster)

  totals_all <- data %>%
    summarise(
      Cluster = "all",
      total_barcodes_count = n()
    )

  per_prefix_all <- data %>%
    count(prefix) %>%
    mutate(prefix = paste0(prefix, "_barcodes_count")) %>%
    pivot_wider(names_from = prefix, values_from = n)

  bind_rows(
    totals_all %>% bind_cols(per_prefix_all),
    summary
  )
}

main <- function() {
  option_list <- list(
    make_option("--input-csv", type = "character", help = "Barcode metadata CSV"),
    make_option("--output-csv", type = "character", help = "Destination summary CSV"),
    make_option("--prefixes", type = "character", default = NULL,
                help = "Comma-separated barcode prefixes to include"),
    make_option("--delimiter", type = "character", default = "_",
                help = "Delimiter separating prefix from barcode [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "summarize_barcodes_by_cluster_and_sample.R --input-csv FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  input_csv <- args$`input-csv`
  output_csv <- args$`output-csv`
  prefixes <- args$prefixes

  if (is.null(input_csv) || !nzchar(input_csv)) stop("--input-csv is required")
  if (is.null(output_csv) || !nzchar(output_csv)) stop("--output-csv is required")
  if (!file.exists(input_csv)) stop("Input CSV not found: ", input_csv)

  data <- readr::read_csv(input_csv, show_col_types = FALSE)
  prefix_vec <- if (is.null(prefixes) || !nzchar(prefixes)) NULL else trimws(strsplit(prefixes, ",")[[1]])

  summary <- build_summary(data, prefix_vec, args$delimiter)
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary, output_csv)
  message("Summary written to ", output_csv)
}

if (sys.nframe() == 0L) {
  main()
}
