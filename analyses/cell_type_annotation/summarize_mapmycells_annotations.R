#!/usr/bin/env Rscript

# Summarise MapMyCells hierarchical annotations per UMAP cluster.
#
# Inputs
#   --mapmy-csv: CSV with columns `barcode`, `class_name`, `subclass_name`, `supertype_name`, `cluster_name` (required).
#   --umap-csv: CSV linking barcodes to cluster IDs (required).
#   --output-csv: destination summary CSV (required).
#   --barcode-column: column in UMAP CSV containing barcodes (default: `Barcode`).
#   --cluster-column: column in UMAP CSV containing cluster IDs (default: `Cluster`).
#
# Output
#   Long-format table with counts and proportions per cluster and annotation level.
#
# Dependencies: optparse, readr, dplyr, tidyr

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
})

main <- function() {
  option_list <- list(
    make_option("--mapmy-csv", type = "character", help = "MapMyCells annotation CSV"),
    make_option("--umap-csv", type = "character", help = "Barcode-to-cluster CSV"),
    make_option("--output-csv", type = "character", help = "Destination summary CSV"),
    make_option("--barcode-column", type = "character", default = "Barcode",
                help = "Barcode column in UMAP CSV [default %default]"),
    make_option("--cluster-column", type = "character", default = "Cluster",
                help = "Cluster column in UMAP CSV [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "summarize_mapmycells_annotations.R --mapmy-csv FILE --umap-csv FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`mapmy-csv`, args$`umap-csv`, args$`output-csv`)))) {
    stop("--mapmy-csv, --umap-csv, and --output-csv are required")
  }
  if (!file.exists(args$`mapmy-csv`)) stop("MapMyCells CSV not found: ", args$`mapmy-csv`)
  if (!file.exists(args$`umap-csv`)) stop("UMAP CSV not found: ", args$`umap-csv`)

  map_df <- readr::read_csv(args$`mapmy-csv`, show_col_types = FALSE)
  umap_df <- readr::read_csv(args$`umap-csv`, show_col_types = FALSE)
  names(map_df) <- sub("^\ufeff", "", names(map_df))
  names(umap_df) <- sub("^\ufeff", "", names(umap_df))

  required_levels <- c("class_name", "subclass_name", "supertype_name", "cluster_name")
  if (!all(required_levels %in% names(map_df))) {
    stop("MapMyCells CSV missing columns: ", paste(setdiff(required_levels, names(map_df)), collapse = ", "))
  }

  if (!all(c(args$`barcode-column`, args$`cluster-column`) %in% names(umap_df))) {
    stop("UMAP CSV must contain columns ", args$`barcode-column`, " and ", args$`cluster-column`)
  }
  map_barcode_col <- if ("barcode" %in% names(map_df)) "barcode" else if ("cell_id" %in% names(map_df)) "cell_id" else NA_character_
  if (is.na(map_barcode_col)) stop("MapMyCells CSV must contain either 'barcode' or 'cell_id'")

  summary_df <- umap_df %>%
    transmute(barcode = .data[[args$`barcode-column`]], Cluster = .data[[args$`cluster-column`]]) %>%
    left_join(map_df %>% rename(barcode = all_of(map_barcode_col)), by = "barcode") %>%
    select(Cluster, all_of(required_levels)) %>%
    pivot_longer(cols = required_levels, names_to = "level", values_to = "designation") %>%
    group_by(Cluster, level, designation) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Cluster, level) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary_df, args$`output-csv`)
  message("Summary written to ", args$`output-csv`)
}

if (sys.nframe() == 0L) {
  main()
}
