#!/usr/bin/env Rscript

# Remove barcodes classified as doublets by scDblFinder and write a filtered 10x HDF5 matrix.
#
# Inputs
#   --input-h5: 10x `filtered_feature_bc_matrix.h5` (or similar) (required).
#   --barcode-csv: CSV with columns `barcode` and `doublet_class` (required).
#   --output-h5: path for the filtered 10x HDF5 file (required).
#   --doublet-class: label indicating doublets within `doublet_class` column (default: `doublet`).
#
# Output
#   10x-formatted `.h5` file excluding barcodes labelled as doublets.
#
# Dependencies: optparse, rhdf5, Matrix, readr, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(Matrix)
  library(readr)
  library(dplyr)
})

read_10x_h5 <- function(path) {
  barcodes <- h5read(path, "matrix/barcodes")
  data <- h5read(path, "matrix/data")
  indices <- h5read(path, "matrix/indices")
  indptr <- h5read(path, "matrix/indptr")
  shape <- h5read(path, "matrix/shape")

  count_matrix <- sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = as.integer(shape)
  )
  rownames(count_matrix) <- h5read(path, "matrix/features/name")
  colnames(count_matrix) <- barcodes

  feature_meta <- list(
    feature_type = h5read(path, "matrix/features/feature_type"),
    feature_id = h5read(path, "matrix/features/id"),
    feature_name = h5read(path, "matrix/features/name")
  )

  list(matrix = count_matrix, features = feature_meta)
}

write_10x_h5 <- function(count_matrix, feature_meta, barcodes, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_path)) file.remove(output_path)

  h5createFile(output_path)
  h5createGroup(output_path, "matrix")
  h5createGroup(output_path, "matrix/features")

  h5write(feature_meta$feature_type, output_path, "matrix/features/feature_type")
  h5write(feature_meta$feature_id, output_path, "matrix/features/id")
  h5write(feature_meta$feature_name, output_path, "matrix/features/name")

  h5write(barcodes, output_path, "matrix/barcodes")
  h5write(count_matrix@x, output_path, "matrix/data")
  h5write(count_matrix@i, output_path, "matrix/indices")
  h5write(count_matrix@p, output_path, "matrix/indptr")
  h5write(dim(count_matrix), output_path, "matrix/shape")
}

match_metadata_barcodes <- function(metadata_barcodes, matrix_barcodes) {
  metadata_barcodes <- unique(as.character(metadata_barcodes))
  metadata_barcodes <- metadata_barcodes[!is.na(metadata_barcodes) & nzchar(metadata_barcodes)]
  direct <- intersect(metadata_barcodes, matrix_barcodes)

  stripped <- sub("^.*_", "", metadata_barcodes)
  stripped_matches <- intersect(stripped, matrix_barcodes)

  matched <- unique(c(direct, stripped_matches))
  if (!length(matched) && length(metadata_barcodes)) {
    warning("No metadata barcodes matched H5 barcodes. If metadata is sample-prefixed, use 'Sample_BARCODE' with an underscore separator.")
  } else if (length(stripped_matches) && length(direct) != length(matched)) {
    message("Matched ", length(stripped_matches), " sample-prefixed metadata barcodes after stripping the prefix before '_'.")
  }
  matched
}

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "Path to 10x HDF5 matrix"),
    make_option("--barcode-csv", type = "character", help = "CSV with scDblFinder classifications"),
    make_option("--output-h5", type = "character", help = "Destination HDF5 path"),
    make_option("--doublet-class", type = "character", default = "doublet",
                help = "Label indicating doublets [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "remove_doublets_scdblfinder_10x_h5.R --input-h5 FILE --barcode-csv FILE --output-h5 FILE")
  args <- parse_args(parser)

  input_h5 <- args$`input-h5`
  barcode_csv <- args$`barcode-csv`
  output_h5 <- args$`output-h5`
  doublet_label <- args$`doublet-class`

  if (any(vapply(list(input_h5, barcode_csv, output_h5), function(x) is.null(x) || !nzchar(x), logical(1)))) {
    stop("Arguments --input-h5, --barcode-csv, and --output-h5 are required")
  }
  if (!file.exists(input_h5)) stop("Input HDF5 file not found: ", input_h5)
  if (!file.exists(barcode_csv)) stop("Barcode CSV not found: ", barcode_csv)

  counts <- read_10x_h5(input_h5)
  barcode_table <- readr::read_csv(barcode_csv, show_col_types = FALSE)

  if (!all(c("barcode", "doublet_class") %in% names(barcode_table))) {
    stop("Barcode CSV must contain columns 'barcode' and 'doublet_class'")
  }

  doublet_barcodes <- barcode_table %>%
    filter(doublet_class == doublet_label) %>%
    pull(barcode)
  matched_doublet_barcodes <- match_metadata_barcodes(doublet_barcodes, colnames(counts$matrix))

  retain <- setdiff(colnames(counts$matrix), matched_doublet_barcodes)
  filtered_matrix <- counts$matrix[, retain, drop = FALSE]

  message("Starting barcodes: ", ncol(counts$matrix))
  message("Metadata rows labelled as doublets: ", length(doublet_barcodes))
  message("Matched and removed H5 doublet barcodes: ", length(matched_doublet_barcodes))
  message("Remaining barcodes: ", ncol(filtered_matrix))

  write_10x_h5(filtered_matrix, counts$features, colnames(filtered_matrix), output_h5)
  message("Filtered matrix written to ", output_h5)
}

if (sys.nframe() == 0L) {
  main()
}
