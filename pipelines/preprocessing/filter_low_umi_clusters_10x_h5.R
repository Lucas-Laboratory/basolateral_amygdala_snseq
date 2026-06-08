#!/usr/bin/env Rscript

# Remove clusters with low quality/UMI counts from a 10x HDF5 count matrix.
#
# Inputs
#   --input-h5: 10x HDF5 file to filter (required).
#   --cluster-csv: CSV linking barcodes to cluster assignments (required).
#   --clusters: comma-separated list of cluster labels to remove (required).
#   --output-h5: destination filtered HDF5 file (required).
#   --barcode-column: column containing barcodes in the cluster CSV (default `Barcode`).
#   --cluster-column: column containing cluster identifiers (default `Cluster`).
#
# Output
#   10x-formatted `.h5` file with specified clusters removed.
#
# Dependencies: optparse, rhdf5, Matrix, readr, dplyr, stringr

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(Matrix)
  library(readr)
  library(dplyr)
  library(stringr)
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
    make_option("--cluster-csv", type = "character", help = "CSV mapping barcodes to clusters"),
    make_option("--clusters", type = "character", help = "Comma-separated cluster labels to remove"),
    make_option("--output-h5", type = "character", help = "Destination HDF5 path"),
    make_option("--barcode-column", type = "character", default = "Barcode",
                help = "Barcode column in cluster CSV [default %default]"),
    make_option("--cluster-column", type = "character", default = "Cluster",
                help = "Cluster column in cluster CSV [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "filter_low_umi_clusters_10x_h5.R --input-h5 FILE --cluster-csv FILE --clusters A,B --output-h5 FILE")
  args <- parse_args(parser)

  input_h5 <- args$`input-h5`
  cluster_csv <- args$`cluster-csv`
  clusters <- args$clusters
  output_h5 <- args$`output-h5`
  barcode_col <- args$`barcode-column`
  cluster_col <- args$`cluster-column`

  if (any(vapply(list(input_h5, cluster_csv, clusters, output_h5), function(x) is.null(x) || !nzchar(x), logical(1)))) {
    stop("Arguments --input-h5, --cluster-csv, --clusters, and --output-h5 are required")
  }
  if (!file.exists(input_h5)) stop("Input HDF5 file not found: ", input_h5)
  if (!file.exists(cluster_csv)) stop("Cluster CSV not found: ", cluster_csv)

  remove_clusters <- str_split(clusters, pattern = ",", simplify = TRUE)
  remove_clusters <- as.vector(remove_clusters)
  remove_clusters <- trimws(remove_clusters[remove_clusters != ""])
  if (!length(remove_clusters)) stop("No clusters provided for removal")

  counts <- read_10x_h5(input_h5)
  cluster_table <- readr::read_csv(cluster_csv, show_col_types = FALSE)
  if (!all(c(barcode_col, cluster_col) %in% names(cluster_table))) {
    stop("Cluster CSV must contain columns '", barcode_col, "' and '", cluster_col, "'")
  }

  barcodes_to_remove <- cluster_table %>%
    filter(.data[[cluster_col]] %in% remove_clusters) %>%
    pull(.data[[barcode_col]])
  matched_barcodes_to_remove <- match_metadata_barcodes(barcodes_to_remove, colnames(counts$matrix))

  retain <- setdiff(colnames(counts$matrix), matched_barcodes_to_remove)
  filtered_matrix <- counts$matrix[, retain, drop = FALSE]

  message("Starting barcodes: ", ncol(counts$matrix))
  message("Metadata rows in clusters {", paste(remove_clusters, collapse = ","), "}: ", length(barcodes_to_remove))
  message("Matched and removed H5 barcodes: ", length(matched_barcodes_to_remove))
  message("Remaining barcodes: ", ncol(filtered_matrix))

  write_10x_h5(filtered_matrix, counts$features, colnames(filtered_matrix), output_h5)
  message("Filtered matrix written to ", output_h5)
}

if (sys.nframe() == 0L) {
  main()
}
