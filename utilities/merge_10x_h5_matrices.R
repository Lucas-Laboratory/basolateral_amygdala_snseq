#!/usr/bin/env Rscript

# Combine multiple 10x Genomics-formatted HDF5 count matrices into a single matrix.
#
# Inputs
#   --input-dir: directory containing 10x `filtered_feature_bc_matrix.h5` files (required).
#   --output-path: destination `.h5` filepath (required; parent directories created as needed).
#   --pattern: optional regular expression matcher for input filenames (default: `\\.h5$`).
#   --allow-mismatched-features: continue when feature sets differ by intersecting genes.
#
# Output
#   10x-compatible HDF5 file containing the merged feature/barcode matrix.
#
# Dependencies: optparse, rhdf5, Matrix

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(Matrix)
})

load_counts <- function(file_path) {
  message("Reading ", basename(file_path))
  datasets <- list(
    barcodes = "matrix/barcodes",
    feature_names = "matrix/features/name",
    feature_ids = "matrix/features/id",
    feature_types = "matrix/features/feature_type",
    data = "matrix/data",
    indices = "matrix/indices",
    indptr = "matrix/indptr",
    shape = "matrix/shape"
  )

  h5 <- lapply(datasets, function(path) {
    tryCatch(
      h5read(file_path, path),
      error = function(err) stop("File ", basename(file_path), " is missing dataset '", path, "': ", err$message)
    )
  })

  count_matrix <- sparseMatrix(
    i = as.integer(h5$indices) + 1L,
    p = as.integer(round(h5$indptr)),
    x = as.numeric(h5$data),
    dims = as.integer(h5$shape)
  )
  colnames(count_matrix) <- h5$barcodes
  rownames(count_matrix) <- h5$feature_names

  list(
    matrix = count_matrix,
    feature_ids = h5$feature_ids,
    feature_types = h5$feature_types,
    feature_names = h5$feature_names
  )
}

write_10x_h5 <- function(entry, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_path)) file.remove(output_path)

  h5createFile(output_path)
  h5createGroup(output_path, "matrix")
  h5createGroup(output_path, "matrix/features")

  h5write(entry$feature_ids, output_path, "matrix/features/id")
  h5write(entry$feature_names, output_path, "matrix/features/name")
  h5write(entry$feature_types, output_path, "matrix/features/feature_type")
  h5write(colnames(entry$matrix), output_path, "matrix/barcodes")
  h5write(entry$matrix@x, output_path, "matrix/data")
  h5write(entry$matrix@i, output_path, "matrix/indices")
  h5write(entry$matrix@p, output_path, "matrix/indptr")
  h5write(dim(entry$matrix), output_path, "matrix/shape")
}

merge_matrices <- function(entries, require_identical_features = TRUE) {
  feature_names <- lapply(entries, `[[`, "feature_names")
  shared <- Reduce(intersect, feature_names)

  if (!length(shared)) {
    stop("No overlapping features detected across the provided matrices.")
  }

  if (require_identical_features) {
    identical_sets <- vapply(feature_names, function(x) setequal(x, shared), logical(1))
    if (!all(identical_sets)) {
      stop("Feature lists differ across input files. Use --allow-mismatched-features to continue with their intersection.")
    }
  }

  lapply(entries, function(entry) {
    idx <- match(shared, entry$feature_names)
    entry$matrix <- entry$matrix[idx, , drop = FALSE]
    entry$feature_names <- entry$feature_names[idx]
    entry$feature_ids <- entry$feature_ids[idx]
    entry$feature_types <- entry$feature_types[idx]
    entry
  })
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory of 10x HDF5 files"),
    make_option("--output-path", type = "character", help = "Destination 10x HDF5 file"),
    make_option("--pattern", type = "character", default = "\\.h5$",
                help = "Regular expression used to match input files [default %default]"),
    make_option("--allow-mismatched-features", action = "store_true", default = FALSE,
                help = "Allow differing feature sets across files by intersecting them")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "merge_10x_h5_matrices.R --input-dir DIR --output-path FILE")
  args <- parse_args(parser)

  input_dir <- args$`input-dir`
  output_path <- args$`output-path`
  pattern <- args$pattern
  require_identical <- !args$`allow-mismatched-features`

  if (is.null(input_dir) || !nzchar(input_dir)) stop("--input-dir is required")
  if (is.null(output_path) || !nzchar(output_path)) stop("--output-path is required")
  if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)

  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No files matched pattern '", pattern, "' in ", input_dir)

  entries <- lapply(files, load_counts)
  aligned <- merge_matrices(entries, require_identical_features = require_identical)
  merged_matrix <- do.call(cbind, lapply(aligned, `[[`, "matrix"))

  base_entry <- aligned[[1]]
  base_entry$matrix <- merged_matrix

  write_10x_h5(base_entry, output_path)
  message("Merged ", length(files), " matrices into ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
