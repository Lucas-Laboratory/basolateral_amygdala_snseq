#!/usr/bin/env Rscript

# Count the number of barcodes (columns) in a 10x Genomics `filtered_feature_bc_matrix.h5` file.
#
# Inputs
#   --input-h5: path to a 10x HDF5 file (required).
#
# Output
#   Prints the barcode count to stdout.
#
# Dependencies: optparse, Seurat, Matrix, hdf5r

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(hdf5r)
})

count_barcodes <- function(input_file) {
  message("Inspecting ", input_file)
  expr_matrix <- tryCatch(
    {
      matrix <- Read10X_h5(filename = input_file)
      if (is.list(matrix)) matrix[[1]] else matrix
    },
    error = function(err) {
      message("Failed to parse as 10x HDF5: ", err$message)
      NULL
    }
  )

  if (!is.null(expr_matrix)) {
    cat(ncol(expr_matrix), "\n")
    return(invisible(ncol(expr_matrix)))
  }

  message("Falling back to generic HDF5 inspection")
  h5file <- H5File$new(input_file, mode = "r")
  on.exit(h5file$close_all())
  print(h5file)
  stop("Unable to automatically extract an expression matrix. Inspect the dataset listing above and adapt the script for custom layouts.")
}

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "Path to 10x HDF5 file")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "count_barcodes_from_10x_h5.R --input-h5 FILE")
  args <- parse_args(parser)

  input_file <- args$`input-h5`
  if (is.null(input_file) || !nzchar(input_file)) stop("--input-h5 is required")
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)

  count_barcodes(input_file)
}

if (sys.nframe() == 0L) {
  main()
}
