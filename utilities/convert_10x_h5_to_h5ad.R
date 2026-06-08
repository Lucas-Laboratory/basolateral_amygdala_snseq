#!/usr/bin/env Rscript

# Convert a 10x Genomics-formatted HDF5 matrix (`filtered_feature_bc_matrix.h5`) into an AnnData `.h5ad` file.
#
# Inputs
#   --input-h5: path to the 10x HDF5 matrix (required).
#   --output-h5ad: destination `.h5ad` filepath (required).
#   --overwrite: replace the destination file when it already exists (default: FALSE).
#   --assay-name: optional name for the primary assay exported to AnnData (default: "X").
#
# Output
#   `.h5ad` file compatible with Scanpy/AnnData containing counts in both `X` and `counts` layers.
#
# Dependencies: optparse, Seurat, SeuratObject, zellkonverter, SingleCellExperiment

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

convert_h5_to_h5ad <- function(input_h5, output_h5ad, overwrite = FALSE, assay_name = "X") {
  if (!file.exists(input_h5)) {
    stop("Input file not found: ", input_h5)
  }
  if (file.exists(output_h5ad) && !overwrite) {
    stop("Output file already exists. Pass --overwrite to replace it: ", output_h5ad)
  }

  message("Reading 10x matrix from ", input_h5)
  counts <- Read10X_h5(input_h5)

  message("Creating Seurat object")
  seurat_obj <- CreateSeuratObject(counts = counts)

  message("Converting to SingleCellExperiment")
  sce <- as.SingleCellExperiment(seurat_obj)
  assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")
  assay(sce, assay_name) <- assay(sce, "counts")

  dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_h5ad)) file.remove(output_h5ad)

  message("Writing AnnData file to ", output_h5ad)
  zellkonverter::writeH5AD(sce, output_h5ad, compression = "gzip")
  invisible(output_h5ad)
}

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "Path to 10x HDF5 matrix"),
    make_option("--output-h5ad", type = "character", help = "Destination AnnData file"),
    make_option("--overwrite", action = "store_true", default = FALSE,
                help = "Overwrite output file if it exists"),
    make_option("--assay-name", type = "character", default = "X",
                help = "Assay name to populate in the AnnData object [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "convert_10x_h5_to_h5ad.R --input-h5 FILE --output-h5ad FILE [options]")
  args <- parse_args(parser)

  input_h5 <- args$`input-h5`
  output_h5ad <- args$`output-h5ad`
  overwrite <- args$overwrite
  assay_name <- args$`assay-name`

  if (is.null(input_h5) || !nzchar(input_h5)) stop("--input-h5 is required")
  if (is.null(output_h5ad) || !nzchar(output_h5ad)) stop("--output-h5ad is required")

  convert_h5_to_h5ad(input_h5, output_h5ad, overwrite = overwrite, assay_name = assay_name)
}

if (sys.nframe() == 0L) {
  main()
}
