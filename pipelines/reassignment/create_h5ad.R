#!/usr/bin/env Rscript

# Split a Seurat object into training/target subsets and export as .h5ad files for scANVI workflows.
#
# Inputs
#   --seurat-rds: Seurat object (.rds) containing RNA counts (required).
#   --output-dir: directory to write `seurat_train.h5ad` and `seurat_target.h5ad` (required).
#   --cluster-column: metadata column specifying clusters (default: `seurat_clusters`).
#   --target-clusters: comma-separated list of cluster ids to designate as target (required).
#   --assay: Seurat assay to export (default: `RNA`).
#
# Dependencies: optparse, Seurat, SingleCellExperiment, zellkonverter, Matrix

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Matrix)
})

write_subset <- function(seurat_obj, path) {
  sce <- as.SingleCellExperiment(seurat_obj)
  assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")
  assay(sce, "X") <- assay(sce, "counts")
  zellkonverter::writeH5AD(sce, path, compression = "gzip")
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--output-dir", type = "character", help = "Directory for .h5ad files"),
    make_option("--cluster-column", type = "character", default = "seurat_clusters"),
    make_option("--target-clusters", type = "character", help = "Comma-separated clusters"),
    make_option("--assay", type = "character", default = "RNA")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "create_h5ad.R --seurat-rds FILE --output-dir DIR --target-clusters A,B [options]")
  args <- parse_args(parser)

  required <- c(args$`seurat-rds`, args$`output-dir`, args$`target-clusters`)
  if (any(!nzchar(required))) stop("--seurat-rds, --output-dir, and --target-clusters are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found")

  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)
  seurat_obj <- readRDS(args$`seurat-rds`)
  if (!args$`cluster-column` %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found: ", args$`cluster-column`)
  }
  DefaultAssay(seurat_obj) <- args$assay
  Idents(seurat_obj) <- seurat_obj@meta.data[[args$`cluster-column`]]

  target_clusters <- trimws(strsplit(args$`target-clusters`, ",")[[1]])
  training_clusters <- setdiff(levels(Idents(seurat_obj)), target_clusters)
  if (!length(training_clusters)) stop("No training clusters remain after excluding targets")

  train_obj <- subset(seurat_obj, idents = training_clusters)
  target_obj <- subset(seurat_obj, idents = target_clusters)

  write_subset(train_obj, file.path(args$`output-dir`, "seurat_train.h5ad"))
  write_subset(target_obj, file.path(args$`output-dir`, "seurat_target.h5ad"))
}

if (sys.nframe() == 0L) {
  main()
}
