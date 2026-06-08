#!/usr/bin/env Rscript

# Merge multiple 10x HDF5 subsets, perform SCTransform, and export clustering results.
#
# Inputs
#   --input-h5: comma-separated list of 10x HDF5 files to merge (required).
#   --output-dir: directory for merged Seurat object outputs (required).
#   --assay: assay name for Seurat object (default: `RNA`).
#   --variable-features: number of variable features (default: 3000).
#   --dims: number of PCA dimensions (default: 50).
#   --k-param / --resolution / --metric / --algorithm: clustering parameters (defaults: 30, 0.4, cosine, 2).
#   --umap-min-dist / --umap-neighbors / --umap-repulsion: UMAP settings (defaults: 0.4, 50, 0.9).
#   --future-max-gb: maximum future globals size in GB for SCTransform (default: 8).
#
# Dependencies: optparse, Seurat, dplyr, ggplot2, future, reticulate

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(reticulate)
})

ensure_glm_gampoi <- function() {
  if (requireNamespace("glmGamPoi", quietly = TRUE)) return(invisible(TRUE))
  message("Installing Bioconductor package glmGamPoi for faster SCTransform.")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("glmGamPoi", ask = FALSE, update = FALSE)
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop("glmGamPoi installation failed; install it manually with BiocManager::install('glmGamPoi')")
  }
  invisible(TRUE)
}

ensure_umap_learn <- function() {
  if (exists("py_require", envir = asNamespace("reticulate"), inherits = FALSE)) {
    reticulate::py_require("umap-learn")
  }
  if (!reticulate::py_module_available("umap")) {
    message("Installing Python package umap-learn for Seurat RunUMAP.")
    reticulate::py_install("umap-learn", pip = TRUE)
  }
  if (!reticulate::py_module_available("umap")) {
    stop("umap-learn installation failed; install it manually with python -m pip install umap-learn")
  }
  invisible(TRUE)
}

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "Comma-separated list of 10x HDF5 files"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--assay", type = "character", default = "RNA"),
    make_option("--variable-features", type = "integer", default = 3000),
    make_option("--dims", type = "integer", default = 50),
    make_option("--k-param", type = "integer", default = 30),
    make_option("--resolution", type = "double", default = 0.4),
    make_option("--metric", type = "character", default = "cosine"),
    make_option("--algorithm", type = "integer", default = 2),
    make_option("--umap-min-dist", type = "double", default = 0.4),
    make_option("--umap-neighbors", type = "integer", default = 50),
    make_option("--umap-repulsion", type = "double", default = 0.9),
    make_option("--future-max-gb", type = "double", default = 8,
                help = "Maximum future globals size in GB for SCTransform [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "merged_subsets_sctransform.R --input-h5 file1.h5,file2.h5 --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-h5`, args$`output-dir`)))) {
    stop("--input-h5 and --output-dir are required")
  }
  h5_files <- trimws(strsplit(args$`input-h5`, ",")[[1]])
  if (!all(file.exists(h5_files))) stop("One or more input H5 files not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)
  options(future.globals.maxSize = args$`future-max-gb` * 1024^3)
  ensure_glm_gampoi()

  seurat_list <- lapply(h5_files, function(path) {
    Read10X_h5(path, use.names = TRUE) %>% CreateSeuratObject(assay = args$assay)
  })
  merged <- merge(seurat_list[[1]], seurat_list[-1])

  merged <- SCTransform(merged, variable.features.n = args$`variable-features`, verbose = TRUE)
  dims_vec <- seq_len(args$dims)
  merged <- RunPCA(merged, npcs = max(args$dims, 50), verbose = TRUE)
  merged <- FindNeighbors(merged, dims = dims_vec, k.param = args$`k-param`,
                          nn.method = "annoy", annoy.metric = args$metric)
  merged <- FindClusters(merged, algorithm = args$algorithm, resolution = args$resolution)
  ensure_umap_learn()
  merged <- RunUMAP(merged, dims = dims_vec, min.dist = args$`umap-min-dist`, metric = args$metric,
                    umap.method = "umap-learn", n.neighbors = args$`umap-neighbors`,
                    repulsion.strength = args$`umap-repulsion`)

  base_name <- paste0("merged_", length(h5_files), "subset")
  ggsave(file.path(args$`output-dir`, paste0("plot-UMAP_", base_name, ".pdf")),
         DimPlot(merged, label = TRUE) + theme_minimal())

  umap_coords <- Embeddings(merged, "umap")
  barcode_df <- data.frame(Barcode = rownames(umap_coords), Cluster = Idents(merged),
                           UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2])
  write.csv(barcode_df, file.path(args$`output-dir`, paste0("dataframe-BarcodeClusterUMAP_", base_name, ".csv")),
            row.names = FALSE)

  SaveH5Seurat(merged, filename = file.path(args$`output-dir`, paste0(base_name, ".h5seurat")), overwrite = TRUE)
}

if (sys.nframe() == 0L) {
  main()
}
