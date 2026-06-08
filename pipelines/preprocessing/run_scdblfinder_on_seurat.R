#!/usr/bin/env Rscript

# Batch-run scDblFinder on Seurat objects stored as `.rds` files and export QC summaries.
#
# Inputs
#   --input-dir: directory containing Seurat objects serialized as `.rds` files (required).
#   --output-dir: destination for scDblFinder reports (required).
#   --pattern: regular expression to select input files (default: `\\.rds$`).
#   --assay: Seurat assay to convert for scDblFinder (default: `RNA`).
#   --reduction: dimensional reduction used for plotting (default: `umap`).
#   --tsne-reduction: optional secondary reduction for plotting (default: `tsne`).
#   --clusters: metadata column used for boxplots (default: `seurat_clusters`).
#   --seed: optional random seed to ensure reproducibility.
#
# Outputs (per input object)
#   - scDblFinder doublet metadata (`doublet_metadata_<sample>.csv`).
#   - Histogram of doublet scores.
#   - UMAP/tSNE visualisations coloured by doublet class.
#   - Boxplot of doublet scores per cluster.
#   - Scatterplot of nUMI versus doublet score.
#
# Dependencies: optparse, Seurat, scDblFinder, SingleCellExperiment, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
})

run_scdblfinder <- function(seurat_obj, assay) {
  DefaultAssay(seurat_obj) <- assay
  sce <- as.SingleCellExperiment(seurat_obj)
  scDblFinder(sce, clusters = FALSE)
}

render_histogram <- function(scores, path) {
  pdf(path)
  hist(scores, breaks = 50, main = "scDblFinder Score Distribution",
       xlab = "Doublet Score", col = "lightblue", border = "grey40")
  dev.off()
}

render_dimplot <- function(seurat_obj, reduction, column, title, path) {
  if (!reduction %in% names(seurat_obj@reductions)) {
    warning("Skipping plot: reduction not found - ", reduction)
    return(invisible(NULL))
  }
  pdf(path)
  print(DimPlot(seurat_obj, reduction = reduction, group.by = column) +
          ggtitle(title) + theme_minimal())
  dev.off()
}

render_boxplot <- function(scores, clusters, path) {
  if (is.data.frame(clusters)) clusters <- clusters[[1]]
  if (is.list(clusters)) clusters <- unlist(clusters, use.names = FALSE)
  pdf(path)
  boxplot(scores ~ clusters,
          main = "Doublet Scores per Cluster",
          xlab = "Cluster", ylab = "Doublet Score",
          col = "lightgrey")
  dev.off()
}

render_scatter <- function(numi, scores, path) {
  pdf(path)
  plot(numi, scores,
       pch = 16, col = rgb(0, 0, 0, 0.4),
       xlab = "nUMI", ylab = "Doublet Score",
       main = "Doublet Score vs nUMI")
  dev.off()
}

process_file <- function(file_path, output_dir, assay, reduction, tsne_reduction, cluster_column) {
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  message("Processing ", sample_name)

  seurat_obj <- readRDS(file_path)
  sce <- run_scdblfinder(seurat_obj, assay)

  doublet_df <- data.frame(
    barcode = colnames(sce),
    doublet_score = sce$scDblFinder.score,
    doublet_class = sce$scDblFinder.class,
    stringsAsFactors = FALSE
  )

  csv_path <- file.path(output_dir, paste0("doublet_metadata_", sample_name, ".csv"))
  write.csv(doublet_df, csv_path, row.names = FALSE)

  seurat_obj$doublet_class <- factor(sce$scDblFinder.class, levels = c("singlet", "doublet"))
  seurat_obj$doublet_score <- sce$scDblFinder.score

  render_histogram(sce$scDblFinder.score, file.path(output_dir, paste0("doublet_score_distribution_", sample_name, ".pdf")))
  render_dimplot(seurat_obj, reduction, "doublet_class", "Doublet Classification", file.path(output_dir, paste0(reduction, "_doublet_classification_", sample_name, ".pdf")))
  if (!is.null(tsne_reduction) && nzchar(tsne_reduction)) {
    render_dimplot(seurat_obj, tsne_reduction, "doublet_class", "Doublet Classification", file.path(output_dir, paste0(tsne_reduction, "_doublet_classification_", sample_name, ".pdf")))
  }
  if (cluster_column %in% colnames(seurat_obj@meta.data)) {
    render_boxplot(seurat_obj$doublet_score, seurat_obj[[cluster_column]], file.path(output_dir, paste0("boxplot_doublet_scores_", sample_name, ".pdf")))
  }
  if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    render_scatter(seurat_obj$nCount_RNA, seurat_obj$doublet_score, file.path(output_dir, paste0("scatter_nUMI_vs_doublet_", sample_name, ".pdf")))
  }
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory with Seurat .rds files"),
    make_option("--output-dir", type = "character", help = "Directory for scDblFinder outputs"),
    make_option("--pattern", type = "character", default = "\\.rds$",
                help = "Regex to select input files [default %default]"),
    make_option("--assay", type = "character", default = "RNA",
                help = "Default assay to use when converting to SCE [default %default]"),
    make_option("--reduction", type = "character", default = "umap",
                help = "Dimensional reduction for plotting [default %default]"),
    make_option("--tsne-reduction", type = "character", default = "tsne",
                help = "Optional secondary reduction for plotting"),
    make_option("--clusters", type = "character", default = "seurat_clusters",
                help = "Metadata column for cluster-level summaries [default %default]"),
    make_option("--seed", type = "integer", default = NULL,
                help = "Random seed for reproducibility")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "run_scdblfinder_on_seurat.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  input_dir <- args$`input-dir`
  output_dir <- args$`output-dir`
  if (is.null(input_dir) || !nzchar(input_dir)) stop("--input-dir is required")
  if (is.null(output_dir) || !nzchar(output_dir)) stop("--output-dir is required")
  if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(args$seed)) set.seed(args$seed)

  files <- list.files(input_dir, pattern = args$pattern, full.names = TRUE)
  if (!length(files)) stop("No files matched pattern '", args$pattern, "' in ", input_dir)

  lapply(files, process_file,
         output_dir = output_dir,
         assay = args$assay,
         reduction = args$reduction,
         tsne_reduction = args$`tsne-reduction`,
         cluster_column = args$clusters)

  message("Processed ", length(files), " objects. Outputs located in ", output_dir)
}

if (sys.nframe() == 0L) {
  main()
}
