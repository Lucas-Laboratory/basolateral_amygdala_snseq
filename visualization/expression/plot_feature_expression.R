#!/usr/bin/env Rscript

# Generate Seurat feature plots for specified genes and save as PDFs.
#
# Inputs
#   --seurat-rds: Seurat object (.rds) containing embedding data (required).
#   --features: comma-separated list of gene features to plot (required).
#   --output-dir: directory for PDF plots (required).
#   --reduction: dimensional reduction to use (default: `umap`).
#   --palette: comma-separated colour vector for FeaturePlot (default: `gray95,firebrick`).
#   --width/--height: PDF width/height (default: 6).
#
# Dependencies: optparse, Seurat, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--features", type = "character", help = "Comma-separated gene list"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--reduction", type = "character", default = "umap"),
    make_option("--palette", type = "character", default = "gray95,firebrick"),
    make_option("--width", type = "double", default = 6),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_feature_expression.R --seurat-rds FILE --features GENE1,GENE2 --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`seurat-rds`, args$features, args$`output-dir`)))) {
    stop("--seurat-rds, --features, and --output-dir are required")
  }
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(args$`seurat-rds`)
  cols <- trimws(strsplit(args$palette, ",")[[1]])
  features <- trimws(strsplit(args$features, ",")[[1]])

  for (feature in features) {
    if (!(feature %in% rownames(seurat_obj))) {
      message("Skipping ", feature, " - not in Seurat object")
      next
    }
    plot_obj <- FeaturePlot(seurat_obj, features = feature, reduction = args$reduction,
                            order = TRUE, cols = cols) + theme_minimal()
    ggsave(file.path(args$`output-dir`, paste0("FeaturePlot_", feature, ".pdf")),
           plot_obj, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
