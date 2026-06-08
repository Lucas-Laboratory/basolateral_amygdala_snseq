#!/usr/bin/env Rscript

# Create a ComplexHeatmap summarising geometric-mean SCT expression per sample/cluster for selected genes.
#
# Inputs
#   --seurat-rds: Seurat object (.rds) containing SCT assay (required).
#   --features-csv: CSV with column `feature` listing genes to include (required).
#   --output-pdf: destination PDF path (required).
#   --assay: Seurat assay to use (default: `SCT`).
#   --max-value: cap for expression values (default: 3.5).
#   --width/--height: PDF size (default: 12 x 14).
#   --cluster-column: metadata column for clusters (default: `seurat_clusters`).
#   --sample-column: optional metadata column for sample prefixes; defaults to barcode prefix before `_`.
#
# Dependencies: optparse, Seurat, dplyr, tidyr, Matrix, ComplexHeatmap, circlize

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
})

geom_mean_sparse <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(0)
  exp(mean(log(x)))
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--features-csv", type = "character", help = "CSV with column 'feature'"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--assay", type = "character", default = "SCT"),
    make_option("--max-value", type = "double", default = 3.5),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 14),
    make_option("--cluster-column", type = "character", default = "seurat_clusters"),
    make_option("--sample-column", type = "character", default = NULL)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_variable_features_heatmap.R --seurat-rds FILE --features-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`seurat-rds`, args$`features-csv`, args$`output-pdf`)))) stop("All inputs are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat file not found")
  if (!file.exists(args$`features-csv`)) stop("Feature CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(args$`seurat-rds`)
  features <- unique(read.csv(args$`features-csv`)$feature)
  features <- intersect(features, rownames(seurat_obj))
  if (!length(features)) stop("No features present in Seurat object")

  meta <- seurat_obj@meta.data
  meta$barcode <- colnames(seurat_obj)
  if (!args$`cluster-column` %in% names(meta)) stop("Cluster column not found")
  meta$cluster <- as.character(meta[[args$`cluster-column`]])
  if (is.null(args$`sample-column`) || !nzchar(args$`sample-column`)) {
    meta$sample <- sub("_.*", "", meta$barcode)
  } else {
    if (!args$`sample-column` %in% names(meta)) stop("Sample column not found")
    meta$sample <- meta[[args$`sample-column`]]
  }

  expr_matrix <- GetAssayData(seurat_obj, assay = args$assay, slot = "data")[features, , drop = FALSE]

  agg_expr <- Matrix::summary(expr_matrix) %>%
    as.data.frame() %>% rename(gene_index = i, barcode_index = j, expression = x) %>%
    mutate(gene = rownames(expr_matrix)[gene_index], barcode = colnames(expr_matrix)[barcode_index]) %>%
    select(gene, barcode, expression) %>%
    left_join(meta, by = "barcode") %>%
    group_by(gene, sample, cluster) %>%
    summarise(expression = geom_mean_sparse(expression), .groups = "drop") %>%
    pivot_wider(names_from = c(sample, cluster), values_from = expression)

  heatmap_matrix <- as.matrix(agg_expr[, -1])
  rownames(heatmap_matrix) <- agg_expr$gene
  heatmap_matrix[heatmap_matrix > args$`max-value`] <- args$`max-value`

  if (nrow(heatmap_matrix) < 2 || ncol(heatmap_matrix) < 2) stop("Insufficient data for clustering")

  col_names <- colnames(heatmap_matrix)
  ordered_cols <- col_names[order(sub("_.*", "", col_names), col_names)]
  heatmap_matrix <- heatmap_matrix[, ordered_cols]

  col_fun <- colorRamp2(c(0, args$`max-value`), c("white", "black"))

  pdf(args$`output-pdf`, width = args$width, height = args$height, useDingbats = FALSE)
  Heatmap(heatmap_matrix,
          cluster_rows = hclust(dist(heatmap_matrix)),
          cluster_columns = FALSE,
          col = col_fun,
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 4))
  dev.off()
  message("Heatmap saved to ", args$`output-pdf`)
}

if (sys.nframe() == 0L) {
  main()
}
