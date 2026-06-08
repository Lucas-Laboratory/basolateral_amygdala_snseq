#!/usr/bin/env Rscript

# Compute a dendrogram of clusters based on SCT expression across samples.
#
# Inputs
#   --seurat-rds: Seurat object (.rds) with SCT assay (required).
#   --output-dir: directory for PDF/CSV outputs (required).
#   --assay: assay name to use (default: `SCT`).
#   --sample-column: optional metadata column for sample; defaults to barcode prefix before `_`.
#   --cluster-column: metadata column for cluster IDs (default: `seurat_clusters`).
#   --width/--height: PDF size (default: 8 x 6).
#
# Outputs
#   - `cluster_dendrogram.pdf` and `cluster_order.csv` in output-dir.
#
# Dependencies: optparse, Seurat, dplyr, tidyr, Matrix, ggplot2, ggdendro

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(ggplot2)
  library(ggdendro)
})

geom_mean_sparse <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(0)
  exp(mean(log(x)))
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--assay", type = "character", default = "SCT"),
    make_option("--sample-column", type = "character", default = NULL),
    make_option("--cluster-column", type = "character", default = "seurat_clusters"),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_feature_dendrogram.R --seurat-rds FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`seurat-rds`, args$`output-dir`)))) stop("--seurat-rds and --output-dir are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(args$`seurat-rds`)
  features <- rownames(seurat_obj)

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

  expr_matrix <- GetAssayData(seurat_obj, assay = args$assay, slot = "data")
  agg_expr <- Matrix::summary(expr_matrix) %>%
    as.data.frame() %>% rename(gene_index = i, barcode_index = j, expression = x) %>%
    mutate(barcode = colnames(expr_matrix)[barcode_index]) %>% select(barcode, expression) %>%
    left_join(meta, by = "barcode") %>%
    group_by(cluster, sample) %>% summarise(expression = geom_mean_sparse(expression), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = expression)

  mat <- as.matrix(agg_expr[, -1])
  rownames(mat) <- agg_expr$cluster
  if (nrow(mat) < 2) stop("Insufficient clusters for dendrogram")

  hc <- hclust(dist(mat), method = "complete")
  cluster_order <- data.frame(Rank = seq_along(hc$order), Cluster = hc$labels[hc$order])
  write.csv(cluster_order, file.path(args$`output-dir`, "cluster_order.csv"), row.names = FALSE)

  pdf(file.path(args$`output-dir`, "cluster_dendrogram.pdf"), width = args$width, height = args$height)
  print(ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE) +
          theme_minimal() + theme(axis.text.y = element_text(size = 10)))
  dev.off()
}

if (sys.nframe() == 0L) {
  main()
}
