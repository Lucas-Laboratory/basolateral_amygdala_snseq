#!/usr/bin/env Rscript

# Create a ComplexHeatmap from FindAllMarkers log2FC matrix filtered by a marker list.
#
# Inputs
#   --logfc-csv: CSV of FindAllMarkers (columns `cluster`, `gene`, `avg_log2FC`) (required).
#   --markers-csv: CSV with `gene` and optional `cell_type`/`rank` columns (required).
#   --output-dir: directory for outputs (required).
#   --assay-name: optional label for heatmap legend (default: `avg_log2FC`).
#   --clusters: optional comma-separated cluster list/order (default: from data).
#   --width/--height: heatmap PDF size (default: 18 x 28).
#
# Outputs
#   - filtered/reordered matrix CSV
#   - dendrogram PDF
#   - ComplexHeatmap PDF
#
# Dependencies: optparse, data.table, dplyr, tidyr, ComplexHeatmap, circlize, dendextend

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
})

geom_mean_sparse <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(0)
  exp(mean(log(x)))
}

main <- function() {
  option_list <- list(
    make_option("--logfc-csv", type = "character", help = "FindAllMarkers CSV"),
    make_option("--markers-csv", type = "character", help = "Marker list CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--assay-name", type = "character", default = "avg_log2FC"),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--width", type = "double", default = 18),
    make_option("--height", type = "double", default = 28)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_log2fc_complex_heatmap.R --logfc-csv FILE --markers-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`logfc-csv`, args$`markers-csv`, args$`output-dir`)))) stop("All arguments required")
  if (!file.exists(args$`logfc-csv`)) stop("FindAllMarkers CSV not found")
  if (!file.exists(args$`markers-csv`)) stop("Marker CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  df_long <- fread(args$`logfc-csv`) %>% as_tibble() %>% mutate(gene = tolower(gene))
  markers <- fread(args$`markers-csv`) %>% as_tibble()
  names(markers) <- sub("^\ufeff", "", names(markers))
  markers <- markers %>% mutate(gene = tolower(gene))
  if (!"rank" %in% names(markers)) markers$rank <- seq_len(nrow(markers))
  if (!"cell_type" %in% names(markers)) markers$cell_type <- "marker_genes"
  markers <- markers %>% arrange(rank)

  df_wide <- df_long %>%
    select(cluster, gene, avg_log2FC) %>%
    pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = list(avg_log2FC = 0)) %>%
    column_to_rownames("gene")

  subset_data <- df_wide[rownames(df_wide) %in% markers$gene, , drop = FALSE]
  subset_data <- subset_data[match(markers$gene, rownames(subset_data)), ]
  subset_data <- subset_data %>% rownames_to_column("gene") %>%
    left_join(select(markers, gene, cell_type), by = "gene") %>%
    column_to_rownames("gene")

  write.csv(subset_data, file.path(args$`output-dir`, "filtered_reordered_matrix.csv"), row.names = TRUE)

  numeric_data <- subset_data[, setdiff(colnames(subset_data), "cell_type"), drop = FALSE]
  hc <- hclust(dist(t(numeric_data)), method = "ward.D2")
  pdf(file.path(args$`output-dir`, "cluster_dendrogram.pdf"), width = 10, height = 6)
  plot(as.dendrogram(hc), main = "Cluster Dendrogram")
  dev.off()

  ordered_clusters <- if (!is.null(args$clusters)) {
    trimws(strsplit(args$clusters, ",")[[1]])
  } else {
    hc$labels[hc$order]
  }
  numeric_data <- numeric_data[, ordered_clusters, drop = FALSE]
  subset_data <- subset_data[, c(ordered_clusters, "cell_type")]
  subset_data$cell_type <- factor(subset_data$cell_type, levels = unique(subset_data$cell_type))

  ordered_genes <- c()
  for (ct in levels(subset_data$cell_type)) {
    idx <- which(subset_data$cell_type == ct)
    sub_mat <- numeric_data[idx, , drop = FALSE]
    if (nrow(sub_mat) > 1) {
      hc_rows <- hclust(dist(sub_mat), method = "ward.D2")
      ordered_genes <- c(ordered_genes, rownames(sub_mat)[hc_rows$order])
    } else {
      ordered_genes <- c(ordered_genes, rownames(sub_mat))
    }
  }
  numeric_data <- numeric_data[ordered_genes, , drop = FALSE]
  subset_data <- subset_data[ordered_genes, ]

  unique_cell_types <- levels(subset_data$cell_type)
  ann_colors <- setNames(rainbow(length(unique_cell_types), s = 0.75, v = 0.85), unique_cell_types)
  row_ha <- rowAnnotation(cell_type = subset_data$cell_type,
                          col = list(cell_type = ann_colors),
                          show_annotation_name = TRUE,
                          annotation_name_side = "top")

  col_fun <- colorRamp2(c(-3.5, 0, 3.5), c("dodgerblue", "white", "firebrick"))
  pdf(file.path(args$`output-dir`, paste0("complex_heatmap_", args$`assay-name`, ".pdf")),
      width = args$width, height = args$height)
  draw(
    Heatmap(as.matrix(numeric_data), name = args$`assay-name`, col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE,
            row_split = subset_data$cell_type, row_gap = unit(10, "mm"),
            column_dend_height = unit(2, "cm"), show_row_names = TRUE,
            show_column_names = TRUE, row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6), left_annotation = row_ha,
            use_raster = FALSE)
  )
  dev.off()
}

if (sys.nframe() == 0L) {
  main()
}
