#!/usr/bin/env Rscript

# Compute local inverse Simpson's index (LISI) scores from a Seurat object and
# generate diagnostic plots.
#
# Inputs
#   --seurat-rds: Seurat object saved via saveRDS (required).
#   --output-dir: directory for CSV and PDF outputs (default: alongside input RDS).
#   --cluster-column: metadata column holding cluster labels (default: `seurat_clusters`).
#   --sample-column: metadata column defining batches/samples for LISI (default: `orig.ident`).
#   --group-column: optional metadata column for colouring a group UMAP (omit to skip).
#   --cluster-colors: optional CSV with columns `cluster` and `hex_color`.
#   --group-colors: optional CSV with columns `group` and `hex_color`.
#   --reduction: Seurat reduction to use for embeddings (default: `umap`).
#   --k: number of nearest neighbours for LISI (default: 30).
#   --seed: RNG seed used when shuffling points before plotting (default: 42).
#   --point-size: base point size for scatter plots (default: 0.5).
#   --width/--height: output dimensions in inches (default: 8 x 6).
#
# Outputs
#   - CSV of per-barcode LISI scores and embedding coordinates.
#   - PDF boxplot of LISI by cluster.
#   - PDF UMAP coloured by LISI.
#   - Optional PDF UMAP coloured by group column when requested.
#
# Dependencies: optparse, Seurat, FNN, dplyr, readr, ggplot2, scales

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(FNN)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tools)
})

compute_lisi <- function(knn_index, sample_ids) {
  neighbors <- knn_index$nn.index
  scores <- numeric(nrow(neighbors))
  for (i in seq_len(nrow(neighbors))) {
    neigh_ids <- sample_ids[neighbors[i, ]]
    freq <- table(neigh_ids)
    probs <- freq / sum(freq)
    scores[i] <- 1 / sum(probs^2)
  }
  scores
}

read_palette <- function(path, key_col, value_col) {
  df <- read_csv(path, show_col_types = FALSE)
  if (!all(c(key_col, value_col) %in% names(df))) {
    stop(sprintf("%s must contain columns '%s' and '%s'", path, key_col, value_col))
  }
  pal <- df[[value_col]]
  names(pal) <- as.character(df[[key_col]])
  pal
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Path to Seurat RDS"),
    make_option("--output-dir", type = "character", default = NULL,
                help = "Directory for outputs (default: alongside RDS)"),
    make_option("--cluster-column", type = "character", default = "seurat_clusters"),
    make_option("--sample-column", type = "character", default = "orig.ident"),
    make_option("--group-column", type = "character", default = NULL),
    make_option("--cluster-colors", type = "character", default = NULL,
                help = "Optional CSV mapping clusters to hex colours"),
    make_option("--group-colors", type = "character", default = NULL,
                help = "Optional CSV mapping groups to hex colours"),
    make_option("--reduction", type = "character", default = "umap"),
    make_option("--k", type = "integer", default = 30),
    make_option("--seed", type = "integer", default = 42),
    make_option("--point-size", type = "double", default = 0.5),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "compute_lisi_scores.R --seurat-rds FILE [options]")
  args <- parse_args(parser)

  if (!nzchar(args$`seurat-rds`)) stop("--seurat-rds is required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat RDS not found")

  output_dir <- args$`output-dir` %||% dirname(normalizePath(args$`seurat-rds`))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading Seurat object from ", args$`seurat-rds`)
  seu <- readRDS(args$`seurat-rds`)
  if (!inherits(seu, "Seurat")) stop("Input RDS must contain a Seurat object")

  meta <- seu@meta.data
  cluster_col <- args$`cluster-column`
  sample_col <- args$`sample-column`

  if (!cluster_col %in% colnames(meta)) {
    stop("Cluster column ", cluster_col, " not found in Seurat metadata")
  }
  if (!sample_col %in% colnames(meta)) {
    stop("Sample column ", sample_col, " not found in Seurat metadata")
  }

  group_col <- args$`group-column`
  if (!is.null(group_col) && !nzchar(group_col)) group_col <- NULL
  if (!is.null(group_col) && !group_col %in% colnames(meta)) {
    stop("Group column ", group_col, " not found in Seurat metadata")
  }

  embeddings <- tryCatch({
    Embeddings(seu, reduction = args$reduction)
  }, error = function(err) NULL)
  if (is.null(embeddings)) stop("Reduction ", args$reduction, " not found in Seurat object")
  embeddings <- as.matrix(embeddings)

  if (nrow(embeddings) != nrow(meta)) {
    stop("Embedding row count does not match metadata rows")
  }

  k <- args$k
  if (k <= 0 || k >= nrow(embeddings)) {
    stop("--k must be between 1 and number of cells - 1")
  }

  message("Computing LISI scores with k=", k)
  knn <- FNN::get.knn(embeddings, k = k)
  lisi_values <- compute_lisi(knn, meta[[sample_col]])

  coords <- as.data.frame(embeddings[, 1:min(ncol(embeddings), 2), drop = FALSE])
  colnames(coords)[1:2] <- c("UMAP_1", "UMAP_2")

  set.seed(args$seed)
  df <- meta %>%
    mutate(barcode = rownames(meta),
           cluster = .data[[cluster_col]],
           sample = .data[[sample_col]],
           lisi_score = lisi_values,
           UMAP_1 = coords$UMAP_1,
           UMAP_2 = coords$UMAP_2) %>%
    sample_frac(size = 1)

  base_name <- file_path_sans_ext(basename(args$`seurat-rds`))
  csv_path <- file.path(output_dir, paste0("lisi_scores_", base_name, ".csv"))
  write_csv(df, csv_path)
  message("Wrote LISI scores to ", csv_path)

  cluster_palette <- if (!is.null(args$`cluster-colors`) && nzchar(args$`cluster-colors`)) {
    read_palette(args$`cluster-colors`, "cluster", "hex_color")
  } else {
    pal <- hue_pal()(length(unique(df$cluster)))
    names(pal) <- sort(unique(df$cluster))
    pal
  }

  boxplot_path <- file.path(output_dir, paste0("plot_lisi_by_cluster_", base_name, ".pdf"))
  p_box <- ggplot(df, aes(x = factor(cluster), y = lisi_score, fill = factor(cluster))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.3, alpha = 0.2, size = 0.6) +
    scale_fill_manual(values = cluster_palette) +
    theme_minimal(base_size = 11) +
    labs(title = "LISI by cluster", x = "Cluster", y = "LISI score", fill = "Cluster")
  ggsave(boxplot_path, p_box, width = args$width, height = args$height)
  message("Wrote LISI boxplot to ", boxplot_path)

  umap_lisi_path <- file.path(output_dir, paste0("plot_umap_lisi_", base_name, ".pdf"))
  p_umap <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = lisi_score)) +
    geom_point(size = args$`point-size`, alpha = 0.8) +
    scale_colour_viridis_c(option = "magma") +
    theme_minimal(base_size = 11) +
    labs(title = "UMAP coloured by LISI", x = "UMAP 1", y = "UMAP 2", colour = "LISI")
  ggsave(umap_lisi_path, p_umap, width = args$width, height = args$height)
  message("Wrote UMAP LISI plot to ", umap_lisi_path)

  if (!is.null(group_col)) {
    group_palette <- if (!is.null(args$`group-colors`) && nzchar(args$`group-colors`)) {
      read_palette(args$`group-colors`, "group", "hex_color")
    } else {
      pal <- hue_pal()(length(unique(df[[group_col]])))
      names(pal) <- sort(unique(df[[group_col]]))
      pal
    }
    group_path <- file.path(output_dir, paste0("plot_umap_", group_col, "_", base_name, ".pdf"))
    p_group <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = .data[[group_col]])) +
      geom_point(size = args$`point-size`, alpha = 0.8) +
      scale_colour_manual(values = group_palette) +
      theme_minimal(base_size = 11) +
      labs(title = paste("UMAP coloured by", group_col), x = "UMAP 1", y = "UMAP 2", colour = group_col)
    ggsave(group_path, p_group, width = args$width, height = args$height)
    message("Wrote UMAP group plot to ", group_path)
  }

  message("Done computing LISI outputs")
}

`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

if (sys.nframe() == 0L) {
  main()
}
