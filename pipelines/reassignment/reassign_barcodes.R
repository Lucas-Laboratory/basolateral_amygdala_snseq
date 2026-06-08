#!/usr/bin/env Rscript

# Reassign Seurat cluster identities using external predictions.
#
# Inputs
#   --barcode-csv: CSV with columns `Barcode` and `Cluster` (required).
#   --prediction-csv: CSV containing `barcode` and predicted cluster column (default `predicted_cluster`).
#   --seurat-rds: Seurat object (.rds) to update (required).
#   --output-rds: destination for updated Seurat object (default: `<input>_reintegrated.rds`).
#   --output-barcode-csv: optional path to write the merged barcode table.
#   --output-dir: optional directory to receive default RDS and merged barcode CSV outputs.
#   --prediction-column: column name holding predicted clusters [default `predicted_cluster`].
#
# Dependencies: optparse, Seurat, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

main <- function() {
  option_list <- list(
    make_option("--barcode-csv", type = "character", help = "Barcode cluster CSV"),
    make_option("--prediction-csv", type = "character", help = "Prediction CSV"),
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--output-rds", type = "character", default = NULL,
                help = "Output path for reintegrated Seurat object"),
    make_option("--output-barcode-csv", type = "character", default = NULL,
                help = "Optional merged barcode CSV"),
    make_option("--output-dir", type = "character", default = NULL,
                help = "Optional output directory for default outputs"),
    make_option("--prediction-column", type = "character", default = "predicted_cluster",
                help = "Prediction column name [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "reassign_barcodes.R --barcode-csv FILE --prediction-csv FILE --seurat-rds FILE [options]")
  args <- parse_args(parser)

  required <- c(args$`barcode-csv`, args$`prediction-csv`, args$`seurat-rds`)
  if (any(!nzchar(required))) stop("--barcode-csv, --prediction-csv, and --seurat-rds are required")
  if (!file.exists(args$`barcode-csv`)) stop("Barcode CSV not found")
  if (!file.exists(args$`prediction-csv`)) stop("Prediction CSV not found")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found")

  barcode_df <- read.csv(args$`barcode-csv`, stringsAsFactors = FALSE)
  if (!all(c("Barcode", "Cluster") %in% names(barcode_df))) {
    stop("Barcode CSV must contain 'Barcode' and 'Cluster' columns")
  }

  prediction_df <- read.csv(args$`prediction-csv`, stringsAsFactors = FALSE)
  if (!"barcode" %in% names(prediction_df)) {
    stop("Prediction CSV must contain a 'barcode' column")
  }
  if (!args$`prediction-column` %in% names(prediction_df)) {
    stop("Prediction CSV missing column: ", args$`prediction-column`)
  }
  prediction_df$barcode <- sub("-target$", "", prediction_df$barcode)
  prediction_df <- prediction_df %>% select(barcode, !!args$`prediction-column`)
  names(prediction_df) <- c("Barcode", "predicted_cluster")

  merged_bc <- barcode_df %>%
    left_join(prediction_df, by = "Barcode") %>%
    mutate(
      old_cluster = Cluster,
      new_cluster = ifelse(!is.na(predicted_cluster), predicted_cluster, Cluster)
    )

  seurat_obj <- readRDS(args$`seurat-rds`)
  orig_counts <- table(Idents(seurat_obj))

  mapping <- setNames(merged_bc$new_cluster, merged_bc$Barcode)
  present <- intersect(names(mapping), colnames(seurat_obj))
  if (!length(present)) stop("No overlapping barcodes between CSV and Seurat object")

  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(seurat_clusters = as.character(seurat_clusters))
  seurat_obj@meta.data[present, "seurat_clusters"] <- mapping[present]
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  new_counts <- table(Idents(seurat_obj))

  info_df <- data.frame(
    cluster = sort(unique(c(names(orig_counts), names(new_counts)))),
    before = as.integer(orig_counts[sort(unique(c(names(orig_counts), names(new_counts))))]),
    after = as.integer(new_counts[sort(unique(c(names(orig_counts), names(new_counts))))])
  )
  info_df[is.na(info_df)] <- 0L
  message("Cluster counts before/after reassignment:")
  print(info_df)

  output_dir <- args$`output-dir`
  if (!is.null(output_dir) && nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(args$`output-rds`) || !nzchar(args$`output-rds`)) {
      args$`output-rds` <- file.path(output_dir, "seurat_reintegrated.rds")
    }
    if (is.null(args$`output-barcode-csv`) || !nzchar(args$`output-barcode-csv`)) {
      args$`output-barcode-csv` <- file.path(output_dir, "reassigned_barcodes.csv")
    }
  }

  output_rds <- args$`output-rds`
  if (is.null(output_rds) || !nzchar(output_rds)) {
    output_rds <- file.path(dirname(args$`seurat-rds`),
                            paste0(tools::file_path_sans_ext(basename(args$`seurat-rds`)), "_reintegrated.rds"))
  }
  saveRDS(seurat_obj, output_rds)
  message("Saved updated Seurat object to ", output_rds)

  if (!is.null(args$`output-barcode-csv`) && nzchar(args$`output-barcode-csv`)) {
    dir.create(dirname(args$`output-barcode-csv`), recursive = TRUE, showWarnings = FALSE)
    write.csv(merged_bc, args$`output-barcode-csv`, row.names = FALSE)
  }
}

if (sys.nframe() == 0L) {
  main()
}
