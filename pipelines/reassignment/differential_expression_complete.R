#!/usr/bin/env Rscript

# Convenience wrapper around cluster-wise differential expression comparisons for VAE reassignment workflows.
#
# Inputs
#   --seurat-rds: reintegrated Seurat object (required).
#   --output-dir: directory for DEG tables (required).
#   --sample-column: metadata column with sample labels (default: inferred from barcode prefix).
#   --comparisons: comma-separated list like `GroupA:GroupB,GroupC:GroupD` (default matches legacy pairs).
#   --min-cells: minimum cells per cluster per group (default: 3).
#   --min-pct: Seurat `min.pct` parameter (default: 0.1).
#   --logfc-threshold: Seurat `logfc.threshold` (default: 0).
#
# Dependencies: optparse, Seurat, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

parse_comparisons <- function(spec) {
  comps <- trimws(strsplit(spec, ",")[[1]])
  lapply(comps, function(pair) {
    parts <- trimws(strsplit(pair, ":")[[1]])
    if (length(parts) != 2) stop("Invalid comparison: ", pair)
    list(a = parts[1], b = parts[2])
  })
}

run_clusterwise_deg <- function(seurat_obj, cells_a, cells_b, name_a, name_b, min_cells, min_pct, logfc_threshold, output_dir) {
  deg_dir <- file.path(output_dir, "DEG")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  Idents(seurat_obj) <- "seurat_clusters"
  results <- list()

  for (cl in levels(Idents(seurat_obj))) {
    sel_a <- intersect(cells_a, WhichCells(seurat_obj, idents = cl))
    sel_b <- intersect(cells_b, WhichCells(seurat_obj, idents = cl))
    if (length(sel_a) >= min_cells && length(sel_b) >= min_cells) {
      markers <- FindMarkers(seurat_obj,
                             ident.1 = sel_a,
                             ident.2 = sel_b,
                             min.pct = min_pct,
                             logfc.threshold = logfc_threshold,
                             verbose = FALSE)
      if (nrow(markers)) {
        markers$cluster <- cl
        markers$comparison <- paste(name_a, "vs", name_b)
        results[[cl]] <- markers
      }
    }
  }

  if (length(results)) {
    combined <- do.call(rbind, results)
    out_path <- file.path(deg_dir, paste0("merged_clusterwise_DEG_", name_a, "_vs_", name_b, ".csv"))
    write.csv(combined, out_path)
    message("Wrote ", basename(out_path))
  } else {
    message("No DEGs for ", name_a, " vs ", name_b)
  }
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--sample-column", type = "character", default = NULL),
    make_option("--comparisons", type = "character",
                default = paste(c("Female-Diestrus-Naive:Male-Naive",
                                  "Female-Proestrus-Naive:Male-Naive",
                                  "Female-Proestrus-Naive:Female-Diestrus-Naive",
                                  "Male-Naive:Female-Diestrus-Naive",
                                  "Male-Naive:Female-Proestrus-Naive",
                                  "Female-Diestrus-Naive:Female-Proestrus-Naive"), collapse = ",")),
    make_option("--min-cells", type = "integer", default = 3),
    make_option("--min-pct", type = "double", default = 0.1),
    make_option("--logfc-threshold", type = "double", default = 0)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "differential_expression_complete.R --seurat-rds FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`seurat-rds`, args$`output-dir`)))) stop("--seurat-rds and --output-dir are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(args$`seurat-rds`)
  if (is.null(args$`sample-column`) || !nzchar(args$`sample-column`)) {
    seurat_obj$sample <- sapply(strsplit(Cells(seurat_obj), "_"), `[`, 1)
  } else {
    if (!args$`sample-column` %in% colnames(seurat_obj@meta.data)) {
      stop("Sample column not found: ", args$`sample-column`)
    }
    seurat_obj$sample <- seurat_obj@meta.data[[args$`sample-column`]]
  }

  comparisons <- parse_comparisons(args$comparisons)
  for (cmp in comparisons) {
    cells_a <- WhichCells(seurat_obj, expression = sample == cmp$a)
    cells_b <- WhichCells(seurat_obj, expression = sample == cmp$b)
    run_clusterwise_deg(seurat_obj, cells_a, cells_b, cmp$a, cmp$b,
                        args$`min-cells`, args$`min-pct`, args$`logfc-threshold`, args$`output-dir`)
  }
}

if (sys.nframe() == 0L) {
  main()
}
