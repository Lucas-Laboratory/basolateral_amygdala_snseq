#!/usr/bin/env Rscript

# Perform cluster-wise differential expression for specified group comparisons on a Seurat object.
#
# Inputs
#   --seurat-rds: Seurat object saved as `.rds` (required).
#   --sample-column: metadata column containing sample/group labels (default: inferred from barcode prefix).
#   --output-dir: directory for differential expression results (required).
#   --comparisons: comma-separated list of comparisons in the form `GroupA:GroupB` (required).
#   --min-cells: minimum cells per cluster per group to perform the test (default: 3).
#   --test: statistical test passed to `FindMarkers` (default: `wilcox`).
#
# Output
#   One CSV per comparison with cluster-wise differential expression statistics.
#
# Dependencies: optparse, Seurat, dplyr, stringr

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(stringr)
})

annotate_samples <- function(seurat_obj, sample_column) {
  if (!is.null(sample_column) && sample_column %in% colnames(seurat_obj@meta.data)) {
    return(seurat_obj)
  }

  inferred <- str_split_fixed(Cells(seurat_obj), "_", 2)[, 1]
  seurat_obj$sample <- inferred
  message("Sample column not provided or absent; inferred labels from barcode prefixes into metadata column 'sample'.")
  seurat_obj
}

parse_comparisons <- function(comparison_string) {
  pairs <- str_split(comparison_string, pattern = ",", simplify = TRUE)
  pairs <- pairs[pairs != ""]
  if (!length(pairs)) stop("No comparisons provided.")

  lapply(pairs, function(entry) {
    groups <- str_split(entry, pattern = ":", simplify = TRUE)
    groups <- str_trim(groups)
    if (ncol(groups) < 2 || any(groups == "")) {
      stop("Malformed comparison: ", entry, " (expected format GroupA:GroupB)")
    }
    list(group1 = groups[1], group2 = groups[2])
  })
}

run_clusterwise_deg <- function(seurat_obj, group1, group2, sample_column, min_cells, test_use) {
  clusters <- levels(Idents(seurat_obj))
  results <- list()

  group1_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[sample_column]] == group1] 
  group2_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[sample_column]] == group2]

  for (cluster in clusters) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster)
    g1_cluster <- intersect(group1_cells, cluster_cells)
    g2_cluster <- intersect(group2_cells, cluster_cells)

    if (length(g1_cluster) >= min_cells && length(g2_cluster) >= min_cells) {
      markers <- FindMarkers(
        object = seurat_obj,
        ident.1 = g1_cluster,
        ident.2 = g2_cluster,
        min.pct = 0,
        min.diff.pct = 0,
        test.use = test_use,
        logfc.threshold = 0,
        only.pos = FALSE,
        verbose = FALSE
      )
      if (nrow(markers)) {
        markers$cluster <- cluster
        markers$comparison <- paste(group1, "vs", group2)
        markers$group1_cells <- length(g1_cluster)
        markers$group2_cells <- length(g2_cluster)
        results[[cluster]] <- markers
      }
    }
  }

  if (!length(results)) return(NULL)
  dplyr::bind_rows(results, .id = NULL)
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Path to Seurat .rds object"),
    make_option("--sample-column", type = "character", default = NULL,
                help = "Metadata column containing sample labels"),
    make_option("--output-dir", type = "character", help = "Directory for DEG outputs"),
    make_option("--comparisons", type = "character", help = "Comma-separated comparisons (GroupA:GroupB)"),
    make_option("--min-cells", type = "integer", default = 3,
                help = "Minimum cells per cluster per group [default %default]"),
    make_option("--test", type = "character", default = "wilcox",
                help = "Test passed to FindMarkers [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "run_clusterwise_differential_expression.R --seurat-rds FILE --output-dir DIR --comparisons A:B,C:D [options]")
  args <- parse_args(parser)

  seurat_path <- args$`seurat-rds`
  sample_column <- args$`sample-column`
  output_dir <- args$`output-dir`
  comparisons <- args$comparisons

  if (any(vapply(list(seurat_path, output_dir, comparisons), function(x) is.null(x) || !nzchar(x), logical(1)))) {
    stop("Arguments --seurat-rds, --output-dir, and --comparisons are required")
  }
  if (!file.exists(seurat_path)) stop("Seurat object not found: ", seurat_path)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(seurat_path)
  if (!is.null(sample_column) && !sample_column %in% colnames(seurat_obj@meta.data)) {
    stop("Specified sample column not present in metadata: ", sample_column)
  }
  seurat_obj <- annotate_samples(seurat_obj, sample_column)
  if (is.null(sample_column)) sample_column <- "sample"

  Idents(seurat_obj) <- factor(Idents(seurat_obj))

  parsed <- parse_comparisons(comparisons)

  for (cmp in parsed) {
    message("Running cluster-wise DEG for ", cmp$group1, " vs ", cmp$group2)
    deg <- run_clusterwise_deg(seurat_obj, cmp$group1, cmp$group2, sample_column, args$`min-cells`, args$test)
    if (is.null(deg)) {
      message("No qualifying clusters for ", cmp$group1, " vs ", cmp$group2)
      next
    }
    output_file <- file.path(output_dir, paste0("clusterwise_deg_", cmp$group1, "_vs_", cmp$group2, ".csv"))
    write.csv(deg, output_file, row.names = TRUE)
    message("Wrote results to ", output_file)
  }
}

if (sys.nframe() == 0L) {
  main()
}
