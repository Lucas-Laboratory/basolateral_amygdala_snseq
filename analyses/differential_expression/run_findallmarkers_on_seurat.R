#!/usr/bin/env Rscript

# Run Seurat's FindAllMarkers on a serialized Seurat object and export the full table.
#
# Inputs
#   --seurat-rds: path to Seurat `.rds` object (required).
#   --output-csv: destination CSV for marker table (required).
#   --test: statistical test to pass to `FindAllMarkers` (default: `wilcox`).
#   --min-pct: minimum expression fraction threshold (default: 0).
#   --logfc-threshold: minimum log fold-change (default: 0).
#   --only-pos: return positive markers only (default: FALSE).
#
# Output
#   CSV containing the results of `FindAllMarkers`.
#
# Dependencies: optparse, Seurat

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Path to Seurat object"),
    make_option("--output-csv", type = "character", help = "Destination CSV"),
    make_option("--test", type = "character", default = "wilcox",
                help = "Test used by FindAllMarkers [default %default]"),
    make_option("--min-pct", type = "double", default = 0,
                help = "Minimum fraction of cells expressing a gene [default %default]"),
    make_option("--logfc-threshold", type = "double", default = 0,
                help = "Minimum log fold-change [default %default]"),
    make_option("--only-pos", action = "store_true", default = FALSE,
                help = "Return only positive markers")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "run_findallmarkers_on_seurat.R --seurat-rds FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  if (is.null(args$`seurat-rds`) || !nzchar(args$`seurat-rds`)) stop("--seurat-rds is required")
  if (is.null(args$`output-csv`) || !nzchar(args$`output-csv`)) stop("--output-csv is required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found: ", args$`seurat-rds`)

  seurat_obj <- readRDS(args$`seurat-rds`)

  markers <- FindAllMarkers(
    object = seurat_obj,
    test.use = args$test,
    min.pct = args$`min-pct`,
    logfc.threshold = args$`logfc-threshold`,
    only.pos = args$`only-pos`,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    return.thresh = 1,
    verbose = FALSE
  )

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(markers, file = args$`output-csv`, row.names = FALSE)
  message("FindAllMarkers table written to ", args$`output-csv`)
}

if (sys.nframe() == 0L) {
  main()
}
