#!/usr/bin/env Rscript

# Perform classical multidimensional scaling on a symmetric similarity matrix and export coordinates.
#
# Inputs
#   --weights-csv: CSV containing a symmetric matrix with matching row/column names (required).
#   --output-csv: path to write coordinates (required).
#   --dimension: number of MDS dimensions to retain (default: 1).
#
# Output
#   CSV with columns `gene` and `position` (or `dim_1`... for higher dimensions).
#
# Dependencies: optparse, readr, tibble

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(tibble)
})

main <- function() {
  option_list <- list(
    make_option("--weights-csv", type = "character", help = "Symmetric weight matrix CSV"),
    make_option("--output-csv", type = "character", help = "Destination CSV"),
    make_option("--dimension", type = "integer", default = 1,
                help = "Number of MDS dimensions [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step03_cmds_sort_by_connection.R --weights-csv FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`weights-csv`, args$`output-csv`)))) {
    stop("--weights-csv and --output-csv are required")
  }
  if (!file.exists(args$`weights-csv`)) stop("Weights CSV not found: ", args$`weights-csv`)
  if (args$dimension < 1) stop("--dimension must be >= 1")

  mat <- as.matrix(read.csv(args$`weights-csv`, row.names = 1, check.names = FALSE))

  if (!all(rownames(mat) == colnames(mat))) {
    stop("Row names and column names must match")
  }
  if (!isTRUE(all.equal(mat, t(mat)))) {
    stop("Matrix is not symmetric")
  }

  max_sim <- max(mat, na.rm = TRUE)
  diss_mat <- max_sim - mat
  coords <- cmdscale(as.dist(diss_mat), k = args$dimension, eig = TRUE)$points

  if (args$dimension == 1) {
    result <- tibble(gene = rownames(coords), position = coords[, 1])
    result <- result[order(result$position), , drop = FALSE]
  } else {
    coord_df <- as.data.frame(coords)
    names(coord_df) <- paste0("dim_", seq_len(ncol(coord_df)))
    coord_df$gene <- rownames(coords)
    result <- coord_df[, c("gene", names(coord_df)[names(coord_df) != "gene"])]
  }

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  write_csv(result, args$`output-csv`)
}

if (sys.nframe() == 0L) {
  main()
}
