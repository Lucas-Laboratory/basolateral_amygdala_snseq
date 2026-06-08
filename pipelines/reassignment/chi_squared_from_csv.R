#!/usr/bin/env Rscript

# Render a heatmap from a chi-squared residual CSV file.
#
# Inputs
#   --input-csv: path to CSV with rownames (required).
#   --output-pdf: destination PDF (default: replace .csv with .pdf).
#   --width: PDF width (default: 6).
#   --height: PDF height (default: 8).
#   --low/--mid/--high: colours for gradient (defaults: dodgerblue3/white/firebrick3).
#
# Dependencies: optparse, readr, ggplot2, reshape2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
  library(reshape2)
})

main <- function() {
  option_list <- list(
    make_option("--input-csv", type = "character", help = "Residual CSV"),
    make_option("--output-pdf", type = "character", default = NULL),
    make_option("--width", type = "double", default = 6),
    make_option("--height", type = "double", default = 8),
    make_option("--low", type = "character", default = "dodgerblue3"),
    make_option("--mid", type = "character", default = "white"),
    make_option("--high", type = "character", default = "firebrick3")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "chi_squared_from_csv.R --input-csv FILE [options]")
  args <- parse_args(parser)

  if (is.null(args$`input-csv`) || !nzchar(args$`input-csv`)) stop("--input-csv is required")
  if (!file.exists(args$`input-csv`)) stop("Input CSV not found")

  output_pdf <- args$`output-pdf`
  if (is.null(output_pdf) || !nzchar(output_pdf)) {
    output_pdf <- sub("\\.csv$", ".pdf", args$`input-csv`)
  }
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

  raw_data <- read.csv(args$`input-csv`, row.names = 1, check.names = FALSE)
  melted <- melt(as.matrix(raw_data), varnames = c("Row", "Column"), value.name = "Value")

  p <- ggplot(melted, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = args$low, mid = args$mid, high = args$high, midpoint = 0, na.value = "grey90") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_blank(), panel.grid = element_blank())

  ggsave(output_pdf, p, width = args$width, height = args$height)
}

if (sys.nframe() == 0L) {
  main()
}
