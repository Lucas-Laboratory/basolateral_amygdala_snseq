#!/usr/bin/env Rscript

# Colour UMAP embeddings using a cluster-to-colour lookup table.
#
# Inputs
#   --umap-csv: CSV with UMAP coordinates and a cluster column (required).
#   --color-csv: CSV mapping cluster identifiers to hex colours (required).
#   --output-pdf: destination PDF path (required).
#   --x-column / --y-column: coordinate column names (default: `umap_1` / `umap_2`).
#   --cluster-column: column containing cluster identifiers (default: `cluster`).
#   --cluster-map-column: column in colour CSV matching cluster ids (default: `list`).
#   --color-column: column in colour CSV with hex values (default: `color`).
#   --point-size: point size for scatter plot (default: 0.5).
#   --legend: include legend (default: FALSE).
#   --width/--height: PDF size (default: 5 x 5).
#
# Dependencies: optparse, ggplot2, data.table, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

main <- function() {
  option_list <- list(
    make_option("--umap-csv", type = "character", help = "UMAP coordinate CSV"),
    make_option("--color-csv", type = "character", help = "Cluster colour CSV"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--x-column", type = "character", default = "umap_1"),
    make_option("--y-column", type = "character", default = "umap_2"),
    make_option("--cluster-column", type = "character", default = "cluster"),
    make_option("--cluster-map-column", type = "character", default = "list"),
    make_option("--color-column", type = "character", default = "color"),
    make_option("--point-size", type = "double", default = 0.5),
    make_option("--legend", action = "store_true", default = FALSE),
    make_option("--width", type = "double", default = 5),
    make_option("--height", type = "double", default = 5)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_cluster_colored_umap.R --umap-csv FILE --color-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`umap-csv`, args$`color-csv`, args$`output-pdf`)))) {
    stop("--umap-csv, --color-csv, and --output-pdf are required")
  }
  if (!file.exists(args$`umap-csv`)) stop("UMAP CSV not found")
  if (!file.exists(args$`color-csv`)) stop("Color CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  umap_df <- fread(args$`umap-csv`)
  color_df <- fread(args$`color-csv`)
  names(umap_df) <- sub("^\ufeff", "", names(umap_df))
  names(color_df) <- sub("^\ufeff", "", names(color_df))

  if (!args$`x-column` %in% names(umap_df) && "UMAP_1" %in% names(umap_df)) args$`x-column` <- "UMAP_1"
  if (!args$`y-column` %in% names(umap_df) && "UMAP_2" %in% names(umap_df)) args$`y-column` <- "UMAP_2"
  if (!args$`cluster-column` %in% names(umap_df) && "Cluster" %in% names(umap_df)) args$`cluster-column` <- "Cluster"
  if (!args$`cluster-map-column` %in% names(color_df) && "cluster" %in% names(color_df)) args$`cluster-map-column` <- "cluster"
  if (!args$`color-column` %in% names(color_df) && "custom_hex" %in% names(color_df)) args$`color-column` <- "custom_hex"

  required_umap <- c(args$`x-column`, args$`y-column`, args$`cluster-column`)
  if (!all(required_umap %in% names(umap_df))) {
    stop("UMAP CSV missing columns: ", paste(setdiff(required_umap, names(umap_df)), collapse = ", "))
  }
  required_color <- c(args$`cluster-map-column`, args$`color-column`)
  if (!all(required_color %in% names(color_df))) {
    stop("Color CSV missing columns: ", paste(setdiff(required_color, names(color_df)), collapse = ", "))
  }

  umap_df[[args$`cluster-column`]] <- as.character(umap_df[[args$`cluster-column`]])
  color_df[[args$`cluster-map-column`]] <- as.character(color_df[[args$`cluster-map-column`]])

  color_map <- setNames(color_df[[args$`color-column`]], color_df[[args$`cluster-map-column`]])
  umap_df$color <- ifelse(umap_df[[args$`cluster-column`]] %in% names(color_map),
                          color_map[umap_df[[args$`cluster-column`]]], "#CCCCCC")

  if (args$legend) {
    umap_df$cluster_factor <- factor(umap_df[[args$`cluster-column`]], levels = names(color_map))
    p <- ggplot(umap_df, aes_string(x = args$`x-column`, y = args$`y-column`, colour = "cluster_factor")) +
      geom_point(size = args$`point-size`) +
      scale_color_manual(values = color_map, na.value = "#CCCCCC", name = "Cluster")
  } else {
    p <- ggplot(umap_df, aes_string(x = args$`x-column`, y = args$`y-column`, colour = "color")) +
      geom_point(size = args$`point-size`) +
      scale_color_identity()
  }

  p <- p + theme_minimal() + labs(title = "UMAP Projection", x = args$`x-column`, y = args$`y-column`)

  ggsave(args$`output-pdf`, p, width = args$width, height = args$height)
  message("UMAP plot saved to ", args$`output-pdf`)
}

if (sys.nframe() == 0L) {
  main()
}
