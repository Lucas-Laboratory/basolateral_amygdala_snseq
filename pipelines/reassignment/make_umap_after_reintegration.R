#!/usr/bin/env Rscript

# Plot UMAP embeddings coloured by cluster with highlighted reassigned barcodes.
#
# Inputs
#   --umap-csv: barcode-level UMAP coordinates with columns `Barcode`, `UMAP_1`, `UMAP_2`, `new_cluster` (required).
#   --hex-csv: cluster-to-hex colour mapping (columns `cluster`, `custom_hex`) (required).
#   --outline-csv: optional CSV listing barcodes to highlight (column `barcode`).
#   --output-pdf: destination PDF (required).
#   --point-size: base point size (default: 1).
#   --highlight-source: `auto`, `reassigned`, `outline`, or `none` (default: auto).
#   --highlight-point-size: highlighted point size (default: derived from point size).
#   --outline-size: stroke width for highlighted points (default: 0.8).
#   --default-color: fallback hex colour (default: `#000000`).
#
# Dependencies: optparse, readr, dplyr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

normalise_barcode <- function(value) {
  gsub("-target$", "", trimws(as.character(value)))
}

main <- function() {
  option_list <- list(
    make_option("--umap-csv", type = "character", help = "UMAP coordinate CSV"),
    make_option("--hex-csv", type = "character", help = "Cluster colour CSV"),
    make_option("--outline-csv", type = "character", default = NULL,
                help = "Optional barcode outline CSV"),
    make_option("--output-pdf", type = "character", help = "Output PDF path"),
    make_option("--point-size", type = "double", default = 1),
    make_option("--highlight-source", type = "character", default = "auto",
                help = "Highlight source: auto, reassigned, outline, or none"),
    make_option("--highlight-point-size", type = "double", default = NA,
                help = "Point size for highlighted/reassigned cells"),
    make_option("--outline-size", type = "double", default = 0.8),
    make_option("--default-color", type = "character", default = "#000000")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "make_umap_after_reintegration.R --umap-csv FILE --hex-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  required <- c(args$`umap-csv`, args$`hex-csv`, args$`output-pdf`)
  if (any(!nzchar(required))) stop("--umap-csv, --hex-csv, and --output-pdf are required")
  if (!file.exists(args$`umap-csv`)) stop("UMAP CSV not found")
  if (!file.exists(args$`hex-csv`)) stop("Hex colour CSV not found")
  if (!is.null(args$`outline-csv`) && !file.exists(args$`outline-csv`)) stop("Outline CSV not found")
  highlight_source <- tolower(args$`highlight-source`)
  allowed_highlight_sources <- c("auto", "reassigned", "outline", "none")
  if (!highlight_source %in% allowed_highlight_sources) {
    stop("--highlight-source must be one of: ", paste(allowed_highlight_sources, collapse = ", "))
  }

  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  hex_df <- read_csv(args$`hex-csv`, show_col_types = FALSE)
  if (!all(c("cluster", "custom_hex") %in% names(hex_df))) stop("Hex CSV must contain 'cluster' and 'custom_hex'")
  hex_colors <- setNames(hex_df$custom_hex, as.character(hex_df$cluster))

  umap_df <- read_csv(args$`umap-csv`, show_col_types = FALSE)
  if (!all(c("Barcode", "UMAP_1", "UMAP_2", "new_cluster") %in% names(umap_df))) {
    stop("UMAP CSV must contain columns: Barcode, UMAP_1, UMAP_2, new_cluster")
  }

  if (!is.null(args$`outline-csv`)) {
    outline_df <- read_csv(args$`outline-csv`, show_col_types = FALSE)
    if (!"barcode" %in% names(outline_df)) stop("Outline CSV must contain 'barcode'")
    outline_barcodes <- normalise_barcode(outline_df$barcode)
  } else {
    outline_barcodes <- character()
  }

  barcode_key <- normalise_barcode(umap_df$Barcode)
  outline_from_csv <- barcode_key %in% outline_barcodes
  reassigned <- rep(FALSE, nrow(umap_df))
  if ("predicted_cluster" %in% names(umap_df)) {
    reassigned <- reassigned | !is.na(umap_df$predicted_cluster)
  }
  if (all(c("old_cluster", "new_cluster") %in% names(umap_df))) {
    reassigned <- reassigned | (!is.na(umap_df$old_cluster) &
                                  !is.na(umap_df$new_cluster) &
                                  as.character(umap_df$old_cluster) != as.character(umap_df$new_cluster))
  }
  if (highlight_source == "reassigned" && !any(reassigned)) {
    warning("No reassigned/newly predicted cells detected from predicted_cluster or old_cluster/new_cluster columns")
  }

  highlight <- switch(highlight_source,
                      auto = outline_from_csv | reassigned,
                      reassigned = reassigned,
                      outline = outline_from_csv,
                      none = rep(FALSE, nrow(umap_df)))
  highlight_point_size <- args$`highlight-point-size`
  if (is.na(highlight_point_size)) {
    highlight_point_size <- max(args$`point-size` + 0.9, args$`point-size` * 1.8)
  }

  umap_df <- umap_df %>%
    mutate(custom_hex = hex_colors[as.character(new_cluster)],
           custom_hex = if_else(is.na(custom_hex), args$`default-color`, custom_hex),
           highlight = highlight)

  missing_clusters <- setdiff(unique(umap_df$new_cluster), names(hex_colors))
  if (length(missing_clusters)) {
    warning("Missing hex codes for clusters: ", paste(missing_clusters, collapse = ", "))
  }

  base_df <- filter(umap_df, !highlight)
  highlight_df_plot <- filter(umap_df, highlight)
  message("Highlighted cells: ", nrow(highlight_df_plot))

  p <- ggplot() +
    geom_point(data = base_df,
               aes(x = UMAP_1, y = UMAP_2, colour = factor(new_cluster)),
               size = args$`point-size`, alpha = 1, shape = 16) +
    geom_point(data = highlight_df_plot,
               aes(x = UMAP_1, y = UMAP_2, fill = factor(new_cluster)),
               colour = "black", size = highlight_point_size, stroke = args$`outline-size`, shape = 21) +
    scale_color_manual(name = "Cluster", values = hex_colors, na.value = args$`default-color`) +
    scale_fill_manual(values = hex_colors, guide = "none", na.value = args$`default-color`) +
    coord_fixed() +
    theme_classic() +
    labs(x = "UMAP 1", y = "UMAP 2")

  ggsave(args$`output-pdf`, p, width = 6, height = 6)
}

if (sys.nframe() == 0L) {
  main()
}
