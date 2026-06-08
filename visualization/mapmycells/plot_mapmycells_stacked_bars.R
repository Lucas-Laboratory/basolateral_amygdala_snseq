#!/usr/bin/env Rscript

# Generate percentage stacked bar plots for MapMyCells classifications.
#
# Inputs
#   --counts-csv: CSV with columns for cluster, count, and class (required).
#   --colors-csv: CSV mapping unique class names to hex colors (required).
#   --output-pdf: destination PDF (required).
#   --cluster-column: column in counts CSV for clusters (default: `cluster`).
#   --count-column: column for counts (default: `count`).
#   --class-column: column for class labels (default: `supertype_name`).
#   --color-class-column: column in colors CSV for class names (default: `unique_supertype_name`).
#   --color-hex-column: column in colors CSV for hex colours (default: `hex_color`).
#   --width/--height: output size in inches (default: 8 x 6).
#
# Dependencies: optparse, readr, dplyr, ggplot2, scales, rlang

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(rlang)
})

main <- function() {
  option_list <- list(
    make_option("--counts-csv", type = "character", help = "Counts CSV"),
    make_option("--colors-csv", type = "character", help = "Colour mapping CSV"),
    make_option("--output-pdf", type = "character", help = "Output PDF"),
    make_option("--cluster-column", type = "character", default = "cluster"),
    make_option("--count-column", type = "character", default = "count"),
    make_option("--class-column", type = "character", default = "supertype_name"),
    make_option("--color-class-column", type = "character", default = "unique_supertype_name"),
    make_option("--color-hex-column", type = "character", default = "hex_color"),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_mapmycells_stacked_bars.R --counts-csv FILE --colors-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`counts-csv`, args$`colors-csv`, args$`output-pdf`)))) stop("All inputs required")
  if (!file.exists(args$`counts-csv`)) stop("Counts CSV not found")
  if (!file.exists(args$`colors-csv`)) stop("Colour CSV not found")
  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)

  df_counts <- read_csv(args$`counts-csv`, show_col_types = FALSE)
  df_colors <- read_csv(args$`colors-csv`, show_col_types = FALSE)
  names(df_counts) <- sub("^\ufeff", "", names(df_counts))
  names(df_colors) <- sub("^\ufeff", "", names(df_colors))

  df_counts <- df_counts %>%
    group_by(.data[[args$`cluster-column`]]) %>%
    mutate(percentage = 100 * .data[[args$`count-column`]] / sum(.data[[args$`count-column`]])) %>%
    ungroup()

  df_colors_join <- df_colors %>%
    select(all_of(c(args$`color-class-column`, args$`color-hex-column`))) %>%
    distinct()
  df_counts <- df_counts %>%
    left_join(df_colors_join, by = setNames(args$`color-class-column`, args$`class-column`))

  color_mapping <- setNames(df_colors[[args$`color-hex-column`]], df_colors[[args$`color-class-column`]])

  p <- ggplot(df_counts, aes(x = .data[[args$`cluster-column`]], y = percentage,
                             fill = .data[[args$`class-column`]])) +
    geom_bar(stat = "identity") +
    labs(x = "Cluster", y = "Percentage (%)", fill = "Class") +
    theme_minimal() +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_fill_manual(values = color_mapping)

  ggsave(args$`output-pdf`, plot = p, width = args$width, height = args$height)
  message("Stacked bar plot saved to ", args$`output-pdf`)
}

if (sys.nframe() == 0L) {
  main()
}
