#!/usr/bin/env Rscript

# Generate per-slice MERFISH scatter plots coloured by custom cluster hex palette.
#
# Inputs
#   --tangram-csv: expression matrix with predicted cluster ids (required).
#   --metadata-csv: cell metadata CSV containing CCF coordinates (required).
#   --hex-csv: CSV mapping clusters to hex colours (required).
#   --output-dir: base directory for outputs (required).
#   --x-axis: coordinate column for x-axis (default: `x`).
#   --y-axis: coordinate column for y-axis (default: `y`).
#   --invert-x/--invert-y: invert axes (default: FALSE/TRUE).
#   --dot-size: point size for PNG scatter plots (default: 5).
#   --pdf-dot-size: point size for PDF scatter plots (default: 0.1).
#   --hide-unassigned: omit cells without a predicted cluster/hex colour (default: keep them as black anatomy).
#   --png-size: pixel size for PNG outputs (default: 3000).
#
# Output structure:
#   <output-dir>/PDF, /PNG, and /Legend directories.
#
# Dependencies: optparse, readr, dplyr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

main <- function() {
  option_list <- list(
    make_option("--tangram-csv", type = "character", help = "Tangram prediction CSV"),
    make_option("--metadata-csv", type = "character", help = "Cell metadata CSV"),
    make_option("--hex-csv", type = "character", help = "Cluster hex colour CSV"),
    make_option("--output-dir", type = "character", help = "Output directory root"),
    make_option("--x-axis", type = "character", default = "x"),
    make_option("--y-axis", type = "character", default = "y"),
    make_option("--invert-x", action = "store_true", default = FALSE),
    make_option("--invert-y", action = "store_true", default = TRUE),
    make_option("--dot-size", type = "double", default = 5),
    make_option("--pdf-dot-size", type = "double", default = 0.1),
    make_option("--hide-unassigned", action = "store_true", default = FALSE,
                help = "Hide black anatomy/background cells without predicted cluster colours"),
    make_option("--show-unassigned", action = "store_true", default = FALSE,
                help = "Deprecated compatibility flag; anatomy/background cells are shown by default"),
    make_option("--png-size", type = "integer", default = 3000)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "custom_hex_merfish_slice_plots.R --tangram-csv FILE --metadata-csv FILE --hex-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  required <- c(args$`tangram-csv`, args$`metadata-csv`, args$`hex-csv`, args$`output-dir`)
  if (any(!nzchar(required))) stop("All input paths and --output-dir are required")
  if (!file.exists(args$`tangram-csv`)) stop("Tangram CSV not found")
  if (!file.exists(args$`metadata-csv`)) stop("Metadata CSV not found")
  if (!file.exists(args$`hex-csv`)) stop("Hex colour CSV not found")

  pdf_dir <- file.path(args$`output-dir`, "PDF")
  png_dir <- file.path(args$`output-dir`, "PNG")
  legend_dir <- file.path(args$`output-dir`, "Legend")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

  tangram <- read_csv(args$`tangram-csv`, show_col_types = FALSE)
  metadata <- read_csv(args$`metadata-csv`, show_col_types = FALSE)
  hexcol <- read_csv(args$`hex-csv`, show_col_types = FALSE)

  names(metadata)[names(metadata) == "cell_label"] <- "cell_id"
  df <- metadata %>%
    left_join(tangram %>% select(cell_id, cluster_id), by = "cell_id") %>%
    left_join(hexcol %>% rename(cluster_id = cluster), by = "cluster_id") %>%
    mutate(is_assigned = !is.na(cluster_id) & !is.na(custom_hex),
           color = ifelse(is_assigned, custom_hex, "#000000"),
           slice = sub(".*\\.", "", brain_section_label))

  show_unassigned <- !args$`hide-unassigned` || args$`show-unassigned`
  if (!show_unassigned) {
    dropped <- sum(!df$is_assigned)
    message("Dropping unassigned/missing-colour cells from plots because --hide-unassigned was set: ", dropped)
    df <- filter(df, is_assigned)
  }
  if (!nrow(df)) {
    stop("No cells available to plot after filtering")
  }

  if (!all(c(args$`x-axis`, args$`y-axis`) %in% names(df))) {
    stop("Specified axis columns not found in metadata")
  }

  df <- df %>%
    mutate(plot_x = (if (args$`invert-x`) -1 else 1) * .data[[args$`x-axis`]],
           plot_y = (if (args$`invert-y`) -1 else 1) * .data[[args$`y-axis`]])

  expand_range <- function(values) {
    rng <- range(values, na.rm = TRUE)
    rng + c(-0.2, 0.2)
  }

  xlim <- expand_range(df$plot_x)
  ylim <- expand_range(df$plot_y)

  base_plot <- function(data, dot_size) {
    unassigned <- filter(data, !is_assigned)
    assigned <- filter(data, is_assigned)

    ggplot() +
      geom_point(data = unassigned, aes(x = plot_x, y = plot_y),
                 colour = "#000000", size = dot_size) +
      geom_point(data = assigned, aes(x = plot_x, y = plot_y, colour = color),
                 size = dot_size) +
      scale_color_identity() +
      coord_fixed(xlim = xlim, ylim = ylim) +
      labs(x = paste0(toupper(args$`x-axis`), " (mm)"),
           y = paste0(toupper(args$`y-axis`), " (mm)")) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(), legend.position = "none")
  }

  for (sl in unique(df$slice)) {
    subdat <- filter(df, slice == sl)
    p_pdf <- base_plot(subdat, args$`pdf-dot-size`) + labs(title = paste("Slice", sl))
    p_png <- base_plot(subdat, args$`dot-size`) + labs(title = paste("Slice", sl))

    pdf(file.path(pdf_dir, paste0("slice_", sl, ".pdf")), width = 6, height = 6)
    print(p_pdf)
    dev.off()

    png(file.path(png_dir, paste0("slice_", sl, ".png")), width = args$`png-size`, height = args$`png-size`, bg = "transparent")
    print(p_png + theme(
      plot.title = element_text(color = "transparent"),
      axis.title = element_text(color = "transparent"),
      axis.text = element_text(color = "transparent"),
      axis.ticks = element_line(color = "transparent"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ))
    dev.off()
  }

  legend_df <- hexcol
  p_leg <- ggplot(legend_df, aes(x = factor(cluster), y = 1, fill = factor(cluster))) +
    geom_tile() +
    scale_fill_manual(name = "Cluster", values = setNames(legend_df$custom_hex, legend_df$cluster),
                      guide = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 10))

  pdf(file.path(legend_dir, "cluster_legend.pdf"), width = 24, height = 2)
  print(p_leg)
  dev.off()

  png(file.path(legend_dir, "cluster_legend.png"), width = 2500, height = 200, bg = "transparent")
  print(p_leg)
  dev.off()
}

if (sys.nframe() == 0L) {
  main()
}
