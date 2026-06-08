#!/usr/bin/env Rscript

# Compare MapMyCells heterotypic doublet estimates with DoubletFinder scores by
# plotting ranked probabilities, binned DoubletFinder averages, and doublet rug.
#
# Inputs
#   --input-csv: CSV containing columns `heterotypic_doublet_estimated_probability`,
#                `doublet_score`, and `doublet_class` (required).
#   --output-dir: directory to save the PDF plot (required).
#   --output-prefix: optional prefix for output filenames.
#   --bin-size: number of ranked cells per DoubletFinder average bin (default: 100).
#   --width/--height: PDF size in inches (default: 8 x 5).
#
# Output
#   - PDF combining LISI-like heterotypic doublet probability curve with
#     DoubletFinder histogram-style bars and a rug showing called doublets.
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
    make_option("--input-csv", type = "character", help = "Merged doublet metrics CSV"),
    make_option("--output-dir", type = "character", help = "Directory for plot output"),
    make_option("--output-prefix", type = "character", default = "mapmycells_doublet"),
    make_option("--bin-size", type = "integer", default = 100,
                help = "Rank bin size for DoubletFinder averages"),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 5)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "mapmycells_doublet_comparison_plot.R --input-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-csv`, args$`output-dir`)))) {
    stop("--input-csv and --output-dir are required")
  }
  if (!file.exists(args$`input-csv`)) stop("Input CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  df <- read_csv(args$`input-csv`, show_col_types = FALSE)
  required <- c("heterotypic_doublet_estimated_probability", "doublet_score", "doublet_class")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Input CSV missing columns: ", paste(missing, collapse = ", "))
  }

  df_plot <- df %>%
    transmute(rank = row_number(),
              heterotypic_prob = heterotypic_doublet_estimated_probability,
              doublet_score = doublet_score,
              doublet_class = as.character(doublet_class))

  bin_size <- args$`bin-size`
  if (bin_size <= 0) stop("--bin-size must be positive")

  df_hist <- df_plot %>%
    mutate(bin_start = ((rank - 1) %/% bin_size) * bin_size + 1,
           bin_mid = bin_start + bin_size / 2) %>%
    group_by(bin_mid) %>%
    summarise(avg_doublet_score = mean(doublet_score, na.rm = TRUE), .groups = "drop")

  output_pdf <- file.path(args$`output-dir`, paste0(args$`output-prefix`, "_comparison.pdf"))

  p <- ggplot(df_plot, aes(x = rank)) +
    geom_col(data = df_hist, aes(x = bin_mid, y = avg_doublet_score), fill = "grey80",
             width = bin_size * 0.9, inherit.aes = FALSE) +
    geom_line(aes(y = heterotypic_prob), colour = "firebrick4", linewidth = 0.6) +
    geom_rug(data = df_plot %>% filter(tolower(doublet_class) == "doublet"),
             sides = "b", alpha = 0.4, colour = "dodgerblue4") +
    scale_y_continuous(name = "Heterotypic doublet estimate",
                       sec.axis = sec_axis(~ ., name = "Avg DoubletFinder score")) +
    scale_x_continuous(name = "Cell rank (sorted input order)") +
    theme_minimal(base_size = 11) +
    labs(title = "MapMyCells vs DoubletFinder doublet comparison",
         caption = sprintf("Bin size: %s", bin_size))

  ggsave(output_pdf, p, width = args$width, height = args$height)
  message("Saved plot to ", output_pdf)
}

if (sys.nframe() == 0L) {
  main()
}
