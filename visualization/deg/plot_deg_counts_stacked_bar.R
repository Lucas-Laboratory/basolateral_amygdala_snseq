#!/usr/bin/env Rscript

# Plot stacked DEG counts (unique/overlap) for comparison pairs from a Venn summary CSV.
#
# Inputs
#   --venn-csv: Venn summary CSV containing columns `grouping`, `set_identifier`, `comp1`, `comp2`, `unique1`, `unique2`, `overlap`, `sign`, `cluster` (required).
#   --output-dir: directory to write PDFs (required).
#   --grouping: grouping filter to use (default: `last`).
#   --clusters: optional comma-separated cluster IDs (default: 0-29).
#   --colors-positive: comma-separated colours for unique1, overlap, unique2 (default: maroon4,slategray,turquoise4).
#   --colors-negative: comma-separated colours for unique1, overlap, unique2 on the negative side (default: black,darkslategray,grey40).
#   --width/--height: PDF size (default: 6 x 4).
#
# Dependencies: optparse, ggplot2, dplyr, tidyr, readr

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

assign_colors <- function(df, colors_pos, colors_neg) {
  df %>% mutate(fill_color = case_when(
    sign == "positive" & segment == "unique1" ~ colors_pos[1],
    sign == "positive" & segment == "overlap" ~ colors_pos[2],
    sign == "positive" & segment == "unique2" ~ colors_pos[3],
    sign == "negative" & segment == "unique1" ~ colors_neg[1],
    sign == "negative" & segment == "overlap" ~ colors_neg[2],
    sign == "negative" & segment == "unique2" ~ colors_neg[3],
    TRUE ~ "grey"
  ))
}

main <- function() {
  option_list <- list(
    make_option("--venn-csv", type = "character", help = "Venn summary CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--grouping", type = "character", default = "last"),
    make_option("--clusters", type = "character", default = "0:29"),
    make_option("--colors-positive", type = "character", default = "maroon4,slategray,turquoise4"),
    make_option("--colors-negative", type = "character", default = "black,darkslategray,grey40"),
    make_option("--width", type = "double", default = 6),
    make_option("--height", type = "double", default = 4)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_deg_counts_stacked_bar.R --venn-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`venn-csv`, args$`output-dir`)))) stop("--venn-csv and --output-dir are required")
  if (!file.exists(args$`venn-csv`)) stop("Venn CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  data_all <- read_csv(args$`venn-csv`, show_col_types = FALSE) %>% filter(grouping == args$grouping)
  set_ids <- unique(data_all$set_identifier)
  colors_pos <- trimws(strsplit(args$`colors-positive`, ",")[[1]])
  colors_neg <- trimws(strsplit(args$`colors-negative`, ",")[[1]])
  cluster_levels <- if (grepl(":", args$clusters)) {
    bounds <- as.integer(strsplit(args$clusters, ":")[[1]])
    as.character(seq(bounds[1], bounds[2]))
  } else {
    trimws(strsplit(args$clusters, ",")[[1]])
  }

  for (setid in set_ids) {
    df <- data_all %>% filter(set_identifier == setid)
    df <- df %>% mutate(comp1_clean = sub(" vs .*", "", comp1), comp2_clean = sub(" vs .*", "", comp2))
    df_long <- df %>% mutate(row_id = row_number()) %>%
      pivot_longer(cols = c(unique1, unique2, overlap), names_to = "segment", values_to = "count") %>%
      mutate(count = as.numeric(count), count = ifelse(sign == "negative", -count, count))

    df_long$cluster <- factor(df_long$cluster, levels = cluster_levels)
    df_long$segment <- factor(df_long$segment, levels = c("unique1", "overlap", "unique2"))
    df_long <- assign_colors(df_long, colors_pos, colors_neg)

    p <- ggplot(df_long, aes(x = cluster, y = count, fill = fill_color)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_hline(yintercept = 0, colour = "black") +
      ggtitle(paste("Compared to", setid)) +
      labs(x = "Cluster", y = "Number of DEGs") +
      theme_bw() +
      scale_fill_identity()

    filename <- file.path(args$`output-dir`, paste0("stacked_bar_compared_to_", setid, ".pdf"))
    ggsave(filename, plot = p, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
