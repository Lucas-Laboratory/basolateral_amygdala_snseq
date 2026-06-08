#!/usr/bin/env Rscript

# Create bar plots summarising ORA term counts across clusters and comparisons.
#
# Inputs
#   --comparison-dirs: comma-separated key=dir strings (e.g. `DvM=path,PvM=path`) (required).
#   --categories: comma-separated GO categories (default: `BP,CC,MF`).
#   --output-dir: directory for output PDFs (required).
#   --clusters: cluster range/list (default: `0:29`).
#   --width/--height: PDF size (default: 12 x 6).
#
# Each directory should contain ORA CSV files with columns `cluster`, `direction`, and optionally others.
#
# Dependencies: optparse, readr, dplyr, tidyr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

parse_comparison_dirs <- function(spec) {
  entries <- trimws(strsplit(spec, ",")[[1]])
  comps <- lapply(entries, function(entry) {
    kv <- trimws(strsplit(entry, "=")[[1]])
    if (length(kv) != 2) stop("Invalid comparison entry: ", entry)
    kv
  })
  setNames(vapply(comps, `[`, character(1), 2), vapply(comps, `[`, character(1), 1))
}

main <- function() {
  option_list <- list(
    make_option("--comparison-dirs", type = "character", help = "Comma-separated key=dir list"),
    make_option("--categories", type = "character", default = "BP,CC,MF"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--clusters", type = "character", default = "0:29"),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_ora_term_count_histograms.R --comparison-dirs key=dir,... --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`comparison-dirs`, args$`output-dir`)))) stop("Required arguments missing")
  comp_dirs <- parse_comparison_dirs(args$`comparison-dirs`)
  if (!all(dir.exists(comp_dirs))) stop("One or more comparison directories not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  categories <- trimws(strsplit(args$categories, ",")[[1]])
  cluster_levels <- if (grepl(":", args$clusters)) {
    bounds <- as.integer(strsplit(args$clusters, ":")[[1]])
    as.character(seq(bounds[1], bounds[2]))
  } else {
    trimws(strsplit(args$clusters, ",")[[1]])
  }

  for (cat in categories) {
    all_counts <- lapply(names(comp_dirs), function(comp) {
      files <- list.files(comp_dirs[[comp]], pattern = paste0("GO_", cat, ".*\\.csv$"), full.names = TRUE)
      if (!length(files)) {
        expand_grid(cluster = cluster_levels, direction = c("positive", "negative")) %>% mutate(count = 0L)
      } else {
        df <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))
        df %>% mutate(cluster = as.character(cluster)) %>%
          count(cluster, direction, name = "count") %>%
          right_join(expand_grid(cluster = cluster_levels, direction = c("positive", "negative")), by = c("cluster", "direction")) %>%
          replace_na(list(count = 0L))
      }
    })
    names(all_counts) <- names(comp_dirs)
    plot_data <- bind_rows(all_counts, .id = "comparison") %>%
      mutate(comparison = factor(comparison, levels = names(comp_dirs)),
             cluster = factor(cluster, levels = cluster_levels),
             count = ifelse(direction == "negative", -count, count))

    p <- ggplot(plot_data, aes(x = cluster, y = count, fill = comparison)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_hline(yintercept = 0) +
      labs(title = paste0("ORA term counts (", cat, ")"), x = "Cluster", y = "Enriched (+) / De-enriched (–) term count") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                         plot.title = element_text(hjust = 0.5))

    ggsave(file.path(args$`output-dir`, paste0("ORA_Counts_", cat, ".pdf")), p, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
