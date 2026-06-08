#!/usr/bin/env Rscript

# Create bar plots summarising enriched/de-enriched GSEA term counts across clusters.
#
# Inputs
#   --comparison-dirs: comma-separated key=path pairs (e.g. `DvM=path1,PvM=path2`) (required).
#   --categories: comma-separated GO categories to process (default: `CC,MF,BP`).
#   --output-dir: directory for PDF outputs (required).
#   --clusters: cluster range or list (default: `0:29`).
#   --width/--height: PDF size (default: 12 x 6).
#
# Each directory should contain CSV files with columns `cluster` and `NES` for each GO category.
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
    if (length(kv) != 2) stop("Invalid --comparison-dirs entry: ", entry)
    kv
  })
  setNames(vapply(comps, `[`, character(1), 2), vapply(comps, `[`, character(1), 1))
}

main <- function() {
  option_list <- list(
    make_option("--comparison-dirs", type = "character", help = "Comma-separated key=dir list"),
    make_option("--categories", type = "character", default = "CC,MF,BP"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--clusters", type = "character", default = "0:29"),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_gsea_term_count_histograms.R --comparison-dirs key=path,... --output-dir DIR [options]")
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
      files <- list.files(comp_dirs[[comp]], pattern = paste0(cat, ".*\\.csv$"), full.names = TRUE)
      counts <- if (!length(files)) {
        tibble(cluster = cluster_levels, enriched_count = 0L, de_enriched_count = 0L)
      } else {
        read_csv(files[1], show_col_types = FALSE) %>%
          mutate(cluster = as.character(cluster), NES = as.numeric(NES)) %>%
          group_by(cluster) %>%
          summarise(enriched_count = sum(NES > 0, na.rm = TRUE),
                    de_enriched_count = sum(NES < 0, na.rm = TRUE), .groups = "drop") %>%
          right_join(tibble(cluster = cluster_levels), by = "cluster") %>%
          replace_na(list(enriched_count = 0L, de_enriched_count = 0L))
      }
      mutate(counts, comparison = comp)
    })
    all_counts <- bind_rows(all_counts)

    plot_data <- all_counts %>%
      pivot_longer(cols = c(enriched_count, de_enriched_count), names_to = "direction", values_to = "count") %>%
      mutate(count = ifelse(direction == "de_enriched_count", -count, count),
             comparison = factor(comparison, levels = names(comp_dirs)),
             cluster = factor(cluster, levels = cluster_levels))

    p <- ggplot(plot_data, aes(x = cluster, y = count, fill = comparison)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_hline(yintercept = 0) +
      labs(title = paste0("GSEA term counts (", cat, ")"), x = "Cluster", y = "Enriched (+) / De-enriched (–) term count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5))

    ggsave(file.path(args$`output-dir`, paste0("GSEA_Counts_", cat, ".pdf")), p, width = args$width, height = args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
