#!/usr/bin/env Rscript

# Generate wide-format GSEA summary tables plus clustered heatmaps for NES (or
# another metric) across clusters.
#
# Inputs
#   --input-dir: directory containing unique_GSEA CSVs and comparison subfolders (required).
#   --output-dir: where CSVs and PDFs will be written (required).
#   --metric-col: column to visualise (default: NES).
#   --padj-cutoff: adjusted p-value threshold used to keep terms (default: 1e-4).
#   --clusters: comma-separated list and/or ranges (e.g. 0:26,28,29).
#   --font-size-x/--font-size-y: axis label sizes (default: 8 / 6).
#   --width/--height: heatmap page size in inches (default: 8 x 12).
#   --unique-pattern: regex to match unique term CSVs (default matches BP/CC/MF DvM/PvM/PvD files).
#   --full-template: sprintf template for full GSEA CSV within comparison subdir (default: GSEA_%s_reint_full_mod_DEG_%s.csv).
#
# Dependencies: optparse, dplyr, readr, tidyr, stringr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

parse_clusters <- function(spec) {
  pieces <- trimws(strsplit(spec, ",")[[1]])
  clusters <- integer(0)
  for (piece in pieces) {
    if (!nzchar(piece)) next
    if (grepl(":", piece, fixed = TRUE)) {
      bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1]])
      if (length(bounds) != 2 || any(is.na(bounds))) stop("Invalid cluster range: ", piece)
      clusters <- c(clusters, seq(bounds[1], bounds[2]))
    } else {
      val <- suppressWarnings(as.integer(piece))
      if (is.na(val)) stop("Invalid cluster value: ", piece)
      clusters <- c(clusters, val)
    }
  }
  unique(clusters)
}

format_param <- function(x) {
  val <- formatC(x, format = "fg", digits = 10, flag = "#")
  gsub(" ", "", val, fixed = TRUE)
}

build_paths <- function(input_dir, comparison, ontology, template) {
  csv_name <- sprintf(template, ontology, comparison)
  file.path(input_dir, comparison, csv_name)
}

required_columns <- c("ID", "Description", "cluster")

read_full_gsea <- function(path, metric_col, clusters) {
  df <- read_csv(path, show_col_types = FALSE)
  missing <- setdiff(c(required_columns, metric_col, "p.adjust"), names(df))
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  df %>%
    mutate(cluster = as.integer(cluster)) %>%
    filter(cluster %in% clusters)
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Base directory for GSEA results"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--metric-col", type = "character", default = "NES"),
    make_option("--padj-cutoff", type = "double", default = 1e-4),
    make_option("--clusters", type = "character", default = "0:26,28,29"),
    make_option("--font-size-x", type = "double", default = 8),
    make_option("--font-size-y", type = "double", default = 6),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 12),
    make_option("--unique-pattern", type = "character",
                default = "^unique_GSEA_(BP|CC|MF)_reint_full_mod_DEG_(DvM|PvM|PvD)\\.csv$"),
    make_option("--full-template", type = "character",
                default = "GSEA_%s_reint_full_mod_DEG_%s.csv")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_gsea_cluster_heatmap.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) stop("--input-dir and --output-dir are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  clusters <- parse_clusters(args$clusters)
  if (!length(clusters)) stop("No clusters parsed from --clusters")
  cluster_cols <- as.character(clusters)

  unique_files <- list.files(args$`input-dir`, pattern = args$`unique-pattern`, full.names = TRUE)
  if (!length(unique_files)) stop("No files matched --unique-pattern in input directory")

  for (unique_fp in unique_files) {
    matches <- str_match(basename(unique_fp), args$`unique-pattern`)
    if (all(is.na(matches))) {
      message("Skipping unmatched file: ", unique_fp)
      next
    }
    ontology <- matches[2]
    comparison <- matches[3]

    unique_df_all <- read_csv(unique_fp, show_col_types = FALSE)
    if (!all(c("ID", "Description") %in% names(unique_df_all))) {
      stop("unique CSV missing ID/Description columns: ", unique_fp)
    }

    gsea_fp <- build_paths(args$`input-dir`, comparison, ontology, args$`full-template`)
    if (!file.exists(gsea_fp)) {
      stop("Full GSEA file not found: ", gsea_fp)
    }

    full_gsea_df <- read_full_gsea(gsea_fp, args$`metric-col`, clusters)
    sig_terms <- full_gsea_df %>%
      filter(p.adjust < args$`padj-cutoff`) %>%
      distinct(ID, Description)
    if (!nrow(sig_terms)) {
      message("No significant terms for ", ontology, " ", comparison, " below padj cutoff; skipping.")
      next
    }

    unique_df <- unique_df_all %>% semi_join(sig_terms, by = c("ID", "Description"))

    dir_summ <- full_gsea_df %>%
      filter(ID %in% sig_terms$ID) %>%
      transmute(ID, Description, cluster = as.character(cluster),
                direction = case_when(
                  .data[[args$`metric-col`]] > 0 ~ "UP",
                  .data[[args$`metric-col`]] < 0 ~ "DOWN",
                  TRUE ~ NA_character_
                )) %>%
      distinct(ID, Description, cluster, direction)

    dir_wide <- dir_summ %>%
      pivot_wider(id_cols = c(ID, Description), names_from = cluster, values_from = direction) %>%
      right_join(unique_df, by = c("ID", "Description")) %>%
      relocate(ID, Description)

    missing_cols <- setdiff(cluster_cols, names(dir_wide))
    for (mc in missing_cols) dir_wide[[mc]] <- NA_character_
    dir_wide <- dir_wide %>% select(ID, Description, all_of(cluster_cols))

    padj_label <- paste0("padj", gsub("\\\\.", "_", format_param(args$`padj-cutoff`)))
    cluster_label <- paste0(min(clusters), "-", max(clusters))
    metric_label <- args$`metric-col`

    dir_fname <- paste("clusters_direction_GSEA", metric_label, padj_label,
                      paste0("cl", cluster_label), ontology, comparison, sep = "_")
    write_csv(dir_wide, file.path(args$`output-dir`, paste0(dir_fname, ".csv")))

    metric_summ <- full_gsea_df %>%
      filter(ID %in% sig_terms$ID) %>%
      transmute(ID, Description, cluster = as.character(cluster),
                value = .data[[args$`metric-col`]]) %>%
      distinct(ID, Description, cluster, value)

    template <- unique_df %>% tidyr::crossing(cluster = cluster_cols)
    full_metric <- template %>%
      left_join(metric_summ, by = c("ID", "Description", "cluster")) %>%
      mutate(metric_final = if_else(is.na(value), 0, value))

    metric_wide <- full_metric %>%
      select(ID, Description, cluster, metric_final) %>%
      pivot_wider(id_cols = c(ID, Description), names_from = cluster, values_from = metric_final)

    met_fname <- paste("clusters", tolower(metric_label), "GSEA", padj_label,
                       paste0("cl", cluster_label), ontology, comparison, sep = "_")
    write_csv(metric_wide, file.path(args$`output-dir`, paste0(met_fname, ".csv")))

    mat_df <- metric_wide %>% select(all_of(cluster_cols))
    mat_df[is.na(mat_df)] <- 0
    mat <- as.matrix(mat_df)
    rownames(mat) <- metric_wide$Description

    if (nrow(mat) >= 2) {
      dist_rows <- dist(mat)
      hc_rows <- hclust(dist_rows, method = "average")
      term_order <- hc_rows$labels[hc_rows$order]
    } else {
      term_order <- metric_wide$Description
    }

    heat_df <- full_metric %>%
      mutate(cluster = as.integer(cluster),
             Description = factor(Description, levels = term_order))

    heat_fname <- paste("heatmap", tolower(metric_label), "GSEA", padj_label,
                        paste0("cl", cluster_label), ontology, comparison, "clustered", sep = "_")

    p <- ggplot(heat_df, aes(x = cluster, y = Description, fill = metric_final)) +
      geom_tile() +
      scale_x_continuous(breaks = clusters, labels = clusters, expand = c(0, 0)) +
      scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0,
                           na.value = "grey90") +
      labs(x = "Cluster", fill = args$`metric-col`,
           title = paste(ontology, comparison,
                         sprintf("(p.adjust < %s; clusters %s)", args$`padj-cutoff`, cluster_label))) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = args$`font-size-x`, angle = 90, vjust = 0.5),
            axis.text.y = element_text(size = args$`font-size-y`),
            axis.ticks.y = element_blank())

    ggsave(filename = file.path(args$`output-dir`, paste0(heat_fname, ".pdf")),
           plot = p, width = args$width, height = args$height)

    message("Processed ", ontology, " ", comparison, " → ", heat_fname)
  }
}

if (sys.nframe() == 0L) {
  main()
}
