#!/usr/bin/env Rscript

# Build a Venn-style summary of overlapping DEGs across pairwise comparisons.
#
# Inputs
#   --manifest: TSV/CSV with columns `label` and `path` pointing to DEG tables (required).
#               Labels must follow `GroupA vs GroupB` (configurable via --label-separator).
#   --output-csv: destination for the overlap summary (required).
#   --label-separator: separator between group names inside labels (default: " vs ").
#   --positive-threshold: minimum log2FC considered positive (default: >0).
#   --negative-threshold: maximum log2FC considered negative (default: <0).
#
# Output
#   CSV with overlap counts grouped by `grouping` (first/last term), `set_identifier`, `cluster`, and `sign`.
#   Columns include `comp1`, `comp2`, `n1`, `n2`, `unique1`, `unique2`, `overlap`.
#
# Dependencies: optparse, readr, dplyr, stringr, tidyr

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

load_deg_tables <- function(manifest_path) {
  manifest <- readr::read_delim(manifest_path, delim = if (grepl("\\.tsv$", manifest_path, ignore.case = TRUE)) "\t" else ",",
                                show_col_types = FALSE)
  if (!all(c("label", "path") %in% names(manifest))) {
    stop("Manifest must contain columns 'label' and 'path'")
  }

  manifest$path <- trimws(manifest$path)
  manifest$label <- trimws(manifest$label)

  lapply(seq_len(nrow(manifest)), function(i) {
    label <- manifest$label[i]
    path <- manifest$path[i]
    if (!file.exists(path)) stop("File not found for label ", label, ": ", path)
    df <- readr::read_csv(path, show_col_types = FALSE)
    if (!all(c("gene", "cluster", "avg_log2FC") %in% names(df))) {
      stop("DEG file for label ", label, " must contain gene, cluster, avg_log2FC columns")
    }
    df$label <- label
    df
  })
}

extract_pairs <- function(labels, separator) {
  parsed <- str_split(labels, pattern = separator, n = 2)
  if (any(vapply(parsed, length, integer(1)) < 2)) {
    stop("Each label must contain separator '", separator, "'")
  }
  data.frame(
    label = labels,
    first = vapply(parsed, `[`, character(1), 1),
    last = vapply(parsed, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

summarise_overlap <- function(df_list, pairs, separator, pos_threshold, neg_threshold) {
  combined <- bind_rows(df_list)
  stopifnot("label" %in% names(combined))

  result_rows <- list()

  add_rows <- function(grouping, identifier, comp1, comp2, cluster_id, genes1, genes2, sign_label) {
    overlap_genes <- intersect(genes1, genes2)
    result_rows[[length(result_rows) + 1]] <<- data.frame(
      grouping = grouping,
      set_identifier = identifier,
      cluster = cluster_id,
      sign = sign_label,
      comp1 = comp1,
      comp2 = comp2,
      n1 = length(genes1),
      n2 = length(genes2),
      unique1 = length(setdiff(genes1, genes2)),
      unique2 = length(setdiff(genes2, genes1)),
      overlap = length(overlap_genes),
      stringsAsFactors = FALSE
    )
  }

  for (group_type in c("first", "last")) {
    group_col <- if (group_type == "first") "first" else "last"
    split_labels <- split(pairs$label, pairs[[group_col]])
    for (identifier in names(split_labels)) {
      labels_in_group <- split_labels[[identifier]]
      if (length(labels_in_group) != 2) next
      comp1 <- labels_in_group[1]
      comp2 <- labels_in_group[2]

      clusters1 <- combined %>% filter(label == comp1) %>% pull(cluster) %>% unique()
      clusters2 <- combined %>% filter(label == comp2) %>% pull(cluster) %>% unique()
      shared_clusters <- intersect(clusters1, clusters2)

      for (cluster_id in shared_clusters) {
        comp1_df <- combined %>% filter(label == comp1, cluster == cluster_id)
        comp2_df <- combined %>% filter(label == comp2, cluster == cluster_id)

        genes1_pos <- comp1_df %>% filter(avg_log2FC > pos_threshold) %>% pull(gene) %>% unique()
        genes2_pos <- comp2_df %>% filter(avg_log2FC > pos_threshold) %>% pull(gene) %>% unique()
        add_rows(group_type, identifier, comp1, comp2, cluster_id, genes1_pos, genes2_pos, "positive")

        genes1_neg <- comp1_df %>% filter(avg_log2FC < neg_threshold) %>% pull(gene) %>% unique()
        genes2_neg <- comp2_df %>% filter(avg_log2FC < neg_threshold) %>% pull(gene) %>% unique()
        add_rows(group_type, identifier, comp1, comp2, cluster_id, genes1_neg, genes2_neg, "negative")
      }
    }
  }

  bind_rows(result_rows)
}

main <- function() {
  option_list <- list(
    make_option("--manifest", type = "character", help = "Manifest mapping labels to DEG CSVs"),
    make_option("--output-csv", type = "character", help = "Destination summary CSV"),
    make_option("--label-separator", type = "character", default = " vs ",
                help = "Separator used between group names [default %default]"),
    make_option("--positive-threshold", type = "double", default = 0,
                help = "Lower bound for positive log2FC [default > %default]"),
    make_option("--negative-threshold", type = "double", default = 0,
                help = "Upper bound for negative log2FC [default < %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "summarize_deg_overlap_counts.R --manifest FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  if (is.null(args$manifest) || !nzchar(args$manifest)) stop("--manifest is required")
  if (is.null(args$`output-csv`) || !nzchar(args$`output-csv`)) stop("--output-csv is required")
  if (!file.exists(args$manifest)) stop("Manifest not found: ", args$manifest)

  deg_tables <- load_deg_tables(args$manifest)
  labels <- vapply(deg_tables, function(df) unique(df$label), character(1))
  pairs <- extract_pairs(labels, args$`label-separator`)

  summary <- summarise_overlap(
    df_list = deg_tables,
    pairs = pairs,
    separator = args$`label-separator`,
    pos_threshold = args$`positive-threshold`,
    neg_threshold = args$`negative-threshold`
  )

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary, args$`output-csv`)
  message("Venn summary written to ", args$`output-csv`)
}

if (sys.nframe() == 0L) {
  main()
}
