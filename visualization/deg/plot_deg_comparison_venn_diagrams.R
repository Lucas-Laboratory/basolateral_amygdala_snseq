#!/usr/bin/env Rscript

# Create scaled two-set Venn diagrams comparing DEG constituents across
# comparisons for each shared cluster, after filtering by adjusted p-value,
# percentage difference, and log2 fold-change thresholds. Outputs per-comparison
# filtered tables, count summaries, and PDFs with positive/negative overlaps.
#
# Inputs
#   --input-dir: directory containing DEG CSVs with columns `gene`, `cluster`,
#                `comparison`, `avg_log2FC`, `pct_diff`, `p_val_adj` (required).
#   --output-dir: where filtered tables and PDFs will be saved (required).
#   --max-p-adj: maximum adjusted p-value to keep (default: 0.01).
#   --min-abs-pct-diff: minimum |pct_diff| required (default: 0.1).
#   --min-abs-logfc: minimum |avg_log2FC| required (default: 0.5).
#   --clusters: optional comma-separated list/ranges limiting clusters processed.
#   --width/--height: PDF dimensions for Venn diagrams (default: 8 x 8).
#
# Dependencies: optparse, readr, dplyr, stringr, grid

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(grid)
})

format_param <- function(x) {
  val <- format(x, scientific = FALSE, trim = TRUE)
  gsub(" ", "", val, fixed = TRUE)
}

parse_clusters <- function(spec) {
  if (is.null(spec) || !nzchar(spec)) return(NULL)
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

required_cols <- c("gene", "cluster", "comparison", "avg_log2FC", "pct_diff", "p_val_adj")

overlap_area <- function(d, r1, r2) {
  if (d >= r1 + r2) {
    0
  } else if (d <= abs(r1 - r2)) {
    pi * min(r1, r2)^2
  } else {
    part1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
    part2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
    part3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    part1 + part2 - part3
  }
}

draw_scaled_venn <- function(set1, set2, comp_names, cluster_num, sign_type) {
  n1 <- length(set1)
  n2 <- length(set2)
  n12 <- length(intersect(set1, set2))

  r1 <- sqrt(n1 / pi)
  r2 <- sqrt(n2 / pi)

  if (n12 <= 0) {
    d <- r1 + r2 + 0.1
  } else if (n12 >= min(n1, n2)) {
    d <- abs(r1 - r2) + 0.1
  } else {
    f <- function(dist) overlap_area(dist, r1, r2) - n12
    d <- uniroot(f, lower = abs(r1 - r2), upper = r1 + r2)$root
  }

  x_min <- min(-r1, d - r2)
  x_max <- max(r1, d + r2)
  y_min <- -max(r1, r2)
  y_max <- max(r1, r2)

  x_margin <- 0.1 * (x_max - x_min)
  y_margin <- 0.1 * (y_max - y_min)
  x_min <- x_min - x_margin
  x_max <- x_max + x_margin
  y_min <- y_min - y_margin
  y_max <- y_max + y_margin

  trans_x <- function(x) (x - x_min) / (x_max - x_min)
  trans_y <- function(y) (y - y_min) / (y_max - y_min)

  center1 <- c(trans_x(0), trans_y(0))
  center2 <- c(trans_x(d), trans_y(0))
  r1_npc <- r1 / (x_max - x_min)
  r2_npc <- r2 / (x_max - x_min)

  grid.newpage()
  grid.circle(x = center1[1], y = center1[2], r = r1_npc,
              gp = gpar(fill = rgb(1, 0, 0, 0.5), col = "black", lwd = 2, fontfamily = "Helvetica"))
  grid.circle(x = center2[1], y = center2[2], r = r2_npc,
              gp = gpar(fill = rgb(0, 0, 1, 0.5), col = "black", lwd = 2, fontfamily = "Helvetica"))

  grid.text(label = n1, x = center1[1], y = center1[2],
            gp = gpar(fontsize = 14, fontfamily = "Helvetica"))
  grid.text(label = n2, x = center2[1], y = center2[2],
            gp = gpar(fontsize = 14, fontfamily = "Helvetica"))

  if (n12 > 0) {
    a <- (r1^2 - r2^2 + d^2) / (2 * d)
    inter_x <- trans_x(a)
    inter_y <- trans_y(0)
    grid.text(label = n12, x = inter_x, y = inter_y,
              gp = gpar(fontsize = 14, fontfamily = "Helvetica", col = "darkgreen"))
  }

  title_str <- sprintf("%s FC genes — cluster %s\n%s vs %s", sign_type, cluster_num, comp_names[1], comp_names[2])
  grid.text(label = title_str, x = 0.5, y = 0.95,
            gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Helvetica"))
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Input directory of DEG CSVs"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--max-p-adj", type = "double", default = 0.01),
    make_option("--min-abs-pct-diff", type = "double", default = 0.1),
    make_option("--min-abs-logfc", type = "double", default = 0.5),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 8)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_deg_comparison_venn_diagrams.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) stop("--input-dir and --output-dir are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  cluster_filter <- parse_clusters(args$clusters)

  files <- list.files(args$`input-dir`, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No CSV files found in input directory")

  filtered_dfs <- list()
  comparison_labels <- character()

  param_suffix <- paste0(
    "mpva-", format_param(args$`max-p-adj`),
    "_mapd-", format_param(args$`min-abs-pct-diff`),
    "_mafc-", format_param(args$`min-abs-logfc`)
  )

  for (file in files) {
    df <- read_csv(file, show_col_types = FALSE)
    missing <- setdiff(required_cols, names(df))
    if (length(missing)) {
      stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
    }

    df <- df %>%
      mutate(cluster = as.integer(cluster),
             abs_pct_diff = abs(pct_diff),
             abs_avg_log2FC = abs(avg_log2FC))

    if (!is.null(cluster_filter)) {
      df <- df %>% filter(cluster %in% cluster_filter)
    }

    if (!nrow(df)) {
      message("No rows remain in ", basename(file), " after filtering clusters; skipping.")
      next
    }

    clusters <- sort(unique(df$cluster))
    comparison_name <- df$comparison[1]
    if (!nzchar(comparison_name)) {
      comparison_name <- tools::file_path_sans_ext(basename(file))
    }
    comparison_labels <- c(comparison_labels, comparison_name)

    all_filtered <- tibble()
    counts <- tibble()

    for (clust in clusters) {
      df_clust <- df %>% filter(cluster == clust)
      df_filtered <- df_clust %>%
        filter(p_val_adj < args$`max-p-adj`,
               abs_pct_diff > args$`min-abs-pct-diff`,
               abs_avg_log2FC > args$`min-abs-logfc`)

      counts <- bind_rows(counts, tibble(comparison = comparison_name,
                                         cluster = clust,
                                         DEG_count = nrow(df_filtered)))
      all_filtered <- bind_rows(all_filtered, df_filtered)
    }

    write_csv(counts, file.path(args$`output-dir`, paste0(comparison_name, "_DEG_counts_", param_suffix, ".csv")))
    write_csv(all_filtered, file.path(args$`output-dir`, paste0(comparison_name, "_DEG_subset_", param_suffix, ".csv")))

    filtered_dfs[[comparison_name]] <- all_filtered
  }

  if (!length(filtered_dfs)) {
    stop("No DEG tables processed; check filters and inputs")
  }

  comparison_first_terms <- sapply(comparison_labels, function(x) str_split(x, " vs ")[[1]][1])
  comparison_sets <- split(comparison_labels, comparison_first_terms)

  for (first_term in names(comparison_sets)) {
    set_comps <- comparison_sets[[first_term]]
    if (length(set_comps) != 2) next

    comp1 <- set_comps[1]
    comp2 <- set_comps[2]

    if (is.null(filtered_dfs[[comp1]]) || is.null(filtered_dfs[[comp2]])) next

    shared_clusters <- intersect(unique(filtered_dfs[[comp1]]$cluster),
                                 unique(filtered_dfs[[comp2]]$cluster))
    shared_clusters <- sort(shared_clusters)
    if (!length(shared_clusters)) next

    pdf_name <- file.path(args$`output-dir`, paste0(first_term, "_VennOverlap_", param_suffix, ".pdf"))
    pdf(file = pdf_name, width = args$width, height = args$height)
    tryCatch({
      for (clust in shared_clusters) {
        comp1_sub <- filtered_dfs[[comp1]] %>% filter(cluster == clust)
        comp2_sub <- filtered_dfs[[comp2]] %>% filter(cluster == clust)

        pos_genes1 <- comp1_sub %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
        pos_genes2 <- comp2_sub %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
        neg_genes1 <- comp1_sub %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()
        neg_genes2 <- comp2_sub %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()

        draw_scaled_venn(pos_genes1, pos_genes2, c(comp1, comp2), clust, "positive")
        draw_scaled_venn(neg_genes1, neg_genes2, c(comp1, comp2), clust, "negative")
      }
    }, finally = dev.off())
    message("Wrote Venn diagrams for ", first_term, " (", comp1, " vs ", comp2, ") to ", pdf_name)
  }
}

if (sys.nframe() == 0L) {
  main()
}
