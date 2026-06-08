#!/usr/bin/env Rscript

# Build bar plots summarising tangram assignments across Allen Brain Atlas substructures.
#
# Inputs
#   --input-dir: directory containing expression, CCF coordinate, and lookup CSVs (required).
#   --output-dir: directory to write merged tables and plots (required).
#   --ccf-pattern: regex matching coordinate files (default: `^ccf_coordinates_.*\\.csv$`).
#   --expression-pattern: regex matching predicted cluster CSVs (default: `^expression-matrix-log2_.*predicted-cluster_id.*\\.csv$`).
#   --parcellation-csv: CSV mapping parcellation_index to substructure (required).
#   --substructure-csv: CSV listing substructures to retain (required, column `substructure`).
#   --subtype-groups: semicolon-separated cluster ranges (default: `NN=0:7;GABA=8:19;Glut=20:29`).
#   --pdf-width: width for output PDFs (default: 20).
#   --pdf-height: height for output PDFs (default: 6).
#
# Dependencies: optparse, readr, dplyr, tidyr, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

parse_groups <- function(spec) {
  if (!nzchar(spec)) return(list())
  pieces <- strsplit(spec, ";")[[1]]
  groups <- lapply(pieces, function(piece) {
    kv <- strsplit(piece, "=")[[1]]
    name <- kv[1]
    range <- kv[2]
    clusters <- if (grepl(":", range)) {
      bounds <- as.integer(strsplit(range, ":")[[1]])
      seq(bounds[1], bounds[2])
    } else {
      as.integer(strsplit(range, ",")[[1]])
    }
    setNames(list(clusters), name)
  })
  Reduce(function(x, y) c(x, y), groups, init = list())
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory with tangram outputs"),
    make_option("--output-dir", type = "character", help = "Directory for plots"),
    make_option("--ccf-pattern", type = "character", default = "^ccf_coordinates_.*\\.csv$"),
    make_option("--expression-pattern", type = "character",
                default = "^expression-matrix-log2_.*predicted-cluster_id.*\\.csv$"),
    make_option("--parcellation-csv", type = "character", help = "Parcellation lookup CSV"),
    make_option("--substructure-csv", type = "character", help = "Substructure list CSV"),
    make_option("--subtype-groups", type = "character", default = "NN=0:7;GABA=8:19;Glut=20:29"),
    make_option("--pdf-width", type = "double", default = 20),
    make_option("--pdf-height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "allen_region_barplot.R --input-dir DIR --output-dir DIR --parcellation-csv FILE --substructure-csv FILE [options]")
  args <- parse_args(parser)

  required <- c(args$`input-dir`, args$`output-dir`, args$`parcellation-csv`, args$`substructure-csv`)
  if (any(!nzchar(required))) stop("--input-dir, --output-dir, --parcellation-csv, --substructure-csv are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)
  if (!file.exists(args$`parcellation-csv`)) stop("Parcellation CSV not found")
  if (!file.exists(args$`substructure-csv`)) stop("Substructure CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  parcellation <- read_csv(args$`parcellation-csv`, show_col_types = FALSE)
  substructures <- read_csv(args$`substructure-csv`, show_col_types = FALSE)$substructure

  ccf_files <- list.files(args$`input-dir`, pattern = args$`ccf-pattern`, full.names = TRUE)
  expr_files <- list.files(args$`input-dir`, pattern = args$`expression-pattern`, full.names = TRUE)
  if (!length(expr_files)) stop("No expression CSVs matched pattern")
  if (length(ccf_files) != length(expr_files)) stop("Mismatched number of expression and CCF files")

  counts_list <- list()
  for (expr_path in expr_files) {
    sample_id <- sub("^expression-matrix-log2_(.*)_predicted.*$", "\\1", basename(expr_path))
    ccf_match <- ccf_files[grepl(sample_id, basename(ccf_files))]
    if (!length(ccf_match)) stop("Missing CCF file for sample ", sample_id)

    expr_df <- read_csv(expr_path, show_col_types = FALSE)
    ccf_df <- read_csv(ccf_match[1], show_col_types = FALSE) %>% select(cell_label, parcellation_index)

    merged <- expr_df %>%
      rename(cell_label = cell_id) %>%
      left_join(ccf_df, by = "cell_label") %>%
      left_join(parcellation, by = "parcellation_index") %>%
      filter(substructure %in% substructures)

    write_csv(merged, file.path(args$`output-dir`, paste0(sample_id, "_merged.csv")))

    counts <- merged %>%
      count(cluster_id, substructure) %>%
      pivot_wider(names_from = substructure, values_from = n, values_fill = 0) %>%
      arrange(cluster_id)
    write_csv(counts, file.path(args$`output-dir`, paste0(sample_id, "_counts.csv")))
    counts_list[[sample_id]] <- counts
  }

  combined_counts <- bind_rows(counts_list) %>%
    pivot_longer(-cluster_id, names_to = "substructure", values_to = "count") %>%
    group_by(cluster_id, substructure) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = substructure, values_from = count, values_fill = 0) %>%
    arrange(cluster_id)
  write_csv(combined_counts, file.path(args$`output-dir`, "combined_counts.csv"))

  observed_substructures <- setdiff(names(combined_counts), "cluster_id")
  fill_colors <- setNames(
    grDevices::hcl.colors(length(observed_substructures), palette = "Dark 3"),
    observed_substructures
  )

  plot_bar <- function(data, y_col, filename, y_label, title = NULL) {
    pdf(file.path(args$`output-dir`, filename), width = args$`pdf-width`, height = args$`pdf-height`)
    print(
      data %>%
        ggplot(aes(x = factor(cluster_id), y = .data[[y_col]], fill = substructure)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = fill_colors) +
        labs(x = "Cluster", y = y_label, title = title) +
        theme_minimal(base_size = 14)
    )
    dev.off()
  }

  long_abs <- combined_counts %>% pivot_longer(-cluster_id, names_to = "substructure", values_to = "count")
  plot_bar(long_abs, "count", "combined_barplot_absolute.pdf", "Cell ROI count")

  cluster_norm <- long_abs %>% group_by(cluster_id) %>% mutate(cluster_prop = count / sum(count)) %>% ungroup()
  plot_bar(cluster_norm, "cluster_prop", "combined_barplot_cluster_normalized.pdf", "Proportion within cluster")

  subtype_norm <- long_abs %>% group_by(substructure) %>% mutate(proportion = count / sum(count)) %>% ungroup()
  plot_bar(subtype_norm, "proportion", "combined_barplot_normalized.pdf", "Proportion of ROIs")

  double_norm <- long_abs %>%
    group_by(substructure) %>% mutate(prop_sub = count / sum(count)) %>% ungroup() %>%
    group_by(cluster_id) %>% mutate(double_norm = prop_sub / max(prop_sub)) %>% ungroup()
  plot_bar(double_norm, "double_norm", "combined_barplot_double_normalized.pdf", "Arbitrary units")

  subtype_groups <- parse_groups(args$`subtype-groups`)
  for (name in names(subtype_groups)) {
    clusters <- subtype_groups[[name]]
    subset_abs <- long_abs %>% filter(cluster_id %in% clusters)
    plot_bar(subset_abs, "count", paste0("combined_barplot_", name, "_absolute.pdf"), "Cell ROI count", paste(name, "Absolute"))

    subset_cluster_norm <- cluster_norm %>% filter(cluster_id %in% clusters)
    plot_bar(subset_cluster_norm, "cluster_prop", paste0("combined_barplot_", name, "_cluster_normalized.pdf"),
             "Proportion within cluster", paste(name, "Cluster Normalized"))

    subset_norm <- subtype_norm %>% filter(cluster_id %in% clusters)
    plot_bar(subset_norm, "proportion", paste0("combined_barplot_", name, "_normalized.pdf"),
             "Proportion of ROIs", paste(name, "Normalized"))

    subset_double <- double_norm %>% filter(cluster_id %in% clusters)
    plot_bar(subset_double, "double_norm", paste0("combined_barplot_", name, "_double_normalized.pdf"),
             "Arbitrary units", paste(name, "Double Normalized"))
  }
}

if (sys.nframe() == 0L) {
  main()
}
