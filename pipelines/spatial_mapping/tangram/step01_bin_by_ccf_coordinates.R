#!/usr/bin/env Rscript

# Bin tangram spatial coordinates and aggregate cell counts by cluster.
#
# Inputs
#   --input-dir: directory containing expression matrices and CCF coordinate CSVs (required).
#   --output-dir: directory for merged and binned outputs (default: `<input>/../Output_binning`).
#   --cluster-range: cluster identifiers to include (e.g. `0:29`) [default: 0:29].
#   --bin-width: bin width in mm for each axis (default: 0.1).
#
# Assumes matching filename stems, e.g. `expression-matrix-log2_ABCA-1.csv` and `ccf_coordinates_ABCA-1.csv`.
#
# Dependencies: optparse, readr, dplyr, tidyr, stringr

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

bin_and_count <- function(df, coord, breaks, clusters) {
  df %>%
    filter(!is.na(.data[[coord]])) %>%
    mutate(bin = cut(.data[[coord]], breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
    group_by(bin, cluster_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(bin, cluster_id = clusters, fill = list(count = 0)) %>%
    pivot_wider(names_from = cluster_id, values_from = count)
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory with expression and CCF CSVs"),
    make_option("--output-dir", type = "character", default = NULL,
                help = "Directory for binned outputs"),
    make_option("--cluster-range", type = "character", default = "0:29",
                help = "Cluster range (e.g. 0:29 or 0,2,4)"),
    make_option("--bin-width", type = "double", default = 0.1,
                help = "Bin width in mm [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step01_bin_by_ccf_coordinates.R --input-dir DIR [options]")
  args <- parse_args(parser)

  if (is.null(args$`input-dir`) || !nzchar(args$`input-dir`)) stop("--input-dir is required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found: ", args$`input-dir`)

  clusters <- if (grepl(":", args$`cluster-range`)) {
    bounds <- as.integer(strsplit(args$`cluster-range`, ":")[[1]])
    seq(bounds[1], bounds[2])
  } else {
    as.integer(strsplit(args$`cluster-range`, ",")[[1]])
  }

  output_dir <- args$`output-dir`
  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- sub("Input_?binning$", "Output_binning", args$`input-dir`, ignore.case = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  expr_files <- list.files(args$`input-dir`, pattern = "^expression-matrix-log2_.*\\.csv$", full.names = TRUE)
  ccf_files <- list.files(args$`input-dir`, pattern = "^ccf_coordinates_.*\\.csv$", full.names = TRUE)
  if (!length(expr_files) || !length(ccf_files)) {
    stop("Input directory must contain expression and CCF coordinate CSVs")
  }

  merged_list <- lapply(expr_files, function(expr_file) {
    abca <- str_extract(basename(expr_file), "ABCA-\\d+")
    ccf_file <- ccf_files[str_detect(basename(ccf_files), paste0("ccf_coordinates_", abca, "\\.csv$"))]
    if (!length(ccf_file)) stop("Missing CCF file for ", abca)

    expr_df <- read_csv(expr_file, show_col_types = FALSE) %>%
      select(cell_id, brain_section_label, cluster_id)
    ccf_df <- read_csv(ccf_file, show_col_types = FALSE)

    expr_df %>%
      left_join(ccf_df %>% select(cell_label, x, y, z), by = c("cell_id" = "cell_label")) %>%
      mutate(abca = abca)
  })

  all_merged <- bind_rows(merged_list)
  axis_breaks <- lapply(c("x", "y", "z"), function(axis) {
    vals <- all_merged[[axis]]
    seq(floor(min(vals, na.rm = TRUE)), ceiling(max(vals, na.rm = TRUE)), by = args$`bin-width`)
  })
  names(axis_breaks) <- c("x", "y", "z")

  binned_list <- list(x = list(), y = list(), z = list())

  for (df in merged_list) {
    abca <- unique(df$abca)
    write_csv(df %>% select(-abca), file.path(output_dir, paste0("merged_", abca, ".csv")))
    for (axis in names(axis_breaks)) {
      binned_wide <- bin_and_count(df, axis, axis_breaks[[axis]], clusters)
      write_csv(binned_wide, file.path(output_dir, paste0("binned_", abca, "_", axis, ".csv")))
      binned_list[[axis]][[abca]] <- binned_wide
    }
  }

  for (axis in names(axis_breaks)) {
    combined <- bind_rows(binned_list[[axis]]) %>%
      group_by(bin) %>%
      summarise(across(everything(), sum), .groups = "drop")
    write_csv(combined, file.path(output_dir, paste0("combined_binned_", axis, ".csv")))
  }
}

if (sys.nframe() == 0L) {
  main()
}
