#!/usr/bin/env Rscript

# Create stacked bar plots and chi-square summaries from a cluster count table.
#
# Inputs
#   --count-csv: CSV with rows as samples (including a `Total` row) and columns as clusters (required).
#   --output-dir: directory for plots and tables (required).
#   --sample-labels: comma-separated row names representing samples (default: Assigned groups).
#   --total-label: row name representing the total row (default: `Total`).
#   --color-map: optional CSV with columns `Sample` and `Color`.
#   --monte-carlo: number of Monte Carlo simulations (default: 10000).
#
# Dependencies: optparse, readr, dplyr, tidyr, ggplot2, reshape2

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(tibble)
})

run_chi_square <- function(count_table, sample_rows, total_row, simulations) {
  sample_totals <- count_table[sample_rows, "Total"]
  grand_total <- count_table[total_row, "Total"]
  expected_props <- sample_totals / grand_total
  clusters <- setdiff(colnames(count_table), "Total")

  chi_sq <- setNames(rep(NA_real_, length(clusters)), clusters)
  p_sim <- setNames(rep(NA_real_, length(clusters)), clusters)
  p_asym <- setNames(rep(NA_real_, length(clusters)), clusters)
  residuals_list <- list()

  for (cluster in clusters) {
    observed <- count_table[sample_rows, cluster]
    names(observed) <- sample_rows
    cells <- sum(observed, na.rm = TRUE)
    if (cells < 5) {
      chi_sq[cluster] <- NA
      p_sim[cluster] <- NA
      p_asym[cluster] <- NA
      residuals_list[[cluster]] <- setNames(rep(NA_real_, length(sample_rows)), sample_rows)
      next
    }

    expected <- cells * expected_props
    sim_test <- suppressWarnings(chisq.test(x = observed, p = expected_props, rescale.p = TRUE,
                                            simulate.p.value = TRUE, B = simulations))
    asym_test <- suppressWarnings(chisq.test(x = observed, p = expected_props, rescale.p = TRUE))

    chi_sq[cluster] <- sim_test$statistic
    p_sim[cluster] <- sim_test$p.value
    p_asym[cluster] <- asym_test$p.value
    residuals_list[[cluster]] <- setNames((observed - expected) / sqrt(expected), sample_rows)
  }

  chi_sq_df <- as.data.frame(t(chi_sq))
  p_sim_df <- as.data.frame(t(p_sim))
  p_asym_df <- as.data.frame(t(p_asym))
  residuals_mat <- matrix(NA_real_, nrow = length(sample_rows), ncol = length(residuals_list),
                          dimnames = list(sample_rows, names(residuals_list)))
  for (cluster in names(residuals_list)) {
    residuals_mat[, cluster] <- residuals_list[[cluster]]
  }
  residuals_df <- as.data.frame(residuals_mat)

  list(chi_sq = chi_sq_df, p_sim = p_sim_df, p_asym = p_asym_df, residuals = residuals_df)
}

main <- function() {
  option_list <- list(
    make_option("--count-csv", type = "character", help = "Cluster count table"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--sample-labels", type = "character",
                default = "Assigned: Male,Assigned: Female Diestrus,Assigned: Female Proestrus"),
    make_option("--total-label", type = "character", default = "Total"),
    make_option("--color-map", type = "character", default = NULL,
                help = "CSV mapping Sample to Color"),
    make_option("--monte-carlo", type = "integer", default = 10000,
                help = "Monte Carlo simulations [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "stacked_bar_chi_squared_from_count_table.R --count-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`count-csv`, args$`output-dir`)))) stop("--count-csv and --output-dir are required")
  if (!file.exists(args$`count-csv`)) stop("Count CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  sample_rows <- trimws(strsplit(args$`sample-labels`, ",")[[1]])
  count_df <- read.csv(args$`count-csv`, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(sample_rows %in% rownames(count_df))) stop("Sample rows missing from table")

  count_df[] <- lapply(count_df, as.numeric)
  if (!args$`total-label` %in% rownames(count_df)) {
    total_row <- colSums(count_df[sample_rows, , drop = FALSE], na.rm = TRUE)
    count_df <- rbind(count_df, Total = total_row)
    rownames(count_df)[nrow(count_df)] <- args$`total-label`
  }

  numeric_clusters <- setdiff(colnames(count_df), "Total")

  long_df <- count_df[sample_rows, , drop = FALSE] %>%
    rownames_to_column("Sample") %>%
    select(-Total) %>%
    pivot_longer(-Sample, names_to = "Cluster", values_to = "Count")

  if (!is.null(args$`color-map`) && nzchar(args$`color-map`)) {
    color_df <- read_csv(args$`color-map`, show_col_types = FALSE)
    if (!all(c("Sample", "Color") %in% names(color_df))) stop("Color map must contain Sample and Color columns")
    color_map <- setNames(color_df$Color, color_df$Sample)
  } else {
    palette <- c("black", "turquoise4", "maroon4", "goldenrod3", "gray50")
    color_map <- setNames(rep(palette, length.out = length(sample_rows)), sample_rows)
  }

  stacked_plot <- ggplot(long_df, aes(x = Cluster, y = Count, fill = Sample)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_map, guide = guide_legend(title = "Sample")) +
    theme_minimal() +
    labs(x = "Integrated Cluster", y = "Number of Barcodes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1)

  chi_results <- run_chi_square(count_df, sample_rows = sample_rows, total_row = args$`total-label`,
                                simulations = args$`monte-carlo`)

  residuals_long <- chi_results$residuals %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Cluster", values_to = "Residual")
  residual_heatmap <- ggplot(residuals_long, aes(x = Cluster, y = Sample, fill = Residual)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0,
                          limits = c(-4, 5.5), name = "Standardised\nResidual") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

  ggsave(file.path(args$`output-dir`, "stacked_bar_counts.pdf"), stacked_plot, width = 6, height = 6)
  ggsave(file.path(args$`output-dir`, "chi_squared_residuals_heatmap.pdf"), residual_heatmap, width = 6, height = 6)

  write_csv(chi_results$chi_sq, file.path(args$`output-dir`, "chi_squared_statistics.csv"))
  write_csv(chi_results$p_sim, file.path(args$`output-dir`, "pvalues_monte_carlo.csv"))
  write_csv(chi_results$p_asym, file.path(args$`output-dir`, "pvalues_asymptotic.csv"))
  write_csv(chi_results$residuals %>% rownames_to_column("Sample"),
            file.path(args$`output-dir`, "chi_squared_residuals.csv"))
}

if (sys.nframe() == 0L) {
  main()
}
