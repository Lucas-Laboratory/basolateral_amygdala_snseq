#!/usr/bin/env Rscript

# Train a Garnett classifier on selected Seurat clusters and classify target cells.
#
# Required inputs
#   --seurat-rds: serialized Seurat object containing RNA assay counts.
#   --marker-file: Garnett marker definition file.
#   --output-dir: directory to receive classifier and classification outputs.
#   --target-clusters: comma-separated cluster ids to classify; the remainder are used for training.
#
# Optional arguments
#   --training-clusters: comma-separated cluster ids to use for training (overrides default).
#   --cluster-id-column: metadata column identifying clusters (default: `seurat_clusters`).
#   --rank-prob-ratio: probability ratio used by `classify_cells` (default: 1.3).
#   --sample-prefixes: comma-separated list of sample prefixes for barcode grouping (e.g. `SampleA,SampleB`).
#   --chi-square-simulations: Monte Carlo draws for chi-square tests (default: 10000).
#   --skip-marker-plot: disable marker quality plotting.
#   --organism-db: Bioconductor OrgDb package for gene mapping (default: `org.Mm.eg.db`).
#
# Outputs
#   - Garnett classifier RDS, prediction tables, optional marker quality plot, and chi-square summaries.
#
# Dependencies: optparse, Seurat, SeuratWrappers, SingleCellExperiment, garnett, monocle,
#               Matrix, Biobase, dplyr, tibble, readr, ggplot2, tidyr, reshape2

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratWrappers)
  library(SingleCellExperiment)
  library(garnett)
  library(monocle)
  library(Matrix)
  library(Biobase)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(reshape2)
})

prepare_cds <- function(seurat_obj, cluster_id_column) {
  expr_matrix <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "counts")
  pheno <- seurat_obj@meta.data
  pheno[[cluster_id_column]] <- as.character(pheno[[cluster_id_column]])
  fd <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))
  cds <- newCellDataSet(
    expr_matrix,
    phenoData = new("AnnotatedDataFrame", data = pheno),
    featureData = new("AnnotatedDataFrame", data = fd),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )
  cds <- estimateSizeFactors(cds)
  cds@featureData@data$gene_id <- rownames(cds)
  cds@featureData@data$gene_short_name <- cds@featureData@data$gene_short_name
  cds@phenoData@data$cluster <- as.character(seurat_obj@meta.data[[cluster_id_column]])
  cds
}

run_chi_square <- function(count_table, sample_labels, simulations) {
  sample_totals <- count_table[sample_labels, "Total"]
  grand_total <- count_table["Total", "Total"]
  expected_props <- sample_totals / grand_total
  clusters <- setdiff(colnames(count_table), "Total")

  chi_sq <- setNames(numeric(), character())
  p_sim <- setNames(numeric(), character())
  p_asym <- setNames(numeric(), character())
  residuals_list <- list()

  for (cluster in clusters) {
    observed <- count_table[sample_labels, cluster]
    names(observed) <- sample_labels
    cells_in_cluster <- sum(observed, na.rm = TRUE)
    if (cells_in_cluster < 5) {
      chi_sq[cluster] <- NA
      p_sim[cluster] <- NA
      p_asym[cluster] <- NA
      residuals_list[[cluster]] <- setNames(rep(NA_real_, length(sample_labels)), sample_labels)
      next
    }

    expected <- cells_in_cluster * expected_props
    sim_test <- suppressWarnings(chisq.test(x = observed, p = expected_props, rescale.p = TRUE,
                                            simulate.p.value = TRUE, B = simulations))
    asym_test <- suppressWarnings(chisq.test(x = observed, p = expected_props, rescale.p = TRUE))

    chi_sq[cluster] <- sim_test$statistic
    p_sim[cluster] <- sim_test$p.value
    p_asym[cluster] <- asym_test$p.value
    residuals_list[[cluster]] <- setNames((observed - expected) / sqrt(expected), sample_labels)
  }

  chi_sq_df <- as.data.frame(t(chi_sq))
  p_sim_df <- as.data.frame(t(p_sim))
  p_asym_df <- as.data.frame(t(p_asym))
  residuals_mat <- matrix(NA_real_, nrow = length(sample_labels), ncol = length(residuals_list),
                          dimnames = list(sample_labels, names(residuals_list)))
  for (cluster in names(residuals_list)) {
    residuals_mat[, cluster] <- residuals_list[[cluster]]
  }
  residuals_df <- as.data.frame(residuals_mat)

  list(chi_sq = chi_sq_df, p_sim = p_sim_df, p_asym = p_asym_df, residuals = residuals_df)
}

plot_residual_heatmap <- function(residuals_df) {
  melted <- melt(as.matrix(residuals_df), varnames = c("Sample", "Cluster"), value.name = "Residual")
  ggplot(melted, aes(x = Cluster, y = Sample, fill = Residual)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0,
                         name = "Standardised\nResidual") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Path to Seurat .rds"),
    make_option("--marker-file", type = "character", help = "Garnett marker file"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--target-clusters", type = "character", help = "Comma-separated cluster ids to classify"),
    make_option("--training-clusters", type = "character", default = NULL,
                help = "Comma-separated cluster ids for training (optional)"),
    make_option("--cluster-id-column", type = "character", default = "seurat_clusters",
                help = "Metadata column containing cluster ids [default %default]"),
    make_option("--rank-prob-ratio", type = "double", default = 1.3,
                help = "rank_prob_ratio passed to classify_cells [default %default]"),
    make_option("--sample-prefixes", type = "character", default = NULL,
                help = "Comma-separated barcode prefixes for sample grouping"),
    make_option("--chi-square-simulations", type = "integer", default = 10000,
                help = "Monte Carlo draws for chi-square tests [default %default]"),
    make_option("--skip-marker-plot", action = "store_true", default = FALSE,
                help = "Skip marker quality plotting"),
    make_option("--organism-db", type = "character", default = "org.Mm.eg.db",
                help = "OrgDb package for gene mapping [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "run_garnett_cluster_classifier.R --seurat-rds FILE --marker-file FILE --output-dir DIR --target-clusters A,B [options]")
  args <- parse_args(parser)

  required <- c(args$`seurat-rds`, args$`marker-file`, args$`output-dir`, args$`target-clusters`)
  if (any(!nzchar(required))) stop("Arguments --seurat-rds, --marker-file, --output-dir, and --target-clusters are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat object not found: ", args$`seurat-rds`)
  if (!file.exists(args$`marker-file`)) stop("Marker file not found: ", args$`marker-file`)

  if (!requireNamespace(args$`organism-db`, quietly = TRUE)) {
    stop("OrgDb package not installed: ", args$`organism-db`)
  }
  org_db <- get(args$`organism-db`)

  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  seurat_obj <- readRDS(args$`seurat-rds`)
  if (!args$`cluster-id-column` %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster metadata column not found: ", args$`cluster-id-column`)
  }

  Idents(seurat_obj) <- seurat_obj@meta.data[[args$`cluster-id-column`]]
  target_clusters <- trimws(strsplit(args$`target-clusters`, ",")[[1]])
  training_clusters <- args$`training-clusters`
  if (is.null(training_clusters) || !nzchar(training_clusters)) {
    training_clusters <- setdiff(levels(Idents(seurat_obj)), target_clusters)
  } else {
    training_clusters <- trimws(strsplit(training_clusters, ",")[[1]])
  }

  if (!length(training_clusters)) stop("No training clusters available")

  training_subset <- subset(seurat_obj, idents = training_clusters)
  target_subset <- subset(seurat_obj, idents = target_clusters)

  cds_train <- prepare_cds(training_subset, args$`cluster-id-column`)
  cds_target <- prepare_cds(target_subset, args$`cluster-id-column`)

  classifier <- train_cell_classifier(
    cds = cds_train,
    marker_file = args$`marker-file`,
    db = org_db,
    cds_gene_id_type = "SYMBOL"
  )
  saveRDS(classifier, file = file.path(args$`output-dir`, "garnett_classifier.rds"))

  if (!args$`skip-marker-plot`) {
    marker_check <- check_markers(
      cds = cds_train,
      marker_file = args$`marker-file`,
      db = org_db,
      cds_gene_id_type = "SYMBOL",
      marker_file_gene_id_type = "SYMBOL"
    )
    marker_plot <- plot_markers(marker_check)
    ggsave(file.path(args$`output-dir`, "marker_quality.pdf"), marker_plot, width = 12, height = 18, units = "in")
  }

  classified <- classify_cells(
    cds = cds_target,
    classifier = classifier,
    db = org_db,
    cluster_extend = TRUE,
    cds_gene_id_type = "SYMBOL",
    rank_prob_ratio = args$`rank-prob-ratio`
  )

  predictions <- tibble(
    cell_barcode = colnames(classified),
    predicted_cell_type = classified@phenoData@data$cell_type
  )
  write_csv(predictions, file.path(args$`output-dir`, "garnett_predictions.csv"))

  wide_predictions <- predictions %>%
    group_by(predicted_cell_type) %>%
    summarise(barcodes = list(cell_barcode), .groups = "drop") %>%
    unnest_wider(barcodes, names_sep = "_")
  write_csv(wide_predictions, file.path(args$`output-dir`, "garnett_predictions_wide.csv"))

  prefixes <- args$`sample-prefixes`
  if (!is.null(prefixes) && nzchar(prefixes)) {
    prefix_vec <- trimws(strsplit(prefixes, ",")[[1]])
    match_mat <- sapply(prefix_vec, function(prefix) startsWith(predictions$cell_barcode, paste0(prefix, "_")))
    if (!is.matrix(match_mat)) match_mat <- matrix(match_mat, ncol = length(prefix_vec))
    colnames(match_mat) <- prefix_vec
    sample_assignments <- apply(match_mat, 1, function(row) {
      if (!any(row)) return("Other")
      paste0("Assigned: ", colnames(match_mat)[which(row)[1]])
    })
    predictions$sample_group <- sample_assignments

    count_table <- predictions %>%
      filter(startsWith(sample_group, "Assigned:")) %>%
      count(sample_group, predicted_cell_type, name = "count") %>%
      pivot_wider(names_from = predicted_cell_type, values_from = count, values_fill = 0) %>%
      arrange(sample_group)

    if (!nrow(count_table)) {
      warning("No barcodes matched the supplied sample prefixes; skipping chi-square summaries.")
    } else {
      numeric_cols <- names(count_table)[sapply(count_table, is.numeric)]
      count_table$Total <- rowSums(count_table[numeric_cols])
      total_row <- count_table[1, ]
      total_row$sample_group <- "Total"
      total_row[1, numeric_cols] <- colSums(count_table[numeric_cols])
      count_table <- bind_rows(count_table, total_row)

      write_csv(count_table, file.path(args$`output-dir`, "garnett_counts_by_sample.csv"))

      sample_rows <- count_table$sample_group[grepl("^Assigned:", count_table$sample_group)]
      chi_results <- run_chi_square(count_table, sample_labels = sample_rows,
                                    simulations = args$`chi-square-simulations`)

      write_csv(chi_results$chi_sq, file.path(args$`output-dir`, "garnett_chi_square_statistics.csv"))
      write_csv(chi_results$p_sim, file.path(args$`output-dir`, "garnett_pvalues_monte_carlo.csv"))
      write_csv(chi_results$p_asym, file.path(args$`output-dir`, "garnett_pvalues_asymptotic.csv"))
      write_csv(chi_results$residuals, file.path(args$`output-dir`, "garnett_residuals.csv"))

      heatmap <- plot_residual_heatmap(chi_results$residuals)
      ggsave(file.path(args$`output-dir`, "garnett_residual_heatmap.pdf"), heatmap, width = 6, height = 4)
    }
  }
}

if (sys.nframe() == 0L) {
  main()
}
