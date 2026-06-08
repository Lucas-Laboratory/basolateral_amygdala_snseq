
# Load required packages
library(optparse)
library(Seurat)
library(dplyr)
library(reticulate)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)

ensure_glm_gampoi <- function() {
  if (requireNamespace("glmGamPoi", quietly = TRUE)) return(invisible(TRUE))
  message("Installing Bioconductor package glmGamPoi for faster SCTransform.")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("glmGamPoi", ask = FALSE, update = FALSE)
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop("glmGamPoi installation failed; install it manually with BiocManager::install('glmGamPoi')")
  }
  invisible(TRUE)
}

ensure_umap_learn <- function() {
  if (exists("py_require", envir = asNamespace("reticulate"), inherits = FALSE)) {
    reticulate::py_require("umap-learn")
  }
  if (!reticulate::py_module_available("umap")) {
    message("Installing Python package umap-learn for Seurat RunUMAP.")
    reticulate::py_install("umap-learn", pip = TRUE)
  }
  if (!reticulate::py_module_available("umap")) {
    stop("umap-learn installation failed; install it manually with python -m pip install umap-learn")
  }
  invisible(TRUE)
}

option_list <- list(
  make_option("--seurat-rds", type = "character",
              help = "Input Seurat object saved as .rds"),
  make_option("--varfeat-csv", type = "character",
              help = "CSV containing first-pass variable features"),
  make_option("--output-dir", type = "character",
              help = "Directory for integration outputs"),
  make_option("--cluster-id", type = "character", default = "7",
              help = "Cluster identity to integrate separately (default: %default)"),
  make_option("--feature-column", type = "character", default = "x",
              help = "Column in --varfeat-csv containing feature names (default: %default)")
)

parser <- OptionParser(
  option_list = option_list,
  usage = "seurat_cluster_integration_template.R --seurat-rds FILE --varfeat-csv FILE --output-dir DIR [options]"
)
args <- parse_args(parser)

required_args <- c(args$`seurat-rds`, args$`varfeat-csv`, args$`output-dir`)
if (!all(nzchar(required_args))) {
  stop("--seurat-rds, --varfeat-csv, and --output-dir are required")
}
if (!file.exists(args$`seurat-rds`)) stop("Seurat RDS not found: ", args$`seurat-rds`)
if (!file.exists(args$`varfeat-csv`)) stop("Variable feature CSV not found: ", args$`varfeat-csv`)

input_rds <- args$`seurat-rds`
output_dir <- args$`output-dir`
var_feats_first_pass_path <- args$`varfeat-csv`
cluster_id <- args$`cluster-id`
feature_column <- args$`feature-column`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load in seurat object and check idents, load in variable feat df, cross compare and make vector
seurat_obj <- readRDS(input_rds)
ident_check <- Idents(seurat_obj)
print(head(ident_check))
feature_list <- read.csv(var_feats_first_pass_path, stringsAsFactors = FALSE)
if (!feature_column %in% colnames(feature_list)) {
    stop("Feature column not found in --varfeat-csv: ", feature_column)
}
variable_features <- feature_list[[feature_column]]
variable_features <- intersect(variable_features, rownames(seurat_obj))

# Recover cluster of interest barcodes, subset the object, create lists of barcodes in and out of the cluster of interest
cluster_barcodes <- WhichCells(seurat_obj, idents = cluster_id)
seurat_cluster <- subset(seurat_obj, cells = cluster_barcodes)
seurat_rest <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), cluster_barcodes))
seurat_cluster$integration_group <- paste0("cluster", cluster_id)
seurat_rest$integration_group <- "rest"
seurat_list <- list(seurat_rest, seurat_cluster)

# Set global memory limit to accomodate SCT matrix manipulation
options(future.globals.maxSize= 3000*1024^2)

# Apply SCTransform to all members of the seurat object
total_cells <- sum(sapply(seurat_list, ncol))
ensure_glm_gampoi()
seurat_list <- lapply(
    seurat_list, 
    function(obj) {
        SCTransform(
            obj,
            ncells = total_cells,
            residual.features = variable_features,
            do.scale = TRUE,
            do.center = TRUE,
            verbose = TRUE
        )
    }
)

# Select and prep integration features
valid_features <- Reduce(intersect, lapply(seurat_list, function(x) rownames(x[["SCT"]]@scale.data)))
features <- intersect(variable_features, valid_features)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# Find anchors
anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    normalization.method = "SCT",
    anchor.features = features
)

# Perform the integration (cluster of interest onto rest)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Define clustering parameters
dims_set <- 1:75 
k_param_set <- 75 
metric_set <- "cosine" 
algo_set <- 2 
resolution_set <- 0.3 
optimization_set <- 0.4 
umap_neighbors_set <- 50L 
repulsion_set <- 0.6
dims_set_cleaned <- paste0(min(dims_set), "-", max(dims_set))

# Run the PCA, will spit warnings, visually confirm PCA with Scree and heatmaps
integrated <- RunPCA(integrated, assay = "SCT", npcs = 100, features = variable_features, verbose = TRUE)
    scree_plot_integrated <- ElbowPlot(integrated, ndims = 100)
    dim_heatmap_integrated <- DimHeatmap(
        integrated,
        dims = 1:5,
        cells = 500,
        balanced = TRUE
    )
print(scree_plot_integrated)
print(dim_heatmap_integrated)

# Find nearest neighbors
integrated <- FindNeighbors(integrated, 
    dims = dims_set, 
    k.param = k_param_set, 
    nn.method = "annoy", 
    annoy.metric = metric_set, 
    verbose = TRUE
)

# Confirm presence of nearest neighbor graph "snn"
names(integrated@graphs)

# Perform clustering
integrated <- FindClusters(integrated,
    graph.name = "SCT_snn",
    algorithm = algo_set, 
    resolution = resolution_set, 
    group.singletons = TRUE, 
    verbose = TRUE
)

# Assign the metadata from the cluster of interest and save table of new cluster constituency
cluster_new_idents <- Idents(integrated)[cluster_barcodes]
cluster_assignments <- split(names(cluster_new_idents), cluster_new_idents)
max_len <- max(sapply(cluster_assignments, length))
padded_df <- as.data.frame(lapply(cluster_assignments, function(x) {
    length(x) <- max_len
    return(x)
}), stringsAsFactors = FALSE)
colnames(padded_df) <- paste0("cluster_", names(cluster_assignments))

# Run UMAP via umap-learn via reticulate
ensure_umap_learn()
integrated <- RunUMAP(integrated, 
      dims = dims_set, 
      min.dist = optimization_set, 
      metric = metric_set,
      umap.method = "umap-learn",
      n.neighbors = umap_neighbors_set, 
      repulsion.strength = repulsion_set, 
      verbose = TRUE,
      densmap = FALSE
      )

# Create and print integrated UMAP plot
umap_plot_integrated <- DimPlot(
    integrated, 
    label = TRUE, 
    shuffle = TRUE, 
    label.size = 2, 
    label.box = TRUE, 
    repel = FALSE, 
    reduction = "umap"
) + 
theme_minimal() + 
theme(
    aspect.ratio = 1,           
    legend.position = "right",  
    legend.box = "vertical"     
)
print(umap_plot_integrated)

# Create and original umap plot
umap_plot_original <- DimPlot(
    seurat_obj, 
    label = TRUE, 
    shuffle = TRUE, 
    label.size = 2, 
    label.box = TRUE, 
    repel = FALSE, 
    reduction = "umap"
) + 
theme_minimal() + 
theme(
    aspect.ratio = 1,           
    legend.position = "right",  
    legend.box = "vertical"     
)

# Create a UMAP from integrated object colored by COI and sample
cluster_barcodes_integrated <- intersect(cluster_barcodes, colnames(integrated))
umap_data <- Embeddings(integrated, "umap") %>%
    as.data.frame() %>%
    setNames(c("UMAP_1", "UMAP_2")) %>%
    mutate(cell = rownames(.))

umap_data$label <- "Integrator"
umap_data$label[umap_data$cell %in% cluster_barcodes_integrated & grepl("^Male-Naive_", umap_data$cell)] <- "Integrand: Male"
umap_data$label[umap_data$cell %in% cluster_barcodes_integrated & grepl("^Female-Proestrus-Naive_", umap_data$cell)] <- "Integrand: Female Proestrus"
umap_data$label[umap_data$cell %in% cluster_barcodes_integrated & grepl("^Female-Diestrus-Naive_", umap_data$cell)] <- "Integrand: Female Diestrus"


umap_data$label <- factor(umap_data$label, levels = c(
    "Integrator",
    "Integrand: Male",
    "Integrand: Female Proestrus",
    "Integrand: Female Diestrus"
))
color_map <- c(
    "Integrator" = "#CCCCCC",
    "Integrand: Male" = "#000000",
    "Integrand: Female Proestrus" = "#8B1D62",
    "Integrand: Female Diestrus" = "#348085"
)

sample_color_umap_integrated <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
    geom_point(size = 0.6) +
scale_color_manual(values = color_map, name = "Source") +
theme_minimal() +
theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.box = "vertical"
)

print(sample_color_umap_integrated)

# Create a count table
Idents(integrated) <- "seurat_clusters"
final_clusters <- Idents(integrated)
cluster_barcodes_integrated <- intersect(cluster_barcodes, colnames(integrated))
get_sample_group <- function(barcode) {
  if (grepl("^Female-Proestrus-Naive_", barcode)) {
    return("Integrand: Female Proestrus")
  } else if (grepl("^Female-Diestrus-Naive_", barcode)) {
    return("Integrand: Female Diestrus")
  } else if (grepl("^Male-Naive_", barcode)) {
    return("Integrand: Male")
  } else {
    return(NA_character_)
  }
}

integrand_meta <- data.frame(
  barcode = cluster_barcodes_integrated,
  sample_group = sapply(cluster_barcodes_integrated, get_sample_group),
  final_cluster = as.character(final_clusters[cluster_barcodes_integrated]),
  stringsAsFactors = FALSE
)

integrand_meta <- integrand_meta[!is.na(integrand_meta$sample_group), ]
summary_table <- table(integrand_meta$sample_group, integrand_meta$final_cluster)
summary_df <- as.data.frame.matrix(summary_table)
summary_df["Total", ] <- colSums(summary_df)
summary_df_integrated_counts <- summary_df[c(
  "Integrand: Female Proestrus",
  "Integrand: Female Diestrus",
  "Integrand: Male",
  "Total"
), , drop = FALSE]
summary_df_integrated_counts$Total <- rowSums(summary_df_integrated_counts)

print(summary_df_integrated_counts)

# Create stacked bar plot of counts
censored_summary_df_integrated_counts <- summary_df_integrated_counts
censored_summary_df_integrated_counts$Total <- NULL
long_df <- censored_summary_df_integrated_counts[1:3, ] %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Cluster", values_to = "Count")

color_map <- c(
    "Integrand: Male" = "black",
    "Integrand: Female Diestrus" = "turquoise4",
    "Integrand: Female Proestrus" = "maroon4"
)

stacked_counts_plot <- ggplot(long_df, aes(x = Cluster, y = Count, fill = Sample)) +
geom_bar(stat = "identity") +
scale_fill_manual(values = color_map) +
theme_minimal() +
ylab("Number of Barcodes") +
xlab("Integrated Cluster") +
theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 1,
    legend.position = "right",
    legend.box = "vertical"
)

print(stacked_counts_plot)

# Run Chi-Squared analysis to look for disproportionate distribution of COI barcodes
chi_sq_table <- setNames(numeric(), character())
p_value_sim_table <- setNames(numeric(), character())
p_value_asymptotic_table <- setNames(numeric(), character())
residuals_sim_list <- list()
residuals_asymptotic_list <- list()

sample_totals <- summary_df_integrated_counts[1:3, "Total"]
grand_total <- summary_df_integrated_counts["Total", "Total"]
expected_proportions <- sample_totals / grand_total
clusters <- setdiff(colnames(summary_df_integrated_counts), "Total")

for (cluster in clusters) {
    observed <- summary_df_integrated_counts[1:3, cluster]
    names(observed) <- rownames(summary_df_integrated_counts)[1:3]
    
    n_drawn <- sum(observed, na.rm = TRUE)
    
    if (n_drawn < 5) {
        chi_sq_table[cluster] <- NA
        p_value_sim_table[cluster] <- NA
        p_value_asymptotic_table[cluster] <- NA
        residuals_sim_list[[cluster]] <- setNames(rep(NA, 3), rownames(summary_df_integrated_counts)[1:3])
        residuals_asymptotic_list[[cluster]] <- setNames(rep(NA, 3), rownames(summary_df_integrated_counts)[1:3])
        next
    }

    expected <- n_drawn * expected_proportions

    test_sim <- suppressWarnings(
        chisq.test(x = observed, p = expected_proportions, rescale.p = TRUE, simulate.p.value = TRUE, B = 10000)
    )

    test_asymptotic <- suppressWarnings(
        chisq.test(x = observed, p = expected_proportions, rescale.p = TRUE)
    )

    chi_sq_table[cluster] <- test_sim$statistic
    p_value_sim_table[cluster] <- test_sim$p.value
    p_value_asymptotic_table[cluster] <- test_asymptotic$p.value

    residuals <- (observed - expected) / sqrt(expected)
    residuals_sim_list[[cluster]] <- setNames(residuals, names(observed))
    residuals_asymptotic_list[[cluster]] <- setNames(residuals, names(observed))
}

chi_sq_df <- as.data.frame(t(chi_sq_table))
p_val_sim_df <- as.data.frame(t(p_value_sim_table))
p_val_asymptotic_df <- as.data.frame(t(p_value_asymptotic_table))
residual_row_names <- rownames(summary_df_integrated_counts)[1:3]
residuals_df <- matrix(NA, nrow = 3, ncol = length(residuals_sim_list),
    dimnames = list(residual_row_names, names(residuals_sim_list)))
for (cluster in names(residuals_sim_list)) {
    residuals_df[, cluster] <- residuals_sim_list[[cluster]]
}
residuals_df <- as.data.frame(residuals_df)

cat("Chi-squared statistics:\n")
print(chi_sq_df)
cat("\nP-values (Monte Carlo simulation):\n")
print(p_val_sim_df)
cat("\nP-values (Asymptotic):\n")
print(p_val_asymptotic_df)
cat("\nStandardized residuals:\n")
print(residuals_df)

# Create a heatmap of chi squared residuals
residuals_long <- reshape2::melt(as.matrix(residuals_df), varnames = c("Sample", "Cluster"), value.name = "Residual")

residuals_heatmap <- ggplot(residuals_long, aes(x = Cluster, y = Sample, fill = Residual)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient2(
        low = "dodgerblue",
        mid = "white",
        high = "firebrick",
        midpoint = 0,
        limits = c(-4, 5),
        name = "Standardized\nResidual"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = (1/7)
    )

print(residuals_heatmap)

# Run FindAllMarkers on integrated object
FindAll_markers <- FindAllMarkers(
    object = integrated,
    test.use = "wilcox",
    min.pct = -Inf,
    min.diff.pct = -Inf,
    verbose = TRUE,
    only.pos = FALSE,
    max.cells.per.ident = Inf,
    logfc.threshold = -Inf
)

# Create dot plots for markers of integrated clusters
create_dot_plot <- function(data, markers, max_dot_size) {
    data_filtered <- data %>%
        filter(gene %in% markers$gene) %>%
        mutate(cluster = as.numeric(as.character(cluster)),
            dot_outline = ifelse(p_val_adj < 0.01, "black", "white")) %>%
        arrange(factor(gene, levels = markers$gene), cluster)
    
    data_filtered <- data_filtered %>%
        mutate(avg_log2FC = pmax(pmin(avg_log2FC, 3.5), -3.5))

    marker_dot_plot <- ggplot(data_filtered, aes(x = factor(cluster), y = factor(gene, levels = markers$gene))) +
        geom_hline(aes(yintercept = as.numeric(factor(gene, levels = markers$gene))), 
            color = "gray70", linetype = "solid", linewidth = 0.3) +
        geom_point(aes(size = pct.1, fill = avg_log2FC, color = dot_outline), shape = 21) +
        scale_size_continuous(range = c(1, max_dot_size), limits = c(0.1, 1)) +
        scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-3.5, 3.5)) +
        scale_color_identity() +
        labs(x = "Cluster", y = "Gene", size = "Pct.1", fill = "Avg Log2FC", color = "P-Value < 0.05 & Avg_Log2FC > 0") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            legend.position = "right",
            legend.box = "vertical",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(marker_dot_plot)
}

marker_genes <- data.frame(
    gene = c("Nnat", "Tshz2", "Ccnjl", "Rspo2", "Doc2b", "Nova1", "Tafa2", "Cdh23", "Chrm3"),
    rank = 1:9
)
max_dot_size <- 2.5
returned_cellmarker_dotplot <- create_dot_plot(FindAll_markers, marker_genes, max_dot_size)
print(returned_cellmarker_dotplot)

# Make feature plots and violin plots for integrated seurat object 
featplot_feats <- c("Nnat", "Tshz2", "Ccnjl", "Rspo2", "Ndst3", "Nova1", "Tafa2", "Cdh23")
featplot_colors <- c("lightgray", "firebrick")
featplot <- FeaturePlot(integrated, featplot_feats, cols = featplot_colors, order = TRUE)
print(featplot)
violinplot <- VlnPlot(integrated, featplot_feats)
print(violinplot)

# Make select feature plots
featplot_feats <- c("Nnat", "Rspo2")
featplot_colors <- c("lightgray", "firebrick")
featplot_select <- FeaturePlot(integrated, featplot_feats, cols = featplot_colors, order = TRUE)
print(featplot_select)

# Make activated feature plots and violin plots
featplot_feats <- c("Fosb", "Junb", "Bdnf", "Fosl2", "Homer1", "Nr4a2", "Crem", "Arc", "Ntrk2", "Btg2")
featplot_colors <- c("lightgray", "firebrick")
featplot_active <- FeaturePlot(integrated, featplot_feats, cols = featplot_colors, order = TRUE)
print(featplot_active)
violinplot_active <- VlnPlot(integrated, featplot_feats)
print(violinplot_active)


# Save plots and data tables
write.csv(padded_df, file.path(output_dir, "integrand-barcode_constits-new-clusters.csv"), row.names = FALSE)

write.csv(summary_df_integrated_counts, file = file.path(output_dir, "integrated_cluster_counts.csv"))

pdf(file.path(output_dir, "integrated_scree.pdf"))
print(scree_plot_integrated)
dev.off()

pdf(file.path(output_dir, "integrated_pca_dim_heatmaps.pdf"))
par(mfrow = c(1, 1)) # Ensure one plot per page
for (dim_idx in 1:5) {
    DimHeatmap(
        integrated,
        dims = dim_idx,
        cells = 500,
        balanced = TRUE
    )
}
dev.off()

anchor_genes <- unique(anchors@anchor.features)
write.csv(anchor_genes, file = file.path(output_dir, "used_anchor_genes.csv"), row.names = FALSE)

pdf(file.path(output_dir, "integrated_UMAP_seurat-colors.pdf"))
print(umap_plot_integrated)
dev.off()

pdf(file.path(output_dir, "original_UMAP_seurat-colors.pdf"))
print(umap_plot_original)
dev.off()

pdf(file.path(output_dir, "integrated_UMAP_COI-sample-colors.pdf"))
print(sample_color_umap_integrated)
dev.off()

pdf(file.path(output_dir, "integrated_stacked_barplot_counts-per-cluster.pdf"))
print(stacked_counts_plot)
dev.off()

write.csv(chi_sq_df, file = file.path(output_dir, "chi_squared_statistic_dataframe.csv"))

write.csv(p_val_sim_df, file = file.path(output_dir, "monte_carlo_B10k_simulated_p-values_dataframe.csv"))

write.csv(p_val_asymptotic_df, file = file.path(output_dir, "asymtotic_p-values_dataframe.csv"))

write.csv(residuals_df, file = file.path(output_dir, "chi-squared_standardized-pearson-residuals_dataframe.csv"), row.names = TRUE)

pdf(file.path(output_dir, "chi-squared_standardized-pearson-residuals_heatmap.pdf"))
print(residuals_heatmap)
dev.off()

FindAll_file_name <- paste0(
    "dataframe-FindAllMarkers_integrated_",
    "_dims-", dims_set_cleaned,
    "_kparam-", k_param_set,
    "_metric-", metric_set,
    "_algo-", algo_set,
    "_res-", resolution_set,
    ".csv"
)
output_file <- file.path(output_dir, FindAll_file_name)
write.csv(FindAll_markers, file = output_file, row.names = FALSE)

pdf(file.path(output_dir, "integrated_cell-marker_dotplot.pdf"), width = 4, height = 2)
print(returned_cellmarker_dotplot)
dev.off()

pdf(file.path(output_dir, "integrated_cell-marker_featplot.pdf"), width = 9, height = 7.5)
print(featplot)
dev.off()

pdf(file.path(output_dir, "integrated_cell-marker_featplot-select.pdf"), width = 6, height = 2.75)
print(featplot_select)
dev.off()

pdf(file.path(output_dir, "integrated_marker_featplot-active.pdf"), width = 12, height = 8)
print(featplot_active)
dev.off()

pdf(file.path(output_dir, "integrated_cell-marker_violinplot.pdf"), width = 12, height = 8)
print(violinplot)
dev.off()

pdf(file.path(output_dir, "integrated_marker_active_violinplot.pdf"), width = 12, height = 8)
print(violinplot_active)
dev.off()

saveRDS(integrated, file = file.path(output_dir, "cluster_integrated_glut_seurat_obj_consistent-variable-features.rds"))
