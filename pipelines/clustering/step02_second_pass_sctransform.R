#!/usr/bin/env Rscript

# Second-pass SCTransform workflow using a variable feature list, followed by PCA, UMAP, clustering, and marker export.
#
# Inputs
#   --input-h5: 10x HDF5 matrix (required).
#   --varfeat-csv: CSV of variable features from first pass (required).
#   --markers-csv: optional markers CSV for dot plot generation.
#   --output-dir: output directory for clustering results (required).
#   --dims: number of PCA dimensions (default: 75).
#   --k-param: neighbor parameter for FindNeighbors (default: 30).
#   --resolution: clustering resolution (default: 0.1).
#   --metric: distance metric for neighbors/UMAP (default: `cosine`).
#   --algorithm: Louvain/Leiden algorithm code (default: 2).
#   --umap-neighbors: number of neighbors for UMAP (default: 50).
#   --umap-repulsion: UMAP repulsion strength (default: 0.9).
#   --umap-min-dist: UMAP min.dist parameter (default: 0.4).
#   --future-max-gb: maximum future globals size in GB for SCTransform (default: 8).
#
# Dependencies: optparse, Seurat, dplyr, ggplot2, ggrepel, pheatmap, patchwork, tibble, tidyr, cluster, reticulate

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(future)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(cluster)
  library(reticulate)
})

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

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "10x HDF5 matrix"),
    make_option("--varfeat-csv", type = "character", help = "Variable feature CSV from first pass"),
    make_option("--markers-csv", type = "character", default = NULL,
                help = "Optional marker gene CSV for dot plots"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--dims", type = "integer", default = 75),
    make_option("--k-param", type = "integer", default = 30),
    make_option("--resolution", type = "double", default = 0.1),
    make_option("--metric", type = "character", default = "cosine"),
    make_option("--algorithm", type = "integer", default = 2),
    make_option("--umap-neighbors", type = "integer", default = 50),
    make_option("--umap-repulsion", type = "double", default = 0.9),
    make_option("--umap-min-dist", type = "double", default = 0.4),
    make_option("--future-max-gb", type = "double", default = 8,
                help = "Maximum future globals size in GB for SCTransform [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step02_second_pass_sctransform.R --input-h5 FILE --varfeat-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-h5`, args$`varfeat-csv`, args$`output-dir`)))) {
    stop("--input-h5, --varfeat-csv, and --output-dir are required")
  }
  if (!file.exists(args$`input-h5`)) stop("Input H5 file not found")
  if (!file.exists(args$`varfeat-csv`)) stop("Variable feature CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)
  options(future.globals.maxSize = args$`future-max-gb` * 1024^3)
  ensure_glm_gampoi()

  h5_object <- Read10X_h5(args$`input-h5`, use.names = TRUE)
  seurat_obj <- CreateSeuratObject(h5_object)
  variable_features <- read.csv(args$`varfeat-csv`, stringsAsFactors = FALSE)[[1]]
  variable_features <- intersect(variable_features, rownames(seurat_obj))

  seurat_obj <- SCTransform(seurat_obj, ncells = ncol(seurat_obj), residual.features = variable_features, verbose = TRUE)
  variable_features <- intersect(variable_features, rownames(seurat_obj))

  feature_data <- HVFInfo(seurat_obj, assay = "SCT") %>% mutate(gene = rownames(.))
  filtered_features <- feature_data[variable_features, , drop = FALSE] %>% arrange(desc(residual_variance))
  top10 <- head(filtered_features, 10)
  base_name <- tools::file_path_sans_ext(basename(args$`input-h5`))

  varfeat_plot <- VariableFeaturePlot(seurat_obj, assay = "SCT") +
    theme_minimal() +
    ggtitle("Variable Feature Plot - Second Pass") +
    geom_text_repel(data = top10, aes(x = gmean, y = residual_variance, label = gene), size = 3, box.padding = 0.5)
  ggsave(file.path(args$`output-dir`, paste0("plot_VarFeat_SecondPass_", base_name, ".pdf")), varfeat_plot)

  varfeat_log_plot <- VariableFeaturePlot(seurat_obj, assay = "SCT") +
    theme_minimal() +
    ggtitle("Variable Feature Plot - Second Pass (log10)") +
    scale_y_continuous(trans = "log10") +
    geom_text_repel(data = top10, aes(x = gmean, y = residual_variance, label = gene), size = 3, box.padding = 0.5)
  ggsave(file.path(args$`output-dir`, paste0("plot_VarFeat_log10_SecondPass_", base_name, ".pdf")), varfeat_log_plot)

  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = max(args$dims, 50), features = variable_features, verbose = TRUE)
  ggsave(file.path(args$`output-dir`, paste0("plot-Scree_", base_name, ".pdf")), ElbowPlot(seurat_obj, ndims = max(args$dims, 50)))

  pdf(file.path(args$`output-dir`, paste0("plot-PCALoadingHeatmap_", base_name, ".pdf")))
  for (dim_idx in seq_len(min(50, args$dims))) {
    DimHeatmap(seurat_obj, dims = dim_idx, cells = 500, balanced = TRUE)
  }
  dev.off()

  dims_vec <- seq_len(args$dims)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_vec, k.param = args$`k-param`,
                              nn.method = "annoy", annoy.metric = args$metric, verbose = TRUE)
  seurat_obj <- FindClusters(seurat_obj, algorithm = args$algorithm, resolution = args$resolution,
                             group.singletons = TRUE, verbose = TRUE)

  ensure_umap_learn()
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_vec, min.dist = args$`umap-min-dist`, metric = args$metric,
                        umap.method = "umap-learn", n.neighbors = args$`umap-neighbors`,
                        repulsion.strength = args$`umap-repulsion`, verbose = TRUE, densmap = FALSE)

  dims_label <- paste0(min(dims_vec), "-", max(dims_vec))
  clustering_dir <- file.path(args$`output-dir`, paste0("k", args$`k-param`, "_res", args$resolution))
  dir.create(clustering_dir, recursive = TRUE, showWarnings = FALSE)

  umap_plot <- DimPlot(seurat_obj, label = TRUE, shuffle = TRUE, label.size = 2, label.box = TRUE,
                       reduction = "umap") +
    theme_minimal() + theme(aspect.ratio = 1, legend.position = "right", legend.box = "vertical")
  ggsave(file.path(clustering_dir,
                   paste0("plot-UMAP_SecondPass_dims-", dims_label,
                          "_kparam-", args$`k-param`, "_metric-", args$metric,
                          "_algo-", args$algorithm, "_res-", args$resolution,
                          "_opt-", args$`umap-min-dist`, "_umapneighbors-", args$`umap-neighbors`,
                          "_umaprepel-", args$`umap-repulsion`, "_", base_name, ".pdf")),
         umap_plot)

  umap_coords <- Embeddings(seurat_obj, "umap")
  barcode_data <- data.frame(Barcode = rownames(umap_coords), Cluster = Idents(seurat_obj),
                              UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2])
  write.csv(barcode_data, file.path(clustering_dir,
                                    paste0("dataframe-BarcodeClusterUMAP_SecondPass_dims-", dims_label,
                                           "_kparam-", args$`k-param`, "_metric-", args$metric,
                                           "_algo-", args$algorithm, "_res-", args$resolution,
                                           "_", base_name, ".csv")), row.names = FALSE)

  markers <- FindAllMarkers(seurat_obj, test.use = "wilcox", min.pct = -Inf, min.diff.pct = -Inf,
                            verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf,
                            logfc.threshold = -Inf)
  write.csv(markers, file.path(clustering_dir,
                               paste0("dataframe-FindAllMarkers_SecondPass_dims-", dims_label,
                                      "_kparam-", args$`k-param`, "_metric-", args$metric,
                                      "_algo-", args$algorithm, "_res-", args$resolution,
                                      "_", base_name, ".csv")), row.names = FALSE)

  if (!is.null(args$`markers-csv`) && nzchar(args$`markers-csv`) && file.exists(args$`markers-csv`)) {
    marker_genes <- read.csv(args$`markers-csv`)
    create_dot_plot <- function(data, markers, output_file, max_bubble_size = 5) {
      data_filtered <- data %>%
        filter(gene %in% markers$gene) %>%
        mutate(cluster = as.numeric(as.character(cluster)),
               bubble_outline = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "black", "white"))
      data_filtered <- data_filtered %>%
        mutate(avg_log2FC = pmax(pmin(avg_log2FC, 3.5), -3.5))
      p <- ggplot(data_filtered, aes(x = factor(cluster), y = factor(gene, levels = markers$gene))) +
        geom_hline(aes(yintercept = as.numeric(factor(gene, levels = markers$gene))),
                   color = "gray70", linetype = "solid", linewidth = 0.3) +
        geom_point(aes(size = pct.1, fill = avg_log2FC, color = bubble_outline), shape = 21) +
        scale_size_continuous(range = c(1, max_bubble_size), limits = c(0.1, 1)) +
        scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick", midpoint = 0,
                             limits = c(-3.5, 3.5)) +
        scale_color_identity() +
        labs(x = "Cluster", y = "Gene", size = "Pct.1", fill = "Avg Log2FC") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              legend.position = "right", legend.box = "vertical", panel.grid = element_blank())
      ggsave(output_file, plot = p, width = 8, height = 7)
    }

    create_dot_plot(markers, marker_genes,
                    file.path(clustering_dir,
                              paste0("dotplot_markers_SecondPass_dims-", dims_label,
                                     "_kparam-", args$`k-param`, "_metric-", args$metric,
                                     "_algo-", args$algorithm, "_res-", args$resolution,
                                     "_", base_name, ".pdf")))
  }
}

if (sys.nframe() == 0L) {
  main()
}
