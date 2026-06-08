#!/usr/bin/env Rscript

# Perform manual SoupX-based background correction for 10x Genomics data.
#
# Required inputs
#   --filtered-h5: 10x `filtered_feature_bc_matrix.h5` file containing cell-associated droplets.
#   --raw-h5: 10x `raw_feature_bc_matrix.h5` file containing all droplets.
#   --output-dir: directory to write diagnostic plots, intermediate objects, and corrected outputs.
#   --sample-label: identifier used to prefix output filenames.
#
# Optional inputs
#   --annotation-csv: table with external annotations to merge (must contain barcode column).
#   --annotation-barcode-column: column in the annotation table containing barcodes (default `cell_id`).
#   --metadata-barcode-column: Seurat metadata column holding barcodes after conversion (default `.cell`).
#   --cell-type-column: metadata column describing cell types used to build SoupX exclusion sets (default `cell_type`).
#   --non-expressed-config: CSV with columns `cell_type` and `gene` describing gene sets that should be absent in each cell type.
#   --marker-csv: CSV listing genes to visualise in SoupX change/marker plots (single column expected).
#   --umi-threshold: knee point (UMI count) highlighted on the droplet rank plot (default 500).
#   --dims: number of PCA dimensions to use for UMAP/neighbor graph (default 50).
#   --resolution: clustering resolution passed to `FindClusters` (default 0.6).
#
# Outputs
#   - Seurat object (`seurat-object_<sample>.rds`).
#   - Droplet rank PDF and UMAP PDFs.
#   - SoupX background statistics and diagnostic plots.
#   - Corrected count matrix as 10x HDF5 (`soupx-corrected-counts_<sample>.h5`) and RDS (`SoupX-corrected-data_<sample>.rds`).
#   - Serialized `SoupChannel` object and marker tables.
#
# Dependencies: optparse, Seurat, SoupX, DropletUtils, ggplot2, dplyr, tibble, Matrix, hdf5r, future

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SoupX)
  library(DropletUtils)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(hdf5r)
  library(future)
  library(readr)
  library(stringr)
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

options(future.globals.maxSize = 8 * 1024^3)

read_marker_table <- function(path) {
  if (is.null(path) || !nzchar(path)) return(character())
  tbl <- readr::read_csv(path, col_types = cols(.default = "c"))
  unique(na.omit(unlist(tbl)))
}

read_non_expressed_sets <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(list(
      Oligodendrocyte = c("Mbp", "Mog"),
      Microglia = c("Tmem119", "Itgam"),
      Astrocyte = c("Gfap"),
      Glutamatergic = c("Slc17a6"),
      Gabaergic = c("Slc32a1")
    ))
  }

  cfg <- readr::read_csv(path, col_types = cols(cell_type = col_character(), gene = col_character()))
  split(cfg$gene, cfg$cell_type)
}

prepare_seurat <- function(filtered_matrix, dims, resolution) {
  seurat_obj <- CreateSeuratObject(counts = filtered_matrix)
  ensure_glm_gampoi()
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = seq_len(dims), verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = seq_len(dims), verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  seurat_obj
}

merge_annotation <- function(seurat_obj, annotation_csv, metadata_barcode_column, annotation_barcode_column) {
  if (is.null(annotation_csv) || !nzchar(annotation_csv)) {
    return(seurat_obj)
  }

  annot <- readr::read_csv(annotation_csv, show_col_types = FALSE)
  if (!annotation_barcode_column %in% names(annot)) {
    stop("Annotation file missing column: ", annotation_barcode_column)
  }

  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column(var = metadata_barcode_column)

  merged <- dplyr::left_join(meta, annot, by = setNames(annotation_barcode_column, metadata_barcode_column))
  rownames(merged) <- merged[[metadata_barcode_column]]
  merged[[metadata_barcode_column]] <- NULL
  seurat_obj@meta.data <- merged
  seurat_obj
}

plot_umaps <- function(seurat_obj, output_dir, sample_label, cluster_column, cell_type_column) {
  cluster_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = cluster_column,
                          label = TRUE, repel = TRUE) +
    ggtitle(paste("UMAP by", cluster_column)) +
    theme_minimal()

  cluster_path <- file.path(output_dir, paste0("UMAP_by_", cluster_column, "_", sample_label, ".pdf"))
  ggsave(cluster_path, plot = cluster_plot, width = 6, height = 5, useDingbats = FALSE)

  if (!is.null(cell_type_column) && cell_type_column %in% colnames(seurat_obj@meta.data)) {
    type_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = cell_type_column,
                         label = TRUE, repel = TRUE) +
      ggtitle(paste("UMAP by", cell_type_column)) +
      theme_minimal()
    type_path <- file.path(output_dir, paste0("UMAP_by_", cell_type_column, "_", sample_label, ".pdf"))
    ggsave(type_path, plot = type_plot, width = 6, height = 5, useDingbats = FALSE)
  }
}

plot_droplet_rank <- function(raw_matrix, filtered_matrix, output_dir, sample_label, umi_threshold) {
  total_raw <- colSums(raw_matrix)
  total_filtered <- colSums(filtered_matrix)

  droplet_df <- data.frame(
    barcode = names(total_raw),
    total_umis = total_raw
  ) %>%
    filter(total_umis > 0) %>%
    arrange(desc(total_umis)) %>%
    mutate(
      droplet_index = row_number(),
      status = ifelse(barcode %in% names(total_filtered), "Filtered", "Background")
    )

  droplet_df$status <- factor(droplet_df$status, levels = c("Filtered", "Background"))
  log_breaks <- c(5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1e6, 5e6, 1e7)

  p <- ggplot(droplet_df, aes(x = droplet_index, y = total_umis, color = status)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = umi_threshold, linetype = "dotted", color = "firebrick3") +
    annotate("text", x = max(droplet_df$droplet_index) * 0.2, y = umi_threshold * 1.2,
             label = paste("Threshold =", umi_threshold, "UMIs"), color = "firebrick3", size = 3) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(breaks = log_breaks, labels = scales::comma) +
    scale_color_manual(values = c("Filtered" = "dodgerblue3", "Background" = "black")) +
    annotation_logticks(sides = "lb") +
    labs(title = "Droplet Rank vs UMI Count", x = "Droplets", y = "Total UMIs") +
    theme_minimal()

  out_path <- file.path(output_dir, paste0("droplet_rank_vs_UMI_", sample_label, ".pdf"))
  ggsave(out_path, plot = p, width = 7, height = 5, useDingbats = FALSE)
}

create_use_to_est <- function(seurat_obj, cell_type_column, gene_sets) {
  meta <- seurat_obj@meta.data
  if (!cell_type_column %in% colnames(meta)) {
    stop("Metadata column not found for cell types: ", cell_type_column)
  }

  barcodes <- rownames(meta)
  cell_types <- meta[[cell_type_column]]

  matrix(FALSE, nrow = length(barcodes), ncol = length(gene_sets),
         dimnames = list(barcodes, names(gene_sets))) %>%
    {use_matrix <- .;
     for (i in seq_along(barcodes)) {
       cell_type <- cell_types[i]
       if (is.na(cell_type)) next
       excluded <- setdiff(names(gene_sets), cell_type)
       use_matrix[i, excluded] <- TRUE
     }
     use_matrix
    }
}

write_10x_h5 <- function(counts, output_path) {
  counts <- as(counts, "dgCMatrix")
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_path)) file.remove(output_path)

  h5 <- H5File$new(output_path, mode = "w")
  on.exit(h5$close_all())

  h5$create_group("matrix")
  h5$create_group("matrix/features")

  h5[["matrix/data"]] <- counts@x
  h5[["matrix/indices"]] <- counts@i
  h5[["matrix/indptr"]] <- counts@p
  h5[["matrix/shape"]] <- c(nrow(counts), ncol(counts))
  h5[["matrix/features/id"]] <- rownames(counts)
  h5[["matrix/features/name"]] <- rownames(counts)
  h5[["matrix/barcodes"]] <- colnames(counts)
}

run_soupx_pipeline <- function(filtered_matrix, raw_matrix, seurat_obj, cell_type_column,
                               non_expressed_sets, marker_genes, output_dir, sample_label) {
  message("Initialising SoupChannel")
  sc <- SoupChannel(raw_matrix, filtered_matrix, calcSoupProfile = FALSE)
  sc <- estimateSoup(sc, soupRange = c(0, 500), keepDroplets = FALSE)

  umap_embeddings <- Embeddings(seurat_obj, "umap")
  sc <- setDR(sc, umap_embeddings)
  sc <- setClusters(sc, setNames(as.character(Idents(seurat_obj)), colnames(seurat_obj)))

  if (!cell_type_column %in% colnames(seurat_obj@meta.data)) {
    stop("Cell type column not present in metadata: ", cell_type_column)
  }

  use_to_est <- create_use_to_est(seurat_obj, cell_type_column, non_expressed_sets)
  sc <- calculateContaminationFraction(sc, non_expressed_sets, use_to_est)

  if (!"rho" %in% colnames(sc$metaData)) {
    stop("SoupX contamination fraction was not estimated; 'rho' missing from metadata")
  }

  markers <- quickMarkers(sc$toc, sc$metaData$clusters, N = 10)
  write.csv(markers, file.path(output_dir, paste0("quick_markers_", sample_label, ".csv")), row.names = TRUE)

  sc <- autoEstCont(sc, tfidfMin = 1, topMarkers = NULL, doPlot = TRUE, verbose = TRUE)

  corrected_counts <- adjustCounts(sc, verbose = 1, clusters = NULL, method = "subtraction", roundToInt = TRUE)

  soup_profile <- sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ]
  write.csv(soup_profile, file.path(output_dir, paste0("background_gene_counts_", sample_label, ".csv")), row.names = TRUE)

  umap_df <- data.frame(
    x = sc$metaData$umap_1,
    y = sc$metaData$umap_2
  )
  rownames(umap_df) <- rownames(sc$metaData)

  top_background <- head(markers$gene, 50)
  change_map_path <- file.path(output_dir, paste0("plotAll_top50_ChangeMap_", sample_label, ".pdf"))
  pdf(change_map_path)
  for (gene_id in top_background) {
    p <- plotChangeMap(sc, cleanedMatrix = corrected_counts, geneSet = gene_id, DR = umap_df) +
      ggtitle(paste(gene_id, "-", sample_label))
    print(p)
  }
  dev.off()

  if (length(marker_genes)) {
    custom_change_map <- file.path(output_dir, paste0("plotAll_markers_ChangeMap_", sample_label, ".pdf"))
    pdf(custom_change_map)
    for (gene_id in marker_genes[marker_genes %in% rownames(sc$toc)]) {
      p <- plotChangeMap(sc, cleanedMatrix = corrected_counts, geneSet = gene_id, DR = umap_df) +
        ggtitle(paste(gene_id, "-", sample_label))
      print(p)
    }
    dev.off()
  }

  marker_distribution_path <- file.path(output_dir, paste0("plotMarkerDistribution_", sample_label, ".pdf"))
  pdf(marker_distribution_path)
  plotMarkerDistribution(sc, non_expressed_sets, maxCells = 150)
  dev.off()

  if (length(marker_genes)) {
    marker_map_path <- file.path(output_dir, paste0("plotMarkerMap_markers_", sample_label, ".pdf"))
    pdf(marker_map_path)
    for (gene_id in marker_genes[marker_genes %in% rownames(sc$toc)]) {
      p <- plotMarkerMap(sc, geneSet = gene_id, DR = umap_df) +
        ggtitle(paste(gene_id, "-", sample_label))
      print(p)
    }
    dev.off()
  }

  sink(file.path(output_dir, paste0("SoupChannel_summary_", sample_label, ".txt")))
  print(str(sc))
  sink()

  list(sc = sc, corrected_counts = corrected_counts)
}

main <- function() {
  option_list <- list(
    make_option("--filtered-h5", type = "character", help = "Filtered 10x HDF5 file"),
    make_option("--raw-h5", type = "character", help = "Raw 10x HDF5 file"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--sample-label", type = "character", help = "Identifier used in output filenames"),
    make_option("--annotation-csv", type = "character", default = NULL,
                help = "Optional CSV with external annotations"),
    make_option("--annotation-barcode-column", type = "character", default = "cell_id",
                help = "Barcode column in the annotation CSV [default %default]"),
    make_option("--metadata-barcode-column", type = "character", default = "cell_id",
                help = "Temporary column name for Seurat metadata barcodes [default %default]"),
    make_option("--cell-type-column", type = "character", default = "cell_type",
                help = "Metadata column holding cell type labels used for SoupX gene sets"),
    make_option("--non-expressed-config", type = "character", default = NULL,
                help = "CSV with `cell_type,gene` columns defining SoupX exclusion gene sets"),
    make_option("--marker-csv", type = "character", default = NULL,
                help = "CSV listing marker genes to visualise in diagnostic plots"),
    make_option("--umi-threshold", type = "double", default = 500,
                help = "UMI count threshold annotated on droplet rank plot [default %default]"),
    make_option("--dims", type = "integer", default = 50,
                help = "Number of PCA dimensions for neighbors/UMAP [default %default]"),
    make_option("--resolution", type = "double", default = 0.6,
                help = "Seurat clustering resolution [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "soupx_manual_correction.R --filtered-h5 FILE --raw-h5 FILE --output-dir DIR --sample-label LABEL [options]")
  args <- parse_args(parser)

  filtered_h5 <- args$`filtered-h5`
  raw_h5 <- args$`raw-h5`
  output_dir <- args$`output-dir`
  sample_label <- args$`sample-label`

  if (any(vapply(list(filtered_h5, raw_h5, output_dir, sample_label), function(x) is.null(x) || !nzchar(x), logical(1)))) {
    stop("Required arguments: --filtered-h5, --raw-h5, --output-dir, --sample-label")
  }
  if (!file.exists(filtered_h5)) stop("Filtered matrix not found: ", filtered_h5)
  if (!file.exists(raw_h5)) stop("Raw matrix not found: ", raw_h5)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  filtered_matrix <- Read10X_h5(filtered_h5, use.names = TRUE)
  raw_matrix <- Read10X_h5(raw_h5, use.names = TRUE)

  seurat_obj <- prepare_seurat(filtered_matrix, dims = args$dims, resolution = args$resolution)
  seurat_obj <- merge_annotation(seurat_obj, args$`annotation-csv`, args$`metadata-barcode-column`, args$`annotation-barcode-column`)
  if (!args$`cell-type-column` %in% colnames(seurat_obj@meta.data) &&
      "seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    message("Cell type column '", args$`cell-type-column`, "' not found; using seurat_clusters for SoupX estimation.")
    args$`cell-type-column` <- "seurat_clusters"
  }

  plot_umaps(seurat_obj, output_dir, sample_label, "seurat_clusters", args$`cell-type-column`)
  plot_droplet_rank(raw_matrix, filtered_matrix, output_dir, sample_label, args$`umi-threshold`)

  seurat_rds <- file.path(output_dir, paste0("seurat-object_", sample_label, ".rds"))
  saveRDS(seurat_obj, seurat_rds)

  marker_genes <- read_marker_table(args$`marker-csv`)
  non_expressed_sets <- read_non_expressed_sets(args$`non-expressed-config`)

  soupx_results <- run_soupx_pipeline(
    filtered_matrix = filtered_matrix,
    raw_matrix = raw_matrix,
    seurat_obj = seurat_obj,
    cell_type_column = args$`cell-type-column`,
    non_expressed_sets = non_expressed_sets,
    marker_genes = marker_genes,
    output_dir = output_dir,
    sample_label = sample_label
  )

  corrected_counts <- soupx_results$corrected_counts
  sc <- soupx_results$sc

  corrected_h5 <- file.path(output_dir, paste0("soupx-corrected-counts_", sample_label, ".h5"))
  write_10x_h5(corrected_counts, corrected_h5)

  corrected_rds <- file.path(output_dir, paste0("SoupX-corrected-data_", sample_label, ".rds"))
  saveRDS(corrected_counts, corrected_rds)

  sc_rds <- file.path(output_dir, paste0("SoupChannel_", sample_label, ".rds"))
  saveRDS(sc, sc_rds)
}

if (sys.nframe() == 0L) {
  main()
}
