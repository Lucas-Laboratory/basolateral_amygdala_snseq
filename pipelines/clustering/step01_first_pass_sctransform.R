#!/usr/bin/env Rscript

# First-pass SCTransform workflow to generate variable features and diagnostic plots.
#
# Inputs
#   --input-h5: 10x HDF5 matrix (required).
#   --output-dir: directory for plots and CSV outputs (required).
#   --variable-features: number of variable features to retain (default: 5000).
#   --min-cells: minimum cells expressing a gene (default: 30).
#   --patterns-to-remove: optional regex patterns (comma-separated) to exclude genes.
#   --assay: Seurat assay to create (default: `RNA`).
#   --future-max-gb: maximum future globals size in GB for SCTransform (default: 8).
#
# Outputs
#   - Variable feature CSV and plots in the specified output directory.
#
# Dependencies: optparse, Seurat, dplyr, ggplot2, ggrepel, pheatmap, future

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(future)
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

main <- function() {
  option_list <- list(
    make_option("--input-h5", type = "character", help = "10x HDF5 matrix"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--variable-features", type = "integer", default = 5000),
    make_option("--min-cells", type = "integer", default = 30),
    make_option("--patterns-to-remove", type = "character", default = NULL,
                help = "Comma-separated regex patterns for gene exclusion"),
    make_option("--assay", type = "character", default = "RNA"),
    make_option("--future-max-gb", type = "double", default = 8,
                help = "Maximum future globals size in GB for SCTransform [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "01_first_pass_sctransform.R --input-h5 FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-h5`, args$`output-dir`)))) {
    stop("--input-h5 and --output-dir are required")
  }
  if (!file.exists(args$`input-h5`)) stop("Input H5 file not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)
  options(future.globals.maxSize = args$`future-max-gb` * 1024^3)
  ensure_glm_gampoi()

  h5_object <- Read10X_h5(args$`input-h5`, use.names = TRUE)
  seurat_obj <- CreateSeuratObject(h5_object, assay = args$assay)

  patterns_default <- c("^mt-", "^Rpl", "^Mrpl", "^Rps", "^Nduf", "\\d+Rik$", "^Gm\\d+",
                        "os$", "os\\d+$", "^[0-9]+$", "^[0-9]{4,}$", "^[A-Z][0-9]{4,}$",
                        "^[A-Z]{2}[0-9]{4,}(\\\\.[0-9]+)?$", "^D[0-9]{1,2}E.*$", "^[0-9]{5}$")
  patterns <- patterns_default
  if (!is.null(args$`patterns-to-remove`) && nzchar(args$`patterns-to-remove`)) {
    patterns <- unique(c(patterns_default, trimws(strsplit(args$`patterns-to-remove`, ",")[[1]])))
  }

  features_to_remove <- unique(rownames(seurat_obj)[grepl(paste(patterns, collapse = "|"), rownames(seurat_obj))])
  seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), features_to_remove))

  gene_counts <- rowSums(GetAssayData(seurat_obj, layer = "counts") > 0)
  filtered_genes <- names(gene_counts[gene_counts >= args$`min-cells`])
  seurat_obj <- subset(seurat_obj, features = filtered_genes)

  total_cells <- ncol(seurat_obj)
  seurat_obj <- SCTransform(seurat_obj,
                            ncells = total_cells,
                            variable.features.n = args$`variable-features`,
                            verbose = TRUE)
  variable_features <- VariableFeatures(seurat_obj, assay = "SCT")
  feature_data <- HVFInfo(seurat_obj, assay = "SCT") %>% mutate(gene = rownames(.))
  filtered_features <- feature_data[variable_features, , drop = FALSE] %>%
    arrange(desc(residual_variance))
  top10 <- head(filtered_features, 10)

  base_name <- tools::file_path_sans_ext(basename(args$`input-h5`))
  varfeat_plot <- VariableFeaturePlot(seurat_obj, assay = "SCT") +
    theme_minimal() +
    ggtitle("Variable Feature Plot - First Pass") +
    geom_text_repel(data = top10, aes(x = gmean, y = residual_variance, label = gene),
                    size = 3, box.padding = 0.5)
  ggsave(file.path(args$`output-dir`, paste0("plot_VarFeat_FirstPass_", base_name, ".pdf")), varfeat_plot)

  varfeat_log_plot <- VariableFeaturePlot(seurat_obj, assay = "SCT") +
    theme_minimal() +
    ggtitle("Variable Feature Plot - First Pass (log10)") +
    scale_y_continuous(trans = "log10") +
    geom_text_repel(data = top10, aes(x = gmean, y = residual_variance, label = gene),
                    size = 3, box.padding = 0.5)
  ggsave(file.path(args$`output-dir`, paste0("plot_VarFeat_FirstPass_log10_", base_name, ".pdf")), varfeat_log_plot)

  write.csv(variable_features,
            file.path(args$`output-dir`, paste0("dataframe_VarFeat_FirstPass_", base_name, ".csv")),
            row.names = FALSE)
}

if (sys.nframe() == 0L) {
  main()
}
