#!/usr/bin/env Rscript

# Perform cluster-wise GO ORA with optional semantic similarity simplification,
# splitting results by direction of change (positive/negative logFC).
#
# Inputs
#   --input-csv: merged DEG CSV with columns `gene`, `avg_log2FC`, `cluster`,
#                `p_val_adj` (required).
#   --output-dir: directory where per-ontology CSVs will be written (required).
#   --comparison-label: label embedded in output filenames (default: derived from
#                       input filename).
#   --ontologies: comma-separated GO ontologies (default: BP,MF,CC).
#   --orgdb: OrgDb package name (default: org.Mm.eg.db).
#   --p-adj-cutoff: adjusted p-value threshold for selecting genes (default: 0.05).
#   --min-logfc: minimum absolute log2 fold-change to retain (default: 0).
#   --simplify-cutoff: semantic similarity cutoff for `simplify` (default: 0.7).
#   --min-gene-count: minimum Entrez IDs required to run enrichment (default: 5).
#
# Outputs
#   - Raw ORA CSVs per ontology/direction.
#   - Simplified ORA CSVs per ontology/direction.
#
# Dependencies: optparse, readr, dplyr, clusterProfiler, GOSemSim

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(GOSemSim)
})

parse_vector <- function(spec) {
  vals <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

load_orgdb <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("OrgDb package not installed: ", pkg)
  }
  get(pkg, envir = asNamespace(pkg))
}

run_enrich <- function(entrez_ids, ont, orgdb, padj_cutoff) {
  if (length(entrez_ids) == 0) return(NULL)
  tryCatch({
    enrichGO(gene          = entrez_ids,
             OrgDb         = orgdb,
             ont           = ont,
             pAdjustMethod = "BH",
             pvalueCutoff  = padj_cutoff,
             qvalueCutoff  = 1,
             readable      = TRUE)
  }, error = function(e) {
    warning(sprintf("enrichGO failed for ontology %s: %s", ont, e$message))
    NULL
  })
}

main <- function() {
  option_list <- list(
    make_option("--input-csv", type = "character", help = "Merged DEG CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--comparison-label", type = "character", default = NULL),
    make_option("--ontologies", type = "character", default = "BP,MF,CC"),
    make_option("--orgdb", type = "character", default = "org.Mm.eg.db"),
    make_option("--p-adj-cutoff", type = "double", default = 0.05),
    make_option("--min-logfc", type = "double", default = 0),
    make_option("--simplify-cutoff", type = "double", default = 0.7),
    make_option("--min-gene-count", type = "integer", default = 5)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "by_cluster_ora_wang_similarity.R --input-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-csv`, args$`output-dir`)))) {
    stop("--input-csv and --output-dir are required")
  }
  if (!file.exists(args$`input-csv`)) stop("Input CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  comparison_label <- args$`comparison-label`
  if (is.null(comparison_label) || !nzchar(comparison_label)) {
    comparison_label <- tools::file_path_sans_ext(basename(args$`input-csv`))
  }

  ontologies <- parse_vector(args$ontologies)
  if (!length(ontologies)) stop("--ontologies produced zero entries")

  df <- read_csv(args$`input-csv`, show_col_types = FALSE)
  required_cols <- c("gene", "avg_log2FC", "cluster", "p_val_adj")
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  deg_filt <- df %>% filter(p_val_adj <= args$`p-adj-cutoff`, abs(avg_log2FC) >= args$`min-logfc`)
  if (!nrow(deg_filt)) stop("No genes remain after filtering")

  clusters <- sort(unique(deg_filt$cluster))
  orgdb <- load_orgdb(args$orgdb)

  for (ont in ontologies) {
    message(sprintf("Processing ontology %s", ont))
    raw_results <- list()
    simp_results <- list()

    for (cluster_id in clusters) {
      for (direction in c("positive", "negative")) {
        deg_sub <- deg_filt %>%
          filter(cluster == cluster_id,
                 if (direction == "positive") avg_log2FC > 0 else avg_log2FC < 0)
        if (!nrow(deg_sub)) next

        entrez_ids <- bitr(deg_sub$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)$ENTREZID
        entrez_ids <- unique(na.omit(entrez_ids))
        if (length(entrez_ids) < args$`min-gene-count`) {
          message(sprintf("Cluster %s (%s): insufficient genes (%d)", cluster_id, direction, length(entrez_ids)))
          next
        }

        ego <- run_enrich(entrez_ids, ont, orgdb, args$`p-adj-cutoff`)
        if (is.null(ego)) next
        df_raw <- as.data.frame(ego)
        if (!nrow(df_raw)) next
        df_raw$cluster <- cluster_id
        df_raw$direction <- direction
        df_raw$category <- ont
        df_raw$simplified <- FALSE
        key <- paste(cluster_id, direction, sep = "_")
        raw_results[[key]] <- df_raw

        ego_simp <- tryCatch(
          simplify(ego,
                   cutoff     = args$`simplify-cutoff`,
                   by         = "p.adjust",
                   select_fun = min,
                   measure    = "Wang"),
          error = function(e) {
            warning(sprintf("simplify failed for cluster %s (%s): %s", cluster_id, direction, e$message))
            NULL
          }
        )
        if (!is.null(ego_simp)) {
          df_simp <- as.data.frame(ego_simp)
          if (nrow(df_simp)) {
            df_simp$cluster <- cluster_id
            df_simp$direction <- direction
            df_simp$category <- ont
            df_simp$simplified <- TRUE
            simp_results[[key]] <- df_simp
          }
        }
      }
    }

    combined_raw <- bind_rows(raw_results)
    combined_simp <- bind_rows(simp_results)

    for (direction in c("positive", "negative")) {
      out_raw <- file.path(args$`output-dir`,
                           sprintf("GO_%s_%s_%s_raw.csv", ont, direction, comparison_label))
      write_csv(filter(combined_raw, direction == !!direction), out_raw)

      out_simp <- file.path(args$`output-dir`,
                            sprintf("GO_%s_%s_%s_simplified.csv", ont, direction, comparison_label))
      write_csv(filter(combined_simp, direction == !!direction), out_simp)

      message(sprintf("Saved %s results for %s (%s)", ont, comparison_label, direction))
    }
  }

  message("ORA processing complete")
}

if (sys.nframe() == 0L) {
  main()
}
