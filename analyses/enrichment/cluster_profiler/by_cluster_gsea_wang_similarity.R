#!/usr/bin/env Rscript

# Run cluster-wise GO GSEA using clusterProfiler and export combined results for
# each ontology. Genes are pre-filtered by minimum detection percentage and the
# ranked list is derived from average log2 fold-change values.
#
# Inputs
#   --input-csv: merged DEG CSV with columns `gene`, `avg_log2FC`, `cluster`,
#                `pct.1`, `pct.2` (required).
#   --output-dir: directory where ontology-level CSVs are written (required).
#   --comparison-label: label embedded in output filenames (default: derived from
#                       input filename).
#   --ontologies: comma-separated GO ontologies to evaluate (default: BP,MF,CC).
#   --orgdb: name of the Bioconductor OrgDb package (default: org.Mm.eg.db).
#   --min-pct: minimum detection fraction in either group to retain a gene (default: 0.1).
#   --min-gs-size / --max-gs-size: gene-set size bounds for GSEA (defaults: 10 / 500).
#   --pvalue-cutoff: p-value cutoff passed to `gseGO` (default: 0.05).
#   --jitter: uniform noise range added to rankings to break ties (default: 1e-6).
#
# Outputs
#   - CSV per ontology combining all cluster results.
#
# Dependencies: optparse, readr, dplyr, clusterProfiler

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(clusterProfiler)
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

run_gsea <- function(ranked_genes, ont, orgdb, min_gs, max_gs, cutoff) {
  if (!length(ranked_genes)) return(NULL)
  tryCatch({
    gseGO(geneList     = ranked_genes,
          OrgDb        = orgdb,
          ont          = ont,
          keyType      = "ENTREZID",
          minGSSize    = min_gs,
          maxGSSize    = max_gs,
          pvalueCutoff = cutoff,
          verbose      = FALSE)
  }, error = function(e) {
    warning(sprintf("gseGO failed for ontology %s: %s", ont, e$message))
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
    make_option("--min-pct", type = "double", default = 0.1),
    make_option("--min-gs-size", type = "integer", default = 10),
    make_option("--max-gs-size", type = "integer", default = 500),
    make_option("--pvalue-cutoff", type = "double", default = 0.05),
    make_option("--jitter", type = "double", default = 1e-6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "by_cluster_gsea_wang_similarity.R --input-csv FILE --output-dir DIR [options]")
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
  required_cols <- c("gene", "avg_log2FC", "cluster", "pct.1", "pct.2")
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  deg_filt <- df %>% filter(!(pct.1 < args$`min-pct` & pct.2 < args$`min-pct`))
  if (!nrow(deg_filt)) stop("No genes remain after min-pct filtering")

  clusters <- sort(unique(deg_filt$cluster))
  orgdb <- load_orgdb(args$orgdb)

  for (ont in ontologies) {
    message(sprintf("Processing ontology %s", ont))
    results_list <- vector("list", length(clusters))
    names(results_list) <- clusters

    for (cluster_id in clusters) {
      deg_sub <- deg_filt %>% filter(cluster == cluster_id)
      if (!nrow(deg_sub)) next

      gene_map <- tryCatch(
        bitr(deg_sub$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb),
        error = function(e) NULL
      )
      if (is.null(gene_map) || !nrow(gene_map)) {
        message(sprintf("Cluster %s: no Entrez IDs mapped", cluster_id))
        next
      }
      ranked <- deg_sub %>% inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
        transmute(ENTREZID, score = avg_log2FC + runif(n(), -args$jitter, args$jitter)) %>%
        arrange(desc(score))
      if (!nrow(ranked)) next
      gene_list <- ranked$score
      names(gene_list) <- ranked$ENTREZID

      gsea_res <- run_gsea(gene_list, ont, orgdb, args$`min-gs-size`, args$`max-gs-size`, args$`pvalue-cutoff`)
      if (is.null(gsea_res)) next

      gsea_res <- setReadable(gsea_res, OrgDb = orgdb)
      df_res <- as.data.frame(gsea_res@result)
      if (!nrow(df_res)) next
      df_res$cluster <- cluster_id
      df_res$category <- ont
      results_list[[as.character(cluster_id)]] <- df_res
      message(sprintf("Cluster %s: %d terms", cluster_id, nrow(df_res)))
    }

    combined <- bind_rows(results_list)
    if (!nrow(combined)) {
      message(sprintf("Ontology %s returned no results", ont))
      next
    }

    out_file <- file.path(args$`output-dir`, sprintf("GSEA_%s_%s.csv", ont, comparison_label))
    write_csv(combined, out_file)
    message(sprintf("Saved %s", out_file))
  }

  message("GSEA processing complete")
}

if (sys.nframe() == 0L) {
  main()
}
