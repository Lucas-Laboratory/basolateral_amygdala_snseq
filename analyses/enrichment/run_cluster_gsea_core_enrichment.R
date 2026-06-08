#!/usr/bin/env Rscript

# Perform cluster-wise GO GSEA across BP/MF/CC ontologies using clusterProfiler.
#
# Inputs
#   --deg-csv: CSV containing differential expression results with columns `gene`, `avg_log2FC`, `cluster`, `pct.1`, `pct.2` (required).
#   --output-dir: directory where GSEA tables will be written (required).
#   --organism-db: Bioconductor OrgDb package to use for ID conversion (default: `org.Mm.eg.db`).
#   --pct-threshold: remove genes expressed in fewer than this fraction in both groups (default: 0.1).
#   --min-gs-size: minimum gene set size for GSEA (default: 10).
#   --max-gs-size: maximum gene set size for GSEA (default: 500).
#   --pvalue-cutoff: p-value cutoff passed to `gseGO` (default: 0.05).
#   --ontologies: comma-separated ontologies to evaluate, e.g. `BP,MF,CC` (default: all three).
#
# Output
#   Per-ontology CSVs capturing the cluster-wise GSEA results with readable gene symbols.
#
# Dependencies: optparse, clusterProfiler, GOSemSim, enrichplot, OrgDb package

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(GOSemSim)
  library(enrichplot)
})

run_cluster_gsea <- function(deg_table, org_db, ontologies, pct_threshold, min_gs_size, max_gs_size, pvalue_cutoff, output_dir) {
  if (!all(c("gene", "avg_log2FC", "cluster", "pct.1", "pct.2") %in% names(deg_table))) {
    stop("Input DEG table must include columns: gene, avg_log2FC, cluster, pct.1, pct.2")
  }

  deg_filtered <- subset(deg_table, !(pct.1 < pct_threshold & pct.2 < pct_threshold))
  clusters <- unique(deg_filtered$cluster)

  for (ont in ontologies) {
    message("Running GO GSEA for ontology ", ont)
    ontology_results <- lapply(clusters, function(cluster_id) {
      sub_df <- subset(deg_filtered, cluster == cluster_id)
      if (!nrow(sub_df)) return(NULL)

      gene_map <- suppressMessages(clusterProfiler::bitr(sub_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db))
      ranked <- merge(sub_df, gene_map, by.x = "gene", by.y = "SYMBOL")
      if (!nrow(ranked)) return(NULL)

      ranking <- ranked$avg_log2FC
      names(ranking) <- ranked$ENTREZID
      ranking <- ranking + runif(length(ranking), -1e-6, 1e-6)
      ranking <- sort(ranking, decreasing = TRUE)

      gsea <- tryCatch(
        gseGO(
          geneList = ranking,
          OrgDb = org_db,
          ont = ont,
          keyType = "ENTREZID",
          minGSSize = min_gs_size,
          maxGSSize = max_gs_size,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE
        ),
        error = function(err) {
          warning("gseGO failed for cluster ", cluster_id, " (", err$message, ")")
          NULL
        }
      )
      if (is.null(gsea)) return(NULL)

      gsea <- setReadable(gsea, OrgDb = org_db)
      result <- as.data.frame(gsea@result)
      if (!nrow(result)) return(NULL)

      result$cluster <- cluster_id
      result$category <- ont
      result
    })

    ontology_results <- do.call(rbind, ontology_results[!vapply(ontology_results, is.null, logical(1))])
    if (is.null(ontology_results) || !nrow(ontology_results)) {
      message("No significant terms for ontology ", ont)
      next
    }

    output_path <- file.path(output_dir, paste0("gsea_", tolower(ont), "_by_cluster.csv"))
    utils::write.csv(ontology_results, output_path, row.names = FALSE)
    message("Saved ontology ", ont, " results to ", output_path)
  }
}

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "Path to DEG CSV"),
    make_option("--output-dir", type = "character", help = "Directory for output tables"),
    make_option("--organism-db", type = "character", default = "org.Mm.eg.db",
                help = "OrgDb package for ID conversion [default %default]"),
    make_option("--pct-threshold", type = "double", default = 0.1,
                help = "Minimum detection fraction retained per group [default %default]"),
    make_option("--min-gs-size", type = "integer", default = 10,
                help = "Minimum genes per set for GSEA [default %default]"),
    make_option("--max-gs-size", type = "integer", default = 500,
                help = "Maximum genes per set for GSEA [default %default]"),
    make_option("--pvalue-cutoff", type = "double", default = 0.05,
                help = "P-value cutoff for gseGO [default %default]"),
    make_option("--ontologies", type = "character", default = "BP,MF,CC",
                help = "Comma-separated GO ontologies [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "run_cluster_gsea_core_enrichment.R --deg-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (is.null(args$`deg-csv`) || !nzchar(args$`deg-csv`)) stop("--deg-csv is required")
  if (is.null(args$`output-dir`) || !nzchar(args$`output-dir`)) stop("--output-dir is required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found: ", args$`deg-csv`)

  if (!requireNamespace(args$`organism-db`, quietly = TRUE)) {
    stop("OrgDb package not installed: ", args$`organism-db`)
  }
  org_db <- get(args$`organism-db`)

  deg_table <- utils::read.csv(args$`deg-csv`, stringsAsFactors = FALSE)

  ontologies <- strsplit(args$ontologies, ",", fixed = TRUE)[[1]]
  ontologies <- trimws(ontologies)
  ontologies <- ontologies[ontologies != ""]
  if (!length(ontologies)) stop("No ontologies provided")

  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  run_cluster_gsea(
    deg_table = deg_table,
    org_db = org_db,
    ontologies = ontologies,
    pct_threshold = args$`pct-threshold`,
    min_gs_size = args$`min-gs-size`,
    max_gs_size = args$`max-gs-size`,
    pvalue_cutoff = args$`pvalue-cutoff`,
    output_dir = args$`output-dir`
  )
}

if (sys.nframe() == 0L) {
  main()
}
