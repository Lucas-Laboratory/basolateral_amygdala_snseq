#!/usr/bin/env Rscript

# Perform GO over-representation analysis with optional Wang similarity simplification.
#
# Inputs
#   --genes-csv: CSV with a `gene` column containing symbols (required).
#   --output-dir: directory for enrichment tables (required).
#   --organism-db: OrgDb package name for mapping (default: `org.Mm.eg.db`).
#   --ontologies: comma-separated list of ontologies (default: `BP,MF,CC`).
#   --pvalue-cutoff: p-value cutoff for enrichGO (default: 0.05).
#   --qvalue-cutoff: q-value cutoff for enrichGO (default: 0.2).
#   --simplify-cutoff: semantic similarity cutoff for `simplify` (default: 0.7).
#
# Output
#   For each ontology, writes `<ontology>_enrichment.csv` and `<ontology>_enrichment_simplified.csv`.
#
# Dependencies: optparse, clusterProfiler, GOSemSim, readr

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(GOSemSim)
  library(readr)
})

main <- function() {
  option_list <- list(
    make_option("--genes-csv", type = "character", help = "CSV containing gene symbols"),
    make_option("--output-dir", type = "character", help = "Directory for enrichment outputs"),
    make_option("--organism-db", type = "character", default = "org.Mm.eg.db",
                help = "OrgDb package for gene mapping [default %default]"),
    make_option("--ontologies", type = "character", default = "BP,MF,CC",
                help = "Comma-separated GO ontologies [default %default]"),
    make_option("--pvalue-cutoff", type = "double", default = 0.05,
                help = "p-value cutoff for enrichGO [default %default]"),
    make_option("--qvalue-cutoff", type = "double", default = 0.2,
                help = "q-value cutoff for enrichGO [default %default]"),
    make_option("--simplify-cutoff", type = "double", default = 0.7,
                help = "Wang similarity cutoff for simplify() [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step05_basic_ora_wang_similarity.R --genes-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`genes-csv`, args$`output-dir`)))) {
    stop("--genes-csv and --output-dir are required")
  }
  if (!file.exists(args$`genes-csv`)) stop("Gene CSV not found: ", args$`genes-csv`)

  if (!requireNamespace(args$`organism-db`, quietly = TRUE)) {
    stop("OrgDb package not installed: ", args$`organism-db`)
  }
  org_db <- get(args$`organism-db`)

  genes <- read_csv(args$`genes-csv`, show_col_types = FALSE)$gene
  entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  entrez_ids <- entrez_map$ENTREZID
  if (!length(entrez_ids)) stop("No genes mapped to Entrez IDs")

  ontologies <- trimws(strsplit(args$ontologies, ",")[[1]])
  ontologies <- ontologies[ontologies != ""]
  if (!length(ontologies)) stop("No ontologies specified")

  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  for (ont in ontologies) {
    enrich_res <- enrichGO(gene = entrez_ids,
                           OrgDb = org_db,
                           ont = ont,
                           pAdjustMethod = "BH",
                           pvalueCutoff = args$`pvalue-cutoff`,
                           qvalueCutoff = args$`qvalue-cutoff`,
                           readable = TRUE)

    enrich_df <- as.data.frame(enrich_res)
    write.csv(enrich_df,
              file = file.path(args$`output-dir`, paste0("GO_", ont, "_enrichment.csv")),
              row.names = FALSE)

    if (!is.null(enrich_res) && nrow(enrich_df)) {
      simplified <- simplify(enrich_res,
                             cutoff = args$`simplify-cutoff`,
                             by = "p.adjust",
                             select_fun = min,
                             measure = "Wang")
      write.csv(as.data.frame(simplified),
                file = file.path(args$`output-dir`, paste0("GO_", ont, "_enrichment_simplified.csv")),
                row.names = FALSE)
    }
  }
}

if (sys.nframe() == 0L) {
  main()
}
