#!/usr/bin/env Rscript

# Annotate differential expression result tables with NCBI Gene descriptions.
#
# Inputs
#   --input-dir: directory containing CSV files with a gene symbol column (required).
#   --output-dir: destination directory for annotated CSV files (default: `<input-dir>/annotated`).
#   --gene-column: column name holding gene symbols (default: `gene`).
#   --organism-db: Bioconductor organism package providing SYMBOL→ENTREZID mappings (default: `org.Mm.eg.db`).
#   --email: optional email address passed to NCBI (recommended for rentrez).
#   --chunk-size: number of Entrez IDs to request per API call (default: 300).
#
# Output
#   Annotated CSVs with `entrez_id`, `brief_desc`, and `long_desc` columns appended.
#
# Dependencies: optparse, AnnotationDbi, rentrez, organism-specific `org.*.eg.db`

suppressPackageStartupMessages({
  library(optparse)
  library(AnnotationDbi)
  library(rentrez)
})

fetch_ncbi_summaries <- function(entrez_ids, chunk_size) {
  valid_ids <- unique(na.omit(entrez_ids))
  if (!length(valid_ids)) {
    return(data.frame(entrez_id = character(), description = character(), summary = character(),
                      stringsAsFactors = FALSE))
  }

  chunks <- split(valid_ids, ceiling(seq_along(valid_ids) / chunk_size))
  pieces <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    message("Requesting annotations for chunk ", i, "/", length(chunks), " (", length(ids), " IDs)")
    api_results <- rentrez::entrez_summary(db = "gene", id = ids)
    pieces[[i]] <- do.call(rbind, lapply(names(api_results), function(key) {
      entry <- api_results[[key]]
      data.frame(
        entrez_id = key,
        description = entry$description,
        summary = entry$summary,
        stringsAsFactors = FALSE
      )
    }))
    Sys.sleep(0.34)
  }

  do.call(rbind, pieces)
}

annotate_tables <- function(input_dir, output_dir, gene_column, organism_db, chunk_size) {
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) {
    stop("No CSV files found in ", input_dir)
  }

  if (!requireNamespace(organism_db, quietly = TRUE)) {
    stop("Bioconductor package not installed: ", organism_db)
  }
  org_db <- get(organism_db)

  all_genes <- unique(unlist(lapply(csv_files, function(path) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE)
    if (!gene_column %in% names(df)) {
      stop("Column '", gene_column, "' not found in ", basename(path))
    }
    df[[gene_column]]
  })))

  message("Mapping ", length(all_genes), " gene symbols to Entrez IDs")
  entrez_ids <- AnnotationDbi::mapIds(
    x = org_db,
    keys = all_genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  summaries <- fetch_ncbi_summaries(unname(entrez_ids), chunk_size)

  for (file in csv_files) {
    df <- utils::read.csv(file, stringsAsFactors = FALSE)
    if (!gene_column %in% names(df)) {
      stop("Column '", gene_column, "' not found in ", basename(file))
    }

    gene_entrez <- unname(entrez_ids[df[[gene_column]]])
    desc_idx <- match(gene_entrez, summaries$entrez_id)

    df$entrez_id <- gene_entrez
    df$brief_desc <- summaries$description[desc_idx]
    df$long_desc <- summaries$summary[desc_idx]

    output_path <- file.path(output_dir, paste0("annotated_", basename(file)))
    utils::write.csv(df, output_path, row.names = FALSE)
    message("Wrote ", basename(output_path))
  }
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory containing DEG CSV files"),
    make_option("--output-dir", type = "character", default = NULL,
                help = "Directory for annotated CSV files [default: <input-dir>/annotated]"),
    make_option("--gene-column", type = "character", default = "gene",
                help = "Column containing gene symbols [default %default]"),
    make_option("--organism-db", type = "character", default = "org.Mm.eg.db",
                help = "Organism annotation package providing SYMBOL→ENTREZID mapping"),
    make_option("--email", type = "character", default = NULL,
                help = "Email address for NCBI API etiquette"),
    make_option("--chunk-size", type = "integer", default = 300,
                help = "Number of Entrez IDs per rentrez request [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "annotate_deg.R --input-dir DIR [options]")
  args <- parse_args(parser)

  input_dir <- args$`input-dir`
  output_dir <- args$`output-dir`
  gene_column <- args$`gene-column`
  organism_db <- args$`organism-db`
  email <- args$email
  chunk_size <- args$`chunk-size`

  if (is.null(input_dir) || !nzchar(input_dir)) stop("--input-dir is required")
  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- file.path(input_dir, "annotated")
  }

  if (!is.null(email)) {
    options(rentrez.email = email)
  }

  annotate_tables(input_dir, output_dir, gene_column, organism_db, chunk_size)
}

if (sys.nframe() == 0L) {
  main()
}
