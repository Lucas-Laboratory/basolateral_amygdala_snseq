#!/usr/bin/env Rscript

# Convert STRING identifiers (with species prefixes) to preferred gene symbols.
#
# Inputs
#   --input-csv: CSV produced by MDS step containing columns `gene` and `position` (required).
#   --output-csv: destination CSV with gene symbols and positions (required).
#   --species: STRING taxonomy id (default: 10090 for mouse).
#   --string-version: STRING database version (default: 11).
#   --id-column: column holding STRING ids or prefixed peptide ids (default: `gene`).
#
# Output
#   CSV with columns `gene_symbol` and `position`.
#
# Dependencies: optparse, STRINGdb, readr, dplyr

suppressPackageStartupMessages({
  library(optparse)
  library(STRINGdb)
  library(readr)
  library(dplyr)
})

sentence_case <- function(x) {
  ifelse(is.na(x) | x == "", NA_character_, paste0(toupper(substr(tolower(x), 1, 1)), substr(tolower(x), 2, nchar(x))))
}

main <- function() {
  option_list <- list(
    make_option("--input-csv", type = "character", help = "Input CSV with STRING ids"),
    make_option("--output-csv", type = "character", help = "Destination CSV"),
    make_option("--species", type = "integer", default = 10090,
                help = "STRING species id [default %default]"),
    make_option("--string-version", type = "character", default = "11",
                help = "STRING database version [default %default]"),
    make_option("--id-column", type = "character", default = "gene",
                help = "Column containing STRING identifiers [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step04_convert_string_id_to_gene_symbol.R --input-csv FILE --output-csv FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-csv`, args$`output-csv`)))) {
    stop("--input-csv and --output-csv are required")
  }
  if (!file.exists(args$`input-csv`)) stop("Input CSV not found: ", args$`input-csv`)

  df <- read_csv(args$`input-csv`, show_col_types = FALSE)
  if (!args$`id-column` %in% names(df)) {
    stop("Column not found: ", args$`id-column`)
  }

  string_db <- STRINGdb$new(version = args$`string-version`, species = args$species, score_threshold = 0)

  df$peptide_id <- sub("^[0-9]+\\.", "", df[[args$`id-column`]])
  mapped <- string_db$map(df, "peptide_id", removeUnmappedRows = FALSE)
  annotated <- string_db$add_proteins_description(mapped)
  annotated$gene_symbol <- sentence_case(annotated$preferred_name)

  result <- annotated %>%
    transmute(gene_symbol = gene_symbol, position = position)

  dir.create(dirname(args$`output-csv`), recursive = TRUE, showWarnings = FALSE)
  write_csv(result, args$`output-csv`)
}

if (sys.nframe() == 0L) {
  main()
}
