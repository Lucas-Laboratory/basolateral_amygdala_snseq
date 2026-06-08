#!/usr/bin/env Rscript

# Retrieve STRING interactions for a gene set and export a weighted network.
#
# Inputs
#   --genes-csv: CSV containing a `gene` column (required).
#   --output-pdf: network visualisation written via igraph (required).
#   --weights-csv: adjacency matrix of STRING combined scores (required).
#   --species: STRING taxonomy id (default: 10090 for mouse).
#   --string-version: STRING database version (default: 11).
#   --score-threshold: combined_score threshold to display edges (default: 700).
#
# Output
#   A PDF network plot and an adjacency matrix CSV.
#
# Dependencies: optparse, STRINGdb, igraph, readr

suppressPackageStartupMessages({
  library(optparse)
  library(STRINGdb)
  library(igraph)
  library(readr)
})

main <- function() {
  option_list <- list(
    make_option("--genes-csv", type = "character", help = "CSV with gene column"),
    make_option("--output-pdf", type = "character", help = "Destination network PDF"),
    make_option("--weights-csv", type = "character", help = "Destination adjacency CSV"),
    make_option("--species", type = "integer", default = 10090,
                help = "STRING species id [default %default]"),
    make_option("--string-version", type = "character", default = "11",
                help = "STRING database version [default %default]"),
    make_option("--score-threshold", type = "integer", default = 700,
                help = "Combined score threshold for plotting [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "step02_string_network_plot.R --genes-csv FILE --output-pdf FILE --weights-csv FILE [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`genes-csv`, args$`output-pdf`, args$`weights-csv`)))) {
    stop("--genes-csv, --output-pdf, and --weights-csv are required")
  }
  if (!file.exists(args$`genes-csv`)) stop("Gene list CSV not found: ", args$`genes-csv`)

  genes_df <- read_csv(args$`genes-csv`, show_col_types = FALSE)
  if (!"gene" %in% names(genes_df)) stop("CSV must include a 'gene' column")

  string_db <- STRINGdb$new(version = args$`string-version`, species = args$species, score_threshold = 0)
  mapped <- string_db$map(data.frame(gene = genes_df$gene), "gene", removeUnmappedRows = TRUE)
  if (!nrow(mapped)) stop("No genes mapped to STRING ids")

  interactions <- string_db$get_interactions(mapped$STRING_id)
  if (!nrow(interactions)) stop("No STRING interactions retrieved")

  filtered <- subset(interactions, from %in% mapped$STRING_id & to %in% mapped$STRING_id)
  graph_edges <- data.frame(
    from = filtered$from,
    to = filtered$to,
    weight = filtered$combined_score,
    stringsAsFactors = FALSE
  )

  g <- graph_from_data_frame(graph_edges, directed = FALSE)
  V(g)$label <- mapped$gene[match(V(g)$name, mapped$STRING_id)]

  layout_coords <- layout_with_fr(g, niter = 2000, area = vcount(g)^2, repulserad = vcount(g)^3)

  dir.create(dirname(args$`output-pdf`), recursive = TRUE, showWarnings = FALSE)
  pdf(args$`output-pdf`, width = 12, height = 12)
  on.exit(dev.off(), add = TRUE)
  set.seed(42)
  plot(g,
       layout = layout_coords,
       vertex.size = 5,
       vertex.label = NA,
       vertex.color = rgb(173/255, 216/255, 230/255, 0.5),
       edge.color = ifelse(E(g)$weight >= args$`score-threshold`, "#00000066", "#00000022"),
       edge.width = pmax(E(g)$weight / 200, 0.2),
       edge.curved = 0)

  adj_mat <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  dir.create(dirname(args$`weights-csv`), recursive = TRUE, showWarnings = FALSE)
  write.csv(adj_mat, file = args$`weights-csv`, row.names = TRUE)
}

if (sys.nframe() == 0L) {
  main()
}
