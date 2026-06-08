#!/usr/bin/env Rscript

# Create a colour-coded STRING network using GO molecular function annotations.
#
# Inputs
#   --weights-csv: adjacency matrix of STRING combined scores (required).
#   --embedding-csv: CSV with columns `position` and STRING identifiers (required).
#   --symbol-csv: CSV mapping positions to gene symbols (required).
#   --gomf-csv: long-format GO MF table containing `gene`, `gomf_description`, `hex_color` (required).
#   --output-pdf: path to save the network plot (required).
#   --default-color: fallback hex colour for nodes without annotations (default: `#808080`).
#   --alpha: node colour transparency (default: 0.25).
#
# Output
#   PDF of the coloured network.
#
# Dependencies: optparse, igraph, readr, dplyr, grDevices

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
  library(readr)
  library(dplyr)
  library(grDevices)
})

build_network <- function(weights_path, embed_path, symbol_path, gomf_path, output_pdf, default_color, alpha) {
  weights <- as.matrix(read.csv(weights_path, row.names = 1, check.names = FALSE))
  embed <- read_csv(embed_path, show_col_types = FALSE)
  symbols <- read_csv(symbol_path, show_col_types = FALSE)
  gomf <- read_csv(gomf_path, show_col_types = FALSE)

  if (!all(c("position", "gene") %in% names(embed))) stop("Embedding CSV must contain 'position' and 'gene' columns")
  if (!all(c("position", "gene_symbol") %in% names(symbols))) {
    if (!all(c("position", "gene") %in% names(symbols))) {
      stop("Symbol CSV must contain columns 'position' and either 'gene_symbol' or 'gene'")
    }
    symbols <- symbols %>% rename(gene_symbol = gene)
  }
  if (!all(c("gene", "hex_color", "gomf_description") %in% names(gomf))) {
    stop("GO MF CSV must include 'gene', 'hex_color', and 'gomf_description'")
  }

  mapping <- embed %>%
    select(position, string_id = gene) %>%
    inner_join(symbols %>% select(position, gene_symbol), by = "position")

  valid_ids <- mapping$string_id
  sub_mat <- weights[rownames(weights) %in% valid_ids, colnames(weights) %in% valid_ids, drop = FALSE]
  if (!nrow(sub_mat) || !ncol(sub_mat)) stop("No overlapping STRING identifiers between weights and mapping")

  rownames(sub_mat) <- mapping$gene_symbol[match(rownames(sub_mat), mapping$string_id)]
  colnames(sub_mat) <- mapping$gene_symbol[match(colnames(sub_mat), mapping$string_id)]

  g <- graph_from_adjacency_matrix(sub_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  layout_coords <- layout_with_fr(g, weights = E(g)$weight)

  gomf_unique <- gomf %>% distinct(gene, .keep_all = TRUE)
  vertex_colors <- sapply(V(g)$name, function(gene) {
    col <- gomf_unique$hex_color[match(gene, gomf_unique$gene)]
    if (is.na(col) || !nzchar(col)) col <- default_color
    adjustcolor(col, alpha.f = alpha)
  })

  legend_palette <- gomf %>% distinct(gomf_description, hex_color)

  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(output_pdf, width = 8, height = 8, bg = "white")
  on.exit(dev.off(), add = TRUE)

  plot(g,
       layout = layout_coords,
       vertex.size = 5,
       vertex.color = vertex_colors,
       vertex.label = NA,
       edge.curved = 0,
       edge.color = adjustcolor("darkgray", alpha.f = 0.25),
       edge.width = 1)

  if (nrow(legend_palette)) {
    legend("topright",
           legend = legend_palette$gomf_description,
           pch = 21,
           pt.bg = adjustcolor(legend_palette$hex_color, alpha.f = alpha),
           pt.cex = 1.5,
           title = "GOMF")
  }
}

main <- function() {
  option_list <- list(
    make_option("--weights-csv", type = "character", help = "STRING adjacency matrix CSV"),
    make_option("--embedding-csv", type = "character", help = "MDS embedding CSV"),
    make_option("--symbol-csv", type = "character", help = "Symbol mapping CSV"),
    make_option("--gomf-csv", type = "character", help = "GO MF annotation CSV"),
    make_option("--output-pdf", type = "character", help = "Destination PDF"),
    make_option("--default-color", type = "character", default = "#808080",
                help = "Fallback node colour [default %default]"),
    make_option("--alpha", type = "double", default = 0.25,
                help = "Node colour transparency [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "make_network_graph_colorized.R --weights-csv FILE --embedding-csv FILE --symbol-csv FILE --gomf-csv FILE --output-pdf FILE [options]")
  args <- parse_args(parser)

  required <- c(args$`weights-csv`, args$`embedding-csv`, args$`symbol-csv`, args$`gomf-csv`, args$`output-pdf`)
  if (any(!nzchar(required))) stop("All input paths and --output-pdf are required")
  if (!file.exists(args$`weights-csv`)) stop("Weights CSV not found: ", args$`weights-csv`)
  if (!file.exists(args$`embedding-csv`)) stop("Embedding CSV not found: ", args$`embedding-csv`)
  if (!file.exists(args$`symbol-csv`)) stop("Symbol CSV not found: ", args$`symbol-csv`)
  if (!file.exists(args$`gomf-csv`)) stop("GO MF CSV not found: ", args$`gomf-csv`)

  build_network(args$`weights-csv`, args$`embedding-csv`, args$`symbol-csv`, args$`gomf-csv`,
                args$`output-pdf`, args$`default-color`, args$`alpha`)
}

if (sys.nframe() == 0L) {
  main()
}
