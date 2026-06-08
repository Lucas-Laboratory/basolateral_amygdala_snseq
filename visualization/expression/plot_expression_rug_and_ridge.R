#!/usr/bin/env Rscript

# Generate rug and ridge plots for ORA gene sets against an indexed gene list.
#
# Inputs
#   --ora-dir: directory containing ORA CSV files with columns `Description`, `geneID`, `p.adjust` (required).
#   --sort-list-csv: CSV with columns `gene`, `index` (required).
#   --output-dir: directory for outputs (required).
#   --p-adjust-threshold: adjusted p-value cutoff (default: 0.001).
#   --width/--height: plot dimensions (default: 6 x 2).
#
# Dependencies: optparse, readr, dplyr, ggplot2, ggridges, circlize

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(circlize)
})

main <- function() {
  option_list <- list(
    make_option("--ora-dir", type = "character", help = "Directory of ORA CSVs"),
    make_option("--sort-list-csv", type = "character", help = "Gene index CSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--p-adjust-threshold", type = "double", default = 0.001),
    make_option("--width", type = "double", default = 6),
    make_option("--height", type = "double", default = 2)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_expression_rug_and_ridge.R --ora-dir DIR --sort-list-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`ora-dir`, args$`sort-list-csv`, args$`output-dir`)))) stop("All inputs required")
  if (!dir.exists(args$`ora-dir`)) stop("ORA directory not found")
  if (!file.exists(args$`sort-list-csv`)) stop("Sort list CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  sort_list <- read_csv(args$`sort-list-csv`, show_col_types = FALSE)
  max_index <- max(sort_list$index, na.rm = TRUE)
  ora_files <- list.files(args$`ora-dir`, pattern = "\\.csv$", full.names = TRUE)

  for (ora_file in ora_files) {
    ora_df <- read_csv(ora_file, show_col_types = FALSE) %>% filter(p.adjust < args$`p-adjust-threshold`)
    if (!nrow(ora_df)) next

    filename <- basename(ora_file)
    ontology <- if (grepl("GO_BP", filename)) "BP" else if (grepl("GO_CC", filename)) "CC" else if (grepl("GO_MF", filename)) "MF" else tools::file_path_sans_ext(filename)

    ontology_dir <- file.path(args$`output-dir`, ontology)
    rug_dir <- file.path(ontology_dir, "rug")
    ridge_dir <- file.path(ontology_dir, "ridge")
    dir.create(rug_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ridge_dir, recursive = TRUE, showWarnings = FALSE)

    ranks <- rank(ora_df$p.adjust, ties.method = "first")
    codes <- sprintf("%04d", ranks)

    out_df <- sort_list
    desc_cols <- make.names(ora_df$Description, unique = TRUE)
    for (i in seq_along(desc_cols)) {
      genes_term <- strsplit(ora_df$geneID[i], "/")[[1]]
      out_df[[desc_cols[i]]] <- as.integer(out_df$gene %in% genes_term)
    }
    write_csv(out_df, file.path(ontology_dir, paste0(ontology, "_binary_matrix.csv")))

    term_colors <- circlize::rand_color(length(desc_cols))
    names(term_colors) <- desc_cols

    for (i in seq_along(desc_cols)) {
      col_name <- desc_cols[i]
      code <- codes[i]
      description <- ora_df$Description[i]
      plot_data <- out_df %>% filter(.data[[col_name]] == 1)
      if (!nrow(plot_data)) next

      p_rug <- ggplot(plot_data, aes(x = index)) +
        geom_rug(sides = "b", colour = term_colors[col_name]) +
        scale_x_continuous(limits = c(1, max_index), expand = c(0, 0)) +
        labs(title = description, x = "Index", y = NULL) +
        theme_minimal() +
        theme(panel.grid = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
      ggsave(file.path(rug_dir, paste0(code, "_", ontology, "_", col_name, "_rug.pdf")), p_rug,
             width = args$width, height = args$height)

      p_ridge <- ggplot(plot_data, aes(x = index, y = 1)) +
        geom_density_ridges(fill = term_colors[col_name], colour = NA, scale = 1) +
        scale_x_continuous(limits = c(1, max_index), expand = c(0, 0)) +
        labs(title = description, x = "Index", y = NULL) +
        theme_minimal() +
        theme(panel.grid = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
      ggsave(file.path(ridge_dir, paste0(code, "_", ontology, "_", col_name, "_ridge.pdf")), p_ridge,
             width = args$width, height = args$height)
    }
  }
}

if (sys.nframe() == 0L) {
  main()
}
