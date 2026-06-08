#!/usr/bin/env Rscript

# Produce bar/dot and optional UpSet plots for GO term counts across clusters.
#
# Inputs
#   --input-dir: directory containing GSEA result CSVs (required).
#   --output-dir: directory for PDFs (required).
#   --file-pattern: regex to match files (default: `^GSEA_(MF|BP|CC)_reint_full_mod_DEG_.*\.csv$`).
#   --alpha: adjusted p-value threshold (default: 0.01).
#   --top-n: comma-separated list of term counts to display (default: `10,20`).
#   --clusters: optional cluster list/range for dot plots (default: none).
#   --make-upset: include UpSet plots (default: TRUE).
#   --width/--height: base PDF size (default: 12 x 7).
#
# Dependencies: optparse, ggplot2, dplyr, tidyr, stringr, patchwork

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

parse_name <- function(path) {
  nm <- basename(path)
  m <- regmatches(nm, regexec("^GSEA_(MF|BP|CC)_reint_full_mod_DEG_([^.]+)\\.csv$", nm))[[1]]
  if (length(m) >= 3) list(ontology = m[2], comparison = m[3]) else list(ontology = "UNK", comparison = tools::file_path_sans_ext(nm))
}

pick_column <- function(df, candidates) {
  ix <- which(tolower(names(df)) %in% tolower(candidates))
  if (length(ix)) names(df)[ix[1]] else NA_character_
}

pick_term_column <- function(df) {
  term_cols <- c("Description", "Term", "Name", "Pathway")
  id_cols <- c("ID", "GO", "GO.ID")
  term <- pick_column(df, term_cols)
  id <- pick_column(df, id_cols)
  list(term = if (!is.na(term)) term else if (!is.na(id)) id else names(df)[which(sapply(df, is.character))[1]],
       id = if (!is.na(id)) id else NA_character_)
}

pick_sig_col <- function(df) {
  adj <- c("padj", "p_adj", "p.adjust", "qvalue", "fdr")
  raw <- c("pval", "p_val", "p.value", "p")
  col <- pick_column(df, adj)
  if (!is.na(col)) return(col)
  pick_column(df, raw)
}

pick_cluster_col <- function(df) {
  col <- pick_column(df, c("cluster", "cluster_id"))
  if (is.na(col)) stop("No cluster column found in CSV")
  col
}

read_gsea <- function(path, alpha) {
  df <- read.csv(path, check.names = FALSE)
  if (!nrow(df)) stop("Empty file: ", basename(path))
  parsed <- parse_name(path)
  term_cols <- pick_term_column(df)
  sig_col <- pick_sig_col(df)
  cluster_col <- pick_cluster_col(df)
  nes_col <- pick_column(df, c("nes", "normalized_enrichment_score"))
  sig_vec <- if (!is.na(sig_col)) {
    p <- suppressWarnings(as.numeric(df[[sig_col]]))
    !is.na(p) & p <= alpha
  } else if (!is.na(nes_col)) {
    nes <- suppressWarnings(as.numeric(df[[nes_col]]))
    !is.na(nes) & nes != 0
  } else rep(FALSE, nrow(df))
  tibble(
    term = df[[term_cols$term]],
    cluster = suppressWarnings(as.integer(as.character(df[[cluster_col]]))),
    sig = sig_vec,
    NES = if (!is.na(nes_col)) suppressWarnings(as.numeric(df[[nes_col]])) else NA_real_,
    comparison = parsed$comparison,
    ontology = parsed$ontology
  )
}

make_bar_dot <- function(dat, top_n, clusters_to_show, ont) {
  term_counts <- dat %>%
    group_by(term) %>% summarise(n_clusters = sum(sig, na.rm = TRUE), sum_abs_NES = sum(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    filter(n_clusters > 0) %>% arrange(desc(n_clusters), desc(sum_abs_NES)) %>% slice_head(n = top_n)
  if (!nrow(term_counts)) return(list(plot = NULL, counts = term_counts))
  top_terms <- term_counts$term
  clusters_all <- if (!is.null(clusters_to_show)) as.character(clusters_to_show) else dat %>% filter(term %in% top_terms, sig, !is.na(cluster)) %>% pull(cluster) %>% unique() %>% sort()
  dat_top <- dat %>% filter(term %in% top_terms) %>% left_join(term_counts, by = "term") %>%
    mutate(term = reorder(term, n_clusters), cluster = factor(cluster, levels = clusters_all)) %>% filter(!is.na(cluster))

  p_bar <- ggplot(term_counts, aes(x = reorder(term, n_clusters), y = n_clusters)) +
    geom_col() + geom_text(aes(label = n_clusters), hjust = -0.1, size = 3) + coord_flip(clip = "off") +
    labs(title = paste0("GO:", ont, " — # clusters (Top ", nrow(term_counts), ")"), x = NULL, y = "# clusters") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), plot.margin = margin(10, 20, 10, 10))

  p_dot <- ggplot(dat_top, aes(x = cluster, y = term)) +
    geom_point(aes(alpha = sig), size = 2.2) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
    labs(title = "Clusters per term", x = "cluster", y = NULL) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 9))

  list(plot = (p_dot | p_bar) + plot_layout(widths = c(1.25, 1)), counts = term_counts)
}

maybe_upset <- function(dat, top_terms, outfile, upset_by, clusters_fixed) {
  if (!length(top_terms) || !nzchar(outfile)) return(invisible(NULL))
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("UpSetR not installed; skipping UpSet plot")
    return(invisible(NULL))
  }
  dat_top <- dat %>% filter(term %in% top_terms)
  if (identical(upset_by, "clusters")) {
    cl_all <- clusters_fixed
    sets <- dat_top %>% filter(sig, !is.na(cluster)) %>% mutate(cluster = as.integer(cluster)) %>%
      distinct(term, cluster) %>% mutate(val = 1L) %>% pivot_wider(names_from = cluster, values_from = val, values_fill = 0L)
    missing_cols <- setdiff(as.character(cl_all), names(sets))
    for (mc in missing_cols) sets[[mc]] <- 0L
    sets <- sets[, c(setdiff(names(sets), as.character(cl_all)), as.character(cl_all)), drop = FALSE]
    rownames(sets) <- make.names(sets[[1]], unique = TRUE)
    df <- as.data.frame(sets[, as.character(cl_all), drop = FALSE])
    colnames(df) <- paste0("cl_", colnames(df))
    intersections <- lapply(colnames(df), identity)
    pdf(outfile, width = 12, height = 6)
    tryCatch(
      UpSetR::upset(df, sets = colnames(df), intersections = intersections, keep.order = TRUE,
                    empty.intersections = "on", mb.ratio = c(0.6, 0.4), sets.bar.color = "grey30",
                    mainbar.y.label = "terms per cluster", sets.x.label = "term count per cluster"),
      error = function(e) {
        message("Skipping UpSet plot for ", basename(outfile), ": ", e$message)
        plot.new()
        title(main = "UpSet plot skipped")
        text(0.5, 0.5, e$message)
      }
    )
    dev.off()
  } else {
    message("UpSet by terms not implemented in this CLI version")
  }
}

main_script <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Input directory"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--file-pattern", type = "character", default = "^GSEA_(MF|BP|CC)_reint_full_mod_DEG_.*\\.csv$"),
    make_option("--alpha", type = "double", default = 0.01),
    make_option("--top-n", type = "character", default = "10,20"),
    make_option("--clusters", type = "character", default = NULL),
    make_option("--make-upset", action = "store_true", default = TRUE),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 7),
    make_option("--upset-height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_top_go_terms_upset.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) stop("--input-dir and --output-dir are required")
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(args$`input-dir`, pattern = args$`file-pattern`, full.names = TRUE)
  if (!length(files)) stop("No input files matched pattern")

  clusters_to_show <- if (!is.null(args$clusters) && nzchar(args$clusters)) {
    if (grepl(":", args$clusters)) {
      bounds <- as.integer(strsplit(args$clusters, ":")[[1]])
      as.integer(seq(bounds[1], bounds[2]))
    } else {
      as.integer(trimws(strsplit(args$clusters, ",")[[1]]))
    }
  } else NULL
  clusters_fixed <- if (!is.null(clusters_to_show)) clusters_to_show else 0:29
  top_ns <- as.integer(trimws(strsplit(args$`top-n`, ",")[[1]]))

  summary_lines <- character()
  for (path in files) {
    parsed <- parse_name(path)
    dat <- read_gsea(path, args$alpha)
    for (tn in top_ns) {
      res <- make_bar_dot(dat, tn, clusters_to_show, parsed$ontology)
      if (is.null(res$plot)) next
      pdf_file <- file.path(args$`output-dir`, sprintf("GO%s_top%02d_%s.pdf", parsed$ontology, tn, parsed$comparison))
      n_terms <- nrow(res$counts)
      pdf_height <- max(5, min(14, n_terms * 0.45 + 2))
      pdf(file = pdf_file, width = args$width, height = ifelse(isTRUE(TRUE), pdf_height, args$height))
      print(res$plot)
      dev.off()
      if (isTRUE(args$`make-upset`)) {
        upset_file <- file.path(args$`output-dir`, sprintf("GO%s_top%02d_%s_UpSet.pdf", parsed$ontology, tn, parsed$comparison))
        maybe_upset(dat, res$counts$term, upset_file, "clusters", clusters_fixed)
      }
    }
  }
}

if (sys.nframe() == 0L) {
  main_script()
}
