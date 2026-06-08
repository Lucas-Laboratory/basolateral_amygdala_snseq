#!/usr/bin/env Rscript

# Score transcription factor activity per cluster using decoupleR and DoRothEA regulons.
#
# Inputs
#   --deg-csv: differential expression CSV containing `gene`, `avg_log2FC`, `cluster`, `pct.1`, `pct.2` (required).
#   --output-dir: directory where per-cluster TF scores are written (required).
#   --species: DoRothEA species dataset to load (`human` or `mouse`, default: `mouse`).
#   --confidence: comma-separated DoRothEA confidence levels to retain (default: `A,B,C`).
#   --method: decoupleR scoring method (`wmean` or `mlm`, default: `wmean`).
#   --minsize: minimum regulon size (default: 5).
#   --resamples: number of permutations for wmean significance estimation (default: 10000).
#   --pct-threshold: drop genes with both detection fractions below this value (default: 0.1).
#   --exclude-clusters: comma-separated list of clusters to skip (optional).
#   --output-format: `per-cluster`, `long`, or `both` (default: `per-cluster`).
#   --comparison-label: comparison name to embed in long-format output filename (default: comparison).
#
# Output
#   Per-cluster CSVs and/or a long-format CSV with TF activity scores sorted by p-value.
#
# Dependencies: optparse, decoupleR, dorothea, dplyr, tidyr, tibble

suppressPackageStartupMessages({
  library(optparse)
  library(decoupleR)
  library(dorothea)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

load_regulon <- function(species, confidence) {
  data(dorothea_mm, package = "dorothea")
  data(dorothea_hs, package = "dorothea")
  dataset <- switch(tolower(species),
                    mouse = dorothea_mm,
                    human = dorothea_hs,
                    stop("Unsupported species: ", species))
  dataset %>% filter(.data$confidence %in% .env$confidence)
}

main <- function() {
  option_list <- list(
    make_option("--deg-csv", type = "character", help = "Cluster-wise DEG CSV"),
    make_option("--output-dir", type = "character", help = "Directory for TF score tables"),
    make_option("--species", type = "character", default = "mouse",
                help = "Species dataset for DoRothEA (mouse/human) [default %default]"),
    make_option("--confidence", type = "character", default = "A,B,C",
                help = "Comma-separated confidence levels [default %default]"),
    make_option("--method", type = "character", default = "wmean",
                help = "decoupleR method: wmean or mlm [default %default]"),
    make_option("--minsize", type = "integer", default = 5,
                help = "Minimum regulon size [default %default]"),
    make_option("--resamples", type = "integer", default = 10000,
                help = "Number of permutations for run_wmean [default %default]"),
    make_option("--seed", type = "integer", default = 42,
                help = "Random seed for run_wmean permutations [default %default]"),
    make_option("--pct-threshold", type = "double", default = 0.1,
                help = "Drop genes expressed below this pct in both groups [default %default]"),
    make_option("--exclude-clusters", type = "character", default = NULL,
                help = "Comma-separated cluster IDs to skip"),
    make_option("--output-format", type = "character", default = "per-cluster",
                help = "Output format: per-cluster, long, or both [default %default]"),
    make_option("--comparison-label", type = "character", default = "comparison",
                help = "Comparison label for long-format output [default %default]"),
    make_option("--output-prefix", type = "character", default = "TF_scores",
                help = "Filename prefix for outputs [default %default]")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "score_tf_activity_dorothea_decouplr.R --deg-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (is.null(args$`deg-csv`) || !nzchar(args$`deg-csv`)) stop("--deg-csv is required")
  if (is.null(args$`output-dir`) || !nzchar(args$`output-dir`)) stop("--output-dir is required")
  if (!file.exists(args$`deg-csv`)) stop("DEG CSV not found: ", args$`deg-csv`)

  confidence <- trimws(strsplit(args$confidence, ",")[[1]])
  method <- tolower(args$method)
  if (!(method %in% c("wmean", "mlm"))) stop("--method must be 'wmean' or 'mlm'")
  output_format <- tolower(args$`output-format`)
  if (!(output_format %in% c("per-cluster", "long", "both"))) {
    stop("--output-format must be 'per-cluster', 'long', or 'both'")
  }
  excl <- if (is.null(args$`exclude-clusters`) || !nzchar(args$`exclude-clusters`)) character() else trimws(strsplit(args$`exclude-clusters`, ",")[[1]])

  deg <- readr::read_csv(args$`deg-csv`, show_col_types = FALSE)
  required_cols <- c("gene", "avg_log2FC", "cluster", "pct.1", "pct.2")
  if (!all(required_cols %in% names(deg))) {
    stop("DEG table missing required columns: ", paste(setdiff(required_cols, names(deg)), collapse = ", "))
  }

  regulon <- load_regulon(args$species, confidence)

  deg <- deg %>%
    filter(!(pct.1 < args$`pct-threshold` & pct.2 < args$`pct-threshold`)) %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(!(cluster %in% excl))

  clusters <- sort(unique(deg$cluster))
  if (!length(clusters)) stop("No clusters available after filtering")

  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  results_list <- list()
  for (cluster_id in clusters) {
    message("Scoring TF activity for cluster ", cluster_id, " with ", method)
    cluster_deg <- deg %>%
      filter(cluster == cluster_id) %>%
      select(gene, avg_log2FC) %>%
      distinct(gene, .keep_all = TRUE)

    matrix_input <- cluster_deg %>%
      column_to_rownames("gene") %>%
      as.matrix()

    res <- if (method == "wmean") {
      run_wmean(
        mat = matrix_input,
        network = regulon,
        .source = "tf",
        .target = "target",
        .mor = "mor",
        minsize = args$minsize,
        times = args$resamples,
        seed = args$seed
      )
    } else {
      run_mlm(
        mat = matrix_input,
        network = regulon,
        .source = "tf",
        .target = "target",
        .mor = "mor",
        minsize = args$minsize
      )
    }

    if ("p_value" %in% names(res)) {
      res <- arrange(res, p_value)
    } else if ("statistic" %in% names(res)) {
      res <- arrange(res, desc(abs(statistic)))
    } else if ("score" %in% names(res)) {
      res <- arrange(res, desc(abs(score)))
    }

    res <- res %>%
      mutate(cluster = cluster_id, .after = condition)
    if ("p_value" %in% names(res) && !("p.adjust" %in% names(res))) {
      res <- res %>% mutate(`p.adjust` = p.adjust(p_value, method = "BH"))
    }

    results_list[[cluster_id]] <- res

    if (output_format %in% c("per-cluster", "both")) {
      output_path <- file.path(args$`output-dir`, paste0(args$`output-prefix`, "_cluster_", cluster_id, "_", method, ".csv"))
      readr::write_csv(res, output_path)
    }
  }

  if (output_format %in% c("long", "both")) {
    long_res <- bind_rows(results_list) %>%
      mutate(cluster = as.integer(as.character(cluster)))
    long_path <- file.path(
      args$`output-dir`,
      paste0(args$`output-prefix`, "_long-format_", args$`comparison-label`, "_", method, ".csv")
    )
    readr::write_csv(long_res, long_path)
  }
}

if (sys.nframe() == 0L) {
  main()
}
