#!/usr/bin/env Rscript

# Generate RRHO2 heatmaps comparing cluster-wise DEG rankings across
# pairwise condition contrasts that share a common reference condition.
#
# Inputs
#   --input-dir: directory containing `merged_clusterwise_DEG_*.csv` files (required).
#   --output-dir: destination directory for RRHO2 outputs (required).
#   --plot-prefix: prefix prepended to output filenames (default: RRHO2).
#   --stepsize: RRHO2 grid step size controlling resolution (default: 500).
#   --color-palette: comma-separated list of hex colors for the nolabel heatmap
#                    (default: "#00007F,blue,#007FFF,cyan,#7FFF7F,yellow,#FF7F00,red,#7F0000").
#   --score-modes: comma-separated ranking metrics to compute; choose from
#                  `avg_log2FC`, `standard_score` (default: both).
#   --min-shared-genes: minimum shared genes required to plot a cluster (default: 10).
#
# Outputs
#   For each valid condition/cluster combination:
#     - RRHO2 heatmap PDF with labels.
#     - RRHO2 heatmap PNG with labels.
#     - RRHO2 heatmap PNG without labels, rendered via base image().
#   A debug log is written when RRHO2 encounters recoverable errors.
#
# Dependencies: optparse, RRHO2

suppressPackageStartupMessages({
  library(optparse)
  library(RRHO2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

parse_vector <- function(spec) {
  vals <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "", x)
}

prepare_score <- function(df, method = c("avg_log2FC", "standard_score")) {
  method <- match.arg(method)
  df$gene <- rownames(df)
  df$cluster_name <- sub("\\..*", "", df$gene)
  df$symbol <- sub(".*\\.", "", df$gene)
  split(df, df$cluster_name) |> lapply(function(sub_df) {
    genes <- sub_df$symbol
    if (method == "avg_log2FC") {
      score <- setNames(sub_df$avg_log2FC, genes)
    } else {
      pvals <- sub_df$p_val_adj
      pvals[pvals <= 0] <- 1e-300
      minus_log <- -log10(pvals)
      minus_log[minus_log > 300] <- 300
      score <- setNames(minus_log * sign(sub_df$avg_log2FC), genes)
    }
    score[!is.na(score)]
  })
}

make_pairs <- function(meta_df) {
  out <- lapply(unique(meta_df$condition_B), function(b) {
    subset_b <- meta_df[meta_df$condition_B == b, , drop = FALSE]
    if (nrow(subset_b) < 2) return(NULL)
    pairs <- combn(seq_len(nrow(subset_b)), 2, simplify = FALSE, FUN = function(idx) {
      data.frame(
        file1 = subset_b$file[idx[1]],
        condA1 = subset_b$condition_A[idx[1]],
        file2 = subset_b$file[idx[2]],
        condA2 = subset_b$condition_A[idx[2]],
        condB = b,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, pairs)
  })
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

plot_rrho2_set <- function(list_a, list_b, base_name, mode, outdir, label_x, label_y,
                           stepsize, palette, debug_store) {
  comp_x <- gsub(" ", "_", label_x)
  comp_y <- gsub(" ", "_", label_y)
  comp_type <- paste0(comp_x, "__vs__", comp_y)
  score_type <- mode

  dups_a <- names(list_a)[duplicated(names(list_a))]
  dups_b <- names(list_b)[duplicated(names(list_b))]
  if (length(dups_a)) message("Removed duplicated genes from list A: ", paste(unique(dups_a), collapse = ", "))
  if (length(dups_b)) message("Removed duplicated genes from list B: ", paste(unique(dups_b), collapse = ", "))
  list_a <- list_a[!duplicated(names(list_a))]
  list_b <- list_b[!duplicated(names(list_b))]

  df_a <- data.frame(gene = names(list_a), score = as.numeric(list_a), stringsAsFactors = FALSE)
  df_b <- data.frame(gene = names(list_b), score = as.numeric(list_b), stringsAsFactors = FALSE)

  rrho_result <- tryCatch(
    suppressWarnings(
      RRHO2_initialize(df_a, df_b,
                       labels = c(label_x, label_y),
                       log10.ind = TRUE,
                       stepsize = stepsize)
    ),
    error = function(e) {
      msg <- sprintf("INIT ERROR [%s_%s]: %s", base_name, mode, e$message)
      message(msg)
      debug_store$push(msg)
      return(NULL)
    }
  )
  if (is.null(rrho_result)) return(invisible(NULL))

  hypermat <- rrho_result$hypermat
  if (is.null(hypermat) || nrow(hypermat) < 2 || ncol(hypermat) < 2) {
    msg <- sprintf("Skipping %s_%s: hypermat too small (%s)", base_name, mode,
                   if (is.null(hypermat)) "NULL" else sprintf("%dx%d", nrow(hypermat), ncol(hypermat)))
    message(msg)
    debug_store$push(msg)
    return(invisible(NULL))
  }

  dir_pdf <- file.path(outdir, comp_type, score_type, "pdf")
  dir_png <- file.path(outdir, comp_type, score_type, "png")
  dir_png_nolabel <- file.path(outdir, comp_type, score_type, "png_nolabels")
  dir.create(dir_pdf, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_png, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_png_nolabel, recursive = TRUE, showWarnings = FALSE)

  pdf(file.path(dir_pdf, paste0(base_name, "_", mode, ".pdf")), width = 7, height = 6)
  tryCatch(RRHO2_heatmap(rrho_result, plot.legend = TRUE, plot.labels = TRUE),
           error = function(e) {
             msg <- sprintf("HEATMAP PDF ERROR [%s_%s]: %s", base_name, mode, e$message)
             message(msg)
             debug_store$push(msg)
           })
  dev.off()

  png(file.path(dir_png, paste0(base_name, "_", mode, ".png")), width = 800, height = 700)
  tryCatch(RRHO2_heatmap(rrho_result, plot.legend = TRUE, plot.labels = TRUE),
           error = function(e) {
             msg <- sprintf("HEATMAP PNG ERROR [%s_%s]: %s", base_name, mode, e$message)
             message(msg)
             debug_store$push(msg)
           })
  dev.off()

  png(file.path(dir_png_nolabel, paste0(base_name, "_", mode, "_nolabels.png")), width = 800, height = 700)
  minimum <- min(hypermat, na.rm = TRUE)
  maximum <- max(hypermat, na.rm = TRUE)
  colorGradient <- colorRampPalette(palette)(101)
  breaks <- seq(minimum, maximum, length.out = length(colorGradient) + 1)
  image(hypermat, col = colorGradient, breaks = breaks, axes = FALSE)
  dev.off()
}

main <- function() {
  option_list <- list(
    make_option("--input-dir", type = "character", help = "Directory with DEG CSVs"),
    make_option("--output-dir", type = "character", help = "Directory for outputs"),
    make_option("--plot-prefix", type = "character", default = "RRHO2"),
    make_option("--stepsize", type = "integer", default = 500),
    make_option("--color-palette", type = "character",
                default = "#00007F,blue,#007FFF,cyan,#7FFF7F,yellow,#FF7F00,red,#7F0000"),
    make_option("--score-modes", type = "character", default = "avg_log2FC,standard_score"),
    make_option("--min-shared-genes", type = "integer", default = 10),
    make_option("--pattern", type = "character", default = "^merged_clusterwise_DEG_.*\\.csv$",
                help = "Regex used to select DEG files")
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "rrho2_cluster_overlap_heatmaps.R --input-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`input-dir`, args$`output-dir`)))) {
    stop("--input-dir and --output-dir are required")
  }
  if (!dir.exists(args$`input-dir`)) stop("Input directory not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  palette <- parse_vector(args$`color-palette`)
  if (length(palette) < 3) stop("--color-palette must provide at least three colours")

  modes <- parse_vector(args$`score-modes`)
  allowed_modes <- c("avg_log2FC", "standard_score")
  if (!length(modes)) stop("--score-modes produced zero entries")
  invalid_modes <- setdiff(modes, allowed_modes)
  if (length(invalid_modes)) stop("Invalid score mode(s): ", paste(invalid_modes, collapse = ", "))

  files <- list.files(args$`input-dir`, pattern = args$pattern, full.names = TRUE)
  if (!length(files)) stop("No files matched pattern in input directory")

  file_meta <- data.frame(
    file = files,
    base = basename(files),
    stringsAsFactors = FALSE
  )
  file_meta$condition_pair <- gsub("^merged_clusterwise_DEG_|\\.csv$", "", file_meta$base)
  file_meta$condition_A <- sub("_vs_.*", "", file_meta$condition_pair)
  file_meta$condition_B <- sub(".*_vs_", "", file_meta$condition_pair)

  comparison_pairs <- make_pairs(file_meta)
  if (is.null(comparison_pairs) || !nrow(comparison_pairs)) {
    stop("Found no condition pairs sharing a common reference condition")
  }

  debug_messages <- new.env(parent = emptyenv())
  debug_messages$list <- character()
  debug_messages$push <- function(msg) {
    debug_messages$list <- c(debug_messages$list, msg)
  }

  for (i in seq_len(nrow(comparison_pairs))) {
    row <- comparison_pairs[i, ]
    message("Processing:", row$condA1, " vs ", row$condA2, " (common ", row$condB, ")")

    df1 <- read.csv(row$file1, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    df2 <- read.csv(row$file2, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

    score_lists1 <- setNames(vector("list", length(modes)), modes)
    score_lists2 <- setNames(vector("list", length(modes)), modes)
    for (mode in modes) {
      score_lists1[[mode]] <- prepare_score(df1, method = mode)
      score_lists2[[mode]] <- prepare_score(df2, method = mode)
    }

    common_clusters <- Reduce(intersect, list(names(score_lists1[[1]]), names(score_lists2[[1]])))
    for (clust in common_clusters) {
      vecs <- lapply(modes, function(m) list(score_lists1[[m]][[clust]], score_lists2[[m]][[clust]]))
      vecs <- unlist(vecs, recursive = FALSE)
      if (any(vapply(vecs, function(x) is.null(x) || is.null(names(x)) || !length(x), logical(1)))) {
        message("Skipping cluster ", clust, " due to missing score vectors")
        next
      }
      shared_genes <- Reduce(intersect, lapply(vecs, names))
      if (length(shared_genes) < args$`min-shared-genes`) {
        message("Skipping cluster ", clust, " (", length(shared_genes), " shared genes)")
        next
      }

      for (mode in modes) {
        score_lists1[[mode]][[clust]] <- score_lists1[[mode]][[clust]][shared_genes]
        score_lists2[[mode]][[clust]] <- score_lists2[[mode]][[clust]][shared_genes]
      }

      base_name <- paste(args$`plot-prefix`,
                         sanitize_name(row$condA1),
                         "vs",
                         sanitize_name(row$condA2),
                         "vs",
                         sanitize_name(row$condB),
                         "cluster",
                         clust,
                         sep = "_")

      label_x <- paste(row$condA1, "vs", row$condB)
      label_y <- paste(row$condA2, "vs", row$condB)

      for (mode in modes) {
        plot_rrho2_set(score_lists1[[mode]][[clust]],
                       score_lists2[[mode]][[clust]],
                       base_name,
                       mode,
                       args$`output-dir`,
                       label_x,
                       label_y,
                       args$stepsize,
                       palette,
                       debug_messages)
      }
    }
  }

  if (length(debug_messages$list)) {
    debug_file <- file.path(args$`output-dir`, "debug_log.txt")
    writeLines(debug_messages$list, con = debug_file)
    message("Wrote debug log to ", debug_file)
  }

  message("RRHO2 comparisons complete")
}

if (sys.nframe() == 0L) {
  main()
}
