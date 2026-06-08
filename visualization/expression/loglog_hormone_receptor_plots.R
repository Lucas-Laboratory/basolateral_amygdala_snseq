#!/usr/bin/env Rscript

# Generate log–log CPM scatter plots and differential expression summaries for
# hormone receptor gated subsets using pseudobulk CPM and Wilcoxon tests.
#
# Inputs
#   --seurat-rds: Seurat object with an assay containing count data (required).
#   --metadata-csv: CSV/TSV with columns `Barcode` and `new_cluster` (required).
#   --output-dir: directory for plots, tables, and summaries (required).
#   --assay: Seurat assay supplying counts (default: `SCT`).
#   --gene-gates: comma-separated receptor genes (default: `Esr1,Esr2,Ar,Pgr`).
#   --gaba-clusters: numeric list/ranges defining GABA clusters (default: `8:19`).
#   --glut-clusters: numeric list/ranges defining glutamatergic clusters (default: `20:26,28,29`).
#   --conditions: comma-separated ordering of experimental groups (default: `Male-Naive,Female-Proestrus-Naive,Female-Diestrus-Naive`).
#   --thresholds: comma-separated expression thresholds aligned with `gene-gates` (default: all 1).
#   --pseudocount: value added before log10 CPM (default: 1e-6).
#   --padj-cutoff: adjusted p-value threshold for significance (default: 0.05).
#   --min-pct: Seurat `min.pct` for FindMarkers (default: 0).
#   --logfc-threshold: Seurat `logfc.threshold` (default: 0).
#
# Outputs
#   - Receptor count histograms PDF.
#   - DEG result CSVs per subset and comparison.
#   - DEG summary table, cell counts per gate, and threshold/sanity logs.
#   - Log–log scatter PDFs for GABA and Glut receptor-positive subsets.
#
# Dependencies: optparse, Seurat, Matrix, dplyr, ggplot2, readr, tidyr

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

parse_numeric_list <- function(spec) {
  if (!nzchar(spec)) return(integer())
  pieces <- trimws(strsplit(spec, ",")[[1]])
  vals <- integer(0)
  for (piece in pieces) {
    if (grepl(":", piece, fixed = TRUE)) {
      bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1]])
      if (length(bounds) != 2 || any(is.na(bounds))) {
        stop("Invalid numeric range: ", piece)
      }
      vals <- c(vals, seq(bounds[1], bounds[2]))
    } else {
      v <- suppressWarnings(as.integer(piece))
      if (is.na(v)) stop("Invalid numeric value: ", piece)
      vals <- c(vals, v)
    }
  }
  unique(vals)
}

parse_thresholds <- function(spec, genes) {
  parts <- as.numeric(trimws(strsplit(spec, ",")[[1]]))
  if (length(parts) == 1L) {
    parts <- rep(parts, length(genes))
  }
  if (length(parts) != length(genes)) {
    stop("--thresholds must supply one value per gene or a single value")
  }
  names(parts) <- genes
  parts
}

smart_read_table <- function(path) {
  header <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) "")
  n_tab <- length(gregexpr("\t", header, fixed = TRUE)[[1]])
  n_com <- length(gregexpr(",", header, fixed = TRUE)[[1]])
  use_tab <- is.finite(n_tab) && is.finite(n_com) && (n_tab > n_com)
  df <- if (use_tab) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  names(df) <- trimws(names(df))
  lower <- tolower(names(df))
  if ("barcode" %in% lower) names(df)[which(lower == "barcode")] <- "Barcode"
  if ("new_cluster" %in% lower) names(df)[which(lower == "new_cluster")] <- "new_cluster"
  df
}

pseudobulk_cpm <- function(cts, cell_mask, conditions, barcode_df) {
  if (!any(cell_mask)) {
    out <- Matrix(0, nrow = nrow(cts), ncol = length(conditions), sparse = TRUE)
    colnames(out) <- conditions
    rownames(out) <- rownames(cts)
    return(out)
  }
  cond <- factor(barcode_df$condition[cell_mask], levels = conditions)
  design <- Matrix::sparse.model.matrix(~ 0 + cond)
  colnames(design) <- levels(cond)
  pb <- cts[, which(cell_mask), drop = FALSE] %*% design
  missing_cols <- setdiff(conditions, colnames(pb))
  if (length(missing_cols)) {
    add <- Matrix(0, nrow = nrow(pb), ncol = length(missing_cols), sparse = TRUE)
    colnames(add) <- missing_cols
    rownames(add) <- rownames(pb)
    pb <- cbind(pb, add)[, conditions, drop = FALSE]
  } else {
    pb <- pb[, conditions, drop = FALSE]
  }
  library_sizes <- Matrix::colSums(pb)
  t(t(pb) / pmax(library_sizes, 1)) * 1e6
}

run_de_test <- function(seu, barcode_df, cells_mask, id1, id2, conditions, assay_name, min_pct, logfc_threshold) {
  if (!any(cells_mask)) return(NULL)
  sub <- subset(seu, cells = barcode_df$Barcode[cells_mask])
  sub$condition <- barcode_df$condition[match(colnames(sub), barcode_df$Barcode)]
  Idents(sub) <- factor(sub$condition, levels = conditions)
  if (sum(Idents(sub) == id1) < 3 || sum(Idents(sub) == id2) < 3) return(NULL)
  de <- FindMarkers(
    sub,
    ident.1 = id1,
    ident.2 = id2,
    test.use = "wilcox",
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    assay = assay_name,
    layer = "counts",
    features = rownames(sub@assays[[assay_name]])
  )
  de$gene <- rownames(de)
  de$test_used <- "wilcox"
  de
}

annotation_counts <- function(df, state_label) {
  df %>%
    filter(state == state_label) %>%
    count(comparison, receptor, name = "n")
}

plot_loglog <- function(df, title_txt, outfile) {
  if (!nrow(df)) return()
  counts_M <- annotation_counts(df, "up_in_M")
  counts_FR <- annotation_counts(df, "up_in_FR")
  counts_FNR <- annotation_counts(df, "up_in_FNR")

  p <- ggplot(df, aes(x = lx, y = ly, colour = state)) +
    geom_point(size = 0.3, alpha = 0.85) +
    scale_colour_manual(values = c(ns = "grey80", up_in_M = "black",
                                   up_in_FR = "maroon4", up_in_FNR = "turquoise4"),
                        guide = "none") +
    coord_fixed() +
    facet_grid(comparison ~ receptor) +
    labs(title = title_txt, x = "log10 CPM (Male-Naive)", y = "log10 CPM (comparison group)") +
    theme_classic(base_size = 10) +
    geom_text(data = counts_M, aes(x = -Inf, y = Inf, label = n),
              inherit.aes = FALSE, hjust = -0.2, vjust = 1.1, size = 3, colour = "black") +
    geom_text(data = counts_FR, aes(x = Inf, y = -Inf, label = n),
              inherit.aes = FALSE, hjust = 1.1, vjust = -0.2, size = 3, colour = "maroon4") +
    geom_text(data = counts_FNR, aes(x = Inf, y = -Inf, label = n),
              inherit.aes = FALSE, hjust = 1.1, vjust = -0.2, size = 3, colour = "turquoise4")

  ggsave(outfile, p, width = 9, height = 5.5)
  message("Wrote ", outfile)
}

main <- function() {
  option_list <- list(
    make_option("--seurat-rds", type = "character", help = "Seurat object (RDS)"),
    make_option("--metadata-csv", type = "character", help = "Barcode annotation CSV/TSV"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--assay", type = "character", default = "SCT"),
    make_option("--gene-gates", type = "character", default = "Esr1,Esr2,Ar,Pgr"),
    make_option("--gaba-clusters", type = "character", default = "8:19"),
    make_option("--glut-clusters", type = "character", default = "20:26,28,29"),
    make_option("--conditions", type = "character", default = "Male-Naive,Female-Proestrus-Naive,Female-Diestrus-Naive"),
    make_option("--thresholds", type = "character", default = "1"),
    make_option("--pseudocount", type = "double", default = 1e-6),
    make_option("--padj-cutoff", type = "double", default = 0.05),
    make_option("--min-pct", type = "double", default = 0),
    make_option("--logfc-threshold", type = "double", default = 0)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "loglog_hormone_receptor_plots.R --seurat-rds FILE --metadata-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  required <- c(args$`seurat-rds`, args$`metadata-csv`, args$`output-dir`)
  if (any(!nzchar(required))) stop("--seurat-rds, --metadata-csv, and --output-dir are required")
  if (!file.exists(args$`seurat-rds`)) stop("Seurat RDS not found")
  if (!file.exists(args$`metadata-csv`)) stop("Metadata CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  genes_gate <- trimws(strsplit(args$`gene-gates`, ",")[[1]])
  genes_gate <- genes_gate[nzchar(genes_gate)]
  if (!length(genes_gate)) stop("--gene-gates must list at least one gene")
  thr <- parse_thresholds(args$thresholds, genes_gate)

  gaba_clusters <- parse_numeric_list(args$`gaba-clusters`)
  glut_clusters <- parse_numeric_list(args$`glut-clusters`)
  if (!length(gaba_clusters) && !length(glut_clusters)) stop("At least one cluster must be provided for GABA or glut groups")

  conditions_order <- trimws(strsplit(args$conditions, ",")[[1]])
  conditions_order <- conditions_order[nzchar(conditions_order)]
  if (length(conditions_order) != 3) stop('--conditions must supply exactly three levels (Male, Female-Proestrus, Female-Diestrus by default)')

  message("Loading Seurat object")
  seu <- readRDS(args$`seurat-rds`)
  if (!inherits(seu, "Seurat")) stop("Input RDS must contain a Seurat object")
  assay_name <- args$assay
  if (!(assay_name %in% names(seu@assays))) stop("Assay ", assay_name, " not found in Seurat object")
  DefaultAssay(seu) <- assay_name
  cts <- GetAssayData(seu, layer = "counts")
  if (!inherits(cts, "dgCMatrix")) cts <- as(cts, "dgCMatrix")

  barcode_df <- smart_read_table(args$`metadata-csv`)
  if (!all(c("Barcode", "new_cluster") %in% colnames(barcode_df))) {
    stop("Metadata CSV must contain 'Barcode' and 'new_cluster' columns")
  }
  barcode_df$new_cluster <- suppressWarnings(as.integer(barcode_df$new_cluster))
  barcode_df <- barcode_df %>%
    mutate(condition = sub("^(.*)_.+$", "\\1", Barcode))

  common_cells <- intersect(colnames(seu), barcode_df$Barcode)
  if (!length(common_cells)) stop("No overlapping barcodes between Seurat object and metadata")
  seu <- seu[, common_cells, drop = FALSE]
  cts <- cts[, common_cells, drop = FALSE]
  barcode_df <- barcode_df[match(common_cells, barcode_df$Barcode), , drop = FALSE]

  hist_df <- lapply(genes_gate, function(g) {
    counts <- if (g %in% rownames(cts)) as.numeric(cts[g, ]) else {
      message("Gene not found, treating as zero: ", g)
      rep(0, ncol(cts))
    }
    data.frame(gene = g, count = counts, stringsAsFactors = FALSE)
  }) %>% bind_rows()
  hist_pdf <- file.path(args$`output-dir`, "receptor_count_histograms.pdf")
  ggsave(hist_pdf,
         ggplot(hist_df, aes(x = count)) +
           geom_histogram(bins = 50) +
           facet_wrap(~ gene, scales = "free_y") +
           labs(title = "Raw count distributions for receptor genes",
                x = "Counts (RNA assay, counts layer)", y = "Cells") +
           theme_classic(base_size = 12),
         width = 8, height = 6)
  message("Wrote ", hist_pdf)

  thr_path <- file.path(args$`output-dir`, "thresholds_used.txt")
  write_lines(
    c("Thresholds (counts >=) for receptor gating",
      sprintf("%s: %s", names(thr), format(thr, trim = TRUE))),
    thr_path
  )

  barcode_df$class <- ifelse(barcode_df$new_cluster %in% gaba_clusters, "GABA",
                             ifelse(barcode_df$new_cluster %in% glut_clusters, "Glut", NA_character_))
  keep <- !is.na(barcode_df$class)
  if (!any(keep)) stop("No cells fall into specified GABA or Glut cluster sets")
  barcode_df <- barcode_df[keep, , drop = FALSE]
  cts <- cts[, barcode_df$Barcode, drop = FALSE]

  for (g in genes_gate) {
    counts <- if (g %in% rownames(cts)) as.numeric(cts[g, ]) else rep(0, ncol(cts))
    barcode_df[[paste0(g, "_pos")]] <- counts >= thr[g]
  }
  barcode_df$None_pos <- !(Reduce(`|`, lapply(genes_gate, function(g) barcode_df[[paste0(g, "_pos")]])))
  genes_gate_ext <- c(genes_gate, "None")

  sanity_path <- file.path(args$`output-dir`, "sanity_checks.txt")
  sanity_lines <- c(
    paste0("assay_name: ", assay_name),
    paste0("assay layers: ", paste0(tryCatch(Layers(seu[[assay_name]]), error = function(e) character()), collapse = ", ")),
    paste0("n_common_cells: ", length(common_cells)),
    "",
    capture.output(print(table(barcode_df$class, useNA = "ifany"))),
    "",
    capture.output(print(table(barcode_df$condition, barcode_df$class)))
  )
  write_lines(sanity_lines, sanity_path)

  pseudocount <- args$pseudocount
  padj_cutoff <- args$`padj-cutoff`
  min_pct <- args$`min-pct`
  logfc_threshold <- args$`logfc-threshold`

  build_plot_df <- function(class_label, include_none = TRUE) {
    idx_class <- barcode_df$class == class_label
    receptors <- if (include_none) genes_gate_ext else genes_gate
    out <- vector("list", length(receptors))
    names(out) <- receptors

    for (g in receptors) {
      gate_mask <- if (g == "None") barcode_df$None_pos else barcode_df[[paste0(g, "_pos")]]
      idx <- idx_class & gate_mask
      if (!any(idx)) next
      cpm <- pseudobulk_cpm(cts, idx, conditions_order, barcode_df)
      male <- as.numeric(cpm[, conditions_order[1]])
      pro  <- if (length(conditions_order) >= 2) as.numeric(cpm[, conditions_order[2]]) else NA
      dies <- if (length(conditions_order) >= 3) as.numeric(cpm[, conditions_order[3]]) else NA

      base_df <- data.frame(gene = rownames(cts), receptor = g, stringsAsFactors = FALSE)

      comparisons <- list(
        list(name = paste0(conditions_order[1], ' vs ', conditions_order[2]),
             x = male, y = pro, id1 = conditions_order[1], id2 = conditions_order[2],
             state_labels = c('up_in_M', 'up_in_FR')),
        list(name = paste0(conditions_order[1], ' vs ', conditions_order[3]),
             x = male, y = dies, id1 = conditions_order[1], id2 = conditions_order[3],
             state_labels = c('up_in_M', 'up_in_FNR')),
        list(name = paste0(conditions_order[2], ' vs ', conditions_order[3]),
             x = pro, y = dies, id1 = conditions_order[2], id2 = conditions_order[3],
             state_labels = c('up_in_FR', 'up_in_FNR'))
      )

      comparison_frames <- list()
      for (cmp in comparisons) {
        if (any(is.na(c(cmp$x, cmp$y)))) next
        de <- run_de_test(seu, barcode_df, idx, cmp$id1, cmp$id2,
                          conditions_order, assay_name, min_pct, logfc_threshold)
        if (!is.null(de)) {
          out_csv <- file.path(args$`output-dir`, paste0(
            'deg_stats_', tolower(class_label), '_', g, '_',
            gsub(' ', '_', cmp$name, fixed = TRUE), '_wilcox.csv'))
          write_csv(de, out_csv)
        } else {
          de <- data.frame(gene = rownames(cts), avg_log2FC = NA_real_, p_val_adj = NA_real_)
        }
        df_cmp <- base_df
        df_cmp$comparison <- cmp$name
        df_cmp$lx <- log10(cmp$x + pseudocount)
        df_cmp$ly <- log10(cmp$y + pseudocount)
        df_cmp <- left_join(df_cmp, de[, c('gene', 'avg_log2FC', 'p_val_adj')], by = 'gene') %>%
          mutate(state = case_when(
            !is.na(p_val_adj) & p_val_adj < padj_cutoff & avg_log2FC > 0 ~ cmp$state_labels[1],
            !is.na(p_val_adj) & p_val_adj < padj_cutoff & avg_log2FC < 0 ~ cmp$state_labels[2],
            TRUE ~ 'ns'))
        comparison_frames[[cmp$name]] <- df_cmp
      }
      out[[g]] <- bind_rows(comparison_frames)
    }
    bind_rows(out)
  }

  plot_gaba <- build_plot_df("GABA", include_none = TRUE)
  plot_glut <- build_plot_df("Glut", include_none = TRUE)

  all_plot_df <- bind_rows(plot_gaba %>% mutate(class = "GABA"),
                           plot_glut %>% mutate(class = "Glut"))

  if (nrow(all_plot_df)) {
    deg_summary <- all_plot_df %>%
      group_by(class, receptor, comparison) %>%
      summarise(n_tested = n(),
                n_sig = sum(!is.na(p_val_adj) & p_val_adj < padj_cutoff),
                up_in_M = sum(state == 'up_in_M'),
                up_in_FR = sum(state == 'up_in_FR'),
                up_in_FNR = sum(state == 'up_in_FNR'),
                .groups = 'drop')

    cells_by_subset <- lapply(genes_gate_ext, function(g) {
      gate_mask <- if (g == 'None') barcode_df$None_pos else barcode_df[[paste0(g, '_pos')]]
      data.frame(class = barcode_df$class[gate_mask], receptor = g, stringsAsFactors = FALSE)
    }) %>% bind_rows() %>% count(class, receptor, name = 'n_cells_subset')

    deg_summary <- left_join(deg_summary, cells_by_subset, by = c('class', 'receptor')) %>%
      mutate(n_cells_subset = coalesce(n_cells_subset, 0L))

    write_csv(deg_summary, file.path(args$`output-dir`, 'deg_summary_by_subset_wilcox.csv'))
  }

  if (nrow(plot_gaba)) {
    plot_loglog(plot_gaba %>% mutate(class = 'GABA'),
                'Log–log CPM scatter (GABA; receptor-gated) — wilcox',
                file.path(args$`output-dir`, 'loglog_scatter_gaba_wilcox.pdf'))
  }
  if (nrow(plot_glut)) {
    plot_loglog(plot_glut %>% mutate(class = 'Glut'),
                'Log–log CPM scatter (Glut; receptor-gated) — wilcox',
                file.path(args$`output-dir`, 'loglog_scatter_glut_wilcox.pdf'))
  }

  cells_summary <- barcode_df %>%
    pivot_longer(cols = paste0(genes_gate, '_pos'), names_to = 'gate', values_to = 'pos') %>%
    mutate(gate = sub('_pos$', '', gate)) %>%
    group_by(class, gate, condition) %>%
    summarise(n_cells = sum(pos), .groups = 'drop') %>%
    arrange(class, gate, condition)
  write_csv(cells_summary, file.path(args$`output-dir`, 'cells_per_gate.csv'))

  message('Done')
}

if (sys.nframe() == 0L) {
  main()
}
