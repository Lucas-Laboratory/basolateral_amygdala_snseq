#!/usr/bin/env Rscript

# Generate per-cluster pie charts and UMAP plots for MapMyCells annotations.
#
# Inputs
#   --umap-csv: CSV containing barcode, cluster, UMAP coordinates (required).
#   --annotations-csv: CSV mapping barcodes to MapMyCells annotations (required).
#   --name-tiers: comma-separated annotation column names to iterate (default: `class_name,subclass_name,supertype_name,cluster_name`).
#   --output-dir: directory for outputs (required).
#   --sample-name: label used in filenames (default derived from UMAP CSV).
#   --top-n: number of categories per cluster to keep before collapsing to "other" (default: 5).
#   --width/--height: PDF size for pie charts/UMAPs (default: 8 x 6).
#
# Dependencies: optparse, ggplot2, dplyr, readr, colorspace, ggpubr

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(colorspace)
})

make_colors <- function(categories) {
  n <- length(categories)
  if (n == 0L) return(setNames(character(), character()))
  combined <- c(
    colorspace::sequential_hcl(ceiling(n / 3), palette = "Hawaii"),
    colorspace::sequential_hcl(ceiling(n / 3), palette = "Batlow"),
    colorspace::sequential_hcl(ceiling(n / 3), palette = "Turku")
  )[seq_len(n)]
  setNames(combined, categories)
}

strip_bom <- function(df) {
  names(df) <- sub("^\ufeff", "", names(df))
  df
}

find_column <- function(df, candidates, label, required = TRUE) {
  exact_match <- candidates[candidates %in% names(df)]
  if (length(exact_match) > 0L) return(exact_match[1])

  case_match <- match(tolower(candidates), tolower(names(df)), nomatch = 0L)
  case_match <- case_match[case_match > 0L]
  if (length(case_match) > 0L) return(names(df)[case_match[1]])

  if (required) {
    stop(label, " column not found. Tried: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

clean_ids <- function(ids) {
  trimws(as.character(ids))
}

drop_replicate_from_prefix <- function(ids) {
  sub("-Rep[0-9]+_", "_", ids)
}

barcode_suffix <- function(ids) {
  sub("^.*_", "", ids)
}

pick_join_strategy <- function(cluster_ids, annotation_ids) {
  cluster_ids <- clean_ids(cluster_ids)
  annotation_ids <- clean_ids(annotation_ids)

  candidates <- list(
    exact = list(cluster = cluster_ids, annotation = annotation_ids),
    replicate_normalized = list(
      cluster = drop_replicate_from_prefix(cluster_ids),
      annotation = drop_replicate_from_prefix(annotation_ids)
    ),
    barcode_suffix = list(
      cluster = barcode_suffix(cluster_ids),
      annotation = barcode_suffix(annotation_ids)
    )
  )

  for (strategy_name in names(candidates)) {
    annotation_key <- candidates[[strategy_name]]$annotation
    valid_annotation_key <- annotation_key[!is.na(annotation_key) & nzchar(annotation_key)]
    if (anyDuplicated(valid_annotation_key) > 0L) next

    cluster_key <- candidates[[strategy_name]]$cluster
    n_matches <- sum(cluster_key %in% valid_annotation_key, na.rm = TRUE)
    if (n_matches > 0L) {
      return(list(
        name = strategy_name,
        cluster_key = cluster_key,
        annotation_key = annotation_key,
        n_matches = n_matches
      ))
    }
  }

  list(name = NA_character_, n_matches = 0L)
}

standardize_umap_data <- function(cluster_data) {
  barcode_col <- find_column(cluster_data, c("Barcode", "barcode", "cell_id", "cellid"), "UMAP barcode")
  cluster_col <- find_column(cluster_data, c("Cluster", "cluster", "seurat_cluster", "seurat_clusters", "cluster_id"), "UMAP cluster")
  x_col <- find_column(cluster_data, c("UMAP_1", "umap_1", "UMAP1", "umap1", "x"), "UMAP x")
  y_col <- find_column(cluster_data, c("UMAP_2", "umap_2", "UMAP2", "umap2", "y"), "UMAP y")

  cluster_data %>%
    transmute(
      Barcode = clean_ids(.data[[barcode_col]]),
      Cluster = .data[[cluster_col]],
      UMAP_1 = .data[[x_col]],
      UMAP_2 = .data[[y_col]]
    )
}

standardize_annotation_plot_data <- function(annotation_data, annotation_barcode_col) {
  annotation_cluster_col <- find_column(
    annotation_data,
    c("Cluster", "cluster", "seurat_cluster", "seurat_clusters", "cluster_id"),
    "Annotation cluster",
    required = FALSE
  )
  annotation_x_col <- find_column(
    annotation_data,
    c("UMAP_1", "umap_1", "UMAP1", "umap1", "x"),
    "Annotation UMAP x",
    required = FALSE
  )
  annotation_y_col <- find_column(
    annotation_data,
    c("UMAP_2", "umap_2", "UMAP2", "umap2", "y"),
    "Annotation UMAP y",
    required = FALSE
  )

  if (any(is.na(c(annotation_cluster_col, annotation_x_col, annotation_y_col)))) {
    stop(
      "No overlapping barcode/cell_id values found between UMAP and annotation CSVs, ",
      "and the annotation CSV does not contain fallback cluster/UMAP columns."
    )
  }

  annotation_data %>%
    mutate(
      Barcode = clean_ids(.data[[annotation_barcode_col]]),
      Cluster = .data[[annotation_cluster_col]],
      UMAP_1 = .data[[annotation_x_col]],
      UMAP_2 = .data[[annotation_y_col]]
    )
}

merge_umap_annotations <- function(cluster_data, annotation_data) {
  annotation_barcode_col <- find_column(
    annotation_data,
    c("Barcode", "barcode", "cell_id", "cellid"),
    "Annotation barcode"
  )

  cluster_data <- standardize_umap_data(cluster_data)
  annotation_ids <- annotation_data[[annotation_barcode_col]]
  join_strategy <- pick_join_strategy(cluster_data$Barcode, annotation_ids)

  if (!is.na(join_strategy$name)) {
    message(
      "Joining UMAP and annotations using ", join_strategy$name,
      " IDs; matched ", join_strategy$n_matches, " of ", nrow(cluster_data), " UMAP rows."
    )
    annotation_join_data <- annotation_data %>%
      mutate(
        annotation_cell_id = clean_ids(.data[[annotation_barcode_col]]),
        .join_key = join_strategy$annotation_key
      ) %>%
      select(-any_of(c("Barcode", "Cluster", "UMAP_1", "UMAP_2")))

    return(cluster_data %>%
      mutate(.join_key = join_strategy$cluster_key) %>%
      left_join(annotation_join_data, by = ".join_key") %>%
      select(-.join_key))
  }

  message(
    "No overlapping barcode/cell_id values found after exact, replicate-normalized, ",
    "or suffix matching. Using cluster and UMAP columns from the annotation CSV."
  )
  standardize_annotation_plot_data(annotation_data, annotation_barcode_col)
}

process_tier <- function(merged_data, tier, output_dir, sample_name, top_n, width, height) {
  categories <- sort(unique(na.omit(merged_data[[tier]])))
  if (length(categories) == 0L) {
    stop("No non-missing annotations found for tier: ", tier)
  }
  colors <- make_colors(categories)
  colors["other"] <- "gray"

  tier_data <- merged_data %>% mutate(color = colors[.data[[tier]]])
  write_csv(tier_data, file.path(output_dir, paste0("MapMyAlignments_", tier, "_", sample_name, ".csv")))

  # Pie charts per cluster
  plot_list <- tier_data %>%
    filter(!is.na(.data[[tier]])) %>%
    group_by(Cluster, .data[[tier]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(rank = rank(-n, ties.method = "first"), category_group = ifelse(rank <= top_n, .data[[tier]], "other")) %>%
    group_by(Cluster, category_group) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    split(.$Cluster) %>%
    lapply(function(df) {
      cluster_id <- unique(df$Cluster)
      ggplot(df, aes(x = "", y = n, fill = category_group)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        scale_fill_manual(values = colors) +
        theme_void() +
        ggtitle(paste("Cluster", cluster_id)) +
        geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3)
    })

  pdf(file.path(output_dir, paste0("MapMyPieCharts_", tier, "_", sample_name, ".pdf")), width = width, height = height)
  for (plot in plot_list) print(plot)
  dev.off()

  # Export top categories per cluster with colors
  top_n_df <- tier_data %>%
    filter(!is.na(.data[[tier]])) %>%
    group_by(Cluster, .data[[tier]]) %>% summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster) %>% mutate(rank = rank(-n, ties.method = "first")) %>%
    filter(rank <= top_n) %>% distinct(.data[[tier]]) %>% mutate(hex_color = colors[.data[[tier]]])
  write_csv(top_n_df, file.path(output_dir, paste0("top", top_n, "Color_", tier, "_", sample_name, ".csv")))

  # UMAP plot
  umap_plot <- ggplot(tier_data, aes(x = UMAP_1, y = UMAP_2, colour = .data[[tier]])) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = colors, na.value = "gray80") +
    theme_minimal() + labs(title = paste("UMAP coloured by", tier), colour = tier)
  ggsave(file.path(output_dir, paste0("UMAP_", tier, "_", sample_name, ".pdf")), umap_plot, width = width, height = height)
}

main <- function() {
  option_list <- list(
    make_option("--umap-csv", type = "character", help = "UMAP coordinate CSV"),
    make_option("--annotations-csv", type = "character", help = "MapMyCells annotation CSV"),
    make_option("--name-tiers", type = "character", default = "class_name,subclass_name,supertype_name,cluster_name"),
    make_option("--output-dir", type = "character", help = "Output directory"),
    make_option("--sample-name", type = "character", default = NULL),
    make_option("--top-n", type = "integer", default = 5),
    make_option("--width", type = "double", default = 8),
    make_option("--height", type = "double", default = 6)
  )

  parser <- OptionParser(option_list = option_list,
                         usage = "plot_mapmycells_pies_and_umap.R --umap-csv FILE --annotations-csv FILE --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`umap-csv`, args$`annotations-csv`, args$`output-dir`)))) stop("All inputs required")
  if (!file.exists(args$`umap-csv`)) stop("UMAP CSV not found")
  if (!file.exists(args$`annotations-csv`)) stop("Annotation CSV not found")
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  sample_name <- args$`sample-name`
  if (is.null(sample_name) || !nzchar(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(args$`umap-csv`))
  }

  cluster_data <- strip_bom(read_csv(args$`umap-csv`, show_col_types = FALSE))
  annotation_data <- strip_bom(read_csv(args$`annotations-csv`, show_col_types = FALSE))

  merged_data <- merge_umap_annotations(cluster_data, annotation_data)
  name_tiers <- trimws(strsplit(args$`name-tiers`, ",")[[1]])

  for (tier in name_tiers) {
    if (!tier %in% names(merged_data)) {
      warning("Skipping tier without column: ", tier)
      next
    }
    process_tier(merged_data, tier, args$`output-dir`, sample_name, args$`top-n`, args$width, args$height)
  }
}

if (sys.nframe() == 0L) {
  main()
}
