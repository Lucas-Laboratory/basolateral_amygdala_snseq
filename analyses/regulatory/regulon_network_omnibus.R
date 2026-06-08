#!/usr/bin/env Rscript

# Build transcription factor → target networks using DoRothEA regulons, differential
# expression results, and optional STRING protein–protein interactions.
#
# Inputs
#   --tf-dir: directory containing decoupleR TF score CSVs (required).
#   --deg-dir: directory containing merged clusterwise DEG CSVs (required).
#   --output-dir: directory for CSV and PDF outputs (required).
#   --conditions: comma-separated TF score condition codes (default: PvD,PvM,DvM).
#   --condition-labels: comma-separated code=description pairs for DEG comparisons
#                       (default matches PvD/PvM/DvM).
#   --clusters: clusters to render (comma-separated integers and/or ranges; use
#               `all` to process every available cluster).
#   --species: organism code for DoRothEA/OrgDb lookups (`mm` or `hs`; default `mm`).
#   --dorothea-levels: DoRothEA confidence levels to retain (default: A,B).
#   --tf-p-cutoff / --deg-padj-cutoff: significance thresholds.
#   --direction-metrics: preference order for decoupleR statistics to determine TF
#                        directionality.
#   --tf-direction: restrict to `up`, `down`, or `any` TF direction (default: any).
#   --deg-gene-id-system: `auto`, `symbol`, or `ensembl` for DEG gene IDs (default: auto).
#   --string-layer: include STRING expansion (off by default).
#   --string-only-deg: limit STRING neighbours to DEG targets (otherwise allow
#                      high-confidence interactors).
#   --string-species / --string-score-min / --string-cache-dir: STRING options.
#   --max-a-per-tf / --max-b-per-tf / --max-string-per-target: per-layer caps.
#   --layout: ggraph layout algorithm (default: fr).
#   --pdf-size: square page size in inches (default: 6).
#   --draw-symbols: draw node points in addition to text labels.
#   --custom-prefix / --custom-suffix: appended to output filenames.
#   --log-file: optional path for run log (default: timestamped file in output dir).
#
# Dependencies: optparse, dplyr, tidyr, stringr, igraph, ggraph, ggplot2,
#               AnnotationDbi, dorothea, STRINGdb, grid, scales

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(AnnotationDbi)
})

suppressPackageStartupMessages({
  library(grid)
  library(scales)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a)) && !(is.character(a) && !nzchar(a[1]))) a else b

parse_vector <- function(spec) {
  vals <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

parse_clusters <- function(spec) {
  if (is.null(spec) || !nzchar(spec) || tolower(spec) == "all") return(NULL)
  parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  out <- integer(0)
  for (p in parts) {
    if (!nzchar(p)) next
    if (grepl(":", p, fixed = TRUE)) {
      bounds <- suppressWarnings(as.integer(strsplit(p, ":", fixed = TRUE)[[1]]))
      if (length(bounds) != 2 || any(is.na(bounds))) stop("Invalid cluster range: ", p)
      out <- c(out, seq(bounds[1], bounds[2]))
    } else {
      val <- suppressWarnings(as.integer(p))
      if (is.na(val)) stop("Invalid cluster value: ", p)
      out <- c(out, val)
    }
  }
  sort(unique(out))
}

parse_map <- function(spec, defaults) {
  mapping <- defaults
  if (is.null(spec) || !nzchar(spec)) return(mapping)
  parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  for (part in parts) {
    kv <- trimws(strsplit(part, "=", fixed = TRUE)[[1]])
    if (length(kv) != 2) stop("Invalid condition-label entry: ", part)
    mapping[kv[1]] <- kv[2]
  }
  mapping
}

option_list <- list(
  make_option("--tf-dir", type = "character", help = "Directory with TF score CSVs"),
  make_option("--deg-dir", type = "character", help = "Directory with DEG CSVs"),
  make_option("--output-dir", type = "character", help = "Directory for outputs"),
  make_option("--conditions", type = "character", default = "PvD,PvM,DvM"),
  make_option("--condition-labels", type = "character", default = NULL,
              help = "Comma-separated code=description pairs"),
  make_option("--clusters", type = "character", default = "8:29"),
  make_option("--species", type = "character", default = "mm"),
  make_option("--dorothea-levels", type = "character", default = "A,B"),
  make_option("--tf-p-cutoff", type = "double", default = 0.05),
  make_option("--deg-padj-cutoff", type = "double", default = 0.05),
  make_option("--direction-metrics", type = "character",
              default = "mlm,wmean_norm,norm_wmean,wmean,wmean_corr,corr_wmean"),
  make_option("--tf-direction", type = "character", default = "any"),
  make_option("--deg-gene-id-system", type = "character", default = "auto"),
  make_option("--include-arrowheads", action = "store_true", default = FALSE),
  make_option("--no-label-sign", action = "store_true", default = FALSE,
              help = "Disable colour-coding TF labels by direction"),
  make_option("--pos-label", type = "character", default = "firebrick4"),
  make_option("--neg-label", type = "character", default = "dodgerblue4"),
  make_option("--neutral-label", type = "character", default = "black"),
  make_option("--tf-fill", type = "character", default = "#FDD0A2"),
  make_option("--a-fill", type = "character", default = "#C6DBEF"),
  make_option("--b-fill", type = "character", default = "#C7E9C0"),
  make_option("--string-fill", type = "character", default = "#BDBDBD"),
  make_option("--tf-text-size", type = "double", default = 6.0),
  make_option("--a-text-size", type = "double", default = 5.0),
  make_option("--b-text-size", type = "double", default = 4.2),
  make_option("--string-text-size", type = "double", default = 3.6),
  make_option("--pdf-size", type = "double", default = 6.0),
  make_option("--layout", type = "character", default = "fr"),
  make_option("--custom-prefix", type = "character", default = ""),
  make_option("--custom-suffix", type = "character", default = "p05"),
  make_option("--max-a-per-tf", type = "integer", default = 25),
  make_option("--max-b-per-tf", type = "integer", default = 15),
  make_option("--max-string-per-target", type = "integer", default = 10),
  make_option("--string-layer", action = "store_true", default = FALSE,
              help = "Enable STRING third-layer expansion"),
  make_option("--string-species", type = "integer", default = 10090),
  make_option("--string-score-min", type = "integer", default = 900),
  make_option("--string-cache-dir", type = "character", default = NULL),
  make_option("--string-only-deg", action = "store_true", default = FALSE,
              help = "Restrict STRING interactors to DEG targets"),
  make_option("--disable-alias-mapping", action = "store_true", default = FALSE),
  make_option("--draw-symbols", action = "store_true", default = FALSE),
  make_option("--seed", type = "integer", default = 42),
  make_option("--log-file", type = "character", default = NULL)
)

# ---------- helper functions copied/adapted from legacy script ----------

pick_col <- function(nms, candidates) {
  ix <- which(tolower(nms) %in% tolower(candidates))
  if (length(ix)) nms[ix[1]] else NA_character_
}

count_fixed <- function(text, pat) {
  pos <- gregexpr(pat, text, fixed = TRUE)[[1]]
  if (length(pos) == 1L && pos[1] == -1L) 0L else length(pos)
}

smart_read_table <- function(path) {
  ln <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) "")
  n_tab <- count_fixed(ln, "\t"); n_com <- count_fixed(ln, ",")
  use_tab <- n_tab > n_com
  df <- if (use_tab) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  else              read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) >= 2) {
    nm1 <- names(df)[1]; v1 <- as.character(df[[1]])
    looks_rowname <- (is.null(nm1) || nm1 == "" || grepl("^X(\\n|$)?", nm1)) &&
      mean(grepl("^[0-9]+\\.[A-Za-z0-9_.:-]+$", v1)) > 0.5
    if (looks_rowname) df <- df[, -1, drop = FALSE]
  }
  df
}

read_tf_file <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)
  nms <- names(df)
  tf_col      <- pick_col(nms, c("source","tf","TF","Source"))
  stat_col    <- pick_col(nms, c("statistic","method","metric"))
  score_col   <- pick_col(nms, c("wmean_norm","norm_wmean","wmean","wmean_corr","corr_wmean","score","activity","estimate","value"))
  pval_col    <- pick_col(nms, c("p_value","pvalue","p.val","p_val","PValue"))
  cluster_col <- pick_col(nms, c("cluster","Cluster","cluster_id"))
  if (any(is.na(c(tf_col, stat_col, score_col, pval_col, cluster_col)))) return(NULL)
  out <- data.frame(
    TF        = as.character(df[[tf_col]]),
    statistic = as.character(df[[stat_col]]),
    score     = suppressWarnings(as.numeric(df[[score_col]])),
    p_value   = suppressWarnings(as.numeric(df[[pval_col]])),
    cluster   = suppressWarnings(as.integer(as.character(df[[cluster_col]]))),
    stringsAsFactors = FALSE
  )
  nm <- basename(path)
  m <- regmatches(nm, regexec("^TF_scores_long-format_([A-Za-z]+).*\\.csv$", nm))[[1]]
  out$cond_code <- if (length(m) >= 2) m[2] else sub("\\.csv$", "", nm)
  out
}

read_deg_file <- function(path) {
  df <- smart_read_table(path)
  if (!nrow(df)) return(NULL)
  chr_ix <- vapply(df, function(x) is.character(x) || is.factor(x), logical(1))
  if (any(chr_ix)) for (j in which(chr_ix)) df[[j]] <- trimws(as.character(df[[j]]))
  nms <- names(df)
  p_col     <- pick_col(nms, c("p_val","p.value","pval","p"))
  padj_col  <- pick_col(nms, c("p_val_adj","padj","adj.P.Val","qvalue","FDR","p.adjust","p_adj","padjust"))
  lfc_col   <- pick_col(nms, c("avg_log2FC","avg_log2fc","log2fc","logFC","avg_diff"))
  cl_col    <- pick_col(nms, c("cluster","Cluster","cluster_id"))
  comp_col  <- pick_col(nms, c("comparison","Comparison"))
  gene_col  <- pick_col(nms, c("gene_symbol","gene","symbol","GeneName","genes","target"))
  if (is.na(gene_col)) {
    char_cols <- nms[vapply(nms, function(x) is.character(df[[x]]) || is.factor(df[[x]]), logical(1))]
    excl <- tolower(c("comparison","cluster","p_val","p.value","pval","p_val_adj","padj","adj.P.Val","qvalue","fdr","p.adjust","p_adj","padjust","avg_log2FC","avg_log2fc","log2fc","logFC","avg_diff"))
    char_cols <- char_cols[!tolower(char_cols) %in% excl]
    if (length(char_cols)) {
      score <- vapply(char_cols, function(x) {
        v <- as.character(df[[x]]); v <- v[!is.na(v)]
        if (!length(v)) return(0)
        mean(grepl("^[A-Za-z0-9_.:-]+$", v))
      }, numeric(1))
      if (max(score) > 0.5) gene_col <- char_cols[which.max(score)]
    }
  }
  if (is.na(comp_col)) {
    bn <- basename(path)
    m <- regmatches(bn, regexec("^(mod_)?merged_clusterwise_DEG_(.*)\\.csv$", bn))[[1]]
    comp_inferred <- if (length(m) >= 3) m[3] else NA_character_
    if (!is.na(comp_inferred)) {
      df$comparison <- comp_inferred
      comp_col <- "comparison"
      nms <- names(df)
    }
  }
  to_num <- function(x) suppressWarnings(as.numeric(x))
  if (!is.na(p_col))    { tmp <- to_num(df[[p_col]]);    if (all(is.na(tmp)) && any(nchar(df[[p_col]])>0))    stop("Non-numeric column: ", p_col, " in ", basename(path));   df[[p_col]] <- tmp }
  if (!is.na(padj_col)) { tmp <- to_num(df[[padj_col]]); if (all(is.na(tmp)) && any(nchar(df[[padj_col]])>0)) stop("Non-numeric column: ", padj_col, " in ", basename(path)); df[[padj_col]] <- tmp }
  if (!is.na(lfc_col))  { tmp <- to_num(df[[lfc_col]]);  if (all(is.na(tmp)) && any(nchar(df[[lfc_col]])>0))  stop("Non-numeric column: ", lfc_col, " in ", basename(path));  df[[lfc_col]] <- tmp }
  req_missing <- c()
  if (is.na(cl_col))   req_missing <- c(req_missing, "cluster")
  if (is.na(comp_col)) req_missing <- c(req_missing, "comparison")
  if (is.na(gene_col)) req_missing <- c(req_missing, "gene")
  if (length(req_missing)) {
    stop("DEG file '", basename(path), "' missing: ", paste(req_missing, collapse = ", "),
         ". Available: ", paste(nms, collapse = ", "))
  }
  if (!is.na(comp_col)) {
    comp_vals <- as.character(df[[comp_col]])
    gene_vals <- as.character(df[[gene_col]])
    if (length(gene_vals) && length(comp_vals) && all(gene_vals %in% unique(comp_vals))) {
      stop("Gene column '", gene_col, "' likely mis-parsed as 'comparison'. Columns: ", paste(nms, collapse = ", "))
    }
  }
  data.frame(
    gene        = as.character(df[[gene_col]]),
    comparison  = as.character(df[[comp_col]]),
    cluster     = suppressWarnings(as.integer(as.character(df[[cl_col]]))),
    p_val       = if (!is.na(p_col))    df[[p_col]]    else NA_real_,
    p_val_adj   = if (!is.na(padj_col)) df[[padj_col]] else NA_real_,
    avg_log2FC  = if (!is.na(lfc_col))  df[[lfc_col]]  else NA_real_,
    stringsAsFactors = FALSE
  )
}

load_all_tf <- function(tf_dir, cond_codes) {
  pats <- paste0("^TF_scores_long-format_(", paste(cond_codes, collapse="|"), ").*\\.csv$")
  files <- list.files(tf_dir, pattern = pats, full.names = TRUE)
  if (!length(files)) stop("No TF files found in: ", tf_dir)
  out <- do.call(rbind, lapply(files, read_tf_file))
  if (is.null(out) || !nrow(out)) stop("No TF rows read.")
  out$statistic <- tolower(out$statistic)
  out
}

load_all_deg <- function(deg_dir) {
  files <- list.files(deg_dir, pattern = "^(mod_)?merged_clusterwise_DEG_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No DEG files matched in: ", deg_dir)
  out_list <- lapply(files, read_deg_file)
  out_list <- out_list[!vapply(out_list, is.null, logical(1))]
  out <- do.call(rbind, out_list)
  if (is.null(out) || !nrow(out)) stop("No DEG rows read from matched files in ", deg_dir)
  out
}

looks_like_ensembl <- function(x) {
  x <- toupper(sub("\\..*$", "", x))
  mean(grepl("^ENS", x, perl = TRUE), na.rm = TRUE) > 0.3
}

map_to_symbol <- function(ids, species_db) {
  ids <- sub("\\..*$", "", ids)
  keys <- unique(ids[!is.na(ids)])
  if (!length(keys)) return(data.frame(ENSEMBL = character(0), SYMBOL = character(0)))
  res <- suppressMessages(AnnotationDbi::select(species_db, keys = keys, columns = "SYMBOL", keytype = "ENSEMBL"))
  if (!nrow(res)) return(data.frame(ENSEMBL = character(0), SYMBOL = character(0)))
  if (!"ENSEMBL" %in% names(res)) {
    keycol <- setdiff(names(res), c("SYMBOL")); if (length(keycol) == 1L) names(res)[names(res) == keycol] <- "ENSEMBL"
  }
  res[!is.na(res$SYMBOL), c("ENSEMBL","SYMBOL"), drop = FALSE]
}

map_symbol_to_ens <- function(symbols, species_db) {
  syms <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(syms)) return(data.frame(SYMBOL = character(0), ENSEMBL = character(0)))
  ens <- suppressMessages(AnnotationDbi::mapIds(species_db, keys = syms, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first"))
  data.frame(SYMBOL = names(ens), ENSEMBL = unname(ens), stringsAsFactors = FALSE)
}

map_alias_to_ens <- function(aliases, species_db) {
  als <- unique(aliases[!is.na(aliases) & nzchar(aliases)])
  if (!length(als)) return(data.frame(INPUT = character(0), SYMBOL = character(0), ENSEMBL = character(0)))
  res <- tryCatch(suppressMessages(AnnotationDbi::select(species_db, keys = als, columns = c("SYMBOL","ENSEMBL"), keytype = "ALIAS")),
                  error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(data.frame(INPUT = character(0), SYMBOL = character(0), ENSEMBL = character(0)))
  names(res)[names(res) == "ALIAS"] <- "INPUT"
  if (any(is.na(res$ENSEMBL) & !is.na(res$SYMBOL))) {
    add <- map_symbol_to_ens(res$SYMBOL[is.na(res$ENSEMBL)], species_db)
    j <- match(res$SYMBOL, add$SYMBOL)
    res$ENSEMBL[is.na(res$ENSEMBL) & !is.na(j)] <- add$ENSEMBL[!is.na(j)]
  }
  res <- res[!is.na(res$SYMBOL) | !is.na(res$ENSEMBL), , drop = FALSE]
  res[!duplicated(res$INPUT), , drop = FALSE]
}

get_dorothea_regulon <- function(species = "mm", levels_keep = c("A","B")) {
  if (!requireNamespace("dorothea", quietly = TRUE)) stop("Package 'dorothea' not installed.")
  obj_name <- if (tolower(species) == "hs") "dorothea_hs" else "dorothea_mm"
  reg <- getExportedValue("dorothea", obj_name)
  if (is.function(reg)) reg <- reg()
  if (is.null(reg)) stop("Could not load DoRothEA object '", obj_name, "'.")
  keep_cols <- intersect(c("tf","target","mor","likelihood","confidence"), names(reg))
  reg <- reg[, keep_cols, drop = FALSE]
  if (!is.null(levels_keep)) reg <- reg[is.na(reg$confidence) | reg$confidence %in% levels_keep, , drop = FALSE]
  reg
}

string_map_symbols <- function(string_obj, symbols) {
  df <- data.frame(symbol = unique(symbols), stringsAsFactors = FALSE)
  mapped <- string_obj$map(df, "symbol", removeUnmappedRows = TRUE)
  data.frame(symbol = mapped$symbol, STRING_id = mapped$STRING_id, stringsAsFactors = FALSE)
}

node_color <- function(type) {
  if (type == "TF") return(col_TF)
  if (type == "STRING") return(col_STRING)
  if (type %in% names(conf_palette)) return(conf_palette[[type]])
  "#BBBBBB"
}

node_shape <- function(type) {
  if (type == "STRING") return(as.integer(shape_STRING))
  if (type %in% dorothea_levels) {
    idx <- match(type, dorothea_levels)
    conf_shapes[((idx - 1L) %% length(conf_shapes)) + 1L]
  } else 21L
}

node_size <- function(type) {
  if (type == "TF") return(tf_text_size)
  if (type == "STRING") return(STRING_text_size)
  if (type %in% names(conf_text_size)) return(conf_text_size[[type]])
  default_conf_text_size
}

run_pipeline <- function() {
  set.seed(run_seed)
  log_info("Start run. Output dir: ", output_dir)
  log_info("Log file: ", log_file)

  log_info("Loading TF files from: ", tf_dir)
  TF_all  <- load_all_tf(tf_dir, use_conditions)
  available_methods <- sort(unique(na.omit(tolower(TF_all$statistic))))
  log_info("TF rows loaded: ", nrow(TF_all), " | methods detected: ", paste(available_methods, collapse=", "))

  chosen_method <- intersect(tolower(direction_metric_order), available_methods)
  chosen_method <- if (length(chosen_method)) chosen_method[1] else if (length(available_methods)) available_methods[1] else NA_character_
  if (is.na(chosen_method)) stop("No usable TF scoring method found in TF files.")
  log_info("Using method for direction: ", chosen_method)
  direction_metrics <- chosen_method

  log_info("Loading DEG files from: ", deg_dir)
  DEG_all <- load_all_deg(deg_dir)
  log_info("DEG rows loaded: ", nrow(DEG_all))
  log_info("Comparisons detected: ", paste(unique(DEG_all$comparison), collapse = " | "))

  species_db <- if (tolower(species) == "mm") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Package 'org.Mm.eg.db' not installed.")
    getNamespace("org.Mm.eg.db")
  } else {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Package 'org.Hs.eg.db' not installed.")
    getNamespace("org.Hs.eg.db")
  }

  if (deg_gene_id_system == "ensembl" || (deg_gene_id_system == "auto" && looks_like_ensembl(DEG_all$gene))) {
    DEG_all$gene_ens <- sub("\\..*$", "", as.character(DEG_all$gene))
    m_sym <- map_to_symbol(DEG_all$gene, species_db)
    names(m_sym) <- c("ENSEMBL","SYMBOL")
    DEG_all <- DEG_all %>%
      left_join(m_sym, by = c("gene_ens" = "ENSEMBL")) %>%
      mutate(gene_symbol = ifelse(is.na(SYMBOL), as.character(gene), SYMBOL)) %>%
      select(-SYMBOL)
    log_info("DEG ID mode: ENSEMBL-like. Mapped to symbols where possible.")
  } else {
    DEG_all$gene_symbol <- as.character(DEG_all$gene)
    map_se <- map_symbol_to_ens(DEG_all$gene_symbol, species_db)
    DEG_all <- DEG_all %>% left_join(map_se, by = c("gene_symbol" = "SYMBOL"))
    DEG_all$gene_ens <- if ("ENSEMBL" %in% names(DEG_all)) DEG_all$ENSEMBL else NA_character_
    DEG_all$ENSEMBL <- NULL
    log_info("DEG ID mode: SYMBOL-like. Mapped to ENSEMBL for joins.")
  }
  if (!"gene_symbol" %in% names(DEG_all)) DEG_all$gene_symbol <- as.character(DEG_all$gene)
  if (!"gene_ens" %in% names(DEG_all))    DEG_all$gene_ens    <- sub("\\..*$", "", as.character(DEG_all$gene))

  reg_all <- get_dorothea_regulon(species = species, levels_keep = dorothea_levels)
  log_info("DoRothEA regulon rows: ", nrow(reg_all),
           " | confidence levels: ", paste(unique(na.omit(reg_all$confidence)), collapse = ", "))

  TF_dir_sig <- TF_all %>%
    group_by(TF, cond_code, cluster) %>%
    summarise(
      p_min     = min(p_value, na.rm = TRUE),
      dir_score = sum(ifelse(statistic %in% tolower(direction_metrics), score, 0), na.rm = TRUE),
      dir       = ifelse(dir_score > 0, "up", ifelse(dir_score < 0, "down", "neutral")),
      .groups = "drop"
    ) %>%
    mutate(sig_tf = is.finite(p_min) & p_min <= tf_p_cutoff)

  DEG_sig <- DEG_all %>% mutate(sig_deg = !is.na(p_val_adj) & p_val_adj <= deg_padj_cutoff)

  if (!requireNamespace("STRINGdb", quietly = TRUE) && include_string_layer) {
    stop("Package 'STRINGdb' not installed but --string-layer requested.")
  }
  if (!dir.exists(string_cache_dir)) dir.create(string_cache_dir, recursive = TRUE, showWarnings = FALSE)
  STRING <- if (include_string_layer) {
    tryCatch(
      STRINGdb::STRINGdb$new(version = "12", species = string_species,
                             score_threshold = string_score_min, input_directory = string_cache_dir),
      error = function(e) { log_info("STRINGdb init failed: ", conditionMessage(e), ". STRING expansion will be skipped."); NULL }
    )
  } else NULL
  log_info("STRING layer: ", ifelse(isTRUE(include_string_layer) && !is.null(STRING), "enabled", "disabled"))

  dir_tag <- tolower(if (tf_direction_keep %in% c("up","down","any")) tf_direction_keep else "any")
  layer_tag <- if (!isTRUE(include_string_layer)) "2layer" else if (isTRUE(string_include_all_sig)) "3layer" else "connectionslayeronly"

  for (cc in use_conditions) {
    comp_long <- unname(cond_map[[cc]]); if (is.na(comp_long)) next
    clust_vec <- if (is.null(clusters_to_plot)) sort(unique(TF_dir_sig$cluster[TF_dir_sig$cond_code == cc])) else clusters_to_plot
    log_info("Condition ", cc, " — clusters: ", paste(clust_vec, collapse = ", "))

    for (cl in clust_vec) {
      TF_union <- TF_dir_sig %>%
        filter(cond_code == cc, cluster == cl, sig_tf) %>%
        { if (identical(tf_direction_keep, "any")) . else filter(., dir == tf_direction_keep) } %>%
        distinct(TF) %>% pull(TF)
      if (!length(TF_union)) { log_info("[", cc, " cl", sprintf("%02d", cl), "] No TFs pass p ≤ ", tf_p_cutoff, " — skip."); next }

      DEG_sub <- DEG_sig %>%
        filter(comparison == comp_long, cluster == cl, sig_deg) %>%
        select(any_of(c("gene_symbol", "gene_ens", "avg_log2FC")))
      if (!nrow(DEG_sub)) { log_info("[", cc, " cl", sprintf("%02d", cl), "] No DEGs pass padj ≤ ", deg_padj_cutoff, " — skip."); next }
      if (!"gene_symbol" %in% names(DEG_sub)) DEG_sub$gene_symbol <- as.character(DEG_sub$gene %||% NA_character_)
      if (!"gene_ens" %in% names(DEG_sub))    DEG_sub$gene_ens    <- sub("\\..*$", "", as.character(DEG_sub$gene %||% NA_character_))
      if (!"avg_log2FC" %in% names(DEG_sub))  DEG_sub$avg_log2FC  <- NA_real_

      ens_to_symbol_lut <- DEG_sub %>% filter(!is.na(gene_ens)) %>% select(gene_ens, gene_symbol) %>% distinct()

      reg_sub <- reg_all %>%
        filter(tf %in% TF_union, is.na(confidence) | confidence %in% dorothea_levels) %>%
        select(any_of(c("tf","target","mor","likelihood","confidence")))

      ali_map <- if (isTRUE(use_alias_mapping)) map_alias_to_ens(reg_sub$target, species_db) else data.frame(INPUT=character(0),SYMBOL=character(0),ENSEMBL=character(0))
      sym_map <- map_symbol_to_ens(reg_sub$target, species_db)

      reg_sub$target_ens <- NA_character_
      if (nrow(ali_map)) {
        idx_a <- match(reg_sub$target, ali_map$INPUT)
        reg_sub$target_ens[!is.na(idx_a)] <- ali_map$ENSEMBL[!is.na(idx_a)]
      }
      if (any(is.na(reg_sub$target_ens)) && nrow(sym_map)) {
        idx_s <- match(reg_sub$target, sym_map$SYMBOL)
        reg_sub$target_ens[is.na(reg_sub$target_ens) & !is.na(idx_s)] <- sym_map$ENSEMBL[na.omit(idx_s[is.na(reg_sub$target_ens) & !is.na(idx_s)])]
      }

      reg_sub$target_symbol_mapped <- NA_character_
      if (nrow(ens_to_symbol_lut)) {
        idx_l <- match(reg_sub$target_ens, ens_to_symbol_lut$gene_ens)
        reg_sub$target_symbol_mapped[!is.na(idx_l)] <- ens_to_symbol_lut$gene_symbol[na.omit(idx_l)]
      }
      if (any(is.na(reg_sub$target_symbol_mapped)) && nrow(ali_map)) {
        idx_as <- match(reg_sub$target, ali_map$INPUT)
        reg_sub$target_symbol_mapped[is.na(reg_sub$target_symbol_mapped) & !is.na(idx_as)] <- ali_map$SYMBOL[na.omit(idx_as[is.na(reg_sub$target_symbol_mapped) & !is.na(idx_as)])]
      }
      if (any(is.na(reg_sub$target_symbol_mapped)) && nrow(sym_map)) {
        idx_ss <- match(reg_sub$target, sym_map$SYMBOL)
        reg_sub$target_symbol_mapped[is.na(reg_sub$target_symbol_mapped) & !is.na(idx_ss)] <- sym_map$SYMBOL[na.omit(idx_ss[is.na(reg_sub$target_symbol_mapped) & !is.na(idx_ss)])]
      }
      reg_sub$target_symbol_mapped[is.na(reg_sub$target_symbol_mapped)] <- reg_sub$target[is.na(reg_sub$target_symbol_mapped)]

      reg_deg <- reg_sub %>%
        inner_join(DEG_sub %>% mutate(target_symbol = gene_symbol), by = c("target_symbol_mapped" = "target_symbol")) %>%
        mutate(target_symbol = target_symbol_mapped)

      if (!nrow(reg_deg)) { log_info("[", cc, " cl", sprintf("%02d", cl), "] No overlaps between TF regulon and DEGs — skip."); next }

      edges_tf <- reg_deg %>%
        mutate(interaction = "TF_target",
               confidence = ifelse(is.na(confidence), "N/A", confidence))

      edges_ppi <- data.frame()
      if (isTRUE(include_string_layer) && !is.null(STRING)) {
        tf_symbols <- unique(edges_tf$tf)
        target_symbols <- unique(edges_tf$target_symbol)
        str_mapped <- string_map_symbols(STRING, c(tf_symbols, target_symbols))
        if (nrow(str_mapped)) {
          tf_map <- str_mapped[str_mapped$symbol %in% tf_symbols, ]
          target_map <- str_mapped[str_mapped$symbol %in% target_symbols, ]
          if (nrow(tf_map) && nrow(target_map)) {
            hits <- STRING$get_interactions(tf_map$STRING_id)
            hits <- hits[hits$to %in% target_map$STRING_id, , drop = FALSE]
            if (!string_include_all_sig) {
              hits <- hits[hits$to %in% target_map$STRING_id[target_map$symbol %in% target_symbols], , drop = FALSE]
            }
            if (nrow(hits)) {
              tf_lookup <- setNames(tf_map$symbol, tf_map$STRING_id)
              target_lookup <- setNames(target_map$symbol, target_map$STRING_id)
              hits$from_symbol <- tf_lookup[hits$from]
              hits$to_symbol <- target_lookup[hits$to]
              hits <- hits[!is.na(hits$from_symbol) & !is.na(hits$to_symbol), , drop = FALSE]
              hits <- hits %>% group_by(from_symbol) %>% arrange(desc(score), .by_group = TRUE) %>%
                mutate(rank = row_number()) %>% filter(rank <= max_STRING_per_target) %>% ungroup()
              edges_ppi <- hits %>%
                transmute(tf = from_symbol,
                          target_symbol = to_symbol,
                          interaction = "STRING_PPI",
                          confidence = as.character(score))
            }
          }
        }
      }

      edges_all <- bind_rows(edges_tf %>% select(tf, target_symbol, interaction, confidence, avg_log2FC),
                             edges_ppi %>% mutate(avg_log2FC = NA_real_))
      edges_all <- edges_all %>%
        mutate(source_symbol = tf,
               edge_color = ifelse(interaction == "TF_target", edge_col_TF_to_deg, edge_col_deg_to_str))

      tf_nodes <- data.frame(id = unique(edges_all$source_symbol[edges_all$interaction == "TF_target" &
                                                                   edges_all$source_symbol %in% TF_union]),
                             type = "TF", stringsAsFactors = FALSE)
      conf_nodes_list <- lapply(dorothea_levels, function(lv) {
        ids <- unique(edges_all$target_symbol[edges_all$interaction == "TF_target" & edges_all$confidence == lv])
        if (!length(ids)) return(data.frame())
        data.frame(id = ids, type = lv, stringsAsFactors = FALSE)
      })
      ids_string <- setdiff(
        unique(edges_all$target_symbol[edges_all$interaction == "STRING_PPI"]),
        unique(edges_all$target_symbol[edges_all$interaction == "TF_target"])
      )
      string_nodes <- if (length(ids_string)) data.frame(id = ids_string, type = "STRING", stringsAsFactors = FALSE)
                      else data.frame(id = character(0), type = character(0), stringsAsFactors = FALSE)

      nodes <- do.call(bind_rows, c(list(tf_nodes), conf_nodes_list, list(string_nodes))) %>%
        distinct(id, .keep_all = TRUE)

      nodes <- nodes %>%
        mutate(
          label       = id,
          fill_col    = vapply(type, node_color, character(1)),
          text_size   = vapply(type, node_size,  numeric(1)),
          pt_shape    = vapply(type, node_shape, integer(1)),
          pt_size     = ifelse(type == "TF", 2.0, 2.5),
          label_color = neutral_label_col
        )

      if (isTRUE(label_sign_when_all) && identical(tf_direction_keep, "any")) {
        tf_dir_map <- TF_dir_sig %>% filter(cond_code == cc, cluster == cl, sig_tf) %>% select(TF, dir)
        if (nrow(tf_dir_map)) {
          tf_mask <- nodes$type == "TF"
          match_idx <- match(nodes$id[tf_mask], tf_dir_map$TF)
          d_vec <- tf_dir_map$dir[match_idx]
          nodes$label_color[tf_mask & !is.na(d_vec) & d_vec == "up"]   <- pos_label_col
          nodes$label_color[tf_mask & !is.na(d_vec) & d_vec == "down"] <- neg_label_col
        }
        lfc_by_symbol <- reg_deg %>% group_by(target_symbol) %>% summarise(lfc_mean = mean(avg_log2FC, na.rm = TRUE), .groups = "drop")
        if (nrow(lfc_by_symbol)) {
          j <- match(nodes$id, lfc_by_symbol$target_symbol)
          s <- sign(lfc_by_symbol$lfc_mean[j])
          nodes$label_color[!is.na(s) & s > 0] <- pos_label_col
          nodes$label_color[!is.na(s) & s < 0] <- neg_label_col
        }
      }

      file_base <- paste0(custom_prefix, cc, "_", dir_tag, "_", layer_tag, "_cl", sprintf("%02d", cl),
                          if (nzchar(custom_suffix)) paste0("_", custom_suffix) else "")
      write.csv(nodes %>% mutate(symbol = id), file.path(output_dir, paste0("network_nodes_", file_base, ".csv")), row.names = FALSE)
      write.csv(edges_all,                file.path(output_dir, paste0("network_edges_", file_base, ".csv")), row.names = FALSE)

      g   <- graph_from_data_frame(edges_all, directed = TRUE, vertices = nodes)
      lay <- create_layout(g, layout = layout_algo)

      pdf_path <- file.path(output_dir, paste0("network_", file_base, ".pdf"))
      pdf(pdf_path, width = pdf_size_in, height = pdf_size_in)
      cap_size_mm <- 2

      p <- ggraph(lay) +
        geom_edge_link(aes(color = edge_color),
                       linewidth = edge_width_fixed,
                       arrow = if (isTRUE(include_arrowheads)) arrow(length = unit(arrow_head_mm, 'mm'), type = 'closed') else NULL,
                       start_cap = circle(cap_size_mm, 'mm'),
                       end_cap   = circle(cap_size_mm, 'mm'),
                       alpha = edge_alpha_fixed) +
        scale_edge_colour_identity()
      if (isTRUE(draw_symbols)) {
        p <- p + geom_node_point(aes(shape = pt_shape, fill = fill_col, size = pt_size),
                                 stroke = 0.3, color = "grey20") +
          scale_shape_identity() + scale_fill_identity()
      }
      p <- p +
        scale_size_identity() +
        geom_node_text(aes(label = label, size = text_size, color = label_color),
                       hjust = 0.5, vjust = 0.5) +
        scale_colour_identity() +
        theme_void() +
        ggtitle(paste0("TF network — ", cc, " (", cond_map[[cc]], "), cluster ", cl))
      print(p)
      dev.off()

      log_info("[", cc, " cl", sprintf("%02d", cl), "] TFs kept: ", length(TF_union),
               " | DEGs kept: ", nrow(DEG_sub),
               " | reg_deg: ", nrow(reg_deg),
               " | edges TF: ", nrow(edges_tf),
               " | edges STRING: ", nrow(edges_ppi),
               " | total edges: ", nrow(edges_all),
               " | wrote CSVs and PDF.")
    }
  }

  log_info("Done.")
}

main <- function() {
  parser <- OptionParser(option_list = option_list,
                         usage = "regulon_network_omnibus.R --tf-dir DIR --deg-dir DIR --output-dir DIR [options]")
  args <- parse_args(parser)

  if (!all(nzchar(c(args$`tf-dir`, args$`deg-dir`, args$`output-dir`)))) {
    stop("--tf-dir, --deg-dir, and --output-dir are required")
  }
  if (!dir.exists(args$`tf-dir`)) stop("TF directory not found: ", args$`tf-dir`)
  if (!dir.exists(args$`deg-dir`)) stop("DEG directory not found: ", args$`deg-dir`)
  dir.create(args$`output-dir`, recursive = TRUE, showWarnings = FALSE)

  use_conditions <<- parse_vector(args$conditions)
  if (!length(use_conditions)) stop("--conditions produced zero entries")
  default_map <- c(
    PvD = "Female-Proestrus-Naive vs Female-Diestrus-Naive",
    PvM = "Female-Proestrus-Naive vs Male-Naive",
    DvM = "Female-Diestrus-Naive vs Male-Naive"
  )
  cond_map_full <- parse_map(args$`condition-labels`, default_map)
  missing_labels <- setdiff(use_conditions, names(cond_map_full))
  if (length(missing_labels)) cond_map_full[missing_labels] <- missing_labels

  tf_dir               <<- args$`tf-dir`
  deg_dir              <<- args$`deg-dir`
  output_dir           <<- args$`output-dir`
  cond_map             <<- cond_map_full
  clusters_to_plot     <<- parse_clusters(args$clusters)
  species              <<- tolower(args$species)
  dorothea_levels      <<- parse_vector(args$`dorothea-levels`)
  tf_p_cutoff          <<- args$`tf-p-cutoff`
  deg_padj_cutoff      <<- args$`deg-padj-cutoff`
  direction_metric_order <<- parse_vector(args$`direction-metrics`)
  tf_direction_keep    <<- match.arg(tolower(args$`tf-direction`), c("any","up","down"))
  deg_gene_id_system   <<- match.arg(tolower(args$`deg-gene-id-system`), c("auto","symbol","ensembl"))
  include_arrowheads   <<- isTRUE(args$`include-arrowheads`)
  label_sign_when_all  <<- !isTRUE(args$`no-label-sign`)
  pos_label_col        <<- args$`pos-label`
  neg_label_col        <<- args$`neg-label`
  neutral_label_col    <<- args$`neutral-label`
  col_TF               <<- args$`tf-fill`
  col_A                <<- args$`a-fill`
  col_B                <<- args$`b-fill`
  col_STRING           <<- args$`string-fill`
  tf_text_size         <<- args$`tf-text-size`
  A_text_size          <<- args$`a-text-size`
  B_text_size          <<- args$`b-text-size`
  STRING_text_size     <<- args$`string-text-size`
  pdf_size_in          <<- args$`pdf-size`
  layout_algo          <<- args$layout
  custom_prefix        <<- args$`custom-prefix`
  custom_suffix        <<- args$`custom-suffix`
  max_A_per_TF         <<- args$`max-a-per-tf`
  max_B_per_TF         <<- args$`max-b-per-tf`
  max_STRING_per_target <<- args$`max-string-per-target`
  include_string_layer <<- isTRUE(args$`string-layer`)
  string_species       <<- args$`string-species`
  string_score_min     <<- args$`string-score-min`
  string_cache_dir     <<- args$`string-cache-dir` %||% file.path(output_dir, "stringdb_cache")
  string_include_all_sig <<- !isTRUE(args$`string-only-deg`)
  use_alias_mapping    <<- !isTRUE(args$`disable-alias-mapping`)
  draw_symbols         <<- isTRUE(args$`draw-symbols`)
  run_seed             <<- args$seed

  max_per_TF_by_level  <<- c(A = max_A_per_TF, B = max_B_per_TF, C = 10L, D = 8L, E = 6L)
  conf_text_size       <<- c(A = A_text_size, B = B_text_size)
  conf_palette         <<- c(A = col_A, B = col_B, C = col_TF, D = "#E5D3F2", E = "#FDE0DD")
  conf_shapes          <<- c(21L, 22L, 24L, 23L, 25L)
  shape_STRING         <<- 24L
  default_conf_text_size <<- 4.0
  arrow_head_mm        <<- 1.25
  edge_width_fixed     <<- 0.3
  edge_alpha_fixed     <<- 1
  edge_col_TF_to_deg   <<- "#444444"
  edge_col_deg_to_str  <<- "#888888"

  log_file_path <- args$`log-file` %||% file.path(output_dir, paste0("regulon_network_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  log_file <<- log_file_path
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)
  log_info <<- function(...) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- paste0("[", ts, "] ", paste(..., collapse = ""))
    cat(line, "\n", file = log_con)
    cat(line, "\n")
  }

  run_pipeline()
}

if (sys.nframe() == 0L) {
  main()
}
