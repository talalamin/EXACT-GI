#!/usr/bin/env Rscript
#Differential Exon Usage (DEU) pipeline using limma-voom
# - Per-gene fits are run on the full exon set but using gene-specific sample subsets
# #TALAL AMIN: BIOBIX (UGENT)

# ---------------- LOAD PACKAGES ---------------------------------------

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(readxl); library(openxlsx)
  library(ggplot2); library(viridis); library(dplyr); library(tidyr)
  library(stringr); library(tools); library(matrixStats); library(scales)
  library(ggtext)
  library(BiocParallel)
})

# ---------------- PATHS & SETTINGS ---------------------------------------
# Replace the following placeholders with your actual file locations
ctrl_rds      <- "path/to/Control_counts_DEU.rds"   # Se/DEU RDS for control samples
tum_rds       <- "path/to/Tumor_counts_DEU.rds"     # Se/DEU RDS for tumor samples
design_dir    <- "path/to/design_dir/"               # per-gene CSV files (Gene.csv)
gene_map_file <- "path/to/gene_ID.xlsx"              # Excel mapping: Gene -> ID (ensembl or similar)

out_dir       <- "path/to/output/DEU"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plots_root    <- file.path(out_dir, "per_contrast_plots"); dir.create(plots_root, recursive = TRUE, showWarnings = FALSE)
qc_dir        <- file.path(out_dir, "QC"); dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
status_dir    <- file.path(out_dir, "status"); dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)

# Candidate genes (example) — edit as required
candidate_genes <- c("IGF2", "ZNF331", "HM13", "MEST") #Candidate genes (as examples givens)

# Contrast definitions: user-friendly name -> vector(level1, level2)
contrast_defs <- list(
  "Tumor vs Control"           = c("Tumor", "Control"),
  "Tumor_LOI vs Tumor_noLOI"   = c("Tumor_LOI", "Tumor_noLOI"),
  "Tumor_LOI vs Control_noLOI" = c("Tumor_LOI", "Control_noLOI")
)

# ---------------- QC / expensive-step FLAGS ----------------------------------
RUN_PVAL_HIST            <- TRUE #Optional
COMPUTE_UNFILTERED_PVALS <- TRUE #Optional
UNFILTERED_SAMPLE_N      <- 20000L #Optional
UNFILTERED_CHUNK_SIZE    <- 5000L #Optional
RUN_BCV_PLOT             <- TRUE #Optional
RUN_PVAL_VIOLIN          <- TRUE #Optional
DROP_WITHIN_PATIENT_REPLICATES <- TRUE #Optional
# ------------------------------------------------------------------------------

# logging
log_file <- file.path(out_dir, "pipeline_log.txt"); writeLines(character(), log_file)
log <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, paste(..., collapse = " ")), file = log_file, append = TRUE)
}
options(error = function() {
  tb <- capture.output(traceback(max.lines = 50))
  cat("UNCAUGHT ERROR\n", file = log_file, append = TRUE)
  cat(paste(tb, collapse = "\n"), file = log_file, append = TRUE)
  q(status = 1)
})

log("START: DEU pipeline (voom/limma only; per-gene fits on full exon set)")

# ---------------- helpers ----------------------------------------------------
safe_sheet_name <- function(s) { s2 <- gsub("[^-A-Za-z0-9_ ]", "_", s); s3 <- gsub("\\s+", " ", trimws(s2)); substr(s3, 1, 31) }
clean_label <- function(s) gsub("_", " ", s)

coerce_pergene <- function(rds_obj) {
  if (is.null(rds_obj)) return(NULL)
  if (is.data.frame(rds_obj)) {
    df <- as.data.frame(rds_obj)
    if (!"contrast" %in% colnames(df)) df$contrast <- "unknown"
    return(df)
  }
  if (is.list(rds_obj)) {
    df_list <- list(); nm <- names(rds_obj); if (is.null(nm)) nm <- paste0("c", seq_along(rds_obj))
    for (i in seq_along(rds_obj)) {
      el <- rds_obj[[i]]; df <- NULL
      if (is.data.frame(el)) df <- el
      else if (is.list(el) && !is.null(el$result) && is.data.frame(el$result)) df <- el$result
      else { cand <- Filter(is.data.frame, el); if (length(cand)) df <- as.data.frame(cand[[1]]) }
      if (!is.null(df)) {
        if (!"contrast" %in% colnames(df)) df$contrast <- if (!is.null(nm[i]) && nzchar(nm[i])) nm[i] else paste0("c", i)
        df_list[[length(df_list) + 1]] <- df
      }
    }
    if (length(df_list) == 0) return(NULL)
    return(bind_rows(df_list))
  }
  NULL
}

extract_exon_short <- function(exon_id) {
  s <- as.character(exon_id)
  s <- sub("^.*:", "", s)
  if (grepl("^E", s)) return(s)
  paste0("E", s)
}

compute_unfiltered_pvalues_subsample <- function(counts_all, sample_group_global, nrows = 20000L, cache = NULL, seed = 1) {
  if (!is.null(cache) && file.exists(cache)) {
    res <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(res)) {
      log("Loaded cached unfiltered p-values from:", cache)
      return(res)
    }
  }
  N_total <- nrow(counts_all)
  use_n <- min(as.integer(nrows), N_total)
  set.seed(seed)
  rows_samp <- if (use_n >= N_total) seq_len(N_total) else sample.int(N_total, use_n)

  grp_all <- sample_group_global$group[match(colnames(counts_all), sample_group_global$sample)]
  counts_s <- counts_all[rows_samp, , drop = FALSE]
  grp_for_cols <- grp_all[match(colnames(counts_s), colnames(counts_all))]
  grp_for_cols <- ifelse(grepl("Tumor", grp_for_cols, ignore.case = TRUE), "Tumor",
                         ifelse(grepl("Control", grp_for_cols, ignore.case = TRUE), "Control", NA_character_))
  uniqg <- unique(na.omit(grp_for_cols))
  if (length(uniqg) < 2) {
    log(" Subsampled columns do not contain both Control and Tumor. Skipping unfiltered p-value computation for this subsample.")
    return(list(rows_samp = rows_samp, p_unfiltered = NULL, method = "subsample_skip_no_both_groups"))
  }

  ytmp <- DGEList(counts = counts_s)
  ytmp <- calcNormFactors(ytmp)
  ytmp$samples$group <- factor(grp_for_cols, levels = c("Control", "Tumor"))

  ytmp <- tryCatch({
    estimateCommonDisp(ytmp)
  }, error = function(e) {
    warning("estimateCommonDisp on subset failed: ", conditionMessage(e)); ytmp
  })

  et <- tryCatch({
    exactTest(ytmp, pair = c("Control", "Tumor"))
  }, error = function(e) {
    stop("exactTest on sampled unfiltered rows failed: ", conditionMessage(e))
  })

  pvals <- as.numeric(topTags(et, n = Inf)$table$PValue)
  names(pvals) <- rownames(counts_s)
  res <- list(rows_samp = rows_samp, p_unfiltered = pvals, method = "subsample_common_disp_exactTest")
  if (!is.null(cache)) tryCatch(saveRDS(res, cache), error = function(e) warning("Failed to save cache: ", conditionMessage(e)))
  log("Computed unfiltered (subsampled) p-values; rows:", length(pvals))
  res
}

# ---------------- 1) Load inputs ------------------------------------------------
log("1) Loading inputs")
stopifnot(file.exists(ctrl_rds), file.exists(tum_rds), file.exists(gene_map_file))
se_c <- readRDS(ctrl_rds); colnames(se_c) <- file_path_sans_ext(basename(colnames(se_c)))
se_t <- readRDS(tum_rds);  colnames(se_t) <- file_path_sans_ext(basename(colnames(se_t)))
gene_map <- read_excel(gene_map_file)
log(sprintf("  Loaded: control=%d, tumor=%d", ncol(se_c), ncol(se_t)))

# transcript metadata
rr <- rowRanges(se_c); mcols_rr <- mcols(rr)
tx_field_candidates <- c("tx_name", "transcript_name", "transcript_id", "tx_id")
tx_field_matches <- intersect(tx_field_candidates, colnames(mcols_rr))
if (length(tx_field_matches) == 0) {
  cl_cols <- colnames(mcols_rr)[sapply(mcols_rr, function(x) inherits(x, "CharacterList"))]
  tx_field <- if (length(cl_cols) >= 1) cl_cols[1] else NA_character_
} else tx_field <- tx_field_matches[1]

if (is.na(tx_field) || !tx_field %in% colnames(mcols_rr)) {
  tx_vec <- rep(NA_character_, length(rr))
} else {
  tx_cl <- tryCatch(mcols_rr[[tx_field]], error = function(e) { log("Error extracting tx field:", conditionMessage(e)); NULL })
  if (is.null(tx_cl)) tx_vec <- rep(NA_character_, length(rr)) else {
    tx_vec <- sapply(tx_cl, function(z) { if (is.null(z) || length(z) == 0) return(NA_character_); paste(as.character(z), collapse = ";") }, USE.NAMES = FALSE)
  }
}
tx_meta_all <- data.frame(exon_id = rownames(se_c), tx_name = tx_vec, stringsAsFactors = FALSE); rownames(tx_meta_all) <- tx_meta_all$exon_id
saveRDS(tx_meta_all, file = file.path(out_dir, "tx_meta_all.rds"))
log("  Saved tx_meta_all")

# ---------------- 2) Merge counts & cache ------------------------------------
counts_file <- file.path(out_dir, "counts_all.rds")
if (file.exists(counts_file)) {
  counts_all <- readRDS(counts_file); log("  Loaded cached counts_all.rds")
} else {
  cnt_c <- assay(se_c, "counts"); cnt_t <- assay(se_t, "counts")
  stopifnot(identical(rownames(cnt_c), rownames(cnt_t)))
  counts_all <- cbind(cnt_c, cnt_t)
  saveRDS(counts_all, counts_file); log("  Built and saved counts_all.rds")
}
all_samples <- colnames(counts_all)

# ---------------- sample grouping: global Tumor/Control -----------------------
sample_group_global <- data.frame(sample = all_samples, group = NA_character_, stringsAsFactors = FALSE)
sample_group_global$group[sample_group_global$sample %in% colnames(se_t)] <- "Tumor"
sample_group_global$group[sample_group_global$sample %in% colnames(se_c)] <- "Control"
sample_group_global$group[is.na(sample_group_global$group)] <- ifelse(grepl("-01", sample_group_global$sample), "Tumor", "Control")
group_global <- factor(sample_group_global$group, levels = c("Control", "Tumor"))
design_global <- model.matrix(~0 + group_global); colnames(design_global) <- c("Control", "Tumor"); rownames(design_global) <- sample_group_global$sample
log(sprintf("  Global groups: Control=%d, Tumor=%d", sum(group_global == "Control"), sum(group_global == "Tumor")))

# ---------------- 3) GLOBAL FILTER & NORM -----------------------------------
filt_rds <- file.path(out_dir, "y_all_filtered.rds")
if (file.exists(filt_rds)) {
  y_all <- readRDS(filt_rds); log("  Loaded cached y_all_filtered.rds")
} else {
  log("  Running global CPM filtering (CPM >= 0.5 in >= 14 samples) + calcNormFactors")
  y0 <- DGEList(counts_all, genes = data.frame(geneid = sub(":.*", "", rownames(counts_all)), exonid = rownames(counts_all), stringsAsFactors = FALSE))
  min_samples_keep <- 14L
  cpm_mat <- cpm(y0)
  keep <- rowSums(cpm_mat >= 0.5) >= min_samples_keep
  y_all <- y0[keep, , keep.lib.sizes = FALSE]
  y_all <- calcNormFactors(y_all)
  saveRDS(y_all, filt_rds); log("  Saved y_all_filtered.rds (rows:", nrow(y_all), ")")
}

# ---------------- 4) GLOBAL VOOM (one-time) ---------------------------------
v_all_rds <- file.path(out_dir, "v_all_voom.rds")
design_for_voom <- design_global
design_for_voom <- design_for_voom[colnames(y_all), , drop = FALSE]
rownames(design_for_voom) <- colnames(y_all)

v_all <- NULL
if (file.exists(v_all_rds)) {
  v_all <- readRDS(v_all_rds)
  if (!is.null(v_all$E) && !is.null(dim(v_all$E)) && !is.null(v_all$weights) && !is.null(dim(v_all$weights))) {
    log("  Loaded cached v_all_voom.rds (aligned).")
  } else {
    v_all <- NULL
  }
}

if (is.null(v_all)) {
  log("  Computing global voom (one-time) on filtered features; samples:", ncol(y_all), "features:", nrow(y_all))
  v_all <- voom(y_all, design_for_voom, plot = FALSE)

  if (is.null(dim(v_all$E))) {
    nr <- nrow(y_all); nc <- ncol(y_all)
    v_all$E <- matrix(v_all$E, nrow = nr, ncol = nc, dimnames = list(rownames(y_all)[seq_len(nr)], colnames(y_all)[seq_len(nc)]))
  } else {
    if (is.null(rownames(v_all$E))) rownames(v_all$E) <- rownames(y_all)[seq_len(nrow(v_all$E))]
    if (is.null(colnames(v_all$E))) colnames(v_all$E) <- colnames(y_all)[seq_len(ncol(v_all$E))]
  }

  if (is.null(v_all$weights) || is.null(dim(v_all$weights))) {
    v_all$weights <- matrix(1, nrow = nrow(v_all$E), ncol = ncol(v_all$E), dimnames = list(rownames(v_all$E), colnames(v_all$E)))
  } else {
    if (is.null(rownames(v_all$weights))) rownames(v_all$weights) <- rownames(v_all$E)[seq_len(nrow(v_all$weights))]
    if (is.null(colnames(v_all$weights))) colnames(v_all$weights) <- colnames(v_all$E)[seq_len(ncol(v_all$weights))]
  }

  saveRDS(v_all, v_all_rds); log("  Saved v_all_voom.rds (rows/cols):", nrow(v_all$E), "/", ncol(v_all$E))
}

# ---------------- 5) P-VALUE HISTOGRAM & UNFILTERED (FULL, CHUNKED) -----------
pval_hist_rds <- file.path(qc_dir, "pvalue_hist_data.rds")
pval_hist_png <- file.path(qc_dir, "pvalue_hist_before_after.png")
pvals_unf_cache_full <- file.path(qc_dir, "pvals_unfiltered_full.rds")
pvals_filt_cache_full <- file.path(qc_dir, "pvals_filtered_full.rds")

p_unfiltered <- NULL
p_filtered <- NULL

if (RUN_PVAL_HIST) {
  log("P-value plots requested (RUN_PVAL_HIST = TRUE)")
  # Optionally compute unfiltered p-values in chunks/subsamples for QC if desired.
} else {
  log("RUN_PVAL_HIST = FALSE -> skipping p-value computations/plots")
}

# ---------------- 6) GLOBAL QC & SAMPLE SUMMARIES --------------------------
log("6) Global QC & sample summaries")
libsizes <- colSums(counts_all)
n_zero <- colSums(counts_all == 0)
pct_zero <- 100 * n_zero / nrow(counts_all)
sample_stats <- data.frame(sample = names(libsizes), libsize = libsizes, n_zero = n_zero, pct_zero = pct_zero, stringsAsFactors = FALSE)
patient_id <- substr(sample_stats$sample, 1, 12)
patient_counts <- table(patient_id)
duplicates_per_patient <- data.frame(patient_id = names(patient_counts), n_samples = as.integer(patient_counts), stringsAsFactors = FALSE)
dupe_summary <- as.data.frame(table(patient_counts)); colnames(dupe_summary) <- c("n_samples_per_patient", "freq_patients")
filter_summary <- data.frame(InitialBins = nrow(counts_all), RetainedBins = nrow(y_all), RemovedBins = nrow(counts_all) - nrow(y_all))
wb_qc <- createWorkbook()
addWorksheet(wb_qc, "sample_stats"); writeData(wb_qc, "sample_stats", sample_stats)
addWorksheet(wb_qc, "duplicates_per_patient"); writeData(wb_qc, "duplicates_per_patient", duplicates_per_patient)
addWorksheet(wb_qc, "dupe_summary"); writeData(wb_qc, "dupe_summary", dupe_summary)
addWorksheet(wb_qc, "filter_summary"); writeData(wb_qc, "filter_summary", filter_summary)
saveWorkbook(wb_qc, file.path(out_dir, "QC_summary.xlsx"), overwrite = TRUE)
log("  Saved QC_summary.xlsx")

# Optionally drop within-patient replicates (keeps first sample per patient/group)
if (DROP_WITHIN_PATIENT_REPLICATES) {
  log("Dropping within-patient replicates (flag TRUE). For each patient and group keep first sample only.")
  pid <- substr(colnames(y_all), 1, 12)
  group_for_cols <- sample_group_global$group[match(colnames(y_all), sample_group_global$sample)]
  tbl <- data.frame(sample = colnames(y_all), patient = pid, group = group_for_cols, stringsAsFactors = FALSE)
  idx_keep <- !duplicated(paste0(tbl$patient, "::", tbl$group))
  keep_samples <- tbl$sample[idx_keep]
  y_all <- y_all[, keep_samples, keep.lib.sizes = TRUE]
  v_all <- voom(y_all, design_global[colnames(y_all), , drop = FALSE], plot = FALSE)
  log("  After dropping within-patient replicates, samples retained:", ncol(y_all))
}

# ---------------- 7) PER-GENE DEU (voom on full exon set but per-gene sample subsets) -
log("7) Running per-gene DEU (voom on full exon set per-gene sample subsets)")

num_workers <- min(length(candidate_genes), 4L)
bp_param <- tryCatch({ MulticoreParam(workers = num_workers) }, error = function(e) { SerialParam() })

process_one_gene <- function(g) {
  log(" Processing gene:", g)
  writeLines(as.character(Sys.time()), file.path(status_dir, paste0(g, ".started")))
  ens <- gene_map$ID[gene_map$Gene == g][1]
  exs_gene <- if (!is.na(ens)) grep(sprintf("^%s:", ens), rownames(y_all), value = TRUE) else character(0)
  if (length(exs_gene) < 1) { log("  No exons retained for", g, "- skipping"); file.create(file.path(status_dir, paste0(g, ".done"))); return(NULL) }

  design_file <- file.path(design_dir, paste0(g, ".csv"))
  design_tbl <- if (file.exists(design_file)) { tmp <- read.csv(design_file, stringsAsFactors = FALSE); tmp$samples <- file_path_sans_ext(basename(tmp$samples)); tmp } else NULL

  results_per_contrast <- list()

  for (sh in names(contrast_defs)) {
    lvl <- contrast_defs[[sh]]
    log("  -- contrast:", sh)
    if (sh == "Tumor vs Control") {
      ds2 <- data.frame(samples = colnames(y_all), group = sample_group_global$group[match(colnames(y_all), sample_group_global$sample)], stringsAsFactors = FALSE)
    } else {
      if (is.null(design_tbl)) { log("   design missing for", g, "- skip contrast", sh); next }
      tmp <- subset(design_tbl, group %in% lvl); ds2 <- tmp
    }

    samp <- intersect(ds2$samples, colnames(y_all))
    if (length(samp) < 2 && !is.null(design_tbl)) {
      s_try <- toupper(trimws(ds2$samples)); samp <- colnames(y_all)[toupper(colnames(y_all)) %in% s_try]
    }
    log("   Samples after intersect/fallback:", length(samp))
    if (length(samp) < 2) { log("   Not enough samples, skipping"); next }

    ds2_rows <- if (sh == "Tumor vs Control") data.frame(samples = samp, group = ifelse(samp %in% colnames(se_t), "Tumor", "Control"), stringsAsFactors = FALSE) else {
      tmp_rows <- ds2[match(samp, ds2$samples), , drop = FALSE]
      if (nrow(tmp_rows) != length(samp)) tmp_rows <- data.frame(samples = samp, group = ifelse(samp %in% colnames(se_t), "Tumor", "Control"), stringsAsFactors = FALSE)
      tmp_rows
    }
    samp_counts_tbl <- as.data.frame(table(factor(ds2_rows$group, levels = lvl))); colnames(samp_counts_tbl) <- c("group","N")
    samples_per_group <- paste0(samp_counts_tbl$group, ":", samp_counts_tbl$N, collapse = ";")
    log("   Per-group counts:", samples_per_group)

    y_sub <- y_all[, samp, keep.lib.sizes = TRUE]

    design_g <- model.matrix(~0 + factor(ds2_rows$group, levels = lvl))
    colnames(design_g) <- paste0("cond", levels(factor(ds2_rows$group, levels = lvl)))
    rownames(design_g) <- ds2_rows$samples
    design_g <- design_g[match(colnames(y_sub), rownames(design_g)), , drop = FALSE]
    rownames(design_g) <- colnames(y_sub)

    reuse_v_all <- FALSE
    if (sh == "Tumor vs Control") {
      if (length(samp) == ncol(y_all) && all(sort(samp) == sort(colnames(y_all)))) reuse_v_all <- TRUE
    }

    if (reuse_v_all) {
      cidx <- match(colnames(y_sub), colnames(v_all$E))
      if (any(is.na(cidx))) {
        log("   Warning: some sample columns not found in v_all; falling back to per-gene voom")
        compute_voom_here <- TRUE
      } else {
        compute_voom_here <- FALSE
        v_sub <- v_all
        v_sub$E <- v_all$E[, cidx, drop = FALSE]
        v_sub$weights <- v_all$weights[, cidx, drop = FALSE]
        colnames(v_sub$E) <- colnames(y_sub)
        colnames(v_sub$weights) <- colnames(y_sub)
        rownames(v_sub$E) <- rownames(y_sub)
        rownames(v_sub$weights) <- rownames(y_sub)
        log("   Reused global v_all for contrast:", sh)
      }
    } else {
      compute_voom_here <- TRUE
    }

    if (compute_voom_here) {
      v_sub <- tryCatch({
        voom(y_sub, design_g, plot = FALSE)
      }, error = function(e) {
        log("   voom error for gene", g, "contrast", sh, ":", conditionMessage(e)); NULL
      })
      if (is.null(v_sub)) next

      if (is.null(dim(v_sub$E))) {
        nr <- nrow(y_sub); nc <- ncol(y_sub)
        v_sub$E <- matrix(v_sub$E, nrow = nr, ncol = nc, dimnames = list(rownames(y_sub)[seq_len(nr)], colnames(y_sub)[seq_len(nc)]))
      } else {
        if (is.null(rownames(v_sub$E))) rownames(v_sub$E) <- rownames(y_sub)[seq_len(nrow(v_sub$E))]
        if (is.null(colnames(v_sub$E))) colnames(v_sub$E) <- colnames(y_sub)[seq_len(ncol(v_sub$E))]
      }
      if (is.null(v_sub$weights) || is.null(dim(v_sub$weights))) {
        v_sub$weights <- matrix(1, nrow = nrow(v_sub$E), ncol = ncol(v_sub$E), dimnames = list(rownames(v_sub$E), colnames(v_sub$E)))
      } else {
        if (is.null(rownames(v_sub$weights))) rownames(v_sub$weights) <- rownames(v_sub$E)[seq_len(nrow(v_sub$weights))]
        if (is.null(colnames(v_sub$weights))) colnames(v_sub$weights) <- colnames(v_sub$E)[seq_len(ncol(v_sub$weights))]
      }
      log("   Computed voom for gene", g, "contrast", sh, "(samples:", ncol(v_sub$E), ")")
    }

    fit1 <- tryCatch({ lmFit(v_sub, design_g) }, error = function(e) { log("   lmFit error:", conditionMessage(e)); NULL })
    if (is.null(fit1)) next

    contrast_expr <- paste0("cond", lvl[1], " - cond", lvl[2])
    contrast_mat <- tryCatch({ makeContrasts(contrasts = contrast_expr, levels = design_g) }, error = function(e) { log("   makeContrasts error:", conditionMessage(e)); NULL })
    if (is.null(contrast_mat)) next

    fit2 <- tryCatch({ contrasts.fit(fit1, contrast_mat) }, error = function(e) { log("   contrasts.fit error:", conditionMessage(e)); NULL })
    if (is.null(fit2)) next

    fit2 <- tryCatch({ eBayes(fit2) }, error = function(e) { log("   eBayes warning:", conditionMessage(e)); fit2 })

    if (!is.null(y_sub$genes) && "geneid" %in% colnames(y_sub$genes)) {
      geneid_vec <- y_sub$genes$geneid
    } else {
      geneid_vec <- sub(":.*", "", rownames(v_sub$E))
    }
    if (length(geneid_vec) != nrow(v_sub$E)) geneid_vec <- sub(":.*", "", rownames(v_sub$E))
    fit2$genes <- data.frame(GeneID = geneid_vec, exonID = rownames(v_sub$E), stringsAsFactors = FALSE)

    ex_all <- tryCatch({ diffSplice(fit2, geneid = "GeneID") }, error = function(e) { log("   diffSplice error:", conditionMessage(e)); NULL })
    if (is.null(ex_all)) { log("   diffSplice returned NULL for", g, sh); next }

    tbl_all <- tryCatch({ topSplice(ex_all, coef = 1, test = "t", number = Inf) }, error = function(e) { log("   topSplice error:", conditionMessage(e)); NULL })
    if (is.null(tbl_all) || nrow(tbl_all) == 0) { log("   topSplice returned no rows for", g, sh); next }
    tbl_all <- as.data.frame(tbl_all)

    # Robust subsetting to gene
    gene_rows <- data.frame()
    gene_col_names <- intersect(c("GeneID", "geneid", "gene"), colnames(tbl_all))
    if (length(gene_col_names) > 0) {
      for (gn in gene_col_names) {
        if (!is.na(ens) && ens %in% tbl_all[[gn]]) {
          gene_rows <- tbl_all[tbl_all[[gn]] == ens, , drop = FALSE]
          break
        }
      }
    }
    if (nrow(gene_rows) == 0 && length(exs_gene) > 0) {
      common <- intersect(rownames(tbl_all), exs_gene)
      if (length(common) > 0) gene_rows <- tbl_all[rownames(tbl_all) %in% common, , drop = FALSE]
    }
    if (nrow(gene_rows) == 0 && !is.na(ens)) {
      pat <- paste0("^", ens, ":")
      matches <- grepl(pat, rownames(tbl_all))
      if (any(matches, na.rm = TRUE)) gene_rows <- tbl_all[matches, , drop = FALSE]
    }
    if (nrow(gene_rows) == 0) {
      exoncols <- intersect(c("exonID", "exonid", "exon", "exon_id"), colnames(tbl_all))
      if (length(exoncols) > 0) {
        for (ec in exoncols) {
          common2 <- intersect(tbl_all[[ec]], exs_gene)
          if (length(common2) > 0) {
            gene_rows <- tbl_all[tbl_all[[ec]] %in% common2, , drop = FALSE]
            break
          }
        }
      }
    }

    if (nrow(gene_rows) == 0) {
      log("   After subsetting, no rows for gene", g, "in contrast", sh)
      next
    }

    exon_col <- if ("exonID" %in% colnames(gene_rows)) "exonID" else if ("exonid" %in% colnames(gene_rows)) "exonid" else if ("exon" %in% colnames(gene_rows)) "exon" else NA_character_
    if (is.na(exon_col)) gene_rows$exon_id <- rownames(gene_rows) else gene_rows$exon_id <- as.character(gene_rows[[exon_col]])

    gene_rows$gene <- g; gene_rows$contrast <- sh
    gene_rows$tx_name <- tx_meta_all$tx_name[match(gene_rows$exon_id, rownames(tx_meta_all))]
    gene_rows$samples_per_group <- samples_per_group

    gene_rds <- file.path(out_dir, sprintf("%s_DEU_results.rds", g))
    if (!file.exists(gene_rds)) {
      saveRDS(list(gene_rows), gene_rds)
    } else {
      existing <- readRDS(gene_rds)
      if (!is.list(existing)) existing <- list(existing)
      existing[[sh]] <- gene_rows
      saveRDS(existing, gene_rds)
    }

    results_per_contrast[[sh]] <- gene_rows
    log(sprintf("   Collected %d rows for %s (%s)", nrow(gene_rows), g, sh))
  } # contrasts

  file.create(file.path(status_dir, paste0(g, ".done")))
  startedfile <- file.path(status_dir, paste0(g, ".started")); if (file.exists(startedfile)) file.remove(startedfile)
  log("  Finished gene:", g)
  results_per_contrast
}

res_list <- bplapply(candidate_genes, process_one_gene, BPPARAM = bp_param)
names(res_list) <- candidate_genes

per_contrast_results <- setNames(vector("list", length(contrast_defs)), names(contrast_defs))
for (i in seq_along(candidate_genes)) {
  g <- candidate_genes[i]; rpg <- res_list[[i]]
  if (is.null(rpg)) next
  for (sh in names(rpg)) {
    if (is.null(per_contrast_results[[sh]])) per_contrast_results[[sh]] <- list()
    per_contrast_results[[sh]][[length(per_contrast_results[[sh]]) + 1]] <- rpg[[sh]]
  }
}

# ---------------- 8) Save per-contrast workbook + dodge plots -----------------
log("8) Writing per-contrast results workbook and dodge plots")
deu_wb <- createWorkbook()
for (sh in names(per_contrast_results)) {
  sheet_name <- safe_sheet_name(sh)
  addWorksheet(deu_wb, sheet_name)
  lst <- per_contrast_results[[sh]]
  if (length(lst) == 0) {
    writeData(deu_wb, sheet_name, data.frame(note = sprintf("No results for contrast: %s", sh)))
    next
  }
  df_all <- bind_rows(lst)
  if (!"exon_id" %in% colnames(df_all)) {
    possible <- intersect(c("exonid", "exonID"), colnames(df_all))
    if (length(possible) > 0) df_all$exon_id <- df_all[[possible[1]]]
  }

  df_plot <- df_all %>% filter(gene %in% candidate_genes & exon_id %in% rownames(y_all))
  if (nrow(df_plot) > 0 && "logFC" %in% colnames(df_plot)) {
    df_plot <- df_plot %>% mutate(
      exon_short = vapply(exon_id, extract_exon_short, FUN.VALUE = character(1)),
      sig_flag = ifelse(!is.na(FDR) & FDR < 0.05 & !is.na(logFC) & abs(logFC) >= 0.5, TRUE, FALSE)
    ) %>% arrange(gene, exon_short)

    df_plot <- df_plot %>% group_by(gene) %>% mutate(
      exon_label_raw = exon_short,
      exon_label = ifelse(sig_flag, paste0(exon_short, "<span style='color:red'>*</span>"), exon_short)
    ) %>% ungroup()

    df_plot <- df_plot %>% mutate(xid = paste0(gene, "___", exon_short))
    df_plot$gene <- factor(df_plot$gene, levels = candidate_genes)
    df_plot <- df_plot %>% arrange(gene, exon_short)
    xid_levels <- unique(df_plot$xid)
    df_plot$xid <- factor(df_plot$xid, levels = xid_levels)

    xid_to_label <- setNames(df_plot$exon_label, df_plot$xid)

    total_exons <- length(xid_levels)
    per_exon_phys_in <- 0.25
    base_width <- 14
    plot_width <- max(base_width, per_exon_phys_in * total_exons)

    p <- ggplot(df_plot, aes(x = xid, y = logFC, fill = gene)) +
         geom_col(width = 0.9, color = NA) +
         labs(title = paste0("Differential exon usage: ", clean_label(sh)),
              x = "Exons", y = "logFC", fill = "Gene") +
         theme_minimal(base_size = 12) +
         theme(
           legend.position = "top",
           legend.direction = "horizontal",
           legend.box = "horizontal",
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           axis.text.x = element_markdown(angle = 45, hjust = 1, size = rel(0.8)),
           axis.title.x = element_text(margin = margin(t = 6)),
           plot.title = element_text(hjust = 0.5)
         ) +
         scale_fill_viridis_d(option = "D") +
         scale_x_discrete(labels = xid_to_label, expand = expansion(add = c(0, 0))) +
         scale_y_continuous(breaks = function(x) {
           rmin <- floor(min(x, na.rm = TRUE) - 1)
           rmax <- ceiling(max(x, na.rm = TRUE) + 1)
           seq(rmin, rmax, by = 0.5)
         })

    plot_file <- file.path(plots_root, paste0("dodge_allgenes_", gsub("\\s+", "_", sheet_name), ".png"))
    ggsave(plot_file, p, width = plot_width, height = 6, dpi = 300)
    log("  Saved dodge (single per-contrast plot, equal-width bars):", plot_file)
  }

  drop_cols <- c("geneid", "exonid", "exon_id", "samples_per_group")
  df_save <- df_all %>% select(-any_of(drop_cols))
  writeData(deu_wb, sheet_name, df_save)
}
saveWorkbook(deu_wb, file.path(out_dir, "DEU_limma_results_per_contrast.xlsx"), overwrite = TRUE)
log("  Wrote DEU_limma_results_per_contrast.xlsx")

# ---------------- 9) Per-gene violin plots (single workbook) -----------------
log("9) Making per-gene violin plots (single workbook output)")
violin_wb <- createWorkbook()
logcpm_all <- cpm(y_all, log = TRUE)

for (g in candidate_genes) {
  gene_rds <- file.path(out_dir, paste0(g, "_DEU_results.rds"))
  if (!file.exists(gene_rds)) { log("  no per-gene RDS for", g, "-> skip violin"); next }
  raw_obj <- readRDS(gene_rds); df <- coerce_pergene(raw_obj)
  if (is.null(df)) { log("  could not coerce", g); next }
  exon_col <- intersect(c("exonid", "exonID", "exon_id"), colnames(df))[1]
  if (!is.na(exon_col)) df$exon_id <- as.character(df[[exon_col]])
  exs <- unique(df$exon_id); exs <- exs[exs %in% rownames(y_all)]
  if (length(exs) == 0) { log("  no exons retained for", g); next }

  design_csv <- file.path(design_dir, paste0(g, ".csv"))
  if (file.exists(design_csv)) {
    dtmp <- read.csv(design_csv, stringsAsFactors = FALSE)
    if (!all(c("samples", "group") %in% colnames(dtmp))) {
      log("  design csv for", g, "missing required columns 'samples'/'group' -> fallback to global grouping")
      samp_use <- colnames(y_all)
      sample_map <- data.frame(sample = samp_use, parent = sample_group_global$group[match(samp_use, sample_group_global$sample)], subgroup = NA_character_, stringsAsFactors = FALSE)
    } else {
      dtmp$samples <- file_path_sans_ext(basename(dtmp$samples))
      samp_use <- intersect(dtmp$samples, colnames(y_all))
      if (length(samp_use) == 0) {
        log("  design csv samples for", g, "not found in filtered samples -> fallback to global grouping")
        samp_use <- colnames(y_all)
        sample_map <- data.frame(sample = samp_use, parent = sample_group_global$group[match(samp_use, sample_group_global$sample)], subgroup = NA_character_, stringsAsFactors = FALSE)
      } else {
        grp_vec <- dtmp$group[match(samp_use, dtmp$samples)]
        parent_vec <- ifelse(grepl("Tumor", grp_vec, ignore.case = TRUE), "Tumor",
                             ifelse(grepl("Control", grp_vec, ignore.case = TRUE), "Control", sample_group_global$group[match(samp_use, sample_group_global$sample)]))
        sample_map <- data.frame(sample = samp_use, parent = parent_vec, subgroup = grp_vec, stringsAsFactors = FALSE)
        log("  Using design.csv for", g, "-> samples:", length(samp_use))
      }
    }
  } else {
    samp_use <- colnames(y_all)
    sample_map <- data.frame(sample = samp_use, parent = sample_group_global$group[match(samp_use, sample_group_global$sample)], subgroup = NA_character_, stringsAsFactors = FALSE)
    log("  No design.csv for", g, "-> using global grouping for all filtered samples (n=", length(samp_use), ")")
  }

  samp_use <- unique(sample_map$sample)
  samp_use <- samp_use[samp_use %in% colnames(logcpm_all)]
  if (length(samp_use) == 0) { log("  after intersection no samples available for", g, "-> skipping"); next }

  logcpm_ex <- logcpm_all[exs, samp_use, drop = FALSE]

  parent_df <- as.data.frame(t(logcpm_ex)); parent_df$sample <- rownames(parent_df)
  parent_long <- pivot_longer(parent_df, -sample, names_to = "exon_id", values_to = "expr")
  parent_long <- left_join(parent_long, sample_map %>% select(sample, parent), by = "sample")
  parent_long$grp_label <- parent_long$parent

  sub_long <- data.frame()
  if (any(!is.na(sample_map$subgroup))) {
    subsamps <- sample_map$sample[!is.na(sample_map$subgroup)]
    subsamps <- subsamps[subsamps %in% colnames(logcpm_ex)]
    if (length(subsamps) > 0) {
      df_sub <- as.data.frame(t(logcpm_ex[, subsamps, drop = FALSE])); df_sub$sample <- rownames(df_sub)
      sub_long <- pivot_longer(df_sub, -sample, names_to = "exon_id", values_to = "expr")
      sub_long <- left_join(sub_long, sample_map %>% select(sample, subgroup), by = "sample")
      sub_long$grp_label <- sub_long$subgroup
    }
  }

  plot_df <- bind_rows(parent_long, sub_long)
  wanted <- c("Control", "Tumor", "Tumor LOI", "Tumor noLOI", "Control noLOI")
  plot_df$grp_label <- gsub("_", " ", plot_df$grp_label)
  plot_df$grp_label <- ifelse(is.na(plot_df$grp_label), NA, plot_df$grp_label)
  plot_df <- plot_df[!is.na(plot_df$grp_label) & plot_df$grp_label %in% wanted, , drop = FALSE]
  if (nrow(plot_df) == 0) {
    parent_df2 <- as.data.frame(t(logcpm_ex)); parent_df2$sample <- rownames(parent_df2)
    parent_long2 <- pivot_longer(parent_df2, -sample, names_to = "exon_id", values_to = "expr")
    parent_long2 <- left_join(parent_long2, sample_map %>% select(sample, parent), by = "sample")
    parent_long2$grp_label <- parent_long2$parent
    parent_long2$grp_label <- gsub("_", " ", parent_long2$grp_label)
    parent_long2 <- parent_long2[parent_long2$grp_label %in% c("Control", "Tumor"), ]
    plot_df <- parent_long2
    if (nrow(plot_df) == 0) { log("  after fallback still no samples for", g); next }
  }

  plot_df$exon_id <- factor(plot_df$exon_id, levels = unique(exs))
  title_main <- paste0(g, " expression among groups (filtered samples / design.csv)")
  p_v <- ggplot(plot_df, aes(x = grp_label, y = expr, fill = grp_label)) +
         geom_violin(trim = FALSE, draw_quantiles = 0.5) +
         geom_jitter(width = 0.12, size = 0.35, alpha = 0.6) +
         facet_wrap(~exon_id, scales = "free_y", ncol = 2) +
         labs(title = title_main, x = "Group", y = "log-CPM") +
         theme_minimal(base_size = 12) +
         theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
         scale_fill_viridis_d(option = "D")
  outdir_gene <- file.path(plots_root, g); dir.create(outdir_gene, recursive = TRUE, showWarnings = FALSE)
  vfile <- file.path(outdir_gene, paste0(g, "_violin.png"))
  ggsave(vfile, p_v, width = 12, height = max(4, 1.5 * ceiling(length(exs) / 2)))
  log("  saved violin image:", vfile)

  summary_df <- plot_df %>% group_by(exon_id, grp_label) %>% summarise(N = n(), Mean = mean(expr, na.rm = TRUE), Median = median(expr, na.rm = TRUE), .groups = "drop")
  addWorksheet(violin_wb, safe_sheet_name(g))
  writeData(violin_wb, sheet = safe_sheet_name(g), summary_df)
}
saveWorkbook(violin_wb, file.path(out_dir, "violin_summary_per_gene.xlsx"), overwrite = TRUE)
log("  Saved violin_summary_per_gene.xlsx")

# ---------------- 10) Final summaries & save -------------------------------
log("10) Saving final summaries and exiting")
cpm_filtered <- cpm(y_all)
gene_cpm <- rowsum(cpm_filtered, y_all$genes$geneid)
summary_rows <- list()
for (g in candidate_genes) {
  ens <- gene_map$ID[gene_map$Gene == g][1]
  gene_vals <- NULL
  if (!is.na(ens) && ens %in% rownames(gene_cpm)) gene_vals <- gene_cpm[ens, , drop = TRUE]
  else if (!is.na(ens)) {
    exs <- grep(sprintf("^%s:", ens), rownames(y_all), value = TRUE)
    if (length(exs) > 0) gene_vals <- colSums(cpm_filtered[exs, , drop = FALSE])
  }
  design_file <- file.path(design_dir, paste0(g, ".csv"))
  if (!file.exists(design_file)) {
    if (!is.null(gene_vals)) {
      for (grp in levels(group_global)) {
        vals <- gene_vals[group_global == grp]
        summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = grp, Mean = mean(vals, na.rm = TRUE), Median = median(vals, na.rm = TRUE), N = sum(!is.na(vals)), stringsAsFactors = FALSE)
      }
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = "Control", Mean = NA_real_, Median = NA_real_, N = 0, stringsAsFactors = FALSE)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = "Tumor", Mean = NA_real_, Median = NA_real_, N = 0, stringsAsFactors = FALSE)
    }
    next
  }
  ds <- read.csv(design_file, stringsAsFactors = FALSE); ds$samples <- file_path_sans_ext(basename(ds$samples))
  for (grp in unique(ds$group)) {
    samp_subset <- intersect(ds$samples[ds$group == grp], colnames(y_all))
    if (length(samp_subset) == 0) {
      up <- toupper(trimws(ds$samples[ds$group == grp])); samp_subset <- colnames(y_all)[toupper(colnames(y_all)) %in% up]
    }
    if (length(samp_subset) == 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = grp, Mean = NA_real_, Median = NA_real_, N = 0, stringsAsFactors = FALSE)
    } else {
      if (is.null(gene_vals)) {
        exs <- grep(sprintf("^%s:", ens), rownames(y_all), value = TRUE)
        if (length(exs) > 0) gene_vals_local <- colSums(cpm_filtered[exs, , drop = FALSE]) else gene_vals_local <- rep(NA_real_, ncol(y_all))
        vals <- gene_vals_local[match(samp_subset, colnames(y_all))]
      } else vals <- gene_vals[match(samp_subset, colnames(y_all))]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = grp, Mean = mean(vals, na.rm = TRUE), Median = median(vals, na.rm = TRUE), N = length(vals), stringsAsFactors = FALSE)
    }
  }
}
summary_df <- if (length(summary_rows) > 0) bind_rows(summary_rows) else data.frame()
wbStats <- createWorkbook()
addWorksheet(wbStats, "PerGene_Subgroup_Summary"); if (nrow(summary_df) > 0) writeData(wbStats, "PerGene_Subgroup_Summary", summary_df)
addWorksheet(wbStats, "FilterInfo"); writeData(wbStats, "FilterInfo", filter_summary)
addWorksheet(wbStats, "SampleStats"); writeData(wbStats, "SampleStats", sample_stats)
saveWorkbook(wbStats, file.path(out_dir, "summary_stats_and_sample_info.xlsx"), overwrite = TRUE)
log("  Saved summary_stats_and_sample_info.xlsx")

log("ALL DONE. Outputs written to:", out_dir)