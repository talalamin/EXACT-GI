#!/usr/bin/env Rscript
# DEU_pipeline_voom_only_globalExonFits_fixed_with_subtypes_final_updated.R
# #TALAL AMIN: BIOBIX (UGENT)
# - Optional subtype-aware contrasts if per-gene design + subtype table provided
# - Minimal changes from the original script; the subtype-vs-control neg-label bug is fixed
#
# USAGE:
#  - Edit the "PATHS & SETTINGS" block below to point to your data files and output directory.
#  - Per-gene design CSVs expected as DESIGN_DIR/<GENE>.csv with columns: samples, group
#  - Optional subtype spreadsheet must contain a sample ID column and a subtype column.


suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tools)
  library(matrixStats)
  library(scales)
  library(ggtext)
  # intentionally avoiding BiocParallel: this script runs genes serially
})

# ---------------- PATHS & SETTINGS ---------------------------------------
# === Replace the placeholders below with your actual file locations ===
ctrl_rds      <- "path/to/Control_counts_DEU.rds"         # RangedSummarizedExperiment or equivalent RDS for controls
tum_rds       <- "path/to/Tumor_counts_DEU.rds"           # RangedSummarizedExperiment or equivalent RDS for tumors
design_dir    <- "path/to/design_dir"                     # directory with per-gene design CSVs (GENE.csv)
gene_map_file <- "path/to/gene_ID.xlsx"                   # Excel mapping: sheet with columns Gene, ID
subtype_xlsx  <- "path/to/TCGA_with_subtypes.xlsx"        # Optional: sample -> subtype mapping (Excel)

out_dir       <- "path/to/output/DEU_subtypes_final_fixed_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plots_root    <- file.path(out_dir, "per_contrast_plots"); dir.create(plots_root, recursive = TRUE, showWarnings = FALSE)
qc_dir        <- file.path(out_dir, "QC"); dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
status_dir    <- file.path(out_dir, "status"); dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)

design_with_subtypes_dir <- file.path(design_dir, "design_with_subtypes"); dir.create(design_with_subtypes_dir, recursive = TRUE, showWarnings = FALSE)

# Candidate genes (examples; edit to your list)
candidate_genes <- c("ZNF331", "HM13")

# Main contrast definitions
contrast_defs <- list(
  "Tumor vs Control"           = c("Tumor", "Control"),
  "Tumor_LOI vs Tumor_noLOI"   = c("Tumor_LOI", "Tumor_noLOI"),
  "Tumor_LOI vs Control_noLOI" = c("Tumor_LOI", "Control_noLOI")
)

# Subtypes of interest (edit to your taxonomy)
SUBTYPES_OF_INTEREST <- c("Luminal A", "Luminal B", "HER2+", "Basal-like")

# If you only want specific exons for some genes, list their exon numbers (strings)
# e.g. KEEP_EXONS <- list(HM13 = c("67"), ZNF331 = c("4","5","10"))
KEEP_EXONS <- list(
  HM13   = c("67"),
  ZNF331 = c("4","5","10")
)

# ---------------- FLAGS & LOGGING ----------------------------------------
DROP_WITHIN_PATIENT_REPLICATES <- TRUE

# protections for subtype covariate usage:
MIN_SUBTYPE_SAMPLES <- 5      # minimum samples per subtype to be considered
MIN_SUBTYPE_LEVELS  <- 2      # require at least this many subtypes with >= MIN_SUBTYPE_SAMPLES

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

log("START: DEU pipeline v2 - subtype-covariate protected usage")

# ---------------- HELPERS -------------------------------------------------
safe_sheet_name <- function(s) {
  s2 <- gsub("[^-A-Za-z0-9_ ]", "_", s)
  s3 <- gsub("\\s+", " ", trimws(s2))
  substr(s3, 1, 31)
}
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

get_exon_number <- function(exon_id) {
  s <- as.character(exon_id)
  s2 <- sub("^.*:", "", s)
  s2 <- sub("^E", "", s2, ignore.case = TRUE)
  s2 <- gsub("^0+", "", s2)
  if (s2 == "") s2 <- "0"
  s2
}
extract_exon_short <- function(exon_id) {
  num <- get_exon_number(exon_id)
  sprintf("E%03s", num)
}

# ---------------- PREPROCESS: augment per-gene design CSVs with subtypes ----
log("Preprocessing per-gene design CSVs to include standardized tumor subtype information...")

if (!file.exists(subtype_xlsx)) {
  log("No subtype spreadsheet found at:", subtype_xlsx, "- script will continue but subtype-aware contrasts will be limited.")
  subtbl <- data.frame()
} else {
  subtbl <- read_excel(subtype_xlsx)
  # detect sample id and subtype columns flexibly
  sname_candidates <- intersect(c("Sample ID", "SampleID", "Sample_ID", "Sample.Id", "sample_id"), colnames(subtbl))
  if (length(sname_candidates) == 0) stop("Cannot find a Sample ID column in subtype xlsx: columns found: ", paste(colnames(subtbl), collapse = ","))
  samp_col <- sname_candidates[1]
  stype_candidates <- intersect(c("Subtype", "subtype", "Tumor Subtype", "Tumor_Type", "Tumor Type"), colnames(subtbl))
  if (length(stype_candidates) == 0) stop("Cannot find a Subtype column in subtype xlsx: columns found: ", paste(colnames(subtbl), collapse = ","))
  stype_col <- stype_candidates[1]

  standardize_subtype <- function(x) {
    if (is.na(x) || !nzchar(as.character(x))) return(NA_character_)
    s <- as.character(x); s <- gsub("^BRCA_", "", s, ignore.case = TRUE); s <- trimws(s); s_low <- tolower(s)
    if (grepl("lum", s_low) && grepl("a", s_low)) return("Luminal A")
    if (grepl("lum", s_low) && grepl("b", s_low)) return("Luminal B")
    if (grepl("her2", s_low) || grepl("her-2", s_low)) return("HER2+")
    if (grepl("basal", s_low)) return("Basal-like")
    if (grepl("normal", s_low)) return("Normal")
    return(s)
  }

  subtbl$Sample_ID <- as.character(subtbl[[samp_col]])
  # canonicalize sample IDs to TCGA-like prefix if present
  subtbl$Sample_ID <- sub("^([^\\-]+\\-[A-Z0-9]{2}\\-[A-Z0-9]{4}\\-[0-9]{2}).*", "\\1", subtbl$Sample_ID)
  subtbl$Subtype_raw <- as.character(subtbl[[stype_col]])
  subtbl$Subtype_std <- vapply(subtbl$Subtype_raw, standardize_subtype, FUN.VALUE = NA_character_, USE.NAMES = FALSE)
  subtbl_small <- subtbl %>% dplyr::select(Sample_ID, Subtype_std) %>% distinct()
}

# Write augmented per-gene designs (adds Subtype and group_subtyped)
for (g in candidate_genes) {
  design_file <- file.path(design_dir, paste0(g, ".csv"))
  if (!file.exists(design_file)) { log("Design missing for gene", g, "- skipping preprocessing"); next }
  d <- read.csv(design_file, stringsAsFactors = FALSE)
  if (!all(c("samples", "group") %in% colnames(d))) { log("Design for", g, "missing 'samples'/'group' -> skipping"); next }
  d$samples <- as.character(d$samples)
  # derive a Sample_ID to match subtype table; keep multiple heuristics for robustness
  d$Sample_ID <- vapply(d$samples, function(x) {
    s <- gsub("\\s+", "", as.character(x))
    if (grepl("-[0-9]{2}[A-Z]$", s)) return(sub("([0-9]{2})[A-Z]$", "\\1", s))
    if (grepl("-[0-9]{2}$", s)) return(s)
    return(s)
  }, FUN.VALUE = "")
  if (exists("subtbl_small") && nrow(subtbl_small) > 0) {
    dm <- left_join(d, subtbl_small, by = c("Sample_ID" = "Sample_ID"))
    dm$Subtype_std[grepl("Control", dm$group, ignore.case = TRUE)] <- "Control"
    missing_subtype_rows <- dm[is.na(dm$Subtype_std) & !grepl("Control", dm$group, ignore.case = TRUE), , drop = FALSE]
    if (nrow(missing_subtype_rows) > 0) {
      log(sprintf("Gene %s: %d tumor sample(s) missing subtype -> removing from design for this gene", g, nrow(missing_subtype_rows)))
      dm <- dm[!(is.na(dm$Subtype_std) & !grepl("Control", dm$group, ignore.case = TRUE)), , drop = FALSE]
    }
    dm <- dm[!duplicated(dm$Sample_ID), , drop = FALSE]
    dm$group_subtyped <- ifelse(is.na(dm$Subtype_std) | dm$Subtype_std %in% c("Control", "Normal"), dm$group, paste0(dm$group, "_", dm$Subtype_std))
    outdf <- dm[, intersect(c("samples", "Sample_ID", "group", "Subtype_std", "group_subtyped"), colnames(dm)), drop = FALSE]
    colnames(outdf)[colnames(outdf) == "Subtype_std"] <- "Subtype"
  } else {
    # no subtype info available; keep original design
    outdf <- d
    outdf$Subtype <- NA_character_
    outdf$group_subtyped <- outdf$group
  }
  out_file <- file.path(design_with_subtypes_dir, paste0(g, "_design_with_subtypes.csv"))
  write.csv(outdf, out_file, row.names = FALSE, quote = TRUE)
  subtype_counts <- as.data.frame(table(outdf$Subtype, useNA = "ifany"), stringsAsFactors = FALSE)
  colnames(subtype_counts) <- c("Subtype", "Count")
  write.table(subtype_counts, file = file.path(design_with_subtypes_dir, paste0(g, "_design_summary.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  if (exists("missing_subtype_rows") && nrow(missing_subtype_rows) > 0) writeLines(missing_subtype_rows$samples, con = file.path(design_with_subtypes_dir, paste0(g, "_missing_samples.txt")))
  log(sprintf(" Wrote augmented design for %s -> %s (kept n=%d rows)", g, out_file, nrow(outdf)))
}

# ---------------- 1) Load inputs ------------------------------------------------
log("1) Loading inputs")
stopifnot(file.exists(ctrl_rds), file.exists(tum_rds), file.exists(gene_map_file))
se_c <- readRDS(ctrl_rds); colnames(se_c) <- file_path_sans_ext(basename(colnames(se_c)))
se_t <- readRDS(tum_rds);  colnames(se_t) <- file_path_sans_ext(basename(colnames(se_t)))
gene_map <- read_excel(gene_map_file)
log(sprintf("  Loaded: control=%d, tumor=%d", ncol(se_c), ncol(se_t)))

# ---- transcript metadata extraction (robust) ----
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

# ---- sample grouping (global Tumor/Control) ----
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
log("  NOTE: Only global filter applied: CPM >= 0.5 in >= 14 samples.")

# ---------------- 4) GLOBAL VOOM (one-time) ---------------------------------
v_all_rds <- file.path(out_dir, "v_all_voom.rds")
design_for_voom <- design_global
design_for_voom <- design_for_voom[colnames(y_all), , drop = FALSE]
rownames(design_for_voom) <- colnames(y_all)

v_all <- NULL
if (file.exists(v_all_rds)) {
  v_all <- readRDS(v_all_rds)
  if (is.null(v_all$E) || is.null(v_all$weights)) v_all <- NULL
}

if (is.null(v_all)) {
  log("  Computing global voom (one-time) on filtered features; samples:", ncol(y_all), "features:", nrow(y_all))
  v_all <- voom(y_all, design_for_voom, plot = FALSE)
  if (is.null(dim(v_all$E))) {
    nr <- nrow(y_all); nc <- ncol(y_all)
    v_all$E <- matrix(v_all$E, nrow = nr, ncol = nc, dimnames = list(rownames(y_all)[seq_len(nr)], colnames(y_all)[seq_len(nc)]))
  }
  if (is.null(v_all$weights) || is.null(dim(v_all$weights))) {
    v_all$weights <- matrix(1, nrow = nrow(v_all$E), ncol = ncol(v_all$E), dimnames = list(rownames(v_all$E), colnames(v_all$E)))
  }
  saveRDS(v_all, v_all_rds); log("  Saved v_all_voom.rds (rows/cols):", nrow(v_all$E), "/", ncol(v_all$E))
}

# ---------------- optional: drop within-patient replicates --------------------
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

# ---------------- robust sample resolution helper -------------------------
resolve_sample_name <- function(sample_str, sample_id, ycolnames) {
  if (!is.na(sample_str) && sample_str %in% ycolnames) return(sample_str)
  if (!is.na(sample_id) && sample_id %in% ycolnames) return(sample_id)
  if (!is.na(sample_id)) {
    pref_matches <- ycolnames[startsWith(ycolnames, sample_id)]
    if (length(pref_matches) == 1) return(pref_matches[1])
    if (length(pref_matches) > 1) {
      pref1 <- pref_matches[grepl("-01[A-Z]?$", pref_matches)]
      if (length(pref1) == 1) return(pref1[1])
      pref2 <- pref_matches[grepl("-01", pref_matches)]
      if (length(pref2) >= 1) return(pref2[1])
      return(pref_matches[1])
    }
  }
  s_no_letter <- sub("([A-Z])$","",sample_str)
  if (!is.na(s_no_letter) && s_no_letter != sample_str) {
    pref_matches2 <- ycolnames[startsWith(ycolnames, s_no_letter)]
    if (length(pref_matches2) >= 1) return(pref_matches2[1])
  }
  ci <- which(tolower(ycolnames) == tolower(sample_str))
  if (length(ci) == 1) return(ycolnames[ci])
  ci2 <- which(tolower(ycolnames) == tolower(sample_id))
  if (length(ci2) == 1) return(ycolnames[ci2])
  return(NA_character_)
}

# ---------------- 7) PER-GENE DEU (serial processing) ---------------------
log("7) Running per-gene DEU (voom on full exon set per-gene sample subsets)")

process_one_gene <- function(g) {
  log(" Processing gene:", g)
  writeLines(as.character(Sys.time()), file.path(status_dir, paste0(g, ".started")))
  ens <- gene_map$ID[gene_map$Gene == g][1]
  exs_gene <- if (!is.na(ens)) grep(sprintf("^%s:", ens), rownames(y_all), value = TRUE) else character(0)
  if (length(exs_gene) < 1) { log("  No exons retained for", g, "- skipping"); file.create(file.path(status_dir, paste0(g, ".done"))); return(NULL) }

  design_file <- file.path(design_with_subtypes_dir, paste0(g, "_design_with_subtypes.csv"))
  if (!file.exists(design_file)) design_file <- file.path(design_dir, paste0(g, ".csv"))
  design_tbl <- if (file.exists(design_file)) { tmp <- read.csv(design_file, stringsAsFactors = FALSE); tmp } else NULL
  if (!is.null(design_tbl) && !"samples" %in% colnames(design_tbl)) { log("  design for", g, "missing 'samples' - ignoring design"); design_tbl <- NULL }
  if (!is.null(design_tbl)) {
    design_tbl$samples <- as.character(design_tbl$samples)
    design_tbl$Sample_ID <- if ("Sample_ID" %in% colnames(design_tbl)) as.character(design_tbl$Sample_ID) else vapply(design_tbl$samples, function(x) { s <- as.character(x); if (grepl("-[0-9]{2}[A-Z]$", s)) return(sub("([0-9]{2})[A-Z]$", "\\1", s)); if (grepl("-[0-9]{2}$", s)) return(s); return(s) }, FUN.VALUE = "")
    design_tbl$Subtype <- if ("Subtype" %in% colnames(design_tbl)) as.character(design_tbl$Subtype) else NA_character_
    design_tbl$group_subtyped <- if ("group_subtyped" %in% colnames(design_tbl)) as.character(design_tbl$group_subtyped) else NA_character_
    ycols <- colnames(y_all)
    design_tbl$matched <- mapply(resolve_sample_name, design_tbl$samples, design_tbl$Sample_ID, MoreArgs = list(ycolnames = ycols), USE.NAMES = FALSE, SIMPLIFY = TRUE)
    n_unmatched <- sum(is.na(design_tbl$matched))
    if (n_unmatched > 0) {
      log(sprintf("  Gene %s: %d design rows could not be matched to filtered samples (they will be ignored): %s", g, n_unmatched, paste(head(design_tbl$samples[is.na(design_tbl$matched)], 10), collapse = ",")))
    }
    design_tbl <- design_tbl[!is.na(design_tbl$matched), , drop = FALSE]
  }

  results_per_contrast <- list()

  # Build contrasts to run (global always included)
  contrasts_to_run <- list()
  contrasts_to_run[["Tumor vs Control"]] <- list(samples = colnames(y_all), group_labels = ifelse(colnames(y_all) %in% colnames(se_t), "Tumor", "Control"), add_subtype = TRUE)

  if (!is.null(design_tbl) && nrow(design_tbl) > 0) {
    loi_rows <- design_tbl$matched[design_tbl$group == "Tumor_LOI"]
    nloi_rows <- design_tbl$matched[design_tbl$group == "Tumor_noLOI"]
    ctrl_no_rows <- design_tbl$matched[design_tbl$group == "Control_noLOI"]
    if (length(loi_rows) >= 2 && length(nloi_rows) >= 2) {
      contrasts_to_run[["Tumor_LOI vs Tumor_noLOI"]] <- list(samples = c(loi_rows, nloi_rows),
                                                            group_labels = c(rep("Tumor_LOI", length(loi_rows)), rep("Tumor_noLOI", length(nloi_rows))),
                                                            add_subtype = TRUE)
      log(sprintf("  Added contrast: Tumor_LOI vs Tumor_noLOI (noLOI=%d, LOI=%d)", length(nloi_rows), length(loi_rows)))
    } else log(sprintf("  Skipping gene-level Tumor_LOI vs Tumor_noLOI (LOI=%d noLOI=%d)", length(loi_rows), length(nloi_rows)))
    if (length(loi_rows) >= 2 && length(ctrl_no_rows) >= 2) {
      contrasts_to_run[["Tumor_LOI vs Control_noLOI"]] <- list(samples = c(loi_rows, ctrl_no_rows),
                                                               group_labels = c(rep("Tumor_LOI", length(loi_rows)), rep("Control_noLOI", length(ctrl_no_rows))),
                                                               add_subtype = TRUE)
      log(sprintf("  Added contrast: Tumor_LOI vs Control_noLOI (LOI=%d ctrl_no=%d)", length(loi_rows), length(ctrl_no_rows)))
    } else log(sprintf("  Skipping gene-level Tumor_LOI vs Control_noLOI (LOI=%d ctrl_no=%d)", length(loi_rows), length(ctrl_no_rows)))
  }

  # Subtype-specific contrasts (no subtype covariate)
  if (!is.null(design_tbl) && nrow(design_tbl) > 0) {
    design_tbl$Subtype[grepl("Control", design_tbl$group, ignore.case = TRUE)] <- "Control"
    present_subtypes <- unique(design_tbl$Subtype[grepl("^Tumor", design_tbl$group, ignore.case = TRUE) & !is.na(design_tbl$Subtype)])
    present_subtypes <- intersect(present_subtypes, SUBTYPES_OF_INTEREST)
    for (st in present_subtypes) {
      tumor_subtype_samples <- design_tbl$matched[grepl("^Tumor", design_tbl$group, ignore.case = TRUE) & design_tbl$Subtype == st]
      tumor_loi_subtype <- design_tbl$matched[design_tbl$group == "Tumor_LOI" & design_tbl$Subtype == st]
      tumor_nloi_subtype <- design_tbl$matched[design_tbl$group == "Tumor_noLOI" & design_tbl$Subtype == st]
      control_all <- intersect(colnames(se_c), colnames(y_all))
      control_no <- design_tbl$matched[design_tbl$group == "Control_noLOI"]

      cname <- paste0("Tumor (", st, ") vs Control")
      if (length(tumor_subtype_samples) >= 2 && length(control_all) >= 2) {
        contrasts_to_run[[cname]] <- list(samples = c(tumor_subtype_samples, control_all),
                                          group_labels = c(rep(paste0("Tumor_", st), length(tumor_subtype_samples)), rep("Control", length(control_all))),
                                          add_subtype = FALSE)
        log(sprintf("  Will attempt contrast: %s -> Tumor n=%d Control n=%d", cname, length(tumor_subtype_samples), length(control_all)))
      } else log(sprintf("  Skip contrast %s (not enough samples) Tumor n=%d Control n=%d", cname, length(tumor_subtype_samples), length(control_all)))

      cname2 <- paste0("Tumor LOI (", st, ") vs Tumor noLOI (", st, ")")
      if (length(tumor_loi_subtype) >= 2 && length(tumor_nloi_subtype) >= 2) {
        contrasts_to_run[[cname2]] <- list(samples = c(tumor_loi_subtype, tumor_nloi_subtype),
                                           group_labels = c(rep(paste0("Tumor_LOI_", st), length(tumor_loi_subtype)), rep(paste0("Tumor_noLOI_", st), length(tumor_nloi_subtype))),
                                           add_subtype = FALSE)
        log(sprintf("  Will attempt contrast: %s -> LOI n=%d noLOI n=%d", cname2, length(tumor_loi_subtype), length(tumor_nloi_subtype)))
      } else log(sprintf("  Skip contrast %s (not enough samples) LOI n=%d noLOI n=%d", cname2, length(tumor_loi_subtype), length(tumor_nloi_subtype)))

      cname3 <- paste0("Tumor LOI (", st, ") vs Control_noLOI")
      if (length(tumor_loi_subtype) >= 2 && length(control_no) >= 2) {
        contrasts_to_run[[cname3]] <- list(samples = c(tumor_loi_subtype, control_no),
                                           group_labels = c(rep(paste0("Tumor_LOI_", st), length(tumor_loi_subtype)), rep("Control_noLOI", length(control_no))),
                                           add_subtype = FALSE)
        log(sprintf("  Will attempt contrast: %s -> LOI n=%d Control_noLOI n=%d", cname3, length(tumor_loi_subtype), length(control_no)))
      } else log(sprintf("  Skip contrast %s (not enough samples) LOI n=%d Control_noLOI n=%d", cname3, length(tumor_loi_subtype), length(control_no)))
    }
  } else {
    log(sprintf("  No design or design had no matched samples for gene %s - will only run global Tumor vs Control", g))
  }

  # Run contrasts
  for (cname in names(contrasts_to_run)) {
    info <- contrasts_to_run[[cname]]
    samp_all <- unique(info$samples)
    samp <- intersect(samp_all, colnames(y_all))
    if (length(samp) < 2) { log("   Skipping contrast (after intersection <2):", cname); next }

    mapping_df <- data.frame(original = info$samples, group_label = info$group_labels, stringsAsFactors = FALSE)
    mapping_df <- mapping_df[mapping_df$original %in% samp, , drop = FALSE]
    if (nrow(mapping_df) < 2) { log("   Not enough mapped rows for contrast:", cname); next }
    ds_rows <- data.frame(samples = mapping_df$original, group = mapping_df$group_label, stringsAsFactors = FALSE)

    ds_rows$Subtype <- NA_character_
    if (!is.null(design_tbl) && nrow(design_tbl) > 0) {
      idx <- match(ds_rows$samples, design_tbl$matched)
      ds_rows$Subtype <- design_tbl$Subtype[idx]
      ds_rows$Subtype[grepl("^Control", ds_rows$group, ignore.case = TRUE)] <- "Control"
    } else {
      ds_rows$Subtype[grepl("^Control", ds_rows$group, ignore.case = TRUE)] <- "Control"
    }

    include_subtype <- FALSE
    if (isTRUE(info$add_subtype) && !is.null(design_tbl) && nrow(design_tbl) > 0) {
      tumor_mask <- grepl("^Tumor", ds_rows$group, ignore.case = TRUE)
      if (any(tumor_mask)) {
        sub_counts <- table(ds_rows$Subtype[tumor_mask], useNA = "no")
        sub_counts <- sub_counts[names(sub_counts) %in% SUBTYPES_OF_INTEREST]
        if (length(sub_counts) > 0) {
          n_levels_ge <- sum(sub_counts >= MIN_SUBTYPE_SAMPLES)
          if (n_levels_ge >= MIN_SUBTYPE_LEVELS) include_subtype <- TRUE
        }
      }
    }

    ds_rows$group <- as.character(ds_rows$group)
    ds_rows$group <- factor(ds_rows$group, levels = unique(ds_rows$group))
    if (include_subtype) {
      ds_rows$Subtype_f <- ifelse(is.na(ds_rows$Subtype), "Unknown", ds_rows$Subtype)
      ds_rows$Subtype_f <- factor(ds_rows$Subtype_f, levels = unique(ds_rows$Subtype_f))
      design_g <- model.matrix(~0 + group + Subtype_f, data = ds_rows)
    } else {
      design_g <- model.matrix(~0 + group, data = ds_rows)
    }
    colnames(design_g) <- make.names(colnames(design_g))
    rownames(design_g) <- ds_rows$samples
    design_g <- design_g[match(samp, rownames(design_g)), , drop = FALSE]
    rownames(design_g) <- samp

    # Determine pos/neg labels (several pre-set patterns handled, including subtype variants)
    pos_label <- NULL; neg_label <- NULL
    if (cname == "Tumor vs Control") { pos_label <- "Tumor"; neg_label <- "Control" }
    else if (cname == "Tumor_LOI vs Tumor_noLOI") { pos_label <- "Tumor_LOI"; neg_label <- "Tumor_noLOI" }
    else if (cname == "Tumor_LOI vs Control_noLOI") { pos_label <- "Tumor_LOI"; neg_label <- "Control_noLOI" }
    else if (grepl("^Tumor \\(", cname)) {
      st <- sub("^Tumor \\(([^)]+)\\).*", "\\1", cname)
      pos_label <- paste0("Tumor_", st); neg_label <- "Control"
    } else if (grepl("^Tumor LOI \\(", cname)) {
      st <- sub("^Tumor LOI \\(([^)]+)\\).*", "\\1", cname)
      pos_label <- paste0("Tumor_LOI_", st)
      # Bugfix: if contrast mentions "vs Control_noLOI", ensure neg_label becomes Control_noLOI
      if (grepl("Control_noLOI", cname, fixed = TRUE)) {
        neg_label <- "Control_noLOI"
      } else {
        neg_label <- paste0("Tumor_noLOI_", st)
      }
    } else {
      gls <- as.character(levels(ds_rows$group))
      if (length(gls) >= 2) { pos_label <- gls[1]; neg_label <- gls[2] } else { log("   Not enough groups for contrast:", cname); next }
    }

    pos_col_candidates <- grep(paste0("^group", make.names(pos_label)), colnames(design_g), value = TRUE)
    neg_col_candidates <- grep(paste0("^group", make.names(neg_label)), colnames(design_g), value = TRUE)
    if (length(pos_col_candidates) == 0) pos_col_candidates <- grep(make.names(pos_label), colnames(design_g), value = TRUE)
    if (length(neg_col_candidates) == 0) neg_col_candidates <- grep(make.names(neg_label), colnames(design_g), value = TRUE)
    if (length(pos_col_candidates) == 0 || length(neg_col_candidates) == 0) {
      log("   Could not detect pos/neg design columns for contrast:", cname, "-> skipping (pos candidates:", paste(pos_col_candidates, collapse = ","), " neg candidates:", paste(neg_col_candidates, collapse = ","),")")
      next
    }
    gp_pos_col <- pos_col_candidates[1]; gp_neg_col <- neg_col_candidates[1]

    # subset y_all to samp (keep all exons)
    y_sub <- y_all[, samp, keep.lib.sizes = TRUE]

    reuse_v_all <- FALSE
    if (length(samp) == ncol(y_all) && all(sort(samp) == sort(colnames(y_all)))) reuse_v_all <- TRUE
    if (reuse_v_all) {
      cidx <- match(colnames(y_sub), colnames(v_all$E))
      if (any(is.na(cidx))) compute_voom_here <- TRUE else {
        compute_voom_here <- FALSE
        v_sub <- v_all
        v_sub$E <- v_all$E[, cidx, drop = FALSE]
        v_sub$weights <- v_all$weights[, cidx, drop = FALSE]
        colnames(v_sub$E) <- colnames(y_sub); colnames(v_sub$weights) <- colnames(y_sub)
        rownames(v_sub$E) <- rownames(y_sub); rownames(v_sub$weights) <- rownames(y_sub)
        log("   Reused global v_all for contrast:", cname)
      }
    } else compute_voom_here <- TRUE

    if (compute_voom_here) {
      v_sub <- tryCatch({ voom(y_sub, design_g, plot = FALSE) }, error = function(e) { log("   voom error for gene", g, "contrast", cname, ":", conditionMessage(e)); NULL })
      if (is.null(v_sub)) next
      if (is.null(dim(v_sub$E))) {
        nr <- nrow(y_sub); nc <- ncol(y_sub)
        v_sub$E <- matrix(v_sub$E, nrow = nr, ncol = nc, dimnames = list(rownames(y_sub)[seq_len(nr)], colnames(y_sub)[seq_len(nc)]))
      }
      if (is.null(v_sub$weights) || is.null(dim(v_sub$weights))) v_sub$weights <- matrix(1, nrow = nrow(v_sub$E), ncol = ncol(v_sub$E), dimnames = list(rownames(v_sub$E), colnames(v_sub$E)))
      log("   Computed voom for gene", g, "contrast", cname, "(samples:", ncol(v_sub$E), ")")
    }

    fit1 <- tryCatch({ lmFit(v_sub, design_g) }, error = function(e) { log("   lmFit error:", conditionMessage(e)); NULL })
    if (is.null(fit1)) next

    contrast_expr <- paste0(gp_pos_col, " - ", gp_neg_col)
    contrast_mat <- tryCatch({ makeContrasts(contrasts = contrast_expr, levels = design_g) }, error = function(e) { log("   makeContrasts error:", conditionMessage(e)); NULL })
    if (is.null(contrast_mat)) { log("   Failed to build contrast matrix for", cname); next }
    fit2 <- tryCatch({ contrasts.fit(fit1, contrast_mat) }, error = function(e) { log("   contrasts.fit error:", conditionMessage(e)); NULL })
    if (is.null(fit2)) next
    fit2 <- tryCatch({ eBayes(fit2) }, error = function(e) { log("   eBayes warning:", conditionMessage(e)); fit2 })

    geneid_vec <- if (!is.null(y_sub$genes) && "geneid" %in% colnames(y_sub$genes)) y_sub$genes$geneid else sub(":.*", "", rownames(v_sub$E))
    if (length(geneid_vec) != nrow(v_sub$E)) geneid_vec <- sub(":.*", "", rownames(v_sub$E))
    fit2$genes <- data.frame(GeneID = geneid_vec, exonID = rownames(v_sub$E), stringsAsFactors = FALSE)

    ex_all <- tryCatch({ diffSplice(fit2, geneid = "GeneID") }, error = function(e) { log("   diffSplice error:", conditionMessage(e)); NULL })
    if (is.null(ex_all)) { log("   diffSplice returned NULL for", g, cname); next }
    tbl_all <- tryCatch({ topSplice(ex_all, coef = 1, test = "t", number = Inf) }, error = function(e) { log("   topSplice error:", conditionMessage(e)); NULL })
    if (is.null(tbl_all) || nrow(tbl_all) == 0) { log("   topSplice returned no rows for", g, cname); next }
    tbl_all <- as.data.frame(tbl_all)

    # Robust subsetting to gene
    gene_rows <- data.frame()
    gene_col_names <- intersect(c("GeneID", "geneid", "gene"), colnames(tbl_all))
    if (length(gene_col_names) > 0 && !is.na(ens)) {
      for (gn in gene_col_names) {
        if (ens %in% tbl_all[[gn]]) { gene_rows <- tbl_all[tbl_all[[gn]] == ens, , drop = FALSE]; break }
      }
    }
    if (nrow(gene_rows) == 0 && length(exs_gene) > 0) {
      common <- intersect(rownames(tbl_all), exs_gene)
      if (length(common) > 0) gene_rows <- tbl_all[rownames(tbl_all) %in% common, , drop = FALSE]
    }
    if (nrow(gene_rows) == 0 && !is.na(ens)) {
      pat <- paste0("^", ens, ":"); matches <- grepl(pat, rownames(tbl_all))
      if (any(matches, na.rm = TRUE)) gene_rows <- tbl_all[matches, , drop = FALSE]
    }
    if (nrow(gene_rows) == 0) {
      exoncols <- intersect(c("exonID", "exonid", "exon", "exon_id"), colnames(tbl_all))
      if (length(exoncols) > 0) {
        for (ec in exoncols) {
          common2 <- intersect(tbl_all[[ec]], exs_gene)
          if (length(common2) > 0) { gene_rows <- tbl_all[tbl_all[[ec]] %in% common2, , drop = FALSE]; break }
        }
      }
    }
    if (nrow(gene_rows) == 0) { log("   After subsetting, no rows for gene", g, "in contrast", cname); next }

    exon_col <- if ("exonID" %in% colnames(gene_rows)) "exonID" else if ("exonid" %in% colnames(gene_rows)) "exonid" else if ("exon" %in% colnames(gene_rows)) "exon" else NA_character_
    if (is.na(exon_col)) gene_rows$exon_id <- rownames(gene_rows) else gene_rows$exon_id <- as.character(gene_rows[[exon_col]])

    # Filter for requested exon numbers if provided
    keep_nums <- KEEP_EXONS[[g]]
    if (!is.null(keep_nums)) {
      gene_rows$exon_num <- vapply(gene_rows$exon_id, get_exon_number, FUN.VALUE = character(1))
      gene_rows <- gene_rows[gene_rows$exon_num %in% keep_nums, , drop = FALSE]
      if (nrow(gene_rows) == 0) { log("   After exon filtering, no requested exons for", g, "in contrast", cname); next }
    }

    gene_rows$gene <- g; gene_rows$contrast <- cname
    gene_rows$tx_name <- tx_meta_all$tx_name[match(gene_rows$exon_id, rownames(tx_meta_all))]

    counts_tbl <- as.data.frame(table(ds_rows$group), stringsAsFactors = FALSE)
    rownames(counts_tbl) <- counts_tbl$Var1
    pos_label_use <- pos_label; neg_label_use <- neg_label
    if (!(pos_label_use %in% ds_rows$group)) {
      poss <- ds_rows$group[grepl(make.names(pos_label), make.names(ds_rows$group), fixed = TRUE)]
      if (length(poss) > 0) pos_label_use <- poss[1]
    }
    if (!(neg_label_use %in% ds_rows$group)) {
      negs <- ds_rows$group[grepl(make.names(neg_label), make.names(ds_rows$group), fixed = TRUE)]
      if (length(negs) > 0) neg_label_use <- negs[1]
    }
    groups_present <- unique(ds_rows$group)
    order_groups <- c()
    if (!is.null(pos_label_use) && pos_label_use %in% groups_present) order_groups <- c(order_groups, pos_label_use)
    if (!is.null(neg_label_use) && neg_label_use %in% groups_present && !(neg_label_use %in% order_groups)) order_groups <- c(order_groups, neg_label_use)
    order_groups <- unique(c(order_groups, groups_present))
    samp_counts_vec <- sapply(order_groups, function(x) { if (x %in% rownames(counts_tbl)) counts_tbl[x, "Freq"] else 0 })
    gene_rows$samples_per_group <- paste0(paste0(order_groups, ":", samp_counts_vec), collapse = ";")

    gene_rds <- file.path(out_dir, sprintf("%s_DEU_results.rds", g))
    if (!file.exists(gene_rds)) saveRDS(list(gene_rows), gene_rds) else {
      existing <- readRDS(gene_rds); if (!is.list(existing)) existing <- list(existing)
      existing[[cname]] <- gene_rows
      saveRDS(existing, gene_rds)
    }

    results_per_contrast[[cname]] <- gene_rows
    log(sprintf("   Collected %d requested exon rows for %s (%s)", nrow(gene_rows), g, cname))
  } # contrasts

  file.create(file.path(status_dir, paste0(g, ".done")))
  startedfile <- file.path(status_dir, paste0(g, ".started")); if (file.exists(startedfile)) file.remove(startedfile)
  log("  Finished gene:", g)
  results_per_contrast
}

# run serially to avoid parallel env issues
res_list <- lapply(candidate_genes, process_one_gene)
names(res_list) <- candidate_genes

# Merge into per_contrast_results
per_contrast_results <- list()
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
  df_all$exon_num <- vapply(df_all$exon_id, get_exon_number, FUN.VALUE = character(1))
  df_all <- df_all %>% rowwise() %>% filter({
    allowed <- KEEP_EXONS[[gene]]
    if (is.null(allowed)) TRUE else exon_num %in% allowed
  }) %>% ungroup()

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
    base_width <- 8
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
           axis.text.x = element_markdown(angle = 45, hjust = 1, size = rel(0.9)),
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
    ggsave(plot_file, p, width = plot_width, height = 4.5, dpi = 300)
    log("  Saved dodge (single per-contrast plot, equal-width bars):", plot_file)
  }

  df_save <- df_all %>% select(-any_of(c("geneid", "exonid", "exon_num")))
  writeData(deu_wb, sheet_name, df_save)
}
saveWorkbook(deu_wb, file.path(out_dir, "DEU_limma_results_per_contrast.xlsx"), overwrite = TRUE)
log("  Wrote DEU_limma_results_per_contrast.xlsx")

# ---------------- 9) Violin plots + final summaries ------------------------
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

  design_csv_aug <- file.path(design_with_subtypes_dir, paste0(g, "_design_with_subtypes.csv"))
  if (file.exists(design_csv_aug)) dtmp <- read.csv(design_csv_aug, stringsAsFactors = FALSE) else dtmp <- read.csv(file.path(design_dir, paste0(g, ".csv")), stringsAsFactors = FALSE)

  if (!all(c("samples", "group") %in% colnames(dtmp))) {
    samp_use <- colnames(y_all)
    sample_map <- data.frame(sample = samp_use, parent = sample_group_global$group[match(samp_use, sample_group_global$sample)], subgroup = NA_character_, stringsAsFactors = FALSE)
  } else {
    dtmp$samples <- as.character(dtmp$samples)
    dtmp$Sample_ID <- if ("Sample_ID" %in% colnames(dtmp)) as.character(dtmp$Sample_ID) else vapply(dtmp$samples, function(x) { s <- as.character(x); if (grepl("-[0-9]{2}[A-Z]$", s)) return(sub("([0-9]{2})[A-Z]$", "\\1", s)); if (grepl("-[0-9]{2}$", s)) return(s); return(s) }, FUN.VALUE = "")
    ycols <- colnames(y_all)
    dtmp$matched <- mapply(resolve_sample_name, dtmp$samples, dtmp$Sample_ID, MoreArgs = list(ycolnames = ycols), USE.NAMES = FALSE, SIMPLIFY = TRUE)
    samp_use <- intersect(dtmp$matched, colnames(y_all))
    if (length(samp_use) == 0) {
      log("  design csv samples for", g, "not found in filtered samples -> fallback to global grouping")
      samp_use <- colnames(y_all)
      sample_map <- data.frame(sample = samp_use, parent = sample_group_global$group[match(samp_use, sample_group_global$sample)], subgroup = NA_character_, stringsAsFactors = FALSE)
    } else {
      grp_vec <- if ("group_subtyped" %in% colnames(dtmp)) dtmp$group_subtyped[match(samp_use, dtmp$matched)] else dtmp$group[match(samp_use, dtmp$matched)]
      parent_vec <- ifelse(grepl("Tumor", grp_vec, ignore.case = TRUE), "Tumor", ifelse(grepl("Control", grp_vec, ignore.case = TRUE), "Control", sample_group_global$group[match(samp_use, sample_group_global$sample)]))
      sample_map <- data.frame(sample = samp_use, parent = parent_vec, subgroup = grp_vec, stringsAsFactors = FALSE)
      log("  Using design.csv for", g, "-> samples:", length(samp_use))
    }
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
    if (length(exs) > 0) gene_vals <- rowSums(cpm_filtered[exs, , drop = FALSE])
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
        if (length(exs) > 0) gene_vals_local <- rowSums(cpm_filtered[exs, , drop = FALSE]) else gene_vals_local <- rep(NA_real_, ncol(y_all))
        vals <- gene_vals_local[match(samp_subset, colnames(y_all))]
      } else vals <- gene_vals[match(samp_subset, colnames(y_all))]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(Gene = g, Subgroup = grp, Mean = mean(vals, na.rm = TRUE), Median = median(vals, na.rm = TRUE), N = length(vals), stringsAsFactors = FALSE)
    }
  }
}
summary_df <- if (length(summary_rows) > 0) bind_rows(summary_rows) else data.frame()
wbStats <- createWorkbook()
addWorksheet(wbStats, "PerGene_Subgroup_Summary"); if (nrow(summary_df) > 0) writeData(wbStats, "PerGene_Subgroup_Summary", summary_df)
addWorksheet(wbStats, "FilterInfo"); writeData(wbStats, "FilterInfo", data.frame(InitialBins = nrow(counts_all), RetainedBins = nrow(y_all), RemovedBins = nrow(counts_all) - nrow(y_all)))
addWorksheet(wbStats, "SampleStats"); writeData(wbStats, "SampleStats", data.frame(sample = colnames(y_all)))
saveWorkbook(wbStats, file.path(out_dir, "summary_stats_and_sample_info.xlsx"), overwrite = TRUE)
log("  Saved summary_stats_and_sample_info.xlsx")

log("ALL DONE. Outputs written to:", out_dir)
