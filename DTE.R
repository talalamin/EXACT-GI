#!/usr/bin/env Rscript
# Transcript-level DGE pipeline for candidate imprinted genes
# #TALAL AMIN: BIOBIX (UGENT)
# Edit the path variables below before running.

# ------------------------- LOAD LIBRARIES -----------------------------------
suppressPackageStartupMessages({
  library(tximport)
  library(edgeR)
  library(limma)
  library(openxlsx)
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})
message("Libraries loaded")

# ------------------------- SETTINGS -----------------------------------------
# Toggleable behavior
RUN_PVAL_HIST        <- TRUE   # compute & save transcript-level p-value histograms #Optional
RUN_VOOMPLOTS        <- TRUE   # save voom diagnostic plots per gene #Optional
RUN_PER_GENE_PLOTS   <- TRUE   # save per-gene violin & histogram plots #Optional
USE_GLOBAL_CPM_FILTER <- TRUE  # Use single global CPM filter (if FALSE, falls back to filterByExpr) #Optional 
GLOBAL_CPM           <- 1      # CPM threshold for global filter (if USE_GLOBAL_CPM_FILTER = TRUE)
GLOBAL_MIN_SAMPLES   <- 14     # minimum number of samples with CPM >= GLOBAL_CPM to keep feature
RUN_DUPLICATE_REMOVAL <- TRUE  # Remove per-patient within-group replicates keeping highest libsize #Optional

# ------------------------- I/O / PATHS -------------------------------------
# Replace these placeholders with your actual directories or relative paths.
tumor_dir   <- "path/to/kallisto/tumor/"
control_dir <- "path/to/kallisto/control/"
design_dir  <- "path/to/design_csvs/"    # per-gene design CSVs: Gene.csv with columns 'samples' and 'group'
out_dir     <- "path/to/output/Transcript_DGE/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Path to tx2gene mapping (RDS). Should be an R data.frame with columns: tx, gene_id, gene_name
tx2g_rds_path <- "path/to/tx2gene.rds"

# ------------------------- CANDIDATE GENES & CONTRASTS -----------------------
GI_genes       <- c("HM13", "MEST", "IGF2", "ZNF331")  #Provide Candidate gene names (like provided examples)
contrast_names <- c("Tumor_LOI_vs_Tumor_noLOI",
                    "Tumor_LOI_vs_Control_noLOI",
                    "Tumor_vs_Control")

clean_label <- function(lbl) gsub("_", " ", lbl)

# ------------------------- FILES & CACHING ---------------------------------
rds_dir <- file.path(out_dir, "RDS_checkpoints")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
txi_rds   <- file.path(rds_dir, "txi_scaledTPM_transcript.rds")
yfull_rds <- file.path(rds_dir, "yFull_transcript.rds")

# ------------------------- STEP 1: list kallisto files -----------------------
tumor_files   <- list.files(tumor_dir,   "abundance.h5$", full.names = TRUE, recursive = TRUE)
control_files <- list.files(control_dir, "abundance.h5$", full.names = TRUE, recursive = TRUE)
all_files     <- c(tumor_files, control_files)
samples       <- basename(dirname(all_files))
names(all_files) <- samples

# ------------------------- STEP 2: tx2gene ----------------------------------
if (!file.exists(tx2g_rds_path)) stop("tx2gene RDS not found at: ", tx2g_rds_path)
tx2g <- readRDS(tx2g_rds_path)
if (!all(c("tx", "gene_id", "gene_name") %in% colnames(tx2g))) {
  # try flexible handling if columns have different names
  colnames(tx2g)[1:3] <- c("tx", "gene_id", "gene_name")
}
tx2g$tx <- sub("\\..*$", "", tx2g$tx)   # strip version suffix if present
tx2g <- tx2g[!duplicated(tx2g$tx), ]

# ------------------------- STEP 3: tximport / build DGEList ------------------
if (file.exists(yfull_rds)) {
  message("Loading cached yFull from ", yfull_rds)
  yFull <- readRDS(yfull_rds)
} else {
  if (file.exists(txi_rds)) {
    message("Loading cached txi from ", txi_rds)
    txi <- readRDS(txi_rds)
  } else {
    message("Running tximport with scaledTPM")
    txi <- tximport(
      files = all_files,
      type = "kallisto",
      txIn = TRUE,
      txOut = TRUE,
      countsFromAbundance = "scaledTPM",
      ignoreTxVersion = TRUE
    )
    saveRDS(txi, txi_rds)
    message("Saved txi to ", txi_rds)
  }

  yFull <- DGEList(counts = txi$counts)
  tr_ids <- rownames(yFull$counts)
  tr_nov <- sub("\\..*$", "", tr_ids)
  idx <- match(tr_nov, tx2g$tx)
  yFull$genes <- data.frame(
    TranscriptID = tr_ids,
    GeneID       = tx2g$gene_id[idx],
    GeneName     = tx2g$gene_name[idx],
    stringsAsFactors = FALSE
  )
  yFull <- calcNormFactors(yFull)
  saveRDS(yFull, yfull_rds)
  message("Saved yFull to ", yfull_rds)
}
message("Built transcript-level DGEList")

# ------------------------- STEP 4: Derive grouping Tumor vs Control -----------
samples_all <- colnames(yFull)
tumor_samples   <- unique(basename(dirname(tumor_files)))
control_samples <- unique(basename(dirname(control_files)))
group2_vec <- ifelse(samples_all %in% tumor_samples, "Tumor",
                     ifelse(samples_all %in% control_samples, "Control", NA))
if (any(is.na(group2_vec))) {
  warning("Dropping samples not found in tumor/control dirs; verify your directory layout")
}
keep_samps <- samples_all[!is.na(group2_vec)]
group2_vec  <- group2_vec[!is.na(group2_vec)]

yFull <- yFull[, keep_samps, keep.lib.sizes = FALSE]
yFull$samples <- yFull$samples[keep_samps, , drop = FALSE]
yFull$samples$group2 <- factor(group2_vec, levels = c("Control", "Tumor"))
yFull <- calcNormFactors(yFull)
message("Defined Tumor vs Control groups")

# ------------------------- STEP 5: Global filter -----------------------------
if (USE_GLOBAL_CPM_FILTER) {
  message("Applying global CPM filter: CPM >= ", GLOBAL_CPM, " in >= ", GLOBAL_MIN_SAMPLES, " samples")
  cpm_full <- cpm(yFull, log = FALSE)
  keep_filter <- rowSums(cpm_full >= GLOBAL_CPM) >= GLOBAL_MIN_SAMPLES
  message("Global CPM filter kept ", sum(keep_filter), " of ", nrow(yFull), " transcripts")
} else {
  message("Using filterByExpr fallback")
  design_all <- model.matrix(~0 + group2, data = yFull$samples)
  colnames(design_all) <- levels(yFull$samples$group2)
  keep_filter <- filterByExpr(yFull, design_all)
  message("filterByExpr kept ", sum(keep_filter), " of ", nrow(yFull), " transcripts")
}

yFilt <- yFull[keep_filter, , keep.lib.sizes = FALSE]
yFilt$samples <- yFull$samples
yFilt$genes   <- yFull$genes[keep_filter, ]
yFilt <- calcNormFactors(yFilt)

# ------------------------- STEP 6: Global p-value histograms -----------------
if (RUN_PVAL_HIST) {
  message("Computing transcript-level p-values before & after filtering (Tumor vs Control)")

  design_before <- model.matrix(~0 + group2, data = yFull$samples)
  colnames(design_before) <- levels(yFull$samples$group2)

  vFull <- voom(yFull, design_before, plot = FALSE)
  fitFull <- lmFit(vFull, design_before)
  cm_all <- makeContrasts(Tumor - Control, levels = design_before)
  fitFull2 <- contrasts.fit(fitFull, cm_all) %>% eBayes()
  pvals_before <- fitFull2$p.value[, 1]

  design_after <- design_before[colnames(yFilt), , drop = FALSE]
  vFilt <- voom(yFilt, design_after, plot = FALSE)
  fitFilt <- lmFit(vFilt, design_after)
  fitFilt2 <- contrasts.fit(fitFilt, cm_all) %>% eBayes()
  pvals_after <- fitFilt2$p.value[, 1]

  df_p <- data.frame(
    pvalue = c(pvals_before, pvals_after),
    status = rep(c("Before filter", "After filter"),
                 times = c(length(pvals_before), length(pvals_after)))
  )

  # log-scale y histogram
  p_hist_log <- ggplot(df_p, aes(x = pvalue, fill = status)) +
    geom_histogram(data = subset(df_p, status == "Before filter"), bins = 50, alpha = 0.45, position = "identity", boundary = 0) +
    geom_histogram(data = subset(df_p, status == "After filter"),  bins = 50, alpha = 0.45, position = "identity", boundary = 0) +
    scale_fill_viridis_d() +
    scale_y_log10() +
    labs(title = "Transcript-level p-value distributions (Tumor vs Control)",
         subtitle = "Before and after global CPM filter",
         x = "p-value", y = "Count (log10 scale)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14), legend.position = "top")

  png(file.path(out_dir, "transcript_pvalue_hist_before_after_filter_logY.png"),
      width = 1000, height = 600, res = 150)
  print(p_hist_log); dev.off()
  message("Saved transcript-level p-value histogram (log Y)")

  # linear y histogram
  p_hist_lin <- ggplot(df_p, aes(x = pvalue, fill = status)) +
    geom_histogram(data = subset(df_p, status == "Before filter"), bins = 50, alpha = 0.45, position = "identity", boundary = 0) +
    geom_histogram(data = subset(df_p, status == "After filter"),  bins = 50, alpha = 0.45, position = "identity", boundary = 0) +
    scale_fill_viridis_d() +
    labs(title = "Transcript-level p-value distributions (Tumor vs Control)",
         subtitle = "Before and after global CPM filter",
         x = "p-value", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14), legend.position = "top")

  png(file.path(out_dir, "transcript_pvalue_hist_before_after_filter_linearY.png"),
      width = 1000, height = 600, res = 150)
  print(p_hist_lin); dev.off()
  message("Saved transcript-level p-value histogram (linear Y)")
}

# ------------------------- STEP 7: Read design CSVs per gene ----------------
design_lists <- lapply(GI_genes, function(g) {
  f <- file.path(design_dir, paste0(g, ".csv"))
  if (!file.exists(f)) return(NULL)
  si <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("samples", "group") %in% colnames(si))) return(NULL)
  si <- subset(si, samples %in% colnames(yFilt))
  if (nrow(si) < 2) return(NULL)
  valid_lvls <- c("Tumor_LOI", "Tumor_noLOI", "Control_LOI", "Control_noLOI")
  si <- subset(si, group %in% valid_lvls)
  if (nrow(si) < 2) return(NULL)
  list(gene = g, samp_info = si)
})
design_lists <- Filter(Negate(is.null), design_lists)
message("Valid genes for transcript-level analysis: ", paste(sapply(design_lists, `[[`, "gene"), collapse = ", "))

# ------------------------- STEP 8: Per-gene transcript-level analysis ------
results_transcripts <- setNames(lapply(contrast_names, function(x) list()), contrast_names)
summary_stats_trans <- list()

# Helper: format p-values with superscript exponent when < 0.01
format_pval_sup <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (!is.numeric(p)) return(as.character(p))
  if (p < 0.01) {
    s <- formatC(p, format = "e", digits = 2)
    parts <- strsplit(s, "e")[[1]]
    mant <- parts[1]
    expn <- parts[2]
    if (!grepl("^[+-]", expn)) expn <- paste0("+", expn)
    sign <- substr(expn, 1, 1)
    digits <- substr(expn, 2, nchar(expn))
    digits <- sprintf("%02d", as.integer(digits))
    sup_map <- c("0" = "\u2070", "1" = "\u00B9", "2" = "\u00B2", "3" = "\u00B3",
                 "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077",
                 "8" = "\u2078", "9" = "\u2079")
    sup_minus <- "\u207B"
    sup_digits <- paste0(sapply(strsplit(digits, "")[[1]], function(d) sup_map[[d]]), collapse = "")
    sup_sign <- if (sign == "-") sup_minus else ""
    out <- paste0(mant, " x 10", sup_sign, sup_digits)
    return(out)
  } else {
    return(formatC(p, format = "f", digits = 3))
  }
}

for (dl in design_lists) {
  g <- dl$gene
  message("Processing gene ", g)
  si <- dl$samp_info
  vs_all <- intersect(si$samples, colnames(yFilt))
  if (length(vs_all) < 2) {
    message("  Skipping ", g, ": fewer than 2 samples")
    next
  }

  # Duplicate removal by highest global libsize (if requested)
  valid_samples <- vs_all
  if (RUN_DUPLICATE_REMOVAL) {
    pid <- substr(valid_samples, 1, 12)
    df_dup <- data.frame(sample = valid_samples,
                         patient = pid,
                         subgroup = si$group[match(valid_samples, si$samples)],
                         stringsAsFactors = FALSE)
    keys <- paste0(df_dup$patient, "::", df_dup$subgroup)
    dup_keys <- names(table(keys))[table(keys) > 1]
    if (length(dup_keys) > 0) {
      global_libsizes <- yFilt$samples$lib.size
      names(global_libsizes) <- rownames(yFilt$samples)
      keep_flags <- rep(TRUE, length(valid_samples)); names(keep_flags) <- valid_samples
      for (k in dup_keys) {
        dup_samples <- df_dup$sample[keys == k]
        glib <- global_libsizes[dup_samples]
        glib[is.na(glib)] <- 0
        keep_sample <- names(glib)[which.max(glib)]
        remove_samples <- setdiff(dup_samples, keep_sample)
        if (length(remove_samples) > 0) keep_flags[remove_samples] <- FALSE
      }
      kept <- names(keep_flags)[keep_flags]
      dropped_n <- sum(!keep_flags)
      if (dropped_n > 0) message("  Dropped ", dropped_n, " within-patient replicate(s) for gene ", g)
      valid_samples <- kept
    }
  }

  vs_all <- valid_samples
  if (length(vs_all) < 2) {
    message("  Skipping ", g, ": <2 samples after deduplication")
    next
  }

  # Subset filtered DGEList to gene samples
  yG <- yFilt[, vs_all]
  yG$samples <- yFilt$samples[vs_all, , drop = FALSE]
  yG$samples$subgroup <- factor(si$group[match(vs_all, si$samples)])
  grp_sub <- yG$samples$subgroup

  # Raw CPM matrix
  raw_cpm_mat <- cpm(yG, normalized.lib.sizes = TRUE)

  # Subgroup design
  design_sub <- model.matrix(~0 + subgroup, data = yG$samples)
  colnames(design_sub) <- levels(grp_sub)

  # Filter transcripts with zero total counts across selected samples
  keep_rows <- rowSums(yG$counts) > 0
  yG2 <- yG[keep_rows, , keep.lib.sizes = FALSE]
  yG2 <- calcNormFactors(yG2)
  raw_cpm_mat2 <- raw_cpm_mat[keep_rows, , drop = FALSE]

  # Identify transcripts belonging to current gene
  gene_names_vec <- yG2$genes$GeneName
  trans_ids_all  <- yG2$genes$TranscriptID
  is_trans_g <- !is.na(gene_names_vec) & (gene_names_vec == g)
  trans_ids <- trans_ids_all[is_trans_g]
  trans_ids <- intersect(trans_ids, rownames(raw_cpm_mat2))
  if (length(trans_ids) == 0) {
    message("  No transcripts of ", g, " after filtering; skipping")
    next
  }
  if (nrow(yG2) < 10) {
    message("  Too few transcripts overall for voom for ", g, "; skipping")
    next
  }

  # Voom diagnostic plot (subgroup)
  voom_sub_plot <- file.path(out_dir, paste0(g, "_transcript_voom_subgroup.png"))
  if (RUN_VOOMPLOTS) {
    png(voom_sub_plot, width = 1200, height = 900, res = 150)
    vGsub <- tryCatch({
      voom(yG2, design_sub, plot = TRUE)
    }, error = function(e) {
      dev.off(); message("  Voom error for ", g, ": ", e$message); NULL
    })
    dev.off()
    if (is.null(vGsub)) {
      message("  Skipping ", g, " due to voom error"); next
    }
  } else {
    vGsub <- tryCatch({ voom(yG2, design_sub, plot = FALSE) }, error = function(e) { message("  voom error: ", e$message); NULL })
    if (is.null(vGsub)) next
  }
  message("  Created voom (subgroup) for ", g)
  voom_mat <- vGsub$E

  # Subset matrices to transcripts of this gene
  raw_cpm_sub <- raw_cpm_mat2[trans_ids, , drop = FALSE]
  logcpm_yG2 <- cpm(yG2, log = TRUE)
  logcpm_sub <- logcpm_yG2[trans_ids, , drop = FALSE]
  log2_raw_sub <- log2(raw_cpm_sub + 1)
  voom_mat_sub <- voom_mat[rownames(voom_mat) %in% trans_ids, , drop = FALSE]
  voom_mat_sub <- voom_mat_sub[trans_ids, , drop = FALSE]

  # Summary stats per transcript x subgroup
  df_vals_all <- tibble::tibble(
    Gene         = g,
    TranscriptID = rep(trans_ids, times = ncol(raw_cpm_sub)),
    Sample       = rep(colnames(raw_cpm_sub), each = length(trans_ids)),
    Subgroup     = rep(yG2$samples$subgroup, each = length(trans_ids)),
    raw_CPM      = as.vector(raw_cpm_sub),
    log2CPM_raw  = as.vector(log2_raw_sub),
    voom_log2CPM = as.vector(voom_mat_sub)
  )
  df_stats <- df_vals_all %>%
    group_by(Gene, TranscriptID, Subgroup) %>%
    summarize(
      N = n(),
      mean_raw_CPM       = mean(raw_CPM, na.rm = TRUE),
      median_raw_CPM     = median(raw_CPM, na.rm = TRUE),
      mean_log2CPM_raw   = mean(log2CPM_raw, na.rm = TRUE),
      median_log2CPM_raw = median(log2CPM_raw, na.rm = TRUE),
      mean_voom_log2CPM  = mean(voom_log2CPM, na.rm = TRUE),
      median_voom_log2CPM= median(voom_log2CPM, na.rm = TRUE),
      .groups = "drop"
    )
  summary_stats_trans[[g]] <- df_stats

  # ----------------- Per-gene plotting (all transcripts in one figure) ----
  if (RUN_PER_GENE_PLOTS) {
    outdir_gene <- file.path(out_dir, g); dir.create(outdir_gene, recursive = TRUE, showWarnings = FALSE)

    # Violin plots of log-CPM by subgroup per transcript
    parent_df <- as.data.frame(t(logcpm_sub))
    parent_df$sample <- rownames(parent_df)
    parent_long <- pivot_longer(parent_df, -sample, names_to = "TranscriptID", values_to = "expr")
    samp_map <- data.frame(sample = colnames(logcpm_sub), subgroup = yG2$samples$subgroup, stringsAsFactors = FALSE)
    parent_long <- left_join(parent_long, samp_map, by = "sample")
    parent_long$grp_label <- parent_long$subgroup
    wanted <- c("Control", "Tumor", "Tumor LOI", "Tumor noLOI", "Control noLOI")
    parent_long$grp_label <- factor(gsub("_", " ", parent_long$grp_label), levels = wanted)
    plot_df <- parent_long[!is.na(parent_long$grp_label) & parent_long$grp_label %in% wanted, , drop = FALSE]
    if (nrow(plot_df) > 0) {
      plot_df$TranscriptID <- factor(plot_df$TranscriptID, levels = unique(plot_df$TranscriptID))
      title_main <- paste0(g, " transcript expression (log-CPM) by subgroup")
      p_v <- ggplot(plot_df, aes(x = grp_label, y = expr, fill = grp_label)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.12, size = 0.4, alpha = 0.6) +
        facet_wrap(~TranscriptID, scales = "free_y", ncol = 2) +
        labs(title = title_main, x = "Group", y = "log-CPM") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
        scale_fill_viridis_d(option = "D")
      vfile <- file.path(outdir_gene, paste0(g, "_violin_all_transcripts.png"))
      ggsave(vfile, p_v, width = 12, height = max(4, 1.5 * ceiling(length(trans_ids) / 2)))
      message("  Saved per-gene violin: ", vfile)
    } else {
      message("  No valid samples for per-gene violin for ", g)
    }

    # Histogram of per-transcript mean expression (mean log-CPM)
    mean_by_tx <- rowMeans(logcpm_sub, na.rm = TRUE)
    if (length(mean_by_tx) > 0) {
      df_hist_tx <- data.frame(mean_expr = mean_by_tx, tx = names(mean_by_tx), stringsAsFactors = FALSE)
      h <- hist(df_hist_tx$mean_expr, breaks = 30, plot = FALSE)
      mids <- h$mids; counts <- h$counts; med <- median(df_hist_tx$mean_expr, na.rm = TRUE); br <- h$breaks
      med_idx <- findInterval(med, br, rightmost.closed = TRUE)
      df_h <- data.frame(mid = mids, count = counts, is_med = seq_along(mids) == med_idx)
      p_hist_tx <- ggplot(df_h, aes(x = mid, y = count)) +
        geom_col(width = diff(br)[1], fill = viridis(1), color = ifelse(df_h$is_med, "red", "black"), linewidth = ifelse(df_h$is_med, 1, 0.2)) +
        labs(title = paste0(g, " per-transcript mean expression (log-CPM)"),
             x = "Mean log-CPM (per transcript)", y = "Count") +
        theme_minimal() +
        annotate("text", x = med, y = max(df_h$count, na.rm = TRUE) * 1.05,
                 label = paste0("Median=", round(med, 2)), vjust = 0, size = 3)
      hfile <- file.path(outdir_gene, paste0(g, "_histogram_mean_by_transcript.png"))
      ggsave(hfile, p_hist_tx, width = 8, height = 5)
      message("  Saved per-gene histogram: ", hfile)
    } else {
      message("  No transcripts to build per-gene histogram for ", g)
    }
  } # end per-gene plots

  # ----------------- Fit contrasts on voom-subgroup object -------------------
  fG_sub <- lmFit(vGsub, design_sub)

  # 1) Tumor_LOI vs Tumor_noLOI
  if (all(c("Tumor_LOI", "Tumor_noLOI") %in% grp_sub)) {
    cm1 <- makeContrasts(Tumor_LOI_vs_Tumor_noLOI = Tumor_LOI - Tumor_noLOI, levels = design_sub)
    f1 <- contrasts.fit(fG_sub, cm1); f1 <- eBayes(f1)
    tt1 <- topTable(f1, coef = "Tumor_LOI_vs_Tumor_noLOI", number = Inf, adjust.method = "BH", sort.by = "none")
    is_g <- rownames(tt1) %in% trans_ids
    res_g_trans <- tt1[is_g, , drop = FALSE]
    if (nrow(res_g_trans) > 0) {
      results_transcripts[["Tumor_LOI_vs_Tumor_noLOI"]][[g]] <- data.frame(Gene = g, TranscriptID = rownames(res_g_trans), res_g_trans, stringsAsFactors = FALSE)
    }
  }

  # 2) Tumor_LOI vs Control_noLOI
  if (all(c("Tumor_LOI", "Control_noLOI") %in% grp_sub)) {
    cm2 <- makeContrasts(Tumor_LOI_vs_Control_noLOI = Tumor_LOI - Control_noLOI, levels = design_sub)
    f2 <- contrasts.fit(fG_sub, cm2); f2 <- eBayes(f2)
    tt2 <- topTable(f2, coef = "Tumor_LOI_vs_Control_noLOI", number = Inf, adjust.method = "BH", sort.by = "none")
    is_g <- rownames(tt2) %in% trans_ids
    res_g_trans2 <- tt2[is_g, , drop = FALSE]
    if (nrow(res_g_trans2) > 0) {
      results_transcripts[["Tumor_LOI_vs_Control_noLOI"]][[g]] <- data.frame(Gene = g, TranscriptID = rownames(res_g_trans2), res_g_trans2, stringsAsFactors = FALSE)
    }
  }

  # 3) Tumor vs Control (group-level)
  yG$samples$group2 <- factor(ifelse(yG$samples$subgroup %in% c("Tumor_LOI", "Tumor_noLOI"), "Tumor", "Control"))
  design2 <- model.matrix(~0 + group2, data = yG$samples)
  colnames(design2) <- levels(yG$samples$group2)

  voom_group_plot <- file.path(out_dir, paste0(g, "_transcript_voom_group.png"))
  png(voom_group_plot, width = 1200, height = 900, res = 150)
  vG2 <- tryCatch({
    voom(yG2, design2, plot = TRUE)
  }, error = function(e) {
    dev.off(); message("  voom(group) error for ", g, ": ", e$message); NULL
  })
  dev.off()
  if (!is.null(vG2)) {
    fG2 <- lmFit(vG2, design2)
    if (all(c("Tumor", "Control") %in% colnames(design2))) {
      cm3 <- makeContrasts(Tumor_vs_Control = Tumor - Control, levels = design2)
      f3 <- contrasts.fit(fG2, cm3); f3 <- eBayes(f3)
      tt3 <- topTable(f3, coef = "Tumor_vs_Control", number = Inf, adjust.method = "BH", sort.by = "none")
      is_g <- rownames(tt3) %in% trans_ids
      res_g_trans3 <- tt3[is_g, , drop = FALSE]
      if (nrow(res_g_trans3) > 0) {
        results_transcripts[["Tumor_vs_Control"]][[g]] <- data.frame(Gene = g, TranscriptID = rownames(res_g_trans3), res_g_trans3, stringsAsFactors = FALSE)
      }
    }
  }
} # end per-gene loop

# ------------------------- STEP 9: Save DGE workbook -------------------------
wb <- createWorkbook()
for (cn in contrast_names) {
  lst <- results_transcripts[[cn]]
  sheet_name <- gsub("_", " ", cn)
  if (length(lst) == 0) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data.frame(Note = paste("No results for", sheet_name)))
  } else {
    df_all <- bind_rows(lst)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df_all)
    message("Sheet ", sheet_name, ": ", nrow(df_all), " rows")
  }
}
fn <- file.path(out_dir, "Transcript_DGE_CandidateGenes.xlsx")
saveWorkbook(wb, fn, overwrite = TRUE)
message("Saved DGE workbook at: ", fn)

# --------------------- STEP 9b: Formatted Excel per contrast ----------------
for (cn in contrast_names) {
  lst <- results_transcripts[[cn]]
  if (length(lst) == 0) { message("No results to format for contrast: ", cn); next }
  df_all <- bind_rows(lst)
  if (nrow(df_all) == 0) next

  num_cols <- c("logFC", "AveExpr", "B", "t")
  for (nc in num_cols) if (nc %in% colnames(df_all)) df_all[[nc]] <- round(df_all[[nc]], 2)

  pcol <- NULL
  if ("P.Value" %in% colnames(df_all)) pcol <- "P.Value"
  if (is.null(pcol) && "PValue" %in% colnames(df_all)) pcol <- "PValue"
  if (is.null(pcol) && "p.value" %in% colnames(df_all)) pcol <- "p.value"
  if (!is.null(pcol)) df_all[[pcol]] <- sapply(df_all[[pcol]], format_pval_sup)

  wb2 <- createWorkbook()
  addWorksheet(wb2, gsub("_", " ", cn))
  writeData(wb2, gsub("_", " ", cn), df_all)
  fn2 <- file.path(out_dir, paste0("Transcript_DGE_CandidateGenes_formatted_", cn, ".xlsx"))
  saveWorkbook(wb2, fn2, overwrite = TRUE)
  message("Saved formatted contrast workbook: ", fn2)
}

# ------------------------- STEP 10: Save summary stats -----------------------
if (length(summary_stats_trans) > 0) {
  df_summary_all <- bind_rows(summary_stats_trans)
  wb_sum <- createWorkbook()
  addWorksheet(wb_sum, "Transcript_SummaryStats")
  writeData(wb_sum, "Transcript_SummaryStats", df_summary_all)
  fn_sum <- file.path(out_dir, "Transcript_SummaryStats_CandidateGenes.xlsx")
  saveWorkbook(wb_sum, fn_sum, overwrite = TRUE)
  message("Saved summary stats workbook: ", fn_sum)
}

# ------------------------- STEP 11: Dodge bar plots --------------------------
for (cn in contrast_names) {
  lst <- results_transcripts[[cn]]
  if (length(lst) == 0) { message("No results for dodge plot: ", cn); next }
  df_res <- bind_rows(lst) %>%
    dplyr::select(Gene, TranscriptID, logFC) %>%
    mutate(xid = paste0(Gene, "___", TranscriptID)) %>%
    arrange(Gene, TranscriptID)

  xid_levels <- unique(df_res$xid)
  df_res$xid <- factor(df_res$xid, levels = xid_levels)
  xid_to_label <- setNames(sub(".*___", "", xid_levels), xid_levels)

  total_trans <- length(xid_levels)
  per_trans_phys_in <- 0.20
  base_width <- 14
  plot_width <- max(base_width, per_trans_phys_in * total_trans)

  p_dodge <- ggplot(df_res, aes(x = xid, y = logFC, fill = Gene)) +
    geom_col(width = 0.9, color = NA) +
    labs(title = paste0("Transcript logFC for contrast: ", clean_label(cn)),
         x = "Transcript ID", y = "logFC", fill = "Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "top") +
    scale_fill_viridis_d() +
    scale_x_discrete(labels = xid_to_label, expand = expansion(add = c(0, 0)))

  fn_dodge <- file.path(out_dir, paste0("dodge_plot_", cn, ".png"))
  ggsave(fn_dodge, p_dodge, width = plot_width, height = 6, dpi = 150)
  message("Saved dodge bar plot for contrast: ", cn)
}

message("Transcript-level pipeline completed successfully")
