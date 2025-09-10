###############################################################################
# DGE Pipeline — per-gene serial processing + Tumor vs Control + plotting
# #TALAL AMIN: BIOBIX (UGENT)
# Usage: edit the path variables below, then run.
###############################################################################

# ---- User settings (edit these before running) ----
RUN_DUPLICATE_REMOVAL <- TRUE   # remove within-patient duplicates (keeps sample with largest global libsize)
NUM_WORKERS           <- 4      # placeholder (serial mode used)
RUN_DODGE_PLOTS       <- TRUE   # generate bar (dodge) plots
GLOBAL_MIN_SAMPLES    <- 14     # minimum number of samples in which gene CPM must exceed threshold (global)
CPM_THRESHOLD_GLOBAL  <- 1      # CPM threshold for global filter
PCT_REQUIRED          <- 0.15   # fraction of samples in each group required to pass CPM filter (Optional /Additional)
GI_genes              <- c("HM13", "MEST", "IGF2", "ZNF331")  # example genes of interest

# ---- Libraries ----
suppressPackageStartupMessages({
  library(biomaRt)
  library(tximport)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(viridis)
  library(openxlsx)
  library(tools)
  library(ggtext)
  library(dplyr)
})

# ---- Paths (REPLACE these placeholders with your actual directories) ----
# Provide directories that contain kallisto "abundance.h5" outputs for tumor and control samples
tumor_dir   <- "path/to/tumor/"     # <-- set to your tumor kallisto output root
control_dir <- "path/to/control/"   # <-- set to your control kallisto output root
design_dir  <- "path/to/design/"    # <-- set to directory containing per-gene design CSVs (Gene.csv)
out_dir     <- "path/to/output/"    # <-- set to desired output directory

# Prepare output directories
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Gather files ----
tumor_files   <- list.files(tumor_dir,   pattern = "abundance.h5$", full.names = TRUE, recursive = TRUE)
control_files <- list.files(control_dir, pattern = "abundance.h5$", full.names = TRUE, recursive = TRUE)
all_files     <- c(tumor_files, control_files)
# Name samples by the directory name that contains abundance.h5 (commonly sample folder name)
names(all_files) <- basename(dirname(all_files))

# ---- tx2gene mapping via biomaRt ----
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                 mart = ensembl)
colnames(tx2gene) <- c("tx_id", "gene_id", "external_gene_name")

# ---- Cache filenames and per-gene directory ----
global_txi_rds   <- file.path(output_dir, "global_txi.rds")
global_libs_rds  <- file.path(output_dir, "global_libsizes.rds")
kept_genes_rds   <- file.path(output_dir, "kept_genes_global_CPMfiltered.rds")
per_gene_rds_dir <- file.path(output_dir, "per_gene_rds")
if (!dir.exists(per_gene_rds_dir)) dir.create(per_gene_rds_dir, recursive = TRUE)

# -----------------------------------------------------------------
# 1) Global tximport + CPM-based global filter (cached)
# -----------------------------------------------------------------
if (file.exists(global_txi_rds) && file.exists(global_libs_rds) && file.exists(kept_genes_rds)) {
  message("Loading cached global objects...")
  global_txi <- readRDS(global_txi_rds)
  global_libsizes <- readRDS(global_libs_rds)
  kept_genes_global <- readRDS(kept_genes_rds)
  message("Loaded cached global objects.")
} else {
  message("Running global tximport (no cache found).")
  files_for_global <- all_files
  global_txi <- tximport(files = files_for_global,
                         type = "kallisto",
                         tx2gene = tx2gene[, c("tx_id", "gene_id")],
                         countsFromAbundance = "scaledTPM",
                         ignoreTxVersion = TRUE)
  saveRDS(global_txi, file = global_txi_rds)

  global_y <- DGEList(counts = global_txi$counts)
  global_libsizes <- colSums(global_txi$counts, na.rm = TRUE)
  saveRDS(global_libsizes, file = global_libs_rds)

  global_cpm <- cpm(global_y)
  global_keep_mask <- rowSums(global_cpm > CPM_THRESHOLD_GLOBAL) >= GLOBAL_MIN_SAMPLES
  kept_genes_global <- rownames(global_y)[global_keep_mask]
  saveRDS(kept_genes_global, file = kept_genes_rds)
  message("Global filtering complete; genes kept:", length(kept_genes_global))
}

# -----------------------------------------------------------------
# Helper functions: annotation, plotting, formatting
# -----------------------------------------------------------------

# Annotate topTable results with gene symbols and attach contrast label
topTable_with_annotations <- function(fit, contrast, tx2gene) {
  tt <- topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
  if (nrow(tt) == 0) return(tt)
  tt$ID <- rownames(tt)
  annotated_table <- merge(tt, unique(tx2gene[, c("gene_id", "external_gene_name")]),
                           by.x = "ID", by.y = "gene_id", all.x = TRUE)
  annotated_table$Gene <- annotated_table$external_gene_name
  annotated_table <- annotated_table[, !(colnames(annotated_table) == "external_gene_name")]
  annotated_table$Contrast <- contrast
  return(annotated_table)
}

# Single-contrast bar plot for genes of interest (underscores replaced by spaces)
create_GI_gene_regulation_plot <- function(contrast_name, combined_df_for_contrast, GI_genes, output_dir) {
  if (!RUN_DODGE_PLOTS) return(NULL)
  df <- combined_df_for_contrast
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  if (!"logFC" %in% colnames(df)) return(NULL)

  df$PlotGene <- ifelse(is.na(df$Gene) | df$Gene == "", df$ID, df$Gene)
  df$PlotGene <- gsub("_", " ", df$PlotGene)
  if (!"QueriedGene" %in% colnames(df)) df$QueriedGene <- df$PlotGene
  df$QueriedGeneDisplay <- gsub("_", " ", df$QueriedGene)
  df$xid <- factor(df$QueriedGeneDisplay, levels = unique(df$QueriedGeneDisplay))

  title_text <- gsub("_", " ", gsub("_vs_", " vs ", contrast_name))

  unique_genes <- unique(df$PlotGene)
  colors <- viridis::viridis(length(unique_genes))
  names(colors) <- unique_genes

  p <- ggplot(df, aes(x = xid, y = logFC, fill = PlotGene)) +
    geom_col(width = 0.7, color = NA) +
    scale_fill_manual(values = colors) +
    labs(title = paste0("Differential Gene Expression : ", title_text),
         x = "", y = "Log2 Fold Change (logFC)", fill = "") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  out_file <- file.path(output_dir, paste0(gsub("\\s+","_", contrast_name), "_combined_dodge_plot.png"))
  ggsave(filename = out_file, plot = p, width = 8, height = 5, dpi = 300)
  return(out_file)
}

# Three-panel combined plot for three main contrasts
create_three_panel_plot <- function(combined_results, GI_genes, output_dir) {
  if (!RUN_DODGE_PLOTS) return(NULL)
  cn_map <- c(
    Tumor_vs_Control = "A: Tumor vs Control",
    Tumor_LOI_vs_Tumor_noLOI = "B: Tumor LOI vs Tumor noLOI",
    Tumor_LOI_vs_Control_noLOI = "C: Tumor LOI vs Control noLOI"
  )

  df_list <- list()
  for (cn in names(cn_map)) {
    df <- combined_results[[cn]]
    if (is.null(df) || nrow(df) == 0) next
    df$Panel <- cn_map[[cn]]
    df$ContrastName <- cn
    df_list[[cn]] <- df
  }
  if (length(df_list) == 0) return(NULL)
  df_all <- bind_rows(df_list)
  if (!"logFC" %in% colnames(df_all)) return(NULL)
  df_all <- df_all[!is.na(df_all$logFC), , drop = FALSE]
  if (nrow(df_all) == 0) return(NULL)

  df_all$PlotGene <- ifelse(is.na(df_all$Gene) | df_all$Gene == "", df_all$ID, df_all$Gene)
  df_all$PlotGene <- gsub("_", " ", df_all$PlotGene)
  if (!"QueriedGene" %in% colnames(df_all)) df_all$QueriedGene <- df_all$PlotGene
  df_all$QueriedGeneDisplay <- gsub("_", " ", df_all$QueriedGene)

  # Keep only requested GI genes
  df_all <- df_all[df_all$QueriedGene %in% GI_genes, , drop = FALSE]
  if (nrow(df_all) == 0) return(NULL)

  df_all$xid <- factor(df_all$QueriedGeneDisplay, levels = unique(df_all$QueriedGeneDisplay))

  unique_genes <- unique(df_all$PlotGene)
  colors <- viridis::viridis(length(unique_genes))
  names(colors) <- unique_genes

  p <- ggplot(df_all, aes(x = xid, y = logFC, fill = PlotGene)) +
    geom_col(width = 0.7, color = NA) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ Panel, ncol = 1, scales = "free_y", strip.position = "top") +
    labs(title = "Differential expression of imprinted genes among conditions",
         x = "", y = "Log2 Fold Change (logFC)", fill = "") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  out_file <- file.path(output_dir, "DGE_combined_three_panel.png")
  ggsave(filename = out_file, plot = p, width = 10, height = 14, dpi = 300)
  return(out_file)
}

# Formatting p-values and rounding numeric columns for human-friendly Excel
format_pval <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (!is.numeric(x)) return(as.character(x))
  if (x < 0.10) {
    return(formatC(x, format = "e", digits = 2))
  } else {
    return(formatC(round(x, 2), format = "f", digits = 2))
  }
}
round_and_format_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, 2))
  pval_cols <- intersect(c("P.Value", "PValue", "adj.P.Val", "FDR", "P.value"), colnames(df))
  for (pc in pval_cols) {
    df[[pc]] <- vapply(df[[pc]], format_pval, FUN.VALUE = character(1))
  }
  df
}

# -----------------------------------------------------------------
# Per-gene processing (serial) with caching
# -----------------------------------------------------------------
contrast_names <- c("Tumor_LOI_vs_Tumor_noLOI", "Tumor_LOI_vs_Control_noLOI", "Tumor_vs_Control")
combined_results <- setNames(vector("list", length(contrast_names)), contrast_names)
for (cn in contrast_names) combined_results[[cn]] <- data.frame()

process_one_gene <- function(gene_name) {
  message("Processing gene: ", gene_name)
  gene_outdir <- file.path(output_dir, gene_name)
  if (!dir.exists(gene_outdir)) dir.create(gene_outdir, recursive = TRUE)

  # Expect per-gene design CSV: columns at minimum 'samples' and 'group'
  design_file <- file.path(design_dir, paste0(gene_name, ".csv"))
  if (!file.exists(design_file)) {
    message("Design file missing for ", gene_name, " -> skipping")
    return(NULL)
  }
  sample_info <- read.csv(design_file, stringsAsFactors = FALSE)
  design_samples <- file_path_sans_ext(basename(sample_info$samples))
  valid_samples <- intersect(design_samples, names(all_files))
  if (length(valid_samples) == 0) {
    message("No valid samples for ", gene_name)
    return(NULL)
  }

  # Duplicate removal within (patient, group) using global library sizes
  if (RUN_DUPLICATE_REMOVAL) {
    pid <- substr(valid_samples, 1, 12) # default heuristic: first 12 chars as patient id
    df_dup <- data.frame(sample = valid_samples,
                         patient = pid,
                         group = sample_info$group[match(valid_samples, design_samples)],
                         stringsAsFactors = FALSE)
    keys <- paste0(df_dup$patient, "::", df_dup$group)
    dup_keys <- names(table(keys))[table(keys) > 1]
    if (length(dup_keys) > 0) {
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
      if (dropped_n > 0) message("Dropped ", dropped_n, " within-patient replicate(s) for gene ", gene_name)
      valid_samples <- kept
    }
  }

  # Build sample table for this gene
  sampleTable <- data.frame(sample = valid_samples,
                            file_path = all_files[valid_samples],
                            group = sample_info$group[match(valid_samples, design_samples)],
                            stringsAsFactors = FALSE)
  sampleTable <- sampleTable[!is.na(sampleTable$group), ]
  if (nrow(sampleTable) == 0) { message("No samples with valid group for ", gene_name); return(NULL) }

  per_gene_rds <- file.path(per_gene_rds_dir, paste0(gene_name, "_top_tables.rds"))
  top_tables <- NULL

  if (file.exists(per_gene_rds)) {
    message("Loading cached per-gene results for ", gene_name)
    top_tables <- readRDS(per_gene_rds)
  } else {
    # tximport for this gene's samples
    valid_files <- sampleTable$file_path
    names(valid_files) <- sampleTable$sample
    txi <- tximport(files = valid_files, type = "kallisto",
                    tx2gene = tx2gene[, c("tx_id", "gene_id")],
                    countsFromAbundance = "scaledTPM", ignoreTxVersion = TRUE)

    # Restrict to genes that passed global filter
    keep_rows <- intersect(rownames(txi$counts), kept_genes_global)
    if (length(keep_rows) == 0) {
      message("No genes remain after global filter for gene ", gene_name)
      return(NULL)
    }

    y <- DGEList(counts = txi$counts[keep_rows, , drop = FALSE])
    rownames(y$samples) <- colnames(y$counts)
    y$samples$group <- factor(sampleTable$group[match(colnames(y$counts), sampleTable$sample)])
    y <- calcNormFactors(y, method = "TMM")

    # subgroup design (e.g., Tumor_LOI, Tumor_noLOI, Control_LOI, Control_noLOI)
    design <- model.matrix(~ 0 + group, data = y$samples)
    colnames(design) <- levels(y$samples$group)

    # Build contrasts only for subgroup comparisons present in the design
    contrasts_to_build <- list()
    if (all(c("Tumor_LOI", "Tumor_noLOI") %in% colnames(design))) {
      contrasts_to_build[["Tumor_LOI_vs_Tumor_noLOI"]] <- "Tumor_LOI - Tumor_noLOI"
    }
    if (all(c("Tumor_LOI", "Control_noLOI") %in% colnames(design))) {
      contrasts_to_build[["Tumor_LOI_vs_Control_noLOI"]] <- "Tumor_LOI - Control_noLOI"
    }

    # Filtering with group-aware CPM criteria
    keep_expr <- filterByExpr(y, design)
    cpm_values <- cpm(y)
    groups <- levels(y$samples$group)
    keep_cpm_per_group <- rep(TRUE, nrow(y))
    for (i in seq_along(groups)) {
      g <- groups[i]
      idx <- which(y$samples$group == g)
      required_count <- ceiling(PCT_REQUIRED * length(idx))
      if (required_count < 1) next
      mask <- rowSums(cpm_values[, idx, drop = FALSE] > 1) >= required_count
      keep_cpm_per_group <- keep_cpm_per_group & mask
    }
    if (any(ceiling(PCT_REQUIRED * as.integer(table(y$samples$group))) < 1)) {
      keep_final <- keep_expr
      message("Small group -> using filterByExpr only for ", gene_name)
    } else {
      keep_final <- keep_expr & keep_cpm_per_group
    }

    y <- y[keep_final, , keep.lib.sizes = FALSE]

    # Voom + linear modeling for subgroup contrasts
    v <- voom(y, design, plot = FALSE)
    fit <- lmFit(v, design)
    top_tables <- list()

    if (length(contrasts_to_build) > 0) {
      contr_args <- c(contrasts_to_build, list(levels = design))
      cm_sub <- tryCatch(do.call(makeContrasts, contr_args), error = function(e) {
        message("makeContrasts failed for subgroups in ", gene_name, ": ", e$message); NULL
      })
      if (!is.null(cm_sub)) {
        fit2 <- tryCatch(contrasts.fit(fit, cm_sub), error = function(e) {
          message("contrasts.fit error (subgroup) for ", gene_name, ": ", e$message); NULL
        })
        if (!is.null(fit2)) {
          fit2 <- eBayes(fit2)
          for (cn in names(contrasts_to_build)) {
            tt_cn <- tryCatch(topTable_with_annotations(fit2, cn, tx2gene), error = function(e) NULL)
            if (is.null(tt_cn)) tt_cn <- data.frame()
            top_tables[[cn]] <- tt_cn
          }
        } else {
          for (cn in names(contrasts_to_build)) top_tables[[cn]] <- data.frame()
        }
      } else {
        for (cn in names(contrasts_to_build)) top_tables[[cn]] <- data.frame()
      }
    } else {
      # Ensure keys exist even if no subgroup contrasts were built
      top_tables[["Tumor_LOI_vs_Tumor_noLOI"]] <- data.frame()
      top_tables[["Tumor_LOI_vs_Control_noLOI"]] <- data.frame()
    }

    # ----------------------------------------------------------------
    # Tumor vs Control comparison (origin-based grouping)
    # ----------------------------------------------------------------
    samples_all <- colnames(y$counts)
    tumor_samples_ids   <- unique(basename(dirname(tumor_files)))
    control_samples_ids <- unique(basename(dirname(control_files)))
    group2_vec <- ifelse(samples_all %in% tumor_samples_ids, "Tumor",
                         ifelse(samples_all %in% control_samples_ids, "Control", NA))
    if (any(is.na(group2_vec))) {
      warning("Dropping samples not in tumor/control dirs for gene ", gene_name)
    }
    keep_samps <- samples_all[!is.na(group2_vec)]
    group2_vec  <- group2_vec[!is.na(group2_vec)]

    if (length(keep_samps) > 1) {
      y_group <- y[, keep_samps, keep.lib.sizes = FALSE]
      rownames(y_group$samples) <- colnames(y_group)
      y_group$samples$group2 <- factor(group2_vec, levels = c("Control", "Tumor"))

      if (length(unique(y_group$samples$group2)) == 2) {
        design2 <- model.matrix(~ 0 + group2, data = y_group$samples)
        colnames(design2) <- levels(y_group$samples$group2)
        vG <- tryCatch(voom(y_group, design2, plot = FALSE), error = function(e) {
          message("voom(group) error for ", gene_name, ": ", e$message); NULL
        })
        if (!is.null(vG)) {
          fG <- lmFit(vG, design2)
          if (all(c("Tumor", "Control") %in% colnames(design2))) {
            cm3 <- tryCatch(makeContrasts(Tumor_vs_Control = Tumor - Control, levels = design2),
                            error = function(e) { message("makeContrasts group failed for ", gene_name, ": ", e$message); NULL })
            if (!is.null(cm3)) {
              f3 <- contrasts.fit(fG, cm3)
              f3 <- eBayes(f3)
              tt3 <- topTable_with_annotations(f3, "Tumor_vs_Control", tx2gene)
              top_tables[["Tumor_vs_Control"]] <- tt3
            } else top_tables[["Tumor_vs_Control"]] <- data.frame()
          } else top_tables[["Tumor_vs_Control"]] <- data.frame()
        } else top_tables[["Tumor_vs_Control"]] <- data.frame()
      } else top_tables[["Tumor_vs_Control"]] <- data.frame()
    } else top_tables[["Tumor_vs_Control"]] <- data.frame()

    # Ensure that each contrast key exists
    for (cn in contrast_names) if (is.null(top_tables[[cn]])) top_tables[[cn]] <- data.frame()

    # Cache per-gene results
    saveRDS(top_tables, file = per_gene_rds)
    message("Saved per-gene results for ", gene_name)
  } # end per-gene caching else

  # Map gene identifiers & generate placeholder rows if no row is present for requested gene
  gene_ids <- unique(tx2gene$gene_id[tx2gene$external_gene_name == gene_name])
  if (length(gene_ids) == 0) gene_ids <- character(0)

  per_gene_rows <- list()
  for (cn in contrast_names) {
    tt <- top_tables[[cn]]
    if (is.null(tt) || nrow(tt) == 0) { per_gene_rows[[cn]] <- NULL; next }
    row_match <- NULL
    if (length(gene_ids) > 0) row_match <- tt[tt$ID %in% gene_ids, , drop = FALSE]
    if (nrow(row_match) == 0 && "Gene" %in% colnames(tt)) row_match <- tt[tt$Gene == gene_name, , drop = FALSE]
    if (nrow(row_match) == 0) {
      placeholder <- data.frame(ID = NA_character_, Gene = gene_name, Contrast = cn, QueriedGene = gene_name, stringsAsFactors = FALSE)
      per_gene_rows[[cn]] <- placeholder
    } else {
      row_match$QueriedGene <- gene_name
      per_gene_rows[[cn]] <- row_match
    }
  }

  return(list(per_gene_rows = per_gene_rows, gene = gene_name))
}

# Run per-gene serially (lapply avoids parallel serialization issues)
res_list <- lapply(GI_genes, process_one_gene)
names(res_list) <- GI_genes

# Combine per-contrast results into tables
for (res in res_list) {
  if (is.null(res)) next
  gene_name <- res$gene
  per_gene_rows <- res$per_gene_rows
  for (cn in contrast_names) {
    row_df <- per_gene_rows[[cn]]
    if (is.null(row_df)) next
    if (!"QueriedGene" %in% colnames(row_df)) row_df$QueriedGene <- gene_name
    combined_results[[cn]] <- bind_rows(combined_results[[cn]], row_df)
  }
}

# -----------------------------------------------------------------
# Write final combined workbook: one sheet per contrast
# -----------------------------------------------------------------
final_wb <- createWorkbook()
for (cn in contrast_names) {
  df <- combined_results[[cn]]
  if (is.null(df) || nrow(df) == 0) {
    addWorksheet(final_wb, sheetName = cn)
    writeData(final_wb, sheet = cn, x = data.frame(note = paste("No results for contrast", cn)))
    next
  }
  if ("QueriedGene" %in% colnames(df)) df <- df %>% select(QueriedGene, everything())
  addWorksheet(final_wb, sheetName = cn)
  writeData(final_wb, sheet = cn, x = df)
}
final_xlsx <- file.path(output_dir, "DGE_combined_results_per_contrast.xlsx")
saveWorkbook(final_wb, file = final_xlsx, overwrite = TRUE)
message("Wrote final combined results workbook: ", final_xlsx)

# Write rounded workbook for human-friendly display
rounded_wb <- createWorkbook()
for (cn in contrast_names) {
  df <- combined_results[[cn]]
  if (is.null(df) || nrow(df) == 0) {
    addWorksheet(rounded_wb, sheetName = cn)
    writeData(rounded_wb, sheet = cn, x = data.frame(note = paste("No results for contrast", cn)))
    next
  }
  df2 <- round_and_format_df(df)
  addWorksheet(rounded_wb, sheetName = cn)
  writeData(rounded_wb, sheet = cn, x = df2)
}
rounded_xlsx <- file.path(output_dir, "DGE_combined_results_per_contrast_rounded.xlsx")
saveWorkbook(rounded_wb, file = rounded_xlsx, overwrite = TRUE)
message("Wrote rounded combined workbook: ", rounded_xlsx)

# Save combined results RDS
saveRDS(combined_results, file = file.path(output_dir, "combined_results_per_contrast.rds"))
message("Saved combined results RDS.")

# Generate plots (one per contrast) and 3-panel plot
for (cn in contrast_names) {
  df <- combined_results[[cn]]
  if (is.null(df) || nrow(df) == 0) next
  if (!"logFC" %in% colnames(df)) next
  create_GI_gene_regulation_plot(cn, df, GI_genes, output_dir)
}
three_panel_png <- create_three_panel_plot(combined_results, GI_genes, output_dir)
if (!is.null(three_panel_png)) message("Wrote 3-panel combined dodge plot: ", three_panel_png)

# Final summary messages
message("All done. Outputs written to: ", output_dir)
message("- Final Excel:", final_xlsx)
message("- Rounded Excel:", rounded_xlsx)
message("- Combined RDS:", file.path(output_dir, "combined_results_per_contrast.rds"))
message("- Per-gene RDS directory:", per_gene_rds_dir)
message("- Global kept genes RDS:", kept_genes_rds)

###############################################################################
# End of script
###############################################################################