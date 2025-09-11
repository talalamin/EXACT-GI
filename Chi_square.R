# R script for LOI (Loss of Imprinting) sample-level calls
# - Consolidated duplicated steps and added robust checks
# - Parameterized input/output paths and optional files

# Libraries ------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(stringr)

# ------------------------------ Parameters ----------------------------------
# EDIT THESE before running -------------------------------------------------
input_excel <- "path/to/TCGA_DI_SNPs.xlsx"   # Excel file with sheets: Tumor_LOI, Tumor_noLOI, Control_LOI, Control_noLOI
output_dir  <- "path/to/output_dir"          # Directory where result files will be written
# Optional files (set to NULL if not used)
imprinted_snps_file <- NULL                   # path/to/imprinted_SNPs.txt (one rsID per line) or NULL
low_cov_file         <- NULL                   # path/to/Gene_Low_Coverage_Samples.csv or NULL

# Analysis settings ----------------------------------------------------------
# Range of LAF thresholds to test (adjust if desired)
threshold_range <- seq(0.01, 0.5, by = 0.001)
# If you already have a SNP list to subset to, provide it here (character vector) or set to NULL to keep all SNPs
snps_of_interest <- NULL

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------- Helpers ------------------------------------
# Safe read of a sheet (returns NULL if missing)
read_sheet_safe <- function(path, sheet) {
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
  if (!(sheet %in% sheets)) return(NULL)
  readxl::read_excel(path, sheet = sheet)
}

# Convert specified columns to numeric if they exist
ensure_numeric_cols <- function(df, cols = c("Ref_Count", "Alt_Count", "LAF")) {
  if (is.null(df)) return(NULL)
  cols_present <- intersect(cols, colnames(df))
  if (length(cols_present) > 0) {
    df <- df %>% mutate(across(all_of(cols_present), ~ as.numeric(.x)))
  }
  df
}

# Compute preliminary union call: for each Gene x Sample, if any SNP has min(Ref,Alt) >= min_biallelic then call LOI
compute_union_call <- function(df, min_biallelic = 2) {
  df %>%
    group_by(Gene, Sample_Name) %>%
    summarise(
      Union_Call = ifelse(any(pmin(Ref_Count, Alt_Count) >= min_biallelic), "LOI", "noLOI"),
      .groups = "drop"
    ) %>%
    left_join(df, by = c("Gene", "Sample_Name"))
}

# Compute sample-level union calls at a given LAF threshold (per gene)
compute_union_call_threshold <- function(df, threshold) {
  df %>%
    group_by(Sample_Name, Group) %>%
    summarise(
      Union_Call = ifelse(any(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2), "LOI", "noLOI"),
      SNPs_Passing = paste(unique(SNP_ID[(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2)]), collapse = ";"),
      Count_Passed = sum(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2),
      Total_SNPs = n(),
      .groups = "drop"
    )
}

# Modified final-call function that optionally checks for imprinted SNP presence and low-coverage samples
compute_union_call_threshold_modified <- function(df, threshold, imprinted_snps = NULL, low_cov_samples = character(0)) {
  union_calls <- df %>%
    group_by(Sample_Name, Group) %>%
    summarise(
      LOI_call = ifelse(any(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2), "LOI", "noLOI"),
      SNPs_Passing = paste(unique(SNP_ID[(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2)]), collapse = ";"),
      Count_Passed = sum(LAF >= threshold & pmin(Ref_Count, Alt_Count) >= 2),
      Total_SNPs = n(),
      .groups = "drop"
    )

  extra_info <- df %>%
    group_by(Sample_Name) %>%
    summarise(
      Imprinted_SNP_Found = any(SNP_ID %in% (imprinted_snps %||% character(0))),
      .groups = "drop"
    )

  union_calls <- left_join(union_calls, extra_info, by = "Sample_Name")

  union_calls %>%
    mutate(
      Final_Call = case_when(
        LOI_call == "LOI" ~ "LOI",
        LOI_call == "noLOI" & Imprinted_SNP_Found ~ "noLOI",
        Sample_Name %in% low_cov_samples ~ "Excluded",
        TRUE ~ "Excluded"
      )
    )
}

# Helper: safe infix for NULL default
`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------ Step 1 --------------------------------------
# Read the four expected sheets (Tumor_LOI, Tumor_noLOI, Control_LOI, Control_noLOI)
message("Reading input sheets...")
Tumor_LOI     <- read_sheet_safe(input_excel, "Tumor_LOI")
Tumor_noLOI   <- read_sheet_safe(input_excel, "Tumor_noLOI")
Control_LOI   <- read_sheet_safe(input_excel, "Control_LOI")
Control_noLOI <- read_sheet_safe(input_excel, "Control_noLOI")

# Verify that at least one tumor and one control sheet are present
if (all(sapply(list(Tumor_LOI, Tumor_noLOI, Control_LOI, Control_noLOI), is.null))) {
  stop("No sheets found in input file. Verify 'input_excel' and sheet names.")
}

# Convert numeric-like columns to numeric for all loaded tables
Tumor_LOI     <- ensure_numeric_cols(Tumor_LOI)
Tumor_noLOI   <- ensure_numeric_cols(Tumor_noLOI)
Control_LOI   <- ensure_numeric_cols(Control_LOI)
Control_noLOI <- ensure_numeric_cols(Control_noLOI)

# Combine tumor and control sets separately; keep only first 12 columns if present
select_first_n <- function(df, n = 12) {
  if (is.null(df)) return(NULL)
  df %>% select(1:min(n, ncol(.)))
}

tumor_data    <- bind_rows(select_first_n(Tumor_LOI), select_first_n(Tumor_noLOI))
control_data  <- bind_rows(select_first_n(Control_LOI), select_first_n(Control_noLOI))

# Basic column checks and required column names
required_cols <- c("Chromosome", "Position", "Ref_Allele", "Alt_Allele", "SNP_ID", "Gene",
                   "Ref_Count", "Alt_Count", "Sample_Name", "expression", "Goovaerts_et_al", "LAF")

# If column names differ, attempt best-effort mapping (common alternatives)
col_map <- function(df) {
  if (is.null(df)) return(NULL)
  names(df) <- make.names(names(df))
  df
}

tumor_data   <- col_map(tumor_data)
control_data <- col_map(control_data)

# Add Group column and combine tumor+control for downstream steps
tumor_data   <- tumor_data %>% mutate(Group = "Tumor")
control_data <- control_data %>% mutate(Group = "Control")
all_data     <- bind_rows(tumor_data, control_data)

# Optional SNP subsetting
if (!is.null(snps_of_interest)) {
  all_data <- all_data %>% filter(SNP_ID %in% snps_of_interest)
}

# ------------------------------ Step 2 --------------------------------------
# For each gene, test a range of LAF thresholds and choose an optimal one using chi-square on Tumor vs Control
message("Computing threshold summaries per gene...")

# Read optional auxiliary files if paths provided
imprinted_snps <- if (!is.null(imprinted_snps_file) && file.exists(imprinted_snps_file)) readLines(imprinted_snps_file) else NULL
low_cov_df      <- if (!is.null(low_cov_file) && file.exists(low_cov_file)) read.csv(low_cov_file, stringsAsFactors = FALSE) else NULL

# Prepare outputs
gene_list <- unique(na.omit(all_data$Gene))
if (length(gene_list) == 0) stop("No genes found in the combined input data. Check column names and content.")

gene_threshold_summary_list <- list()
final_cutoff_summary_list   <- list()
final_union_list            <- list()

for (g in gene_list) {
  message("Processing gene: ", g)
  gene_data <- all_data %>% filter(Gene == g)
  if (nrow(gene_data) == 0) next

  gene_summary <- vector("list", length(threshold_range))

  idx <- 1
  for (th in threshold_range) {
    union_calls <- compute_union_call_threshold(gene_data, th)

    # Simpler robust counts
    tLOI   <- sum(union_calls$Union_Call == "LOI" & union_calls$Group == "Tumor", na.rm = TRUE)
    tnoLOI <- sum(union_calls$Union_Call == "noLOI" & union_calls$Group == "Tumor", na.rm = TRUE)
    cLOI   <- sum(union_calls$Union_Call == "LOI" & union_calls$Group == "Control", na.rm = TRUE)
    cnoLOI <- sum(union_calls$Union_Call == "noLOI" & union_calls$Group == "Control", na.rm = TRUE)

    # chi-square test if there are observations in both groups
    p_val <- NA
    if ((tLOI + tnoLOI) > 0 && (cLOI + cnoLOI) > 0) {
      chi_table <- matrix(c(tLOI, tnoLOI, cLOI, cnoLOI), nrow = 2, byrow = TRUE)
      p_val <- tryCatch(chisq.test(chi_table)$p.value, error = function(e) NA)
    }

    gene_summary[[idx]] <- data.frame(
      Gene = g,
      LAF_Threshold = th,
      Tumor_LOI = tLOI,
      Tumor_noLOI = tnoLOI,
      Control_LOI = cLOI,
      Control_noLOI = cnoLOI,
      Chi_Square_P = p_val,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }

  gene_summary_df <- do.call(rbind, gene_summary)
  gene_threshold_summary_list[[g]] <- gene_summary_df

  # Select best threshold: prefer smallest (most significant) chi-square p-value, otherwise fallback to max separation
  if (any(!is.na(gene_summary_df$Chi_Square_P))) {
    best_row <- gene_summary_df %>% filter(!is.na(Chi_Square_P)) %>% arrange(Chi_Square_P) %>% slice(1)
  } else {
    # fallback: pick threshold maximizing absolute difference between tumor LOI and noLOI counts
    gene_summary_df <- gene_summary_df %>% mutate(AbsDiff = abs(Tumor_LOI - Tumor_noLOI) + abs(Control_LOI - Control_noLOI))
    best_row <- gene_summary_df %>% arrange(desc(AbsDiff)) %>% slice(1)
  }

  best_threshold <- best_row$LAF_Threshold
  best_pvalue    <- best_row$Chi_Square_P

  final_cutoff_summary_list[[g]] <- data.frame(
    Gene = g,
    Best_LAF_Threshold = best_threshold,
    Chi_Square_P = best_pvalue,
    Tumor_LOI = best_row$Tumor_LOI,
    Tumor_noLOI = best_row$Tumor_noLOI,
    Control_LOI = best_row$Control_LOI,
    Control_noLOI = best_row$Control_noLOI,
    stringsAsFactors = FALSE
  )

  # Low coverage sample list for this gene (if provided in a compatible format)
  low_cov_samples <- character(0)
  if (!is.null(low_cov_df) && "Gene" %in% colnames(low_cov_df) && "sample_names" %in% colnames(low_cov_df)) {
    if (g %in% low_cov_df$Gene) {
      low_cov_row <- low_cov_df %>% filter(Gene == g)
      low_cov_samples <- unlist(strsplit(low_cov_row$sample_names[1], ";")) %>% str_trim()
    }
  }

  final_union_calls <- compute_union_call_threshold_modified(gene_data, best_threshold, imprinted_snps, low_cov_samples) %>%
    mutate(Gene = g, Best_LAF_Threshold = best_threshold)

  final_union_list[[g]] <- final_union_calls
}

# Combine lists into data frames
gene_threshold_summary_df <- do.call(rbind, gene_threshold_summary_list)
final_cutoff_summary_df   <- do.call(rbind, final_cutoff_summary_list)
final_union_calls_df      <- do.call(rbind, final_union_list)

# Save outputs
write_xlsx(gene_threshold_summary_df, file.path(output_dir, "DI_Gene_Threshold_Summary_All.xlsx"))
write_xlsx(final_cutoff_summary_df,   file.path(output_dir, "DI_Final_Gene_Cutoff_Summary.xlsx"))
write_xlsx(final_union_calls_df,      file.path(output_dir, "DI_Gene_Sample_LOI_Final.xlsx"))

message("Step 2 complete: Threshold summaries and final union calls saved to: ", output_dir)

# ------------------------------ Step 3 --------------------------------------
# Sample retention and noLOI filtering summary (optional)
message("Computing sample retention and noLOI filtering summary...")

sample_retention_summary <- final_union_calls_df %>%
  group_by(Gene) %>%
  summarise(
    Total_Samples = n(),
    Retained = sum(Final_Call %in% c("LOI", "noLOI")),
    Excluded = sum(Final_Call == "Excluded"),
    Perc_Retained = (Retained / Total_Samples) * 100,
    Perc_Excluded = (Excluded / Total_Samples) * 100,
    .groups = "drop"
  )

# noLOI extra filtering: check per-sample coverage thresholds and LAF relative to optimal cutoff
cov_thresholds <- c(4, 10)
noLOI_filter_summary <- list()

for (g in unique(final_union_calls_df$Gene)) {
  opt_thresh <- final_cutoff_summary_df %>% filter(Gene == g) %>% pull(Best_LAF_Threshold)
  gene_noLOI <- final_union_calls_df %>% filter(Gene == g, Final_Call == "noLOI")
  orig_noLOI_count <- nrow(gene_noLOI)

  retained_counts <- sapply(cov_thresholds, function(cov_thresh) {
    retained <- 0
    for (samp in gene_noLOI$Sample_Name) {
      sample_data <- all_data %>% filter(Gene == g, Sample_Name == samp)
      if (nrow(sample_data) == 0) next
      sample_data <- sample_data %>% mutate(Total_Count = pmax(0, Ref_Count) + pmax(0, Alt_Count))
      pass_cov <- all(sample_data$Total_Count >= cov_thresh)
      pass_laf <- all(sample_data$LAF <= (opt_thresh / 2))
      if (pass_cov && pass_laf) retained <- retained + 1
    }
    retained
  })

  removed <- orig_noLOI_count - retained_counts
  perc_removed <- if (orig_noLOI_count > 0) (removed / orig_noLOI_count) * 100 else NA

  noLOI_filter_summary[[g]] <- data.frame(
    Gene = g,
    Original_noLOI = orig_noLOI_count,
    Retained_noLOI_cov4 = retained_counts[1],
    Removed_noLOI_cov4 = removed[1],
    Perc_Removed_noLOI_cov4 = perc_removed[1],
    Retained_noLOI_cov10 = retained_counts[2],
    Removed_noLOI_cov10 = removed[2],
    Perc_Removed_noLOI_cov10 = perc_removed[2],
    stringsAsFactors = FALSE
  )
}

noLOI_filter_summary_df <- do.call(rbind, noLOI_filter_summary)
final_sample_summary <- left_join(sample_retention_summary, noLOI_filter_summary_df, by = "Gene")

write_xlsx(final_sample_summary, file.path(output_dir, "Sample_Retention_Summary.xlsx"))
message("Step 3 complete: Sample retention and noLOI filtering summary saved.")

# ----------------------------------------------------------------------------
# Script complete. 
# Edit parameters at the top (input_excel, output_dir, optional files)
# ----------------------------------------------------------------------------