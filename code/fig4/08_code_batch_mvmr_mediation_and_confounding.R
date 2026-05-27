rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(MVMR)
})

required_pkgs <- c("data.table", "dplyr", "TwoSampleMR", "MVMR")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

pval_threshold <- 5e-6
mv_multiple_pval_threshold <- 5e-6
clump_r2 <- 0.001
clump_kb <- 10000
harmonise_strictness <- 2

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])

  if (length(script_path) > 0) {
    return(normalizePath(script_path[[1]], winslash = "/", mustWork = TRUE))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) {
      return(frame$ofile)
    }
    ""
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

prepare_sumstats_no_frq <- function(path, phenotype, id) {
  dt <- fread(file = path, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
  dt <- unique(dt, by = "SNP")
  dt[["SNP"]] <- as.character(dt[["SNP"]])
  dt[["CHR"]] <- coerce_numeric(dt[["CHR"]])
  dt[["BP"]] <- coerce_numeric(dt[["BP"]])
  dt[["A1"]] <- toupper(as.character(dt[["A1"]]))
  dt[["A2"]] <- toupper(as.character(dt[["A2"]]))
  dt[["BETA"]] <- coerce_numeric(dt[["BETA"]])
  dt[["SE"]] <- coerce_numeric(dt[["SE"]])
  dt[["P"]] <- coerce_numeric(dt[["P"]])
  dt[["N"]] <- coerce_numeric(dt[["N"]])
  dt[["EAF"]] <- NA_real_
  dt[["Phenotype"]] <- phenotype
  dt[["id"]] <- id
  as.data.frame(dt)
}

prepare_sumstats_with_frq <- function(path, phenotype, id) {
  dt <- fread(file = path, select = c("SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"))
  dt <- unique(dt, by = "SNP")
  dt[["SNP"]] <- as.character(dt[["SNP"]])
  dt[["CHR"]] <- coerce_numeric(dt[["CHR"]])
  dt[["BP"]] <- coerce_numeric(dt[["BP"]])
  dt[["A1"]] <- toupper(as.character(dt[["A1"]]))
  dt[["A2"]] <- toupper(as.character(dt[["A2"]]))
  dt[["FRQ"]] <- coerce_numeric(dt[["FRQ"]])
  dt[["BETA"]] <- coerce_numeric(dt[["BETA"]])
  dt[["SE"]] <- coerce_numeric(dt[["SE"]])
  dt[["P"]] <- coerce_numeric(dt[["P"]])
  dt[["N"]] <- coerce_numeric(dt[["N"]])
  dt[["EAF"]] <- dt[["FRQ"]]
  dt[["Phenotype"]] <- phenotype
  dt[["id"]] <- id
  as.data.frame(dt)
}

run_standard_mvmr <- function(exposure_list, outcome_df, outcome_label, output_prefix, output_dir) {
  message("Preparing multivariable exposures for ", output_prefix, " ...")
  
  snp_tracking <- data.frame(Step = character(), Count = integer(), stringsAsFactors = FALSE)
  
  # Ensure SCZ follows the strict 5e-8 threshold if it exists in the exposure_list
  for (i in seq_along(exposure_list)) {
    if (unique(exposure_list[[i]]$Phenotype) == "Schizophrenia") {
      message("  Applying strict 5e-8 P-value threshold for Schizophrenia exposure selection...")
      # Temporarily inflate P values > 5e-8 to 1 so they aren't picked as primary IVs for SCZ
      # but are still available as proxies/background for other traits.
      exposure_list[[i]]$P <- ifelse(exposure_list[[i]]$P > 5e-8, 1, exposure_list[[i]]$P)
    }
  }

  mv_exposure <- TwoSampleMR::mv_extract_exposures_local(
    filenames_exposure = exposure_list,
    phenotype_col = rep("Phenotype", length(exposure_list)),
    snp_col = rep("SNP", length(exposure_list)),
    beta_col = rep("BETA", length(exposure_list)),
    se_col = rep("SE", length(exposure_list)),
    eaf_col = rep("EAF", length(exposure_list)),
    effect_allele_col = rep("A1", length(exposure_list)),
    other_allele_col = rep("A2", length(exposure_list)),
    pval_col = rep("P", length(exposure_list)),
    samplesize_col = rep("N", length(exposure_list)),
    id_col = rep("id", length(exposure_list)),
    pval_threshold = pval_threshold,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    harmonise_strictness = harmonise_strictness
  )

  snp_tracking <- rbind(snp_tracking, data.frame(Step = "After clumping & extract", Count = length(unique(mv_exposure$SNP))))

  # [CRITICAL FIX]: "Exposure Intersection"
  expected_exposure_count <- length(exposure_list)
  snp_counts <- table(mv_exposure$SNP)
  valid_snps <- names(snp_counts)[snp_counts == expected_exposure_count]
  
  mv_exposure <- mv_exposure[mv_exposure$SNP %in% valid_snps, ]
  mv_exposure <- mv_exposure[!is.na(mv_exposure$beta.exposure), ]

  snp_tracking <- rbind(snp_tracking, data.frame(Step = "After enforcing Exposure Intersection (no NAs)", Count = length(unique(mv_exposure$SNP))))

  exposure_map <- as.data.frame(unique(mv_exposure[, c("id.exposure", "exposure")]), stringsAsFactors = FALSE)

  message("Formatting outcome for ", output_prefix, " ...")
  mv_outcome <- TwoSampleMR::format_data(
    dat = outcome_df,
    type = "outcome",
    snps = mv_exposure$SNP,
    snp_col = "SNP",
    chr_col = "CHR",
    pos_col = "BP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "FRQ",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    samplesize_col = "N"
  )

  message("Harmonising data for ", output_prefix, " ...")
  mv_dat <- TwoSampleMR::mv_harmonise_data(
    exposure_dat = mv_exposure,
    outcome_dat = mv_outcome
  )

  valid_harmonised <- nrow(mv_dat[["exposure_beta"]])
  snp_tracking <- rbind(snp_tracking, data.frame(Step = "After full MVMR harmonisation", Count = valid_harmonised))

  write.csv(snp_tracking, file.path(output_dir, paste0(output_prefix, "_snp_tracking.csv")), row.names = FALSE)

  if (valid_harmonised < 3) {
    message("Warning: Too few SNPs (< 3) after harmonising. Skipping.")
    return(list(
      results = data.frame(exposure=character(), outcome=character(), method=character()),
      summary = data.frame()
    ))
  }

  message("Running standard IVW-MVMR for ", output_prefix, " ...")
  mv_fit <- TwoSampleMR::mv_multiple(
    mv_dat,
    pval_threshold = mv_multiple_pval_threshold
  )

  result_df <- mv_fit$result %>%
    as.data.frame() %>%
    left_join(exposure_map, by = "id.exposure")

  result_df$exposure <- result_df$exposure.y
  result_df$outcome <- outcome_label
  result_df$method <- "IVW-MVMR"
  result_df <- result_df %>%
    select(-any_of(c("exposure.x", "exposure.y"))) %>%
    relocate(all_of(c("outcome", "method")), .after = "exposure")

  result_or <- TwoSampleMR::generate_odds_ratios(result_df)

  mvmr_input <- MVMR::format_mvmr(
    BXGs = mv_dat[["exposure_beta"]],
    BYG = mv_dat[["outcome_beta"]],
    seBXGs = mv_dat[["exposure_se"]],
    seBYG = mv_dat[["outcome_se"]],
    RSID = rownames(mv_dat[["exposure_beta"]])
  )

  strength_df <- MVMR::strength_mvmr(mvmr_input)
  pleiotropy_df <- MVMR::pleiotropy_mvmr(mvmr_input)

  strength_out <- data.frame(
    exposure = result_df$exposure,
    conditional_F = as.numeric(strength_df[1, seq_len(nrow(result_df))]),
    stringsAsFactors = FALSE
  )

  pleiotropy_out <- as.data.frame(pleiotropy_df)
  instrument_out <- data.frame(
    SNP = rownames(mv_dat[["exposure_beta"]]),
    stringsAsFactors = FALSE
  )

  write.csv(result_or, file.path(output_dir, paste0(output_prefix, "_results.csv")), row.names = FALSE)
  write.csv(strength_out, file.path(output_dir, paste0(output_prefix, "_strength.csv")), row.names = FALSE)
  write.csv(pleiotropy_out, file.path(output_dir, paste0(output_prefix, "_pleiotropy.csv")), row.names = TRUE)
  write.csv(instrument_out, file.path(output_dir, paste0(output_prefix, "_instruments.csv")), row.names = FALSE)

  summary_out <- result_or %>%
    left_join(strength_out, by = "exposure") %>%
    mutate(instrument_count = nrow(instrument_out))

  list(
    results = result_or,
    strength = strength_out,
    pleiotropy = pleiotropy_out,
    instruments = instrument_out,
    summary = summary_out
  )
}

script_path <- get_script_path()
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!is.na(script_path)) {
  script_dir <- dirname(script_path)
} else {
  script_dir <- file.path(project_root, "MVMR_SCZ_PTSD_3tasks")
}

output_dir <- file.path(project_root, "MVMR_SCZ_PTSD_3tasks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

task_configs <- data.frame(
  task = c("PairMatchingA2", "symDigital", "TowerRearrange"),
  phenotype = c("PairMatchingA2", "symDigital", "TowerRearrange"),
  id = c("pairmatching_a2_fixed", "symdigital_fixed", "tower_rearrange_fixed"),
  source_file = c("Memory_extracted_full.txt", "SymbolDigit_extracted_full.txt", "Tower_extracted_full.txt"),
  stringsAsFactors = FALSE
)

manifest <- bind_rows(
  task_configs %>%
    transmute(
      analysis = paste0("SCZ_", task, "_lefthands"),
      context = "SCZ + task -> Left-handedness",
      exposure_1 = "Schizophrenia",
      exposure_2 = phenotype,
      outcome = "Left-handedness",
      task,
      source_file
    ),
  task_configs %>%
    transmute(
      analysis = paste0("lefthands_", task, "_PTSD"),
      context = "Left-handedness + task -> PTSD",
      exposure_1 = "Left-handedness",
      exposure_2 = phenotype,
      outcome = "PTSD",
      task,
      source_file
    )
)

write.csv(manifest, file.path(output_dir, "analysis_manifest.csv"), row.names = FALSE)

scz_df <- prepare_sumstats_no_frq(
  file.path(project_root, "SCZ_lastest_tb.txt"),
  phenotype = "Schizophrenia",
  id = "SCZ_fixed"
)

lefthands_df <- prepare_sumstats_with_frq(
  file.path(project_root, "lefthands_tb.txt"),
  phenotype = "Lefthands",
  id = "lefthands_fixed"
)

ptsd_df <- prepare_sumstats_with_frq(
  file.path(project_root, "PTSD_tb.txt"),
  phenotype = "PTSD",
  id = "PTSD_fixed"
)

all_summaries <- list()

for (row_index in seq_len(nrow(task_configs))) {
  task_row <- task_configs[row_index, ]
  task_df <- prepare_sumstats_no_frq(
    file.path(project_root, task_row$source_file),
    phenotype = task_row$phenotype,
    id = task_row$id
  )

  scz_prefix <- paste0("strict5e8_MVMR_SCZ_", task_row$task, "_lefthands")
  scz_fit <- run_standard_mvmr(
    exposure_list = list(scz_df, task_df),
    outcome_df = lefthands_df,
    outcome_label = "Left-handedness",
    output_prefix = scz_prefix,
    output_dir = output_dir
  )

  if (nrow(scz_fit$summary) > 0) {
    all_summaries[[length(all_summaries) + 1]] <- scz_fit$summary %>%
      mutate(
        analysis = scz_prefix,
        context = "SCZ + task -> Left-handedness",
        task = task_row$task,
        task_source_file = task_row$source_file
      )
  }

  ptsd_prefix <- paste0("strict5e8_MVMR_lefthands_", task_row$task, "_PTSD")
  ptsd_fit <- run_standard_mvmr(
    exposure_list = list(lefthands_df, task_df),
    outcome_df = ptsd_df,
    outcome_label = "PTSD",
    output_prefix = ptsd_prefix,
    output_dir = output_dir
  )

  if (nrow(ptsd_fit$summary) > 0) {
    all_summaries[[length(all_summaries) + 1]] <- ptsd_fit$summary %>%
      mutate(
        analysis = ptsd_prefix,
        context = "Left-handedness + task -> PTSD",
        task = task_row$task,
        task_source_file = task_row$source_file
      )
  }
}

if (length(all_summaries) > 0) {
  combined_summary <- bind_rows(all_summaries)
  
  write.csv(combined_summary, file.path(output_dir, "combined_mvmr_results_strictSCZ.csv"), row.names = FALSE)
  message("Batch completed with strict SCZ threshold. Combined result saved to: combined_mvmr_results_strictSCZ.csv")
} else {
  message("No valid results to combine for any task.")
}

  if (nrow(ptsd_fit$summary) > 0) {
    all_summaries[[length(all_summaries) + 1]] <- ptsd_fit$summary %>%
      mutate(
        analysis = ptsd_prefix,
        context = "Left-handedness + task -> PTSD",
        task = task_row$task,
        task_source_file = task_row$source_file
      )
  }
}

combined_results <- bind_rows(all_summaries) %>%
  relocate(analysis, context, task, task_source_file, exposure, outcome, method)

write.csv(combined_results, file.path(output_dir, "combined_mvmr_results.csv"), row.names = FALSE)

message("All six MVMR analyses finished. Outputs written to: ", output_dir)