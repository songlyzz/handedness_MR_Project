rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
})


# ─── Config ────────────────────────────────────────────────────────────────────
BASE_DIR   <- "./UKBbehavior"
OUTPUT_DIR <- file.path(BASE_DIR, "MR_Cognition_PTSD")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

pval_threshold   <- 5e-6
clump_r2         <- 0.001
clump_kb         <- 10000
f_stat_threshold <- 10     # Weak instrument threshold: F < 10 are removed

coerce_numeric <- function(x) suppressWarnings(as.numeric(x))

# Safely parse p-values that may come back as strings like "<0.001" or ">0.999"
parse_pval <- function(x) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) return(NA_real_)
  x_chr <- as.character(x)
  cleaned <- sub("^[<>]=?\\s*", "", x_chr)
  suppressWarnings(as.numeric(cleaned))
}

# ─── Loaders ───────────────────────────────────────────────────────────────────
load_no_frq <- function(path, phenotype, id) {
  dt <- fread(file = path, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
  dt <- unique(dt, by = "SNP")
  dt[, SNP  := as.character(SNP)]
  dt[, CHR  := coerce_numeric(CHR)]
  dt[, BP   := coerce_numeric(BP)]
  dt[, A1   := toupper(A1)]
  dt[, A2   := toupper(A2)]
  dt[, BETA := coerce_numeric(BETA)]
  dt[, SE   := coerce_numeric(SE)]
  dt[, P    := coerce_numeric(P)]
  dt[, N    := coerce_numeric(N)]
  dt[, EAF  := NA_real_]
  dt[, Phenotype := phenotype]
  dt[, id   := id]
  as.data.frame(dt)
}

load_with_frq <- function(path, phenotype, id) {
  dt <- fread(file = path, select = c("SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"))
  dt <- unique(dt, by = "SNP")
  dt[, SNP  := as.character(SNP)]
  dt[, CHR  := coerce_numeric(CHR)]
  dt[, BP   := coerce_numeric(BP)]
  dt[, A1   := toupper(A1)]
  dt[, A2   := toupper(A2)]
  dt[, FRQ  := coerce_numeric(FRQ)]
  dt[, BETA := coerce_numeric(BETA)]
  dt[, SE   := coerce_numeric(SE)]
  dt[, P    := coerce_numeric(P)]
  dt[, N    := coerce_numeric(N)]
  dt[, EAF  := FRQ]
  dt[, Phenotype := phenotype]
  dt[, id   := id]
  as.data.frame(dt)
}

# ─── Single-direction MR ───────────────────────────────────────────────────────
run_mr <- function(exposure_raw, outcome_raw, exposure_name, outcome_name) {
  message("  [MR] ", exposure_name, " -> ", outcome_name)

  # 1. Format exposure
  exp_fmt <- format_data(
    dat               = exposure_raw,
    type              = "exposure",
    snp_col           = "SNP",
    chr_col           = "CHR",
    pos_col           = "BP",
    beta_col          = "BETA",
    se_col            = "SE",
    eaf_col           = "EAF",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    pval_col          = "P",
    samplesize_col    = "N"
  )
  exp_fmt <- exp_fmt[exp_fmt$pval.exposure <= pval_threshold, ]

  n_sig <- nrow(exp_fmt)
  message("    Significant SNPs before clumping: ", n_sig)
  if (n_sig == 0) {
    message("    No significant SNPs. Skipping.")
    return(NULL)
  }

  # 2. Clump
  exp_clumped <- tryCatch(
    clump_data(exp_fmt, clump_r2 = clump_r2, clump_kb = clump_kb, pop = "EUR"),
    error = function(e) {
      message("    Clumping error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(exp_clumped) || nrow(exp_clumped) == 0) {
    message("    No SNPs after clumping. Skipping.")
    return(NULL)
  }
  n_iv_preFilt <- nrow(exp_clumped)
  message("    Instruments after clumping: ", n_iv_preFilt)

  # F-statistic per SNP: F = (beta/se)^2
  exp_clumped$F_stat <- (exp_clumped$beta.exposure / exp_clumped$se.exposure)^2
  mean_F_all   <- mean(exp_clumped$F_stat, na.rm = TRUE)
  n_weak       <- sum(exp_clumped$F_stat < f_stat_threshold, na.rm = TRUE)
  exp_clumped  <- exp_clumped[!is.na(exp_clumped$F_stat) & exp_clumped$F_stat >= f_stat_threshold, ]
  n_iv         <- nrow(exp_clumped)
  mean_F_final <- if (n_iv > 0) mean(exp_clumped$F_stat, na.rm = TRUE) else NA_real_
  message("    After F-stat filter (F>=10): ", n_iv, " instruments kept, ", n_weak, " removed  |  mean F (all) = ", round(mean_F_all, 1), ",  mean F (kept) = ", round(mean_F_final, 1))

  if (n_iv == 0) {
    message("    No instruments remain after F-stat filtering. Skipping.")
    return(NULL)
  }

  # 3. Format outcome
  out_fmt <- format_data(
    dat               = outcome_raw,
    type              = "outcome",
    snps              = exp_clumped$SNP,
    snp_col           = "SNP",
    chr_col           = "CHR",
    pos_col           = "BP",
    beta_col          = "BETA",
    se_col            = "SE",
    eaf_col           = "EAF",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    pval_col          = "P",
    samplesize_col    = "N"
  )
  if (is.null(out_fmt) || nrow(out_fmt) == 0) {
    message("    No matching SNPs found in outcome. Skipping.")
    return(NULL)
  }

  # 4. Harmonise
  dat <- harmonise_data(exp_clumped, out_fmt, action = 2)
  dat <- dat[dat$mr_keep == TRUE, ]
  n_harm <- nrow(dat)
  message("    SNPs after harmonisation: ", n_harm)

  if (n_harm < 3) {
    message("    Too few SNPs (<3). Skipping.")
    return(NULL)
  }

  # 5. MR methods: IVW + Egger + Weighted Median + Weighted Mode
  res <- mr(
    dat,
    method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  )
  res_or <- generate_odds_ratios(res)

  # Leave-One-Out sensitivity analysis
  loo_res <- tryCatch(
    {
      loo <- mr_leaveoneout(dat)
      loo$direction <- paste0(exposure_name, " -> ", outcome_name)
      loo
    },
    error = function(e) {
      message("    LOO error: ", conditionMessage(e))
      NULL
    }
  )

  # Egger intercept (pleiotropy test)
  egger_int <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  egger_intercept <- if (!is.null(egger_int) && nrow(egger_int) > 0) egger_int$egger_intercept[1] else NA_real_
  egger_int_pval  <- if (!is.null(egger_int) && nrow(egger_int) > 0) egger_int$pval[1]             else NA_real_

  # Heterogeneity Q
  het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  q_ivw        <- if (!is.null(het) && any(het$method == "Inverse variance weighted")) het$Q[het$method == "Inverse variance weighted"][1]      else NA_real_
  q_ivw_pval   <- if (!is.null(het) && any(het$method == "Inverse variance weighted")) het$Q_pval[het$method == "Inverse variance weighted"][1] else NA_real_
  q_egger      <- if (!is.null(het) && any(het$method == "MR Egger"))                 het$Q[het$method == "MR Egger"][1]                        else NA_real_
  q_egger_pval <- if (!is.null(het) && any(het$method == "MR Egger"))                 het$Q_pval[het$method == "MR Egger"][1]                   else NA_real_

  # 6. MR-PRESSO
  presso_global_pval    <- NA_real_
  presso_corrected_b    <- NA_real_
  presso_corrected_se   <- NA_real_
  presso_corrected_pval <- NA_real_
  presso_n_outliers     <- NA_integer_

  if (has_presso && n_harm >= 4) {
    presso_res <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome     = "beta.outcome",
        BetaExposure    = "beta.exposure",
        SdOutcome       = "se.outcome",
        SdExposure      = "se.exposure",
        OUTLIERtest     = TRUE,
        DISTORTIONtest  = TRUE,
        data            = dat,
        NbDistribution  = 10000,
        SignifThreshold = 0.05
      ),
      error = function(e) {
        message("    MR-PRESSO error: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(presso_res)) {
      presso_global_pval <- tryCatch(parse_pval(presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]), error = function(e) NA_real_)
      raw_res <- tryCatch(presso_res[["MR-PRESSO results"]][["Main MR results"]], error = function(e) NULL)
      if (!is.null(raw_res) && is.data.frame(raw_res)) {
        if ("MR Analysis" %in% names(raw_res)) {
          corr_idx <- grepl("Outlier-corrected", raw_res[["MR Analysis"]])
        } else {
          corr_idx <- grepl("Outlier-corrected", rownames(raw_res))
        }
        corr_row <- raw_res[corr_idx, , drop = FALSE]
        if (is.data.frame(corr_row) && nrow(corr_row) > 0) {
          presso_corrected_b    <- tryCatch(as.numeric(corr_row[["Causal Estimate"]][1]), error = function(e) NA_real_)
          presso_corrected_se   <- tryCatch(as.numeric(corr_row[["Sd"]][1]),              error = function(e) NA_real_)
          presso_corrected_pval <- tryCatch(parse_pval(corr_row[["P-value"]][1]),         error = function(e) NA_real_)
        }
      }
      outliers <- tryCatch(presso_res[["MR-PRESSO results"]][["Outlier Test"]], error = function(e) NULL)
      if (!is.null(outliers) && is.data.frame(outliers) && "Pvalue" %in% names(outliers))
        presso_n_outliers <- sum(sapply(outliers[["Pvalue"]], parse_pval) < 0.05, na.rm = TRUE)
    }
  }

  res_or$direction               <- paste0(exposure_name, " -> ", outcome_name)
  res_or$exposure_label          <- exposure_name
  res_or$outcome_label           <- outcome_name
  res_or$n_sig_snps              <- n_sig
  res_or$n_instruments_preFilter <- n_iv_preFilt
  res_or$n_weak_instruments      <- n_weak
  res_or$n_instruments           <- n_iv
  res_or$mean_F                  <- round(mean_F_final, 2)
  res_or$n_snps_harmonised       <- n_harm
  res_or$egger_intercept         <- egger_intercept
  res_or$egger_intercept_pval    <- egger_int_pval
  res_or$Q_ivw                   <- q_ivw
  res_or$Q_ivw_pval              <- q_ivw_pval
  res_or$Q_egger                 <- q_egger
  res_or$Q_egger_pval            <- q_egger_pval
  res_or$presso_global_pval      <- presso_global_pval
  res_or$presso_n_outliers       <- presso_n_outliers
  res_or$presso_corrected_b      <- presso_corrected_b
  res_or$presso_corrected_se     <- presso_corrected_se
  res_or$presso_corrected_pval   <- presso_corrected_pval

  list(results = res_or, loo = loo_res)
}

# ─── Task definitions ─────────────────────────────────────────────────────────
task_configs <- list(
  list(name = "Matrix",      file = "Matrix_extracted_full.txt",      id = "matrix_cog"),
  list(name = "Memory",      file = "Memory_extracted_full.txt",      id = "memory_cog"),
  list(name = "RT",          file = "RT_extracted_full.txt",          id = "rt_cog"),
  list(name = "SymbolDigit", file = "SymbolDigit_extracted_full.txt", id = "symboldigit_cog"),
  list(name = "TMTB",        file = "TMTB_extracted_full.txt",        id = "tmtb_cog"),
  list(name = "Tower",       file = "Tower_extracted_full.txt",       id = "tower_cog"),
  list(name = "VNR",         file = "VNR_extracted_full.txt",         id = "vnr_cog")
)

# ─── Load PTSD ────────────────────────────────────────────────────────────────
message("Loading PTSD sumstats...")
ptsd_df <- load_with_frq(
  file.path(BASE_DIR, "PTSD_tb.txt"),
  phenotype = "PTSD",
  id        = "ptsd"
)
message("  Loaded ", nrow(ptsd_df), " variants for PTSD.")

# ─── Main loop: 7 cognitive tasks -> PTSD (7 analyses) ────────────────────────
all_results <- list()all_loo     <- list()
for (task in task_configs) {
  message("\n=== ", task$name, " -> PTSD ===")

  cog_df <- load_no_frq(
    file.path(BASE_DIR, task$file),
    phenotype = task$name,
    id        = task$id
  )
  message("  Loaded ", nrow(cog_df), " variants for ", task$name)

  run_res <- run_mr(cog_df, ptsd_df, task$name, "PTSD")
  if (!is.null(run_res)) {
    all_results[[length(all_results) + 1]] <- run_res$results
    if (!is.null(run_res$loo)) all_loo[[length(all_loo) + 1]] <- run_res$loo
  }

  rm(cog_df)
  gc()
}

# ─── Combine & save ───────────────────────────────────────────────────────────
if (length(all_results) > 0) {
  combined <- bind_rows(all_results) %>%
    select(
      direction, exposure_label, outcome_label, method,
      n_sig_snps, n_instruments, n_snps_harmonised, nsnp,
      b, se, pval, lo_ci, up_ci,
      or, or_lci95, or_uci95,
      egger_intercept, egger_intercept_pval,
      Q_ivw, Q_ivw_pval, Q_egger, Q_egger_pval,
      everything()
    )

  out_path <- file.path(OUTPUT_DIR, "mr_cognition_ptsd.csv")
  write.csv(combined, out_path, row.names = FALSE)

  if (length(all_loo) > 0) {
    combined_loo <- bind_rows(all_loo)
    loo_path <- file.path(OUTPUT_DIR, "mr_cognition_ptsd_loo.csv")
    write.csv(combined_loo, loo_path, row.names = FALSE)
    message("LOO saved to: ", loo_path)
  }

  message("\n=== Completed! ===")
  message("Analyses: ", length(unique(combined$direction)))
  message("Total rows (analysis x method): ", nrow(combined))
  message("Saved to: ", out_path)
} else {
  message("No valid MR results were produced.")
}
