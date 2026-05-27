# ------ packages ------
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(MendelianRandomization)
library(dplyr)

safe_numeric <- function(x) {
  if (length(x) == 0 || is.null(x)) {
    return(NA_real_)
  }

  if (is.numeric(x)) {
    return(as.numeric(x[[1]]))
  }

  cleaned <- gsub("^<", "", as.character(x[[1]]))
  suppressWarnings(as.numeric(cleaned))
}

# ------ helper: read & format exposure ------
prepare_exposure <- function(dat, has_eaf = TRUE) {
  p_thresh <- 5e-6
  sub <- subset(dat, dat$P < p_thresh)

  fmt_args <- list(
    dat               = sub,
    type              = "exposure",
    header            = TRUE,
    snp_col           = "SNP",
    chr_col           = "CHR",
    pos_col           = "BP",
    beta_col          = "BETA",
    se_col            = "SE",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    pval_col          = "P",
    samplesize_col    = "N"
  )
  if (has_eaf) fmt_args$eaf_col <- "FRQ"

  exp <- do.call(format_data, fmt_args)
  exp <- clump_data(exp, clump_r2 = 0.001)

  exp$R2   <- get_r_from_bsen(exp$beta.exposure,
                               exp$se.exposure,
                               exp$samplesize.exposure)^2
  exp$Fval <- (exp$samplesize.exposure - 2) * exp$R2 / (1 - exp$R2)

  # Fix: filter on Fval (not the non-existent column F)
  exp <- exp[exp$Fval > 10, ]
  exp
}

# ------ helper: format outcome ------
prepare_outcome <- function(dat, snps, has_eaf = TRUE) {
  fmt_args <- list(
    dat               = dat,
    type              = "outcome",
    snps              = snps,
    header            = TRUE,
    snp_col           = "SNP",
    chr_col           = "CHR",
    pos_col           = "BP",
    beta_col          = "BETA",
    se_col            = "SE",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    pval_col          = "P",
    samplesize_col    = "N"
  )
  if (has_eaf) fmt_args$eaf_col <- "FRQ"
  do.call(format_data, fmt_args)
}

# ------ helper: run full MR pipeline for one pair ------
run_mr_pair <- function(exp_dat, out_dat, exp_name, out_name) {
  harm <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  harm$id.exposure <- exp_name
  harm$id.outcome  <- out_name
  harm <- harm[harm$mr_keep == "TRUE", ]

  if ("Fval" %in% colnames(exp_dat)) {
    cols_to_add <- intersect(c("SNP", "Fval", "R2"), colnames(exp_dat))
    harm <- merge(harm, exp_dat[, cols_to_add, drop = FALSE], by = "SNP", all.x = TRUE)
  }

  if (nrow(harm) < 3) {
    message(sprintf("Skipping %s -> %s: fewer than 3 valid IVs after harmonisation",
                    exp_name, out_name))
    return(NULL)
  }

  # MRPRESSO
  presso <- tryCatch(
    mr_presso(
      BetaOutcome    = "beta.outcome",
      BetaExposure   = "beta.exposure",
      SdOutcome      = "se.outcome",
      SdExposure     = "se.exposure",
      data           = harm,
      OUTLIERtest    = TRUE,
      DISTORTIONtest = TRUE,
      NbDistribution = 10000,
      SignifThreshold = 0.05
    ),
    error = function(e) {
      message("MR-PRESSO failed: ", e$message); NULL
    }
  )

  het    <- mr_heterogeneity(harm)
  het$I2 <- abs(100 * (het$Q - het$Q_df) / het$Q)

  pleio  <- mr_pleiotropy_test(harm)
  direct <- directionality_test(harm)

  mr_res <- mr(harm)
  or_res <- generate_odds_ratios(mr_res)

  # Extract all 4 main methods by name (robust to row ordering)
  ivw     <- or_res %>% dplyr::filter(method == "Inverse variance weighted") %>% mutate(method_label = "IVW")
  mregger <- or_res %>% dplyr::filter(method == "MR Egger")                  %>% mutate(method_label = "MR-Egger")
  wme     <- or_res %>% dplyr::filter(method == "Weighted median")            %>% mutate(method_label = "Weighted_Median")
  wmode   <- or_res %>% dplyr::filter(method == "Weighted mode")              %>% mutate(method_label = "Weighted_Mode")

  # Leave-One-Out sensitivity analysis
  loo_res <- tryCatch(
    {
      loo <- mr_leaveoneout(harm)
      loo$pair <- paste0(exp_name, " -> ", out_name)
      loo
    },
    error = function(e) { message("LOO failed: ", e$message); NULL }
  )

  # Compile PRESSO row
  presso_row <- data.frame(
    presso_estimate = NA_real_,
    presso_sd = NA_real_,
    presso_tstat = NA_real_,
    presso_pval = NA_real_,
    presso_global_pval = NA_real_,
    presso_distortion_pval = NA_real_,
    exposure = exp_name,
    outcome = out_name,
    stringsAsFactors = FALSE
  )

  if (!is.null(presso)) {
    presso_main <- as.data.frame(presso$`Main MR results`)
    if (nrow(presso_main) >= 2) {
      presso_row$presso_estimate <- safe_numeric(presso_main$`Causal Estimate`[2])
      presso_row$presso_sd <- safe_numeric(presso_main$Sd[2])
      presso_row$presso_tstat <- safe_numeric(presso_main$`T-stat`[2])
      presso_row$presso_pval <- safe_numeric(presso_main$`P-value`[2])
    }

    presso_row$presso_global_pval <- safe_numeric(
      presso$`MR-PRESSO results`$`Global Test`$Pvalue
    )
    presso_row$presso_distortion_pval <- safe_numeric(
      presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
    )
  }

  # Pleiotropy row
  pleio_row <- pleio[, c("egger_intercept", "se", "pval")]
  colnames(pleio_row) <- c("egger_intercept", "egger_se", "egger_intercept_pval")

  # Heterogeneity (IVW row)
  het_ivw <- het[het$method == "Inverse variance weighted", c("Q", "Q_df", "Q_pval", "I2"), drop = FALSE]

  # Directionality
  direct_row <- direct[, c("correct_causal_direction", "steiger_pval")]

  list(
    harm      = harm,
    ivw       = ivw,
    mregger   = mregger,
    wme       = wme,
    wmode     = wmode,
    loo       = loo_res,
    het       = het_ivw,
    pleio     = pleio_row,
    direct    = direct_row,
    presso    = presso_row
  )
}

# ------ load raw data ------
ADHD         <- read.table("ADHD_tb.txt",        header = TRUE)
PTSD         <- read.table("PTSD_tb.txt",        header = TRUE)
Schizophrenia <- read.table("SCZ_lastest_tb.txt", header = TRUE)
lefthands    <- read.table("lefthands_tb.txt",   header = TRUE)

# ------ prepare exposures ------
# has_eaf: ADHD=TRUE, PTSD=TRUE, SCZ=FALSE (no FRQ col), lefthands=TRUE
ADHD_exposure  <- prepare_exposure(ADHD,          has_eaf = TRUE)
PTSD_exposure  <- prepare_exposure(PTSD,          has_eaf = TRUE)
SCZ_exposure   <- prepare_exposure(Schizophrenia, has_eaf = FALSE)
LH_exposure    <- prepare_exposure(lefthands,     has_eaf = TRUE)

disease_exposures <- list(
  ADHD          = ADHD_exposure,
  PTSD          = PTSD_exposure,
  Schizophrenia = SCZ_exposure
)

# ------ Direction 1: Disease → Handedness ------
message("=== Direction 1: Disease -> Handedness ===")

fwd_results <- list()

for (dis in names(disease_exposures)) {
  exp    <- disease_exposures[[dis]]
  out    <- prepare_outcome(lefthands, snps = exp$SNP, has_eaf = TRUE)
  result <- run_mr_pair(exp, out, exp_name = dis, out_name = "lefthands")
  if (!is.null(result)) fwd_results[[dis]] <- result
}

# ------ Direction 2: Handedness → Disease ------
message("=== Direction 2: Handedness -> Disease ===")

rev_results <- list()

has_eaf_map <- list(ADHD = TRUE, PTSD = TRUE, Schizophrenia = FALSE)
raw_data_map <- list(ADHD = ADHD, PTSD = PTSD, Schizophrenia = Schizophrenia)

for (dis in names(raw_data_map)) {
  out    <- prepare_outcome(raw_data_map[[dis]],
                            snps    = LH_exposure$SNP,
                            has_eaf = has_eaf_map[[dis]])
  result <- run_mr_pair(LH_exposure, out, exp_name = "lefthands", out_name = dis)
  if (!is.null(result)) rev_results[[dis]] <- result
}

# ------ assemble output tables ------
collect_table <- function(results_list, field) {
  rows <- lapply(names(results_list), function(nm) {
    r <- results_list[[nm]][[field]]
    if (is.null(r)) return(NULL)
    r$pair_key <- nm
    r
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

# IVW / MR-Egger / Weighted Median
make_mr_table <- function(results_list) {
  rows <- lapply(names(results_list), function(nm) {
    r <- results_list[[nm]]
    bind_rows(r$ivw, r$mregger, r$wme, r$wmode)
  })
  bind_rows(rows)
}

fwd_mr   <- make_mr_table(fwd_results)
rev_mr   <- make_mr_table(rev_results)
all_mr   <- bind_rows(
  fwd_mr %>% mutate(direction = "disease_to_handedness"),
  rev_mr %>% mutate(direction = "handedness_to_disease")
)

# Heterogeneity
fwd_het  <- collect_table(fwd_results, "het")  %>% mutate(direction = "disease_to_handedness")
rev_het  <- collect_table(rev_results, "het")  %>% mutate(direction = "handedness_to_disease")
all_het  <- bind_rows(fwd_het, rev_het)

# Pleiotropy
fwd_pleio <- collect_table(fwd_results, "pleio") %>% mutate(direction = "disease_to_handedness")
rev_pleio <- collect_table(rev_results, "pleio") %>% mutate(direction = "handedness_to_disease")
all_pleio <- bind_rows(fwd_pleio, rev_pleio)

# PRESSO
fwd_presso <- collect_table(fwd_results, "presso") %>%
  mutate(direction = "disease_to_handedness")
rev_presso <- collect_table(rev_results, "presso") %>%
  mutate(direction = "handedness_to_disease")
all_presso <- bind_rows(fwd_presso, rev_presso)

# Directionality
fwd_direct <- collect_table(fwd_results, "direct") %>% mutate(direction = "disease_to_handedness")
rev_direct <- collect_table(rev_results, "direct") %>% mutate(direction = "handedness_to_disease")
all_direct <- bind_rows(fwd_direct, rev_direct)

# SNPs / Instruments (with F values)
fwd_harm <- collect_table(fwd_results, "harm") %>% mutate(direction = "disease_to_handedness")
rev_harm <- collect_table(rev_results, "harm") %>% mutate(direction = "handedness_to_disease")
all_harm <- bind_rows(fwd_harm, rev_harm)

# ------ write CSV outputs ------
write.csv(all_mr,    "MR_diseases_handedness_main.csv",          row.names = FALSE)
write.csv(all_het,   "MR_diseases_handedness_heterogeneity.csv", row.names = FALSE)
write.csv(all_pleio, "MR_diseases_handedness_pleiotropy.csv",    row.names = FALSE)
write.csv(all_presso,"MR_diseases_handedness_presso.csv",        row.names = FALSE)
write.csv(all_direct,"MR_diseases_handedness_direction.csv",     row.names = FALSE)
write.csv(all_harm,  "MR_diseases_handedness_SNPs_Fval.csv",     row.names = FALSE)

# LOO
collect_loo <- function(results_list) {
  rows <- lapply(names(results_list), function(nm) results_list[[nm]]$loo)
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  bind_rows(rows)
}
fwd_loo <- collect_loo(fwd_results)
rev_loo <- collect_loo(rev_results)
all_loo <- bind_rows(
  if (!is.null(fwd_loo)) fwd_loo %>% mutate(direction = "disease_to_handedness"),
  if (!is.null(rev_loo)) rev_loo %>% mutate(direction = "handedness_to_disease")
)
if (!is.null(all_loo) && nrow(all_loo) > 0)
  write.csv(all_loo, "MR_diseases_handedness_loo.csv", row.names = FALSE)

message("Done. Results written to working directory.")
