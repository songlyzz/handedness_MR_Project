
# ------ package ------
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(cowplot)
library(gprofiler2)
library(MVMR)
library(ieugwasr)
library(gwasglue)
library(MendelianRandomization)
library(openxlsx)

# ------ run -----
# ----step1: psychiatric disease to handedness----
# XXXX is exposure abbreviation
XXXX <- read.table("./exposure/summary/XXXX_tb.txt", header = T)
XXXX_exposure <- subset(XXXX, ADHD$P<5e-6)
XXXX_exposure <- format_data(
  dat = XXXX_exposure,
  type = "exposure",
  header = TRUE,
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
XXXX_exposure <- clump_data(XXXX_exposure,
                            clump_r2 = 0.001)
XXXX_exposure$R2 <- get_r_from_bsen(XXXX_exposure$beta.exposure, XXXX_exposure$se.exposure, XXXX_exposure$samplesize.exposure)^2
XXXX_exposure$Fval <- (XXXX_exposure$samplesize.exposure - 2) * XXXX_exposure$R2 / (1 - XXXX_exposure$R2)
# filter F > 10 instrument
XXXX_exposure <- XXXX_exposure[XXXX_exposure$F > 10,]
summary(XXXX_exposure$Fval)

# OOOO is outcome abbreviation
OOOO <- read.table("./outcome/summary/OOOO_tb.txt", header = T)
OOOO_outcome <- format_data(
    dat = OOOO,
    type = "outcome",
    snps = XXXX_exposure$SNP,
    header = TRUE,
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

# ----step2: MR ----
SD_mr_data <- harmonise_data(exposure_dat = XXXX_exposure,
                              outcome_dat = XXXX_outcome )
SD_mr_data <- SD_mr_data[SD_mr_data$mr_keep=="TRUE",]

SD_presso <- mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data=SD_mr_data, 
                      OUTLIERtest = TRUE,
                      DISTORTIONtest = TRUE,
                      NbDistribution = 10000,
                      SignifThreshold = 0.05)
SD_heterogeneity <- mr_heterogeneity(SD_mr_data)
SD_pleiotropy <- mr_pleiotropy_test(SD_mr_data)
SD_direct <- directionality_test(SD_mr_data)
    
SD_mr_results <- mr(SD_mr_data)
SD_or_results <- generate_odds_ratios(SD_mr_results)

