
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

# ----step1: extract exposure1 and exposure2 ----
# xx1 and xx2 represent exposure variable 1 and 2
MV_exposure <- mv_extract_exposures_local(
    c(paste("./exposure/summary/","xx1","_tb.txt",sep = ""),paste("./exposure/summary/","xx2","_tb.txt",sep = "")),
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    ncase_col = "N",
    pval_threshold = 5e-06
  )
  
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

# ----step2: MVMR ----
  mvdata <- mv_harmonise_data(MV_exposure, OOOO_outcome)
  mvmr_ivw <-mv_multiple(mvdata)
  mvmr_ivw <- mvmr_ivw$result %>% 
    as.data.frame() %>%
    as_tibble() %>% 
    janitor::clean_names() %>%
    mutate(
      exposure = ifelse(exposure == 'exposure1', "xx1", "xx2"), 
      outcome = 'OOOO',
      method = 'IVW-MVMR'
    ) %>%
    relocate(outcome, method, .after = exposure)
  
  mvmr_data <- format_mvmr(BXGs = mvdata[["exposure_beta"]],
                           BYG = mvdata[["outcome_beta"]],
                           seBXGs = mvdata[["exposure_se"]],
                           seBYG = mvdata[["outcome_se"]],
                           RSID = rownames(mvdata[["exposure_beta"]]))
  
  mvmr_strength <- strength_mvmr(mvmr_data)
  mvmr_pleiotropy <- pleiotropy_mvmr(mvmr_data)
  
