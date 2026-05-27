library(GenomicSEM)

base_dir <- "c:/Working files/01 Experiment Data/Publicated Dataset/UKB behavior"
output_dir <- file.path(base_dir, "LDSC_ptsd_cognition")
ld_dir <- file.path(base_dir, "eur_w_ld_chr")
cognition_munged_dir <- file.path(base_dir, "LDSC_cognition")
disease_munged_dir <- file.path(base_dir, "LDSC_disease")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

ptsd_munged_file <- file.path(disease_munged_dir, "PTSD.sumstats.gz")
cog_names <- c("Matrix", "Memory", "RT", "SymbolDigit", "TMTB", "Tower", "VNR") 
cog_munged_files <- file.path(cognition_munged_dir, paste0(cog_names, ".sumstats.gz"))
all_names <- c("PTSD", cog_names)
all_munged_files <- c(ptsd_munged_file, cog_munged_files)

cat("Using existing munged sumstats\n")
cat("Traits:", paste(all_names, collapse = ", "), "\n")

ldsc_out <- ldsc(
  traits = all_munged_files,
  sample.prev = rep(NA, length(all_names)),
  population.prev = rep(NA, length(all_names)),
  ld = ld_dir,
  wld = ld_dir,
  trait.names = all_names,
  ldsc.log = "ldsc_ptsd_cognition_log"
)

S <- ldsc_out$S
V <- ldsc_out$V
I <- ldsc_out$I

colnames(S) <- rownames(S) <- all_names
colnames(I) <- rownames(I) <- all_names

Rg <- cov2cor(S)
colnames(Rg) <- rownames(Rg) <- all_names

vech_idx <- function(i, j, n) {
  if (i < j) {
    tmp <- i
    i <- j
    j <- tmp
  }
  (j - 1) * n - j * (j - 1) / 2 + i
}

n <- length(all_names)
se_rg <- matrix(NA_real_, n, n)
colnames(se_rg) <- rownames(se_rg) <- all_names

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) {
      next
    }
    idx_ij <- vech_idx(i, j, n)
    se_rg[i, j] <- sqrt(V[idx_ij, idx_ij]) / sqrt(abs(S[i, i] * S[j, j]))       
  }
}

ptsd_idx <- which(all_names == "PTSD")
cog_idx <- which(all_names != "PTSD")

rg_vec <- Rg[ptsd_idx, cog_idx]
se_vec <- se_rg[ptsd_idx, cog_idx]
z_vec <- rg_vec / se_vec
p_vec <- 2 * pnorm(-abs(z_vec))
h2_ptsd <- S[ptsd_idx, ptsd_idx]
h2_cog <- diag(S)[cog_idx]
intercept_vec <- I[ptsd_idx, cog_idx]

results <- data.frame(
  Exposure = "PTSD",
  Outcome = all_names[cog_idx],
  rg = round(rg_vec, 4),
  SE = round(se_vec, 4),
  Z = round(z_vec, 3),
  P = signif(p_vec, 4),
  h2_PTSD = round(h2_ptsd, 4),
  h2_Cognition = round(h2_cog, 4),
  CrossIntercept = round(intercept_vec, 4),
  row.names = NULL
)

cat("\nResults LDSC:\n")
print(results, digits = 4)

write.csv(results, "ldsc_ptsd_cognition_rg.csv", row.names = FALSE)
write.csv(as.data.frame(Rg), "ldsc_ptsd_cognition_Rg_full_matrix.csv")
write.csv(as.data.frame(S), "ldsc_ptsd_cognition_S_covariance.csv")
write.csv(as.data.frame(I), "ldsc_ptsd_cognition_intercept.csv")
cat("\nDone\n")