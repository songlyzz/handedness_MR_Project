library(GenomicSEM)
library(data.table)

base_dir <- "./UKBbehavior"
output_dir <- file.path(base_dir, "LDSC_schizophrenia_cognition")
ld_dir <- file.path(base_dir, "eur_w_ld_chr")
hm3 <- file.path(base_dir, "w_hm3.snplist.gz")
cognition_munged_dir <- file.path(base_dir, "LDSC_cognition")
disease_munged_dir <- file.path(base_dir, "LDSC_disease")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

scz_raw_file <- file.path(base_dir, "SCZ_lastest_tb.txt")
scz_munged_file <- file.path(disease_munged_dir, "SCZ.sumstats.gz")

cog_files <- c(
  "Matrix_extracted_full.txt",
  "Memory_extracted_full.txt",
  "RT_extracted_full.txt",
  "SymbolDigit_extracted_full.txt",
  "TMTB_extracted_full.txt",
  "Tower_extracted_full.txt",
  "VNR_extracted_full.txt"
)

cog_names <- c("Matrix", "Memory", "RT", "SymbolDigit", "TMTB", "Tower", "VNR")
cog_Ns    <- c(11356, 331679, 330024, 87741, 78547, 11263, 171304)

local_cog_files <- file.path(output_dir, cog_files)
source_cog_files <- file.path(base_dir, cog_files)
# We don't want to copy the old broken munged files anymore, so clear the shortcut
source_cog_munged_files <- file.path(cognition_munged_dir, paste0(cog_names, "_INVALID.sumstats.gz"))

all_names <- c("SCZ", cog_names)
all_Ns    <- c(501319, cog_Ns)

local_munged_files <- file.path(output_dir, paste0(all_names, ".sumstats.gz"))
source_munged_files <- c(scz_munged_file, source_cog_munged_files)

if (all(file.exists(source_munged_files))) {
  cat("=== munged sumstats ===\n")
  for (idx in seq_along(source_munged_files)) {
    src <- source_munged_files[idx]
    dst <- local_munged_files[idx]
    ok <- file.copy(src, dst, overwrite = TRUE)
    if (!ok) {
      stop("Failed to copy munged file: ", src)
    }
    cat("  copy:", basename(dst), "\n")
  }
} else {
  cat("=== munged error， raw munge ===\n")
  for (idx in seq_along(source_cog_files)) {
    src <- source_cog_files[idx]
    dst <- local_cog_files[idx]
    dt <- fread(src)
    setnames(dt, old = c("beta", "se", "p"), new = c("BETA", "SE", "P"), skip_absent = TRUE)
    fwrite(dt, dst, sep = "\t", na = "NA", quote = FALSE)
    cat("  已写入:", basename(dst), "\n")
  }

  cat("\n=== munge： sumstats ===\n")
  munge(
    files = c(scz_raw_file, local_cog_files),
    hm3 = hm3,
    trait.names = all_names,
    N = all_Ns,
    info.filter = 0,
    maf.filter = 0
  )
}

cat("\n=== LDSC ===\n")
munged_files <- paste0(all_names, ".sumstats.gz")

ldsc_out <- ldsc(
  traits = munged_files,
  sample.prev = rep(NA, length(all_names)),
  population.prev = rep(NA, length(all_names)),
  ld = ld_dir,
  wld = ld_dir,
  trait.names = all_names,
  ldsc.log = "ldsc_schizophrenia_cognition"
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

scz_idx <- which(all_names == "SCZ")
cog_idx <- which(all_names != "SCZ")

rg_vec <- Rg[scz_idx, cog_idx]
se_vec <- se_rg[scz_idx, cog_idx]
z_vec <- rg_vec / se_vec
p_vec <- 2 * pnorm(-abs(z_vec))
h2_scz <- S[scz_idx, scz_idx]
h2_cog <- diag(S)[cog_idx]
intercept_vec <- I[scz_idx, cog_idx]

results <- data.frame(
  Exposure = "SCZ",
  Outcome = all_names[cog_idx],
  rg = round(rg_vec, 4),
  SE = round(se_vec, 4),
  Z = round(z_vec, 3),
  P = signif(p_vec, 4),
  h2_SCZ = round(h2_scz, 4),
  h2_Cognition = round(h2_cog, 4),
  CrossIntercept = round(intercept_vec, 4),
  row.names = NULL
)

write.csv(results, "ldsc_schizophrenia_cognition_rg.csv", row.names = FALSE)
write.csv(as.data.frame(Rg), "ldsc_schizophrenia_cognition_Rg_full_matrix.csv")
write.csv(as.data.frame(S), "ldsc_schizophrenia_cognition_S_covariance.csv")
write.csv(as.data.frame(I), "ldsc_schizophrenia_cognition_intercept.csv")

