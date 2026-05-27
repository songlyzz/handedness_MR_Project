# =============================================================================
# ldsc_handedness_disease.R
# GenomicSEM 
# =============================================================================

library(GenomicSEM)
library(data.table)

# ── Input files ──────────────────────────────────────────────────────────────
traits <- c(
  "Handedness.sumstats.gz",
  "ADHD.sumstats.gz",
  "PTSD.sumstats.gz",
  "SCZ.sumstats.gz"
)

trait_names <- c(
  "Handedness",
  "ADHD",
  "PTSD",
  "SCZ"
)

# LD score reference files
ld_dir <- "eur_w_ld_chr"

# ── Run LDSC ────────────────────────────────────────────────────────────────
ldsc_out <- ldsc(
  traits          = traits,
  sample.prev     = rep(NA, length(trait_names)),
  population.prev = rep(NA, length(trait_names)),
  ld              = ld_dir,
  wld             = ld_dir,
  trait.names     = trait_names
)

# ── Extract matrices ─────────────────────────────────────────────────────────
S <- ldsc_out$S
V <- ldsc_out$V
I <- ldsc_out$I

colnames(S) <- rownames(S) <- all_names
colnames(I) <- rownames(I) <- all_names

Rg <- cov2cor(S)
colnames(Rg) <- rownames(Rg) <- all_names

vech_idx <- function(i, j, n) {
  if (i < j) { tmp <- i; i <- j; j <- tmp }
  (j - 1) * n - j * (j - 1) / 2 + i
}

n     <- length(all_names)
se_rg <- matrix(NA_real_, n, n)
colnames(se_rg) <- rownames(se_rg) <- all_names

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) next
    idx_ij <- vech_idx(i, j, n)
    se_rg[i, j] <- sqrt(V[idx_ij, idx_ij]) / sqrt(abs(S[i, i] * S[j, j]))
  }
}

# ── Handedness vs SD ────────────────────────────────────────────────
hand_idx <- which(all_names == "Handedness")
dis_idx  <- which(all_names != "Handedness")

rg_vec        <- Rg[hand_idx, dis_idx]
se_vec        <- se_rg[hand_idx, dis_idx]
z_vec         <- rg_vec / se_vec
p_vec         <- 2 * pnorm(-abs(z_vec))
h2_hand       <- S[hand_idx, hand_idx]
h2_dis        <- diag(S)[dis_idx]
intercept_vec <- I[hand_idx, dis_idx]

results <- data.frame(
  Exposure       = "Handedness",
  Outcome        = all_names[dis_idx],
  rg             = round(rg_vec, 4),
  SE             = round(se_vec, 4),
  Z              = round(z_vec,  3),
  P              = signif(p_vec, 4),
  h2_Handedness  = round(h2_hand, 4),
  h2_Disease     = round(h2_dis,  4),
  CrossIntercept = round(intercept_vec, 4),
  row.names      = NULL
)

print(results, digits = 4)

# ─────────────────────────────────────────────────────────────────────
write.csv(results,           "ldsc_handedness_disease_rg.csv",    row.names = FALSE)
write.csv(as.data.frame(Rg), "ldsc_disease_Rg_full_matrix.csv")
write.csv(as.data.frame(S),  "ldsc_disease_S_covariance.csv")
write.csv(as.data.frame(I),  "ldsc_disease_intercept.csv")
