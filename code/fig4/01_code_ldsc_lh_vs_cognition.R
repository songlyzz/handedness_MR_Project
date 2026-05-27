library(GenomicSEM)
library(data.table)

base_dir <- "./UKBbehavior"
ld_dir <- "./eur_w_ld_chr"
hm3 <- "./w_hm3.snplist.gz"

setwd(base_dir)

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


cog_Ns <- c(11356, 331679, 330024, 87741, 78547, 11263, 171304)

all_files <- c("lefthands_tb.txt", cog_files)
all_names <- c("Handedness",       cog_names)
all_Ns    <- c(580148,             cog_Ns)

for (f in cog_files) {
  dt <- fread(f)
  setnames(dt, old = c("beta", "se", "p"), new = c("BETA", "SE", "P"), skip_absent = TRUE)
  fwrite(dt, f, sep = "\t", na = "NA", quote = FALSE)
  cat("  deal:", f, "\n")
}

munge(
  files       = all_files,
  hm3         = hm3,
  trait.names = all_names,
  N           = all_Ns,
  info.filter = 0,
  maf.filter  = 0
)

munged_files <- paste0(all_names, ".sumstats.gz")

ldsc_out <- ldsc(
  traits          = munged_files,
  sample.prev     = rep(NA, length(all_names)),
  population.prev = rep(NA, length(all_names)),
  ld              = ld_dir,
  wld             = ld_dir,
  trait.names     = all_names,
  ldsc.log        = "ldsc_handedness_cognition"
)

