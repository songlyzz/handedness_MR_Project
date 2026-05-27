#!/usr/bin/env Rscript
# extract_g_factor.R
# Single-factor CFA for general cognitive ability (g)

pkgs <- c("lavaan", "dplyr", "readr", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(lavaan)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ── Paths ─────────────────────────────────────────────
BASE_DIR <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
DATA_FILE <- file.path(BASE_DIR, "analysis_variables_final.csv")

df <- read_csv(DATA_FILE, show_col_types = FALSE) %>%
  filter(handedness_I0 %in% c(1, 2)) %>%
  mutate(
    left_handed = as.integer(handedness_I0 == 2),
    sex_c = as.character(as.integer(sex))
  )

pc_cols <- paste0("PC", 1:10)
pc_cols <- pc_cols[pc_cols %in% names(df)]

# ── Tasks ─────────────────────────────────────────────
TASKS <- list(
  list(label="RT",     col="RT_I0",           age="age_I0",      ctr="centre_I0",     dir=-1),
  list(label="PM_err", col="PM_err_R1_I0_a2", age="age_I0",      ctr="centre_I0",     dir=-1),
  list(label="FI",     col="FI",              age="FI_age",      ctr="FI_centre",     dir= 1),
  list(label="SymDig", col="SymDig_correct",  age="SymDig_age",  ctr="SymDig_centre", dir= 1),
  list(label="Tower",  col="Tower_correct",   age="Tower_age",   ctr="Tower_centre",  dir= 1),
  list(label="Matrix", col="Matrix_correct",  age="Matrix_age",  ctr="Matrix_centre", dir= 1),
  list(label="TMT_B",  col="TMT_B_dur",       age="TMT_B_age",   ctr="TMT_B_centre",  dir=-1)
)

# ── Residualize tasks ─────────────────────────────────
resid_list <- list()

for (tk in TASKS) {
  
  sub <- df %>%
    mutate(
      .row = row_number(),
      y    = as.numeric(.data[[tk$col]]),
      age  = as.numeric(.data[[tk$age]]),
      ctr  = as.character(as.integer(as.numeric(.data[[tk$ctr]])))
    ) %>%
    select(.row, y, age, ctr, sex_c, all_of(pc_cols)) %>%
    drop_na()
  
  formula_str <- paste0(
    "y ~ age + sex_c + factor(ctr)",
    if (length(pc_cols) > 0) paste0(" + ", paste(pc_cols, collapse=" + ")) else ""
  )
  
  fit <- lm(as.formula(formula_str), data = sub)
  
  resid_z <- scale(residuals(fit))[,1] * tk$dir
  
  resid_list[[tk$label]] <- setNames(resid_z, sub$.row)
}

# ── Common sample ─────────────────────────────────────
common_rows <- Reduce(intersect, lapply(resid_list, function(x) as.integer(names(x))))

g_data <- data.frame(row_id = common_rows)

for (nm in names(resid_list)) {
  g_data[[nm]] <- as.numeric(resid_list[[nm]][as.character(common_rows)])
}

g_data <- g_data %>%
  left_join(
    df %>% mutate(.row = row_number()) %>%
      select(.row, left_handed),
    by = c("row_id" = ".row")
  )

# ── CFA ───────────────────────────────────────────────
indicators <- names(resid_list)

model <- paste0(
  "g =~ ",
  paste(indicators, collapse = " + ")
)

fit <- cfa(
  model,
  data      = g_data[, indicators],
  estimator = "MLR",
  std.lv    = TRUE
)

# ── Loadings ──────────────────────────────────────────
params <- parameterEstimates(fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  transmute(
    indicator   = rhs,
    loading_std = std.all,
    p_value     = pvalue
  )

# ── Fit indices ───────────────────────────────────────
fit_df <- tibble(
  metric = c("CFI", "TLI", "RMSEA", "SRMR"),
  value  = c(
    fitMeasures(fit, "cfi.robust"),
    fitMeasures(fit, "tli.robust"),
    fitMeasures(fit, "rmsea.robust"),
    fitMeasures(fit, "srmr")
  )
)

# ── Factor scores ─────────────────────────────────────
g_data$g_score <- lavPredict(fit)[, "g"]

# ── Handedness effect ─────────────────────────────────
lm_g <- lm(g_score ~ left_handed, data = g_data)

coef_g <- summary(lm_g)$coefficients["left_handed", ]
ci_g   <- confint(lm_g)["left_handed", ]

result <- tibble(
  beta  = coef_g["Estimate"],
  se    = coef_g["Std. Error"],
  t     = coef_g["t value"],
  p     = coef_g["Pr(>|t|)"],
  ci_lo = ci_g[1],
  ci_hi = ci_g[2]
)

# ── Save ──────────────────────────────────────────────
write_csv(params,  file.path(BASE_DIR, "g_factor_loadings.csv"))
write_csv(fit_df, file.path(BASE_DIR, "g_factor_fit.csv"))
write_csv(result, file.path(BASE_DIR, "g_handedness_result.csv"))
write_csv(
  g_data %>% select(row_id, left_handed, g_score),
  file.path(BASE_DIR, "g_factor_scores.csv")
)

# ── Plot ──────────────────────────────────────────────
p <- ggplot(params, aes(x = loading_std, y = reorder(indicator, loading_std))) +
  geom_point(size = 3, color = "#2c6e9d") +
  theme_classic() +
  labs(
    x = "Standardized Loading",
    y = NULL,
    title = "General Cognitive Factor (g)"
  )

ggsave(file.path(BASE_DIR, "g_factor_loadings.png"), p, width = 6, height = 4)
ggsave(file.path(BASE_DIR, "g_factor_loadings.pdf"), p, width = 6, height = 4)

cat("Done.\n")