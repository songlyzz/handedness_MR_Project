"""
run_ols_analysis.py
────────────────────────────────────────────────────────────────────────
OLS regression

test：
      M1: outcome ~ left_handed + age + sex + C(centre)
      M2: outcome ~ left_handed + age + sex + C(centre) + C(geno_batch) + PC1…PC10

outcome：
  ols_results.csv        
  ols_results_M1.csv     — Model 1 
  ols_results_M2.csv     — Model 2 
"""

from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# BASE_DIR  = Path(r"")
DATA_FILE = BASE_DIR / "analysis_variables_final.csv"
OUT_ALL   = BASE_DIR / "ols_results.csv"
OUT_M1    = BASE_DIR / "ols_results_M1.csv"
OUT_M2    = BASE_DIR / "ols_results_M2.csv"

# ── ──────────────────────────────────────────────────────────
TASKS = [
    ("RT",             "RT_I0",            "age_I0",        "centre_I0"),
    ("PM_err_a1",      "PM_err_R1_I0_a1",  "age_I0",        "centre_I0"),
    ("PM_err_a2",      "PM_err_R1_I0_a2",  "age_I0",        "centre_I0"),
    ("PM_time_a1",     "PM_time_R1_I0_a1", "age_I0",        "centre_I0"),
    ("PM_time_a2",     "PM_time_R1_I0_a2", "age_I0",        "centre_I0"),
    ("FI",             "FI",               "FI_age",         "FI_centre"),
    ("ProspMem",       "ProspMem_firsttry", "ProspMem_age",   "ProspMem_centre"),
    ("DigitSpan",      "DigitSpan",        "DigitSpan_age",  "DigitSpan_centre"),
    ("TMT_A_dur",      "TMT_A_dur",        "TMT_A_age",      "TMT_A_centre"),
    ("TMT_A_err",      "TMT_A_err",        "TMT_A_age",      "TMT_A_centre"),
    ("TMT_B_dur",      "TMT_B_dur",        "TMT_B_age",      "TMT_B_centre"),
    ("TMT_B_err",      "TMT_B_err",        "TMT_B_age",      "TMT_B_centre"),
    ("Matrix_correct", "Matrix_correct",   "Matrix_age",     "Matrix_centre"),
    ("Matrix_timeper", "Matrix_timeper",   "Matrix_age",     "Matrix_centre"),
    ("PAL_correct",    "PAL_correct",      "PAL_age",        "PAL_centre"),
    ("Tower_correct",  "Tower_correct",    "Tower_age",      "Tower_centre"),
    ("SymDig_correct", "SymDig_correct",   "SymDig_age",     "SymDig_centre"),
    ("BookLetter",     "BookLetter",       "Book_age",       "Book_centre"),
    ("Vocab_level",    "Vocab_level",      "Vocab_age",      "Vocab_centre"),
    ("Vocab_SCA",      "Vocab_SCA",        "Vocab_age",      "Vocab_centre"),
]

PC_COLS = [f"PC{i}" for i in range(1, 11)]


def to_cat_str(series: pd.Series) -> pd.Series:
    def _conv(x):
        if pd.isna(x):
            return np.nan
        try:
            return str(int(float(x)))
        except (ValueError, TypeError):
            return np.nan
    return series.apply(_conv)


# ── 1. ───────────────────────────────────────────────────────

df = pd.read_csv(DATA_FILE, low_memory=False)

df = df[df["handedness_I0"].isin([1.0, 2.0])].copy()
n_right = int((df["handedness_I0"] == 1.0).sum())
n_left  = int((df["handedness_I0"] == 2.0).sum())


df["left_handed"] = (df["handedness_I0"] == 2.0).astype(int)
if "PAL_correct" not in df.columns and "PAL_errors" in df.columns:
    df["PAL_correct"] = pd.to_numeric(df["PAL_errors"], errors="coerce")
pm_raw = pd.to_numeric(df.get("ProspMem"), errors="coerce")
df["ProspMem_firsttry"] = np.where(pm_raw.isna(), np.nan, (pm_raw == 1.0).astype(float))

df["sex"] = to_cat_str(df["sex"])

CAT_COLS = [
    "centre_I0", "centre_I2",
    "FI_centre", "ProspMem_centre", "DigitSpan_centre",
    "TMT_A_centre", "TMT_B_centre", "Matrix_centre",
    "PAL_centre", "Tower_centre", "SymDig_centre",
    "Book_centre", "Vocab_centre",
    "geno_batch",
]
for c in CAT_COLS:
    if c in df.columns:
        df[c] = to_cat_str(df[c])

# ──────────────────────────────────────────────────────

results: list[dict] = []

for task_label, out_col, age_col, centre_col in TASKS:
    if out_col not in df.columns:
        print(f"  [{task_label}] pass")
        continue
    if centre_col not in df.columns:
        print(f"  [{task_label}] centre pass")
        continue

    print(f"  [{task_label}]", end="", flush=True)

    out_safe = out_col.replace("-", "_").replace(".", "_")
    age_safe = age_col.replace("-", "_").replace(".", "_")
    ctr_safe = centre_col.replace("-", "_").replace(".", "_")

    for model_id in (1, 2):
        base_need  = [out_col, age_col, centre_col, "left_handed", "sex"]
        geno_need  = PC_COLS
        needed_raw = base_need + (geno_need if model_id == 2 else [])

        needed_exist = [c for c in needed_raw if c in df.columns]
        sub = df[needed_exist].copy().dropna()

        n_total = len(sub)
        n_r = int((sub["left_handed"] == 0).sum())
        n_l = int((sub["left_handed"] == 1).sum())

        if n_l < 20 or n_r < 20:
            print(f" M{model_id}", end="")
            continue

        rename_map = {out_col: out_safe, age_col: age_safe, centre_col: ctr_safe}
        sub = sub.rename(columns=rename_map)

        covariates = [age_safe, "sex", f"C({ctr_safe})"]
        if model_id == 2:
            covariates += [p for p in PC_COLS if p in sub.columns]

        formula = f"{out_safe} ~ left_handed + " + " + ".join(covariates)

        try:
            fit = smf.ols(formula, data=sub).fit()
        except Exception as e:
            print(f" M{model_id}[erorr:{e}]", end="")
            continue

        beta  = fit.params["left_handed"]
        se    = fit.bse["left_handed"]
        tstat = fit.tvalues["left_handed"]
        pval  = fit.pvalues["left_handed"]
        ci    = fit.conf_int().loc["left_handed"]
        r2    = fit.rsquared

        sd_x = float(sub["left_handed"].std(ddof=1))
        sd_y = float(pd.to_numeric(sub[out_safe], errors="coerce").std(ddof=1))
        std_beta = beta * sd_x / sd_y if sd_y > 0 else np.nan

        results.append({
            "task":        task_label,
            "outcome_col": out_col,
            "model":       f"M{model_id}",
            "n_total":     n_total,
            "n_right":     n_r,
            "n_left":      n_l,
            "beta":        round(beta,    6),
            "se":          round(se,      6),
            "t":           round(tstat,   4),
            "p_raw":       pval,
            "ci_lo":       round(float(ci.iloc[0]), 6),
            "ci_hi":       round(float(ci.iloc[1]), 6),
            "std_beta":    round(std_beta, 6),
            "r2":          round(r2,       6),
        })
        print(f" M{model_id}✓", end="", flush=True)

    print()

# ─────────────────────────────────────────────────────────
print("\n[3/4] Benjamini-Hochberg FDR ...")
res_df = pd.DataFrame(results)

for m in ("M1", "M2"):
    mask = res_df["model"] == m
    if mask.sum() < 2:
        res_df.loc[mask, "p_fdr"] = res_df.loc[mask, "p_raw"]
        continue
    _, p_fdr, _, _ = multipletests(res_df.loc[mask, "p_raw"], method="fdr_bh")
    res_df.loc[mask, "p_fdr"] = p_fdr
    print(f"  {m}: {mask.sum()} ")

res_df["p_raw"] = res_df["p_raw"].round(8)
res_df["p_fdr"] = res_df["p_fdr"].round(8)

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""
res_df["sig_fdr"] = res_df["p_fdr"].apply(sig_label)

# ─────────────────────────────────────────────────────────
col_order = ["task","outcome_col","model",
             "n_total","n_right","n_left",
             "beta","se","t","p_raw","p_fdr","sig_fdr",
             "ci_lo","ci_hi","std_beta","r2"]
res_df = res_df[[c for c in col_order if c in res_df.columns]]

res_df.to_csv(OUT_ALL, index=False, encoding="utf-8-sig")
res_df[res_df["model"]=="M1"].to_csv(OUT_M1, index=False, encoding="utf-8-sig")
res_df[res_df["model"]=="M2"].to_csv(OUT_M2, index=False, encoding="utf-8-sig")
print(f"  → {OUT_ALL}")
print(f"  → {OUT_M1}")
print(f"  → {OUT_M2}")

# ────────────────────────────────────────────────────────
for model_name in ("M1", "M2"):
    label = "model 1" if model_name == "M1" else "model 2"
    sub = res_df[res_df["model"] == model_name]
    print(f"\n{'='*100}")
    print(f"  {label}  ({model_name})")
    print(f"  cog ~ left_handed + age + sex + C(centre)" +
          ("" if model_name == "M1" else " + C(geno_batch) + PC1..10"))
    print(f"{'='*100}")
    hdr = f"  {'task':<16} {'Nt':>7} {'NR':>7} {'NL':>7} {'Beta':>9} {'SE':>7} "
    hdr += f"{'t':>7} {'p_raw':>10} {'p_fdr':>10} {'std_β':>8}  sig"
    print(hdr)
    print("  " + "-"*96)
    for _, r in sub.iterrows():
        print(
            f"  {r['task']:<16} {r['n_total']:>7,} {r['n_right']:>7,} {r['n_left']:>7,} "
            f"{r['beta']:>9.4f} {r['se']:>7.4f} {r['t']:>7.3f} "
            f"{r['p_raw']:>10.2e} {r['p_fdr']:>10.2e} {r['std_beta']:>8.4f}  {r['sig_fdr']}"
        )

