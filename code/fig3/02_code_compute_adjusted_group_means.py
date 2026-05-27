from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "analysis_variables_final.csv"
OUT_FILE = BASE_DIR / "group_means_adjusted.csv"

TASKS = [
    ("RT", "RT_I0", "age_I0", "centre_I0"),
    ("PM_err_a1", "PM_err_R1_I0_a1", "age_I0", "centre_I0"),
    ("PM_err_a2", "PM_err_R1_I0_a2", "age_I0", "centre_I0"),
    ("PM_time_a1", "PM_time_R1_I0_a1", "age_I0", "centre_I0"),
    ("PM_time_a2", "PM_time_R1_I0_a2", "age_I0", "centre_I0"),
    ("FI", "FI", "FI_age", "FI_centre"),
    ("ProspMem", "ProspMem_firsttry", "ProspMem_age", "ProspMem_centre"),
    ("DigitSpan", "DigitSpan", "DigitSpan_age", "DigitSpan_centre"),
    ("TMT_A_dur", "TMT_A_dur", "TMT_A_age", "TMT_A_centre"),
    ("TMT_A_err", "TMT_A_err", "TMT_A_age", "TMT_A_centre"),
    ("TMT_B_dur", "TMT_B_dur", "TMT_B_age", "TMT_B_centre"),
    ("TMT_B_err", "TMT_B_err", "TMT_B_age", "TMT_B_centre"),
    ("Matrix_correct", "Matrix_correct", "Matrix_age", "Matrix_centre"),
    ("Matrix_timeper", "Matrix_timeper", "Matrix_age", "Matrix_centre"),
    ("PAL_correct", "PAL_correct", "PAL_age", "PAL_centre"),
    ("Tower_correct", "Tower_correct", "Tower_age", "Tower_centre"),
    ("SymDig_correct", "SymDig_correct", "SymDig_age", "SymDig_centre"),
    ("BookLetter", "BookLetter", "Book_age", "Book_centre"),
    ("Vocab_level", "Vocab_level", "Vocab_age", "Vocab_centre"),
    ("Vocab_SCA", "Vocab_SCA", "Vocab_age", "Vocab_centre"),
]

PC_COLS = [f"PC{i}" for i in range(1, 11)]

NEGATIVE_DIRECTION = {
    "RT_I0",
    "PM_err_R1_I0_a1", "PM_err_R1_I0_a2",
    "PM_time_R1_I0_a1", "PM_time_R1_I0_a2",
    "TMT_A_dur", "TMT_A_err",
    "TMT_B_dur", "TMT_B_err",
    "Matrix_timeper",
}


def to_cat_str(series):
    return series.apply(
        lambda x: str(int(float(x)))
        if pd.notna(x) else np.nan
    )


# Load data
df = pd.read_csv(DATA_FILE, low_memory=False)

df = df[df["handedness_I0"].isin([1.0, 2.0])].copy()
df["left_handed"] = (df["handedness_I0"] == 2.0).astype(int)

# Recode variables
if "PAL_correct" not in df.columns and "PAL_errors" in df.columns:
    df["PAL_correct"] = pd.to_numeric(df["PAL_errors"], errors="coerce")

pm_raw = pd.to_numeric(df.get("ProspMem"), errors="coerce")
df["ProspMem_firsttry"] = np.where(
    pm_raw.isna(),
    np.nan,
    (pm_raw == 1.0).astype(float)
)

# Convert categorical variables
CAT_COLS = [
    "sex",
    "centre_I0", "centre_I2",
    "FI_centre", "ProspMem_centre", "DigitSpan_centre",
    "TMT_A_centre", "TMT_B_centre", "Matrix_centre",
    "PAL_centre", "Tower_centre", "SymDig_centre",
    "Book_centre", "Vocab_centre",
]

for c in CAT_COLS:
    if c in df.columns:
        df[c] = to_cat_str(df[c])

rows = []

for task_label, out_col, age_col, centre_col in TASKS:

    if out_col not in df.columns:
        continue

    out_safe = out_col.replace("-", "_").replace(".", "_")
    age_safe = age_col.replace("-", "_").replace(".", "_")
    ctr_safe = centre_col.replace("-", "_").replace(".", "_")

    needed = (
        [out_col, age_col, centre_col, "left_handed", "sex"]
        + [p for p in PC_COLS if p in df.columns]
    )

    sub = df[[c for c in needed if c in df.columns]].dropna().copy()

    n_right = int((sub["left_handed"] == 0).sum())
    n_left = int((sub["left_handed"] == 1).sum())

    if n_left < 20 or n_right < 20:
        continue

    sub = sub.rename(columns={
        out_col: out_safe,
        age_col: age_safe,
        centre_col: ctr_safe
    })

    # Covariate-only model
    covariates = (
        [age_safe, "sex", f"C({ctr_safe})"]
        + [p for p in PC_COLS if p in sub.columns]
    )

    formula = f"{out_safe} ~ " + " + ".join(covariates)

    try:
        fit = smf.ols(formula, data=sub).fit()
    except Exception:
        continue

    # Adjusted scores
    grand_mean = float(sub[out_safe].mean())
    adjusted = fit.resid.values + grand_mean

    # Standardization
    z = (adjusted - adjusted.mean()) / adjusted.std(ddof=1)

    # Reverse negative-direction measures
    if out_col in NEGATIVE_DIRECTION:
        z = -z

    mask_r = sub["left_handed"].values == 0
    mask_l = sub["left_handed"].values == 1

    z_r = z[mask_r]
    z_l = z[mask_l]

    se_r = z_r.std(ddof=1) / np.sqrt(len(z_r))
    se_l = z_l.std(ddof=1) / np.sqrt(len(z_l))

    diff = float(z_l.mean() - z_r.mean())
    se_diff = np.sqrt(se_r**2 + se_l**2)

    rows.append({
        "task": task_label,
        "n_right": n_right,
        "n_left": n_left,
        "z_mean_right": round(float(z_r.mean()), 5),
        "z_mean_left": round(float(z_l.mean()), 5),
        "effect_size_d": round(diff, 5),
        "ci_lower": round(diff - 1.96 * se_diff, 5),
        "ci_upper": round(diff + 1.96 * se_diff, 5),
    })

# Save results
out_df = pd.DataFrame(rows)
out_df.to_csv(OUT_FILE, index=False)

print(out_df)
print(f"\nSaved to: {OUT_FILE}")