from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parent

MEANS_FILE = BASE / "group_means_adjusted.csv"
RESULT_FILE = BASE / "mixed_model_results_M2.csv"

# Tasks included in FDR correction
SELECTED = [
    "RT",
    "PM_err_a1", "PM_err_a2",
    "FI", "ProspMem", "DigitSpan",
    "TMT_A_dur", "TMT_B_dur", "TMT_BA_dur",
    "Matrix_correct", "Matrix_timeper",
    "PAL_correct", "Tower_correct", "SymDig_correct",
    "BookLetter", "Vocab_SCA",
]

# Load data
df = pd.read_csv(MEANS_FILE)
res = pd.read_csv(RESULT_FILE)

df["task"] = df["task"].replace({"PAL_errors": "PAL_correct"})
res["task"] = res["task"].replace({"PAL_errors": "PAL_correct"})

# Recalculate BH-FDR
res_sel = res[res["task"].isin(SELECTED)].copy()

_, p_fdr, _, _ = multipletests(res_sel["p_raw"], method="fdr_bh")
res_sel["p_fdr"] = p_fdr

def sig_label(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""

res_sel["sig_fdr"] = res_sel["p_fdr"].apply(sig_label)

sig_map = res_sel.set_index("task")["sig_fdr"].to_dict()
df["sig_fdr"] = df["task"].map(sig_map).fillna("")

# Labels
LABELS = {
    "RT": "Reaction Time",
    "PM_err_a1": "PM Error (1)",
    "PM_err_a2": "PM Error (2)",
    "FI": "Fluid Intelligence",
    "ProspMem": "Prospective Memory",
    "DigitSpan": "Digit Span",
    "TMT_A_dur": "TMT-A",
    "TMT_B_dur": "TMT-B",
    "TMT_BA_dur": "TMT B-A",
    "Matrix_correct": "Matrix Accuracy",
    "Matrix_timeper": "Matrix Time",
    "PAL_correct": "PAL",
    "Tower_correct": "Tower",
    "SymDig_correct": "Symbol Digit",
    "BookLetter": "Letter Reading",
    "Vocab_SCA": "Vocabulary",
}

GROUPS = [
    ("Processing Speed / Memory",
     ["RT", "PM_err_a1", "PM_err_a2"]),

    ("General Ability",
     ["FI", "ProspMem", "DigitSpan"]),

    ("Executive Function",
     ["TMT_A_dur", "TMT_B_dur", "TMT_BA_dur"]),

    ("Reasoning / Learning",
     ["Matrix_correct", "Matrix_timeper",
      "PAL_correct", "Tower_correct",
      "SymDig_correct"]),

    ("Language",
     ["BookLetter", "Vocab_SCA"]),
]

# Build plotting order
plot_items = []
y = 0

for group_name, tasks in reversed(GROUPS):

    plot_items.append((y, None, group_name, True))
    y += 0.8

    for task in tasks:
        if task not in df["task"].values:
            continue

        plot_items.append(
            (y, task, LABELS.get(task, task), False)
        )
        y += 1

    y += 0.4

# Colors
C_SIG_POS = "#1f77b4"
C_SIG_NEG = "#d62728"
C_NS = "#888888"

# Marker types
METHOD_MARKER = {
    "OLS-residual-d": "s",
    "logit-d": "^",
    "tobit-d": "D",
}

# Figure
fig_h = max(7, len(plot_items) * 0.45)

fig, ax = plt.subplots(figsize=(9, fig_h))

ax.axvline(0, linestyle="--", linewidth=1, color="black", alpha=0.7)

y_ticks = []
y_labels = []

# Draw points
for y_pos, task, label, is_header in plot_items:

    if is_header:
        y_ticks.append(y_pos)
        y_labels.append(label)
        continue

    row = df[df["task"] == task].iloc[0]

    d = float(row["cohen_d"])
    lo = float(row["cohen_d_ci_lo"])
    hi = float(row["cohen_d_ci_hi"])

    sig = str(row["sig_fdr"])
    method = str(row["effect_method"])

    if sig:
        color = C_SIG_POS if d > 0 else C_SIG_NEG
        lw = 2
        ms = 7
    else:
        color = C_NS
        lw = 1.2
        ms = 5

    marker = METHOD_MARKER.get(method, "o")

    # CI line
    ax.plot(
        [lo, hi],
        [y_pos, y_pos],
        color=color,
        linewidth=lw,
        zorder=2
    )

    # Point estimate
    ax.plot(
        d,
        y_pos,
        marker=marker,
        color=color,
        markersize=ms,
        zorder=3
    )

    y_ticks.append(y_pos)
    y_labels.append(label)

# Axis settings
ax.set_yticks(y_ticks)

tick_objs = ax.set_yticklabels(y_labels, fontsize=8)

headers = {
    label for (_, _, label, is_h) in plot_items if is_h
}

for t in tick_objs:
    if t.get_text() in headers:
        t.set_fontweight("bold")

ax.set_xlabel("Cohen's d (left − right)")
ax.set_title(
    "Handedness and Cognitive Performance",
    fontsize=11,
    fontweight="bold"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Legend
legend_elements = [
    Line2D(
        [0], [0],
        marker="s",
        color="w",
        markerfacecolor=C_SIG_POS,
        markersize=8,
        label="Significant left-hander advantage"
    ),
    Line2D(
        [0], [0],
        marker="s",
        color="w",
        markerfacecolor=C_SIG_NEG,
        markersize=8,
        label="Significant right-hander advantage"
    ),
    Line2D(
        [0], [0],
        marker="s",
        color="w",
        markerfacecolor=C_NS,
        markersize=8,
        label="Non-significant"
    ),
]

ax.legend(
    handles=legend_elements,
    loc="lower right",
    fontsize=7
)

plt.tight_layout()

# Save figure
for ext in ["png", "pdf", "svg"]:
    fig.savefig(
        BASE / f"forest_cohend.{ext}",
        dpi=300,
        bbox_inches="tight"
    )

plt.show()