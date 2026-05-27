"""
plot_forest_gfactor.py  (v2 — g-factor structure)
────────────────────────────────────────────────────────────────────────
Forest plot with three-tier layout:

  上半部分: 7 个 g-factor 指标任务          （浅蓝背景）
  中间行:   General Factor g (CFA score)   （紫色菱形，特殊标注）
  下半部分: 其余 7 个认知任务（无强制分组）
            ✗ 已排除：TMT_BA_dur、Matrix_timeper

  FDR 矫正: 全部 15 个任务统一 BH-FDR

副作用:
  • 将 g_factor 行（cohen_d / CI / N）合并到 group_means_adjusted.csv
  • 将 g_factor 行（p_raw 等）合并到 mixed_model_results_M2.csv
  • 输出: forest_gfactor.{svg / pdf / png}
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as scipy_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# ── 路径 ────────────────────────────────────────────────────────────────
BASE         = Path(__file__).resolve().parent
MEANS_CSV    = BASE / "group_means_adjusted.csv"
M2_CSV       = BASE / "mixed_model_results_M2.csv"
G_SCORES_CSV = BASE / "g_factor_scores.csv"

# ── 任务分组 ─────────────────────────────────────────────────────────────
G_TASKS = [
    "RT", "PM_err_a2", "FI",
    "SymDig_correct", "Tower_correct", "Matrix_correct", "TMT_B_dur",
]
OTHER_TASKS = [
    "PM_err_a1", "ProspMem", "DigitSpan",
    "TMT_A_dur", "PAL_correct",
    "BookLetter", "Vocab_SCA",
]
ALL_TASKS = G_TASKS + ["g_factor"] + OTHER_TASKS   # 15 个

# ── 显示标签 ─────────────────────────────────────────────────────────────
LABELS = {
    "RT":             "Reaction Time",
    "PM_err_a2":      "Pairs Matching — Errors (attempt 2)",
    "FI":             "Fluid Intelligence",
    "SymDig_correct": "Symbol-Digit Coding",
    "Tower_correct":  "Tower Rearranging",
    "Matrix_correct": "Matrix Reasoning (accuracy)",
    "TMT_B_dur":      "Trail Making B — Duration",
    "g_factor":       "General Factor  g  (CFA, 7 indicators)",
    "PM_err_a1":      "Pairs Matching — Errors (attempt 1)",
    "ProspMem":       "Prospective Memory †",
    "DigitSpan":      "Digit Span",
    "TMT_A_dur":      "Trail Making A — Duration",
    "PAL_correct":    "Paired Associates Learning — Correct",
    "BookLetter":     "Letter Reading (NART proxy)",
    "Vocab_SCA":      "Vocabulary (SCA, UKB recommended)",
}

# ════════════════════════════════════════════════════════════════════════
# 1. 计算 g_factor Cohen's d（来自 CFA factor scores，已残差化）
# ════════════════════════════════════════════════════════════════════════
gs   = pd.read_csv(G_SCORES_CSV)
gs_r = gs[gs["left_handed"] == 0]["g_score"].dropna().values
gs_l = gs[gs["left_handed"] == 1]["g_score"].dropna().values

nr, nl   = len(gs_r), len(gs_l)
mr, ml   = gs_r.mean(), gs_l.mean()
sr, sl   = gs_r.std(ddof=1), gs_l.std(ddof=1)
s_pool   = np.sqrt(((nr - 1)*sr**2 + (nl - 1)*sl**2) / (nr + nl - 2))
d_g      = (ml - mr) / s_pool
se_d_g   = np.sqrt((nr + nl) / (nr * nl) + d_g**2 / (2*(nr + nl)))
d_g_lo   = d_g - 1.96 * se_d_g
d_g_hi   = d_g + 1.96 * se_d_g
t_g, p_g = scipy_stats.ttest_ind(gs_l, gs_r)

print(f"\ng_factor: d={d_g:+.4f}  95% CI [{d_g_lo:+.4f}, {d_g_hi:+.4f}]"
      f"  p={p_g:.4e}  N_right={nr}  N_left={nl}")

# ════════════════════════════════════════════════════════════════════════
# 2. 加载 CSV，合并 g_factor 行
# ════════════════════════════════════════════════════════════════════════
df = pd.read_csv(MEANS_CSV)
m2 = pd.read_csv(M2_CSV)

df["task"] = df["task"].replace({"PAL_errors": "PAL_correct"})
m2["task"] = m2["task"].replace({"PAL_errors": "PAL_correct"})

g_gmrow = {
    "task": "g_factor", "outcome_col": "g_factor",
    "direction": "higher_better", "reversed_for_plot": False,
    "special_note": "CFA g-score (7 indicators, MLR)",
    "n_right": nr, "n_left": nl,
    "raw_mean_right": mr, "raw_sd_right": sr,
    "raw_mean_left": ml, "raw_sd_left": sl,
    "effect_method": "OLS-residual-d",
    "cohen_d": d_g, "cohen_d_ci_lo": d_g_lo, "cohen_d_ci_hi": d_g_hi,
}
if "g_factor" not in df["task"].values:
    df = pd.concat([df, pd.DataFrame([g_gmrow])], ignore_index=True)
    print(f"  g_factor added → {MEANS_CSV.name}")
else:
    for col, val in g_gmrow.items():
        df.loc[df["task"] == "g_factor", col] = val
    print(f"  g_factor updated → {MEANS_CSV.name}")

g_m2row = {
    "task": "g_factor", "outcome_col": "g_factor",
    "method": "OLS-residual-d",
    "n_total": nr + nl, "n_right": nr, "n_left": nl,
    "beta": d_g * s_pool, "se": se_d_g * s_pool,
    "stat": t_g, "p_raw": p_g,
    "ci_lo": d_g_lo, "ci_hi": d_g_hi,
    "std_beta": d_g, "note": "CFA g-factor score (7 indicators)",
}
if "g_factor" not in m2["task"].values:
    m2 = pd.concat([m2, pd.DataFrame([g_m2row])], ignore_index=True)
else:
    for col, val in g_m2row.items():
        m2.loc[m2["task"] == "g_factor", col] = val

# 持久化写入
df_save = df.copy()
df_save.to_csv(MEANS_CSV, index=False)
m2.to_csv(M2_CSV, index=False)
print(f"  Saved {MEANS_CSV.name}  ({len(df_save)} rows)")
print(f"  Saved {M2_CSV.name}  ({len(m2)} rows)")

# ════════════════════════════════════════════════════════════════════════
# 3. BH-FDR 矫正（全部 15 个任务统一）
# ════════════════════════════════════════════════════════════════════════
m2_plot = m2[m2["task"].isin(ALL_TASKS)].copy()
m2_plot["_ord"] = m2_plot["task"].map({t: i for i, t in enumerate(ALL_TASKS)})
m2_plot = m2_plot.sort_values("_ord").reset_index(drop=True)

_, p_fdr_arr, _, _ = multipletests(m2_plot["p_raw"].astype(float), method="fdr_bh")
m2_plot["p_fdr"] = p_fdr_arr

def _sig(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""

m2_plot["sig_fdr"] = m2_plot["p_fdr"].apply(_sig)

print("\n=== 15-task unified BH-FDR results ===")
for _, r in m2_plot.iterrows():
    s = r["sig_fdr"] if r["sig_fdr"] else "ns"
    print(f"  {r['task']:<24} p_raw={float(r['p_raw']):.3e}  p_fdr={r['p_fdr']:.4f}  {s}")

# 合并 sig_fdr 到 df（仅内存用于绘图）
sig_map = m2_plot.set_index("task")["sig_fdr"].fillna("").to_dict()
df["sig_fdr"] = df["task"].map(sig_map).fillna("")

# ════════════════════════════════════════════════════════════════════════
# 4. 构建绘图条目（top→bottom, invert y-axis）
# ════════════════════════════════════════════════════════════════════════
# item dict: task | label | is_header | is_gfactor | is_spacer | y
items = []

items.append(dict(task=None, label="g-factor Indicators  (7 tasks)",
                  is_header=True, is_gfactor=False, is_spacer=False))
for t in G_TASKS:
    if df[df["task"] == t].shape[0] == 0:
        continue
    items.append(dict(task=t, label=LABELS[t],
                      is_header=False, is_gfactor=False, is_spacer=False))

items.append(dict(task=None, label="", is_header=False,
                  is_gfactor=False, is_spacer=True))
items.append(dict(task="g_factor", label=LABELS["g_factor"],
                  is_header=False, is_gfactor=True, is_spacer=False))
items.append(dict(task=None, label="", is_header=False,
                  is_gfactor=False, is_spacer=True))

items.append(dict(task=None, label="Other Cognitive Tasks  (not in g-factor)",
                  is_header=True, is_gfactor=False, is_spacer=False))
for t in OTHER_TASKS:
    if df[df["task"] == t].shape[0] == 0:
        continue
    items.append(dict(task=t, label=LABELS[t],
                      is_header=False, is_gfactor=False, is_spacer=False))

# y 坐标（0 = top，向下增大；最终 invert_yaxis）
y = 0.0
for it in items:
    it["y"] = y
    y += 0.55 if it["is_spacer"] else 1.0
total_y = y

# ════════════════════════════════════════════════════════════════════════
# 5. 绘图
# ════════════════════════════════════════════════════════════════════════
C_SIG_LEFT  = "#1a6faf"
C_SIG_RIGHT = "#c0392b"
C_NS        = "#888888"
C_GFACTOR   = "#6a0dad"   # 深紫 — General Factor 专用

METHOD_MARKER = {"OLS-residual-d": "s", "logit-d": "^", "tobit-d": "D"}

def get_color(is_gf, sig, d):
    if is_gf: return C_GFACTOR
    if sig in ("*","**","***"):
        return C_SIG_LEFT if d > 0 else C_SIG_RIGHT
    return C_NS

# x 范围
data_items = [it for it in items
              if it.get("task") and not it["is_header"] and not it["is_spacer"]]
lo_vals = [float(df[df["task"]==it["task"]].iloc[0]["cohen_d_ci_lo"])
           for it in data_items if df[df["task"]==it["task"]].shape[0]]
hi_vals = [float(df[df["task"]==it["task"]].iloc[0]["cohen_d_ci_hi"])
           for it in data_items if df[df["task"]==it["task"]].shape[0]]
x_lo = min(lo_vals) - 0.04
x_hi = max(hi_vals) + 0.04

fig_h = max(9.0, total_y * 0.52 + 2.5)
fig, ax = plt.subplots(figsize=(10.5, fig_h))

# 背景色
gf_ys = [it["y"] for it in items if it.get("task") in G_TASKS]
if gf_ys:
    ax.axhspan(min(gf_ys) - 0.45, max(gf_ys) + 0.45, color="#e8f0fb", zorder=0)
gf_row = next((it for it in items if it.get("task") == "g_factor"), None)
if gf_row:
    ax.axhspan(gf_row["y"] - 0.48, gf_row["y"] + 0.48, color="#f0e6ff", zorder=0)

ax.axvline(0, color="#333", lw=1.0, ls="--", alpha=0.6, zorder=1)

y_ticks, y_labels = [], []

for it in items:
    yp = it["y"]
    if it["is_spacer"]:
        continue
    if it["is_header"]:
        y_ticks.append(yp)
        y_labels.append(it["label"])
        continue

    task = it["task"]
    row  = df[df["task"] == task]
    if row.shape[0] == 0:
        continue
    row = row.iloc[0]

    d      = float(row["cohen_d"])
    lo     = float(row["cohen_d_ci_lo"])
    hi     = float(row["cohen_d_ci_hi"])
    sig    = str(row.get("sig_fdr", "")) if pd.notna(row.get("sig_fdr", "")) else ""
    method = str(row.get("effect_method", "OLS-residual-d"))
    is_gf  = it["is_gfactor"]

    color  = get_color(is_gf, sig, d)
    lw     = 2.8 if is_gf else (2.0 if sig else 1.2)
    ms     = 10  if is_gf else (7  if sig else 5)
    marker = "D" if is_gf else METHOD_MARKER.get(method, "s")

    ax.plot([lo, hi], [yp, yp], color=color, lw=lw,
            solid_capstyle="round", zorder=2)
    ax.plot(d, yp, marker=marker, color=color, ms=ms,
            mew=0.8, mec="white", zorder=3)

    y_ticks.append(yp)
    y_labels.append(it["label"])

# 右侧数值标注
for it in data_items:
    task = it["task"]
    row  = df[df["task"] == task]
    if row.shape[0] == 0:
        continue
    row     = row.iloc[0]
    d       = float(row["cohen_d"])
    lo      = float(row["cohen_d_ci_lo"])
    hi      = float(row["cohen_d_ci_hi"])
    sig     = str(row.get("sig_fdr", "")) if pd.notna(row.get("sig_fdr", "")) else ""
    is_gf   = it["is_gfactor"]
    color   = get_color(is_gf, sig, d)
    sig_str = sig if sig else "ns"

    ax.text(
        x_hi + 0.005, it["y"],
        f"{d:+.3f} [{lo:+.3f}, {hi:+.3f}]  {sig_str}",
        va="center", ha="left", fontsize=6.5, color=color,
        fontfamily="monospace", zorder=4, clip_on=False,
        fontweight="bold" if is_gf else "normal",
    )

# 轴
ax.set_xlim(x_lo, x_hi)
ax.set_ylim(total_y + 0.8, -0.8)   # 小 y 在上（已 invert）
ax.set_yticks(y_ticks)

tick_objs = ax.set_yticklabels(y_labels, fontsize=8.5)
headers_set = {it["label"] for it in items if it["is_header"]}
for txt in tick_objs:
    lbl = txt.get_text()
    if lbl in headers_set:
        txt.set_fontweight("bold"); txt.set_fontstyle("italic")
        txt.set_fontsize(8.5);     txt.set_color("#334466")
    elif lbl == LABELS["g_factor"]:
        txt.set_fontweight("bold")
        txt.set_color(C_GFACTOR)
        txt.set_fontsize(9.0)

ax.set_xlabel("Cohen's $d$  (left − right-handers, covariate-adjusted)",
              fontsize=10)
ax.set_title(
    "Handedness and Cognition: Effect Size Forest Plot\n"
    "(BH-FDR corrected across all 15 tasks;  positive = left-hander advantage)",
    fontsize=11, fontweight="bold", pad=10,
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)

# 图例
legend_elements = [
    Line2D([0],[0], marker="s", color="w", markerfacecolor=C_SIG_LEFT,
           markersize=8, label="Sig. left-hander advantage  (FDR < 0.05)"),
    Line2D([0],[0], marker="s", color="w", markerfacecolor=C_SIG_RIGHT,
           markersize=8, label="Sig. right-hander advantage (FDR < 0.05)"),
    Line2D([0],[0], marker="s", color="w", markerfacecolor=C_NS,
           markersize=8, label="Non-significant"),
    Line2D([0],[0], marker="D", color="w", markerfacecolor=C_GFACTOR,
           markersize=9, label="General Factor  g  (CFA composite, ◆)"),
    Line2D([0],[0], marker="^", color="w", markerfacecolor="#555",
           markersize=8, label="Logit-d  (▲ ProspMem †)"),
]
ax.legend(handles=legend_elements, loc="lower right",
          fontsize=7.5, framealpha=0.9, edgecolor="#cccccc")

fig.text(
    0.01, 0.005,
    "† ProspMem: ordered logistic d = log-OR / (π/√3).\n"
    "g-factor: CFA factor score from 7 age/centre-residualized indicators (MLR estimator).",
    fontsize=6.5, color="#777777", va="bottom",
)

plt.tight_layout(rect=[0, 0.04, 0.77, 1.0])

# ════════════════════════════════════════════════════════════════════════
# 6. 保存图像
# ════════════════════════════════════════════════════════════════════════
for ext in ("svg", "pdf", "png"):
    out = BASE / f"forest_gfactor.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  saved → {out}")

plt.close()
print("\n完成。")
