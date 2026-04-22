"""
03_abandonment_by_tx_month.py
==============================
Analyzes whether the treatment month at which abandonment occurs predicts
subsequent risk of re-notification and death.

Month of abandonment = ceil((end_date - tx_start) / 30.44 days).
Groups: Month 1, 2, 3, 4, 5, 6, 7+ (months 8+ have small N and are pooled).

Uses the same cohort construction as 02_abandonment_stratified.py:
  - First abandonment per individual
  - Administrative censoring: December 31, 2024
  - Competing-risks (Aalen-Johansen) CIF for re-notification and death

Outputs
-------
  figures/cif_renotification_by_tx_month.png
  figures/cif_death_by_tx_month.png
  figures/tx_month_distribution.png
  data/abandonment_cohort_txmonth.csv
  data/abandonment_risk_table_txmonth.csv
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR = Path("data")
FIG_DIR  = Path("figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

CSV_CLEANED = DATA_DIR / "Final_table_cleaned.csv"
COHORT_OUT  = DATA_DIR / "abandonment_cohort_txmonth.csv"
TABLE_OUT   = DATA_DIR / "abandonment_risk_table_txmonth.csv"

CENSOR_DATE      = pd.Timestamp("2024-12-31")
ABANDON_OUTCOMES = {"Abandono", "Abandono Primario"}
FOLLOW_YEARS     = [1, 2, 3]
MAX_PLOT_YEARS   = 5

# ── 1. Load full dataset ───────────────────────────────────────────────────────
print("Loading cleaned dataset…")
df = pd.read_csv(CSV_CLEANED, low_memory=False)
print(f"  {len(df):,} rows  |  {df['sinan_clean'].nunique():,} individuals")

date_cols = ["notification_date", "end_date", "dod", "tx_start"]
for col in date_cols:
    df[col] = pd.to_datetime(df[col], errors="coerce")

df = df.sort_values(["sinan_clean", "notification_date"])

# ── 2. Build abandonment index episodes, keeping tx_start ─────────────────────
covariate_cols = ["hiv", "aids", "address_type", "age_tb", "sex",
                  "dob", "race", "dod", "tx_start"]

aband_raw = (
    df[df["case_outcome"].isin(ABANDON_OUTCOMES)]
    .dropna(subset=["end_date"])
    .sort_values("end_date")
    .groupby("sinan_clean", as_index=False)
    .first()
)[["sinan_clean", "end_date", "notification_date", "case_outcome",
   "case_type"] + covariate_cols]

aband_raw = aband_raw.rename(columns={
    "end_date":          "abandon_end",
    "notification_date": "index_notif_date",
    "case_outcome":      "index_outcome",
    "case_type":         "index_case_type",
    "dod":               "index_dod",
})

# Remove post-censor abandonments
aband_raw = aband_raw[aband_raw["abandon_end"] <= CENSOR_DATE].copy()
print(f"Abandonment cohort (first episode, valid end_date): {len(aband_raw):,}")

# ── 3. Compute treatment duration at abandonment ──────────────────────────────
aband_raw["tx_days"] = (
    (aband_raw["abandon_end"] - aband_raw["tx_start"]).dt.days
)

# Exclude rows with missing/implausible tx_start
n_orig = len(aband_raw)
aband_raw = aband_raw[aband_raw["tx_days"] > 0].copy()
n_valid   = len(aband_raw)
print(f"  Dropped {n_orig - n_valid} rows with tx_days ≤ 0 or missing tx_start")

# Month of abandonment (1-indexed); cap beyond month 7 as "7+"
aband_raw["tx_month_raw"] = np.ceil(aband_raw["tx_days"] / 30.44).astype(int)
aband_raw["tx_month_grp"] = aband_raw["tx_month_raw"].apply(
    lambda m: f"Month {m}" if m <= 6 else "Month 7+"
)

# Distribution
print("\nMonth-of-abandonment distribution:")
dist = aband_raw["tx_month_grp"].value_counts().sort_index()
print(dist.to_string())

# ── 4. Find next notification after abandonment ───────────────────────────────
sub_events = df.merge(aband_raw[["sinan_clean", "abandon_end"]], on="sinan_clean")
sub_events = sub_events[sub_events["notification_date"] > sub_events["abandon_end"]]

next_notif = (
    sub_events.groupby("sinan_clean")["notification_date"].min()
    .reset_index().rename(columns={"notification_date": "next_notif_date"})
)
aband = aband_raw.merge(next_notif, on="sinan_clean", how="left")

# ── 5. Resolve death date ─────────────────────────────────────────────────────
dod_by_indiv = (
    df.dropna(subset=["dod"])
    .groupby("sinan_clean")["dod"].min()
    .reset_index().rename(columns={"dod": "dod_any"})
)
aband = aband.merge(dod_by_indiv, on="sinan_clean", how="left")
aband["death_date"] = pd.to_datetime(
    np.where(aband["dod_any"] > aband["abandon_end"], aband["dod_any"], pd.NaT)
)

# Cap at censor date
aband.loc[aband["next_notif_date"] > CENSOR_DATE, "next_notif_date"] = pd.NaT
aband.loc[aband["death_date"]      > CENSOR_DATE, "death_date"]      = pd.NaT

# ── 6. Assign competing-risks outcome ─────────────────────────────────────────
def assign_outcome(row):
    t_n, t_d = row["next_notif_date"], row["death_date"]
    has_n, has_d = pd.notna(t_n), pd.notna(t_d)
    if has_n and has_d:
        return (1, t_n) if t_n <= t_d else (2, t_d)
    elif has_n: return (1, t_n)
    elif has_d: return (2, t_d)
    else:       return (0, CENSOR_DATE)

outcomes = aband.apply(assign_outcome, axis=1, result_type="expand")
outcomes.columns = ["event", "event_date"]
aband = pd.concat([aband, outcomes], axis=1)
aband["event_date"] = pd.to_datetime(aband["event_date"])
aband["time_years"] = (aband["event_date"] - aband["abandon_end"]).dt.days / 365.25

bad = aband["time_years"] <= 0
print(f"\nDropping {bad.sum()} rows with time ≤ 0")
aband = aband[~bad].copy()

print(f"\nFinal cohort: {len(aband):,}")
print(f"  Censored (0)    : {(aband['event']==0).sum():,}")
print(f"  Re-notified (1) : {(aband['event']==1).sum():,}")
print(f"  Died (2)        : {(aband['event']==2).sum():,}")

aband.to_csv(COHORT_OUT, index=False)
print(f"Cohort saved → {COHORT_OUT}")

# ── 7. Core estimator ─────────────────────────────────────────────────────────
def aalen_johansen(times, events, cause):
    order    = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]
    S, CI    = 1.0, 0.0
    n        = len(t_sorted)
    pts      = [(0.0, 0.0)]
    prev_t   = -np.inf
    for i, t in enumerate(t_sorted):
        if t == prev_t:
            continue
        mask    = t_sorted == t
        d_all   = np.sum(e_sorted[mask] > 0)
        d_cause = np.sum(e_sorted[mask] == cause)
        if n > 0 and d_all > 0:
            CI += S * (d_cause / n)
            S  *= 1 - d_all / n
        n     -= mask.sum()
        prev_t = t
        pts.append((t, CI))
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

def interp_cif(t_arr, ci_arr, year):
    idx = np.searchsorted(t_arr, year, side="right") - 1
    return float(ci_arr[max(0, min(idx, len(ci_arr)-1))])

# ── 8. Month group ordering & palette ────────────────────────────────────────
MONTH_ORDER = ["Month 1","Month 2","Month 3","Month 4","Month 5","Month 6","Month 7+"]
# Gradient from early (yellow-green) to late (deep purple)
PALETTE = ["#f7dc6f","#f0a500","#e05c5c","#b03060","#7b2d8b","#3d5fa0","#1a3a5c"]

CAUSE_LABELS  = {1: "Re-notification", 2: "Death"}
CAUSE_COLORS  = {1: "#E05C5C",         2: "#5C7AE0"}

# ── 9. Distribution figure ────────────────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: bar chart of n by month group
grp_counts = (
    aband.groupby("tx_month_grp")
    .size()
    .reindex(MONTH_ORDER)
    .fillna(0)
    .astype(int)
)
axes[0].bar(grp_counts.index, grp_counts.values,
            color=PALETTE, edgecolor="white", linewidth=0.8)
for i, (grp, n) in enumerate(grp_counts.items()):
    axes[0].text(i, n + 50, f"{n:,}", ha="center", va="bottom", fontsize=9)
axes[0].set_xlabel("Month of treatment at abandonment", fontsize=11)
axes[0].set_ylabel("Number of individuals", fontsize=11)
axes[0].set_title("Distribution of Abandonment\nby Treatment Month", fontsize=12, fontweight="bold")
axes[0].tick_params(axis="x", rotation=25)

# Panel B: % re-notified vs % died by month (1-yr risk)
re_notif_1yr = []
death_1yr    = []
for grp in MONTH_ORDER:
    sub = aband[aband["tx_month_grp"] == grp]
    if len(sub) < 20:
        re_notif_1yr.append(np.nan)
        death_1yr.append(np.nan)
        continue
    t1, c1 = aalen_johansen(sub["time_years"].values, sub["event"].values, 1)
    t2, c2 = aalen_johansen(sub["time_years"].values, sub["event"].values, 2)
    re_notif_1yr.append(interp_cif(t1, c1, 1) * 100)
    death_1yr.append(interp_cif(t2, c2, 1) * 100)

x = np.arange(len(MONTH_ORDER))
width = 0.38
ax2 = axes[1]
bars1 = ax2.bar(x - width/2, re_notif_1yr, width, color="#E05C5C",
                label="Re-notification", edgecolor="white")
bars2 = ax2.bar(x + width/2, death_1yr,    width, color="#5C7AE0",
                label="Death", edgecolor="white")
for bar, val in zip(bars1, re_notif_1yr):
    if not np.isnan(val):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.3,
                 f"{val:.1f}%", ha="center", va="bottom", fontsize=7.5, color="#a00")
ax2.set_xticks(x)
ax2.set_xticklabels(MONTH_ORDER, rotation=25)
ax2.set_ylabel("1-year Cumulative Incidence (%)", fontsize=11)
ax2.set_title("1-year Risk by Month of Abandonment", fontsize=12, fontweight="bold")
ax2.legend(fontsize=10)
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))

plt.suptitle("Treatment Month at Abandonment — Overview\n"
             f"(N = {len(aband):,}; censored Dec 31, 2024)",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR / "tx_month_distribution.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {FIG_DIR / 'tx_month_distribution.png'}")

# ── 10. CIF curves by tx month group ─────────────────────────────────────────
risk_rows = []

for cause in [1, 2]:
    fig, ax = plt.subplots(figsize=(10, 6))

    for grp, color in zip(MONTH_ORDER, PALETTE):
        sub = aband[aband["tx_month_grp"] == grp]
        if len(sub) < 20:
            continue
        t_arr, ci_arr = aalen_johansen(sub["time_years"].values,
                                       sub["event"].values, cause)
        ax.step(t_arr, ci_arr * 100, where="post", lw=2,
                color=color, label=f"{grp}  (n={len(sub):,})")

        row = {"outcome": CAUSE_LABELS[cause], "tx_month_group": grp, "n": len(sub)}
        for yr in FOLLOW_YEARS:
            row[f"risk_{yr}yr_%"] = round(interp_cif(t_arr, ci_arr, yr) * 100, 1)
        risk_rows.append(row)

    for yr in FOLLOW_YEARS:
        ax.axvline(yr, color="grey", lw=0.7, ls=":")

    ax.set_xlim(0, MAX_PLOT_YEARS)
    ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=12)
    ax.set_ylabel(f"Cumulative Incidence of {CAUSE_LABELS[cause]} (%)", fontsize=12)
    ax.set_title(
        f"Risk of {CAUSE_LABELS[cause]} by Treatment Month at Abandonment\n"
        f"(N = {len(aband):,}; censored Dec 31, 2024)",
        fontsize=13, fontweight="bold"
    )
    ax.legend(fontsize=9.5, loc="upper left", framealpha=0.9,
              title="Month of TX at abandonment", title_fontsize=9)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    plt.tight_layout()

    cname = "renotification" if cause == 1 else "death"
    fname = FIG_DIR / f"cif_{cname}_by_tx_month.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fname}")

# ── 11. Risk table ────────────────────────────────────────────────────────────
risk_table = pd.DataFrame(risk_rows)[
    ["outcome", "tx_month_group", "n",
     "risk_1yr_%", "risk_2yr_%", "risk_3yr_%"]
]
risk_table.to_csv(TABLE_OUT, index=False)
print(f"\nRisk table saved → {TABLE_OUT}")

print("\n── Risk Table ───────────────────────────────────────────────────────────")
for outcome in ["Re-notification", "Death"]:
    print(f"\n  {outcome}:")
    sub = risk_table[risk_table["outcome"] == outcome].drop(columns="outcome")
    print(sub.to_string(index=False))

print("\nDone ✓")
