"""
02_abandonment_stratified.py
============================
Stratified competing-risks analysis of outcomes after TB treatment abandonment.

Corrections vs. v1:
  - Administrative censoring fixed to December 31, 2024.
  - Separate CIF figures for re-notification and death.
  - Stratified CIF curves: age group, sex, HIV status, incarceration,
    homelessness (address_type).

Outputs
-------
  figures/cif_renotification_overall.png
  figures/cif_death_overall.png
  figures/cif_renotification_by_<stratum>.png  (5 figures)
  figures/cif_death_by_<stratum>.png           (5 figures)
  data/abandonment_cohort_v2.csv               (cohort with strata)
  data/abandonment_risk_table.csv              (1/2/3-yr risks by subgroup)
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
COHORT_OUT  = DATA_DIR / "abandonment_cohort_v2.csv"
TABLE_OUT   = DATA_DIR / "abandonment_risk_table.csv"

CENSOR_DATE = pd.Timestamp("2024-12-31")
ABANDON_OUTCOMES = {"Abandono", "Abandono Primario"}

# ── 1. Load full dataset ───────────────────────────────────────────────────────
print("Loading cleaned dataset…")
df = pd.read_csv(CSV_CLEANED, low_memory=False)
print(f"  {len(df):,} rows  |  {df['sinan_clean'].nunique():,} individuals")

date_cols = ["notification_date", "end_date", "dod"]
for col in date_cols:
    df[col] = pd.to_datetime(df[col], errors="coerce")

df = df.sort_values(["sinan_clean", "notification_date"])

# ── 2. Build abandonment index episodes ───────────────────────────────────────
# Include key covariates from the index episode
covariate_cols = ["hiv", "aids", "address_type", "age_tb", "sex",
                  "dob", "race", "dod"]

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
print(f"\nAbandonment cohort (first episode, valid end_date): {len(aband_raw):,}")

# ── 3. Find next notification after abandonment ──────────────────────────────
sub_events = df.merge(aband_raw[["sinan_clean", "abandon_end"]], on="sinan_clean")
sub_events = sub_events[sub_events["notification_date"] > sub_events["abandon_end"]]

next_notif = (
    sub_events.groupby("sinan_clean")["notification_date"].min()
    .reset_index().rename(columns={"notification_date": "next_notif_date"})
)
aband = aband_raw.merge(next_notif, on="sinan_clean", how="left")

# ── 4. Resolve death date (any row for that individual) ───────────────────────
dod_by_indiv = (
    df.dropna(subset=["dod"])
    .groupby("sinan_clean")["dod"].min()
    .reset_index().rename(columns={"dod": "dod_any"})
)
aband = aband.merge(dod_by_indiv, on="sinan_clean", how="left")
aband["death_date"] = pd.to_datetime(
    np.where(aband["dod_any"] > aband["abandon_end"], aband["dod_any"], pd.NaT)
)

# ── 5. Apply administrative censoring (Dec 31, 2024) ─────────────────────────
# Cap outcomes at censor date
aband.loc[aband["next_notif_date"] > CENSOR_DATE, "next_notif_date"] = pd.NaT
aband.loc[aband["death_date"]      > CENSOR_DATE, "death_date"]      = pd.NaT

# Exclude individuals who abandoned AFTER the censor date
aband = aband[aband["abandon_end"] <= CENSOR_DATE].copy()
print(f"After excluding post-censor abandonments: {len(aband):,}")

# ── 6. Assign competing-risks outcome ─────────────────────────────────────────
# event = 0 censored | 1 re-notification | 2 death
def assign_outcome(row):
    t_n = row["next_notif_date"]
    t_d = row["death_date"]
    has_n = pd.notna(t_n)
    has_d = pd.notna(t_d)
    if has_n and has_d:
        return (1, t_n) if t_n <= t_d else (2, t_d)
    elif has_n:
        return (1, t_n)
    elif has_d:
        return (2, t_d)
    else:
        return (0, CENSOR_DATE)

outcomes = aband.apply(assign_outcome, axis=1, result_type="expand")
outcomes.columns = ["event", "event_date"]
aband = pd.concat([aband, outcomes], axis=1)
aband["event_date"] = pd.to_datetime(aband["event_date"])

# Time in years from abandonment end
aband["time_years"] = (aband["event_date"] - aband["abandon_end"]).dt.days / 365.25

# Drop implausible rows
bad = aband["time_years"] <= 0
print(f"Dropping {bad.sum()} rows with time ≤ 0")
aband = aband[~bad].copy()

print(f"\nFinal cohort: {len(aband):,}")
print(f"  Censored (0)    : {(aband['event']==0).sum():,}")
print(f"  Re-notified (1) : {(aband['event']==1).sum():,}")
print(f"  Died (2)        : {(aband['event']==2).sum():,}")

# ── 7. Derive stratification variables ───────────────────────────────────────

# Age group
def age_group(a):
    if pd.isna(a):   return "Unknown"
    elif a < 18:     return "<18"
    elif a < 35:     return "18–34"
    elif a < 50:     return "35–49"
    elif a < 65:     return "50–64"
    else:            return "65+"
aband["age_group"] = aband["age_tb"].apply(age_group)

# Sex
aband["sex_grp"] = aband["sex"].map({"M": "Male", "F": "Female"}).fillna("Unknown")

# HIV positive
def hiv_group(row):
    if row["aids"] == "S":         return "AIDS"
    if row["hiv"] == "Pos":        return "HIV+"
    if row["hiv"] == "Neg":        return "HIV-"
    return "Unknown/not tested"
aband["hiv_grp"] = aband.apply(hiv_group, axis=1)

# Incarceration (address_type == "DETENTO")
aband["incarcerated"] = aband["address_type"].apply(
    lambda x: "Incarcerated" if str(x).strip().upper() == "DETENTO" else "Not incarcerated"
)

# Homelessness (address_type == "SEM RESIDENCIA FIXA")
aband["homeless"] = aband["address_type"].apply(
    lambda x: "Homeless" if "SEM RESIDENCIA" in str(x).strip().upper() else "Not homeless"
)

aband.to_csv(COHORT_OUT, index=False)
print(f"\nCohort saved → {COHORT_OUT}")

# ── 8. Core functions ─────────────────────────────────────────────────────────

def aalen_johansen(times, events, cause):
    """Aalen-Johansen CIF estimator for one cause."""
    order   = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]
    unique_t = np.unique(t_sorted)
    S  = 1.0
    CI = 0.0
    n  = len(t_sorted)
    pts = [(0.0, 0.0)]
    i = 0
    for t in unique_t:
        mask = t_sorted == t
        d_all   = np.sum(e_sorted[mask] > 0)
        d_cause = np.sum(e_sorted[mask] == cause)
        if n > 0 and d_all > 0:
            CI += S * (d_cause / n)
            S  *= 1 - d_all / n
        n -= mask.sum()
        pts.append((t, CI))
    t_arr  = np.array([p[0] for p in pts])
    ci_arr = np.array([p[1] for p in pts])
    return t_arr, ci_arr

def interp_cif(t_arr, ci_arr, year):
    idx = np.searchsorted(t_arr, year, side="right") - 1
    return ci_arr[max(0, min(idx, len(ci_arr)-1))]

def cif_for_group(sub, cause):
    return aalen_johansen(sub["time_years"].values, sub["event"].values, cause)

FOLLOW_YEARS = [1, 2, 3]
MAX_PLOT_YEARS = 5

# ── 9. Plot helpers ───────────────────────────────────────────────────────────

# Colour palettes
PALETTES = {
    "age_group":    ["#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#B07AA1"],
    "sex_grp":      ["#4E79A7","#E15759","#888"],
    "hiv_grp":      ["#B07AA1","#E15759","#4E79A7","#aaa"],
    "incarcerated": ["#E15759","#4E79A7"],
    "homeless":     ["#F28E2B","#4E79A7"],
}

ORDER = {
    "age_group":    ["<18","18–34","35–49","50–64","65+","Unknown"],
    "sex_grp":      ["Male","Female","Unknown"],
    "hiv_grp":      ["AIDS","HIV+","HIV-","Unknown/not tested"],
    "incarcerated": ["Incarcerated","Not incarcerated"],
    "homeless":     ["Homeless","Not homeless"],
}

STRATUM_LABELS = {
    "age_group":    "Age group",
    "sex_grp":      "Sex",
    "hiv_grp":      "HIV/AIDS status",
    "incarcerated": "Incarceration",
    "homeless":     "Housing status",
}

CAUSE_LABELS = {1: "Re-notification", 2: "Death"}
CAUSE_COLORS_OVERALL = {1: "#E05C5C", 2: "#5C7AE0"}


def plot_cif(aband, strat_col, cause, ax, palette):
    """
    Draw stratified CIF curves on ax.
    Returns list of (label, n, cif_vals_at_years) for the risk table.
    """
    groups  = [g for g in ORDER[strat_col] if g in aband[strat_col].unique()]
    rows    = []
    for grp, color in zip(groups, palette):
        sub = aband[aband[strat_col] == grp]
        if len(sub) < 20:
            continue
        t_arr, ci_arr = cif_for_group(sub, cause)
        ax.step(t_arr, ci_arr * 100, where="post", lw=2, color=color, label=f"{grp} (n={len(sub):,})")

        row = {"subgroup": grp, "n": len(sub)}
        for yr in FOLLOW_YEARS:
            # only report if enough follow-up
            row[f"risk_{yr}yr_%"] = round(interp_cif(t_arr, ci_arr, yr) * 100, 1)
        rows.append(row)

    for yr in FOLLOW_YEARS:
        ax.axvline(yr, color="grey", lw=0.7, ls=":")

    ax.set_xlim(0, MAX_PLOT_YEARS)
    ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=10)
    ax.set_ylabel(f"Cumulative Incidence of {CAUSE_LABELS[cause]} (%)", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.85)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    return rows


# ── 10. Overall CIF figures ───────────────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")

for cause in [1, 2]:
    fig, ax = plt.subplots(figsize=(8, 5))
    t_arr, ci_arr = cif_for_group(aband, cause)
    ax.step(t_arr, ci_arr * 100, where="post", lw=2.5,
            color=CAUSE_COLORS_OVERALL[cause], label="Overall")
    ax.fill_between(t_arr, ci_arr * 100, step="post",
                    alpha=0.12, color=CAUSE_COLORS_OVERALL[cause])

    for yr in FOLLOW_YEARS:
        val = interp_cif(t_arr, ci_arr, yr) * 100
        ax.axvline(yr, color="grey", lw=0.7, ls=":")
        ax.annotate(f"{val:.1f}%", xy=(yr, val),
                    xytext=(yr + 0.08, val + 0.3),
                    fontsize=10, fontweight="bold",
                    color=CAUSE_COLORS_OVERALL[cause])

    ax.set_xlim(0, MAX_PLOT_YEARS)
    ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=12)
    ax.set_ylabel(f"Cumulative Incidence (%)", fontsize=12)
    ax.set_title(
        f"Risk of {CAUSE_LABELS[cause]} After TB Treatment Abandonment\n"
        f"(N = {len(aband):,}; censored Dec 31, 2024)",
        fontsize=13, fontweight="bold"
    )
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    ax.legend(fontsize=10)
    plt.tight_layout()
    fname = FIG_DIR / f"cif_{'renotification' if cause==1 else 'death'}_overall.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fname}")

# ── 11. Stratified figures ────────────────────────────────────────────────────
strata = ["age_group", "sex_grp", "hiv_grp", "incarcerated", "homeless"]

all_risk_rows = []

for strat_col in strata:
    label = STRATUM_LABELS[strat_col]
    palette = PALETTES[strat_col]

    for cause in [1, 2]:
        fig, ax = plt.subplots(figsize=(9, 5.5))

        rows = plot_cif(aband, strat_col, cause, ax, palette)

        ax.set_title(
            f"Risk of {CAUSE_LABELS[cause]} by {label}\n"
            f"After TB Treatment Abandonment  (N = {len(aband):,}; censored Dec 31, 2024)",
            fontsize=12, fontweight="bold"
        )
        plt.tight_layout()
        cname = "renotification" if cause == 1 else "death"
        fname = FIG_DIR / f"cif_{cname}_by_{strat_col}.png"
        fig.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {fname}")

        for row in rows:
            row["outcome"]  = CAUSE_LABELS[cause]
            row["stratum"]  = label
            row["variable"] = strat_col
        all_risk_rows.extend(rows)

# ── 12. Risk table ────────────────────────────────────────────────────────────
# Add overall rows
for cause in [1, 2]:
    t_arr, ci_arr = cif_for_group(aband, cause)
    row = {"stratum": "Overall", "variable": "overall",
           "subgroup": "All", "n": len(aband),
           "outcome": CAUSE_LABELS[cause]}
    for yr in FOLLOW_YEARS:
        row[f"risk_{yr}yr_%"] = round(interp_cif(t_arr, ci_arr, yr) * 100, 1)
    all_risk_rows.insert(0 if cause == 1 else 1, row)

risk_table = pd.DataFrame(all_risk_rows)[
    ["outcome", "stratum", "subgroup", "n",
     "risk_1yr_%", "risk_2yr_%", "risk_3yr_%"]
]
risk_table.to_csv(TABLE_OUT, index=False)
print(f"\nRisk table saved → {TABLE_OUT}")

print("\n── Risk Table ───────────────────────────────────────────────────────────")
for outcome in ["Re-notification", "Death"]:
    print(f"\n  {outcome}:")
    sub = risk_table[risk_table["outcome"] == outcome]
    print(sub.drop(columns="outcome").to_string(index=False))

print("\nDone ✓")
