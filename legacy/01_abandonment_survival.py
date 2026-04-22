"""
01_abandonment_survival.py
==========================
Competing-risks survival analysis for individuals who abandoned TB treatment.

Index event  : First abandonment episode (case_outcome = "Abandono" or
               "Abandono Primario") per individual; time zero = end_date.
Outcomes     :
  1. Re-notification  – any new notification (notification_date) for the
                        same individual that occurs AFTER the abandonment
                        end_date.  We capture the earliest such event.
  2. Death            – date of death (dod) > abandonment end_date.
     When BOTH a re-notification and death exist, we take whichever is first
     (competing risks framework).
Censoring    : Individuals with neither outcome are censored at the latest
               end_date observed in the whole dataset (administrative
               censoring) or their own last known end_date, whichever is
               earlier.

Outputs
-------
  data/abandonment_cohort.csv       – one row per individual (cohort)
  data/abandonment_survival_table.csv – cumulative incidence table
  figures/abandonment_cuminc.png    – Cumulative Incidence Function plot
  figures/abandonment_km.png        – Kaplan-Meier (event-free survival) plot
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
DATA_DIR    = Path("data")
FIG_DIR     = Path("figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

CSV_IN      = DATA_DIR / "Final_table_cleaned.csv"
COHORT_OUT  = DATA_DIR / "abandonment_cohort.csv"
TABLE_OUT   = DATA_DIR / "abandonment_survival_table.csv"

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("Loading cleaned dataset…")
df = pd.read_csv(CSV_IN, low_memory=False)
print(f"  {len(df):,} rows  |  {df['sinan_clean'].nunique():,} individuals")

# Parse dates
date_cols = ["notification_date", "end_date", "dod"]
for col in date_cols:
    df[col] = pd.to_datetime(df[col], errors="coerce")

df = df.sort_values(["sinan_clean", "notification_date"])

# ── 2. Identify abandonment episodes ─────────────────────────────────────────
ABANDON_OUTCOMES = {"Abandono", "Abandono Primario"}
aband_mask = df["case_outcome"].isin(ABANDON_OUTCOMES)

# Take the FIRST abandonment per individual (earliest end_date among abandonment rows)
aband = (
    df[aband_mask]
    .dropna(subset=["end_date"])                     # need a valid time-zero
    .sort_values("end_date")
    .groupby("sinan_clean", as_index=False)
    .first()                                          # first abandonment episode
    .rename(columns={
        "end_date":           "abandon_end",
        "notification_date":  "index_notif_date",
        "case_outcome":       "index_outcome",
        "case_type":          "index_case_type",
        "dod":                "index_dod",           # dod recorded on that row
    })
    [["sinan_clean", "abandon_end", "index_notif_date",
      "index_outcome", "index_case_type", "index_dod",
      "sex", "dob", "age_tb", "race"]]
    .copy()
)
print(f"\nAbandonment cohort (first abandonment per person, valid end_date):")
print(f"  {len(aband):,} individuals")

# ── 3. Find subsequent events ─────────────────────────────────────────────────
# For each individual in the cohort, look at ALL their rows after abandon_end

# Merge the abandon_end back into the full dataset
sub_events = df.merge(
    aband[["sinan_clean", "abandon_end"]], on="sinan_clean"
)

# Only rows whose notification_date is STRICTLY after abandonment end
sub_events = sub_events[sub_events["notification_date"] > sub_events["abandon_end"]]

# Earliest subsequent notification per individual
next_notif = (
    sub_events.groupby("sinan_clean")["notification_date"]
    .min()
    .reset_index()
    .rename(columns={"notification_date": "next_notif_date"})
)

aband = aband.merge(next_notif, on="sinan_clean", how="left")

# ── 4. Resolve death date ─────────────────────────────────────────────────────
# dod may appear on ANY row for that individual (often duplicated across rows).
# Take the minimum non-null dod across all rows per individual.
dod_by_indiv = (
    df.dropna(subset=["dod"])
    .groupby("sinan_clean")["dod"]
    .min()
    .reset_index()
    .rename(columns={"dod": "dod_any"})
)
aband = aband.merge(dod_by_indiv, on="sinan_clean", how="left")

# Only count death if it is AFTER the abandonment end_date
aband["death_date"] = np.where(
    aband["dod_any"] > aband["abandon_end"],
    aband["dod_any"],
    pd.NaT
)
aband["death_date"] = pd.to_datetime(aband["death_date"])

# ── 5. Build competing-risks outcome ─────────────────────────────────────────
# event = 0  censored
# event = 1  re-notification (first if before death)
# event = 2  death (first if before re-notification)

def assign_outcome(row):
    t_notif = row["next_notif_date"]
    t_death = row["death_date"]
    has_notif = pd.notna(t_notif)
    has_death = pd.notna(t_death)
    if has_notif and has_death:
        return (1, t_notif) if t_notif <= t_death else (2, t_death)
    elif has_notif:
        return (1, t_notif)
    elif has_death:
        return (2, t_death)
    else:
        return (0, pd.NaT)

outcomes = aband.apply(assign_outcome, axis=1, result_type="expand")
outcomes.columns = ["event", "event_date"]
aband = pd.concat([aband, outcomes], axis=1)

# Administrative censoring date = max end_date in entire dataset
censor_date = df["end_date"].max()
print(f"\nAdministrative censoring date: {censor_date.date()}")

# For censored, event_date = administrative censoring
aband.loc[aband["event"] == 0, "event_date"] = censor_date
aband["event_date"] = pd.to_datetime(aband["event_date"])

# ── 6. Compute time-to-event (years) ─────────────────────────────────────────
aband["time_years"] = (
    (aband["event_date"] - aband["abandon_end"]).dt.days / 365.25
)

# Drop rare negative times (data errors: event before abandonment end)
bad = aband["time_years"] <= 0
print(f"Dropping {bad.sum()} rows with time ≤ 0 (data inconsistency)")
aband = aband[~bad].copy()

print(f"\nFinal cohort: {len(aband):,} individuals")
print(f"  Censored (0)      : {(aband['event']==0).sum():,}")
print(f"  Re-notified (1)   : {(aband['event']==1).sum():,}")
print(f"  Died (2)          : {(aband['event']==2).sum():,}")

median_follow = aband["time_years"].median()
print(f"  Median follow-up  : {median_follow:.2f} years")

aband.to_csv(COHORT_OUT, index=False)
print(f"\nCohort saved → {COHORT_OUT}")

# ── 7. Kaplan-Meier (event-free survival) ────────────────────────────────────
# Treats any event (re-notification OR death) as failure

def kaplan_meier(times, events_any):
    """Returns (time_points, survival) arrays."""
    df_km = pd.DataFrame({"t": times, "e": events_any}).sort_values("t")
    times_unique = np.sort(df_km["t"].unique())
    S = 1.0
    results = [(0.0, 1.0)]
    n_risk = len(df_km)
    for t in times_unique:
        at_t   = df_km["t"] == t
        d      = df_km.loc[at_t & (df_km["e"] == 1), "e"].sum()
        n      = n_risk
        if n > 0 and d > 0:
            S *= (1 - d / n)
        n_risk -= at_t.sum()
        results.append((t, S))
    return zip(*results)

km_times, km_surv = kaplan_meier(
    aband["time_years"].values,
    (aband["event"] > 0).astype(int).values
)
km_t = np.array(list(km_times))
km_s = np.array(list(km_surv))

# ── 8. Cumulative Incidence Functions (Aalen-Johansen) ───────────────────────
def aalen_johansen(times, events, cause):
    """
    Simple Aalen-Johansen estimator for a single cause in a two-cause
    competing-risks setting.
    Returns (time_points, CIF).
    """
    df_aj = pd.DataFrame({"t": times, "e": events}).sort_values("t")
    times_unique = np.sort(df_aj["t"].unique())

    S_overall = 1.0  # overall survival (no event of any kind)
    CIF = 0.0
    n_risk = len(df_aj)
    results = [(0.0, 0.0)]

    for t in times_unique:
        at_t   = df_aj["t"] == t
        d_all  = (df_aj.loc[at_t, "e"] > 0).sum()
        d_cause= (df_aj.loc[at_t, "e"] == cause).sum()
        n      = n_risk
        if n > 0 and d_all > 0:
            # Increment in CIF = S(t-) × h_cause(t)
            CIF += S_overall * (d_cause / n)
            S_overall *= (1 - d_all / n)
        n_risk -= at_t.sum()
        results.append((t, CIF))
    return zip(*results)

cif1_t_z, cif1_cif_z = aalen_johansen(aband["time_years"].values, aband["event"].values, cause=1)
cif2_t_z, cif2_cif_z = aalen_johansen(aband["time_years"].values, aband["event"].values, cause=2)
cif1_t = np.array(list(cif1_t_z)); cif1 = np.array(list(cif1_cif_z))
cif2_t = np.array(list(cif2_t_z)); cif2 = np.array(list(cif2_cif_z))

# ── 9. Risk table at 1, 2, 3 years ───────────────────────────────────────────
def interp_at(t_arr, v_arr, years):
    idx = np.searchsorted(t_arr, years, side="right") - 1
    idx = np.clip(idx, 0, len(v_arr) - 1)
    return v_arr[idx]

rows = []
for yr in [1, 2, 3]:
    cif1_val = interp_at(cif1_t, cif1, yr)
    cif2_val = interp_at(cif2_t, cif2, yr)
    km_val   = interp_at(km_t, km_s, yr)
    n_at_risk = (aband["time_years"] >= yr).sum()
    rows.append({
        "year":                  yr,
        "CIF_renotified_%":      round(cif1_val * 100, 1),
        "CIF_death_%":           round(cif2_val * 100, 1),
        "CIF_either_%":          round((1 - km_val) * 100, 1),  # 1 - KM
        "EventFree_Surv_%":      round(km_val * 100, 1),
        "n_at_risk":             n_at_risk,
    })

table = pd.DataFrame(rows)
table.to_csv(TABLE_OUT, index=False)

print("\n── Cumulative Incidence Table ────────────────────────────────────────────")
print(table.to_string(index=False))

# ── 10. Plot: CIF (competing risks) ──────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    "Outcomes After TB Treatment Abandonment\n"
    f"(N = {len(aband):,} individuals; index event = first abandonment)",
    fontsize=14, fontweight="bold", y=1.01
)

## -- Panel A: Cumulative Incidence Functions --
ax = axes[0]
ax.step(cif1_t, cif1 * 100, where="post", color="#E05C5C", lw=2,
        label="Re-notification (any new episode)")
ax.step(cif2_t, cif2 * 100, where="post", color="#5C7AE0", lw=2,
        label="Death")
total_cif = cif1 + cif2
ax.step(cif1_t, np.interp(cif1_t, cif2_t, cif2) * 100 + cif1 * 100,
        where="post", color="#888", lw=1.5, ls="--",
        label="Either (stacked)")

# Annotations at 1, 2, 3 years
for yr in [1, 2, 3]:
    for cif_arr, cif_t_arr, color in [
        (cif1, cif1_t, "#E05C5C"),
        (cif2, cif2_t, "#5C7AE0"),
    ]:
        val = interp_at(cif_t_arr, cif_arr, yr) * 100
        ax.annotate(f"{val:.1f}%",
                    xy=(yr, val), xytext=(yr + 0.05, val + 0.5),
                    fontsize=8, color=color)
    ax.axvline(yr, color="grey", lw=0.7, ls=":")

ax.set_xlim(0, 5)
ax.set_ylim(0, None)
ax.set_xlabel("Years since abandonment", fontsize=11)
ax.set_ylabel("Cumulative Incidence (%)", fontsize=11)
ax.set_title("Cumulative Incidence Functions\n(Competing Risks)", fontsize=12)
ax.legend(fontsize=9, loc="upper left")
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))

## -- Panel B: Kaplan-Meier event-free survival --
ax2 = axes[1]
ax2.step(km_t, km_s * 100, where="post", color="#2C7A3A", lw=2.5,
         label="Event-free survival")
ax2.fill_between(km_t, km_s * 100, step="post", alpha=0.1, color="#2C7A3A")

for yr in [1, 2, 3]:
    val = interp_at(km_t, km_s, yr) * 100
    ax2.axvline(yr, color="grey", lw=0.7, ls=":")
    ax2.annotate(f"{val:.1f}%",
                 xy=(yr, val), xytext=(yr + 0.05, val + 0.3),
                 fontsize=9, color="#2C7A3A", fontweight="bold")

ax2.set_xlim(0, 5)
ax2.set_ylim(0, 100)
ax2.set_xlabel("Years since abandonment", fontsize=11)
ax2.set_ylabel("Event-free Survival (%)", fontsize=11)
ax2.set_title("Kaplan-Meier: Event-free Survival\n(re-notification OR death)", fontsize=12)
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))

# Add risk table under KM plot
risk_years = [0, 1, 2, 3, 4, 5]
n_risk_vals = [int((aband["time_years"] >= yr).sum()) for yr in risk_years]
ax2.annotate(
    "N at risk: " + "   ".join(
        [f"{yr}yr: {n:,}" for yr, n in zip(risk_years[1:], n_risk_vals[1:])]
    ),
    xy=(0.01, 0.04), xycoords="axes fraction", fontsize=8, color="#444",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", alpha=0.8)
)

plt.tight_layout()
fig.savefig(FIG_DIR / "abandonment_cuminc.png", dpi=150, bbox_inches="tight")
print(f"\nFigure saved → {FIG_DIR / 'abandonment_cuminc.png'}")

# ── 11. Breakdown: What type of re-notification? ──────────────────────────────
# Among re-notified, what was the case_type of the next episode?
next_case_type = (
    sub_events.sort_values("notification_date")
    .groupby("sinan_clean")[["case_type", "notification_date"]]
    .first()
    .reset_index()
    .rename(columns={"case_type": "next_case_type",
                     "notification_date": "next_notif_date_check"})
)
renotif = aband[aband["event"] == 1].merge(next_case_type, on="sinan_clean", how="left")
print("\n── Case type of next episode (re-notified individuals) ──")
print(renotif["next_case_type"].value_counts().to_string())

print("\nDone ✓")
