"""
04_abandonment_full_analysis.py
================================
REVISED analysis with comprehensive death ascertainment.

Key methodological changes from previous versions:
  1. Death is now ascertained from TWO sources combined:
       a) 'dod' field (date of death) on any row for that individual
       b) 'case_outcome' = "Obito TB" or "Obito NTB" on any row
          → death date = end_date of that episode
     The earliest date across both sources is used.

  2. Re-notification and death are treated as INDEPENDENT outcomes:
       - Re-notification: time to first new episode after abandonment.
         Death before re-notification is a competing event (Aalen-Johansen).
       - Death: time to death from any source after abandonment.
         Re-notification does NOT censor for death — individuals are
         followed until death or Dec 31, 2024, regardless.
         This is a standard Kaplan-Meier / 1-KM analysis.

All stratified analyses (age, sex, HIV, incarceration, homelessness,
treatment month) are redone with these corrections.

Outputs  (figures/ and data/)
─────────────────────────────
Overall:
  cif_renotification_overall_v2.png
  cif_death_overall_v2.png

Stratified by demographic/clinical factors:
  cif_renotification_by_<stratum>_v2.png   (5 figs)
  cif_death_by_<stratum>_v2.png            (5 figs)

Stratified by treatment month:
  tx_month_distribution_v2.png
  cif_renotification_by_tx_month_v2.png
  cif_death_by_tx_month_v2.png

Tables:
  abandonment_cohort_v3.csv
  abandonment_risk_table_v3.csv
  abandonment_risk_table_txmonth_v3.csv
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

# ── Constants ─────────────────────────────────────────────────────────────────
DATA_DIR         = Path("data")
FIG_DIR          = Path("figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)
CSV_CLEANED      = DATA_DIR / "Final_table_cleaned.csv"
COHORT_OUT       = DATA_DIR / "abandonment_cohort_v3.csv"
TABLE_STRAT_OUT  = DATA_DIR / "abandonment_risk_table_v3.csv"
TABLE_MONTH_OUT  = DATA_DIR / "abandonment_risk_table_txmonth_v3.csv"
CENSOR_DATE      = pd.Timestamp("2024-12-31")
ABANDON_OUTCOMES = {"Abandono", "Abandono Primario"}
OBITO_OUTCOMES   = {"Obito TB", "Obito NTB"}
FOLLOW_YEARS     = [1, 2, 3]
MAX_PLOT_YEARS   = 5

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("Loading cleaned dataset…")
df = pd.read_csv(CSV_CLEANED, low_memory=False)
print(f"  {len(df):,} rows  |  {df['sinan_clean'].nunique():,} individuals")

for col in ["notification_date", "end_date", "dod", "tx_start"]:
    df[col] = pd.to_datetime(df[col], errors="coerce")

df = df.sort_values(["sinan_clean", "notification_date"])

# ── 2. Comprehensive death date per individual ────────────────────────────────
# Source A: dod field (min over all rows)
dod_src_a = (
    df.dropna(subset=["dod"])
    .groupby("sinan_clean")["dod"].min()
    .reset_index().rename(columns={"dod": "death_src_a"})
)

# Source B: Obito case_outcome → use end_date of that row
obito_rows = df[df["case_outcome"].isin(OBITO_OUTCOMES)].dropna(subset=["end_date"])
dod_src_b  = (
    obito_rows.groupby("sinan_clean")["end_date"].min()
    .reset_index().rename(columns={"end_date": "death_src_b"})
)

# Merge and take earliest
dod_combined = dod_src_a.merge(dod_src_b, on="sinan_clean", how="outer")
dod_combined["death_date_all"] = dod_combined[["death_src_a","death_src_b"]].min(axis=1)
print(f"\nDeath ascertainment:")
print(f"  Individuals with dod field     : {len(dod_src_a):,}")
print(f"  Individuals with Obito outcome : {len(dod_src_b):,}")
print(f"  Combined (either)              : {len(dod_combined):,}")

# ── 3. Build abandonment index episodes ──────────────────────────────────────
covariate_cols = ["hiv", "aids", "address_type", "age_tb", "sex",
                  "dob", "race", "tx_start"]

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
})
aband_raw = aband_raw[aband_raw["abandon_end"] <= CENSOR_DATE].copy()
print(f"\nAbandonment cohort: {len(aband_raw):,}")

# ── 4. Next notification after abandonment (for re-notification outcome) ─────
sub_events = df.merge(aband_raw[["sinan_clean","abandon_end"]], on="sinan_clean")
sub_events = sub_events[sub_events["notification_date"] > sub_events["abandon_end"]]
next_notif = (
    sub_events.groupby("sinan_clean")["notification_date"].min()
    .reset_index().rename(columns={"notification_date": "next_notif_date"})
)
aband = aband_raw.merge(next_notif, on="sinan_clean", how="left")

# Cap re-notification at censor date
aband.loc[aband["next_notif_date"] > CENSOR_DATE, "next_notif_date"] = pd.NaT

# ── 5. Death: comprehensive, AFTER abandonment, capped at censor date ────────
aband = aband.merge(dod_combined[["sinan_clean","death_date_all"]],
                    on="sinan_clean", how="left")

# Only count death strictly after the abandonment end_date
aband["death_date"] = pd.to_datetime(
    np.where(aband["death_date_all"] > aband["abandon_end"],
             aband["death_date_all"], pd.NaT)
)
aband.loc[aband["death_date"] > CENSOR_DATE, "death_date"] = pd.NaT

# ── 6. Re-notification outcome (competing risks: death censors if before notif)
# event_rn: 0=censored, 1=re-notified, 2=died before re-notified (competing)
def assign_renotif(row):
    t_n = row["next_notif_date"]
    t_d = row["death_date"]
    has_n, has_d = pd.notna(t_n), pd.notna(t_d)
    if has_n and has_d:
        if t_n <= t_d: return (1, t_n)
        else:          return (2, t_d)   # died first → competing event
    elif has_n: return (1, t_n)
    elif has_d: return (2, t_d)
    else:       return (0, CENSOR_DATE)

rn_out = aband.apply(assign_renotif, axis=1, result_type="expand")
rn_out.columns = ["event_rn", "event_date_rn"]
aband = pd.concat([aband, rn_out], axis=1)
aband["event_date_rn"] = pd.to_datetime(aband["event_date_rn"])
aband["time_rn"] = (aband["event_date_rn"] - aband["abandon_end"]).dt.days / 365.25

# ── 7. Death outcome: independent of re-notification ─────────────────────────
# event_d: 0=alive/censored, 1=died
aband["event_d"] = aband["death_date"].notna().astype(int)
aband["event_date_d"] = aband["death_date"].where(
    aband["death_date"].notna(), CENSOR_DATE
)
aband["event_date_d"] = pd.to_datetime(aband["event_date_d"])
aband["time_d"] = (aband["event_date_d"] - aband["abandon_end"]).dt.days / 365.25

# Drop implausible rows (negative times in either outcome)
bad = (aband["time_rn"] <= 0) | (aband["time_d"] <= 0)
print(f"Dropping {bad.sum()} rows with non-positive time")
aband = aband[~bad].copy()

print(f"\nFinal cohort: {len(aband):,}")
print(f"RE-NOTIFICATION:")
print(f"  Censored (0)          : {(aband['event_rn']==0).sum():,}")
print(f"  Re-notified (1)       : {(aband['event_rn']==1).sum():,}")
print(f"  Died before notif (2) : {(aband['event_rn']==2).sum():,}")
print(f"DEATH (comprehensive, independent):")
print(f"  Alive/censored (0)    : {(aband['event_d']==0).sum():,}")
print(f"  Died (1)              : {(aband['event_d']==1).sum():,}")

# ── 8. Stratification variables ───────────────────────────────────────────────
def age_group(a):
    if pd.isna(a): return "Unknown"
    elif a < 18:   return "<18"
    elif a < 35:   return "18–34"
    elif a < 50:   return "35–49"
    elif a < 65:   return "50–64"
    else:          return "65+"

def hiv_group(row):
    if row["aids"] == "S":   return "AIDS"
    if row["hiv"] == "Pos":  return "HIV+"
    if row["hiv"] == "Neg":  return "HIV-"
    return "Unknown/not tested"

aband["age_group"]   = aband["age_tb"].apply(age_group)
aband["sex_grp"]     = aband["sex"].map({"M":"Male","F":"Female"}).fillna("Unknown")
aband["hiv_grp"]     = aband.apply(hiv_group, axis=1)
aband["incarcerated"]= aband["address_type"].apply(
    lambda x: "Incarcerated" if str(x).strip().upper()=="DETENTO" else "Not incarcerated")
aband["homeless"]    = aband["address_type"].apply(
    lambda x: "Homeless" if "SEM RESIDENCIA" in str(x).strip().upper() else "Not homeless")

# Treatment month
aband["tx_days"] = (aband["abandon_end"] - aband["tx_start"]).dt.days
aband_tx = aband[aband["tx_days"] > 0].copy()
aband_tx["tx_month_raw"] = np.ceil(aband_tx["tx_days"] / 30.44).astype(int)
aband_tx["tx_month_grp"] = aband_tx["tx_month_raw"].apply(
    lambda m: f"Month {m}" if m <= 6 else "Month 7+")

aband.to_csv(COHORT_OUT, index=False)
print(f"\nCohort saved → {COHORT_OUT}")

# ── 9. Estimators ─────────────────────────────────────────────────────────────
def aalen_johansen(times, events, cause):
    """CIF for re-notification (with death as competing event)."""
    order  = np.argsort(times)
    ts, es = times[order], events[order]
    S, CI, n = 1.0, 0.0, len(ts)
    pts = [(0.0, 0.0)]
    prev = -np.inf
    for t in ts:
        if t == prev: continue
        mask  = ts == t
        d_all = np.sum(es[mask] > 0)
        d_c   = np.sum(es[mask] == cause)
        if n > 0 and d_all > 0:
            CI += S * (d_c / n)
            S  *= 1 - d_all / n
        n -= mask.sum(); prev = t
        pts.append((t, CI))
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

def kaplan_meier(times, events):
    """1-KM for death (independent; re-notification does not censor)."""
    order  = np.argsort(times)
    ts, es = times[order], events[order]
    S, n   = 1.0, len(ts)
    pts = [(0.0, 0.0)]
    prev = -np.inf
    for t in ts:
        if t == prev: continue
        mask = ts == t
        d    = np.sum(es[mask] == 1)
        if n > 0 and d > 0:
            S *= 1 - d / n
        n -= mask.sum(); prev = t
        pts.append((t, 1 - S))         # return 1-KM (= cumulative incidence)
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

def interp_val(t_arr, v_arr, year):
    idx = np.searchsorted(t_arr, year, side="right") - 1
    return float(v_arr[max(0, min(idx, len(v_arr)-1))])

# ── 10. Plot helpers ──────────────────────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")

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

MONTH_ORDER = ["Month 1","Month 2","Month 3","Month 4","Month 5","Month 6","Month 7+"]
MONTH_PAL   = ["#f7dc6f","#f0a500","#e05c5c","#b03060","#7b2d8b","#3d5fa0","#1a3a5c"]

def make_cif_fig(data, strat_col, outcome, palette, order, label, suffix=""):
    """Draw stratified CIF (re-notification) figure."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    rows = []
    groups = [g for g in order if g in data[strat_col].unique()]
    for grp, color in zip(groups, palette):
        sub = data[data[strat_col] == grp]
        if len(sub) < 20: continue
        t_arr, ci_arr = aalen_johansen(sub["time_rn"].values, sub["event_rn"].values, 1)
        ax.step(t_arr, ci_arr*100, where="post", lw=2, color=color,
                label=f"{grp}  (n={len(sub):,})")
        row = {"subgroup":grp, "n":len(sub)}
        for yr in FOLLOW_YEARS:
            row[f"risk_{yr}yr_%"] = round(interp_val(t_arr, ci_arr, yr)*100, 1)
        rows.append(row)
    for yr in FOLLOW_YEARS: ax.axvline(yr, color="grey", lw=0.7, ls=":")
    ax.set_xlim(0, MAX_PLOT_YEARS); ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=11)
    ax.set_ylabel("Cumulative Incidence of Re-notification (%)", fontsize=11)
    ax.set_title(f"Risk of Re-notification by {label}\n"
                 f"After TB Treatment Abandonment  (N={len(data):,}; censored Dec 31, 2024)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.85)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    plt.tight_layout()
    fname = FIG_DIR / f"cif_renotification_by_{strat_col}{suffix}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  Saved: {fname}")
    return rows

def make_death_fig(data, strat_col, palette, order, label, suffix=""):
    """Draw stratified 1-KM (death) figure."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    rows = []
    groups = [g for g in order if g in data[strat_col].unique()]
    for grp, color in zip(groups, palette):
        sub = data[data[strat_col] == grp]
        if len(sub) < 20: continue
        t_arr, ci_arr = kaplan_meier(sub["time_d"].values, sub["event_d"].values)
        ax.step(t_arr, ci_arr*100, where="post", lw=2, color=color,
                label=f"{grp}  (n={len(sub):,})")
        row = {"subgroup":grp, "n":len(sub)}
        for yr in FOLLOW_YEARS:
            row[f"risk_{yr}yr_%"] = round(interp_val(t_arr, ci_arr, yr)*100, 1)
        rows.append(row)
    for yr in FOLLOW_YEARS: ax.axvline(yr, color="grey", lw=0.7, ls=":")
    ax.set_xlim(0, MAX_PLOT_YEARS); ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=11)
    ax.set_ylabel("Cumulative Incidence of Death (%)", fontsize=11)
    ax.set_title(f"Risk of Death by {label}\n"
                 f"After TB Treatment Abandonment  (N={len(data):,}; censored Dec 31, 2024)\n"
                 f"[Death from dod field OR Obito case outcome; re-notification does not censor]",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.85)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    plt.tight_layout()
    fname = FIG_DIR / f"cif_death_by_{strat_col}{suffix}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  Saved: {fname}")
    return rows

# ── 11. Overall figures ───────────────────────────────────────────────────────
print("\n── Overall figures ──────────────────────────────────────────────────────")

def overall_fig(data, outcome_label, t_col, e_col, estimator_fn,
                color, y_label, fname):
    fig, ax = plt.subplots(figsize=(8, 5))
    t_arr, ci_arr = estimator_fn(data[t_col].values, data[e_col].values)
    ax.step(t_arr, ci_arr*100, where="post", lw=2.5, color=color, label="Overall")
    ax.fill_between(t_arr, ci_arr*100, step="post", alpha=0.12, color=color)
    for yr in FOLLOW_YEARS:
        val = interp_val(t_arr, ci_arr, yr)*100
        ax.axvline(yr, color="grey", lw=0.7, ls=":")
        ax.annotate(f"{val:.1f}%", xy=(yr, val), xytext=(yr+0.08, val+0.25),
                    fontsize=10, fontweight="bold", color=color)
    ax.set_xlim(0, MAX_PLOT_YEARS); ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f"Risk of {outcome_label} After TB Treatment Abandonment\n"
                 f"(N={len(data):,}; censored Dec 31, 2024)", fontsize=13, fontweight="bold")
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    plt.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  Saved: {FIG_DIR / fname}")
    return t_arr, ci_arr

rn_t_ov, rn_ci_ov = overall_fig(
    aband, "Re-notification", "time_rn", "event_rn",
    lambda t, e: aalen_johansen(t, e, 1),
    "#E05C5C", "Cumulative Incidence (%)",
    "cif_renotification_overall_v2.png"
)
d_t_ov, d_ci_ov = overall_fig(
    aband, "Death (comprehensive)", "time_d", "event_d",
    kaplan_meier,
    "#5C7AE0", "Cumulative Incidence of Death (%)",
    "cif_death_overall_v2.png"
)

# ── 12. Stratified figures & risk table ──────────────────────────────────────
print("\n── Stratified figures (demographic/clinical) ────────────────────────────")
strata = ["age_group","sex_grp","hiv_grp","incarcerated","homeless"]
all_risk_rows = []

# Overall row
for outcome, t_arr, ci_arr in [("Re-notification", rn_t_ov, rn_ci_ov),
                                 ("Death",           d_t_ov,  d_ci_ov)]:
    row = {"outcome":outcome,"stratum":"Overall","variable":"overall","subgroup":"All","n":len(aband)}
    for yr in FOLLOW_YEARS:
        row[f"risk_{yr}yr_%"] = round(interp_val(t_arr, ci_arr, yr)*100, 1)
    all_risk_rows.append(row)

for sc in strata:
    label   = STRATUM_LABELS[sc]
    palette = PALETTES[sc]
    order   = ORDER[sc]

    rn_rows = make_cif_fig(aband, sc, "rn", palette, order, label, "_v2")
    d_rows  = make_death_fig(aband, sc, palette, order, label, "_v2")

    for row in rn_rows:
        row.update({"outcome":"Re-notification","stratum":label,"variable":sc})
    for row in d_rows:
        row.update({"outcome":"Death","stratum":label,"variable":sc})
    all_risk_rows.extend(rn_rows + d_rows)

# Save stratified table
risk_table = pd.DataFrame(all_risk_rows)[
    ["outcome","stratum","subgroup","n","risk_1yr_%","risk_2yr_%","risk_3yr_%"]]
risk_table.to_csv(TABLE_STRAT_OUT, index=False)
print(f"\nStratified risk table → {TABLE_STRAT_OUT}")

print("\n── Risk Table ───────────────────────────────────────────────────────────")
for outcome in ["Re-notification","Death"]:
    print(f"\n  {outcome}:")
    print(risk_table[risk_table["outcome"]==outcome].drop(columns="outcome").to_string(index=False))

# ── 13. Treatment month analysis ─────────────────────────────────────────────
print("\n── Treatment month figures ──────────────────────────────────────────────")

# Distribution + 1yr bar overview
grp_counts = (aband_tx.groupby("tx_month_grp").size()
              .reindex(MONTH_ORDER).fillna(0).astype(int))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].bar(grp_counts.index, grp_counts.values,
            color=MONTH_PAL, edgecolor="white")
for i,(g,n) in enumerate(grp_counts.items()):
    axes[0].text(i, n+30, f"{n:,}", ha="center", va="bottom", fontsize=9)
axes[0].set_xlabel("Month of treatment at abandonment", fontsize=11)
axes[0].set_ylabel("Number of individuals", fontsize=11)
axes[0].set_title("Distribution of Abandonment\nby Treatment Month", fontsize=12, fontweight="bold")
axes[0].tick_params(axis="x", rotation=25)

rn_1yr, d_1yr = [], []
for grp in MONTH_ORDER:
    sub = aband_tx[aband_tx["tx_month_grp"]==grp]
    if len(sub) < 20:
        rn_1yr.append(np.nan); d_1yr.append(np.nan); continue
    t1,c1 = aalen_johansen(sub["time_rn"].values, sub["event_rn"].values, 1)
    t2,c2 = kaplan_meier(sub["time_d"].values, sub["event_d"].values)
    rn_1yr.append(interp_val(t1,c1,1)*100)
    d_1yr.append(interp_val(t2,c2,1)*100)

x = np.arange(len(MONTH_ORDER)); w = 0.38
ax2 = axes[1]
bars1 = ax2.bar(x-w/2, rn_1yr, w, color="#E05C5C", label="Re-notification", edgecolor="white")
bars2 = ax2.bar(x+w/2, d_1yr,  w, color="#5C7AE0", label="Death",           edgecolor="white")
for bar, val in zip(bars1, rn_1yr):
    if not np.isnan(val):
        ax2.text(bar.get_x()+bar.get_width()/2, val+0.3,
                 f"{val:.1f}%", ha="center", va="bottom", fontsize=7.5, color="#a00")
for bar, val in zip(bars2, d_1yr):
    if not np.isnan(val):
        ax2.text(bar.get_x()+bar.get_width()/2, val+0.05,
                 f"{val:.1f}%", ha="center", va="bottom", fontsize=7, color="#334")
ax2.set_xticks(x); ax2.set_xticklabels(MONTH_ORDER, rotation=25)
ax2.set_ylabel("1-year Cumulative Incidence (%)", fontsize=11)
ax2.set_title("1-year Risk by Month of Abandonment\n[Comprehensive death ascertainment]",
              fontsize=12, fontweight="bold")
ax2.legend(fontsize=10)
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))

plt.suptitle(f"Treatment Month at Abandonment — Overview  (N={len(aband_tx):,}; censored Dec 31, 2024)",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR/"tx_month_distribution_v2.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIG_DIR/'tx_month_distribution_v2.png'}")

# CIF curves by tx month
month_risk_rows = []
for cause_type, label_txt, y_label, estimator, event_col, time_col in [
    ("renotification", "Re-notification", "Cumulative Incidence of Re-notification (%)",
     lambda sub: aalen_johansen(sub["time_rn"].values, sub["event_rn"].values, 1),
     "event_rn", "time_rn"),
    ("death", "Death (comprehensive)", "Cumulative Incidence of Death (%)",
     lambda sub: kaplan_meier(sub["time_d"].values, sub["event_d"].values),
     "event_d", "time_d"),
]:
    fig, ax = plt.subplots(figsize=(10, 6))
    for grp, color in zip(MONTH_ORDER, MONTH_PAL):
        sub = aband_tx[aband_tx["tx_month_grp"]==grp]
        if len(sub) < 20: continue
        t_arr, ci_arr = estimator(sub)
        ax.step(t_arr, ci_arr*100, where="post", lw=2, color=color,
                label=f"{grp}  (n={len(sub):,})")
        row = {"outcome":label_txt, "tx_month_group":grp, "n":len(sub)}
        for yr in FOLLOW_YEARS:
            row[f"risk_{yr}yr_%"] = round(interp_val(t_arr, ci_arr, yr)*100, 1)
        month_risk_rows.append(row)

    for yr in FOLLOW_YEARS: ax.axvline(yr, color="grey", lw=0.7, ls=":")
    ax.set_xlim(0, MAX_PLOT_YEARS); ax.set_ylim(0)
    ax.set_xlabel("Years since abandonment", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    death_note = ("\n[dod field + Obito case outcome; re-notification does not censor]"
                  if cause_type=="death" else "")
    ax.set_title(f"Risk of {label_txt} by Treatment Month at Abandonment{death_note}\n"
                 f"(N={len(aband_tx):,}; censored Dec 31, 2024)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9.5, loc="upper left", framealpha=0.9,
              title="Month of TX at abandonment", title_fontsize=9)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.0f%%"))
    plt.tight_layout()
    fname = FIG_DIR / f"cif_{cause_type}_by_tx_month_v2.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  Saved: {fname}")

# Save month table
month_table = pd.DataFrame(month_risk_rows)[
    ["outcome","tx_month_group","n","risk_1yr_%","risk_2yr_%","risk_3yr_%"]]
month_table.to_csv(TABLE_MONTH_OUT, index=False)
print(f"\nMonth risk table → {TABLE_MONTH_OUT}")

print("\n── Treatment Month Risk Table ───────────────────────────────────────────")
for outcome in ["Re-notification","Death (comprehensive)"]:
    print(f"\n  {outcome}:")
    print(month_table[month_table["outcome"]==outcome].drop(columns="outcome").to_string(index=False))

print("\nDone ✓")
