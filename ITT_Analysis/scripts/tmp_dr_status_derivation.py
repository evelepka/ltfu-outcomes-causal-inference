"""
tmp_dr_status_derivation.py
---------------------------
Exploratory derivation of a per-patient drug-resistance status (`dr_status`)
by rolling up four notification-level variables across all rows for each
sinan_clean in Data/Final_table_cleaned.csv:

  - tmr_tb     (Xpert MTB/RIF)
  - resistance (DST overall: SENS / TB R / TB MR / AND)
  - rifasens   (RIF DST: Sens / Resist / ...)
  - isonisens  (INH DST: Sens / Resist / ...)

Hierarchical coding rule (first match wins, applied per sinan_clean across ALL
their rows):

  1. RR/MDR-TB:
       tmr_tb == "Mtb detectado - Rifamp resistente"  OR
       rifasens == "Resist"                           OR
       resistance == "TB MR"
  2. INH-mono resistance: NOT in (1), AND isonisens == "Resist" somewhere.
  3. Sensitive: NOT in (1) or (2), AND any positive evidence of sensitivity:
       tmr_tb == "Mtb detectado - Rifamp sensivel"          OR
       resistance == "SENS"                                  OR
       (rifasens == "Sens" AND isonisens == "Sens")
  4. Not Evaluated: everyone else, including resistance == "TB R" with no
     specific drug-level information.

Output:
  - Lookup CSV at ITT_Analysis/results/resistance/dr_status_lookup.csv
    (sinan_clean, dr_status — one row per unique sinan_clean)
  - Summary report at ITT_Analysis/results/resistance/dr_status_summary.txt

NOTE: This is a side derivation. It does NOT modify itt_cohort.csv and is NOT
a new cohort file (CLAUDE.md prohibits cohort proliferation). Filter at
read-time as needed.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------
# Project root resolution (mirrors 01_itt_cohort_selection.py)
# ----------------------------------------------------------------------------
def _find_project_root():
    if os.environ.get("TB_ABANDONMENT_ROOT"):
        p = Path(os.environ["TB_ABANDONMENT_ROOT"]).expanduser()
        if p.exists():
            return p
    candidates = [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]
    for c in candidates:
        if c.exists():
            return c
    return Path(__file__).resolve().parents[2]


BASE_DIR = _find_project_root()
FINAL_DATA = BASE_DIR / "Data" / "Final_table_cleaned.csv"
COHORT_CSV = BASE_DIR / "ITT_Analysis" / "data" / "itt_cohort.csv"
OUT_DIR = BASE_DIR / "ITT_Analysis" / "results" / "resistance"
OUT_LOOKUP = OUT_DIR / "dr_status_lookup.csv"
OUT_SUMMARY = OUT_DIR / "dr_status_summary.txt"

OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Base Dir: {BASE_DIR}")
if not FINAL_DATA.exists():
    raise FileNotFoundError(
        f"Final_table_cleaned.csv not found at {FINAL_DATA}. "
        f"Set TB_ABANDONMENT_ROOT or mount Google Drive."
    )

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
print("Loading cleaned dataset...")
USECOLS = ["sinan_clean", "tmr_tb", "resistance", "rifasens", "isonisens"]
df = pd.read_csv(FINAL_DATA, usecols=USECOLS, low_memory=False)
print(f"  rows: {len(df):,}   unique sinan_clean: {df['sinan_clean'].nunique():,}")


# Normalize whitespace on the categorical strings (defensive — values may have
# trailing spaces in SINAN exports).
for col in ["tmr_tb", "resistance", "rifasens", "isonisens"]:
    df[col] = df[col].astype("string").str.strip()


# ----------------------------------------------------------------------------
# Per-row evidence flags
# ----------------------------------------------------------------------------
TMR_RIF_R = "Mtb detectado - Rifamp resistente"
TMR_RIF_S = "Mtb detectado - Rifamp sensivel"

def _eq(series: pd.Series, value: str) -> pd.Series:
    """Equality comparison that returns plain bool (NA -> False)."""
    return (series == value).fillna(False).astype(bool)


df["ev_rr_mdr"] = (
    _eq(df["tmr_tb"], TMR_RIF_R)
    | _eq(df["rifasens"], "Resist")
    | _eq(df["resistance"], "TB MR")
)
df["ev_inh_resist"] = _eq(df["isonisens"], "Resist")
df["ev_tmr_sens"] = _eq(df["tmr_tb"], TMR_RIF_S)
df["ev_dst_sens"] = _eq(df["resistance"], "SENS")
df["ev_rif_sens"] = _eq(df["rifasens"], "Sens")
df["ev_inh_sens"] = _eq(df["isonisens"], "Sens")
df["ev_tbR_unspec"] = _eq(df["resistance"], "TB R")


# ----------------------------------------------------------------------------
# Roll up to per-patient
# ----------------------------------------------------------------------------
print("Rolling up evidence per sinan_clean...")
ev_cols = [
    "ev_rr_mdr", "ev_inh_resist", "ev_tmr_sens", "ev_dst_sens",
    "ev_rif_sens", "ev_inh_sens", "ev_tbR_unspec",
]
pt = df.groupby("sinan_clean", as_index=False)[ev_cols].max()  # any-True


# ----------------------------------------------------------------------------
# Apply hierarchy
# ----------------------------------------------------------------------------
# Vectorized hierarchical assignment
pt["dr_status"] = "Not Evaluated"
sens_mask = (
    pt["ev_tmr_sens"]
    | pt["ev_dst_sens"]
    | (pt["ev_rif_sens"] & pt["ev_inh_sens"])
)
pt.loc[sens_mask, "dr_status"] = "Sensitive"
pt.loc[pt["ev_inh_resist"], "dr_status"] = "INH-mono resistance"
pt.loc[pt["ev_rr_mdr"], "dr_status"] = "RR/MDR-TB"


# Anomaly: did any patient end up with both Sens AND Resist evidence across
# rows? (e.g., Sens RIF on one notification, Resist RIF on another — a real
# data-quality concern.)
# Build "Resist-anywhere" flags per patient and align with pt's sinan_clean.
def _any_per_patient(series_eq: pd.Series) -> pd.Series:
    return series_eq.groupby(df["sinan_clean"]).any()


rif_resist_any = _any_per_patient(_eq(df["rifasens"], "Resist")).reindex(pt["sinan_clean"]).values
inh_resist_any = _any_per_patient(_eq(df["isonisens"], "Resist")).reindex(pt["sinan_clean"]).values
tmr_resist_any = _any_per_patient(_eq(df["tmr_tb"], TMR_RIF_R)).reindex(pt["sinan_clean"]).values
dst_resist_any = _any_per_patient(
    _eq(df["resistance"], "TB R") | _eq(df["resistance"], "TB MR")
).reindex(pt["sinan_clean"]).values

pt["anom_rif_both"] = pt["ev_rif_sens"].values & rif_resist_any
pt["anom_inh_both"] = pt["ev_inh_sens"].values & inh_resist_any
pt["anom_tmr_both"] = pt["ev_tmr_sens"].values & tmr_resist_any
pt["anom_dst_both"] = pt["ev_dst_sens"].values & dst_resist_any


# ----------------------------------------------------------------------------
# Write lookup
# ----------------------------------------------------------------------------
lookup = pt[["sinan_clean", "dr_status"]].copy()
lookup.to_csv(OUT_LOOKUP, index=False)
print(f"Wrote lookup: {OUT_LOOKUP}  ({len(lookup):,} rows)")


# ----------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------
DR_ORDER = ["RR/MDR-TB", "INH-mono resistance", "Sensitive", "Not Evaluated"]

overall_counts = (
    pt["dr_status"]
    .value_counts()
    .reindex(DR_ORDER)
    .fillna(0)
    .astype(int)
)

# Of those landing in Not Evaluated, how many had `resistance == "TB R"`
# somewhere (the deliberately-punted subgroup)?
ne_mask = pt["dr_status"] == "Not Evaluated"
tbR_punted = int((ne_mask & pt["ev_tbR_unspec"]).sum())
tbR_total = int(pt["ev_tbR_unspec"].sum())

# Anomaly counts (mixed evidence across rows for same patient)
anom_rif = int(pt["anom_rif_both"].sum())
anom_inh = int(pt["anom_inh_both"].sum())
anom_tmr = int(pt["anom_tmr_both"].sum())
anom_dst = int(pt["anom_dst_both"].sum())
anom_any = int(
    (pt["anom_rif_both"] | pt["anom_inh_both"] | pt["anom_tmr_both"] | pt["anom_dst_both"]).sum()
)


# ----------------------------------------------------------------------------
# Cross-tab vs ITT cohort
# ----------------------------------------------------------------------------
print("Loading ITT cohort for cross-tabs...")
cohort = pd.read_csv(
    COHORT_CSV,
    usecols=["sinan_clean", "resistance_clean", "itt_group"],
    low_memory=False,
)
print(f"  cohort rows: {len(cohort):,}")

merged = cohort.merge(lookup, on="sinan_clean", how="left")
n_unmatched = int(merged["dr_status"].isna().sum())
merged["dr_status"] = merged["dr_status"].fillna("Not Evaluated")

cohort_counts = (
    merged["dr_status"]
    .value_counts()
    .reindex(DR_ORDER)
    .fillna(0)
    .astype(int)
)

xt_resist = pd.crosstab(merged["dr_status"], merged["resistance_clean"], dropna=False).reindex(DR_ORDER)
xt_group = pd.crosstab(merged["dr_status"], merged["itt_group"], dropna=False).reindex(DR_ORDER)


# ----------------------------------------------------------------------------
# Print + write report
# ----------------------------------------------------------------------------
def fmt_section(title: str) -> str:
    return "\n" + "=" * 78 + f"\n  {title}\n" + "=" * 78 + "\n"


lines = []
lines.append("DR-status derivation summary")
lines.append(f"Source:  {FINAL_DATA}")
lines.append(f"Cohort:  {COHORT_CSV}")
lines.append(f"Output:  {OUT_LOOKUP}")
lines.append("")
lines.append(f"Unique sinan_clean in raw file: {df['sinan_clean'].nunique():,}")
lines.append(f"Cohort rows (LTFU + Non-LTFU):  {len(cohort):,}")
lines.append(f"Cohort rows unmatched in lookup (set to Not Evaluated): {n_unmatched:,}")

lines.append(fmt_section("Overall counts (raw file, all sinan_clean)"))
lines.append(overall_counts.to_string())
lines.append("")
lines.append(f"Total: {int(overall_counts.sum()):,}")

lines.append(fmt_section("Counts within ITT cohort"))
lines.append(cohort_counts.to_string())
lines.append("")
lines.append(f"Total: {int(cohort_counts.sum()):,}")

lines.append(fmt_section("Cross-tab: dr_status x resistance_clean (ITT cohort)"))
lines.append(xt_resist.to_string())

lines.append(fmt_section("Cross-tab: dr_status x itt_group (ITT cohort)"))
lines.append(xt_group.to_string())

lines.append(fmt_section("Unresolved 'TB R' (resistance==TB R, no rifa/isoni Resist info)"))
lines.append(f"Patients with resistance=='TB R' anywhere:                  {tbR_total:,}")
lines.append(f"  Of those, classified as 'Not Evaluated' (deliberate punt): {tbR_punted:,}")
lines.append(f"  Of those, captured by hierarchy (RR/MDR or INH-mono):      {tbR_total - tbR_punted:,}")

lines.append(fmt_section("Anomalies: mixed Sens+Resist evidence across rows for same patient"))
lines.append(f"Patients with rifasens both Sens AND Resist:              {anom_rif:,}")
lines.append(f"Patients with isonisens both Sens AND Resist:             {anom_inh:,}")
lines.append(f"Patients with tmr_tb RIF-sens AND RIF-resist rows:        {anom_tmr:,}")
lines.append(f"Patients with resistance SENS AND (TB R or TB MR) rows:   {anom_dst:,}")
lines.append(f"Patients with ANY of the above mixed-evidence patterns:   {anom_any:,}")
lines.append("")
lines.append("(Note: hierarchy treats any RIF/INH 'Resist' as overriding 'Sens' on")
lines.append(" another row, so these patients are classified as resistant.)")

report = "\n".join(lines)
print(report)

with open(OUT_SUMMARY, "w") as fh:
    fh.write(report + "\n")
print(f"\nWrote summary: {OUT_SUMMARY}")
