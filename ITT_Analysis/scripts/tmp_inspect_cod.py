"""Inspect cause_of_death_code in the cohort."""
import os
from pathlib import Path
import pandas as pd


def _root():
    for c in [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]:
        if c.exists():
            return c
    raise FileNotFoundError("project root not found")


BASE = _root()
COHORT = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"
RAW = BASE / "Data" / "Final_table_cleaned.csv"

cohort = pd.read_csv(COHORT, low_memory=False, parse_dates=["death_date"])
print(f"Cohort N: {len(cohort):,}; deaths in cohort: {(cohort['event_d']==1).sum():,}")

# Pull cause_of_death_code from raw, join on sinan_clean using max date_of_death
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=["sinan_clean", "dod", "cause_of_death_code", "case_outcome"])
raw["dod"] = pd.to_datetime(raw["dod"], errors="coerce")

# Some individuals have multiple notifications; pick the cause-of-death record
# with the latest dod (the actual death record)
death_recs = raw.dropna(subset=["dod", "cause_of_death_code"]).copy()
death_recs = death_recs.sort_values("dod").groupby("sinan_clean", as_index=False).last()
print(f"Raw with cause_of_death_code: {len(death_recs):,}")

m = cohort.merge(
    death_recs[["sinan_clean", "cause_of_death_code"]],
    on="sinan_clean", how="left"
)
print(f"Cohort with COD code linked: {m['cause_of_death_code'].notna().sum():,}")

# Of cohort deaths (event_d=1), how many have COD code?
deaths = m[m["event_d"] == 1]
print(f"Cohort deaths (event_d=1): {len(deaths):,}")
print(f"  ... with COD code:        {deaths['cause_of_death_code'].notna().sum():,}")
print(f"  ... without COD code:     {deaths['cause_of_death_code'].isna().sum():,}")

# Look at format and top codes
print("\nTop 30 cause-of-death codes among cohort deaths:")
print(deaths["cause_of_death_code"].value_counts(dropna=False).head(30))

# Count how many start with each ICD-10 letter, focused on TB-relevant
print("\nDistribution of first letter of COD code (cohort deaths):")
deaths_codes = deaths["cause_of_death_code"].dropna().astype(str).str.upper().str.strip()
print(deaths_codes.str[0].value_counts().head(20))

# Specifically TB block: A15-A19, B90, B20
print("\nTB-block prefix tallies:")
for prefix in ["A15", "A16", "A17", "A18", "A19", "B90", "B20", "B21", "B22"]:
    n = deaths_codes.str.startswith(prefix).sum()
    print(f"  starts with {prefix}: {n}")

# Respiratory J block
print("\nRespiratory J prefix:")
for prefix in ["J0", "J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8", "J9"]:
    n = deaths_codes.str.startswith(prefix).sum()
    print(f"  starts with J{prefix[1]}: {n}")
