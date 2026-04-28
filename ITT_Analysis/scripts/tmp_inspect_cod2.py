"""
Build a hybrid TB-attributable cause-of-death indicator combining:
  (1) SIM ICD-10 cause_of_death_code  (where available)
  (2) TBweb case_outcome 'Obito TB' vs 'Obito NTB'  (for the index episode
      itself; relevant only for Non-LTFU arm where the case closed as Obito)

Check missingness by itt_group, and tabulate the resulting hybrid variable.
"""
import os
from pathlib import Path
import pandas as pd
import numpy as np


def _root():
    for c in [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]:
        if c.exists():
            return c
    raise FileNotFoundError


BASE = _root()
COHORT = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"
RAW = BASE / "Data" / "Final_table_cleaned.csv"

cohort = pd.read_csv(COHORT, low_memory=False, parse_dates=["death_date", "best_start", "end_date"])

# Recover index case_outcome (same as my earlier sanity script)
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=["sinan_clean", "case_type", "case_outcome", "end_date",
                           "dod", "cause_of_death_code"])
raw["end_date"] = pd.to_datetime(raw["end_date"], errors="coerce")
raw["dod"] = pd.to_datetime(raw["dod"], errors="coerce")

novo = raw[raw["case_type"].str.strip().str.title() == "Novo"]
novo = novo[novo["case_outcome"].notna()
            & (novo["case_outcome"].str.strip() != "")
            & (novo["case_outcome"] != "Mud Diag")]
TRANSFER = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
novo = novo[~novo["case_outcome"].isin(TRANSFER)]
first = novo.sort_values("end_date").groupby("sinan_clean", as_index=False).first()
first = first[["sinan_clean", "case_outcome"]]
cohort = cohort.merge(first, on="sinan_clean", how="left")

# Pull SIM cause_of_death_code (last record per sinan with both dod and code)
death_recs = raw.dropna(subset=["dod", "cause_of_death_code"]).copy()
death_recs = death_recs.sort_values("dod").groupby("sinan_clean", as_index=False).last()
cohort = cohort.merge(
    death_recs[["sinan_clean", "cause_of_death_code"]],
    on="sinan_clean", how="left"
)

# ----------------------------------------------------------------------------
# Define hybrid TB-attribution
# ----------------------------------------------------------------------------
# TB code prefixes (ICD-10): A15-A19, B90, plus B200 (HIV w/ mycobacterial).
TB_PREFIXES = ("A15", "A16", "A17", "A18", "A19", "B90")
RESP_PREFIXES = ("J",)  # J00-J99
HIV_PREFIXES = ("B20", "B21", "B22", "B23", "B24")

cod = cohort["cause_of_death_code"].astype(str).str.upper().str.strip()
cohort["cod_tb_strict"] = cod.str.startswith(TB_PREFIXES)
cohort["cod_resp"]      = cod.str.startswith(RESP_PREFIXES)
cohort["cod_hiv"]       = cod.str.startswith(HIV_PREFIXES)
# B20.0 specifically (HIV with mycobacterial infection)
cohort["cod_hiv_tb"]    = cod.str.startswith("B200")
cohort["cod_known"]     = cohort["cause_of_death_code"].notna() & (cod != "NAN") & (cod != "")

# TBweb case_outcome attribution
cohort["tbweb_obito_tb"]  = cohort["case_outcome"] == "Obito TB"
cohort["tbweb_obito_ntb"] = cohort["case_outcome"] == "Obito NTB"

# Hybrid: a death is "TB-attributable" if SIM ICD code is in TB-block,
# OR (no SIM code) but TBweb closed the index episode as Obito TB.
# A death is "non-TB" if SIM ICD is non-TB-non-respiratory-non-HIV,
# OR (no SIM code) but TBweb closed as Obito NTB.
# For LTFU patients: Obito case_outcome doesn't apply (they're Abandono);
# rely entirely on SIM.
def attribute(row):
    if row["event_d"] != 1:
        return "no_death"
    # Tier 1: SIM ICD-10 if available
    if row["cod_known"]:
        if row["cod_tb_strict"] or row["cod_hiv_tb"]:
            return "tb_strict"
        elif row["cod_resp"]:
            return "respiratory"
        elif row["cod_hiv"]:
            return "hiv_other"
        else:
            return "non_tb"
    # Tier 2: TBweb attribution (only meaningful for Non-LTFU)
    if row["case_outcome"] == "Obito TB":
        return "tb_via_tbweb"
    if row["case_outcome"] == "Obito NTB":
        return "ntb_via_tbweb"
    return "unknown"

cohort["cod_class"] = cohort.apply(attribute, axis=1)

print("=" * 70)
print("Cause-of-death classification by arm (counts among event_d=1)")
print("=" * 70)
deaths = cohort[cohort["event_d"] == 1].copy()
print(f"Total deaths: {len(deaths):,}")
print()
print(pd.crosstab(deaths["cod_class"], deaths["itt_group"], margins=True))
print()
print("Same as %, by arm:")
ct = pd.crosstab(deaths["cod_class"], deaths["itt_group"], normalize="columns") * 100
print(ct.round(1))
print()

# ----------------------------------------------------------------------------
# Build TB / non-TB binary indicators for the analysis
# ----------------------------------------------------------------------------
# Strict TB: SIM-ICD-10 TB block OR TBweb Obito TB
# Broad TB-related: strict + HIV-related (B20.x) + respiratory (J)
# Non-TB: SIM-ICD-10 non-TB-non-resp-non-HIV OR TBweb Obito NTB
# Unknown: no SIM code AND (no Obito attribution OR Abandono case_outcome)
deaths["tb_strict_event"] = deaths["cod_class"].isin(["tb_strict", "tb_via_tbweb"])
deaths["tb_broad_event"]  = deaths["cod_class"].isin(
    ["tb_strict", "tb_via_tbweb", "respiratory", "hiv_other"]
)
deaths["non_tb_event"]    = deaths["cod_class"].isin(["non_tb", "ntb_via_tbweb"])
deaths["unknown_event"]   = deaths["cod_class"] == "unknown"

print("Binary TB-attribution indicators:")
for arm in ["Loss to follow-up", "Non-LTFU"]:
    d = deaths[deaths["itt_group"] == arm]
    n = len(d)
    print(f"  {arm} (N deaths = {n:,}):")
    print(f"    TB strict:    {d['tb_strict_event'].sum():>5,} ({100*d['tb_strict_event'].sum()/n:.1f}%)")
    print(f"    TB broad:     {d['tb_broad_event'].sum():>5,} ({100*d['tb_broad_event'].sum()/n:.1f}%)")
    print(f"    Non-TB:       {d['non_tb_event'].sum():>5,} ({100*d['non_tb_event'].sum()/n:.1f}%)")
    print(f"    Unknown:      {d['unknown_event'].sum():>5,} ({100*d['unknown_event'].sum()/n:.1f}%)")
    # check: should sum to N
    s = (d['tb_strict_event'].sum() + d['non_tb_event'].sum() +
         (d['cod_class'].isin(['respiratory', 'hiv_other'])).sum() +
         d['unknown_event'].sum())
    print(f"    (sum check:   {s:,} should equal {n:,})")
    print()

# Among the 'unknown' bucket: what's their itt_group, case_outcome, and timing?
unk = deaths[deaths["cod_class"] == "unknown"]
if len(unk):
    print(f"Unknown-bucket deaths (N={len(unk):,}):")
    print("  by itt_group:")
    print(unk["itt_group"].value_counts())
    print("  by case_outcome:")
    print(unk["case_outcome"].value_counts())
