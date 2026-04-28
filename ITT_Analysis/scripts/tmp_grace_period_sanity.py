"""
tmp_grace_period_sanity.py
--------------------------
Calibration counts for the 30-day classification immortal-time issue.

(1) LTFU arm: deaths occurring within 30 days of end_date.
    These are SIM-ascertained deaths that "shouldn't" exist under a strict
    application of the 30-day-no-contact rule, because a patient who died
    within 30 days of last visit would have been classified as
    "Obito on treatment", not abandono. They exist because TBweb already
    closed the case as abandono before SIM caught the death.

(2) Non-LTFU arm: Obito-on-treatment deaths with short tx duration.
    Bounds the population that *might* have been would-be-LTFU but were
    reclassified to died-on-treatment under the 30-day rule.

Reads the cohort and the upstream cleaned table to recover case_outcome
for the non-LTFU subgroup (the cohort itself collapses outcomes to
itt_group).
"""
import os
from pathlib import Path
import pandas as pd
import numpy as np


def _find_root():
    env = os.environ.get("TB_ABANDONMENT_ROOT")
    if env and Path(env).exists():
        return Path(env)
    for c in [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]:
        if c.exists():
            return c
    raise FileNotFoundError("project root not found")


BASE = _find_root()
COHORT_CSV = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"
RAW_CSV = BASE / "Data" / "Final_table_cleaned.csv"

print(f"Cohort: {COHORT_CSV}")
print(f"Raw:    {RAW_CSV}")
print()

cohort = pd.read_csv(COHORT_CSV, low_memory=False, parse_dates=["best_start", "end_date", "death_date"])

# ---------------------------------------------------------------------------
# Recover case_outcome for the cohort by joining back to the raw table on the
# index episode (first Novo). We replicate the cohort selection: first Novo
# per sinan_clean among valid (non-MudDiag, non-empty) outcomes.
# ---------------------------------------------------------------------------
raw = pd.read_csv(RAW_CSV, low_memory=False, usecols=[
    "sinan_clean", "case_type", "case_outcome", "end_date"
])
raw["end_date"] = pd.to_datetime(raw["end_date"], errors="coerce")
novo = raw[raw["case_type"].str.strip().str.title() == "Novo"]
novo = novo[novo["case_outcome"].notna()
            & (novo["case_outcome"].str.strip() != "")
            & (novo["case_outcome"] != "Mud Diag")]
TRANSFER = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
novo = novo[~novo["case_outcome"].isin(TRANSFER)]
first = novo.sort_values("end_date").groupby("sinan_clean", as_index=False).first()
first = first[["sinan_clean", "case_outcome"]]

cohort = cohort.merge(first, on="sinan_clean", how="left")
print(f"Cohort N: {len(cohort):,}")
print(f"With case_outcome recovered: {cohort['case_outcome'].notna().sum():,}")
print()

# ---------------------------------------------------------------------------
# (1) LTFU arm: deaths within 30 days of end_date
# ---------------------------------------------------------------------------
ltfu = cohort[cohort["itt_group"] == "Loss to follow-up"].copy()
ltfu["days_to_death"] = (ltfu["death_date"] - ltfu["end_date"]).dt.days

n_ltfu = len(ltfu)
n_ltfu_died = (ltfu["event_d"] == 1).sum()

# Days from end_date to death — strictly post-end (event_d=1 means death after end_date)
deaths_le_30 = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] <= 30)).sum()
deaths_le_60 = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] <= 60)).sum()
deaths_le_90 = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] <= 90)).sum()
deaths_le_180 = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] <= 180)).sum()
deaths_le_365 = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] <= 365)).sum()

# Negative days_to_death are technically excluded from event_d=1 by cohort-build
# logic (death_date >= end_date), but check anyway.
neg = ((ltfu["event_d"] == 1) & (ltfu["days_to_death"] < 0)).sum()
print("=" * 70)
print("(1) LTFU arm: timing of post-end_date deaths")
print("=" * 70)
print(f"  LTFU N:                          {n_ltfu:,}")
print(f"  LTFU with any post-end death:    {n_ltfu_died:,}  ({100*n_ltfu_died/n_ltfu:.2f}% of LTFU)")
print(f"  ... within 30 days of end_date:  {deaths_le_30:,}  ({100*deaths_le_30/n_ltfu_died:.2f}% of LTFU deaths;  {100*deaths_le_30/n_ltfu:.3f}% of LTFU)")
print(f"  ... within 60 days:              {deaths_le_60:,}  ({100*deaths_le_60/n_ltfu_died:.2f}% of LTFU deaths)")
print(f"  ... within 90 days:              {deaths_le_90:,}  ({100*deaths_le_90/n_ltfu_died:.2f}% of LTFU deaths)")
print(f"  ... within 180 days:             {deaths_le_180:,}  ({100*deaths_le_180/n_ltfu_died:.2f}% of LTFU deaths)")
print(f"  ... within 365 days:             {deaths_le_365:,}  ({100*deaths_le_365/n_ltfu_died:.2f}% of LTFU deaths)")
if neg > 0:
    print(f"  (sanity: {neg} death dates fall *before* end_date — should be 0)")

# Distribution of days-to-death for LTFU deaths
ltfu_died = ltfu[(ltfu["event_d"] == 1) & ltfu["days_to_death"].notna()]
print()
print("  LTFU days-to-death distribution (for those who died):")
print(f"    median: {ltfu_died['days_to_death'].median():.0f} d")
print(f"    p10:    {ltfu_died['days_to_death'].quantile(0.10):.0f} d")
print(f"    p25:    {ltfu_died['days_to_death'].quantile(0.25):.0f} d")
print(f"    p75:    {ltfu_died['days_to_death'].quantile(0.75):.0f} d")
print()

# ---------------------------------------------------------------------------
# (2) Non-LTFU "Obito" deaths with short tx duration
#     case_outcome in {Obito TB, Obito NTB}
#     tx duration = end_date - best_start
# ---------------------------------------------------------------------------
nonltfu = cohort[cohort["itt_group"] == "Non-LTFU"].copy()
nonltfu["tx_days"] = (nonltfu["end_date"] - nonltfu["best_start"]).dt.days

obito_set = {"Obito TB", "Obito NTB"}
obito = nonltfu[nonltfu["case_outcome"].isin(obito_set)].copy()

n_nonltfu = len(nonltfu)
n_obito = len(obito)

print("=" * 70)
print("(2) Non-LTFU arm: Obito-on-treatment timing")
print("=" * 70)
print(f"  Non-LTFU N:                      {n_nonltfu:,}")
print(f"  ... case_outcome in Obito set:   {n_obito:,}  ({100*n_obito/n_nonltfu:.2f}% of non-LTFU)")
print()

# Distribution of tx_days at the time the Obito case was closed
print("  Obito-on-treatment tx_days = (end_date - best_start):")
print(f"    median: {obito['tx_days'].median():.0f} d")
print(f"    p10:    {obito['tx_days'].quantile(0.10):.0f} d")
print(f"    p25:    {obito['tx_days'].quantile(0.25):.0f} d")
print(f"    p75:    {obito['tx_days'].quantile(0.75):.0f} d")
print()

bands = [(0, 7), (0, 14), (0, 30), (0, 60), (0, 90), (0, 180)]
for lo, hi in bands:
    n = ((obito["tx_days"] >= lo) & (obito["tx_days"] <= hi)).sum()
    print(f"    tx_days in [{lo:>3d}, {hi:>3d}]: {n:>6,}  ({100*n/n_obito:5.2f}% of Obito-on-tx)")

print()
print("  Among these, how short was the time from death to end_date?")
print("  (death_date and end_date should be ≈equal for died-on-tx;")
print("   that's how the case_outcome=Obito flips. Confirm:)")
obito["dod_minus_end"] = (obito["death_date"] - obito["end_date"]).dt.days
print(f"    median |death_date - end_date|: {obito['dod_minus_end'].abs().median():.1f} d")
print(f"    pct with death_date == end_date: "
      f"{100*(obito['dod_minus_end'] == 0).sum()/obito['dod_minus_end'].notna().sum():.1f}%")
print()

# How many Obito-on-tx deaths might plausibly have been would-be-LTFU?
# We don't know last-visit date, but a hard upper bound is:
# tx_days < some threshold (e.g., 60d) AND in age/risk profile resembling LTFU
# We just report the count under several thresholds.
print("  Bounds on 'would-be-LTFU' under various tx_days cutoffs:")
for thr in [30, 60, 90, 120, 180]:
    n = (obito["tx_days"] <= thr).sum()
    print(f"    Obito with tx_days ≤ {thr:>3d}: {n:>6,}  "
          f"({100*n/n_nonltfu:5.3f}% of non-LTFU; "
          f"{100*n/(n_nonltfu + n):5.3f}% of an expanded LTFU group)")
print()

# ---------------------------------------------------------------------------
# Quick proportion-LTFU shift if we reassigned all Obito-on-tx with tx_days≤T
# to LTFU
# ---------------------------------------------------------------------------
print("=" * 70)
print("If we reassigned Obito-on-tx deaths to LTFU under a tx_days ≤ T rule:")
print("=" * 70)
print(f"  Current LTFU N:                  {n_ltfu:,}  ({100*n_ltfu/len(cohort):.2f}%)")
for thr in [30, 60, 90]:
    n_add = (obito["tx_days"] <= thr).sum()
    new_ltfu = n_ltfu + n_add
    print(f"  + Obito with tx_days ≤ {thr:>3d}d:    "
          f"+{n_add:,}  →  {new_ltfu:,}  "
          f"({100*new_ltfu/len(cohort):.2f}%)")
