"""
Check empirically how TBweb sets the LTFU end_date.

Two competing definitions:
  (A) Date of last visit. tx_days can be small (or zero — Abandono primario).
  (B) Date 30 days after last visit. Earliest possible tx_days = 30.

Plus we have three LTFU outcome labels in TBweb:
  - Abandono            : interrupted treatment after starting
  - Abandono Primario   : never returned after first prescription
  - Faltoso             : non-attender (older terminology)

If (B) is correct, none of these should have tx_days < 30.

Also break down by: % at exactly 30 days, % < 30 days, % at exactly 60 days.
A spike at 30 d would be a fingerprint of definition (B).
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

cohort = pd.read_csv(COHORT, low_memory=False,
                      parse_dates=["best_start", "end_date"])

# Recover the original case_outcome (Abandono / Abandono Primario / Faltoso)
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=["sinan_clean", "case_type", "case_outcome", "end_date"])
raw["end_date"] = pd.to_datetime(raw["end_date"], errors="coerce")
novo = raw[raw["case_type"].str.strip().str.title() == "Novo"]
novo = novo[novo["case_outcome"].notna()
            & (novo["case_outcome"].str.strip() != "")
            & (novo["case_outcome"] != "Mud Diag")]
TRANSFER = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
novo = novo[~novo["case_outcome"].isin(TRANSFER)]
first = novo.sort_values("end_date").groupby("sinan_clean", as_index=False).first()
cohort = cohort.merge(first[["sinan_clean", "case_outcome"]], on="sinan_clean", how="left")

ltfu = cohort[cohort["itt_group"] == "Loss to follow-up"].copy()
ltfu["tx_days"] = (ltfu["end_date"] - ltfu["best_start"]).dt.days

print(f"LTFU N = {len(ltfu):,}")
print()
print("Distribution of case_outcome among LTFU:")
print(ltfu["case_outcome"].value_counts(dropna=False))
print()

print("=" * 70)
print("Distribution of tx_days (end_date - best_start) by case_outcome")
print("=" * 70)
print()
for label in ["Abandono", "Abandono Primario", "Faltoso"]:
    sub = ltfu[ltfu["case_outcome"] == label]
    if len(sub) == 0:
        continue
    n = len(sub)
    print(f"--- {label} (N = {n:,}) ---")
    quantiles = sub["tx_days"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).to_dict()
    print(f"  median tx_days: {sub['tx_days'].median():.0f}")
    print(f"  p5 / p25 / p50 / p75 / p95: "
          f"{quantiles[0.05]:.0f} / {quantiles[0.25]:.0f} / "
          f"{quantiles[0.5]:.0f} / {quantiles[0.75]:.0f} / {quantiles[0.95]:.0f}")
    # Count near key thresholds
    for thr in [-1, 0, 7, 14, 21, 28, 29, 30, 31, 32, 45, 60, 90]:
        n_at = (sub["tx_days"] == thr).sum()
        n_le = (sub["tx_days"] <= thr).sum()
        if thr in [-1, 0, 7, 14, 28, 29, 30, 31, 32, 60]:
            print(f"    tx_days {thr:>4d}: exact={n_at:>5,}  cum (≤{thr}) = {n_le:>6,} ({100*n_le/n:5.2f}%)")
    print(f"    ANY tx_days < 30:  {(sub['tx_days'] < 30).sum():>6,} ({100*(sub['tx_days']<30).sum()/n:5.2f}%)")
    print(f"    tx_days exactly 30: {(sub['tx_days']==30).sum():>6,}  (would be huge spike under defn B)")
    print(f"    tx_days exactly 0:  {(sub['tx_days']==0).sum():>6,}  (only possible under defn A — Primario)")
    print()

# Histogram of tx_days for all LTFU, top 30 most common values
print("=" * 70)
print("Top 20 most common tx_days values among ALL LTFU:")
print("=" * 70)
print(ltfu["tx_days"].value_counts().head(20).sort_index())
print()

# Granular histogram, days 0-90
print("Frequency by day, 0-90:")
hist = ltfu["tx_days"].value_counts().sort_index()
hist = hist[(hist.index >= -5) & (hist.index <= 90)]
for d, n in hist.items():
    bar = "#" * min(80, int(n / max(hist.values) * 80))
    print(f"  day {int(d):>3d}: {n:>5,}  {bar}")

# For "Abandono" specifically (the most common label, the post-start LTFU)
print()
print("=" * 70)
print("'Abandono' subgroup specifically — fingerprint test")
print("=" * 70)
ab = ltfu[ltfu["case_outcome"] == "Abandono"].copy()
n = len(ab)
print(f"N = {n:,}")
print(f"  tx_days < 30:  {(ab['tx_days']<30).sum()} ({100*(ab['tx_days']<30).sum()/n:.2f}%)")
print(f"  tx_days = 30:  {(ab['tx_days']==30).sum()} ({100*(ab['tx_days']==30).sum()/n:.2f}%)")
print(f"  tx_days = 31:  {(ab['tx_days']==31).sum()}")
print(f"  tx_days = 60:  {(ab['tx_days']==60).sum()}")
print(f"  tx_days = 61:  {(ab['tx_days']==61).sum()}")

# Check: if defn B, the *minimum* tx_days for a real Abandono should be ~30
print(f"\n  min tx_days: {ab['tx_days'].min()}")
print(f"  if defn B holds, min should be ~30")
print(f"  10 smallest tx_days for Abandono:")
print(ab["tx_days"].nsmallest(10).tolist())
