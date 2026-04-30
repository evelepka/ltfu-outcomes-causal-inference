"""
Sanity checks prompted by Sueli's email (2026-04-30) on TBweb closure rules.

Two questions:
  Q1. Among LTFU patients labelled `Abandono Primario`, how many lack a
      recorded tx_start? Per Sueli, these are "diagnosis confirmed, never
      started treatment" — closed >=90d after notification. Including
      them in an LTFU-vs-on-treatment contrast is conceptually a
      category error since they never received the intervention.

  Q2. Is the 2016 introduction of `Abandono Primario` as a closure code
      visible as a step change in per-label LTFU composition? If pre-2016
      `Abandono` absorbs what would later be `Primario`, we expect
      Primario fraction to rise sharply at 2016.

Read-only: queries cohort + raw TBweb pull. Does not modify cohort.
"""
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
    raise FileNotFoundError("Could not locate Abandonment Paper GDrive root")


BASE = _root()
COHORT = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"
RAW = BASE / "Data" / "Final_table_cleaned.csv"

print(f"Reading cohort: {COHORT}")
cohort = pd.read_csv(COHORT, low_memory=False,
                     parse_dates=["best_start", "end_date"])

print(f"Reading raw: {RAW}")
raw = pd.read_csv(
    RAW,
    low_memory=False,
    usecols=["sinan_clean", "case_type", "case_outcome",
             "end_date", "tx_start", "notification_date", "diagnostic_date"],
)
for c in ["end_date", "tx_start", "notification_date", "diagnostic_date"]:
    raw[c] = pd.to_datetime(raw[c], errors="coerce")

# Recover index-episode case_outcome + index tx_start (same logic as 01_)
TRANSFER = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
novo = raw[raw["case_type"].astype(str).str.strip().str.title() == "Novo"]
novo = novo[novo["case_outcome"].notna()
           & (novo["case_outcome"].astype(str).str.strip() != "")
           & (novo["case_outcome"] != "Mud Diag")
           & ~novo["case_outcome"].isin(TRANSFER)]
first = (novo.sort_values("end_date")
              .groupby("sinan_clean", as_index=False)
              .first()[["sinan_clean", "case_outcome", "tx_start",
                        "notification_date", "diagnostic_date"]])

cohort = cohort.merge(first, on="sinan_clean", how="left",
                     suffixes=("", "_raw"))

ltfu = cohort[cohort["itt_group"] == "Loss to follow-up"].copy()
print(f"\nLTFU N = {len(ltfu):,}")
print("LTFU sub-label counts:")
print(ltfu["case_outcome"].value_counts(dropna=False))

# ============================================================
# Q1: Primario without tx_start
# ============================================================
print("\n" + "=" * 70)
print("Q1: Abandono Primario WITHOUT recorded tx_start (never-treated)")
print("=" * 70)
prim = ltfu[ltfu["case_outcome"] == "Abandono Primario"].copy()
n_prim = len(prim)
n_prim_no_tx = prim["tx_start"].isna().sum()
n_prim_with_tx = prim["tx_start"].notna().sum()
print(f"Primario total in LTFU cohort:        {n_prim:>6,}")
print(f"  with tx_start recorded:             {n_prim_with_tx:>6,} "
      f"({100*n_prim_with_tx/n_prim:.1f}%)")
print(f"  WITHOUT tx_start (never-treated):   {n_prim_no_tx:>6,} "
      f"({100*n_prim_no_tx/n_prim:.1f}%)")
print(f"As share of full LTFU cohort:         "
      f"{100*n_prim_no_tx/len(ltfu):.2f}%")
print(f"As share of full ITT cohort:          "
      f"{100*n_prim_no_tx/len(cohort):.2f}%")

# Compare end_date - notification spread for the two groups
prim["dx_to_close"] = (prim["end_date"] - prim["notification_date"]).dt.days
print("\nDays from notification to end_date (Primario):")
for grp, sub in prim.groupby(prim["tx_start"].isna()):
    label = "no tx_start" if grp else "with tx_start"
    q = sub["dx_to_close"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).to_dict()
    print(f"  {label:>14s} (N={len(sub):>5,}): "
          f"p5/p25/p50/p75/p95 = "
          f"{q[0.05]:.0f}/{q[0.25]:.0f}/{q[0.5]:.0f}/"
          f"{q[0.75]:.0f}/{q[0.95]:.0f}")

# Sueli says the no-tx_start case requires >=90d post-notification — verify
no_tx = prim[prim["tx_start"].isna()]
if len(no_tx) > 0:
    n_under_90 = (no_tx["dx_to_close"] < 90).sum()
    print(f"\n  Sueli rule check: >=90d notification->close required.")
    print(f"  N with notif->close < 90d: {n_under_90:,} "
          f"({100*n_under_90/len(no_tx):.1f}% of no-tx Primario) "
          f"— expected near zero post-2016/post-rule-enforcement.")

# What does best_start look like for no-tx Primario? It gets proxied from
# diagnostic_date or notification_date.
no_tx_in_cohort = ltfu[(ltfu["case_outcome"] == "Abandono Primario")
                      & ltfu["tx_start"].isna()].copy()
no_tx_in_cohort["tx_days"] = (no_tx_in_cohort["end_date"]
                              - no_tx_in_cohort["best_start"]).dt.days
print(f"\n  tx_days (end_date - best_start proxy) for no-tx Primario:")
q = no_tx_in_cohort["tx_days"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).to_dict()
print(f"    p5/p25/p50/p75/p95 = "
      f"{q[0.05]:.0f}/{q[0.25]:.0f}/{q[0.5]:.0f}/"
      f"{q[0.75]:.0f}/{q[0.95]:.0f}")
print(f"    These patients enter follow-up at a synthetic origin "
      f"(diag/notif date) but never received TB therapy.")

# ============================================================
# Q2: 2016 step change in Primario share?
# ============================================================
print("\n" + "=" * 70)
print("Q2: Sub-label composition of LTFU by calendar year of end_date")
print("=" * 70)
ltfu["year"] = ltfu["end_date"].dt.year
crosstab = pd.crosstab(ltfu["year"], ltfu["case_outcome"])
crosstab["TOTAL"] = crosstab.sum(axis=1)
for col in ["Abandono", "Abandono Primario", "Faltoso"]:
    if col in crosstab.columns:
        crosstab[f"{col}_pct"] = (100 * crosstab[col] /
                                  crosstab["TOTAL"]).round(1)

cols_show = ["Abandono", "Abandono Primario", "Faltoso", "TOTAL",
             "Abandono_pct", "Abandono Primario_pct", "Faltoso_pct"]
cols_show = [c for c in cols_show if c in crosstab.columns]
print(crosstab[cols_show].to_string())

# Save outputs for the manuscript Appendix
OUT = BASE / "ITT_Analysis" / "results" / "Tables"
OUT.mkdir(parents=True, exist_ok=True)
crosstab[cols_show].to_csv(OUT / "S_ltfu_sublabel_by_year.csv")
print(f"\nWrote: {OUT / 'S_ltfu_sublabel_by_year.csv'}")

# Also dump a concise summary table for Q1
q1 = pd.DataFrame({
    "metric": [
        "LTFU total",
        "Abandono",
        "Faltoso",
        "Abandono Primario (total)",
        "Abandono Primario, with tx_start",
        "Abandono Primario, NO tx_start (never-treated)",
        "Never-treated as % of LTFU",
        "Never-treated as % of ITT cohort",
    ],
    "value": [
        f"{len(ltfu):,}",
        f"{(ltfu['case_outcome'] == 'Abandono').sum():,}",
        f"{(ltfu['case_outcome'] == 'Faltoso').sum():,}",
        f"{n_prim:,}",
        f"{n_prim_with_tx:,}",
        f"{n_prim_no_tx:,}",
        f"{100*n_prim_no_tx/len(ltfu):.2f}%",
        f"{100*n_prim_no_tx/len(cohort):.2f}%",
    ],
})
q1.to_csv(OUT / "S_primario_summary.csv", index=False)
print(f"Wrote: {OUT / 'S_primario_summary.csv'}")

# ============================================================
# Exclusion list for sensitivity analysis: never-treated Primario
# ============================================================
excl = ltfu[(ltfu["case_outcome"] == "Abandono Primario")
            & ltfu["tx_start"].isna()][["sinan_clean"]].copy()
excl_path = (BASE / "ITT_Analysis" / "data"
             / "primario_no_txstart_ids.csv")
excl.to_csv(excl_path, index=False)
print(f"Wrote {len(excl):,} exclusion IDs to: {excl_path}")
