#!/usr/bin/env python3
"""Could the missing-outcome exclusion hide undocumented disengagement?

Appendix 1.1 says that of the individuals excluded for a missing or incomplete
programmatic outcome, only a handful had a usable treatment-closure date and
start date, and fewer still met every other eligibility criterion -- so
re-classifying them as lost to follow-up could not move any estimate. Those
counts were typed once against the old exclusion chain (12,430 excluded) and had
no script behind them. This derives them from the same data and the same filter
definitions as 61_make_figS1_consort.py.

For each person in the missing-outcome group, their earliest first-episode
("Novo") record is taken, and counted if:
  step 1  end date inside 2013-01-01 .. 2023-12-31 and a recorded treatment
          start not after the end date (the 61 date checks), and
  step 2  also aged >= 15, not dead on or before treatment start, and not
          identified only post-mortem (every remaining criterion; transfer
          cannot apply to a record with no outcome).

Usage:  python3 ITT_Analysis/scripts/61b_missing_outcome_reclassification.py
Output: ITT_Analysis/results/missing_outcome_reclassification.csv
"""
import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
OUT = ROOT / "ITT_Analysis" / "results" / "missing_outcome_reclassification.csv"

df = pd.read_csv(ROOT / "Data" / "Final_table_cleaned.csv", low_memory=False)
for c in ("end_date", "dod", "tx_start"):
    df[c] = pd.to_datetime(df[c], errors="coerce")

# comprehensive death date, as in 61
a = df.dropna(subset=["dod"]).groupby("sinan_clean")["dod"].max()
b = (df[df.case_outcome.isin({"Obito TB", "Obito NTB"})].dropna(subset=["end_date"])
       .groupby("sinan_clean")["end_date"].max())
death = pd.concat([a, b], axis=1).max(axis=1).rename("death_date")

novo = df[df.case_type.str.strip().str.title() == "Novo"]
blank = novo.case_outcome.isna() | (novo.case_outcome.astype(str).str.strip() == "")
mud = novo.case_outcome.eq("Mud Diag")
valid_ids = set(novo[~blank & ~mud].sinan_clean)
lost = novo[~novo.sinan_clean.isin(valid_ids)]
missing_ids = set(lost[blank.loc[lost.index]].sinan_clean)   # as n_missing in 61

g = (novo[novo.sinan_clean.isin(missing_ids)].sort_values("end_date")
     .groupby("sinan_clean", as_index=False).first().join(death, on="sinan_clean"))

dates_ok = (g.end_date.notna()
            & (g.end_date >= pd.Timestamp("2013-01-01"))
            & (g.end_date <= pd.Timestamp("2023-12-31"))
            & g.tx_start.notna() & ~(g.end_date < g.tx_start))
rest_ok = ((g.age_tb >= 15)
           & ~(g.death_date.notna() & (g.death_date <= g.tx_start))
           & (g.disease_discovery != "Descob. Apos Obito"))

rows = [
    ("excluded_missing_outcome", len(missing_ids)),
    ("valid_window_and_start_dates", int(dates_ok.sum())),
    ("meet_all_other_criteria", int((dates_ok & rest_ok).sum())),
]
out = pd.DataFrame(rows, columns=["step", "n"])
if OUT.exists():
    raise SystemExit(f"{OUT.name} exists; move it aside to regenerate")
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
print(f"\n  wrote {OUT}")
