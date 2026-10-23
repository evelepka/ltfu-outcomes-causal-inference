#!/usr/bin/env python3
"""What share of ELIGIBLE individuals does the analysed cohort represent?

Reviewer 3 (R2) noted that ~72% of individuals in TBweb are analysed. Most of
the gap is eligibility (not a first new-case episode, reclassified as not TB,
aged < 15, death before treatment, treatment outside the study period,
diagnosed only post-mortem) rather than exclusion of eligible people for data
problems (missing outcome, transfer during treatment, unusable treatment
dates). This classifies every individual once, on their own record, with the
filter definitions of 61_make_figS1_consort.py, and reports the cohort as a
share of the eligible.

Record used per person: the earliest first-episode ("Novo") record with a valid
outcome, as in 61; for people with no valid outcome, their earliest Novo record.

Eligibility is tested before data quality, so a person failing both is counted
as ineligible. People whose record has no end date cannot be placed in the
study period; they are reported separately (no_end_date) and the share is
given with them counted both ways.

Usage:  python3 ITT_Analysis/scripts/61c_eligible_denominator.py
Output: ITT_Analysis/results/cohort_eligible_denominator.csv
"""
import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
OUT = ROOT / "ITT_Analysis" / "results" / "cohort_eligible_denominator.csv"
if OUT.exists():
    raise SystemExit(f"{OUT.name} exists; move it aside to regenerate")

df = pd.read_csv(ROOT / "Data" / "Final_table_cleaned.csv", low_memory=False)
for c in ("end_date", "dod", "tx_start"):
    df[c] = pd.to_datetime(df[c], errors="coerce")
a = df.dropna(subset=["dod"]).groupby("sinan_clean")["dod"].max()
b = (df[df.case_outcome.isin({"Obito TB", "Obito NTB"})].dropna(subset=["end_date"])
       .groupby("sinan_clean")["end_date"].max())
death = pd.concat([a, b], axis=1).max(axis=1).rename("death_date")

n_total = df.sinan_clean.nunique()
novo = df[df.case_type.str.strip().str.title() == "Novo"]
blank = novo.case_outcome.isna() | (novo.case_outcome.astype(str).str.strip() == "")
mud = novo.case_outcome.eq("Mud Diag")
valid = novo[~blank & ~mud]
lost = novo[~novo.sinan_clean.isin(set(valid.sinan_clean))]
s_blank = set(lost[blank.loc[lost.index]].sinan_clean)
s_mud = set(lost[mud.loc[lost.index]].sinan_clean) - s_blank

rec = pd.concat([
    valid.sort_values("end_date").groupby("sinan_clean", as_index=False).first(),
    lost.sort_values("end_date").groupby("sinan_clean", as_index=False).first(),
]).join(death, on="sinan_clean").set_index("sinan_clean")
rec["missing_outcome"] = rec.index.isin(s_blank)
rec["mud_diag"] = rec.index.isin(s_mud)

TR = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
no_end = rec.end_date.isna()
in_window = (rec.end_date >= pd.Timestamp("2013-01-01")) & (rec.end_date <= pd.Timestamp("2023-12-31"))
# ~(age >= 15), not (age < 15): 61 keeps only age >= 15, so a missing age is excluded
inelig = (rec.mud_diag | ~(rec.age_tb >= 15)
          | (rec.death_date.notna() & (rec.death_date <= rec.tx_start))
          | (~no_end & ~in_window)
          | (rec.disease_discovery == "Descob. Apos Obito"))
quality = (rec.missing_outcome | rec.case_outcome.isin(TR)
           | rec.tx_start.isna() | (rec.end_date < rec.tx_start))

cohort = ~inelig & ~no_end & ~quality
n_cohort = int(cohort.sum())
expected = len(pd.read_csv(ROOT / "ITT_Analysis" / "data" / "itt_cohort.csv",
                           usecols=["sinan_clean"], low_memory=False))
assert n_cohort == expected, f"classification gives {n_cohort:,}, cohort file has {expected:,}"

elig_with_noend = int((~inelig).sum())             # no end date counted eligible
elig_without = int((~inelig & ~no_end).sum())      # no end date counted ineligible
rows = [
    ("tbweb_individuals", n_total),
    ("no_first_episode_record", n_total - novo.sinan_clean.nunique()),
    ("first_episode_individuals", len(rec)),
    ("ineligible", int(inelig.sum())),
    ("no_end_date_not_otherwise_ineligible", int((~inelig & no_end).sum())),
    ("eligible_excl_no_end_date", elig_without),
    ("eligible_incl_no_end_date", elig_with_noend),
    ("cohort", n_cohort),
    ("pct_of_tbweb", round(100 * n_cohort / n_total, 2)),
    ("pct_of_eligible_excl_no_end_date", round(100 * n_cohort / elig_without, 2)),
    ("pct_of_eligible_incl_no_end_date", round(100 * n_cohort / elig_with_noend, 2)),
    ("eligible_excluded_missing_outcome", int((~inelig & ~no_end & rec.missing_outcome).sum())),
    ("eligible_excluded_transfer", int((~inelig & ~no_end & ~rec.missing_outcome & rec.case_outcome.isin(TR)).sum())),
    ("eligible_excluded_treatment_dates", int((~inelig & ~no_end & ~rec.missing_outcome & ~rec.case_outcome.isin(TR) & quality).sum())),
]
# the manuscript's wording: adults whose first episode ended within the study
# period, before the remaining eligibility filters
adult_win = (rec.age_tb >= 15) & in_window
pre = rec.death_date.notna() & (rec.death_date <= rec.tx_start)
for name, d in (("adult_first_episode_end_in_window", adult_win),
                ("adult_end_in_window_excl_change_of_diagnosis", adult_win & ~rec.mud_diag),
                ("adult_end_in_window_excl_cod_and_pretreatment_death", adult_win & ~rec.mud_diag & ~pre)):
    rows += [(name, int(d.sum())), ("pct_of_" + name, round(100 * n_cohort / d.sum(), 2))]
out = pd.DataFrame(rows, columns=["quantity", "value"])
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
print(f"\n  wrote {OUT}")
