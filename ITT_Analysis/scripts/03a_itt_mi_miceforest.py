"""
03a_itt_mi_miceforest.py
------------------------
Fast multiple imputation for the LTFU subgroup, replacing the R `mice::mice`
step in 03_itt_multiple_imputation_models.R.

Uses miceforest (LightGBM-backed MICE). Produces m=5 imputed datasets as CSVs
under ITT_Analysis/data/mi/imp_01.csv ... imp_05.csv, which the downstream R
script reads and pools with Rubin's rules via mice::as.mira() / mice::pool().

Variables imputed (must match 03_itt_multiple_imputation_models.R):
    categorical — age_group, sex, race_clean, edu_clean, hiv_aids, diabetes,
                  alcohol, drug_use, incarcerated, homelessness, hosp_admission,
                  clinical_clean, dot_status, tx_month_grp
    event/time  — fg_status, fg_time, event_d, time_d

Reference levels and NA-string handling mirror the R script.

2026-08-31 — added `tx_year`, the CALENDAR YEAR of treatment start, as a
categorical predictor in the imputation models. Calendar year is the strongest
single predictor of missingness in the cohort: race_clean is unrecorded for
10.6% of 2014 treatment starts, rising monotonically to 31.1% of 2023 starts
(edu_clean 23.7% -> 32.6%), a steeper gradient than for any patient
characteristic. The imputation model previously omitted it, so the
missing-at-random assumption was conditional on a set that left out the
dominant driver of the mechanism; conditioning on year brings the imputation
model closer to the plausible MAR set. Note this is about the IMPUTATION model
only — calendar year is a negligible CONFOUNDER of the exposure-mortality
association (adjusting for it moves the estimate by ~0.3%), so it is
deliberately NOT added to the analysis models.

Categorical (one level per year, 2013-2023) rather than a linear or spline
term, because the gradient is a recording-practice step change concentrated in
the last two years, not a smooth trend. `best_start` is fully observed
(0 unparseable of 171,048), so tx_year is a pure predictor and is never itself
imputed. It is carried into the output CSVs as a new column; nothing is
dropped, and downstream scripts select by name.

2026-08-31 (same day, second pass) — added two AUXILIARY predictors alongside
tx_year, at the owner's direction:

  * `mental_health` — the second-strongest predictor of EDUCATION missingness:
    edu_clean is unrecorded for 40.9% of the 2,622 patients with a recorded
    mental health condition against 26.9% of the 168,426 without. Fully
    observed (0 missing), so a pure predictor.
  * `dx_setting_aux` — diagnosis setting, which the manuscript already asserted
    was in the imputation model when it was not. Edu missingness runs from
    20.2% (active finding in the community) to 31.8% (emergency/inpatient).

Both are AUXILIARY: they are NOT among the 13 adjustment covariates and must
NOT be added to the analysis models. A variable can sharpen an imputation model
without belonging in the substantive adjustment set. edu_clean carries the
largest missingness in the study (27.1%), so its predictors matter most.

WHY `dx_setting_aux` AND NOT `diagnosis_setting` DIRECTLY. Contrary to the
brief, diagnosis_setting is NOT fully recorded: 4,178 patients (2.44%) have no
setting. Putting the raw column into the imputation frame would have made
miceforest impute it, silently turning a 4-variable imputation into a
5-variable one and handing downstream code a diagnosis_setting with no
missingness — several scripts (02_make_itt_table1.py, 60_make_manuscript.py,
55c_landmark_smd_psweighted.py) do `.fillna("Missing")` on it and would have
lost that category. So the auxiliary copy carries the missing rows as an
explicit "Unrecorded" level and the raw `diagnosis_setting` column passes
through untouched. This keeps the variable a pure predictor, as intended, and
the unrecorded group is itself the most informative level (37.9% edu
missingness, higher than any recorded setting) rather than something to fill in.

Net effect on the frame: three predictors added (tx_year, mental_health,
dx_setting_aux); the set of variables actually IMPUTED is unchanged at four
(race_clean 13.25%, edu_clean 27.06%, dot_status 6.97%, clinical_clean 1 row).
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import miceforest as mf


# ---------------------------------------------------------------------------
# Project root resolution — mirrors 01_itt_cohort_selection.py
# ---------------------------------------------------------------------------
def _find_project_root() -> Path:
    env = os.environ.get("TB_ABANDONMENT_ROOT")
    if env:
        p = Path(env).expanduser()
        if p.exists():
            return p
    candidates = [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/TB SP 2026/LTFU Paper",
    ]
    for c in candidates:
        if c.exists():
            return c
    return Path(__file__).resolve().parents[2]


BASE_DIR = _find_project_root()
COHORT_CSV = BASE_DIR / "ITT_Analysis" / "data" / "itt_cohort.csv"
OUT_DIR = BASE_DIR / "ITT_Analysis" / "data" / "mi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"[paths] BASE_DIR     = {BASE_DIR}")
print(f"[paths] cohort input = {COHORT_CSV}")
print(f"[paths] MI output    = {OUT_DIR}")

if not COHORT_CSV.exists():
    raise FileNotFoundError(f"itt_cohort.csv not found at {COHORT_CSV}")


# ---------------------------------------------------------------------------
# Load cohort. MI is performed on the FULL cohort (both LTFU and Non-LTFU
# arms). Downstream models in 03_itt_multiple_imputation_models.R filter
# to LTFU after reading the imputed datasets — imputing on the full cohort
# borrows information across arms and yields better conditional
# distributions for the imputation models.
# ---------------------------------------------------------------------------
df = pd.read_csv(COHORT_CSV, low_memory=False)
print(f"Full cohort N = {len(df):,}  (itt_group: {df['itt_group'].value_counts().to_dict()})")

# Match the R script's NA-string harmonization on the same columns
NA_STRINGS = {"Missing", "Ignorado", "Unknown", "", "nan", "NaN"}
NA_COLS = ["race_clean", "edu_clean", "dot_status", "alcohol", "drug_use",
           "diabetes", "hosp_admission", "hiv_aids"]
for c in NA_COLS:
    df[c] = df[c].apply(lambda x: np.nan if pd.isna(x) or str(x) in NA_STRINGS else x)

# Times: floor at 0.001 to avoid Surv() barfing on zeros (matches R script)
for c in ("time_rn", "time_d"):
    df[c] = pd.to_numeric(df[c], errors="coerce").clip(lower=0.001)

# Derive Fine-Gray composite status and time
#   1 = retreatment, 2 = death (pre-retreatment), 0 = censored
event_rn = pd.to_numeric(df["event_rn"], errors="coerce").fillna(0).astype(int)
event_d = pd.to_numeric(df["event_d"], errors="coerce").fillna(0).astype(int)
df["fg_status"] = np.select(
    [event_rn == 1, (event_d == 1) & (event_rn == 0)],
    [1, 2],
    default=0,
)
df["fg_time"] = np.where(event_rn == 1, df["time_rn"], df["time_d"])

# Age group bucket (same cuts as R: <25, 25-44, 45-64, 65+)
df["age_group"] = pd.cut(
    pd.to_numeric(df["age_tb"], errors="coerce"),
    bins=[14, 24, 44, 64, 150],
    labels=["15-24", "25-44", "45-64", "≥65"],
)

# Calendar year of treatment start (2026-08-31). One level per year, NOT a
# linear term — the missingness gradient is a recording-practice step change,
# not a smooth trend. `best_start` is complete, so this must not introduce any
# missingness of its own; assert that rather than trusting it.
_start = pd.to_datetime(df["best_start"], errors="coerce")
if _start.isna().any():
    raise ValueError(
        f"best_start unparseable for {int(_start.isna().sum())} rows; tx_year "
        f"would enter the imputation frame with missingness of its own"
    )
df["tx_year"] = _start.dt.year.astype(int).astype(str)

# Auxiliary predictor: diagnosis setting with its 4,178 unrecorded patients
# (2.44%) held as an explicit level. A separate column so the raw
# `diagnosis_setting` still passes through with its NAs intact — downstream code
# labels those "Missing" and must keep being able to. See the docstring.
df["dx_setting_aux"] = (
    df["diagnosis_setting"]
    .apply(lambda x: np.nan if pd.isna(x) or str(x) in NA_STRINGS else x)
    .fillna("Unrecorded")
    .astype(str)
)

# `mental_health` needs no derivation, but it must be fully observed to stay a
# pure predictor rather than becoming a fifth imputed variable.
_mh = df["mental_health"].apply(
    lambda x: np.nan if pd.isna(x) or str(x) in NA_STRINGS else x)
if _mh.isna().any():
    raise ValueError(
        f"mental_health missing for {int(_mh.isna().sum())} rows; it was added "
        f"as an AUXILIARY predictor on the stated basis that it is fully "
        f"recorded, so handle its missingness explicitly rather than letting "
        f"miceforest impute it"
    )

# ---------------------------------------------------------------------------
# Build the MI frame with correct dtypes. miceforest infers by dtype:
#   category -> categorical imputation; float -> numeric imputation.
# ---------------------------------------------------------------------------
CATEGORICAL = [
    "age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
    "alcohol", "drug_use", "incarcerated", "homelessness", "hosp_admission",
    "clinical_clean", "dot_status", "tx_month_grp",
    # calendar year of treatment start (2026-08-31) — fully observed, included
    # as a predictor because it is the dominant driver of the missingness
    # mechanism for race_clean and edu_clean. `tx_month_grp` above is the month
    # of THERAPY (duration on treatment), which carries no calendar-time
    # information, so year was genuinely absent from the model before this.
    "tx_year",
    # auxiliary predictors (2026-08-31) — fully observed, in the model only to
    # sharpen the conditional distribution of edu_clean, which carries the
    # largest missingness in the study. NOT adjustment covariates: do not add
    # them to the analysis models in the downstream R scripts.
    "mental_health", "dx_setting_aux",
    # exposure — not missing, but included as a predictor in the imputation
    # models so that LTFU vs Non-LTFU heterogeneity informs conditional
    # distributions for the imputed variables.
    "itt_group",
]
# fg_status and event_d are categorical at the R side; keep them so here too
CATEGORICAL += ["fg_status", "event_d"]

NUMERIC = ["fg_time", "time_d"]

MI_COLS = CATEGORICAL + NUMERIC

df_mi = df[MI_COLS].copy().reset_index(drop=True)
for c in CATEGORICAL:
    df_mi[c] = df_mi[c].astype("category")
for c in NUMERIC:
    df_mi[c] = pd.to_numeric(df_mi[c], errors="coerce")

print("\n[pre-MI] missing counts per variable:")
print(df_mi.isna().sum().to_string())

# ---------------------------------------------------------------------------
# Run miceforest
# ---------------------------------------------------------------------------
M = 5
MAXIT = 5
SEED = 42

print(f"\nStarting miceforest imputation (m={M}, iterations={MAXIT}) ...")
# mean_match_candidates=0 tells miceforest to skip KDTree-based PMM and use
# raw LightGBM predictions directly. We were tripping over a KDTree "data
# must be finite" error on some variables; raw-prediction imputation is
# cleaner and plenty accurate at this sample size.
kernel = mf.ImputationKernel(
    df_mi,
    num_datasets=M,
    random_state=SEED,
    mean_match_candidates=0,
)
kernel.mice(iterations=MAXIT, verbose=True)
print("Imputation finished.")

# ---------------------------------------------------------------------------
# Attach the non-imputed columns the downstream R script needs for modeling
# (everything except the 18 MI vars). Write one CSV per imputation.
# ---------------------------------------------------------------------------
PASS_THROUGH = [c for c in df.columns if c not in MI_COLS]

for i in range(M):
    completed = kernel.complete_data(dataset=i)
    out = pd.concat(
        [completed.reset_index(drop=True), df[PASS_THROUGH].reset_index(drop=True)],
        axis=1,
    )
    out_path = OUT_DIR / f"imp_{i + 1:02d}.csv"
    out.to_csv(out_path, index=False)
    print(f"  wrote {out_path}  (rows={len(out):,})")

print("\nDone. Downstream script: 03_itt_multiple_imputation_models.R")
