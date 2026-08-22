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

# ---------------------------------------------------------------------------
# Build the MI frame with correct dtypes. miceforest infers by dtype:
#   category -> categorical imputation; float -> numeric imputation.
# ---------------------------------------------------------------------------
CATEGORICAL = [
    "age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
    "alcohol", "drug_use", "incarcerated", "homelessness", "hosp_admission",
    "clinical_clean", "dot_status", "tx_month_grp",
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
