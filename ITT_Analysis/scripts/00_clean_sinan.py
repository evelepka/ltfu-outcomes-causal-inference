"""
00_clean_sinan.py
=================
Clean the SINAN TB notification dataset.

Problem
-------
The individual identifier 'sinan' is a 7-digit number, but leading zeros
were dropped during data export.  As a result, the same person can appear
as e.g. "00144" and "144".  To confirm two records belong to the same
person, date-of-birth (dob) AND sex must match.

Strategy
--------
1. Zero-pad sinan to 7 digits  →  sinan_padded
2. Detect any padded sinan that corresponds to records with inconsistent
   dob / sex (i.e. genuine different people who happen to share the same
   numeric ID after padding).
3. For conflicting groups, assign a new unique ID so that records sharing
   the SAME (sinan_padded, dob, sex) cluster stay together.
4. Flag every row with audit columns for downstream transparency.

Output
------
Cleaned CSV saved to: tb_recurrence_analysis/data/Final_table_cleaned.csv
A conflict report CSV: tb_recurrence_analysis/data/sinan_conflicts.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
CSV_IN  = (
    "/Users/jasonandrews/Library/CloudStorage/"
    "GoogleDrive-jasonandr@gmail.com/My Drive/"
    "TB Recurrence Predictors/Original data /Final table.csv"
)
OUT_DIR = Path("/Users/jasonandrews/.gemini/antigravity/scratch/"
               "tb_recurrence_analysis/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_OUT      = OUT_DIR / "Final_table_cleaned.csv"
CONFLICT_OUT = OUT_DIR / "sinan_conflicts.csv"

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("Loading CSV (this may take a moment)…")
df = pd.read_csv(CSV_IN, low_memory=False, dtype={"sinan": str})
print(f"  Rows: {len(df):,}  |  Columns: {df.shape[1]}")

# ── 2. Zero-pad sinan to 7 digits ─────────────────────────────────────────────
df["sinan_original"] = df["sinan"]          # keep original for audit
df["sinan_padded"]   = df["sinan"].str.strip().str.zfill(7)

n_before = df["sinan"].nunique()
n_after  = df["sinan_padded"].nunique()
print(f"\nUnique sinan before padding : {n_before:,}")
print(f"Unique sinan after  padding : {n_after:,}")
print(f"  (Records affected by padding: "
      f"{(df['sinan'] != df['sinan_padded']).sum():,} rows)")

# ── 3. Normalise dob and sex for identity matching ────────────────────────────
df["dob_parsed"] = pd.to_datetime(df["dob"], errors="coerce")
df["sex_norm"]   = df["sex"].str.strip().str.upper()

# ── 4. Detect conflicts within each padded sinan ──────────────────────────────
#
# A "conflict" = same sinan_padded, but records disagree on dob or sex.
# These are probably distinct people who share the numeric ID after padding.

def identity_key(row):
    """String key combining dob + sex; NaT/NaN treated as 'UNKNOWN'."""
    dob = row["dob_parsed"].strftime("%Y-%m-%d") if pd.notna(row["dob_parsed"]) else "UNKNOWN_DOB"
    sex = row["sex_norm"] if pd.notna(row["sex_norm"]) else "UNKNOWN_SEX"
    return f"{dob}__{sex}"

print("\nBuilding identity keys…")
df["identity_key"] = df.apply(identity_key, axis=1)

# For each padded sinan, how many distinct identity keys are there?
grp = (
    df.groupby("sinan_padded")["identity_key"]
    .nunique()
    .reset_index()
    .rename(columns={"identity_key": "n_identities"})
)
conflicts = grp[grp["n_identities"] > 1]
print(f"Padded sinan values with identity conflicts: {len(conflicts):,}")

# ── 5. Assign resolved sinan_clean ────────────────────────────────────────────
#
# For non-conflicting sinans: sinan_clean = sinan_padded  (same for everyone
#   in that group — they're the same person even if dob/sex has minor NaN).
#
# For conflicting sinans: sinan_clean = sinan_padded + "_v" + cluster_index
#   where cluster_index ranks the identity sub-groups by earliest notification.

conflict_set = set(conflicts["sinan_padded"])

print("Resolving conflicts and assigning sinan_clean…")

resolved_ids = []

for padded, group in df.groupby("sinan_padded"):
    if padded not in conflict_set:
        # No conflict — everyone in this group is the same person
        resolved_ids.append(
            pd.Series(padded, index=group.index, name="sinan_clean")
        )
    else:
        # Conflict — split by identity key, assign sub-ID
        # Sort keys consistently so v1 = the identity with the most records
        key_counts = group["identity_key"].value_counts()
        key_rank   = {k: i + 1 for i, k in enumerate(key_counts.index)}
        clean_ids  = group["identity_key"].map(
            lambda k, kr=key_rank: f"{padded}_v{kr[k]}"
        )
        resolved_ids.append(clean_ids)

df["sinan_clean"] = pd.concat(resolved_ids).sort_index()

# ── 6. Audit flags ────────────────────────────────────────────────────────────
df["flag_leading_zero_padded"] = df["sinan_original"] != df["sinan_padded"]
df["flag_identity_conflict"]   = df["sinan_padded"].isin(conflict_set)

# ── 7. Save conflict report ───────────────────────────────────────────────────
if len(conflicts) > 0:
    conflict_detail = (
        df[df["flag_identity_conflict"]]
        [["sinan_original", "sinan_padded", "sinan_clean",
          "dob", "dob_parsed", "sex", "identity_key",
          "notification_date", "case_type"]]
        .sort_values(["sinan_padded", "identity_key", "notification_date"])
    )
    conflict_detail.to_csv(CONFLICT_OUT, index=False)
    print(f"\nConflict report saved → {CONFLICT_OUT}")
    print("\nConflict summary (padded sinan | n_rows | identities):")
    for sid in sorted(conflict_set):
        sub = df[df["sinan_padded"] == sid]
        for key, cnt in sub["identity_key"].value_counts().items():
            print(f"  {sid}  |  {cnt} rows  |  {key}  →  sinan_clean={sub.loc[sub['identity_key']==key,'sinan_clean'].iloc[0]}")

# ── 8. Summary statistics ─────────────────────────────────────────────────────
print("\n── Summary ──────────────────────────────────────────────────────────────")
print(f"Total rows                        : {len(df):,}")
print(f"Rows with leading-zero fix        : {df['flag_leading_zero_padded'].sum():,}")
print(f"Rows in conflict groups           : {df['flag_identity_conflict'].sum():,}")
print(f"Unique individuals (sinan_clean)  : {df['sinan_clean'].nunique():,}")
print(f"  vs. unique sinan (raw)          : {df['sinan_original'].nunique():,}")
print(f"  vs. unique sinan (padded only)  : {df['sinan_padded'].nunique():,}")

# ── 9. Save cleaned CSV ───────────────────────────────────────────────────────
# Move sinan_clean to front alongside original columns
cols_front = ["sinan_clean", "sinan_padded", "sinan_original",
              "flag_leading_zero_padded", "flag_identity_conflict"]
other_cols = [c for c in df.columns if c not in cols_front
              and c not in ["sinan", "dob_parsed", "sex_norm",
                             "identity_key", "sinan_padded",
                             "sinan_original", "flag_leading_zero_padded",
                             "flag_identity_conflict"]]
# Keep original sinan col renamed for clarity; drop helper cols
df_out = df[cols_front + other_cols].copy()

print(f"\nSaving cleaned dataset…")
df_out.to_csv(CSV_OUT, index=False)
print(f"Saved → {CSV_OUT}  ({len(df_out):,} rows × {df_out.shape[1]} cols)")
print("\nDone ✓")
