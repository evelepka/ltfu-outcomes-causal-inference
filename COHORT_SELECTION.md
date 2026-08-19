# Cohort Selection: Technical Reference

> **Status:** Authoritative. Last refreshed 2026-04-22. Regenerate whenever
> `ITT_Analysis/scripts/01_itt_cohort_selection.py` or `Data/Final_table_cleaned.csv` changes.

This document describes the single analytic cohort used across all causal
analyses (g-formula, MSM, target trial emulation, landmark survival). **There
is one cohort file, not two**; the "abandon-only" and "abandon-vs-non-abandon"
comparisons are both computed from the same file by filtering at read-time on
the `itt_group` column.

---

## 1. Data flow

```
Raw SINAN + SIM linkage (Brazilian TB notification + death registry)
         │
         │  Data/Final_table_cleaned.csv     (98 MB, static input; regenerated
         │                                   by 00_clean_sinan.py at repo root)
         ▼
ITT_Analysis/scripts/01_itt_cohort_selection.py    ← SINGLE SOURCE OF TRUTH
         │
         ├─► ITT_Analysis/data/itt_cohort.csv            (authoritative cohort)
         ├─► Data/exclusion_flowchart.csv                (auto-generated)
         └─► ITT_Analysis/results/Inclusion_Exclusion_ITT.docx
```

All ≈46 downstream R and Python scripts read only from `itt_cohort.csv`. No
script writes it except `01_itt_cohort_selection.py`.

### File locations

| File | Path |
|---|---|
| Input | `$PROJECT_ROOT/Data/Final_table_cleaned.csv` |
| Output (cohort) | `$PROJECT_ROOT/ITT_Analysis/data/itt_cohort.csv` |
| Output (flowchart) | `$PROJECT_ROOT/Data/exclusion_flowchart.csv` |
| Output (docx) | `$PROJECT_ROOT/ITT_Analysis/results/Inclusion_Exclusion_ITT.docx` |

`$PROJECT_ROOT` resolves in this order (see `_paths.py` / `_paths.R`):
1. `TB_ABANDONMENT_ROOT` environment variable
2. `~/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/TB SP 2026/LTFU Paper`
3. `~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper`
4. `~/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper`

A candidate only counts if it actually contains `ITT_Analysis/data`. There is
**no** script-relative fallback: if nothing resolves, the scripts raise. The old
fallback resolved `PROJECT_ROOT` to `$HOME` on an unmounted machine and every
downstream path silently became `~/ITT_Analysis/...` (ADR-0002).

---

## 2. Unit of observation

One row per **individual** (`sinan_clean`), based on their **first "Novo" (new
case) TB episode** with a valid outcome code.

- `sinan_clean` is the cleaned SINAN identifier produced by `00_clean_sinan.py`:
  zero-padded, with DOB+sex conflicts resolved via `_v1`/`_v2` suffixes.
- "Novo" case_type selects new cases, excluding retreatment / recidiva / missing.
- Valid outcome codes exclude `""`, NA, and `"Mud Diag"` (diagnosis change).
- Among qualifying episodes per individual, the earliest by `end_date` is
  retained.

---

## 3. Inclusion / exclusion flow

Numbers below are auto-written to `Data/exclusion_flowchart.csv` every time the
script runs; treat that CSV as authoritative. They are pinned in
`test/golden/prose_numbers.yaml` and re-checked by `bash test/run_fast.sh`, so
this table cannot silently fall out of date again (it did, for three months and
across two exclusion-criteria changes — see `docs/dead-ends.md`).

| Step | N remaining | Excluded |
|---|---:|---:|
| Initial: First Novo episode per individual (transfers already removed) | 198,409 | — |
| Exclude: Age < 15 years | 192,231 | 6,178 |
| Exclude: Treatment end date outside 2013–2023 | 174,354 | 17,877 |
| Exclude: No recorded treatment start date (Abandono Primario) | 171,576 | 2,778 |
| Exclude: Death on or before treatment start/proxy (immortal-time guard) | 171,070 | 506 |
| Exclude: Invalid dates (end_date < tx_start) | 171,069 | 1 |
| Exclude: Post-mortem identification (discovered after death) | 171,048 | 21 |
| **Final cohort N = 171,048** | | |

Composition:

- **Loss to follow-up (LTFU / abandoned): 20,830** — `case_outcome ∈ {Abandono, Abandono Primario, Faltoso}`
- **Non-LTFU (maintained care): 150,218** — any other valid definitive outcome
- **Excluded: transfers (`Transf Outro Municipio`, `Transf Outro Estado/Pais`)** — dropped before counting, not in the flowchart

---

## 4. Exposure definition

Column: `itt_group` (two levels after transfer exclusion)

| Value | Index episode `case_outcome` |
|---|---|
| `"Loss to follow-up"` | `Abandono`, `Abandono Primario`, `Faltoso` |
| `"Non-LTFU"` | Any other valid outcome (cure, treatment completion, death in care, etc.) |

Note that "Faltoso" (default loss to follow-up) is grouped with "Abandono"
(formal abandonment). This is intentional — both represent treatment
disengagement.

---

## 5. Deriving the two comparison cohorts

Both derive from `itt_cohort.csv` by filtering at read-time.

### 5a. Abandon-vs-non-abandon (main causal comparison)
Use the full file. Models regress outcomes on `itt_group` as a two-level
factor. Reference level = `"Non-LTFU"`. Example:

```r
df <- read.csv("ITT_Analysis/data/itt_cohort.csv")
df$itt_group <- factor(df$itt_group, levels = c("Non-LTFU", "Loss to follow-up"))
```

Used by: `05_itt_g_formula_analysis.R`, `14_itt_landmark_mortality_analysis.R`,
`18_itt_time_dependent_cox.R`, `24_itt_ltfu_msm_ipw.R`, target-trial scripts
(30–32, 41), and all plot scripts (33–41).

### 5b. Abandoners-only (LTFU subgroup analyses)
Filter after read. Example:

```r
df <- read.csv("ITT_Analysis/data/itt_cohort.csv") |>
      dplyr::filter(itt_group == "Loss to follow-up")
```

Used by: `03_itt_multiple_imputation_models.R`, `13_itt_ltfu_outcomes_descriptive.R`,
`15_ltfu_stratified_curves.R`, `16_ltfu_stratified_homelessness.R`,
`17_ltfu_stratified_timing.R`, `19_itt_ltfu_alluvial_plot.R`,
`25_itt_ltfu_gformula_timing.R`, `26_itt_ltfu_counterfactual_returners.R`,
`27`–`29` (abandonment timing/duration/post-dropout).

---

## 6. Time variables — two distinct anchors

The cohort contains **two pairs of time/event variables** for survival analysis.
Choosing the wrong one causes immortal-time bias. **This is the source of the
April 2026 bugfix (commits 8671d40, 92bc70b).**

| Variable pair | Time origin | Use in |
|---|---|---|
| `time_d`, `event_d` | `end_date` (index episode outcome date) | **Landmark / post-treatment** analyses. Script 14 and the LTFU descriptive/stratified scripts. Measures mortality conditional on surviving to the end of the index episode. |
| `time_d_tx`, `time_rn_tx` | `best_start` (treatment initiation proxy) | **G-formula, MSM, target trial** — any causal intent-to-treat estimand. Script 05 and plot scripts 33–41. Anchoring here prevents "phantom early mortality" artifacts in the LTFU arm. |

`best_start` is derived by fallback chain: `tx_start` → `diagnostic_date` →
`notification_date` (first non-null). This rescues cohort members with
missing `tx_start` values.

Any new causal analysis script **must** use `time_d_tx` / `time_rn_tx`.
Any landmark / conditional-on-survival analysis may use `time_d` / `time_rn`.
Don't mix them within a single model.

---

## 7. Outcomes

- **Death**: ascertained from two sources:
  - SIM linkage (`dod` field)
  - SINAN `case_outcome ∈ {Obito TB, Obito NTB}` with `end_date` as death date
  - `death_date_comprehensive` = max of the two sources per individual, censored at 2024-12-31
- **Re-notification** (`event_rn`, `time_rn`): any SINAN notification for the
  same `sinan_clean` after the index `end_date`, censored at 2024-12-31.
- **Competing-risk assignment** (`assign_renotif_itt` / `_tx`):
  - 0 = censored, 1 = re-notification, 2 = death

---

## 8. Covariates (harmonized in script 01)

| Column | Source field | Levels |
|---|---|---|
| `age_group` | `age_tb` | 15–24, 25–44, 45–64, 65+ |
| `sex` | `sex` | Male, Female |
| `race_clean` | `race` | White, Black or Mixed, Other, NA |
| `edu_clean` | `education` | None, ≤7 years, 8–11 years, ≥12 years, NA |
| `hiv_aids` | `hiv` OR `aids` | Positive, Negative, NA |
| `dot_status` | `tx_administration_type` | Yes/No/NA (directly observed therapy) |
| `incarcerated` | `address_type == "DETENTO"` | Yes/No/NA |
| `homelessness` | `address_type == "SEM RESIDENCIA FIXA"` | Yes/No/NA |
| `alcohol`, `drug_use`, `mental_health`, `tobacco_use`, `diabetes`, `other_immuno_condition`, `hosp_admission` | `alcoholism`, `drug_use`, `mental_issue`, etc. | Yes/No/NA |
| `clinical_clean` | `clinical_classif` | Pulmonary / Extrapulmonary / Pulmonary and Extrapulmonary or disseminated |
| `diagnosis_setting` | `disease_discovery` | Outpatient / Emergency / Active finding in institution / Active finding in community / Contact investigation |
| `lab_confirmed_stat` | `lab_confirmed` | Yes/No/NA |
| `bac1_clean` | `bac1` (acid-fast bacilli, smear 1) | Positive / Negative / Not Evaluated |
| `sputum_culture_clean` | `sputum_culture` | Positive / Negative / Not Evaluated |
| `resistance_clean` | `resistance` | Sensitive / Resistant (Any) / Not Evaluated |
| `tx_month_grp` | derived from `best_start → end_date` | <2 months, 2 to <4 months, ≥4 months |

---

## 9. Regenerating the cohort

Any time `Final_table_cleaned.csv` changes, or the cohort script changes,
rerun:

```bash
python3 ITT_Analysis/scripts/01_itt_cohort_selection.py
```

Then run the fast tier, which reports whether the rebuild moved any number this
repo's prose or figures depend on:

```bash
bash test/run_fast.sh
```

Outputs land in the Google Drive paths listed in §1. Commit any script
changes to git.

**Downstream scripts are not auto-rerun.** After regenerating `itt_cohort.csv`,
consider whether to:

1. Re-run `03_itt_multiple_imputation_models.R` (MI dataset is downstream)
2. Re-run `05_itt_g_formula_analysis.R` and target-trial scripts
3. Re-render figures (33–41)

---

## 10. Things that are NOT the authoritative cohort

Explicitly deprecated / do not use:

- `Data/analysis_ready_cohort.csv` — older LTFU-only file (N=21,640) with
  different column naming (`abandon_end`, `index_outcome`, `incarceration`).
  Read only by the `SP-TB-spatial-analyses` repo, which is scheduled as future
  work. Do **not** use for ITT analyses.
- `Data/cohort_with_spatial.csv` — spatial-enriched derivative of
  `analysis_ready_cohort.csv`. Same caveat.
- Anything in `legacy/` — retired pre-ITT crude analysis.

Files previously present that have been deleted (2026-04-22):
`analysis_cohort_21640.csv`, `analysis_ready_cohort_BACKUP.csv`,
`abandonment_cohort.csv`, `cure_survival_cohort.csv`, `eligible_cohort.csv`.

---

## 11. Known caveats / next investigations

- **"21,640" figure elsewhere is NOT the ITT LTFU count.** Older files in
  `Data/` (`analysis_ready_cohort.csv`, the now-deleted
  `analysis_cohort_21640.csv`, the historical `exclusion_flowchart.csv`)
  reference N=21,640. That was the N of an *LTFU-only* cohort constructed
  under pre-ITT logic (the old flowchart started from 38,423 abandonment
  records, not 198,409 Novo episodes). The current ITT cohort has 20,830
  LTFU within a total of 171,048.
- **R script working-directory inconsistency**: Some R scripts read
  `"ITT_Analysis/data/itt_cohort.csv"`, others read
  `"Abandonment Paper/ITT_Analysis/data/itt_cohort.csv"`. Both work if the
  working directory is set appropriately, but this should be harmonized.
