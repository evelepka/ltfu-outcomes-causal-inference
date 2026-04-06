# ITT Cohort Analysis Package: Mortality & Abandonment

This package contains the final, harmonized database and scripts for the analysis of the impact of treatment abandonment (Loss to follow-up) on mortality outcomes.

## 1. Master Database (`ITT_Analysis/data/itt_cohort.csv`)
The database contains 172,463 individuals carefully selected from the SINAN-TB and SIM linked dataset (2013-2023).

### Key Time Variables:
- `best_start`: The proxy treatment start date used as Time Zero. Calculated as: `tx_start` -> `diagnostic_date` -> `notification_date`.
- `end_date`: Date of the index treatment outcome (Abandonment or Cure).
- `death_date`: Date of death (from SIM or SINAN outcome).
- `time_d`: Survival time (years) starting from `end_date` (Landmark view).
- `time_d_tx`: Survival time (years) starting from `best_start` (G-formula view).
- `tx_month_grp`: Duration of the index treatment before its outcome.

### Harmonized Adjusted Variables:
All multivariate models (Cox, G-formula, Landmark) now use these 13 variables consistently:
`age_group, sex, race_clean, edu_clean, hiv_aids, dot_status, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, diabetes`.

## 2. Selection Script (`01_itt_cohort_selection.py`)
This script implements the complex inclusion/exclusion logic:
- Selects the **first** "Novo" episode per individual.
- Filters for Age >= 15 and years 2013-2023.
- **Rigor:** Excludes individuals who died ON or BEFORE the `best_start` to avoid immortal time bias and selection bias.
- Handles proxy dates consistently to avoid dropping patients with missing `tx_start`.

## 3. How to Reuse
To run the analysis with a new database:
1. Ensure the new CSV has the same column names as `Final_table_cleaned.csv`.
2. Run `01_itt_cohort_selection.py` to generate the cohort.
3. Run the subsequent R/Python scripts (02 to 17) in order.

---
*Created on 2026-03-21 to ensure consistency and prevent future errors in cohort formation.*
