# Pipeline reference: ITT_Analysis/scripts/

This is the analysis pipeline for the LTFU-mortality manuscript. Scripts are numbered by phase; later phases generally depend on earlier ones via intermediate CSVs.

## Single source of truth: `data/itt_cohort.csv`

Produced by `01_itt_cohort_selection.py`. Contains 171,069 individuals with the following key columns:

| Column | Description |
|---|---|
| `sinan_clean` | Patient identifier |
| `itt_group` | `Loss to follow-up` or `Non-LTFU` |
| `best_start` | Treatment start (`tx_start`; required) |
| `end_date` | Date of index treatment outcome |
| `death_date` | Date of death (SIM or TBweb) |
| `time_d`, `event_d` | Survival from `end_date` (landmark view) |
| `time_d_tx`, `event_d` | Survival from `best_start` (g-formula view) |
| `time_rn`, `event_rn` | Retreatment indicator |
| `tx_month_grp` | Index treatment duration grouped (<2 mo / 2–4 mo / ≥4 mo) |
| 13 baseline covariates | `age_group`, `sex`, `race_clean`, `edu_clean`, `hiv_aids`, `dot_status`, `alcohol`, `drug_use`, `incarcerated`, `homelessness`, `hosp_admission`, `clinical_clean`, `diabetes` |

All multivariable models use these 13 variables consistently.

## Phase layout

| Phase | Scripts | What they produce |
|---|---|---|
| 00–01 | data cleaning, cohort selection | `itt_cohort.csv` |
| 02 | Table 1 | `Table1_ITT.docx` |
| 03a, 03 | MI + within-LTFU Cox / Fine-Gray | `mi/imp_01..05.csv`, `multivariable_results_mi_cc.csv` |
| 04 | Table 2 (combined multivariable) | `Table2_ITT.docx` |
| 05–12 | g-formula, mediation, retreatment timing, landmark KM | `g_formula_*.csv` |
| 13–19 | within-LTFU descriptive curves, stratified KM, alluvial plot | various |
| 20–22 | manuscript-building helpers, QBA / E-value, comprehensive HTML report | `Analyses_Summary*.html` |
| 23–29 | g-methods triangulation, MSM-IPW, abandonment timing/duration | various |
| 30 series | target-trial pipeline (rolling, MI, early/late, grace, cause-specific, period, defnB, complete-case) | `target_trial_*.csv` |
| 31 | LTFU subgroup interactions | `ltfu_subgroup_effect_modification.csv` |
| 32 series | target-trial subgroup interactions | `target_trial_subgroup_*.csv`, `target_trial_resistance_*.csv`, `target_trial_period_*.csv` |
| 33–41 | HR(t) plots, methodology comparison, time-varying crude, panels, RMST, delayed impact, piecewise | various PNGs/CSVs |
| 50–54 | manuscript figures 1–5 | `Figure_1`–`Figure_5_*.png`/`.pdf` |
| 55–55c | baseline composition at landmark (SMD, Love plot, PS-weighted) | `landmark_*smd*.csv`, `landmark_loveplot.png` |
| 56–56b | within-LTFU return-stratified mortality (all-cause + cause-specific) | `return_stratified_*.csv` |
| 57–57c | IPCW for return-to-care, Bayesian Cox | `ltfu_ipcw_return.csv`, `ipcw_ltfu_vs_ontx_m3.csv`, `bayesian_cox_*.csv` |
| 58–58b | g-computation counterfactual for returners | `counterfactual_returners*.csv` |
| 59 | severity-stratified on-treatment mortality | `ontx_severity_stratified.csv` |
| 60 | competing-risks (death-on-tx vs LTFU) | `competing_risks_by_severity.csv` |
| 60b | manuscript-assembly script (defnB variant) | `Draft_*_v2_defnB.docx` |
| 61 | corrected multi-source DR classifier | `resistance/dr_status_corrected.csv` |

## Outputs

All script outputs land in `$PROJECT_ROOT/ITT_Analysis/results/`. That directory is gitignored.

## Path resolution

R scripts source `_paths.R` which resolves the project root via, in order:
1. `TB_ABANDONMENT_ROOT` environment variable
2. Known Google Drive mounts for current collaborators
3. Script-relative fallback (`../../` from `ITT_Analysis/scripts/`)

Python scripts use `os.path.dirname(__file__)` and walk up to find the project root.

## Data dictionary

`Master_Causal_Analysis/` contains the Portuguese variable dictionary for the source SINAN-TBweb fields.
