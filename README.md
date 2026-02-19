# Outcomes After TB Treatment Abandonment

Survival analysis of re-notification and mortality risk following tuberculosis (TB) treatment abandonment in Brazil, using SINAN surveillance data (2013–2024).

## Research Question

Among individuals who abandoned TB treatment, what is the cumulative risk of:
1. **Re-entering the TB notification system** (re-notification / retreatment)?
2. **Dying** — at 1, 2, and 3 years after abandonment?

And how do these risks vary by age, sex, HIV/AIDS status, incarceration, homelessness, and the treatment month at which abandonment occurred?

## Key Findings

**Overall cohort: N = 29,784 individuals with a first abandonment episode**

| Outcome | 1 year | 2 years | 3 years |
|---|---|---|---|
| Re-notification | 35.2% | 40.1% | 42.6% |
| Death (comprehensive) | 3.2% | 5.1% | 6.4% |

**Highest-risk subgroups:**
- **AIDS patients:** 49% re-notification and 16.5% death at 3 years
- **Homeless individuals:** 47% re-notification and 10.7% death at 3 years  
- **Month 2 abandonees** have the highest re-notification risk (~61% at 3yr)
- **Month 1 abandonees** have the highest death risk (6.6% at 1yr)

## Methods Summary

- **Index event:** First abandonment episode per individual (`case_outcome` = "Abandono" / "Abandono Primario")
- **Time zero:** `end_date` (treatment outcome date)
- **Administrative censoring:** December 31, 2024
- **Re-notification estimator:** Aalen-Johansen CIF (death is a competing event)
- **Death estimator:** 1 − Kaplan-Meier (re-notification does NOT censor; individuals followed through subsequent episodes)
- **Death ascertainment:** `dod` field + `case_outcome` ∈ {"Obito TB", "Obito NTB"} — combined for ~4.6× more deaths than `dod` alone

See [`METHODS.md`](METHODS.md) for full reproducible technical description.

## Repository Structure

```
├── 00_clean_sinan.py                    # Step 1: SINAN ID cleaning & zero-padding
├── 04_abandonment_full_analysis.py      # Step 2+: Main analysis (all figures & tables)
├── METHODS.md                           # Full technical methods document
├── figures/                             # All output figures (PNG)
│   ├── cif_renotification_overall_v2.png
│   ├── cif_death_overall_v2.png
│   ├── cif_renotification_by_age_group_v2.png
│   ├── cif_death_by_age_group_v2.png
│   ├── cif_renotification_by_hiv_grp_v2.png
│   ├── cif_death_by_hiv_grp_v2.png
│   ├── cif_renotification_by_incarcerated_v2.png
│   ├── cif_death_by_incarcerated_v2.png
│   ├── cif_renotification_by_homeless_v2.png
│   ├── cif_death_by_homeless_v2.png
│   ├── cif_renotification_by_sex_grp_v2.png
│   ├── cif_death_by_sex_grp_v2.png
│   ├── cif_renotification_by_tx_month_v2.png
│   ├── cif_death_by_tx_month_v2.png
│   └── tx_month_distribution_v2.png
└── results/                             # Risk tables (CSV)
    ├── abandonment_risk_table_v3.csv    # Risks by demographic/clinical strata
    └── abandonment_risk_table_txmonth_v3.csv  # Risks by treatment month
```

## Data

The raw data (`Final table.csv`) is from SINAN (Brazil's TB notification system) and is **not included** in this repository due to patient privacy. The cleaned cohort files are also excluded.  

To reproduce the analysis, place the raw CSV at the path specified in `00_clean_sinan.py` and run the scripts in order.

## Exploratory / Superseded Scripts

These earlier scripts are retained for transparency but are superseded by `04_abandonment_full_analysis.py`:
- `01_abandonment_survival.py` — initial analysis (death as competing risk only)
- `02_abandonment_stratified.py` — v1 stratified analysis (dod field only)
- `03_abandonment_by_tx_month.py` — tx-month analysis (dod field only)

## Dependencies

```
python >= 3.9
pandas
numpy
matplotlib
```

No external survival analysis packages are required — all estimators are implemented from scratch.

## Citation / Contact

Jason Andrews Lab, Stanford University  
Analysis performed February 2026.
