# Technical Methods: Outcomes After TB Treatment Abandonment

## Overview

This document describes in full reproducible detail the analysis of outcomes (re-notification and death) following tuberculosis (TB) treatment abandonment in a Brazilian TB surveillance dataset (SINAN). The analysis aims to estimate the cumulative risk of re-entering the TB notification system and the risk of death at 1, 2, and 3 years after abandonment, overall and stratified by key clinical and social risk factors.

---

## Data Source

**Dataset:** Final Table — SINAN (Sistema de Informação de Agravos de Notificação), Brazil  
**Location:** `My Drive / TB Recurrence Predictors / Original data / Final table.csv`  
**Format:** CSV (~88 MB) and Excel (.xlsx)  
**Dimensions:** 270,735 rows × 54 columns  
**Unit of observation:** Each row is a TB treatment notification (episode). Individuals may appear multiple times.

### Key Variables Used

| Variable | Description |
|---|---|
| `sinan` | Individual identifier (7-digit numeric; leading zeros were dropped in export) |
| `dob` | Date of birth (text, e.g. "June 22, 1972") |
| `sex` | Sex (M/F) |
| `notification_date` | Date notification was entered in the system |
| `tx_start` | Treatment start date |
| `end_date` | Date of treatment outcome |
| `case_outcome` | Treatment outcome (e.g. "Cura", "Abandono", "Obito TB") |
| `case_type` | Episode type (e.g. "Novo", "Retr Aband", "Recidiva") |
| `dod` | Date of death (where recorded) |
| `hiv` | HIV test result ("Pos", "Neg", "N/realiz", "S/inf") |
| `aids` | AIDS diagnosis ("S"=yes, "N"=no) |
| `address_type` | Address category ("ENDERECO PADRAO", "DETENTO", "SEM RESIDENCIA FIXA") |
| `age_tb` | Age at time of TB notification |
| `lab_confirmed` | Laboratory confirmation (1=yes, 0=no) |

---

## Step 1: Data Cleaning — SINAN ID Standardisation

**Script:** `00_clean_sinan.py`

### Problem
The `sinan` identifier is nominally a 7-digit number, but leading zeros were stripped during data export (e.g., individual `0000144` may appear as `"144"`). This creates apparent duplicates.

### Solution
1. **Zero-padding:** All `sinan` values are left-padded with zeros to exactly 7 digits (`str.zfill(7)`), creating `sinan_padded`.
2. **Identity verification:** To confirm two records belong to the same individual, `date_of_birth` (parsed to datetime) AND `sex` must match within each `sinan_padded` group.
3. **Conflict resolution:** 17 `sinan_padded` values mapped to 2 genuinely different individuals (different DOB/sex). These were split into sub-IDs (e.g. `0047624_v1`, `0047624_v2`), ranked by number of records.
4. **Output column:** `sinan_clean` — the definitive individual identifier used in all downstream analyses.

### Results
- Rows affected by zero-padding: **24,601**
- Identity conflicts resolved: **17 sinan IDs → 34 distinct individuals**
- Final unique individuals: **235,629**

### Audit Columns Added
- `sinan_original` — raw value from source file
- `sinan_padded` — zero-padded version
- `flag_leading_zero_padded` — boolean: was padding applied?
- `flag_identity_conflict` — boolean: was this row in a conflict group?

---

## Step 2: Cohort Definition — Abandonment Index Events

**Script:** `04_abandonment_full_analysis.py`

### Inclusion Criteria
- `case_outcome` ∈ {"Abandono", "Abandono Primario"}
- `end_date` is not missing (required as time zero)
- `end_date` ≤ December 31, 2024 (administrative censoring date)
- Only the **first** abandonment episode per individual is used as the index event (most interpretable epidemiologically)

### Time Zero
`end_date` of the index abandonment episode.

### Final Cohort
**N = 29,784 individuals**  
(1,092 excluded from tx-month analysis due to missing/implausible `tx_start`)

---

## Step 3: Outcome Definitions

### Outcome 1 — Re-notification
**Definition:** The earliest `notification_date` for the same individual (`sinan_clean`) that occurs **strictly after** the abandonment `end_date`.

**Estimator:** Aalen-Johansen Cumulative Incidence Function (CIF), treating death before re-notification as a **competing event**. This prevents overestimating re-notification risk in the presence of mortality.

- event = 1: re-notification
- event = 2 (competing): death before re-notification
- event = 0: censored at December 31, 2024

### Outcome 2 — Death (Comprehensive)
**Definition:** Death ascertained from TWO sources combined:

1. **`dod` field** — date of death recorded on any row for that individual (minimum across all rows)
2. **`case_outcome` ∈ {"Obito TB", "Obito NTB"}** — the `end_date` of any such episode for that individual (minimum across all rows)

The earliest date across both sources is used as the death date. Only deaths occurring **strictly after** the abandonment `end_date` are counted.

**Critical methodological decision:** Re-notification does **NOT** censor for death. Individuals are followed until death or December 31, 2024, regardless of whether they were re-notified. This captures deaths that occur after a subsequent treatment episode.

**Estimator:** Kaplan-Meier (1 − KM), i.e. standard cause-specific cumulative incidence where re-notification is treated as an independent (non-censoring) event.

### Why This Matters
Previous analyses using `dod` alone captured only **496 deaths (1.7%)** in the cohort. Adding `case_outcome` Obito episodes increased this to **2,269 deaths (7.6%)**  — a 4.6× improvement in death ascertainment.

---

## Step 4: Statistical Analysis

All analyses use:
- **Administrative censoring date:** December 31, 2024
- **Time scale:** Years from abandonment `end_date`
- **Risk time points reported:** 1, 2, and 3 years

### Aalen-Johansen Estimator (Re-notification CIF)

For a two-cause competing risks setting (cause 1 = re-notification, cause 2 = death):

$$\hat{F}_1(t) = \sum_{t_j \leq t} \hat{S}(t_{j-1}) \cdot \frac{d_{1j}}{n_j}$$

where:
- $\hat{S}(t_{j-1})$ = overall event-free survival just before time $t_j$
- $d_{1j}$ = number of re-notification events at time $t_j$
- $n_j$ = number at risk just before $t_j$

### Kaplan-Meier (Death, 1−KM)

$$\hat{F}_2(t) = 1 - \prod_{t_j \leq t} \left(1 - \frac{d_{2j}}{n_j}\right)$$

where $d_{2j}$ = deaths at $t_j$ and re-notification events are **not** treated as censoring events.

### Implementation Note
Both estimators are implemented from scratch in Python (NumPy/Pandas) without external survival analysis libraries, to ensure full reproducibility and transparency.

---

## Step 5: Stratified Analyses

### Stratification Variables

| Variable | Groups | Source column(s) |
|---|---|---|
| Age group | <18, 18–34, 35–49, 50–64, 65+, Unknown | `age_tb` at index episode |
| Sex | Male, Female, Unknown | `sex` at index episode |
| HIV/AIDS status | AIDS, HIV+, HIV−, Unknown/not tested | `aids` (S=AIDS), `hiv` |
| Incarceration | Incarcerated, Not incarcerated | `address_type` == "DETENTO" |
| Homelessness | Homeless, Not homeless | `address_type` contains "SEM RESIDENCIA" |
| Treatment month | Month 1–6, Month 7+ | `ceil((end_date − tx_start) / 30.44 days)` |

### Treatment Month at Abandonment
For each individual with a valid `tx_start`, the number of months on treatment before abandonment is computed as:

```
tx_month = ceil((abandon_end - tx_start).days / 30.44)
```

Groups: Month 1, Month 2, Month 3, Month 4, Month 5, Month 6, Month 7+ (capped to maintain sample sizes). Individuals with `tx_days ≤ 0` or missing `tx_start` are excluded from this analysis (n = 1,092 excluded; analytic n = 28,692).

---

## Step 6: Outputs

### Figures (in `figures/`)

| File | Description |
|---|---|
| `cif_renotification_overall_v2.png` | Overall re-notification CIF |
| `cif_death_overall_v2.png` | Overall death CIF (comprehensive) |
| `cif_renotification_by_age_group_v2.png` | Re-notification by age group |
| `cif_death_by_age_group_v2.png` | Death by age group |
| `cif_renotification_by_sex_grp_v2.png` | Re-notification by sex |
| `cif_death_by_sex_grp_v2.png` | Death by sex |
| `cif_renotification_by_hiv_grp_v2.png` | Re-notification by HIV/AIDS status |
| `cif_death_by_hiv_grp_v2.png` | Death by HIV/AIDS status |
| `cif_renotification_by_incarcerated_v2.png` | Re-notification by incarceration |
| `cif_death_by_incarcerated_v2.png` | Death by incarceration |
| `cif_renotification_by_homeless_v2.png` | Re-notification by housing status |
| `cif_death_by_homeless_v2.png` | Death by housing status |
| `tx_month_distribution_v2.png` | Distribution of abandonment by tx month + 1yr bar chart |
| `cif_renotification_by_tx_month_v2.png` | Re-notification by treatment month |
| `cif_death_by_tx_month_v2.png` | Death by treatment month |

### Tables (in `results/`)

| File | Description |
|---|---|
| `abandonment_risk_table_v3.csv` | 1/2/3-yr risks by all demographic/clinical strata |
| `abandonment_risk_table_txmonth_v3.csv` | 1/2/3-yr risks by treatment month group |
| `sinan_conflicts.csv` | Rows with identity conflicts resolved in cleaning |

---

## Key Results Summary

### Overall Cohort (N = 29,784)

| Outcome | 1 year | 2 years | 3 years |
|---|---|---|---|
| Re-notification | 35.2% | 40.1% | 42.6% |
| Death (comprehensive) | 3.2% | 5.1% | 6.4% |

### Notable Subgroup Findings

**Re-notification highest in:**
- AIDS patients: 49.0% at 1 yr, 57.7% at 3 yr
- Homeless individuals: 39.5% at 1 yr, 47.1% at 3 yr
- Month 2 abandonees: 51.1% at 1 yr, 60.8% at 3 yr (peak; early enough to still be symptomatic, not yet cured)

**Death highest in:**
- AIDS patients: 9.2% at 1 yr, 16.5% at 3 yr
- HIV+ (non-AIDS): 7.2% at 1 yr, 12.5% at 3 yr
- Homeless: 5.6% at 1 yr, 10.7% at 3 yr
- 65+: 5.9% at 1 yr, 12.0% at 3 yr
- Month 1 abandonees: 6.6% at 1 yr (highest by tx month; reflects sickest early abandoners)

**Incarcerated individuals** show paradoxically lower re-notification (27.7% at 1 yr) and death (0.9% at 1yr) — likely reflecting release from the prison system and loss to the SINAN database rather than genuinely better outcomes.

---

## Reproducibility Notes

### Environment
- Python 3.x
- Libraries: `pandas`, `numpy`, `matplotlib`
- No external survival analysis packages required (estimators implemented from scratch)

### Running the Analysis
Scripts should be run in order from the project root directory:

```bash
python3 00_clean_sinan.py          # Clean and standardise SINAN IDs
python3 04_abandonment_full_analysis.py  # Main analysis (all stratified figures + tables)
```

### Data Path
The raw data file must be accessible at:
```
/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/
  My Drive/TB Recurrence Predictors/Original data/Final table.csv
```

Update `CSV_CLEANED` in `00_clean_sinan.py` and `CSV_IN` in `04_abandonment_full_analysis.py` if the path differs.

### Intermediate Scripts (exploratory / earlier versions)
- `01_abandonment_survival.py` — initial analysis (death as competing risk only; superseded)
- `02_abandonment_stratified.py` — v1 stratified analysis (dod only; superseded)
- `03_abandonment_by_tx_month.py` — tx-month analysis (dod only; superseded)
- `04_abandonment_full_analysis.py` — **final, definitive analysis** (comprehensive death, all strata)

---

## Limitations

1. **Death undercounting remains possible.** SINAN is a notification database, not a mortality registry. Deaths outside the TB system (e.g., from other causes in the community, not captured as Obito outcome) may be missed.
2. **Administrative censoring assumes non-informative censoring.** Individuals lost before December 31, 2024 for non-random reasons (e.g., emigration) are treated as censored.
3. **Index episode only.** Only the first abandonment per individual is used. Repeat-abandoners may have different outcome profiles.
4. **Covariates from index episode.** HIV status, address type, etc. are taken from the index abandonment episode and may change over time.
5. **`tx_start` quality.** 1,092 individuals (~3.7%) are excluded from the treatment-month analysis due to missing or implausible `tx_start` values (tx_days ≤ 0).
