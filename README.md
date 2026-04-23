# Outcomes After TB Treatment Abandonment

Causal epidemiological analysis of retreatment and long-term mortality following tuberculosis (TB) treatment abandonment in São Paulo, Brazil (2013–2023).

## Research Question

Historically, evaluating the mortality impact of TB treatment abandonment is confounded by immortal time bias (late abandoners appear falsely protected) and symptomatic presentation bias (returning to care falsely appears harmful). This repository utilizes advanced causal inference modeling to answer:
1. What is the true mortality penalty of abandoning TB therapy over the 6-month treatment course? 
2. Does returning to therapy rescue patients, or does it serve as a proxy for clinical deterioration?
3. Which patient subpopulations disproportionately suffer the relative mortality penalties of abandonment?

## Key Findings (Causal Analysis Phase)

**Overall ITT cohort:** 172,463 individuals successfully initiating therapy, of which 21,619 (12.5%) abandoned treatment.
- **The Paradox of Retreatment:** Among those who abandoned, 43.8% re-entered treatment. Using Marginal Structural Models (MSM) with Inverse Probability Weighting (IPW), returning to care carried an **overwhelming 5-fold increase in mortality (aHR ~5.0)**. Returning to care functions not as a protective behavior, but as a distress beacon for severe, acute disease progression.
- **Sequential Target Trials:** To eliminate immortal time bias, we compared abandoners explicitly against mathematically matched patients who remained *actively on therapy*. Patients abandoning as late as **Month 4 still suffered a 2.26-fold relative mortality risk**, proving late abandonment is catastrophic.
- **Competing Risks (Effect Modification):** Subgroup interaction models revealed the mathematical penalty for abandoning is actually *worse* for healthy populations (age 15-24 HR 4.14 vs. age 45-64 HR 1.92). Marginalized populations (homeless, HIV+) suffer such overwhelming absolute baseline mortality irrespective of treatment that their relative penalty for dropping out is mechanically blunted.

## Repository Structure

```
├── ITT_Analysis/                                   # Authoritative causal pipeline
│   ├── scripts/
│   │   ├── _paths.R                                # Shared project-root resolver for R
│   │   ├── 01_itt_cohort_selection.py              # Cohort builder → itt_cohort.csv (single source of truth)
│   │   ├── 02_make_itt_table1.py                   # Table 1 (LTFU vs Non-LTFU baseline)
│   │   ├── 03a_itt_mi_miceforest.py                # Multiple-imputation datasets (miceforest)
│   │   ├── 03_itt_multiple_imputation_models.R     # MI Cox + Fine–Gray pooled via Rubin's rules
│   │   ├── 04_itt_make_table2.py                   # Table 2 (multivariable MI + complete-case)
│   │   ├── 30c_itt_target_trial_early_late.R       # Sequential target-trial HR array (early/late split)
│   │   ├── 32b_itt_target_trial_subgroups_mi.R     # MI-pooled subgroup target-trial interactions
│   │   ├── 50_make_fig1_descriptive.R              # Figure 1 — descriptive panel
│   │   ├── 50a_alluvial_three_column.py            # Figure 1A — 3-column Sankey alluvial
│   │   ├── 51_make_fig2_stratified_retreatment.R   # Figure 2 — stratified CIF of retreatment
│   │   ├── 52_make_fig3_causal_mortality.R         # Figure 3 — causal mortality panels
│   │   ├── make_results_html.py                    # Self-contained results HTML (drops a copy on Desktop)
│   │   └── archive/                                # Superseded scripts from earlier pipeline iterations
│   ├── Master_Causal_Analysis/                     # Data dictionary (PT)
│   └── README_FINAL.md                             # Pipeline README
├── figures/                                        # Current paper figures (Figure_1/2/3 + alluvial)
├── 00_clean_sinan.py                               # SINAN ID cleaning → Final_table_cleaned.csv
├── COHORT_SELECTION.md                             # Technical cohort documentation
├── CLAUDE.md                                       # Guidance for agents/collaborators
├── legacy/                                         # Pre-ITT crude analysis (do not run)
└── README.md                                       # This file
```

All data, intermediate CSVs, and regenerated figures live in Google Drive at `~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper/` — never committed to git. Scripts resolve this path automatically (see `_paths.R` / `01_itt_cohort_selection.py`). See `COHORT_SELECTION.md` for the regeneration workflow.

## Pipeline order

```
01 → 03a → 03 → 04                (tables + MI models)
01 → 30c, 32b                     (target-trial CSVs for Figure 3)
01 → 50a → 50, 51, 52             (figures)
→ make_results_html.py            (bundles everything into one HTML)
```

Running `make_results_html.py` produces `TB_Abandonment_Results_YYYY-MM-DD.html` on the Desktop — a self-contained report with base64-embedded figures.

## Documentation
- `COHORT_SELECTION.md` — cohort definition, inclusion/exclusion, time variables, and regeneration workflow
- `CLAUDE.md` — agent / collaborator guidance (where data lives, common pitfalls, what not to do)
- `ITT_Analysis/README_FINAL.md` — pipeline-level README

## Dependencies
- **R**: `survival`, `dplyr`, `broom`, `survivalROC`
- **Python >= 3.9**: `pandas`, `numpy`, `matplotlib`, `python-docx`

## Citation / Contact
Evelyn Lepka de Lima and Jason Andrews
Analysis finalized April 2026.
