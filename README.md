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
├── ITT_Analysis/                                # Authoritative causal pipeline
│   ├── scripts/
│   │   ├── 01_itt_cohort_selection.py           # Cohort builder → itt_cohort.csv (single source of truth)
│   │   ├── 03_itt_multiple_imputation_models.R  # MI for LTFU subgroup
│   │   ├── 05_itt_g_formula_analysis.R          # Parametric G-computation
│   │   ├── 14_itt_landmark_mortality_analysis.R # Landmark survival (post-treatment)
│   │   ├── 24_itt_ltfu_msm_ipw.R                # Retreatment MSM-IPW
│   │   ├── 30-41_itt_target_trial_*.R           # Sequential & piecewise target trials, RMST, subgroups
│   │   └── 33-37_itt_plot_*.R                   # Main paper figures
│   ├── results/                                 # Walkthrough PDFs, KM datasets, figures
│   └── Master_Causal_Analysis/                  # Data dictionary (PT) + working docs
├── 00_clean_sinan.py                            # SINAN ID cleaning → Final_table_cleaned.csv
├── COHORT_SELECTION.md                          # Technical cohort documentation
├── CLAUDE.md                                    # Guidance for agents/collaborators
├── legacy/                                      # Pre-ITT crude analysis (do not run)
└── README.md                                    # This file
```

Data lives in Google Drive at `~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper/` — never committed to git. See `COHORT_SELECTION.md` for paths and the cohort-regeneration workflow.

## Documentation
- `COHORT_SELECTION.md` — cohort definition, inclusion/exclusion, time variables, and regeneration workflow
- `CLAUDE.md` — agent / collaborator guidance (where data lives, common pitfalls, what not to do)
- `ITT_Analysis/README_FINAL.md` — pipeline-level README

## Abstract Synthesis
For a complete, print-ready synthesis of the advanced causal findings, please see `ITT_Analysis/results/causal_analysis_summary.pdf`.

## Dependencies
- **R**: `survival`, `dplyr`, `broom`, `survivalROC`
- **Python >= 3.9**: `pandas`, `numpy`, `matplotlib`, `python-docx`

## Citation / Contact
Evelyn Lepka de Lima and Jason Andrews
Analysis finalized April 2026.
