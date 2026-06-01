# Outcomes After TB Treatment Loss to Follow-Up

Causal-inference analysis of mortality after tuberculosis treatment loss to follow-up (LTFU) in São Paulo State, Brazil (2013–2023), linking the state TB surveillance system (TBweb) with the national death registry (SIM).

The headline analysis is a **sequential target-trial emulation** with a symmetric 30-day grace-period landmark, comparing patients who became LTFU at each month of therapy with concurrent on-treatment controls under aligned time at risk.

## Cohort

- **N = 171,069** individuals initiating tuberculosis therapy (2013–2023; first episode; recorded treatment start required).
- **20,830 (12.2%)** experienced loss to follow-up during the index treatment.
- **150,239** had alternative outcomes (cure, on-treatment death, treatment failure, regimen change).
- Median post-LTFU follow-up: 5.0 years (IQR 2.5–8.0).

## Key results

- **Late-mortality penalty of LTFU is causal and persistent across months of treatment.** Sequential target-trial adjusted hazard ratios for late mortality (6–24 months post-landmark) range from **1.83 (95% CI 1.27–2.63)** at month 6 to **2.85 (2.33–3.48)** at month 3.
- **The excess is TB-specific.** Cause-specific HRs for TB-attributable death are 1.61–4.36 across months 1–6; non-TB death HRs are essentially null (0.96–1.48). This negative-control contrast supports a causal effect of interrupted TB therapy rather than residual selection.
- **The mortality surrounding return to care is largely attributable to the antecedent disengagement, not to re-entry.** In a risk-set (incidence-density) matched analysis anchored at the time of return, patients who returned had markedly higher 24-month mortality than matched still-LTFU patients (HR 6.93, 95% CI 5.95–8.08) — an association reflecting confounding by indication, since return is typically prompted by recurrent or progressive disease. A g-computation counterfactual estimated that returning patients, whose observed 24-month mortality was 9.3%, would have experienced only 3.5% had they instead completed treatment — indicating the excess is driven by the disengagement itself. Within-LTFU IPCW and Bayesian-Cox analyses (appendix) give the same qualitative conclusion.
- **Early-window aHRs sit at or below unity** (range 0.63–1.16 across months 1–6) and reflect residual survival/selection biases that the target-trial design cannot fully eliminate, not a protective effect of disengagement.
- **Effect modification:** relative penalty is largest in younger and stably housed individuals (where competing-mortality baseline is low) and smallest in homeless individuals and people aged ≥65 (where baseline mortality is already high).

## Repository layout

```
outcomes-after-tb-abandonment/
├── README.md                  (this file)
├── COHORT_SELECTION.md        cohort definition and inclusion/exclusion details
├── CLAUDE.md                  agent / collaborator guidance
├── .gitignore                 excludes data/, results/, and manuscript drafts
├── figures/
│   ├── Figure_1_descriptive.{png,pdf}        post-LTFU trajectories
│   ├── Figure_2_predictors.{png,pdf}         within-LTFU multivariable Cox + Fine-Gray (mortality + retreatment predictors)
│   ├── Figure_3_return_pathway.{png,pdf}      return-to-care: risk-set matched mortality + g-computation counterfactual
│   ├── Figure_4_causal_mortality.{png,pdf}    target-trial HR(t) + monthly aHR array
│   ├── Figure_5_cause_specific.{png,pdf}      TB vs non-TB cause-specific aHRs
│   └── appendix/Figure_S8_cif_retreatment.{png,pdf}  24-mo cumulative incidence of retreatment (appendix figure)
└── ITT_Analysis/
    ├── README.md              pipeline-level documentation
    ├── Master_Causal_Analysis/ data dictionary (PT)
    └── scripts/
        ├── _paths.R                              shared project-root resolver
        ├── 00_clean_sinan.py                     SINAN ID/date cleaning
        ├── 01_itt_cohort_selection.py            cohort builder
        ├── 02_make_itt_table1.py                 baseline characteristics
        ├── 03a_itt_mi_miceforest.py              multiple imputation
        ├── 03_itt_multiple_imputation_models.R   within-LTFU Cox + Fine-Gray
        ├── ...                                   (full pipeline below)
        └── 67_make_fig3_return_pathway.py        return-to-care figure (Figure 3) builder
```

Raw data files and intermediate CSVs are not committed (patient privacy). All script paths resolve via `_paths.R` and `01_itt_cohort_selection.py`'s root resolver; see `COHORT_SELECTION.md` for the regeneration workflow.

The current manuscript (`Draft_*.docx`), appendix (`Appendix_*.docx`), and markdown source live in the project's shared Google Drive folder, not in this repo. Ask a collaborator for access if you need them.

## Pipeline order

The pipeline runs in numbered phases. Phase numbers correspond to file prefixes.

```
00 → 01                                   data cleaning + cohort selection
01 → 02                                   Table 1
01 → 03a → 03                             within-LTFU multivariable (MI + Cox + Fine-Gray)
01 → 05–28                                descriptive + g-formula + mediation + landmark
01 → 30h, 30d, 30i, 30j, 32d, 32f         target-trial: defn-B (primary), grace, cause-specific, period, subgroups, resistance
01 → 30h_incl_no_tx_sensitivity           tx_start-inclusion sensitivity (appendix Table S1)
01 → 55, 55b, 55c                         baseline composition at landmark (SMD + Love + PS-weighted)
01 → 56, 56b                              within-LTFU return-stratified mortality
01 → 57, 57b, 57c                         IPCW + Bayesian Cox for return-to-care mediation
01 → 58, 58b                              g-computation counterfactual for returners
01 → 59                                   severity-stratified on-treatment mortality
01 → 60                                   competing-risks framing
01 → 61                                   corrected multi-source DR classifier
01 → 50, 50a, 51, 52b, 53, 54, 67         figures 1–5 (52b causal mortality; 67 return-to-care)
```

## Key methodological choices

- **Time origin:** sequential target-trial emulation with a symmetric 30-day grace period at each monthly landmark, addressing immortal-time bias from both treatment-start and treatment-completion anchoring.
- **Treatment-start requirement:** the primary cohort requires a recorded treatment-start date; re-including the 1,394 individuals without one (treatment start imputed from the diagnostic, then notification, date) is a sensitivity analysis showing late-window aHRs within ~6% at every monthly trial (appendix Table S1; `30h_incl_no_tx_sensitivity.R`).
- **Confounding adjustment:** 13 baseline covariates (age, sex, race, education, HIV, DOT, alcohol, drug use, incarceration, homelessness, hospitalisation, clinical form, diabetes), with multiple imputation by miceforest (5 datasets) and Rubin's-rules pooling.
- **Return-to-care:** characterized with (i) a risk-set (incidence-density) matched analysis anchored at the time of return — avoiding the immortal-time bias of comparing returners with non-returners from the time of LTFU; (ii) a g-computation counterfactual of returners' mortality had they completed treatment; and (iii) within-LTFU ITT, naive return-censored, and IPCW estimands that bound persistent-disengagement mortality.
- **Cause-of-death attribution:** hybrid SIM ICD-10 + TBweb `case_outcome` for primary; SIM-only for sensitivity.
- **Drug-resistance:** corrected multi-source classifier (phenotypic SINAN + Xpert/rifasens + isoniazid) in `61_build_dr_status_corrected.py`. The historical `dr_status_lookup.csv` was reclassifying culture-positive-but-DST-untested patients as Sensitive and is superseded.

## Reproducibility

All scripts read from `$PROJECT_ROOT` (Google Drive mount in our setup) and write to `$PROJECT_ROOT/ITT_Analysis/results/`. To reproduce on a new machine:
1. Set the environment variable `TB_ABANDONMENT_ROOT` to your project root.
2. Place `Data/Final_table_cleaned.csv` (cleaned SINAN+SIM linkage) and `Data/death_dates.csv` at the expected locations.
3. Run `01_itt_cohort_selection.py` to build the cohort.
4. Run the phase scripts in numbered order. Most phases are independent given the cohort and the 5 MI datasets.

## Dependencies

**R (≥ 4.5):** `survival`, `cmprsk`, `dplyr`, `mice`, `brms`, `ggplot2`, `patchwork`, `scales`, `broom`

**Python (≥ 3.9):** `pandas`, `numpy`, `scipy`, `statsmodels`, `lifelines`, `miceforest`, `matplotlib`, `python-docx`

## Authors and contact

Evelyn Lepka de Lima, Jason Andrews, and colleagues.

For questions: `jasonandr@stanford.edu`.
