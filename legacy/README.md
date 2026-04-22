# Legacy reference files

Preserved for historical reference. **Not part of the active pipeline** — do not run as-is and do not cite in new analyses.

The current authoritative pipeline is `ITT_Analysis/`. See `COHORT_SELECTION.md` in the repo root for the cohort definition.

## Contents

### From predecessor `TB-Abandonment-Analysis` repo
- `MASTER_RUN_ANALYSIS.py` — orchestrator for the old (pre-ITT) analysis pipeline
- `config.py` / `config.R` — path configuration pattern; superseded by the env-var + candidate-path resolution in `ITT_Analysis/scripts/01_itt_cohort_selection.py`
- `setup_sherlock_env.sh` — environment setup for Stanford Sherlock cluster runs
- `README_abandonment.md` — initial project README

### Crude survival analysis (pre-ITT, pre-target-trial)
Did not account for immortal time bias or confounding. Superseded by `ITT_Analysis/scripts/` (g-formula, MSM, target trial emulation).

- `01_abandonment_survival.py` — crude abandonment survival (wrote `abandonment_cohort.csv`, now deleted)
- `02_abandonment_stratified.py` — stratified crude analysis
- `03_abandonment_by_tx_month.py` — treatment-month stratified survival
- `04_abandonment_full_analysis.py` — "final" crude analysis before the causal rewrite
- `METHODS.md` — documentation of the crude analysis methodology
