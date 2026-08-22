# Guidance for agents working on this repo

Read this first. It encodes lessons learned from cleanup (2026-04-22) that
cost the project non-trivial work.

## What this repo is

Causal epidemiological analysis of outcomes (retreatment + mortality) after TB
treatment abandonment in São Paulo, Brazil (2013–2023). Uses SINAN-SIM linked
data from Brazilian public-health surveillance.

- **Authoritative pipeline:** `ITT_Analysis/scripts/` (numbered 01–41, plus
  `tmp_*.R` scratch files). R + Python. Includes g-formula, MSM-IPW, target
  trial emulation, landmark analyses, piecewise trials, RMST.
- **Data cleaning upstream:** top-level `00_clean_sinan.py` produces
  `Data/Final_table_cleaned.csv` (the input to the cohort selection).
- **Legacy (do not run):** `legacy/` — pre-ITT crude survival analysis and
  the predecessor project's orchestrator/configs.

## Where data lives

Data is **never** committed to git. All CSVs, shapefiles, and figures live in
Google Drive at:

```
~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper/
```

The cohort-selection script resolves this path automatically (env var
`TB_ABANDONMENT_ROOT` overrides; otherwise tries known GDrive mounts). See
`ITT_Analysis/scripts/01_itt_cohort_selection.py:_find_project_root`.

A collaborator (Evelyn Lepka de Lima) has her own GDrive mount at
`/Users/evelynlepkadelima/...`; the resolver tries that too.

## The cohort situation

**There is ONE cohort file, not two.** `ITT_Analysis/data/itt_cohort.csv`
contains both `"Loss to follow-up"` and `"Non-LTFU"` individuals (N=171,048).
The "abandoners-only" cohort is `filter(itt_group == "Loss to follow-up")`
applied at read-time. See `COHORT_SELECTION.md` for full detail.

Do not create new cohort files. Do not cache subsets to disk. Filter at
read-time.

## Time variable trap (read this before touching survival code)

The cohort has **two time/event pairs** with different origins. Picking
the wrong one causes immortal-time bias.

- `time_d_tx`, `time_rn_tx` — origin at treatment start (`best_start`).
  **Use for g-formula / MSM / target trial / any causal contrast.**
- `time_d`, `time_rn` — origin at `end_date`.
  **Use for landmark analyses conditional on surviving the index episode.**

This was the subject of the April 2026 bugfix (commits 8671d40, 92bc70b).
Figure 1 was re-rendered because the LTFU arm looked falsely protected when
early dropouts' retrospective exposure time was counted as "at risk under
the treatment arm".

## Paths in R vs Python

- **Python** scripts in `ITT_Analysis/scripts/` use the project-root resolver
  in script 01 as a model. Absolute-ish paths.
- **R** scripts currently use relative paths, inconsistently:
  - some read `"ITT_Analysis/data/itt_cohort.csv"`
  - some read `"Abandonment Paper/ITT_Analysis/data/itt_cohort.csv"`
  Both work if `setwd()` is right. Harmonizing is on the todo list but not
  urgent. Don't casually rewrite them.

## When asked to run the pipeline

1. The cohort script reads from GDrive and writes to GDrive. If GDrive isn't
   mounted locally, fail loudly — don't synthesize substitute paths.
2. Downstream R scripts require the cohort CSV to exist. If you re-run cohort
   selection, MI (`03_itt_multiple_imputation_models.R`) and any downstream
   analyses may need re-running too. Don't silently assume they're current.
3. Plots 33–41 depend on a lot of upstream results. Rerun the chain, not just
   the plot.

## Do NOT do these things

- **Don't create new cohort CSVs.** If you need a subset, filter on read.
  History: five dead cohort files piled up in Data/ (deleted 2026-04-22).
- **Don't read `analysis_ready_cohort.csv` from ITT scripts.** It's a legacy
  LTFU-only file used only by the separate `SP-TB-spatial-analyses` repo.
  Its column naming (`abandon_end`, `index_outcome`, `incarceration`) differs
  from the ITT cohort and WILL produce wrong results if you wire it into an
  ITT model.
- **Don't commit data to git.** `.gitignore` should already block CSVs, but
  the `Data/` directory itself is outside the repo anyway.
- **Don't move scripts without grepping for references first** (`METHODS.md`
  and `README.md` reference 01–04 by name, etc.)
- **Don't touch the spatial repo (`SP-TB-spatial-analyses`) as part of ITT
  work.** It's deliberately a separate project.

## Related repos

- `github.com/evelepka/ltfu-outcomes-causal-inference` — this repo (ITT
  pipeline, **canonical**). It is the `origin` remote and `main` is the
  authoritative branch. Push here; `git pull`/`git push` with no arguments
  target it (local `main` tracks `origin/main`).
- `github.com/jasonandr/outcomes-after-tb-abandonment` — DEPRECATED backup.
  No longer a remote on this clone (removed 2026-06-02); was stale on an old
  `master` branch. Do not treat as authoritative.
- `github.com/jasonandr/SP-TB-spatial-analyses` — private; geocoding pipeline
  and São Paulo spatial analyses. Reads `analysis_ready_cohort.csv` from the
  same GDrive Data/ folder. Future work.

## Documentation files in this repo

- `README.md` — high-level project summary and findings
- `COHORT_SELECTION.md` — technical cohort definition (this doc's companion)
- `CLAUDE.md` — you are here
- `legacy/README.md` — what's in the legacy folder and why
- `ITT_Analysis/README_FINAL.md` — the ITT pipeline's own README
- `ITT_Analysis/Master_Causal_Analysis/Dicionario_de_Dados.md` — data
  dictionary (Portuguese)
