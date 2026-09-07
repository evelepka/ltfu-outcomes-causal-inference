#!/usr/bin/env python3
"""Evidence bearing on the missing-at-random assumption (appendix section 4.2).

WHY THIS EXISTS
---------------
Reviewer 1, comment 1.13: "the authors assume MAR but provide no justification."
MAR is not testable -- the observed data are equally compatible with MAR and
MNAR, because the distinguishing information is exactly what is missing. So this
script does NOT test MAR. It produces the three things that can be shown:

  1. HOW MUCH is missing, and whether missingness differs by exposure arm.
     Establishes that the data are not MCAR and that the imputation model must
     condition on the exposure (it does).

  2. WHAT PREDICTS missingness among the fully observed covariates. MAR requires
     the imputation model to contain the drivers of missingness; this shows which
     they are, so the claim can be checked rather than asserted.

  3. WHETHER the missingness indicator still predicts mortality once the observed
     covariates are conditioned on. This comparison exists only because the
     OUTCOME is complete here and the missingness is confined to covariates; it
     would be unavailable if the outcome were the thing missing.

Point 3 is the closest thing to evidence, and it is still not a test. A residual
association does not falsify MAR: the imputation model includes the outcome, the
exposure and follow-up time, so MAR is assumed conditional on those, and MAR does
not require the missingness indicator to be independent of the outcome. Nor does a
null result establish MAR. Both readings are stated in the appendix text.

The standardisation follows the paper's own idiom: fit a model, predict under each
value of the indicator for every individual, and average over the observed
covariate distribution.

Outputs
  ITT_Analysis/results/missingness_mechanism.csv   one row per covariate
  ITT_Analysis/results/missingness_by_year.csv     recording-practice trend

Usage
  python3 ITT_Analysis/scripts/49c_missingness_mechanism.py
  B=0 python3 ...   skip the bootstrap (point estimates only)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[2]
COHORT = ROOT / "ITT_Analysis" / "data" / "itt_cohort.csv"
OUT_MAIN = ROOT / "ITT_Analysis" / "results" / "missingness_mechanism.csv"
OUT_YEAR = ROOT / "ITT_Analysis" / "results" / "missingness_by_year.csv"

B = int(os.environ.get("B", "500"))
SEED = 20260831

# HORIZON and the follow-up restriction, 2026-08-31.
#
# The first version of this script used a five-year horizon over the whole cohort
# and scored anyone without a recorded death as a five-year survivor. That
# misclassifies administrative censoring as survival, and it does so
# DIFFERENTIALLY here: race is unrecorded for 31.1% of 2023 treatment starts
# against 10.6% of 2014, so the unrecorded group is enriched with people who
# cannot yet have an observed five-year death. The race contrast moved from
# -0.32 to -0.73 pp once that was corrected, which is small but is the difference
# between a number a reviewer can check and one they can dismiss.
#
# Fixed by restricting to individuals with a FULL horizon of potential follow-up.
# Two years is the default rather than five: it admits treatment starts through
# 2022 and so retains the pandemic-era rise in missingness, where a five-year
# window would discard exactly the years that carry the recording-practice
# signal. The horizon is arbitrary for a missingness diagnostic -- it is not an
# effect estimate -- so coverage is worth more than matching the paper's headline.
# Set HORIZON=5 to reproduce the restricted five-year version.
HORIZON = float(os.environ.get("HORIZON", "2"))

# The 14 model covariates, minus geo4 which is complete by construction (an
# unmatched municipality is assigned the urban reference), matching the
# complete-case rule in script 43b.
COVARS = ["age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
          "alcohol", "drug_use", "incarcerated", "homelessness",
          "hosp_admission", "clinical_clean", "dot_status"]

d = pd.read_csv(COHORT, low_memory=False)
d["start"] = pd.to_datetime(d["best_start"], errors="coerce")
d["year"] = d["start"].dt.year
d["ltfu"] = (d["itt_group"] == "Loss to follow-up").astype(int)

# TWO POPULATIONS, on purpose.
#
# The censoring problem affects only the comparison that uses the OUTCOME. How
# much is missing, whether it differs by arm, and what predicts it involve no
# outcome at all, so they are computed on the FULL cohort. Restricting those to a
# follow-up window would discard 2023 -- the year race is unrecorded for 31.1% of
# starts against 10.6% in 2014 -- and so throw away the strongest evidence for the
# recording-practice mechanism for no methodological gain.
#
# The standardised mortality difference IS outcome-based, so it is computed on
# individuals with a full horizon of potential follow-up, where censoring can
# never be scored as survival.
CUT = (d["start"] + pd.to_timedelta(d["time_d"] * 365.25, unit="D")).max()
FULL_FU = d["start"] <= (CUT - pd.Timedelta(days=int(HORIZON * 365.25)))
d["died"] = ((d["event_d"] == 1) & (d["time_d"] <= HORIZON)).astype(int)
print(f"administrative cut-off {CUT.date()}")
print(f"  descriptives          : full cohort, {len(d):,}, "
      f"starts {d.year.min()}-{d.year.max()}")
print(f"  outcome-based contrast: {FULL_FU.sum():,} with a full {HORIZON:g} y of "
      f"potential follow-up, starts {d.loc[FULL_FU,'year'].min()}-"
      f"{d.loc[FULL_FU,'year'].max()}")

# blanks count as missing: the CSV encodes missing categoricals both ways
for c in COVARS:
    if d[c].dtype == object:
        d[c] = d[c].replace(r"^\s*$", np.nan, regex=True)

# conditioning set: covariates with NO missingness of their own, so they can
# condition without themselves needing imputation
OBS = [c for c in COVARS if d[c].isna().sum() == 0]
EXTRA = [c for c in ("tobacco_use", "mental_health", "other_immuno_condition",
                     "resistance_clean", "lab_confirmed_stat")
         if c in d.columns and d[c].isna().sum() == 0]
COND = OBS + EXTRA

print(f"cohort {len(d):,}  |  horizon {HORIZON} y  |  bootstrap B={B}")
print(f"conditioning on {len(COND)} fully observed covariates:\n  {', '.join(COND)}\n")

# design matrix and outcome for the OUTCOME-BASED contrast: full-follow-up only
dr = d.loc[FULL_FU].reset_index(drop=True)
XB = pd.get_dummies(dr[COND].astype(str), drop_first=False).values.astype(float)
Y = dr["died"].values
rng = np.random.default_rng(SEED)


def standardised_diff(flag, idx=None):
    """Standardised difference in two-year mortality, flag=1 vs flag=0."""
    if idx is None:
        idx = np.arange(len(flag))
    X = np.column_stack([XB[idx], flag[idx]])
    m = LogisticRegression(max_iter=2000, C=1e6, solver="lbfgs")
    m.fit(X, Y[idx])
    n = X.shape[0]
    X1 = np.column_stack([XB[idx], np.ones(n)])
    X0 = np.column_stack([XB[idx], np.zeros(n)])
    return (m.predict_proba(X1)[:, 1].mean() - m.predict_proba(X0)[:, 1].mean()) * 100


rows = []
for v in COVARS:
    miss = d[v].isna().values
    n_miss = int(miss.sum())
    rec = {"variable": v, "n_missing": n_miss,
           "pct_missing": round(100 * n_miss / len(d), 2)}
    if n_miss == 0:
        rows.append(rec)
        print(f"  {v:<16} complete")
        continue

    rec["pct_missing_ltfu"] = round(100 * miss[d.ltfu == 1].mean(), 2)
    rec["pct_missing_incare"] = round(100 * miss[d.ltfu == 0].mean(), 2)
    rec["diff_ltfu_incare_pp"] = round(rec["pct_missing_ltfu"]
                                       - rec["pct_missing_incare"], 2)
    mr = dr[v].isna().values
    rec["crude_mort_diff_pp"] = round(
        100 * (Y[mr].mean() - Y[~mr].mean()), 2)

    flag = dr[v].isna().values.astype(float)      # restricted population
    point = standardised_diff(flag)
    rec["std_mort_diff_pp"] = round(point, 2)
    if B > 0 and int(mr.sum()) >= 50:
        draws, failures = [], []
        for _ in range(B):
            # over dr, NOT d: XB and Y are the restricted population. Drawing
            # indices over the full cohort put them out of range, every
            # replicate raised IndexError, and the bare except below swallowed
            # it -- so the CI columns silently vanished from the output.
            idx = rng.integers(0, len(dr), len(dr))
            try:
                draws.append(standardised_diff(flag, idx))
            except Exception as exc:                            # noqa: BLE001
                failures.append(repr(exc))
                continue
        if failures:
            # a handful of singular fits is tolerable; wholesale failure is a bug
            print(f"    {len(failures)}/{B} replicates failed, e.g. {failures[0][:90]}")
            if len(failures) > B // 2:
                raise RuntimeError(
                    f"{v}: {len(failures)} of {B} bootstrap replicates failed -- "
                    f"refusing to emit a CI-less row silently")
        if draws:
            lo, hi = np.percentile(draws, [2.5, 97.5])
            rec["std_ci_low"] = round(lo, 2)
            rec["std_ci_high"] = round(hi, 2)

    # strongest observed predictor of this variable being missing;
    # meaningless when only a handful of values are missing
    best = None
    if n_miss < 50:
        rows.append(rec)
        print(f"  {v:<16} {n_miss:>7,} -- too few for a diagnostic")
        continue
    for c in COND + ["year"]:
        g = pd.DataFrame({"m": miss, "lv": d[c].astype(str)}).groupby("lv")["m"]
        agg = g.agg(["mean", "size"])
        agg = agg[agg["size"] >= 200]
        if len(agg) < 2:
            continue
        spread = (agg["mean"].max() - agg["mean"].min()) * 100
        if best is None or spread > best[1]:
            best = (c, spread, agg["mean"].idxmin(), agg["mean"].min() * 100,
                    agg["mean"].idxmax(), agg["mean"].max() * 100)
    if best:
        rec.update(top_predictor=best[0],
                   top_predictor_spread_pp=round(best[1], 1),
                   top_lo_level=str(best[2]), top_lo_pct=round(best[3], 1),
                   top_hi_level=str(best[4]), top_hi_pct=round(best[5], 1))
    rows.append(rec)
    ci = (f" ({rec.get('std_ci_low')} to {rec.get('std_ci_high')})"
          if "std_ci_low" in rec else "")
    print(f"  {v:<16} {n_miss:>7,} ({rec['pct_missing']:5.2f}%)  "
          f"LTFU {rec['pct_missing_ltfu']:5.2f}% vs {rec['pct_missing_incare']:5.2f}%  "
          f"crude {rec['crude_mort_diff_pp']:+.2f} pp  "
          f"standardised {rec['std_mort_diff_pp']:+.2f} pp{ci}")
    if best:
        print(f"  {'':<16} strongest predictor: {best[0]} "
              f"({best[2]} {best[3]:.1f}% -> {best[4]} {best[5]:.1f}%)")

out = pd.DataFrame(rows)
out.to_csv(OUT_MAIN, index=False)
print(f"\nwrote {OUT_MAIN.relative_to(ROOT)} ({len(out)} rows)")

# ---------------------------------------------------------------------------
# Every predictor level, not just the strongest
# ---------------------------------------------------------------------------
# The appendix quotes several of these (homelessness, mental health,
# hospitalisation for education; year for race), so each needs a row in a source
# file rather than existing only in this script's stdout.
pred_rows = []
for v in [c for c in COVARS if d[c].isna().sum() > 50]:
    miss = d[v].isna().values
    for c in COND + ["year"]:
        agg = (pd.DataFrame({"m": miss, "lv": d[c].astype(str)})
               .groupby("lv")["m"].agg(["mean", "size"]))
        for lv, r in agg.iterrows():
            if r["size"] < 200:
                continue
            pred_rows.append({"variable": v, "predictor": c, "level": lv,
                              "n": int(r["size"]),
                              "pct_missing": round(100 * r["mean"], 2)})
pred = pd.DataFrame(pred_rows)
OUT_PRED = ROOT / "ITT_Analysis" / "results" / "missingness_predictors.csv"
pred.to_csv(OUT_PRED, index=False)
print(f"wrote {OUT_PRED.relative_to(ROOT)} ({len(pred)} rows)")
print("\n  levels quoted in appendix 4.2:")
for v, c, lv in [("edu_clean", "homelessness", "Yes"),
                 ("edu_clean", "homelessness", "No"),
                 ("edu_clean", "mental_health", "Yes"),
                 ("edu_clean", "hosp_admission", "Yes"),
                 ("race_clean", "year", "2014"),
                 ("race_clean", "year", "2023")]:
    sel = pred[(pred.variable == v) & (pred.predictor == c) & (pred.level == lv)]
    if len(sel):
        print(f"    {v:<12} {c:<16} {lv:<6} {sel.iloc[0].pct_missing:5.2f}%")

incomplete = [c for c in COVARS if d[c].isna().sum() > 0]
yr = (d.groupby("year")[incomplete].apply(lambda g: (g.isna().mean() * 100).round(2))
      .join(d.groupby("year").size().rename("n")))
yr.to_csv(OUT_YEAR)
print(f"wrote {OUT_YEAR.relative_to(ROOT)}")
print(f"\n{yr.to_string()}")

# ---------------------------------------------------------------------------
# DOT status: does "unrecorded" behave like "not given" rather than "given"?
# ---------------------------------------------------------------------------
# This one matters for more than the assumption. Imputation under MAR assigns the
# unrecorded to the recorded conditional distribution, which is dominated by
# "Yes". If unrecorded in fact usually means DOT was not given, that misclassifies
# them on an adjustment covariate.
print("\n--- DOT status: profile of unrecorded against recorded ---")
prof = d.assign(dot3=d["dot_status"].fillna("(unrecorded)")).groupby("dot3").agg(
    n=("died", "size"),
    mort5_pct=("died", lambda s: round(100 * s.mean(), 2)),
    ltfu_pct=("ltfu", lambda s: round(100 * s.mean(), 2)),
)
recorded = int(prof.loc[prof.index != "(unrecorded)", "n"].sum())
yes = int(prof.loc["Yes", "n"]) if "Yes" in prof.index else 0
prof["share_of_recorded_pct"] = [
    round(100 * n / recorded, 1) if i != "(unrecorded)" else np.nan
    for i, n in zip(prof.index, prof["n"])]
print(prof.to_string())
print(f"  'Yes' is {round(100 * yes / recorded, 1)}% of recorded values, so a "
      f"missing-at-random imputation sends most unrecorded individuals there.")
prof.to_csv(ROOT / "ITT_Analysis" / "results" / "missingness_dot_profile.csv")
print(f"  wrote ITT_Analysis/results/missingness_dot_profile.csv")

print("\n--- for the appendix text ---")
for v in incomplete:
    r = out.loc[out.variable == v].iloc[0]
    ci = (f" (95% CI {r.std_ci_low} to {r.std_ci_high})"
          if "std_ci_low" in out.columns and pd.notna(r.get("std_ci_low")) else "")
    print(f"  {v}: {int(r.n_missing):,} ({r.pct_missing}%); "
          f"LTFU {r.pct_missing_ltfu}% vs in care {r.pct_missing_incare}%; "
          f"crude {r.crude_mort_diff_pp} pp -> standardised "
          f"{r.std_mort_diff_pp} pp{ci}")
