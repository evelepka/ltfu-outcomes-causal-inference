#!/usr/bin/env python3
"""
=============================================================================
 CLONE-CENSOR-WEIGHT (CCW) ANALYSIS OF LTFU AND MORTALITY  --  v2
=============================================================================

 Replaces ccw_analise_comentada.py / ccw_estimativas_pontuais.py /
 ccw_bootstrap.py, which disagreed with each other and contained three
 errors that invalidated the point estimates:

   (1) Non-adherers in the "disengage" arm were censored at treatment END
       (m_end) rather than at the end of the TARGET WINDOW.  This stretched
       the grace period across the whole treatment course, so on-treatment
       person-time and on-treatment deaths were counted as if the patient
       had followed a "disengage in month 0" strategy.  Because the error
       was worst for early windows, it manufactured the apparent
       "earlier disengagement -> higher mortality" gradient.
   (2) Administrative censoring was absent: everyone got 24 person-months
       regardless of whether 24 months of follow-up existed.  5.3% of the
       cohort has <730 days to 2024-12-31, and it is differential
       (9.5% of LTFU vs 4.8% of non-LTFU).
   (3) Deviation-vs-death ties were resolved at MONTH granularity, so a
       death in the same month as disengagement was credited to the
       "remain engaged" arm.

 v2 also changes the ESTIMAND from disjoint windows to nested strategies
 ("disengage by month T"), which is Hernan's formulation and the one
 Reviewer 3 proposed.  Disjoint windows put almost no one in the exposed
 arm for early T, so positivity fails and the weights explode; nested
 strategies keep the adherent group growing with T.

-----------------------------------------------------------------------------
 TARGET TRIAL
-----------------------------------------------------------------------------
 Time zero            : treatment initiation
 Eligibility          : the analysis cohort (n = 171,048)
 Strategies, for each T in 1..6 (months of therapy, PAPER convention):
     A_T : experience a 30+ day treatment interruption at some point
           during the first T months of therapy
     R   : remain engaged (never a 30+ day interruption; complete therapy)
 Grace period         : months 1..T.  A patient still engaged and alive
                        inside the window is compatible with BOTH strategies,
                        so deaths in that window count in BOTH arms.  This is
                        the mechanism that removes the early-window artifact.
 Deviation / censoring:
     R   : censored when the patient disengages (any month)
     A_T : censored at the END of month T if still engaged
           (adherents -- those who disengaged within the window -- are
            never artificially censored)
 Weighting            : stabilized IPCW, fit ARM-SPECIFICALLY and only on
                        person-months actually at risk of the artificial
                        censoring event
 Outcome              : all-cause death (configurable; see OUTCOME_COL)
 Estimand             : risk difference and risk ratio at 6, 12 and 24 months
                        from treatment initiation
 Inference            : person-level nonparametric bootstrap, refitting the
                        entire pipeline (clone -> censor -> weight -> estimate)

-----------------------------------------------------------------------------
 MONTH CONVENTION -- read this before comparing to anything
-----------------------------------------------------------------------------
 Internally months are 0-indexed (month 0 = days 0-30.4).
 ALL PRINTED / SAVED OUTPUT uses the PAPER convention, month = internal + 1,
 so "month 1" means the first month of therapy in both this script and the
 manuscript.  The v1 scripts printed internal indices, which is why the CCW
 report's "month 0" stratum corresponds to the manuscript's month 1.

 Usage:
   python3 ccw_v2.py                  # point estimates + diagnostics
   python3 ccw_v2.py --bootstrap 200  # add bootstrap CIs
   python3 ccw_v2.py --ltfu-date closure   # R3's sensitivity (no 30-day shift)
=============================================================================
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
ROOT = Path("/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/"
            ".shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/LTFU Paper")
COHORT_CSV = ROOT / "ITT_Analysis/data/itt_cohort.csv"
OUTDIR = ROOT / "CCW_analysis/results_v2"

MONTH_DAYS = 30.4          # days per month
HORIZON_M = 24             # analysis horizon, months (internal index 0..23)
ADMIN_CENSOR = pd.Timestamp("2024-12-31")   # end of mortality follow-up
STRATEGY_MONTHS = [1, 2, 3, 4, 5, 6]        # T, in PAPER months
REPORT_AT = [6, 12, 24]                     # report risk at these months
WEIGHT_TRUNC = (0.01, 0.99)                 # stabilized-weight truncation
OUTCOME_COL = "event_d"                     # all-cause death indicator
OUTCOME_DATE = "death_date"

# Baseline covariates for the censoring (weight) model.
# All measured at treatment initiation, so unaffected by the exposure.
# NOTE: missing values are carried as their own level ("__missing__").  The
# main manuscript uses multiple imputation; this is a deliberate difference
# that must either be reconciled or reported as such.
COVS = ["age_group", "sex", "hiv_aids", "homelessness",
        "hosp_admission", "drug_use", "clinical_clean", "dot_status"]


# ---------------------------------------------------------------------------
# STEP 1 -- PERSON-LEVEL TIMELINE, ALL CLOCKS IN DAYS FROM TREATMENT START
# ---------------------------------------------------------------------------
def load_timeline(path=COHORT_CSV, ltfu_date_mode="shift30", verbose=True):
    """Build one row per person with day-level clocks.

    ltfu_date_mode:
      'shift30' -- disengagement = end_date - 30d (manuscript primary)
      'closure' -- disengagement = end_date       (Reviewer 3's preference:
                   the date the LTFU definition is actually met)
    """
    df = pd.read_csv(path, low_memory=False)
    for c in ["best_start", "end_date", OUTCOME_DATE]:
        df[c] = pd.to_datetime(df[c], errors="coerce")

    df["is_ltfu"] = df["itt_group"].eq("Loss to follow-up").to_numpy()

    # --- disengagement clock -------------------------------------------------
    shift = 30 if ltfu_date_mode == "shift30" else 0
    dis_date = df["end_date"] - pd.Timedelta(days=shift)
    t_dis = (dis_date - df["best_start"]).dt.days.to_numpy(dtype="float64")
    t_dis[~df["is_ltfu"].to_numpy()] = np.inf          # never disengages

    # Back-shifting can push disengagement before treatment start (1,118
    # patients, min -29 d, under 'shift30').  v1 silently floored these to
    # month 0, where they became ~25% of the "primary abandonment" stratum.
    # Here they are floored explicitly and COUNTED.
    n_neg = int(np.sum(np.isfinite(t_dis) & (t_dis < 0)))
    t_dis = np.where(np.isfinite(t_dis), np.maximum(t_dis, 0.0), np.inf)

    # --- other clocks --------------------------------------------------------
    # NB: treatment-end time (end_date - best_start) is deliberately NOT carried
    # forward.  Censoring non-adherents at treatment end was bug (1) in v1; the
    # window end is the correct censoring point and is derived from T alone.
    t_death = (df[OUTCOME_DATE] - df["best_start"]).dt.days.to_numpy(dtype="float64")
    t_death = np.where(df[OUTCOME_COL].to_numpy() == 1, t_death, np.inf)
    t_death = np.where(np.isfinite(t_death), np.maximum(t_death, 0.0), np.inf)

    # THE FIX FOR BUG (2): real administrative censoring, per person.
    t_admin = (ADMIN_CENSOR - df["best_start"]).dt.days.to_numpy(dtype="float64")

    out = pd.DataFrame({
        "is_ltfu": df["is_ltfu"].to_numpy(),
        "t_dis": t_dis,
        "t_death": t_death,
        "t_admin": t_admin,
    })
    for c in COVS:
        out[c] = df[c].astype("object").where(df[c].notna(), "__missing__").astype(str)

    if verbose:
        n = len(out)
        short = int(np.sum(t_admin < HORIZON_M * MONTH_DAYS))
        print(f"  cohort                        : {n:,}")
        print(f"  LTFU                          : {out.is_ltfu.sum():,} "
              f"({out.is_ltfu.mean():.1%})")
        print(f"  deaths (any time)             : {int(np.isfinite(t_death).sum()):,}")
        print(f"  ltfu_date mode                : {ltfu_date_mode}")
        if shift:
            print(f"  disengagement < treat. start  : {n_neg:,} "
                  f"({n_neg/max(out.is_ltfu.sum(),1):.1%} of LTFU) -> floored to month 1")
        print(f"  < {HORIZON_M}mo potential follow-up  : {short:,} ({short/n:.1%})"
              f"  [LTFU {np.mean(t_admin[out.is_ltfu.to_numpy()] < HORIZON_M*MONTH_DAYS):.1%}"
              f" vs non-LTFU {np.mean(t_admin[~out.is_ltfu.to_numpy()] < HORIZON_M*MONTH_DAYS):.1%}]")
    return out


# ---------------------------------------------------------------------------
# STEP 2 -- CLONE AND CENSOR
# ---------------------------------------------------------------------------
def build_arm(tl, arm, T=None):
    """One arm of clones.

    arm='remain'    : deviates when the patient disengages (any time)
    arm='disengage' : adherents = disengaged within paper months 1..T (never
                      censored); everyone else deviates at the END of month T

    Month arithmetic is done in INTEGER months throughout.  Do not reintroduce
    a months->days conversion here: MONTH_DAYS = 30.4 is not binary-exact, so
    floor(24 * 30.4 / 30.4) == 23, which previously truncated follow-up to 23
    months and shifted the disengage-arm boundary a month early.

    Internal months are 0-indexed; paper month = internal + 1.  Strategy T
    therefore admits internal disengagement months 0..T-1.

    The deviation-vs-death tie (bug 3) is resolved at DAY level in the remain
    arm, where deviation happens at an individual instant t_dis.  In the
    disengage arm the boundary falls exactly on a month edge, so the integer
    comparison m_death < T is already unambiguous.
    """
    t_dis = tl["t_dis"].to_numpy()
    t_death = tl["t_death"].to_numpy()
    t_admin = tl["t_admin"].to_numpy()

    # days -> integer months, exactly once
    with np.errstate(invalid="ignore"):
        m_dis = np.where(np.isfinite(t_dis), np.floor(t_dis / MONTH_DAYS), np.inf)
        m_death = np.where(np.isfinite(t_death), np.floor(t_death / MONTH_DAYS), np.inf)
        m_admin = np.floor(t_admin / MONTH_DAYS)

    if arm == "remain":
        m_dev = m_dis                                    # inf if never disengages
        death_first = t_death < t_dis                    # day-level tie break
        at_risk_mode = "all"
    elif arm == "disengage":
        adherent = m_dis <= (T - 1)                      # exact integer test
        # THE FIX FOR BUG (1): non-adherents are censored at the end of the
        # TARGET WINDOW, not at treatment completion.
        m_dev = np.where(adherent, np.inf, float(T))
        death_first = m_death < m_dev                    # month edge, no ties
        at_risk_mode = "window_end"
    else:
        raise ValueError(arm)

    # integer month at which the clone is censored, whichever comes first
    m_cens = np.minimum(m_dev, np.minimum(m_admin, float(HORIZON_M)))

    has_event = death_first & (t_death < t_admin) & (m_death < HORIZON_M)

    # months contributed: through the death month if an event, else up to (not
    # including) the month in which censoring occurs
    nmonths = np.where(has_event, m_death + 1, m_cens)
    nmonths = np.clip(nmonths, 0, HORIZON_M).astype(np.int32)

    # was the stop caused by DEVIATION (weightable) rather than admin/horizon?
    dev_binds = (~has_event) & np.isfinite(m_dev) \
                & (m_dev <= m_admin) & (m_dev < HORIZON_M)

    out = pd.DataFrame({
        "nmonths": nmonths,
        "event_month": np.where(has_event, m_death, np.inf),
        "dev_binds": dev_binds,
        "dev_month": m_dev,
    })
    for c in COVS:
        out[c] = tl[c].to_numpy()
    out["at_risk_mode"] = at_risk_mode
    out["T"] = -1 if T is None else T
    return out[out["nmonths"] > 0].reset_index(drop=True)


# ---------------------------------------------------------------------------
# STEP 3 -- EXPAND TO PERSON-MONTHS
# ---------------------------------------------------------------------------
def expand(clones):
    n = clones["nmonths"].to_numpy()
    idx = np.repeat(np.arange(len(clones)), n)
    L = clones.iloc[idx].reset_index(drop=True)
    L["month"] = np.concatenate([np.arange(k) for k in n]).astype(np.int32)
    L["clone"] = idx.astype(np.int64)

    L["death"] = ((L["event_month"].to_numpy() < np.inf)
                  & (L["month"].to_numpy() == L["event_month"].to_numpy())).astype(np.int8)

    # artificial-censoring indicator: fires in the LAST contributed month of a
    # clone whose stop was caused by deviation
    last = L["month"].to_numpy() == (L["nmonths"].to_numpy() - 1)
    L["dev_next"] = (last & L["dev_binds"].to_numpy()
                     & (L["death"].to_numpy() == 0)).astype(np.int8)

    # rows genuinely AT RISK of the artificial censoring event.  Restricting
    # the weight model to these avoids the separation that arises when months
    # in which nobody can deviate are included as zeros.
    if L["at_risk_mode"].iloc[0] == "all":
        L["at_risk"] = np.int8(1)
    else:                                    # 'window_end': only the last month of the window
        T = int(L["T"].iloc[0])
        L["at_risk"] = (L["month"].to_numpy() == (T - 1)).astype(np.int8)
    return L


# ---------------------------------------------------------------------------
# STEP 4 -- STABILIZED IPCW, ARM-SPECIFIC
# ---------------------------------------------------------------------------
def add_ipcw(L, verbose=False, label=""):
    """sw = prod_k (1 - p_num_k) / (1 - p_den_k) over months at risk."""
    n_rows = len(L)
    p_den = np.zeros(n_rows)
    p_num = np.zeros(n_rows)

    fit_rows = (L["at_risk"].to_numpy() == 1)
    y = L["dev_next"].to_numpy()[fit_rows]

    if fit_rows.sum() > 0 and 0 < y.sum() < len(y):
        sub = L.loc[fit_rows, COVS + ["month"]].copy()
        sub["month"] = sub["month"].astype(str)
        enc_d = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
        Xd = enc_d.fit_transform(sub[COVS + ["month"]])
        enc_n = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
        Xn = enc_n.fit_transform(sub[["month"]])
        # mild regularization: C=1e6 in v1 invited separation on sparse cells
        p_den[fit_rows] = LogisticRegression(max_iter=1000, C=1.0
                                             ).fit(Xd, y).predict_proba(Xd)[:, 1]
        p_num[fit_rows] = LogisticRegression(max_iter=1000, C=1.0
                                             ).fit(Xn, y).predict_proba(Xn)[:, 1]

    eps = 1e-8
    ratio = (1.0 - p_num) / np.clip(1.0 - p_den, eps, None)
    cum = pd.Series(ratio, index=L.index).groupby(L["clone"]).cumprod()
    # Censoring at the END of month k removes the clone from month k+1 onward,
    # so month k must be weighted by the product through k-1 only.  Without
    # this lag the row on which a clone is censored discounts its own censoring
    # and stabilized weights no longer average 1.
    L["sw"] = cum.groupby(L["clone"]).shift(1).fillna(1.0).to_numpy()

    lo, hi = np.quantile(L["sw"].to_numpy(), WEIGHT_TRUNC)
    L["swt"] = np.clip(L["sw"].to_numpy(), lo, hi)

    if verbose:
        sw = L["sw"].to_numpy()
        pct = float(np.mean((sw < lo) | (sw > hi)))
        flag = "" if 0.9 <= sw.mean() <= 1.1 else "   <-- CHECK: far from 1.0"
        print(f"      weights[{label:<14s}] mean={sw.mean():.3f} sd={sw.std():.3f} "
              f"min={sw.min():.3f} max={sw.max():.2f} trunc={pct:.2%}{flag}")
        if fit_rows.sum() > 0 and not (0 < y.sum() < len(y)):
            print(f"      WARNING[{label}]: censoring model not identified "
                  f"(at-risk rows={fit_rows.sum():,}, events={int(y.sum())}); "
                  f"weights set to 1")
    return L


# ---------------------------------------------------------------------------
# STEP 5 -- WEIGHTED CUMULATIVE INCIDENCE
# ---------------------------------------------------------------------------
def weighted_risk(L, label="", strict=True):
    """Discrete-time weighted cumulative incidence, indexed by internal month."""
    w = L["swt"].to_numpy()
    d = L["death"].to_numpy()
    m = L["month"].to_numpy()
    num = np.bincount(m, weights=w * d, minlength=HORIZON_M)
    den = np.bincount(m, weights=w, minlength=HORIZON_M)
    n_rows = np.bincount(m, minlength=HORIZON_M)

    # Self-checks.  A hazard of exactly 1 means some month is populated only by
    # deaths, which is how the 30.4-day rounding bug surfaced (month 23 held 47
    # rows, all deaths -> cumulative incidence 100%).
    if strict:
        thin = np.where((n_rows > 0) & (n_rows < 100))[0]
        if len(thin):
            print(f"      WARNING[{label}]: months with <100 person-months: "
                  f"{thin.tolist()} (n={n_rows[thin].tolist()})")
        if np.any(m >= HORIZON_M):
            raise AssertionError(f"{label}: month index >= HORIZON_M")

    hz = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    if strict and np.any(hz >= 1.0):
        raise AssertionError(
            f"{label}: weighted hazard >= 1 at month(s) "
            f"{np.where(hz >= 1.0)[0].tolist()} -- a month containing only "
            f"deaths. Check month arithmetic / administrative censoring.")
    return 1.0 - np.cumprod(1.0 - hz)          # length HORIZON_M, index = month


def risk_at(ci, month_paper):
    """Risk at a paper-convention month (1-indexed)."""
    return float(ci[min(month_paper, HORIZON_M) - 1])


# ---------------------------------------------------------------------------
# STEP 6 -- ONE CONTRAST
# ---------------------------------------------------------------------------
def reference_arm(tl, verbose=False):
    """'Remain engaged' does not depend on T, so build it once."""
    L = expand(build_arm(tl, "remain"))
    L = add_ipcw(L, verbose=verbose, label="remain")
    return weighted_risk(L, "remain", strict=verbose), L


def contrast(tl, T, ci_ref, verbose=False):
    L = expand(build_arm(tl, "disengage", T=T))
    L = add_ipcw(L, verbose=verbose, label=f"disengage T={T}")
    ci_dis = weighted_risk(L, f"disengage T={T}", strict=verbose)
    res = {"T": T}
    for mo in REPORT_AT:
        r1, r0 = risk_at(ci_dis, mo), risk_at(ci_ref, mo)
        res[f"risk{mo}_dis"] = 100 * r1
        res[f"risk{mo}_rem"] = 100 * r0
        res[f"rr{mo}"] = r1 / r0 if r0 > 0 else np.nan
        res[f"rd{mo}"] = 100 * (r1 - r0)
    m_dis = np.floor(np.where(np.isfinite(tl["t_dis"].to_numpy()),
                              tl["t_dis"].to_numpy(), np.inf) / MONTH_DAYS)
    res["n_adherent"] = int(np.sum(m_dis <= (T - 1)))
    res["deaths_dis_arm"] = int(L["death"].sum())
    res["mean_weight"] = float(L["sw"].mean())
    # Positivity flag: a stabilized weight mean far from 1 means the adherent
    # group cannot be reweighted back to the full cohort on baseline covariates
    # alone -- the estimate is not trustworthy however tight its CI looks.
    res["estimable"] = bool(0.9 <= res["mean_weight"] <= 1.1)
    return res


# ---------------------------------------------------------------------------
# BOOTSTRAP
# ---------------------------------------------------------------------------
def bootstrap(tl, B, seed=2026):
    """Person-level resample; refit clone -> censor -> weight -> estimate.

    Uses the SAME covariate set as the point estimates (v1's bootstrap used a
    reduced set 'for speed', so its CIs did not correspond to its estimates).
    """
    keys = [f"{m}{mo}" for mo in REPORT_AT for m in ("rr", "rd")]
    draws = {T: {k: [] for k in keys} for T in STRATEGY_MONTHS}
    rng = np.random.default_rng(seed)
    n = len(tl)
    t0 = time.time()
    for b in range(B):
        d = tl.iloc[rng.integers(0, n, n)].reset_index(drop=True)
        ci_ref, _ = reference_arm(d)
        for T in STRATEGY_MONTHS:
            r = contrast(d, T, ci_ref)
            for k in keys:
                draws[T][k].append(r[k])
        if (b + 1) % 10 == 0:
            el = time.time() - t0
            print(f"    rep {b+1}/{B}  {el/(b+1):.1f}s/rep  "
                  f"ETA {el/(b+1)*(B-b-1)/60:.1f} min", flush=True)
    return draws


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=0,
                    help="number of bootstrap replicates (0 = point estimates only)")
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"], default="shift30")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("CCW v2 -- nested 'disengage by month T' strategies vs remain engaged")
    print("=" * 78)
    tl = load_timeline(ltfu_date_mode=args.ltfu_date)

    print("\n--- weight diagnostics ---")
    ci_ref, Lref = reference_arm(tl, verbose=True)
    rows = []
    for T in STRATEGY_MONTHS:
        rows.append(contrast(tl, T, ci_ref, verbose=True))

    print(f"\n--- reference arm: remain engaged ---")
    for mo in REPORT_AT:
        print(f"    {mo:>2}-month risk : {100*risk_at(ci_ref, mo):5.2f}%")
    print(f"    deaths contributing: {int(Lref['death'].sum()):,}")

    res = pd.DataFrame(rows)
    print("\n--- effect of disengaging by month T (paper months) ---")
    hdr = (f"{'T':>2} {'adherent':>9} {'6mo RR':>8} {'6mo RD':>8} "
           f"{'24mo RR':>8} {'24mo RD':>8} {'risk24 dis':>11} {'ok':>4}")
    print(hdr)
    print("-" * len(hdr))
    for _, r in res.iterrows():
        print(f"{int(r['T']):>2} {int(r['n_adherent']):>9,} "
              f"{r['rr6']:>8.2f} {r['rd6']:>+8.2f} "
              f"{r['rr24']:>8.2f} {r['rd24']:>+8.2f} "
              f"{r['risk24_dis']:>10.2f}% {'y' if r['estimable'] else 'NO':>4}")
    bad = res.loc[~res["estimable"], "T"].tolist()
    if bad:
        print(f"\n  !! T={bad} failed the positivity check (stabilized weights far "
              f"from 1).\n     Too few patients follow 'disengage this early' for the "
              f"adherent group to be\n     reweighted to the cohort on baseline "
              f"covariates. Do not report these as\n     estimates; report them as "
              f"non-estimable. Note that T=1 also absorbs the\n     "
              f"{'1,118' if args.ltfu_date == 'shift30' else '0'} patients whose "
              f"back-shifted disengagement date precedes treatment start.")

    print("\n  Interpreting the T gradient: these are NESTED strategies, so the\n"
          "  exposed group GROWS with T ('disengage at any point in months 1..T').\n"
          "  RR(T) rising with T therefore does NOT mean later disengagement is\n"
          "  more harmful -- it cannot be read as a timing dose-response. For the\n"
          "  'no safe point to disengage' claim, contrast consecutive T or use\n"
          "  disjoint windows and accept the wider intervals.")

    out = {"config": {"ltfu_date_mode": args.ltfu_date, "horizon_m": HORIZON_M,
                      "month_days": MONTH_DAYS, "covariates": COVS,
                      "weight_truncation": WEIGHT_TRUNC,
                      "admin_censor": str(ADMIN_CENSOR.date()),
                      "month_convention": "paper (1-indexed)"},
           "reference_risk": {str(mo): 100 * risk_at(ci_ref, mo) for mo in REPORT_AT},
           "point": res.to_dict(orient="records")}

    if args.bootstrap:
        print(f"\n--- bootstrap ({args.bootstrap} reps, full covariate set) ---")
        draws = bootstrap(tl, args.bootstrap, seed=args.seed)
        print("\n--- effect of disengaging by month T, with 95% CI ---")
        for T in STRATEGY_MONTHS:
            pt = res.loc[res["T"] == T].iloc[0]
            ci = {}
            for k in draws[T]:
                a = np.asarray(draws[T][k], dtype=float)
                a = a[np.isfinite(a)]
                ci[k] = (float(np.percentile(a, 2.5)),
                         float(np.percentile(a, 97.5))) if len(a) else (np.nan, np.nan)
            out.setdefault("ci", {})[str(T)] = ci
            print(f"  T={T}: 24mo RR {pt['rr24']:.2f} "
                  f"({ci['rr24'][0]:.2f}-{ci['rr24'][1]:.2f})   "
                  f"24mo RD {pt['rd24']:+.2f}pp "
                  f"({ci['rd24'][0]:+.2f} to {ci['rd24'][1]:+.2f})   "
                  f"6mo RR {pt['rr6']:.2f} "
                  f"({ci['rr6'][0]:.2f}-{ci['rr6'][1]:.2f})")

    tag = args.ltfu_date
    (OUTDIR / f"ccw_v2_{tag}.json").write_text(json.dumps(out, indent=2))
    res.to_csv(OUTDIR / f"ccw_v2_{tag}_point.csv", index=False)
    print(f"\nSaved -> {OUTDIR}/ccw_v2_{tag}.json  (+ _point.csv)")

    print("\n" + "=" * 78)
    print("NOT YET IMPLEMENTED -- required before this can be the primary analysis:")
    print("  * cause-specific outcomes (TB vs non-TB death).  The cohort file has")
    print("    no cause-of-death column; merge SIM ICD-10 / TBweb closure cause and")
    print("    set OUTCOME_COL / OUTCOME_DATE.  This is the manuscript's main")
    print("    defense against residual confounding and cloning does NOT fix it.")
    print("  * subgroup contrasts (age, sex, HIV, homelessness, resistance, period)")
    print("  * reconciliation with the manuscript's multiple imputation: weights")
    print("    here carry missingness as its own level ('__missing__')")
    print("=" * 78)


if __name__ == "__main__":
    main()
