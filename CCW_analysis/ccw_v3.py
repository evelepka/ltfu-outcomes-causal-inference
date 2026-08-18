#!/usr/bin/env python3
"""
=============================================================================
 CLONE-CENSOR-WEIGHT (CCW) ANALYSIS OF LTFU AND MORTALITY  --  v3
=============================================================================
 Supersedes ccw_v2.py (and the three v1 scripts, which can be deleted).

 v3 adds the three things v2 flagged as missing before CCW could be primary:

   1. CAUSE-SPECIFIC OUTCOMES (TB-attributable vs non-TB death), replicating
      the attribution scheme in ITT_Analysis/scripts/30e_itt_target_trial_
      cause_specific.R so the CCW figure is comparable to main-text Figure 5.
   2. SUBGROUPS: age, sex, HIV, homelessness, drug resistance, calendar period
      (main-text Figure 4D).
   3. MULTIPLE IMPUTATION: runs across ITT_Analysis/data/mi/imp_*.csv and
      pools, replacing v2's missing-as-its-own-level shortcut.

 Carried over from v2 (do not regress these):
   * Non-adherents are censored at the END OF THE TARGET WINDOW, not at
     treatment completion.  Censoring at treatment end was the v1 error that
     stretched the grace period across the whole course and manufactured the
     "earlier disengagement is worse" gradient.
   * Real administrative censoring per person (5.3% of the cohort has <24
     months of possible follow-up, 9.5% of LTFU vs 4.8% of non-LTFU).
   * All month arithmetic in INTEGER months.  MONTH_DAYS = 30.4 is not
     binary-exact: floor(24 * 30.4 / 30.4) == 23.  Never convert months->days.
   * IPCW lagged one month: censoring at the end of month k discounts months
     k+1 onward, not month k itself.
   * Self-checks that raise on a weighted hazard >= 1 and flag stabilized
     weights whose mean is far from 1 (a positivity failure).

-----------------------------------------------------------------------------
 TARGET TRIAL
-----------------------------------------------------------------------------
 Time zero   : treatment initiation
 Strategies, for each T in 1..6 (PAPER months; internal index = paper - 1):
     A_T : experience a 30+ day treatment interruption at some point during
           the first T months of therapy
     R   : remain engaged (never a 30+ day interruption)
 Grace       : months 1..T.  A patient still engaged and alive inside the
               window is compatible with BOTH strategies, so deaths there are
               counted in BOTH arms.  This is what removes the early-window
               artifact that the landmark design could not.
 Censoring   : R   -> at disengagement (any month)
               A_T -> at the end of month T if still engaged
                      (adherents are never artificially censored)
 Weighting   : stabilized IPCW, arm-specific, fit only on person-months
               genuinely at risk of the artificial censoring event
 Outcomes    : all-cause death; TB-attributable death; non-TB death
               (cause-specific: competing-cause deaths are censored)
 Estimand    : risk difference / risk ratio at 6, 12, 24 months from
               treatment initiation
 Inference   : person-level bootstrap refitting the whole pipeline, with
               imputation uncertainty propagated (see --rubin)

 Usage
   python3 ccw_v3.py                              # point estimates, all MI
   python3 ccw_v3.py --bootstrap 500              # + CIs (cycles imputations)
   python3 ccw_v3.py --bootstrap 200 --rubin      # + CIs (full M x B, Rubin)
   python3 ccw_v3.py --subgroups                  # subgroup contrasts
   python3 ccw_v3.py --ltfu-date closure          # Reviewer 3 timing sensitivity
   python3 ccw_v3.py --cause tb_simonly           # SIM-only attribution
=============================================================================
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack, vstack
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder

# ---------------------------------------------------------------------------
# PATHS / CONFIG
# ---------------------------------------------------------------------------
ROOT = Path("/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/"
            ".shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/LTFU Paper")
MI_DIR = ROOT / "ITT_Analysis/data/mi"
RAW_CSV = ROOT / "Data/Final_table_cleaned.csv"
OUTDIR = ROOT / "CCW_analysis/results_v3"
CACHE = OUTDIR / "cause_lookup_fixedattr.csv"   # ADR-0003

MONTH_DAYS = 30.4
HORIZON_M = 24
ADMIN_CENSOR = pd.Timestamp("2024-12-31")
STRATEGY_MONTHS = [1, 2, 3, 4, 5, 6]
PRIMARY_T = 6                       # strategy used for subgroup contrasts
REPORT_AT = [6, 12, 24]
WEIGHT_TRUNC = (0.01, 0.99)
DISJOINT = False   # set by --disjoint

# `geo4` is a covariate in BOTH designs as of 2026-08-18. It was added to the
# rolling landmark first; keeping CCW's set identical matters because CCW is the
# secondary analysis supporting the response to Reviewer 3, and a reviewer
# comparing the two would otherwise find them adjusted differently.
COVS = ["age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
        "alcohol", "drug_use", "incarcerated", "homelessness",
        "hosp_admission", "clinical_clean", "dot_status", "geo4"]
GEO_REF = "Urbano"      # reference level; custody folds in here

SUBGROUPS = ["age_group", "sex", "hiv_aids", "homelessness",
             "resistance_clean", "period"]

CAUSES = {
    # name          -> (event mask column, human label)
    "all_cause":    ("__any_death__", "all-cause death"),
    "tb_hybrid":    ("tb_hybrid",     "TB-attributable death (hybrid: SIM + TBweb)"),
    "nontb_hybrid": ("nontb_hybrid",  "non-TB death (hybrid: SIM + TBweb)"),
    "tb_simonly":   ("tb_simonly",    "TB-attributable death (SIM ICD-10 only)"),
    "nontb_simonly": ("nontb_simonly", "non-TB death (SIM ICD-10 only)"),
}


# ---------------------------------------------------------------------------
# CAUSE-OF-DEATH ATTRIBUTION
# ---------------------------------------------------------------------------
def build_cause_lookup(force=False, verbose=True):
    """Replicates classify_cod() from 30e_itt_target_trial_cause_specific.R.

    TB-strict   : SIM ICD-10 ^A15-A19, ^B90, ^B200
    respiratory : ^J
    hiv_other   : ^B20-B24 excluding B200
    non_tb      : any other known code
    unknown     : no SIM code -> fall back to TBweb case_outcome
                  'Obito TB'  -> tb_via_tbweb
                  'Obito NTB' -> ntb_via_tbweb

    NOTE, matching 30e: respiratory and hiv_other deaths belong to NEITHER
    hybrid class.  They are therefore censored in both cause-specific
    analyses, as are unknown-cause deaths.
    """
    if CACHE.exists() and not force:
        return pd.read_csv(CACHE)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    if verbose:
        print(f"  building cause lookup from {RAW_CSV.name} (~98 MB, one-off)...")
    raw = pd.read_csv(RAW_CSV, low_memory=False,
                      usecols=["sinan_clean", "case_type", "case_outcome",
                               "end_date", "dod", "cause_of_death_code"])
    raw["end_date"] = pd.to_datetime(raw["end_date"], format="%B %d, %Y", errors="coerce")
    raw["dod"] = pd.to_datetime(raw["dod"], format="%B %d, %Y", errors="coerce")

    TRANSFER = ["Transf Outro Municipio", "Transf Outro Estado/Pais"]
    novo = raw[raw["case_type"].astype(str).str.strip().str.lower().eq("novo")].copy()
    co = novo["case_outcome"].astype(str).str.strip()
    novo = novo[co.ne("") & co.ne("nan") & co.ne("Mud Diag") & ~co.isin(TRANSFER)]
    novo = novo.sort_values("end_date")
    first = novo.drop_duplicates("sinan_clean", keep="first")[["sinan_clean", "case_outcome"]]

    # --- ADR-0003 FIX: Obito outcome from ANY episode -------------------
    # An LTFU patient's index episode closes as `Abandono`, so index-only
    # lookup cannot see their death. 1,058 of 1,668 LTFU deaths (63.4%) are
    # recorded on a retreatment episode (`Retr Aband`, `Recidiva`). Verified
    # same-death: median 0-day lag to the cohort death date.
    _ob = raw[raw["case_outcome"].astype(str).str.strip()
              .isin(["Obito TB", "Obito NTB"])]
    _ob = (_ob.sort_values("end_date")
              .drop_duplicates("sinan_clean", keep="last")
              [["sinan_clean", "case_outcome"]]
              .rename(columns={"case_outcome": "_obito"}))
    first = first.merge(_ob, on="sinan_clean", how="outer")
    first["case_outcome"] = first["_obito"].combine_first(first["case_outcome"])
    first = first.drop(columns=["_obito"])
    if verbose:
        print(f"  [ADR-0003] Obito recovered from any episode: {len(_ob):,}")

    dr = raw[raw["dod"].notna() & raw["cause_of_death_code"].notna()].copy()
    dr["cause_of_death_code"] = dr["cause_of_death_code"].astype(str).str.strip().str.upper()
    dr = dr[dr["cause_of_death_code"].ne("") & dr["cause_of_death_code"].ne("NAN")]
    dr = dr.sort_values("dod").drop_duplicates("sinan_clean", keep="last")
    dr = dr[["sinan_clean", "cause_of_death_code"]]

    a = first.merge(dr, on="sinan_clean", how="outer")
    cod = a["cause_of_death_code"].fillna("")
    known = cod.ne("")
    tb_strict = known & cod.str.match(r"^(A1[5-9]|B90|B200)")
    resp = known & cod.str.match(r"^J\d")
    hiv_other = known & cod.str.match(r"^B2[0-4]") & ~cod.str.match(r"^B200")

    cls = pd.Series("unknown", index=a.index)
    cls[known & ~tb_strict & ~resp & ~hiv_other] = "non_tb"
    cls[hiv_other] = "hiv_other"
    cls[resp] = "respiratory"
    cls[tb_strict] = "tb_strict"
    outc = a["case_outcome"].fillna("")
    cls[(cls == "unknown") & outc.eq("Obito TB")] = "tb_via_tbweb"
    cls[(cls == "unknown") & outc.eq("Obito NTB")] = "ntb_via_tbweb"

    a["cod_class"] = cls
    a["tb_hybrid"] = cls.isin(["tb_strict", "tb_via_tbweb"])
    a["nontb_hybrid"] = cls.isin(["non_tb", "ntb_via_tbweb"])
    a["tb_simonly"] = cls.eq("tb_strict")
    a["nontb_simonly"] = cls.eq("non_tb")

    out = a[["sinan_clean", "cod_class", "tb_hybrid", "nontb_hybrid",
             "tb_simonly", "nontb_simonly"]]
    out.to_csv(CACHE, index=False)
    if verbose:
        print(f"  cause classes: {dict(cls.value_counts())}")
    return out


# ---------------------------------------------------------------------------
# TIMELINE (one imputation)
# ---------------------------------------------------------------------------
_GEO_CACHE = {}


def build_geo_lookup(verbose=True):
    """sinan_clean -> geo4. Mirror of build_geo_lookup() in _rolling.R.

    Keyed on the RAW TBweb `tx_city` string via the crosswalk from
    47_build_ibge_typology.py, so no normalisation happens here. Custody is
    already folded into the reference level by that script, because
    geo_class == "Prison" is 100% collinear with the `incarcerated` covariate.

    City is taken from the INDEX episode. That does not hit invariant 8:
    `tx_city` is present on every episode, so there is no differential
    missingness. A modal/latest rule WOULD hit it (9.86% of LTFU patients have
    >1 city versus 1.35% of non-LTFU).
    """
    if "df" in _GEO_CACHE:
        return _GEO_CACHE["df"]
    xw_path = ROOT / "ITT_Analysis" / "external" / "municipality_typology_sp.csv"
    if not xw_path.exists():
        raise SystemExit("[geo] missing crosswalk; run "
                         "ITT_Analysis/scripts/47_build_ibge_typology.py")
    xw = pd.read_csv(xw_path)[["municipality", "geo4"]].rename(
        columns={"municipality": "tx_city"})

    raw = pd.read_csv(RAW_CSV, low_memory=False,
                      usecols=["sinan_clean", "case_type", "case_outcome",
                               "end_date", "tx_city"])
    raw["end_date"] = pd.to_datetime(raw.end_date, format="%B %d, %Y", errors="coerce")
    TRANSFER = ["Transf Outro Municipio", "Transf Outro Estado/Pais"]
    novo = raw[raw.case_type.astype(str).str.strip().str.lower().eq("novo")]
    novo = novo[novo.case_outcome.notna()
                & novo.case_outcome.astype(str).str.strip().ne("")
                & novo.case_outcome.ne("Mud Diag")
                & ~novo.case_outcome.isin(TRANSFER)]
    first = novo.sort_values("end_date").drop_duplicates("sinan_clean")[
        ["sinan_clean", "tx_city"]]
    anyc = raw.sort_values("end_date").drop_duplicates("sinan_clean")[
        ["sinan_clean", "tx_city"]].rename(columns={"tx_city": "tx_city_any"})
    first = first.merge(anyc, on="sinan_clean", how="outer")
    first["tx_city"] = first.tx_city.fillna(first.tx_city_any)

    out = first[["sinan_clean", "tx_city"]].merge(xw, on="tx_city", how="left")
    cov = out.geo4.notna().mean()
    if cov < 0.999:
        raise SystemExit(f"[geo] crosswalk covers only {cov:.2%} of patients; "
                         f"re-run 47_build_ibge_typology.py")
    if verbose:
        print(f"  [geo] {out.geo4.notna().sum():,} patients mapped ({cov:.2%})")
    out = out[["sinan_clean", "geo4"]]
    _GEO_CACHE["df"] = out
    return out


def load_timeline(imp_path, cause_lookup, ltfu_date_mode="shift30", verbose=False):
    df = pd.read_csv(imp_path, low_memory=False)
    df = df.merge(cause_lookup, on="sinan_clean", how="left")
    # always merged, so no call site can fit without a covariate that is in COVS
    df = df.merge(build_geo_lookup(verbose=False), on="sinan_clean", how="left")
    df["geo4"] = df.geo4.fillna(GEO_REF)
    for c in ["best_start", "end_date", "death_date"]:
        df[c] = pd.to_datetime(df[c], errors="coerce")

    is_ltfu = df["itt_group"].eq("Loss to follow-up").to_numpy()

    shift = 30 if ltfu_date_mode == "shift30" else 0
    t_dis = ((df["end_date"] - pd.Timedelta(days=shift)) - df["best_start"]
             ).dt.days.to_numpy(dtype="float64")
    t_dis[~is_ltfu] = np.inf
    n_neg = int(np.sum(np.isfinite(t_dis) & (t_dis < 0)))
    t_dis = np.where(np.isfinite(t_dis), np.maximum(t_dis, 0.0), np.inf)

    died = df["event_d"].to_numpy() == 1
    t_death = (df["death_date"] - df["best_start"]).dt.days.to_numpy(dtype="float64")
    t_death = np.where(died & np.isfinite(t_death), np.maximum(t_death, 0.0), np.inf)
    t_admin = (ADMIN_CENSOR - df["best_start"]).dt.days.to_numpy(dtype="float64")

    tl = pd.DataFrame({"is_ltfu": is_ltfu, "t_dis": t_dis,
                       "t_death": t_death, "t_admin": t_admin})
    tl["__any_death__"] = died
    for c in ["tb_hybrid", "nontb_hybrid", "tb_simonly", "nontb_simonly"]:
        tl[c] = df[c].fillna(False).to_numpy().astype(bool) & died
    for c in COVS:
        tl[c] = df[c].astype(str).to_numpy()
    # subgroup-only variables
    tl["resistance_clean"] = df["resistance_clean"].astype(str).to_numpy()
    tl["period"] = np.where(df["best_start"].dt.year <= 2019,
                            "2013-2019 (pre-COVID)", "2020-2023")

    if verbose:
        print(f"  {imp_path.name}: n={len(tl):,}  LTFU={is_ltfu.sum():,}  "
              f"deaths={int(died.sum()):,}  neg_disengage={n_neg:,}")
        for k in ["tb_hybrid", "nontb_hybrid", "tb_simonly", "nontb_simonly"]:
            print(f"      {k:<15s} {int(tl[k].sum()):>7,}")
    return tl


# ---------------------------------------------------------------------------
# CLONE AND CENSOR
# ---------------------------------------------------------------------------
def build_arm(tl, arm, T=None, cause="all_cause"):
    """One arm of clones for one outcome.

    Cause-specific handling (mirrors 30e): a death from a competing cause is
    treated as CENSORING at its own event time, as are unknown-cause deaths
    and (in the hybrid scheme) respiratory / HIV-other deaths.
    """
    t_dis = tl["t_dis"].to_numpy()
    t_death = tl["t_death"].to_numpy()
    t_admin = tl["t_admin"].to_numpy()

    ev_col = CAUSES[cause][0]
    is_target = tl[ev_col].to_numpy().astype(bool)
    any_death = tl["__any_death__"].to_numpy().astype(bool)

    # competing event: a death that is not the outcome of interest
    t_compete = np.where(any_death & ~is_target, t_death, np.inf)
    t_event = np.where(is_target, t_death, np.inf)

    with np.errstate(invalid="ignore"):
        m_dis = np.where(np.isfinite(t_dis), np.floor(t_dis / MONTH_DAYS), np.inf)
        m_event = np.where(np.isfinite(t_event), np.floor(t_event / MONTH_DAYS), np.inf)
        m_admin = np.floor(t_admin / MONTH_DAYS)
        m_compete = np.where(np.isfinite(t_compete),
                             np.floor(t_compete / MONTH_DAYS), np.inf)

    if arm == "remain":
        m_dev = m_dis
        death_first = t_event < t_dis                # day-level tie break
        at_risk_mode = "all"
    elif arm == "disengage":
        if DISJOINT:
            # Strategy: disengage DURING month T specifically (disjoint windows).
            # Nested strategies cannot express a timing gradient because the
            # exposed group grows with T; disjoint windows can, at the cost of a
            # much smaller adherent group and weaker positivity.
            adherent = m_dis == (T - 1)
            # people who already disengaged before the window deviated then;
            # people still engaged at the window's end deviate there
            m_dev = np.where(adherent, np.inf,
                             np.where(m_dis < (T - 1), m_dis, float(T)))
        else:
            adherent = m_dis <= (T - 1)
            m_dev = np.where(adherent, np.inf, float(T))
        death_first = m_event < m_dev                # month edge, no ties
        at_risk_mode = "window_end"
    else:
        raise ValueError(arm)

    m_cens = np.minimum(np.minimum(m_dev, m_compete),
                        np.minimum(m_admin, float(HORIZON_M)))
    has_event = (death_first & (t_event < t_admin) & (m_event < HORIZON_M)
                 & (t_event < t_compete))

    nmonths = np.where(has_event, m_event + 1, m_cens)
    nmonths = np.clip(nmonths, 0, HORIZON_M).astype(np.int32)

    dev_binds = ((~has_event) & np.isfinite(m_dev) & (m_dev <= m_admin)
                 & (m_dev <= m_compete) & (m_dev < HORIZON_M))

    keep = nmonths > 0
    return {"nmonths": nmonths[keep],
            "event_month": np.where(has_event, m_event, np.inf)[keep],
            "dev_binds": dev_binds[keep],
            "pat": tl["pat"].to_numpy()[keep],
            "at_risk_mode": at_risk_mode,
            "T": -1 if T is None else T}


def _ranges(n):
    """concatenated arange(k) for k in n, without a Python loop"""
    tot = int(n.sum())
    out = np.ones(tot, dtype=np.int32)
    out[0] = 0
    starts = np.cumsum(n)[:-1]
    out[starts] = 1 - n[:-1]
    return np.cumsum(out, dtype=np.int32)


def expand(arm):
    """Person-months as plain numpy arrays.

    Covariates are NOT expanded.  Every covariate here is measured at baseline
    and so is time-fixed; carrying only the integer pattern code avoids
    materializing millions of duplicated strings.
    """
    n = arm["nmonths"]
    month = _ranges(n)
    clone = np.repeat(np.arange(len(n), dtype=np.int64), n)
    ev = np.repeat(arm["event_month"], n)
    death = (np.isfinite(ev) & (month == ev)).astype(np.int8)
    last = month == (np.repeat(n, n) - 1)
    dev_next = (last & np.repeat(arm["dev_binds"], n) & (death == 0)).astype(np.int8)
    at_risk = (np.ones(len(month), dtype=bool) if arm["at_risk_mode"] == "all"
               else month == (arm["T"] - 1))
    return {"month": month, "clone": clone, "death": death,
            "dev_next": dev_next, "at_risk": at_risk,
            "pat": np.repeat(arm["pat"], n), "counts": n}


# ---------------------------------------------------------------------------
# STABILIZED IPCW
# ---------------------------------------------------------------------------
def _fit_cells(pat, month, y, Xpat):
    """Pooled logistic for the deviation hazard, fit on AGGREGATED cells.

    Every covariate is time-fixed, so a person-month's design row is fully
    determined by (covariate pattern, month).  Aggregating to those cells
    collapses ~3.4M rows to ~110k and is algebraically identical to fitting on
    the full expansion with binomial counts.  Returns per-row probabilities for
    the covariate model and the month-only (numerator) model.
    """
    cell = pat.astype(np.int64) * (HORIZON_M + 1) + month
    codes, uniq = pd.factorize(cell)
    n_cells = len(uniq)
    n_at = np.bincount(codes, minlength=n_cells)
    k_at = np.bincount(codes, weights=y, minlength=n_cells)
    cpat = (uniq // (HORIZON_M + 1)).astype(np.int64)
    cmon = (uniq % (HORIZON_M + 1)).astype(np.int64)

    Bm = csr_matrix((np.ones(n_cells), (np.arange(n_cells), cmon)),
                    shape=(n_cells, HORIZON_M + 1))
    Xd = hstack([Xpat[cpat], Bm], format="csr")

    def fit(X):
        X2 = vstack([X, X], format="csr")
        y2 = np.r_[np.ones(n_cells), np.zeros(n_cells)]
        w2 = np.r_[k_at, n_at - k_at]
        ok = w2 > 0
        m = LogisticRegression(max_iter=200, C=1.0).fit(X2[ok], y2[ok],
                                                        sample_weight=w2[ok])
        return m.predict_proba(X)[:, 1]

    return fit(Xd)[codes], fit(Bm)[codes]


def add_ipcw(rows, Xpat, verbose=False, label=""):
    n = len(rows["month"])
    p_den = np.zeros(n)
    p_num = np.zeros(n)
    fr = rows["at_risk"]
    y = rows["dev_next"][fr].astype(float)
    fitted = False
    if fr.sum() > 0 and 0 < y.sum() < len(y):
        p_den[fr], p_num[fr] = _fit_cells(rows["pat"][fr], rows["month"][fr],
                                          y, Xpat)
        fitted = True

    # Vectorized per-clone cumulative product, LAGGED one month: censoring at
    # the end of month k discounts months k+1 onward, not month k itself.
    lr = np.log((1.0 - p_num) / np.clip(1.0 - p_den, 1e-8, None))
    cs = np.cumsum(lr)
    cnt = rows["counts"]
    starts = np.r_[0, np.cumsum(cnt)[:-1]]
    base = np.repeat(cs[starts] - lr[starts], cnt)   # cumsum before clone start
    sw = np.exp(cs - base - lr)                      # exclusive => already lagged

    lo, hi = np.quantile(sw, WEIGHT_TRUNC)
    swt = np.clip(sw, lo, hi)
    if verbose:
        flag = "" if 0.9 <= sw.mean() <= 1.1 else "  <-- positivity CHECK"
        print(f"      w[{label:<22s}] mean={sw.mean():.3f} sd={sw.std():.3f} "
              f"max={sw.max():.2f} trunc={np.mean((sw < lo) | (sw > hi)):.2%}{flag}")
        if not fitted:
            print(f"      WARNING[{label}]: censoring model not identified "
                  f"(at-risk={int(fr.sum()):,}, events={int(y.sum())}); weights=1")
    rows["sw"], rows["swt"] = sw, swt
    return rows


def weighted_risk(rows, label="", strict=False):
    w, d, m = rows["swt"], rows["death"], rows["month"]
    num = np.bincount(m, weights=w * d, minlength=HORIZON_M)
    den = np.bincount(m, weights=w, minlength=HORIZON_M)
    hz = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    if strict and np.any(hz >= 1.0):
        raise AssertionError(f"{label}: weighted hazard >= 1 at month(s) "
                             f"{np.where(hz >= 1.0)[0].tolist()} -- month "
                             f"containing only events. Check month arithmetic.")
    return 1.0 - np.cumprod(1.0 - hz)


def risk_at(ci, month_paper):
    return float(ci[min(month_paper, HORIZON_M) - 1])


# ---------------------------------------------------------------------------
# CONTRASTS
# ---------------------------------------------------------------------------
def attach_patterns(tl, covs):
    """Factorize the baseline covariate pattern once per dataset.

    Returns (tl with a 'pat' code column, sparse dummy matrix over patterns).
    Both are reusable across strategies, causes AND bootstrap replicates: a
    resample changes pattern COUNTS but not the pattern universe, so the design
    matrix never has to be rebuilt.
    """
    tl = tl.copy()
    key = pd.MultiIndex.from_frame(tl[covs].astype(str))
    codes, uniq = pd.factorize(key)
    tl["pat"] = codes.astype(np.int64)
    Xpat = OneHotEncoder(handle_unknown="ignore", sparse_output=True
                         ).fit_transform(uniq.to_frame(index=False))
    return tl, csr_matrix(Xpat)


def reference(tl, cause, Xpat, verbose=False, strict=False):
    """'Remain engaged' depends on cause and dataset but NOT on T, so build it
    once per (dataset, cause) and reuse across strategies.  Rebuilding it per T
    roughly doubles the cost of every point estimate and every bootstrap rep."""
    r = add_ipcw(expand(build_arm(tl, "remain", cause=cause)), Xpat,
                 verbose=verbose, label=f"remain/{cause}")
    return weighted_risk(r, f"remain/{cause}", strict), int(r["death"].sum())


def one_contrast(tl, T, cause, Xpat, verbose=False, strict=False, ref=None):
    if ref is None:
        ref = reference(tl, cause, Xpat, verbose=verbose, strict=strict)
    ci0, ev_rem = ref
    Ld = add_ipcw(expand(build_arm(tl, "disengage", T=T, cause=cause)), Xpat,
                  verbose=verbose, label=f"T={T}/{cause}")
    ci1 = weighted_risk(Ld, f"T={T}/{cause}", strict)

    res = {"T": T, "cause": cause}
    for mo in REPORT_AT:
        r1, r0 = risk_at(ci1, mo), risk_at(ci0, mo)
        res[f"risk{mo}_dis"] = 100 * r1
        res[f"risk{mo}_rem"] = 100 * r0
        res[f"rr{mo}"] = r1 / r0 if r0 > 0 else np.nan
        res[f"rd{mo}"] = 100 * (r1 - r0)
    res["mean_weight"] = float(Ld["sw"].mean())
    res["estimable"] = bool(0.9 <= res["mean_weight"] <= 1.1)
    res["events_dis"] = int(Ld["death"].sum())
    res["events_rem"] = ev_rem
    return res


def pool_mi(rows):
    """Rubin point estimate across imputations: mean on the log scale for
    ratios, natural scale for risks and differences."""
    out = {k: rows[0][k] for k in ("T", "cause") if k in rows[0]}
    for k in rows[0]:
        if k in ("T", "cause"):
            continue
        v = np.array([r[k] for r in rows], dtype=float)
        if k.startswith("rr"):
            out[k] = float(np.exp(np.nanmean(np.log(v))))
        elif k == "estimable":
            out[k] = bool(np.all(v))
        else:
            out[k] = float(np.nanmean(v))
    return out


# ---------------------------------------------------------------------------
# BOOTSTRAP
# ---------------------------------------------------------------------------
def bootstrap(timelines, xpats, jobs, B, seed, rubin=False):
    """Person-level bootstrap.

    Default ('cycling'): each replicate resamples persons AND uses one
    imputation, cycling b % M.  One pass of B replicates propagates sampling
    and imputation uncertainty together.  Cost: B fits.

    --rubin: B replicates WITHIN each imputation, then Rubin's rules on
    log(RR)/RD -- total variance = mean(within) + (1 + 1/M) * between.
    Statistically cleaner, costs M x B fits.
    """
    keys = [f"{s}{mo}" for mo in REPORT_AT for s in ("rr", "rd")]
    M = len(timelines)
    rng = np.random.default_rng(seed)
    n = len(timelines[0])
    t0 = time.time()

    causes_in_jobs = sorted({c for _, c in jobs})

    if not rubin:
        draws = {j: {k: [] for k in keys} for j in jobs}
        for b in range(B):
            tl, Xp = timelines[b % M], xpats[b % M]
            d = tl.iloc[rng.integers(0, n, n)].reset_index(drop=True)
            refs = {c: reference(d, c, Xp) for c in causes_in_jobs}
            for j in jobs:
                T, cause = j
                r = one_contrast(d, T, cause, Xp, ref=refs[cause])
                for k in keys:
                    draws[j][k].append(r[k])
            if (b + 1) % 10 == 0:
                el = time.time() - t0
                print(f"    rep {b+1}/{B}  {el/(b+1):.1f}s/rep  "
                      f"ETA {el/(b+1)*(B-b-1)/60:.1f} min", flush=True)
        return {j: {k: (float(np.percentile(np.asarray(v)[np.isfinite(v)], 2.5)),
                        float(np.percentile(np.asarray(v)[np.isfinite(v)], 97.5)))
                    for k, v in d0.items()} for j, d0 in draws.items()}

    # ---- full Rubin ----
    per_imp = {j: {k: [] for k in keys} for j in jobs}     # point per imputation
    within = {j: {k: [] for k in keys} for j in jobs}      # bootstrap var per imputation
    for mi, (tl, Xp) in enumerate(zip(timelines, xpats)):
        refs0 = {c: reference(tl, c, Xp) for c in causes_in_jobs}
        pt = {j: one_contrast(tl, j[0], j[1], Xp, ref=refs0[j[1]]) for j in jobs}
        bd = {j: {k: [] for k in keys} for j in jobs}
        for b in range(B):
            d = tl.iloc[rng.integers(0, n, n)].reset_index(drop=True)
            refs = {c: reference(d, c, Xp) for c in causes_in_jobs}
            for j in jobs:
                r = one_contrast(d, j[0], j[1], Xp, ref=refs[j[1]])
                for k in keys:
                    bd[j][k].append(r[k])
            if (b + 1) % 10 == 0:
                el = time.time() - t0
                done = mi * B + b + 1
                tot = M * B
                print(f"    imp {mi+1}/{M} rep {b+1}/{B}  {el/done:.1f}s/rep  "
                      f"ETA {el/done*(tot-done)/60:.1f} min", flush=True)
        for j in jobs:
            for k in keys:
                v = np.asarray(bd[j][k], dtype=float)
                v = v[np.isfinite(v)]
                sc = np.log(v) if k.startswith("rr") else v
                within[j][k].append(float(np.var(sc, ddof=1)))
                p = pt[j][k]
                per_imp[j][k].append(float(np.log(p) if k.startswith("rr") else p))

    out = {}
    for j in jobs:
        out[j] = {}
        for k in keys:
            q = np.array(per_imp[j][k]); ub = np.mean(within[j][k])
            bvar = np.var(q, ddof=1) if len(q) > 1 else 0.0
            tot = ub + (1 + 1 / M) * bvar
            qb, se = float(np.mean(q)), float(np.sqrt(tot))
            lo, hi = qb - 1.96 * se, qb + 1.96 * se
            out[j][k] = (float(np.exp(lo)), float(np.exp(hi))) if k.startswith("rr") \
                else (lo, hi)
    return out


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
# mtime is meaningless in this tree -- Google Drive sync rewrites it on files
# nobody touched (CLAUDE.md invariant 3). So the ONLY reliable staleness signal
# is content. On 2026-08-18 `rolling_landmark_cause.csv` was found reporting a
# pre-primary-abandonment cohort after sitting stale for two days, because
# nothing recorded what had produced it. These hashes make the same failure
# detectable here: `tools/check_ccw_provenance.py` re-hashes the inputs and
# fails if they no longer match what the output claims.
def _sha16(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def build_provenance(imp_files):
    me = Path(__file__).resolve()
    return {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script": me.name,
        "script_sha16": _sha16(me),
        "raw_csv": str(RAW_CSV.relative_to(ROOT)),
        "raw_csv_sha16": _sha16(RAW_CSV) if RAW_CSV.exists() else None,
        "geo_crosswalk": "ITT_Analysis/external/municipality_typology_sp.csv",
        "geo_crosswalk_sha16": _sha16(
            ROOT / "ITT_Analysis" / "external" / "municipality_typology_sp.csv"),
        "cause_lookup": str(CACHE.relative_to(ROOT)),
        "cause_lookup_sha16": _sha16(CACHE) if CACHE.exists() else None,
        "imputations": {Path(f).name: _sha16(f) for f in imp_files},
        "adrs": "ADR-0003 cause-from-any-episode; ADR-0005 no conditional CCW",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=0)
    ap.add_argument("--rubin", action="store_true",
                    help="full M x B bootstrap with Rubin pooling (slower, cleaner)")
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"], default="shift30")
    ap.add_argument("--subgroups", action="store_true")
    ap.add_argument("--disjoint", action="store_true",
                    help="disengage DURING month T (timing gradient) "
                         "instead of nested by-month-T strategies")
    ap.add_argument("--cause", default="all_cause,tb_hybrid,nontb_hybrid",
                    help="comma-separated: " + ",".join(CAUSES))
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    global DISJOINT
    DISJOINT = args.disjoint
    OUTDIR.mkdir(parents=True, exist_ok=True)
    causes = [c.strip() for c in args.cause.split(",")]
    for c in causes:
        assert c in CAUSES, f"unknown cause {c}"

    print("=" * 78)
    print("CCW v3 -- nested 'disengage by month T' vs remain engaged")
    print(f"  ltfu_date={args.ltfu_date}  imputations={args.imputations}  "
          f"causes={causes}")
    print("=" * 78)

    lookup = build_cause_lookup()
    imp_files = sorted(MI_DIR.glob("imp_*.csv"))[:args.imputations]
    assert imp_files, f"no imputations in {MI_DIR}"
    print("\n--- imputed datasets ---")
    tls = [load_timeline(p, lookup, args.ltfu_date, verbose=(i == 0))
           for i, p in enumerate(imp_files)]
    pairs = [attach_patterns(tl, COVS) for tl in tls]
    timelines = [a for a, _ in pairs]
    xpats = [b for _, b in pairs]
    print(f"  covariate patterns: {xpats[0].shape[0]:,} distinct "
          f"({xpats[0].shape[1]} dummies)")

    # ---------------- main contrasts ----------------
    print("\n--- weight diagnostics (imputation 1) ---")
    for cause in causes:
        one_contrast(timelines[0], PRIMARY_T, cause, xpats[0],
                     verbose=True, strict=True)

    rows = []
    for cause in causes:
        refs = [reference(tl, cause, Xp) for tl, Xp in zip(timelines, xpats)]
        for T in STRATEGY_MONTHS:
            rows.append(pool_mi([one_contrast(tl, T, cause, Xp, ref=rf)
                                 for tl, Xp, rf in zip(timelines, xpats, refs)]))
    res = pd.DataFrame(rows)

    for cause in causes:
        sub = res[res["cause"] == cause]
        print(f"\n--- {CAUSES[cause][1]} ---")
        print(f"    reference (remain engaged) 24-mo risk: "
              f"{sub['risk24_rem'].iloc[0]:.2f}%")
        hdr = (f"{'T':>2} {'6mo RR':>8} {'6mo RD':>8} {'24mo RR':>8} "
               f"{'24mo RD':>8} {'risk24':>8} {'events':>8} {'ok':>4}")
        print(hdr); print("-" * len(hdr))
        for _, r in sub.iterrows():
            print(f"{int(r['T']):>2} {r['rr6']:>8.2f} {r['rd6']:>+8.2f} "
                  f"{r['rr24']:>8.2f} {r['rd24']:>+8.2f} "
                  f"{r['risk24_dis']:>7.2f}% {int(r['events_dis']):>8,} "
                  f"{'y' if r['estimable'] else 'NO':>4}")

    out = {"provenance": build_provenance(imp_files),
           "config": {"ltfu_date": args.ltfu_date, "horizon_m": HORIZON_M,
                      "primary_T": PRIMARY_T, "covariates": COVS,
                      "imputations": len(imp_files), "causes": causes,
                      "weight_truncation": WEIGHT_TRUNC,
                      "month_convention": "paper (1-indexed)"},
           "main": res.to_dict(orient="records")}

    # ---------------- subgroups ----------------
    if args.subgroups:
        print(f"\n{'='*78}\nSUBGROUPS (strategy T={PRIMARY_T}, "
              f"all-cause unless noted)\n{'='*78}")
        srows = []
        for sg in SUBGROUPS:
            covs_sg = [c for c in COVS if c != sg]
            levels = sorted(set(timelines[0][sg].tolist()))
            levels = [l for l in levels if l not in ("nan", "", "__missing__")]
            print(f"\n--- {sg} ---")
            hdr = (f"{'level':<24s} {'n':>8} {'24mo RR':>8} {'24mo RD':>8} "
                   f"{'risk24':>8} {'ok':>4}")
            print(hdr); print("-" * len(hdr))
            for lvl in levels:
                per_imp = []
                for tl in timelines:
                    s = tl[tl[sg].to_numpy() == lvl]
                    if len(s) < 500 or int(s["__any_death__"].sum()) < 20:
                        continue
                    s, Xs = attach_patterns(s.reset_index(drop=True), covs_sg)
                    per_imp.append(one_contrast(
                        s, PRIMARY_T, "all_cause", Xs,
                        ref=reference(s, "all_cause", Xs)))
                if not per_imp:
                    print(f"{lvl:<24s} {'--':>8}  (too few for stable estimation)")
                    continue
                r = pool_mi(per_imp)
                r.update({"subgroup": sg, "level": lvl,
                          "n": int((timelines[0][sg].to_numpy() == lvl).sum())})
                srows.append(r)
                print(f"{lvl:<24s} {r['n']:>8,} {r['rr24']:>8.2f} "
                      f"{r['rd24']:>+8.2f} {r['risk24_dis']:>7.2f}% "
                      f"{'y' if r['estimable'] else 'NO':>4}")
        if srows:
            sres = pd.DataFrame(srows)
            sres.to_csv(OUTDIR / f"ccw_v3_{args.ltfu_date}_subgroups.csv", index=False)
            out["subgroups"] = sres.to_dict(orient="records")
            print("\n  Reminder: relative and absolute measures rank subgroups")
            print("  differently. Read RD, not RR, for where the preventable")
            print("  deaths are -- this is the point Reviewer 1 asks you to make.")

    # ---------------- bootstrap ----------------
    if args.bootstrap:
        jobs = [(T, c) for c in causes for T in STRATEGY_MONTHS]
        mode = "Rubin M x B" if args.rubin else "cycling imputations"
        print(f"\n{'='*78}\nBOOTSTRAP: {args.bootstrap} reps, {mode}\n{'='*78}")
        ci = bootstrap(timelines, xpats, jobs, args.bootstrap,
                       args.seed, args.rubin)
        out["ci_mode"] = mode
        out["ci"] = {f"{c}|T{T}": {k: list(v) for k, v in ci[(T, c)].items()}
                     for (T, c) in jobs}
        for cause in causes:
            print(f"\n--- {CAUSES[cause][1]} ---")
            for T in STRATEGY_MONTHS:
                pt = res[(res.cause == cause) & (res["T"] == T)].iloc[0]
                c = ci[(T, cause)]
                print(f"  T={T}: 24mo RR {pt['rr24']:.2f} "
                      f"({c['rr24'][0]:.2f}-{c['rr24'][1]:.2f})   "
                      f"RD {pt['rd24']:+.2f}pp "
                      f"({c['rd24'][0]:+.2f} to {c['rd24'][1]:+.2f})   "
                      f"6mo RR {pt['rr6']:.2f} "
                      f"({c['rr6'][0]:.2f}-{c['rr6'][1]:.2f})")

    tag = args.ltfu_date
    (OUTDIR / f"ccw_v3_{tag}.json").write_text(json.dumps(out, indent=2))
    res.to_csv(OUTDIR / f"ccw_v3_{tag}_main.csv", index=False)
    print(f"\nSaved -> {OUTDIR}/ccw_v3_{tag}.json (+ _main.csv)")


if __name__ == "__main__":
    main()
