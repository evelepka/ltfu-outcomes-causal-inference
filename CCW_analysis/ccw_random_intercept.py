#!/usr/bin/env python3
"""CCW with a municipality random intercept in the censoring-weight model.

WHY THIS EXISTS
---------------
Reviewer 1 (comment 1, R2 round) asks that non-independence of patients within
a municipality be handled directly. `ccw_cluster_bootstrap.py` answers the
VARIANCE half by resampling municipalities; this answers the BIAS half by
putting municipality inside the model that builds the weights.

WHAT IS AND IS NOT A RANDOM INTERCEPT HERE
------------------------------------------
There is no outcome regression in this estimator: risk comes from a weighted
nonparametric hazard (`ccw_v3.weighted_risk`), and the IPCW does all the
adjustment. So "random intercept in the outcome model" has nothing to attach
to; the intercept lives in the censoring-weight model alone. That is the model
that produced RD = 2.4931 pp at 60 months.

The intercept is imposed as an L2 penalty on the municipality indicators only.
`LogisticRegression` applies one C to every column, so the municipality block is
instead SCALED by sigma: scaling a column by s divides the effective penalty on
its coefficient by s^2, and with the estimator's C = 1 that makes the penalty on
a municipality coefficient 1/sigma^2 -- a Gaussian random intercept of variance
sigma^2. sigma^2 = 0.3009 is the ML estimate from
`glmer(is_ltfu ~ 14 covariates + (1|muni))` (lme4 1.1.37, R 4.5.1).

Municipality is taken from the INDEX episode, the rule `build_geo_lookup` uses
for geo4 (invariant 5). NOTE: `ccw_cluster_bootstrap.py` instead takes the first
row per patient; that affects only which block a patient falls in, not any point
estimate, but the two should be reconciled before publication.

BOOTSTRAP UNITS
---------------
  patient    resample patients, as `ccw_main_bootstrap.py`. Correct if the
             estimand is the effect IN Sao Paulo State: the 640 municipalities
             are a census, not a sample.
  cluster    resample municipalities as whole blocks. Conservative; generalises
             beyond the state.
  certainty  Sao Paulo city enters every replicate (its patients resampled with
             replacement) and the other 639 are resampled as whole blocks. The
             capital holds 33.5% of the cohort and there are only ~7.6 effective
             clusters, so a plain cluster bootstrap is largely a lottery over how
             often the capital is drawn; treating it as a certainty unit removes
             that lottery while keeping between-municipality variation.

Usage:
    python3 ccw_random_intercept.py --boot-unit patient --bootstrap 1000
Output: results_v3/<--out>, never overwriting an existing file.
"""
from __future__ import annotations
import argparse, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix, diags
from sklearn.preprocessing import OneHotEncoder

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ccw_v3 as C

SCRIPT_VERSION = "2026-10-16"
SIGMA2_DEFAULT = 0.3009          # glmer ML estimate, see module docstring


def municipality(imp_path: Path, n: int) -> np.ndarray:
    """tx_city from the INDEX episode -- the rule in C.build_geo_lookup."""
    raw = pd.read_csv(C.RAW_CSV, low_memory=False,
                      usecols=["sinan_clean", "case_type", "case_outcome",
                               "end_date", "tx_city"])
    raw["end_date"] = pd.to_datetime(raw.end_date, format="%B %d, %Y",
                                     errors="coerce")
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
    ids = pd.read_csv(imp_path, low_memory=False, usecols=["sinan_clean"])
    city = ids.merge(first[["sinan_clean", "tx_city"]],
                     on="sinan_clean", how="left")["tx_city"]
    assert len(city) == n, "municipality lookup lost rows"
    return city.fillna("__NA__").to_numpy()


def patterns(tl: pd.DataFrame, sigma: float | None, covs=None):
    """Covariate-pattern codes and the design matrix.

    sigma=None reproduces the published model (no municipality). `covs` lets a
    subgroup run drop the stratifying variable, as `ccw_subgroup_bootstrap` does.
    """
    base = list(covs) if covs is not None else list(C.COVS)
    if sigma is None:
        return C.attach_patterns(tl, base)
    covs = base + ["muni"]
    tl = tl.copy()
    key = pd.MultiIndex.from_frame(tl[covs].astype(str))
    codes, uniq = pd.factorize(key)
    tl["pat"] = codes.astype(np.int64)
    enc = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    X = enc.fit_transform(uniq.to_frame(index=False))
    start = sum(len(c) for c in enc.categories_[:-1])   # municipality is last
    s = np.ones(X.shape[1]); s[start:] = sigma
    return tl, csr_matrix(X @ diags(s))


def contrast(tl, Xp, T, H, cause="all_cause"):
    curves = []
    for arm, t in (("disengage", T), ("remain", None)):
        r = C.add_ipcw(C.expand(C.build_arm(tl, arm, T=t, cause=cause)), Xp)
        w, d, m = r["swt"], r["death"], r["month"]
        num = np.bincount(m, weights=w * d, minlength=H)
        den = np.bincount(m, weights=w, minlength=H)
        hz = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        curves.append(1 - np.cumprod(1 - hz))
    ci1, ci0 = curves
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(ci0 > 0, ci1 / ci0, np.nan)
    return 100 * (ci1 - ci0), rr


def blocks(codes: np.ndarray, n_groups: int):
    order = np.argsort(codes, kind="stable")
    bd = np.searchsorted(codes[order], np.arange(n_groups + 1))
    return [order[bd[k]:bd[k + 1]] for k in range(n_groups)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot-unit", choices=["patient", "cluster", "certainty"],
                    default="patient")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--horizon", type=int, default=60)
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"],
                    default="closure")
    ap.add_argument("--sigma2", type=float, default=SIGMA2_DEFAULT)
    ap.add_argument("--no-random-intercept", action="store_true",
                    help="published model; for Run B, where only the "
                         "resampling unit changes")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    H, T, tag = a.horizon, C.PRIMARY_T, a.ltfu_date
    C.HORIZON_M = H
    C.REPORT_AT = sorted(set(C.REPORT_AT + [H]))
    sigma = None if a.no_random_intercept else float(np.sqrt(a.sigma2))
    model = "published" if sigma is None else "muni_random_intercept"

    name = a.out or (f"ccw_re_patient_bootstrap_{tag}_h{H}.csv"
                     if a.boot_unit == "patient" and sigma is not None else
                     f"ccw_{model}_{a.boot_unit}_bootstrap_{tag}_h{H}.csv")
    dest = C.OUTDIR / name
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")

    started = datetime.now(timezone.utc)
    print(f"ccw_random_intercept.py {SCRIPT_VERSION} | start "
          f"{started.isoformat(timespec='seconds')}")
    print(f"  model={model} sigma2={a.sigma2 if sigma else float('nan')} "
          f"unit={a.boot_unit} B={a.bootstrap} M={a.imputations} "
          f"horizon={H} T={T} tag={tag} seed={a.seed}")

    lookup = C.build_cause_lookup()
    imps = sorted(C.MI_DIR.glob("imp_*.csv"))[:a.imputations]
    tls, xps, groups, capital = [], [], [], []
    for i, p in enumerate(imps):
        tl = C.load_timeline(p, lookup, tag, verbose=(i == 0))
        tl["muni"] = municipality(p, len(tl))
        codes, uniq = pd.factorize(tl["muni"])
        tl, Xp = patterns(tl, sigma)
        tls.append(tl); xps.append(Xp)
        gp = blocks(codes, len(uniq))
        groups.append(gp)
        big = int(np.argmax([len(g) for g in gp]))
        capital.append(big)
        if i == 0:
            sh = np.array([len(g) for g in gp], float); sh /= sh.sum()
            print(f"  municipalities {len(uniq)}  effective "
                  f"{1/np.sum(sh**2):.1f}  largest {100*sh.max():.1f}% "
                  f"({uniq[big]})")
            print(f"  design matrix: {Xp.shape[0]:,} patterns x "
                  f"{Xp.shape[1]:,} columns")
    M = len(tls)

    rds, rrs = zip(*[contrast(tl, Xp, T, H) for tl, Xp in zip(tls, xps)])
    rd_pt = np.mean(rds, axis=0)
    rr_pt = np.exp(np.nanmean(np.log(rrs), axis=0))
    print(f"  point {H} m: RD {rd_pt[H-1]:+.4f} pp   RR {rr_pt[H-1]:.4f}",
          flush=True)

    # validation gate: the published model must reproduce the published point
    if sigma is None:
        ref = C.OUTDIR / f"ccw_main_bootstrap_{tag}_h{H}.csv"
        if ref.exists():
            want = pd.read_csv(ref).iloc[H - 1]["rd"]
            if abs(rd_pt[H - 1] - want) > 1e-6:
                raise SystemExit(f"point estimate {rd_pt[H-1]:.6f} does not "
                                 f"reproduce {want:.6f} in {ref.name}")
            print(f"  validation gate OK (reproduces {want:.4f})")

    rng = np.random.default_rng(a.seed)
    dr = np.full((a.bootstrap, H), np.nan)
    br = np.full((a.bootstrap, H), np.nan)
    t0 = time.time()
    for b in range(a.bootstrap):
        k = b % M
        tl, Xp, gp = tls[k], xps[k], groups[k]
        if a.boot_unit == "patient":
            idx = rng.integers(0, len(tl), len(tl))
        elif a.boot_unit == "cluster":
            pick = rng.integers(0, len(gp), len(gp))
            idx = np.concatenate([gp[j] for j in pick])
        else:                                     # certainty unit
            cap = gp[capital[k]]
            others = [g for j, g in enumerate(gp) if j != capital[k]]
            pick = rng.integers(0, len(others), len(others))
            idx = np.concatenate(
                [cap[rng.integers(0, len(cap), len(cap))]]
                + [others[j] for j in pick])
        try:
            dr[b], br[b] = contrast(tl.iloc[idx].reset_index(drop=True), Xp, T, H)
        except Exception as e:                    # noqa: BLE001
            print(f"  rep {b+1} failed: {e}")
        if (b + 1) % 50 == 0:
            el = time.time() - t0
            print(f"  rep {b+1}/{a.bootstrap} | {el/(b+1):.1f}s/rep | "
                  f"ETA {el/(b+1)*(a.bootstrap-b-1)/60:.0f} min", flush=True)

    runtime = time.time() - t0
    ok = ~np.isnan(dr[:, H - 1])
    rd_lo, rd_hi = np.percentile(dr[ok], [2.5, 97.5], axis=0)
    rr_lo, rr_hi = np.nanpercentile(br[ok], [2.5, 97.5], axis=0)
    out = pd.DataFrame({
        "month": np.arange(1, H + 1), "rd": rd_pt, "rd_lo": rd_lo,
        "rd_hi": rd_hi, "rr": rr_pt, "rr_lo": rr_lo, "rr_hi": rr_hi,
        "crosses_zero": (rd_lo <= 0) & (rd_hi >= 0), "B_ok": int(ok.sum()),
        "ltfu_date": tag, "T": T, "cause": "all_cause", "n_imp": M,
        "model": model, "sigma2": (a.sigma2 if sigma is not None else np.nan),
        "boot_unit": a.boot_unit, "script_version": SCRIPT_VERSION,
        "run_started_utc": started.isoformat(timespec="seconds"),
        "runtime_s": round(runtime, 1)})
    out.to_csv(dest, index=False)

    r = out.iloc[H - 1]
    print(f"\n  {model} / {a.boot_unit}  {H} m: "
          f"RD {r.rd:+.4f} ({r.rd_lo:+.4f}, {r.rd_hi:+.4f})  "
          f"width {r.rd_hi - r.rd_lo:.4f}")
    print(f"  RR {r.rr:.4f} ({r.rr_lo:.4f}, {r.rr_hi:.4f})   "
          f"B used {int(r.B_ok)}   runtime {runtime/60:.1f} min")
    print(f"  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
