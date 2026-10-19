#!/usr/bin/env python3
"""CCW subgroup risk differences under the municipality random-intercept model.

Companion to `ccw_subgroup_bootstrap.py`, which produces the published subgroup
numbers. This reruns the same 15 strata with a municipality random intercept in
the censoring-weight model, so the two can be compared cell by cell.

HOW IT REUSES THE EXISTING CODE
-------------------------------
`subgroup_points` and `subgroup_replicate` in `ccw_subgroup_bootstrap` both call
`ccw_v3.attach_patterns(stratum, covs_without_the_stratifier)`. Patching that one
function is therefore enough to put municipality into every stratum's weight
model, and nothing else about the subgroup machinery changes -- there remains one
implementation of the estimator.

VALIDATION GATE
---------------
The published gate checks the main effect against 2.4493. Under this model the
main effect is 2.4931 (ccw_re_patient_bootstrap_closure_h60.csv), so the gate is
re-pointed; a mismatch means the municipality column or the shrinkage is wrong.

Usage:
    python3 ccw_subgroup_random_intercept.py --bootstrap 300
Output: results_v3/ccw_subgroups_re_h60_bootstrap.csv (never overwritten)
"""
from __future__ import annotations
import argparse, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ccw_v3 as C
import ccw_random_intercept as RI
import ccw_subgroup_bootstrap as SB

SCRIPT_VERSION = "2026-10-17"
PUBLISHED_RE = {"rd": 2.4931, "rr": 1.2680}      # Run A point estimate
TOL = 5e-3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=300)
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--horizon", type=int, default=60)
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"], default="closure")
    ap.add_argument("--sigma2", type=float, default=RI.SIGMA2_DEFAULT)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="ccw_subgroups_re_h60_bootstrap.csv")
    a = ap.parse_args()

    H = a.horizon
    SB.HORIZON = H
    SB.set_horizon(H)
    sigma = float(np.sqrt(a.sigma2))
    dest = C.OUTDIR / a.out
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")

    started = datetime.now(timezone.utc)
    print(f"ccw_subgroup_random_intercept.py {SCRIPT_VERSION} | "
          f"start {started.isoformat(timespec='seconds')}")
    print(f"  sigma2={a.sigma2} B={a.bootstrap} M={a.imputations} horizon={H} "
          f"tag={a.ltfu_date} seed={a.seed}")

    lookup = C.build_cause_lookup()
    imps = sorted(C.MI_DIR.glob("imp_*.csv"))[:a.imputations]
    timelines = []
    for i, p in enumerate(imps):
        tl = C.load_timeline(p, lookup, a.ltfu_date, verbose=(i == 0))
        tl["muni"] = RI.municipality(p, len(tl))
        timelines.append(tl)

    # every stratum's weight model now carries the municipality intercept
    C.attach_patterns = lambda tl, covs: RI.patterns(tl, sigma, covs=covs)

    xpats = [RI.patterns(tl, sigma)[1] for tl in timelines]
    timelines = [RI.patterns(tl, sigma)[0] for tl in timelines]

    me = SB.main_effect(timelines, xpats)
    rd, rr = me[f"rd{H}"], me[f"rr{H}"]
    print(f"\n  main effect {H} m: RD {rd:+.4f}  RR {rr:.4f}   "
          f"(expected {PUBLISHED_RE['rd']:+.4f} / {PUBLISHED_RE['rr']:.4f})")
    if abs(rd - PUBLISHED_RE["rd"]) > TOL:
        raise SystemExit("does not reproduce the random-intercept main effect; "
                         "stopping rather than producing untrustworthy subgroups")
    print("  validation gate OK")

    print("\n  subgroup point estimates...")
    pts = SB.subgroup_points(timelines, xpats)
    for r in pts:
        print(f"    {r['subgroup']:<18} {r['level']:<22} "
              f"RD {r[f'rd{H}']:+7.3f}  RR {r[f'rr{H}']:.3f}  n={r['n']:,}")

    sg_levels = [(r["subgroup"], r["level"]) for r in pts]
    rng = np.random.default_rng(a.seed)
    M = len(timelines)
    draws = {k: {"rr": [], "rd": []} for k in sg_levels}
    t0 = time.time()
    for b in range(a.bootstrap):
        res = SB.subgroup_replicate(timelines[b % M], sg_levels, rng)
        for k, (vrr, vrd) in res.items():
            draws[k]["rr"].append(vrr)
            draws[k]["rd"].append(vrd)
        if (b + 1) % 10 == 0:
            el = time.time() - t0
            print(f"    rep {b+1}/{a.bootstrap}  {el/(b+1):.1f}s/rep  "
                  f"ETA {el/(b+1)*(a.bootstrap-b-1)/60:.1f} min", flush=True)
    runtime = time.time() - t0

    rows = []
    for r in pts:
        k = (r["subgroup"], r["level"])
        rdv = np.array(draws[k]["rd"], float)
        rrv = np.array(draws[k]["rr"], float)
        lo, hi = (np.percentile(rdv, [2.5, 97.5]) if len(rdv) else (np.nan, np.nan))
        rlo, rhi = (np.nanpercentile(rrv, [2.5, 97.5]) if len(rrv) else (np.nan, np.nan))
        rows.append({"subgroup": r["subgroup"], "level": r["level"], "n": r["n"],
                     "rd": r[f"rd{H}"], "rr": r[f"rr{H}"],
                     "risk_dis": r[f"risk{H}_dis"], "risk_rem": r[f"risk{H}_rem"],
                     "estimable": r["estimable"], "B_ok": int(len(rdv)),
                     "rd_lo": lo, "rd_hi": hi, "rd_boot_mean": float(np.mean(rdv)) if len(rdv) else np.nan,
                     "rr_lo": rlo, "rr_hi": rhi, "rr_boot_mean": float(np.mean(rrv)) if len(rrv) else np.nan,
                     "model": "muni_random_intercept", "sigma2": a.sigma2,
                     "boot_unit": "patient", "script_version": SCRIPT_VERSION,
                     "run_started_utc": started.isoformat(timespec="seconds"),
                     "runtime_s": round(runtime, 1)})
    pd.DataFrame(rows).to_csv(dest, index=False)
    print(f"\n  wrote {dest}   runtime {runtime/60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
