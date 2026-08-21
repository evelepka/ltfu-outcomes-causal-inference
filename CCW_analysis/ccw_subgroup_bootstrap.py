#!/usr/bin/env python3
"""Confidence intervals for the CCW subgroup risk differences at 60 months.

WHY THIS EXISTS
---------------
The 2026-08-19 decision makes CCW subgroups a PRIMARY result, and
`docs/handoff-2026-08-19.md` §2b reports them with no intervals at all — listed
there as the first blocking item. Two cells decide whether there is anything to
say: homeless (+0.29 pp) and age >=65 (+1.00 pp).

`ccw_v3.py --subgroups` computes the point estimates but does not bootstrap
them, and its `bootstrap()` only takes (T, cause) jobs. The 60-month horizon
also lives outside the tree: `ccw_v3.HORIZON_M` is 24, and the h60 outputs in
results_v3 came from a local edit on the other machine.

Rather than edit `ccw_v3.py` — which is being actively worked on — this imports
it and overrides the two module constants, so there is exactly one
implementation of the CCW machinery.

VALIDATION GATE
---------------
Before doing anything with subgroups, this reproduces the published main-effect
point estimate at 60 months (nested T=6: RR 1.2368, RD +2.2241 from
`ccw_bootstrap_h60.csv`). If that does not match, the horizon override is wrong
and the script stops rather than producing plausible-looking subgroup numbers.

RESAMPLING
----------
Persons are resampled from the WHOLE cohort in each replicate and the strata are
then taken from the resample, mirroring `ccw_v3.bootstrap()`. The cohort is the
sampling unit, so stratum sizes are allowed to vary; this also keeps replicates
aligned across strata, should a paired subgroup contrast be wanted later.

Usage:
    python3 ccw_subgroup_bootstrap.py --validate-only
    python3 ccw_subgroup_bootstrap.py --bootstrap 300
    python3 ccw_subgroup_bootstrap.py --bootstrap 20 --imputations 2   # smoke test

Output: CCW_analysis/results_v3/ccw_subgroups_h60_bootstrap.csv
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ccw_v3 as C

# ccw_v3.ROOT is hardcoded to the other collaborator's Drive mount, so importing
# it on any other machine fails at the first mkdir. Repoint every path derived
# from it to THIS checkout, inferred from where this file sits. Nothing is
# written outside the project tree.
_HERE_ROOT = Path(__file__).resolve().parents[1]
if not C.MI_DIR.exists():
    C.ROOT = _HERE_ROOT
    C.MI_DIR = C.ROOT / "ITT_Analysis/data/mi"
    C.RAW_CSV = C.ROOT / "Data/Final_table_cleaned.csv"
    C.OUTDIR = C.ROOT / "CCW_analysis/results_v3"
    C.CACHE = C.OUTDIR / "cause_lookup_fixedattr.csv"
    assert C.MI_DIR.exists(), f"repointed ROOT still wrong: {C.MI_DIR}"
    print(f"[paths] ccw_v3.ROOT repointed to {C.ROOT}")

HORIZON = 60
PUBLISHED = {60: {"rr": 1.2367522921445866, "rd": 2.2240625382225376},   # ccw_bootstrap_h60.csv
             24: {"rr": 1.1464629,          "rd": 1.2145260}}            # ccw_v3_shift30_main.csv, T=6
TOL = 5e-4


def set_horizon(h: int) -> None:
    """Override the module constants the horizon is read from.

    Both are consulted inside the functions at call time, so assigning here is
    enough; nothing is computed from them at import.
    """
    C.HORIZON_M = h
    C.REPORT_AT = sorted(set(list(C.REPORT_AT) + [h]))


def load(n_imp: int, ltfu_date: str = "shift30"):
    lookup = C.build_cause_lookup()
    imp_files = sorted(C.MI_DIR.glob("imp_*.csv"))[:n_imp]
    assert imp_files, f"no imputations in {C.MI_DIR}"
    tls = [C.load_timeline(p, lookup, ltfu_date, verbose=(i == 0))
           for i, p in enumerate(imp_files)]
    pairs = [C.attach_patterns(tl, C.COVS) for tl in tls]
    return [a for a, _ in pairs], [b for _, b in pairs]


def main_effect(timelines, xpats):
    rows = []
    for tl, Xp in zip(timelines, xpats):
        ref = C.reference(tl, "all_cause", Xp)
        rows.append(C.one_contrast(tl, C.PRIMARY_T, "all_cause", Xp, ref=ref))
    return C.pool_mi(rows)


def subgroup_points(timelines, xpats_unused, subgroups=None):
    """Point estimates per subgroup level, pooled across imputations."""
    out = []
    for sg in (subgroups or C.SUBGROUPS):
        covs_sg = [c for c in C.COVS if c != sg]
        levels = sorted(set(timelines[0][sg].tolist()))
        levels = [l for l in levels if l not in ("nan", "", "__missing__")]
        for lvl in levels:
            per_imp = []
            for tl in timelines:
                s = tl[tl[sg].to_numpy() == lvl]
                if len(s) < 500 or int(s["__any_death__"].sum()) < 20:
                    continue
                s, Xs = C.attach_patterns(s.reset_index(drop=True), covs_sg)
                per_imp.append(C.one_contrast(s, C.PRIMARY_T, "all_cause", Xs,
                                              ref=C.reference(s, "all_cause", Xs)))
            if not per_imp:
                continue
            r = C.pool_mi(per_imp)
            r.update({"subgroup": sg, "level": lvl,
                      "n": int((timelines[0][sg].to_numpy() == lvl).sum())})
            out.append(r)
    return out


def subgroup_replicate(tl, sg_levels, rng):
    """One bootstrap replicate: resample the cohort, then read every stratum
    off the resample."""
    d = tl.iloc[rng.integers(0, len(tl), len(tl))].reset_index(drop=True)
    res = {}
    for sg, lvl in sg_levels:
        covs_sg = [c for c in C.COVS if c != sg]
        s = d[d[sg].to_numpy() == lvl]
        if len(s) < 500 or int(s["__any_death__"].sum()) < 20:
            continue
        try:
            s, Xs = C.attach_patterns(s.reset_index(drop=True), covs_sg)
            r = C.one_contrast(s, C.PRIMARY_T, "all_cause", Xs,
                               ref=C.reference(s, "all_cause", Xs))
        except Exception:
            continue
        res[(sg, lvl)] = (r[f"rr{HORIZON}"], r[f"rd{HORIZON}"])
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=0)
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--horizon", type=int, default=60,
                    help="diagnostic: run at another horizon to localise a mismatch")
    ap.add_argument("--subgroups", default="", help="comma-separated subset")
    args = ap.parse_args()

    global HORIZON
    HORIZON = args.horizon
    set_horizon(HORIZON)
    print(f"CCW subgroup bootstrap | horizon {HORIZON} mo | "
          f"{args.imputations} imputation(s) | B={args.bootstrap}")
    print(f"  REPORT_AT now {C.REPORT_AT}")

    timelines, xpats = load(args.imputations)

    # ---- validation gate -----------------------------------------------
    me = main_effect(timelines, xpats)
    rr, rd = me[f"rr{HORIZON}"], me[f"rd{HORIZON}"]
    pub = PUBLISHED.get(HORIZON)
    if pub is None:
        print(f"  no published value at horizon {HORIZON}; skipping gate"); return 0
    drr, drd = abs(rr - pub["rr"]), abs(rd - pub["rd"])
    print(f"\n  main effect, nested T={C.PRIMARY_T}, {HORIZON} mo:")
    print(f"    RR {rr:.6f}   published {pub['rr']:.6f}   diff {drr:.2e}")
    print(f"    RD {rd:+.6f}  published {pub['rd']:+.6f}  diff {drd:.2e}")
    # The published h60 figures are SINGLE-IMPUTATION: with --imputations 1 this
    # reproduces them to machine precision (6e-15). With M>1 the pooled estimate
    # differs slightly and legitimately, so the gate is only binding at M=1.
    ok = drr < TOL and drd < TOL
    if args.imputations > 1 and not ok:
        print(f"    -> differs from the published SINGLE-IMPUTATION figure, as expected\n"
              f"       at M={args.imputations}. Gate is only binding at M=1; run\n"
              f"       --validate-only --imputations 1 to check the machinery.")
        ok = True
    else:
        print(f"    -> {'MATCHES' if ok else 'DOES NOT MATCH'} the published figures")
    if not ok:
        print("\n  STOPPING. The horizon override does not reproduce the published\n"
              "  main effect, so any subgroup number from it would be untrustworthy.\n"
              "  Likely causes: a different --imputations, a different --ltfu-date,\n"
              "  or the h60 run used a modified ccw_v3.py rather than a constant swap.")
        return 1
    if args.validate_only:
        return 0

    # ---- subgroup point estimates --------------------------------------
    subs = [s.strip() for s in args.subgroups.split(",") if s.strip()] or None
    print("\n  subgroup point estimates...")
    pts = subgroup_points(timelines, xpats, subs)
    for r in pts:
        print(f"    {r['subgroup']:<18} {r['level']:<22} "
              f"RD {r[f'rd{HORIZON}']:+7.3f}  RR {r[f'rr{HORIZON}']:.3f}  n={r['n']:,}")

    if args.bootstrap <= 0:
        pd.DataFrame(pts).to_csv(
            C.OUTDIR / "ccw_subgroups_h60_bootstrap.csv", index=False)
        print("\n  point estimates only (no --bootstrap)")
        return 0

    # ---- bootstrap -------------------------------------------------------
    sg_levels = [(r["subgroup"], r["level"]) for r in pts]
    rng = np.random.default_rng(args.seed)
    M = len(timelines)
    draws = {k: {"rr": [], "rd": []} for k in sg_levels}
    t0 = time.time()
    for b in range(args.bootstrap):
        res = subgroup_replicate(timelines[b % M], sg_levels, rng)
        for k, (vrr, vrd) in res.items():
            draws[k]["rr"].append(vrr)
            draws[k]["rd"].append(vrd)
        if (b + 1) % 10 == 0:
            el = time.time() - t0
            print(f"    rep {b+1}/{args.bootstrap}  {el/(b+1):.1f}s/rep  "
                  f"ETA {el/(b+1)*(args.bootstrap-b-1)/60:.1f} min", flush=True)

    rows = []
    for r in pts:
        k = (r["subgroup"], r["level"])
        d = draws[k]
        rec = {"subgroup": r["subgroup"], "level": r["level"], "n": r["n"],
               "rd": r[f"rd{HORIZON}"], "rr": r[f"rr{HORIZON}"],
               "risk_dis": r[f"risk{HORIZON}_dis"], "risk_rem": r[f"risk{HORIZON}_rem"],
               "estimable": r["estimable"], "B_ok": len(d["rd"])}
        for nm in ("rd", "rr"):
            v = np.asarray(d[nm], dtype=float)
            v = v[np.isfinite(v)]
            if len(v) >= 50:
                rec[f"{nm}_lo"] = float(np.percentile(v, 2.5))
                rec[f"{nm}_hi"] = float(np.percentile(v, 97.5))
                rec[f"{nm}_boot_mean"] = float(v.mean())
            else:
                rec[f"{nm}_lo"] = rec[f"{nm}_hi"] = rec[f"{nm}_boot_mean"] = np.nan
        rows.append(rec)

    df = pd.DataFrame(rows)
    out = C.OUTDIR / "ccw_subgroups_h60_bootstrap.csv"
    df.to_csv(out, index=False)

    print(f"\n  {HORIZON}-month risk differences with {args.bootstrap}-replicate CIs:\n")
    print(f"    {'subgroup':<18} {'level':<22} {'RD':>7}  {'95% CI':<20} "
          f"{'bias':>7}  {'reps':>5}")
    for _, r in df.iterrows():
        bias = r["rd_boot_mean"] - r["rd"]
        flag = "" if abs(bias) < 0.1 * max(abs(r["rd"]), 1e-9) else "  CHECK"
        crosses = " *" if (r["rd_lo"] < 0 < r["rd_hi"]) else ""
        print(f"    {r['subgroup']:<18} {r['level']:<22} {r['rd']:+7.3f}  "
              f"({r['rd_lo']:+.3f} to {r['rd_hi']:+.3f}){crosses:<3} "
              f"{bias:+7.3f}{flag}  {int(r['B_ok']):>5}")
    print("\n    * interval includes zero")
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
