#!/usr/bin/env python3
"""CCW subgroup risk differences with a choice of bootstrap resampling unit.

WHY THIS EXISTS
---------------
`ccw_subgroup_bootstrap.py` and `ccw_subgroup_random_intercept.py` both resample
PATIENTS. Reviewer 1 (comment 1, R2 round) is about dependence within
municipality, and Run A showed the municipality random intercept does not widen
the patient-level interval (0.68 -> 0.71 at 60 months), so the resampling unit
has to carry that dependence instead. This runs the same 15 strata under the
certainty-unit cluster design.

A FAST PATH WAS TRIED AND REJECTED
----------------------------------
Precomputing each stratum's design matrix once and reusing it across replicates
looks safe -- `attach_patterns` says a resample changes pattern counts, not the
pattern universe -- but it is NOT equivalent here: `_fit_cells` fits a RIDGE-
penalised logistic on an over-parameterised one-hot encoding, so adding or
removing columns changes the penalised solution. Measured discrepancy was
2.2e-02 pp on the subgroup RDs, the same order as the effects being compared.
This file therefore calls `attach_patterns` per stratum per replicate, exactly
as the published script does; only the draw differs. `--verify` re-checks that
against the original implementation on a shared seed.

RESAMPLING UNITS
----------------
  patient    each patient an independent draw (the published convention).
  cluster    municipalities resampled as whole blocks.
  certainty  Sao Paulo city in every replicate with its patients resampled
             within it; the other 639 resampled as whole blocks. Removes the
             capital lottery that makes a plain cluster bootstrap wide and
             lopsided, while keeping between-municipality variation.

Usage:
    python3 ccw_subgroup_certainty.py --verify 3 --imputations 2
    python3 ccw_subgroup_certainty.py --boot-unit certainty --bootstrap 300
Output: results_v3/ccw_subgroups_<model>_<unit>_h60_bootstrap.csv
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

SCRIPT_VERSION = "2026-10-17b"
MIN_N, MIN_DEATHS = 500, 20


def strata_of(tl):
    out = []
    for sg in C.SUBGROUPS:
        for lvl in sorted(set(tl[sg].tolist())):
            if lvl in ("nan", "", "__missing__"):
                continue
            if (tl[sg].to_numpy() == lvl).sum() >= MIN_N:
                out.append((sg, lvl))
    return out


def replicate_from_index(tl, sg_levels, idx, H):
    """Byte-for-byte `SB.subgroup_replicate`, except the resample is given.

    Keeping `attach_patterns` inside the loop is deliberate -- see the module
    docstring on the rejected fast path.
    """
    d = tl.iloc[idx].reset_index(drop=True)
    res = {}
    for sg, lvl in sg_levels:
        covs_sg = [c for c in C.COVS if c != sg]
        s = d[d[sg].to_numpy() == lvl]
        if len(s) < MIN_N or int(s["__any_death__"].sum()) < MIN_DEATHS:
            continue
        try:
            s, Xs = C.attach_patterns(s.reset_index(drop=True), covs_sg)
            r = C.one_contrast(s, C.PRIMARY_T, "all_cause", Xs,
                               ref=C.reference(s, "all_cause", Xs))
        except Exception:                                  # noqa: BLE001
            continue
        res[(sg, lvl)] = (r[f"rr{H}"], r[f"rd{H}"])
    return res


def draw(rng, unit, n, gp, cap):
    if unit == "patient":
        return rng.integers(0, n, n)
    if unit == "cluster":
        pick = rng.integers(0, len(gp), len(gp))
        return np.concatenate([gp[j] for j in pick])
    others = [g for j, g in enumerate(gp) if j != cap]
    pick = rng.integers(0, len(others), len(others))
    big = gp[cap]
    return np.concatenate([big[rng.integers(0, len(big), len(big))]]
                          + [others[j] for j in pick])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot-unit", choices=["patient", "cluster", "certainty"],
                    default="certainty")
    ap.add_argument("--bootstrap", type=int, default=300)
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--horizon", type=int, default=60)
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"], default="closure")
    ap.add_argument("--sigma2", type=float, default=RI.SIGMA2_DEFAULT)
    ap.add_argument("--no-random-intercept", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--verify", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    H = a.horizon
    SB.HORIZON = H
    SB.set_horizon(H)
    sigma = None if a.no_random_intercept else float(np.sqrt(a.sigma2))
    model = "published" if sigma is None else "re"
    if sigma is not None:
        C.attach_patterns = lambda tl, covs: RI.patterns(tl, sigma, covs=covs)

    started = datetime.now(timezone.utc)
    print(f"ccw_subgroup_certainty.py {SCRIPT_VERSION} | "
          f"start {started.isoformat(timespec='seconds')}")
    print(f"  model={model} unit={a.boot_unit} B={a.bootstrap} M={a.imputations} "
          f"horizon={H} seed={a.seed}")

    lookup = C.build_cause_lookup()
    imps = sorted(C.MI_DIR.glob("imp_*.csv"))[:a.imputations]
    tls, groups, caps = [], [], []
    levels = None
    for i, p in enumerate(imps):
        tl = C.load_timeline(p, lookup, a.ltfu_date, verbose=(i == 0))
        tl["muni"] = RI.municipality(p, len(tl))
        if sigma is not None:
            tl, _ = RI.patterns(tl, sigma)
        else:
            tl, _ = C.attach_patterns(tl, C.COVS)
        if levels is None:
            levels = strata_of(tl)
            print(f"  strata: {len(levels)}")
        codes, uniq = pd.factorize(tl["muni"])
        gp = RI.blocks(codes, len(uniq))
        tls.append(tl); groups.append(gp)
        caps.append(int(np.argmax([len(g) for g in gp])))
    M = len(tls)

    if a.verify:
        print(f"\n  verifying {a.verify} replicates against "
              f"SB.subgroup_replicate ...")
        worst = 0.0
        for b in range(a.verify):
            k = b % M
            r1 = SB.subgroup_replicate(tls[k], levels,
                                       np.random.default_rng(1000 + b))
            idx = np.random.default_rng(1000 + b).integers(
                0, len(tls[k]), len(tls[k]))
            r2 = replicate_from_index(tls[k], levels, idx, H)
            for key in set(r1) & set(r2):
                worst = max(worst, abs(r1[key][1] - r2[key][1]))
            diff = set(r1) ^ set(r2)
            if diff:
                print(f"    rep {b}: strata differ: {diff}")
        print(f"  largest RD discrepancy: {worst:.3e}")
        print("  OK" if worst < 1e-9 else "  NOT EQUIVALENT -- do not use")
        return 0 if worst < 1e-9 else 1

    name = a.out or f"ccw_subgroups_{model}_{a.boot_unit}_h{H}_bootstrap.csv"
    dest = C.OUTDIR / name
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")

    print("\n  point estimates...")
    pts = SB.subgroup_points(tls, None)
    for r in pts:
        print(f"    {r['subgroup']:<18} {r['level']:<22} "
              f"RD {r[f'rd{H}']:+7.3f}  RR {r[f'rr{H}']:.3f}  n={r['n']:,}")
    sg_levels = [(r["subgroup"], r["level"]) for r in pts]

    rng = np.random.default_rng(a.seed)
    draws = {k: {"rr": [], "rd": []} for k in sg_levels}
    t0 = time.time()
    for b in range(a.bootstrap):
        k = b % M
        idx = draw(rng, a.boot_unit, len(tls[k]), groups[k], caps[k])
        for key, (vrr, vrd) in replicate_from_index(
                tls[k], sg_levels, idx, H).items():
            draws[key]["rr"].append(vrr); draws[key]["rd"].append(vrd)
        if (b + 1) % 10 == 0:
            el = time.time() - t0
            print(f"    rep {b+1}/{a.bootstrap}  {el/(b+1):.1f}s/rep  "
                  f"ETA {el/(b+1)*(a.bootstrap-b-1)/60:.1f} min", flush=True)
    runtime = time.time() - t0

    rows_out = []
    for r in pts:
        key = (r["subgroup"], r["level"])
        rdv = np.array(draws[key]["rd"], float)
        rrv = np.array(draws[key]["rr"], float)
        lo, hi = (np.percentile(rdv, [2.5, 97.5]) if len(rdv) else (np.nan,)*2)
        rlo, rhi = (np.nanpercentile(rrv, [2.5, 97.5]) if len(rrv) else (np.nan,)*2)
        rows_out.append({
            "subgroup": key[0], "level": key[1], "n": r["n"],
            "rd": r[f"rd{H}"], "rr": r[f"rr{H}"],
            "risk_dis": r[f"risk{H}_dis"], "risk_rem": r[f"risk{H}_rem"],
            "estimable": r["estimable"], "B_ok": int(len(rdv)),
            "rd_lo": lo, "rd_hi": hi,
            "rd_boot_mean": float(np.mean(rdv)) if len(rdv) else np.nan,
            "rr_lo": rlo, "rr_hi": rhi,
            "rr_boot_mean": float(np.mean(rrv)) if len(rrv) else np.nan,
            "model": ("muni_random_intercept" if sigma else "published"),
            "sigma2": (a.sigma2 if sigma else np.nan), "boot_unit": a.boot_unit,
            "script_version": SCRIPT_VERSION,
            "run_started_utc": started.isoformat(timespec="seconds"),
            "runtime_s": round(runtime, 1)})
    pd.DataFrame(rows_out).to_csv(dest, index=False)
    print(f"\n  wrote {dest}   runtime {runtime/60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
