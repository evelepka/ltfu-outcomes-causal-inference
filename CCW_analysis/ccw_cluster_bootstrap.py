#!/usr/bin/env python3
"""Municipality cluster bootstrap for the primary CCW contrast.

WHY THIS EXISTS
---------------
Reviewer 1 (comment 1, R2 round) asks that non-independence of patients within a
municipality be addressed by GEE, a mixed model or a bootstrap, and discourages
treating it as a limitation. The published interval resamples PATIENTS, which
assumes each patient is an independent draw.

This resamples MUNICIPALITIES. Each draw takes a whole municipality's patients as
a block, so whatever makes patients in one city alike travels with them, and the
spread across replicates reflects between-city variation rather than
between-patient variation alone.

WHAT TO EXPECT
--------------
Sao Paulo State has 640 municipalities but only 7.6 EFFECTIVE clusters
(Herfindahl 1/sum(share^2)); Sao Paulo city alone holds 33.5% of the cohort and
the top five hold 52%. A cluster bootstrap on this structure is a lottery over
whether the capital is drawn none, once or three times, so the interval may be
wide for reasons of design rather than of evidence. That outcome is itself the
answer to the comment, and `docs/dead-ends.md` records the same structure ruling
out cluster-robust sandwich standard errors.

The point estimate is unchanged by construction; only the interval moves.

Usage:
    python3 ccw_cluster_bootstrap.py --bootstrap 1000
Output: results_v3/ccw_cluster_bootstrap_<tag>_h<H>.csv
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ccw_v3 as C


def curves(tl, Xp, T, cause, H):
    out = []
    for arm, t in (("disengage", T), ("remain", None)):
        r = C.add_ipcw(C.expand(C.build_arm(tl, arm, T=t, cause=cause)), Xp)
        w, d, m = r["swt"], r["death"], r["month"]
        num = np.bincount(m, weights=w * d, minlength=H)
        den = np.bincount(m, weights=w, minlength=H)
        hz = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        out.append(1 - np.cumprod(1 - hz))
    return out


def contrast(tl, Xp, T, cause, H):
    ci1, ci0 = curves(tl, Xp, T, cause, H)
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(ci0 > 0, ci1 / ci0, np.nan)
    return 100 * (ci1 - ci0), rr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--imputations", type=int, default=5)
    ap.add_argument("--horizon", type=int, default=60)
    ap.add_argument("--ltfu-date", choices=["shift30", "closure"], default="closure")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    H, T, tag = a.horizon, C.PRIMARY_T, a.ltfu_date
    C.HORIZON_M = H
    C.REPORT_AT = sorted(set(C.REPORT_AT + [H]))

    lookup = C.build_cause_lookup()
    imps = sorted(C.MI_DIR.glob("imp_*.csv"))[:a.imputations]
    raw = pd.read_csv(C.RAW_CSV, low_memory=False,
                      usecols=["sinan_clean", "tx_city"]).drop_duplicates("sinan_clean")

    tls, xps, groups = [], [], []
    for i, p in enumerate(imps):
        tl = C.load_timeline(p, lookup, tag, verbose=(i == 0))
        ids = pd.read_csv(p, low_memory=False, usecols=["sinan_clean"])["sinan_clean"]
        city = ids.to_frame().merge(raw, on="sinan_clean", how="left")["tx_city"]
        assert len(city) == len(tl)
        codes, uniq = pd.factorize(city)
        tl, Xp = C.attach_patterns(tl, C.COVS)
        tls.append(tl); xps.append(Xp)
        # row indices per municipality, so a draw takes the whole block
        order = np.argsort(codes, kind="stable")
        bounds = np.searchsorted(codes[order], np.arange(len(uniq) + 1))
        groups.append([order[bounds[k]:bounds[k + 1]] for k in range(len(uniq))])
        if i == 0:
            sh = np.array([len(g) for g in groups[0]], dtype=float); sh /= sh.sum()
            print(f"  municipalities {len(uniq)}  effective {1/np.sum(sh**2):.1f}  "
                  f"largest {100*sh.max():.1f}%")
    M = len(tls)

    rds, rrs = zip(*[contrast(tl, Xp, T, "all_cause", H) for tl, Xp in zip(tls, xps)])
    rd_pt = np.mean(rds, axis=0)
    rr_pt = np.exp(np.nanmean(np.log(rrs), axis=0))
    print(f"\ncluster bootstrap | {tag} | M={M} | B={a.bootstrap}")
    print(f"  point {H} m: RD {rd_pt[H-1]:+.4f} pp   RR {rr_pt[H-1]:.4f}")

    rng = np.random.default_rng(a.seed)
    dr = np.full((a.bootstrap, H), np.nan)
    br = np.full((a.bootstrap, H), np.nan)
    sizes = np.zeros(a.bootstrap)
    t0 = time.time()
    for b in range(a.bootstrap):
        k = b % M
        tl, Xp, gp = tls[k], xps[k], groups[k]
        pick = rng.integers(0, len(gp), len(gp))
        idx = np.concatenate([gp[j] for j in pick])
        sizes[b] = len(idx)
        d = tl.iloc[idx].reset_index(drop=True)
        try:
            dr[b], br[b] = contrast(d, Xp, T, "all_cause", H)
        except Exception as e:
            print(f"  rep {b+1} failed: {e}")
        if (b + 1) % 50 == 0:
            el = time.time() - t0
            print(f"  rep {b+1}/{a.bootstrap} | {el/(b+1):.1f}s/rep | "
                  f"ETA {el/(b+1)*(a.bootstrap-b-1)/60:.0f} min", flush=True)

    ok = ~np.isnan(dr[:, H - 1])
    rd_lo, rd_hi = np.percentile(dr[ok], [2.5, 97.5], axis=0)
    rr_lo, rr_hi = np.nanpercentile(br[ok], [2.5, 97.5], axis=0)
    out = pd.DataFrame({"month": np.arange(1, H + 1), "rd": rd_pt,
                        "rd_lo": rd_lo, "rd_hi": rd_hi, "rr": rr_pt,
                        "rr_lo": rr_lo, "rr_hi": rr_hi, "B_ok": int(ok.sum()),
                        "unit": "municipality", "ltfu_date": tag, "n_imp": M})
    dest = C.OUTDIR / f"ccw_cluster_bootstrap_{tag}_h{H}.csv"
    out.to_csv(dest, index=False)

    r = out.iloc[H - 1]
    print(f"\n  resampled cohort size: median {np.median(sizes[ok]):,.0f}  "
          f"range {sizes[ok].min():,.0f}-{sizes[ok].max():,.0f}")
    print(f"  CLUSTER  {H} m: RD {r.rd:+.4f} ({r.rd_lo:+.4f}, {r.rd_hi:+.4f})  "
          f"width {r.rd_hi-r.rd_lo:.4f}")
    pb = C.OUTDIR / f"ccw_main_bootstrap_{tag}_h{H}.csv"
    if pb.exists():
        q = pd.read_csv(pb).iloc[H - 1]
        print(f"  PATIENT  {H} m: RD {q.rd:+.4f} ({q.rd_lo:+.4f}, {q.rd_hi:+.4f})  "
              f"width {q.rd_hi-q.rd_lo:.4f}")
        print(f"  width ratio cluster/patient: {(r.rd_hi-r.rd_lo)/(q.rd_hi-q.rd_lo):.2f}x")
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
