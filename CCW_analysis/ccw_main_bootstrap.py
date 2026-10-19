#!/usr/bin/env python3
"""Main CCW contrast with confidence intervals, month by month.

WHY THIS EXISTS
---------------
`ccw_bootstrap.py` produced the published `ccw_bootstrap_h60.csv` but reads its
input from a scratch directory belonging to a long-gone session, so it cannot be
re-run: the paper's primary interval had no runnable source. This replaces it,
reading the imputations directly and reusing ccw_v3's own machinery so there is
one implementation of the estimator.

It also reports the WHOLE monthly curve rather than three horizons, because
Reviewer 1 (comment 2) and Reviewer 3 both ask about the shape of Figure 4 in
early follow-up, which three points cannot show.

WHAT CHANGED IN THE ESTIMATOR (ADR-0006)
----------------------------------------
  * clones censored inside month 0 are censored, not deleted (`build_arm`)
  * the disengagement date is case closure, not closure minus 30 days, so no
    clone is censored before the 30-day LTFU definition can be met

Both corrections raise the comparator arm's early risk back toward the exposed
arm's and shrink the early negative risk difference by about two thirds; the
residual is healthy-leaver selection that baseline covariates cannot reach.

VALIDATION GATE
---------------
Reproduces the 60-month point estimate in ccw_v3_<tag>_main.csv before
bootstrapping, so a horizon or convention mismatch stops the run rather than
producing plausible-looking intervals.

Usage:
    python3 ccw_main_bootstrap.py --bootstrap 1000
    python3 ccw_main_bootstrap.py --bootstrap 1000 --ltfu-date shift30
Output: results_v3/ccw_main_bootstrap_<tag>_h<H>.csv
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


def curves(tl, Xp, T, cause, H):
    """Weighted CIF per month for both arms."""
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
    ap.add_argument("--cause", default="all_cause")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    H, T, tag = a.horizon, C.PRIMARY_T, a.ltfu_date
    C.HORIZON_M = H
    C.REPORT_AT = sorted(set(C.REPORT_AT + [H]))

    lookup = C.build_cause_lookup()
    imps = sorted(C.MI_DIR.glob("imp_*.csv"))[:a.imputations]
    assert imps, f"no imputations in {C.MI_DIR}"
    pairs = [C.attach_patterns(C.load_timeline(p, lookup, tag, verbose=(i == 0)), C.COVS)
             for i, p in enumerate(imps)]
    tls, xps = [t for t, _ in pairs], [x for _, x in pairs]
    M = len(tls)

    # ---- point estimate, pooled across imputations -----------------------
    rds, rrs = zip(*[contrast(tl, Xp, T, a.cause, H) for tl, Xp in zip(tls, xps)])
    rd_pt = np.mean(rds, axis=0)
    rr_pt = np.exp(np.nanmean(np.log(rrs), axis=0))

    # ---- gate against ccw_v3's own output --------------------------------
    ref = C.OUTDIR / f"ccw_v3_{tag}_main.csv"
    if ref.exists():
        d = pd.read_csv(ref)
        row = d[d["T"].eq(T) & d["cause"].eq(a.cause)]
        if len(row):
            drd = abs(rd_pt[H - 1] - row[f"rd{H}"].iloc[0])
            drr = abs(rr_pt[H - 1] - row[f"rr{H}"].iloc[0])
            print(f"\n[gate] {ref.name}: RD diff {drd:.2e}  RR diff {drr:.2e}")
            if drd > 1e-9 or drr > 1e-9:
                print("[gate] NAO BATE -- parando antes de bootstrapar.")
                return 1
            print("[gate] bate.")

    print(f"\nCCW main bootstrap | {tag} | T={T} | {a.cause} | H={H} | "
          f"M={M} | B={a.bootstrap}")
    print(f"  ponto {H} m: RD {rd_pt[H-1]:+.4f} pp   RR {rr_pt[H-1]:.4f}")

    rng = np.random.default_rng(a.seed)
    dr = np.full((a.bootstrap, H), np.nan)
    br = np.full((a.bootstrap, H), np.nan)
    t0 = time.time()
    for b in range(a.bootstrap):
        tl, Xp = tls[b % M], xps[b % M]
        d = tl.iloc[rng.integers(0, len(tl), len(tl))].reset_index(drop=True)
        try:
            dr[b], br[b] = contrast(d, Xp, T, a.cause, H)
        except Exception as e:
            print(f"  rep {b+1} falhou: {e}")
        if (b + 1) % 50 == 0:
            el = time.time() - t0
            print(f"  rep {b+1}/{a.bootstrap} | {el/(b+1):.1f}s/rep | "
                  f"ETA {el/(b+1)*(a.bootstrap-b-1)/60:.0f} min", flush=True)

    ok = ~np.isnan(dr[:, H - 1])
    rd_lo, rd_hi = np.percentile(dr[ok], [2.5, 97.5], axis=0)
    rr_lo, rr_hi = np.nanpercentile(br[ok], [2.5, 97.5], axis=0)

    out = pd.DataFrame({"month": np.arange(1, H + 1),
                        "rd": rd_pt, "rd_lo": rd_lo, "rd_hi": rd_hi,
                        "rr": rr_pt, "rr_lo": rr_lo, "rr_hi": rr_hi,
                        "crosses_zero": (rd_lo <= 0) & (0 <= rd_hi),
                        "B_ok": int(ok.sum()), "ltfu_date": tag,
                        "T": T, "cause": a.cause, "n_imp": M})
    dest = C.OUTDIR / f"ccw_main_bootstrap_{tag}_h{H}.csv"
    out.to_csv(dest, index=False)

    print(f"\n   mes |    RD (pp)  |       IC 95%        | cruza zero")
    for mo in [m for m in (1, 2, 3, 4, 5, 6, 7, 8, 12, 24, H) if m <= H]:
        r = out.iloc[mo - 1]
        print("   %3d |    %+7.4f  | %+7.4f , %+7.4f  |   %s"
              % (mo, r.rd, r.rd_lo, r.rd_hi, "sim" if r.crosses_zero else "nao"))
    r = out.iloc[H - 1]
    print(f"\n  TITULAR {H} m: RD {r.rd:+.4f} pp ({r.rd_lo:+.4f}, {r.rd_hi:+.4f})  "
          f"RR {r.rr:.4f} ({r.rr_lo:.4f}, {r.rr_hi:.4f})  B_ok={int(r.B_ok)}")
    print(f"escrito {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
