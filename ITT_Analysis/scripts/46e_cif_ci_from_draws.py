#!/usr/bin/env python3
"""Build the 46d confidence intervals from however many bootstrap replicates exist.

46d writes its summary CSV only after the last replicate. At ~190 s per replicate
a full B=180 run takes ~9.5 hours, so anyone who needs the intervals before it
finishes would otherwise be blocked on a job that is already 70% done and whose
draws are sitting on disk.

46d appends to rolling_cause_cif_draws.csv after every replicate precisely so
this is possible. This reads whatever is there and produces the same output the
script would have produced, with the achieved replicate count recorded in the
file rather than assumed.

Safe to run while 46d is still going: read-only on the draws, and it writes to a
different filename than 46d does.

    python3 46e_cif_ci_from_draws.py
    python3 46e_cif_ci_from_draws.py --min-reps 50
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ITT_RESULTS_DIR                     # noqa: E402

KEYS = ["dmon", "cause", "time_y"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-reps", type=int, default=20,
                    help="refuse to write an interval from fewer than this")
    args = ap.parse_args()

    res = Path(ITT_RESULTS_DIR) / "rolling_cause_cif.csv"
    drw = Path(ITT_RESULTS_DIR) / "rolling_cause_cif_draws.csv"
    for f in (res, drw):
        if not f.exists():
            sys.exit(f"missing {f} -- has 46d run?")

    pt = pd.read_csv(res)
    # a run cut mid-write can leave one truncated final line
    bd = pd.read_csv(drw, on_bad_lines="skip")
    n_reps = bd["rep"].nunique()
    print(f"{len(bd):,} draw rows from {n_reps} replicates")
    if n_reps < args.min_reps:
        sys.exit(f"only {n_reps} replicates -- too few for an interval "
                 f"(--min-reps to override)")

    ci = (bd.groupby(KEYS, dropna=False)
            .agg(rd_lo=("rd", lambda s: s.quantile(.025)),
                 rd_hi=("rd", lambda s: s.quantile(.975)),
                 rr_lo=("rr", lambda s: s.quantile(.025)),
                 rr_hi=("rr", lambda s: s.quantile(.975)),
                 n_reps=("rd", "count"))
            .reset_index())
    out = pt.merge(ci, on=KEYS, how="left")
    out["boot_status"] = f"partial_{n_reps}_of_180"

    dst = Path(ITT_RESULTS_DIR) / "rolling_cause_cif_boot_partial.csv"
    out.to_csv(dst, index=False)                        # write before printing
    print(f"wrote {dst}")

    for hz in sorted(out.time_y.unique()):
        s = out[out.dmon.isna() & (out.time_y == hz)]
        if not len(s):
            continue
        print(f"\n--- overall, {hz:g} y, {n_reps} replicates ---")
        for _, r in s.iterrows():
            flag = ""
            if pd.notna(r.rd_lo) and r.rd_lo < 0 < r.rd_hi:
                flag = "  *includes zero*"
            print(f"  {r['cause']:8s} RD {r['rd']:+7.3f} "
                  f"({r['rd_lo']:+7.3f} to {r['rd_hi']:+7.3f}){flag}")

    bym = out[out.dmon.notna()]
    if len(bym):
        print(f"\n--- by month, 5 y, {n_reps} replicates ---")
        for m in sorted(bym.dmon.unique()):
            s = bym[(bym.dmon == m) & (bym.time_y == 5)]
            parts = []
            for c in ("tb", "nontb", "unclass", "all"):
                r = s[s.cause == c]
                if not len(r):
                    continue
                r = r.iloc[0]
                z = "*" if pd.notna(r.rd_lo) and r.rd_lo < 0 < r.rd_hi else " "
                parts.append(f"{c} {r['rd']:+6.3f}({r['rd_lo']:+6.2f},{r['rd_hi']:+6.2f}){z}")
            print(f"  month {int(m)}: " + "  ".join(parts))
        print("\n  * = interval includes zero, i.e. not reportable as an effect")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
