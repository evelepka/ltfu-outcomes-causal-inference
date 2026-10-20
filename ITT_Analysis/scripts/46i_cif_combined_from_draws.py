#!/usr/bin/env python3
"""Combined cause quantities from the 46d bootstrap draws.

The Results report two quantities that are functions of several causes at once:
  not_tb    the excess NOT attributed to tuberculosis = nontb + unclass
  tb_share  the share of the all-cause excess attributed to tuberculosis, %
Adding the per-cause interval limits would be wrong; the interval has to come
from the same replicate's draws, which 46d keeps (one row per rep x cause).

Point estimates are taken from 46d's own point-estimate file, not from the
draws, so the reported number is the pooled estimate and only the interval is
bootstrapped -- the convention 46d and 46e use.

    python3 46i_cif_combined_from_draws.py --suffix _NEWMI_PRIMARY_underlying_certainty

Writes rolling_cause_cif_combined<suffix>.csv; never overwrites.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ITT_RESULTS_DIR                     # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="", help="OUT_SUFFIX of the 46d run")
    ap.add_argument("--min-reps", type=int, default=100)
    a = ap.parse_args()
    res = Path(ITT_RESULTS_DIR)
    pt = pd.read_csv(res / f"rolling_cause_cif{a.suffix}.csv")
    dr = pd.read_csv(res / f"rolling_cause_cif_draws{a.suffix}.csv")
    dst = res / f"rolling_cause_cif_combined{a.suffix}.csv"
    if dst.exists():
        raise SystemExit(f"{dst.name} exists; move it aside to regenerate")
    pt, dr = pt[pt.dmon.isna()], dr[dr.dmon.isna()]
    rows = []
    for t in sorted(pt.time_y.unique()):
        p = pt[pt.time_y == t].set_index("cause").rd
        w = dr[dr.time_y == t].pivot(index="rep", columns="cause", values="rd").dropna()
        if len(w) < a.min_reps:
            raise SystemExit(f"only {len(w)} complete replicates at {t} y")
        for name, point, draws in (
                ("not_tb", p["nontb"] + p["unclass"], w["nontb"] + w["unclass"]),
                ("tb_share", 100 * p["tb"] / p["all"], 100 * w["tb"] / w["all"])):
            rows.append(dict(quantity=name, time_y=t, value=point,
                             lo=draws.quantile(0.025), hi=draws.quantile(0.975),
                             n_reps=len(w)))
    out = pd.DataFrame(rows)
    out.to_csv(dst, index=False)
    print(out.round(3).to_string(index=False))
    print(f"\n  wrote {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
