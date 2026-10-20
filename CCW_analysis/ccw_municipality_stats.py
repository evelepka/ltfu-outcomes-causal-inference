#!/usr/bin/env python3
"""Municipality cluster statistics quoted in Appendix 9 and the response letter.

WHY THIS EXISTS
---------------
The effective number of clusters (7.6), the capital's share of the cohort
(33.5%) and the two intraclass correlations (8.4% for disengagement, 3.4% for
death) were quoted from a docstring, a run-time print and two script
constants, so no results file existed to verify them against. This writes
them to one.

The municipality is the one the analyses use: ccw_random_intercept.municipality(),
the index-episode tx_city rule shared with ccw_v3.build_geo_lookup.

The ICCs are on the latent logistic scale, sigma^2 / (sigma^2 + pi^2/3), from
the random-intercept variances the models actually carry:
  disengagement  SIGMA2_DEFAULT in ccw_random_intercept.py (CCW censoring model)
  death in 5 y   MUNI_SIGMA2 default in ITT_Analysis/scripts/_rolling.R (landmark)
Both are glmer ML estimates; see the docstrings of those two files.

Usage:  python3 CCW_analysis/ccw_municipality_stats.py
Output: results_v3/ccw_municipality_stats.csv, never overwriting an existing file.
"""
from __future__ import annotations
import math, re, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ccw_v3 as C
import ccw_random_intercept as RI

CAPITAL = "SAO PAULO"


def rolling_sigma2() -> float:
    src = (C.ROOT / "ITT_Analysis/scripts/_rolling.R").read_text(encoding="utf-8")
    m = re.search(r'MUNI_SIGMA2\s*<-\s*as\.numeric\(Sys\.getenv\("MUNI_SIGMA2",\s*unset\s*=\s*"([0-9.]+)"\)\)', src)
    if not m:
        raise SystemExit("could not read MUNI_SIGMA2's default from _rolling.R")
    return float(m.group(1))


def main() -> int:
    dest = C.OUTDIR / "ccw_municipality_stats.csv"
    if dest.exists():
        raise SystemExit(f"{dest.name} exists; move it aside to regenerate")
    imp = sorted(C.MI_DIR.glob("imp_*.csv"))[0]
    n = len(pd.read_csv(imp, low_memory=False, usecols=["sinan_clean"]))
    city = pd.Series(RI.municipality(imp, n))
    counts = city.value_counts()
    share = counts / counts.sum()
    if CAPITAL not in counts.index:
        raise SystemExit(f"{CAPITAL!r} not among municipality names; check spelling")
    s_ltfu, s_death = RI.SIGMA2_DEFAULT, rolling_sigma2()
    icc = lambda s: s / (s + math.pi ** 2 / 3)
    rows = [
        ("n_patients", n, "imputed cohort, one row per patient"),
        ("n_municipalities", int((counts.index != "__NA__").sum()), "distinct tx_city, excluding missing"),
        ("n_missing_municipality", int(counts.get("__NA__", 0)), "patients with no tx_city on any episode"),
        ("effective_clusters", round(1 / float(np.sum(share ** 2)), 4), "inverse Herfindahl index of municipality shares"),
        ("capital_share_pct", round(100 * float(share[CAPITAL]), 4), f"share of patients with tx_city {CAPITAL}"),
        ("top5_share_pct", round(100 * float(share.iloc[:5].sum()), 4), "five largest municipalities"),
        ("n_municipalities_lt50", int((counts[counts.index != "__NA__"] < 50).sum()), "municipalities with fewer than 50 patients"),
        ("sigma2_ltfu", s_ltfu, "random-intercept variance, disengagement (CCW censoring model)"),
        ("icc_ltfu_pct", round(100 * icc(s_ltfu), 4), "latent-scale ICC, disengagement"),
        ("sigma2_death", s_death, "random-intercept variance, death within 5 y (landmark)"),
        ("icc_death_pct", round(100 * icc(s_death), 4), "latent-scale ICC, death within 5 y"),
    ]
    out = pd.DataFrame(rows, columns=["quantity", "value", "definition"])
    out["source_imputation"] = imp.name
    out.to_csv(dest, index=False)
    print(out.to_string(index=False))
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
