"""
Exploratory: 1- and 2-year cumulative incidence of mortality and retreatment
among LTFU patients, by year of treatment start (best_start), 2013-2022.

Time origin = best_start (treatment start).  Use time_d_tx / time_rn_tx —
NOT time_d / time_rn (those are anchored at end_date and would mis-align
the K-year windows). See CLAUDE.md "Time variable trap".

2023 dropped: cohort inclusion requires end_date <= 2023-12-31, so 2023
treatment starts are systematically depleted of long episodes (right-
truncation). Censoring date = 2024-12-31, so all 2013-2022 starters have
>=2y of potential follow-up; the naive proportion equals the Aalen-Johansen
estimate at K=1y and K=2y for retreatment with death as a competing event.

Wilson 95% CIs via scipy.stats.binomtest.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def _find_project_root() -> Path:
    if os.environ.get("TB_ABANDONMENT_ROOT"):
        p = Path(os.environ["TB_ABANDONMENT_ROOT"]).expanduser()
        if p.exists():
            return p
    candidates = [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]
    for c in candidates:
        if c.exists():
            return c
    return Path(__file__).resolve().parents[2]


BASE_DIR = _find_project_root()
COHORT_CSV = BASE_DIR / "ITT_Analysis" / "data" / "itt_cohort.csv"
OUT_DIR = BASE_DIR / "ITT_Analysis" / "results" / "year_trends"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Reading cohort: {COHORT_CSV}")
df = pd.read_csv(
    COHORT_CSV,
    usecols=[
        "sinan_clean", "itt_group", "best_start",
        "time_d_tx", "event_d",
        "time_rn_tx", "event_rn",
    ],
    parse_dates=["best_start"],
)

# Restrict to LTFU
df = df[df["itt_group"] == "Loss to follow-up"].copy()
df["year_start"] = df["best_start"].dt.year

# Drop pre-2013 stragglers and 2023 (right-truncated)
pre = df[df["year_start"] < 2013]
if len(pre):
    print(f"Dropping {len(pre)} pre-2013 records (year_start={sorted(pre['year_start'].unique())})")
df = df[(df["year_start"] >= 2013) & (df["year_start"] <= 2022)].copy()

print(f"N LTFU 2013-2022: {len(df)}")


def wilson_ci(k: int, n: int):
    if n == 0:
        return (np.nan, np.nan)
    return stats.binomtest(k, n).proportion_ci(method="wilson")


# Per-year K-year status
rows = []
for yr, g in df.groupby("year_start"):
    n = len(g)
    # Mortality at K years (origin = best_start)
    k1_mort = int(((g["event_d"] == 1) & (g["time_d_tx"] <= 1.0)).sum())
    k2_mort = int(((g["event_d"] == 1) & (g["time_d_tx"] <= 2.0)).sum())
    # Retreatment at K years — naive proportion (= AJ since no censoring before K
    # for this cohort with end_date <= 2023-12-31 and censoring at 2024-12-31)
    k1_rt = int(((g["event_rn"] == 1) & (g["time_rn_tx"] <= 1.0)).sum())
    k2_rt = int(((g["event_rn"] == 1) & (g["time_rn_tx"] <= 2.0)).sum())

    m1_lo, m1_hi = wilson_ci(k1_mort, n)
    m2_lo, m2_hi = wilson_ci(k2_mort, n)
    r1_lo, r1_hi = wilson_ci(k1_rt, n)
    r2_lo, r2_hi = wilson_ci(k2_rt, n)

    rows.append({
        "year_start": int(yr),
        "n_ltfu": n,
        "mort_1y_pct": 100 * k1_mort / n,
        "mort_1y_lo": 100 * m1_lo,
        "mort_1y_hi": 100 * m1_hi,
        "mort_2y_pct": 100 * k2_mort / n,
        "mort_2y_lo": 100 * m2_lo,
        "mort_2y_hi": 100 * m2_hi,
        "retreat_1y_pct": 100 * k1_rt / n,
        "retreat_1y_lo": 100 * r1_lo,
        "retreat_1y_hi": 100 * r1_hi,
        "retreat_2y_pct": 100 * k2_rt / n,
        "retreat_2y_lo": 100 * r2_lo,
        "retreat_2y_hi": 100 * r2_hi,
    })

tab = pd.DataFrame(rows).sort_values("year_start").reset_index(drop=True)

out_csv = OUT_DIR / "outcomes_ltfu_by_year.csv"
tab.to_csv(out_csv, index=False)
print(f"Wrote {out_csv}")
print(tab.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

# Plot — two panels (mortality, retreatment) with 1y and 2y series
fig, (ax_m, ax_r) = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True)

for ax in (ax_m, ax_r):
    ax.axvspan(2019.5, 2021.5, color="#fde0c5", alpha=0.5, zorder=0,
               label="COVID-19 (2020-2021)")
    ax.grid(axis="y", alpha=0.3)

# Mortality
ax_m.errorbar(
    tab["year_start"], tab["mort_1y_pct"],
    yerr=[tab["mort_1y_pct"] - tab["mort_1y_lo"], tab["mort_1y_hi"] - tab["mort_1y_pct"]],
    fmt="o-", color="#3182bd", capsize=3, lw=1.5, ms=5, label="1-year",
)
ax_m.errorbar(
    tab["year_start"], tab["mort_2y_pct"],
    yerr=[tab["mort_2y_pct"] - tab["mort_2y_lo"], tab["mort_2y_hi"] - tab["mort_2y_pct"]],
    fmt="s--", color="#08519c", capsize=3, lw=1.5, ms=5, label="2-year",
)
ax_m.set_ylabel("All-cause mortality from\ntreatment start (%)")
ax_m.set_title("LTFU outcomes by year of treatment start (origin = best_start)")
ax_m.legend(loc="upper left", fontsize=8, frameon=False)

# Retreatment
ax_r.errorbar(
    tab["year_start"], tab["retreat_1y_pct"],
    yerr=[tab["retreat_1y_pct"] - tab["retreat_1y_lo"], tab["retreat_1y_hi"] - tab["retreat_1y_pct"]],
    fmt="o-", color="#e6550d", capsize=3, lw=1.5, ms=5, label="1-year",
)
ax_r.errorbar(
    tab["year_start"], tab["retreat_2y_pct"],
    yerr=[tab["retreat_2y_pct"] - tab["retreat_2y_lo"], tab["retreat_2y_hi"] - tab["retreat_2y_pct"]],
    fmt="s--", color="#a63603", capsize=3, lw=1.5, ms=5, label="2-year",
)
ax_r.set_ylabel("Retreatment cumulative\nincidence (%, death = competing)")
ax_r.set_xlabel("Year of treatment start (best_start)")
ax_r.set_xticks(tab["year_start"])
ax_r.legend(loc="upper left", fontsize=8, frameon=False)

plt.tight_layout()
out_png = OUT_DIR / "outcomes_ltfu_by_year.png"
plt.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"Wrote {out_png}")
