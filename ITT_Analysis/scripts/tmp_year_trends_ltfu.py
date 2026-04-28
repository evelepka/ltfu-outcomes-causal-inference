"""
Exploratory: LTFU rate by year of treatment start (best_start), 2013-2023.

Year axis = year(best_start). LTFU is defined at index-episode end_date from
case_outcome and is therefore fully observed for every cohort member — no
follow-up window is required for this panel.

Caveat (2023 right-truncation): cohort inclusion requires end_date <=
2023-12-31. Non-LTFU (treatment-completion) episodes typically last ~6 months,
so late-2023 starters whose treatment ends in 2024 are excluded. LTFU
episodes end early and are largely retained. This biases the 2023 LTFU% upward
and the 2022 LTFU% modestly upward. Reported here, not corrected.
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
    usecols=["sinan_clean", "itt_group", "best_start", "end_date"],
    parse_dates=["best_start", "end_date"],
)

df["year_start"] = df["best_start"].dt.year
df["is_ltfu"] = (df["itt_group"] == "Loss to follow-up").astype(int)

# Drop the two pre-2013 records (treatment started long before cohort window;
# included only because their index episode end_date fell inside 2013-2023).
pre = df[df["year_start"] < 2013]
if len(pre):
    print(f"Dropping {len(pre)} pre-2013 records (year_start={sorted(pre['year_start'].unique())})")
df = df[df["year_start"] >= 2013].copy()

# Per-year LTFU rate with Wilson 95% CI
rows = []
for yr, g in df.groupby("year_start"):
    n = len(g)
    n_ltfu = int(g["is_ltfu"].sum())
    p = n_ltfu / n
    lo, hi = stats.binomtest(n_ltfu, n).proportion_ci(method="wilson")
    rows.append({
        "year_start": int(yr),
        "n_total": n,
        "n_ltfu": n_ltfu,
        "n_non_ltfu": n - n_ltfu,
        "ltfu_pct": 100 * p,
        "ltfu_ci_lo": 100 * lo,
        "ltfu_ci_hi": 100 * hi,
    })
tab = pd.DataFrame(rows).sort_values("year_start").reset_index(drop=True)

out_csv = OUT_DIR / "ltfu_by_year.csv"
tab.to_csv(out_csv, index=False)
print(f"Wrote {out_csv}")
print(tab.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

# Plot
fig, (ax_n, ax_p) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True,
                                  gridspec_kw={"height_ratios": [1, 1.6]})

# COVID shading (Brazil first case Feb 2020; major impact 2020-2021)
for ax in (ax_n, ax_p):
    ax.axvspan(2019.5, 2021.5, color="#fde0c5", alpha=0.5, zorder=0,
               label="COVID-19 (2020-2021)")

# Top: stacked bars of N
ax_n.bar(tab["year_start"], tab["n_non_ltfu"], color="#9ecae1", label="Non-LTFU")
ax_n.bar(tab["year_start"], tab["n_ltfu"], bottom=tab["n_non_ltfu"],
         color="#d94801", label="LTFU")
ax_n.set_ylabel("Cohort size (N)")
ax_n.legend(loc="upper left", fontsize=8, frameon=False)
ax_n.set_title("Index-episode cohort, by year of treatment start")

# Bottom: LTFU% with CI
ax_p.errorbar(
    tab["year_start"], tab["ltfu_pct"],
    yerr=[tab["ltfu_pct"] - tab["ltfu_ci_lo"], tab["ltfu_ci_hi"] - tab["ltfu_pct"]],
    fmt="o-", color="#a63603", capsize=3, lw=1.5, ms=5,
)
ax_p.set_ylabel("LTFU at index episode end (%)")
ax_p.set_xlabel("Year of treatment start (best_start)")
ax_p.set_xticks(tab["year_start"])
ax_p.grid(axis="y", alpha=0.3)

# Annotate 2023 with caveat marker
last_yr = int(tab["year_start"].max())
last_pct = float(tab.loc[tab["year_start"] == last_yr, "ltfu_pct"].iloc[0])
ax_p.annotate(
    "2023 right-truncated\n(Non-LTFU undercount)",
    xy=(last_yr, last_pct),
    xytext=(last_yr - 2.5, last_pct + 4),
    fontsize=8, ha="left",
    arrowprops=dict(arrowstyle="-", color="grey", lw=0.8),
)

plt.tight_layout()
out_png = OUT_DIR / "ltfu_by_year.png"
plt.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"Wrote {out_png}")
