#!/usr/bin/env python3
"""Figure S7 - annual LTFU incidence (% of cohort) by year of treatment start.

Re-derived on the CURRENT primary cohort (itt_cohort.csv, N=171,069) so the
year totals match the manuscript cohort. 2023 is shown as an open marker
because the index treatment episode is incompletely ascertained for the final
year (administrative censoring at the end of the study window), which inflates
the apparent LTFU proportion.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper"
COHORT = os.path.join(BASE, "ITT_Analysis/data/itt_cohort.csv")
OUT_CSV = os.path.join(BASE, "ITT_Analysis/results/year_trends/ltfu_by_year.csv")
OUT = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S7_ltfu_incidence.png"

d = pd.read_csv(COHORT, low_memory=False,
                usecols=["itt_group", "best_start"])
d["year"] = pd.to_datetime(d["best_start"], errors="coerce").dt.year
d = d[d["year"].between(2013, 2023)]
d["ltfu"] = (d["itt_group"] == "Loss to follow-up").astype(int)

g = d.groupby("year")["ltfu"].agg(["sum", "count"]).reset_index()
g.columns = ["year", "n_ltfu", "n_total"]
g["pct"] = 100 * g["n_ltfu"] / g["n_total"]
# Wilson 95% CI
z = 1.959964
p = g["n_ltfu"] / g["n_total"]; n = g["n_total"]
denom = 1 + z**2 / n
centre = (p + z**2 / (2*n)) / denom
half = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
g["lo"] = 100*(centre - half); g["hi"] = 100*(centre + half)
g.to_csv(OUT_CSV, index=False)
print(g.to_string(index=False))

fig, ax = plt.subplots(figsize=(9.5, 5.2))
# pandemic shading
ax.axvspan(2020, 2022, color="0.88", zorder=0, label="COVID-19 pandemic era")

complete = g[g.year <= 2022]
partial = g[g.year == 2023]

ax.errorbar(complete.year, complete.pct,
            yerr=[complete.pct - complete.lo, complete.hi - complete.pct],
            fmt="o-", color="#c0392b", lw=2, markersize=6, capsize=3,
            markeredgecolor="white", markeredgewidth=0.6)
# connector to 2023 + open marker
ax.plot([2022, 2023], [complete.pct.iloc[-1], partial.pct.iloc[0]],
        color="#c0392b", lw=2, ls=":")
ax.errorbar(partial.year, partial.pct,
            yerr=[partial.pct - partial.lo, partial.hi - partial.pct],
            fmt="o", mfc="white", mec="#c0392b", color="#c0392b",
            markersize=8, capsize=3, lw=2)
ax.annotate("2023: index episode\nincompletely ascertained",
            xy=(2023, partial.pct.iloc[0]),
            xytext=(2020.4, partial.pct.iloc[0] - 0.5),
            fontsize=8.5, color="dimgrey", ha="left", va="top")

ax.set_xticks(range(2013, 2024))
ax.set_xlabel("Year of TB notification", fontsize=11)
ax.set_ylabel("Annual LTFU incidence (% of new TB cohort)", fontsize=11)
ax.set_ylim(0, max(g.hi) * 1.10)
ax.legend(frameon=False, fontsize=10, loc="upper left")
ax.grid(axis="y", color="0.9", lw=0.7)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.1)
print("wrote", OUT)
