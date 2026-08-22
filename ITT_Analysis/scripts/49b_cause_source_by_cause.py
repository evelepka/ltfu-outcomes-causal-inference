#!/usr/bin/env python3
"""Where the cause of death came from, split BY CAUSE as well as by arm.

Reviewer 2 asks for the death-certificate versus programmatic split separately
for tuberculosis and non-tuberculosis deaths. The response letter records this as
outstanding, but nothing new needs computing: 49's
rolling_late_death_sources.csv already carries arm x cause x detection x source,
and this collapses it onto the axis the reviewer asked about.

The point of the table is differential attribution. If the death certificate
supplies the cause for a similar share of deaths in both arms, a difference in
attributed cause is unlikely to be an artefact of where the cause came from.

Output: ITT_Analysis/results/cause_source_by_cause.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ITT_RESULTS_DIR                       # noqa: E402

RES = Path(ITT_RESULTS_DIR)
src = pd.read_csv(RES / "rolling_late_death_sources.csv")

rows = []
for (arm, cause), g in src.groupby(["arm", "cause"], sort=False):
    n = int(g.n.sum())
    cert = int(g.loc[g.cause_src == "SIM ICD code", "n"].sum())
    prog = int(g.loc[g.cause_src == "TBweb outcome", "n"].sum())
    rows.append(dict(arm=arm, cause=cause, n_deaths=n,
                     n_death_certificate=cert,
                     pct_death_certificate=round(100 * cert / n, 1) if n else None,
                     n_programmatic=prog,
                     pct_programmatic=round(100 * prog / n, 1) if n else None))
out = pd.DataFrame(rows)

# share of each arm's late-window deaths that each cause accounts for
tot = out.groupby("arm").n_deaths.transform("sum")
out["pct_of_arm_deaths"] = (100 * out.n_deaths / tot).round(1)

dst = RES / "cause_source_by_cause.csv"
out.to_csv(dst, index=False)                              # write before printing
print(f"wrote {dst}\n")
print(out.to_string(index=False))

print("\nThe comparison the reviewer is asking for:")
for cause in ("TB", "non-TB"):
    s = out[out.cause == cause].set_index("arm")
    if {"LTFU", "in care"} <= set(s.index):
        a, b = s.loc["LTFU", "pct_death_certificate"], s.loc["in care", "pct_death_certificate"]
        verdict = "comparable" if abs(a - b) < 5 else "DIFFERENTIAL"
        print(f"  {cause:7s} cause from a death certificate: "
              f"LTFU {a:.1f}% vs in care {b:.1f}%   -> {verdict}")
