#!/usr/bin/env python3
"""Figure 5 — what the design choice does to the answer.

THE MESSAGE
-----------
In one dataset, holding cohort, covariates and horizon fixed and changing ONLY the
design, the LTFU-mortality association runs from 0.35 (apparently protective) to
4.97 (nearly fivefold harmful). Published estimates sit inside that range. So the
spread in the literature is not a spread in the biology, and large published
estimates should not be read as the expected benefit of retention interventions.

WHY IT IS A FIGURE
------------------
It answers Reviewer 2 major comment 2 (situate our effect sizes against prior
work) and Reviewer 1 pivotal comment 1 through the "clinical research agendas"
clause. It adds NO new analysis -- every number already exists in
literature_design_reproduction.csv and rolling_landmark.csv -- so it does not
offend Reviewer 1 pivotal comment 2 on excessive scope.

TWO RULES THIS FIGURE MUST OBEY
-------------------------------
1. Published estimates and our reproductions are kept in SEPARATE blocks. They are
   not interchangeable: different populations, covariates and outcome definitions.
   The claim is about the reproductions, where only the design differs.
2. SMRs against the general population (Kolappan 2006, Romanowski 2019) are
   deliberately EXCLUDED. They are not LTFU-versus-comparator contrasts and would
   not belong on this axis.

Also on purpose: the primary reported estimand is a 60-month RISK DIFFERENCE
(+2.22 pp), which is not a hazard ratio and is therefore NOT plotted here. A
footnote says so, so nobody reads 2.42 as the headline.

Usage:  python3 ITT_Analysis/scripts/71_make_fig5_design_forest.py
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
OUTDIR = ROOT / "Plos Medicine" / "R1" / "Figures" / "Draft"
RESULTS = ROOT / "ITT_Analysis" / "results"
STEM = "LTFU_design_forest"

BLUE, VERM, GREY, INK = "#0072B2", "#D55E00", "#8A8A8A", "#222222"
AMBER_INK, PURP = "#8A6100", "#7B4FA0"
FS_T, FS_L, FS_N = 8.6, 7.2, 6.5

# --- published estimates: LTFU/default vs a comparator, hazard ratios only ----
PUB = [
    ("García-García 2002, Mexico", 8.9, 3.3, 24.4,
     "TB-specific death, 34 events; outcome definition itself depends on exposure"),
    ("Nájera-Ortiz 2012, Mexico", 5.74, 3.59, 9.18,
     "exposure = treatment <6 months; comparator completed therapy"),
    ("Cunha 2017, Brazil", 3.65, 2.28, 5.83,
     "people living with HIV only; 317 patients"),
]

# --- our data, only the design changes ---------------------------------------
rep = pd.read_csv(RESULTS / "literature_design_reproduction.csv")
r = {row["design"]: row for _, row in rep.iterrows()}
OURS = [
    ("A", "Day-0 origin, cured comparator",
     "the comparator must survive AND complete therapy, so it cannot die for ~6\n"
     "months: immortal time lands on the COMPARATOR and adds to healthy-survivor\n"
     "selection. Both inflate."),
    ("C", "Month-6 landmark, cured comparator",
     "aligning the origins removes the immortal time, but the comparator is still\n"
     "selected for good prognosis."),
    ("B", "Day-0 origin, all-non-LTFU comparator",
     "all 12,336 on-treatment deaths are classified NOT-LTFU, because you cannot be\n"
     "lost to follow-up if you die first. The comparator absorbs the front-loaded\n"
     "acute mortality the exposed arm is definitionally incapable of having."),
    ("D", "Declaration origin, in-care comparator",
     "comparator restricted to patients still in care on that same date. OUR DESIGN."),
]

rows, labels, notes, colors = [], [], [], []
for lab, est, lo, hi, note in PUB:
    rows.append((est, lo, hi)); labels.append(lab); notes.append(note); colors.append(GREY)
for code, lab, note in OURS:
    row = r[code]
    rows.append((row["aHR"], row["CI_L"], row["CI_H"]))
    labels.append(lab)
    notes.append(note)
    colors.append(VERM if code == "D" else BLUE)

n_pub = len(PUB)
y = np.arange(len(rows))[::-1].astype(float)
y[:n_pub] += 0.85          # gap between the two blocks

from matplotlib.transforms import blended_transform_factory as blend

fig = plt.figure(figsize=(7.6, 4.4))
# labels live in a left margin and values in a right margin, OUTSIDE the data
# area; drawing them inside the axes made them collide with the intervals.
ax = fig.add_axes([0.295, 0.215, 0.43, 0.70])
tr = blend(ax.transAxes, ax.transData)     # x in axes fraction, y in data units
ax.set_xscale("log")
ax.set_xlim(0.22, 28)
ax.axvline(1.0, color=INK, lw=0.9, ls=(0, (3, 3)), zorder=1)

for i, ((est, lo, hi), col) in enumerate(zip(rows, colors)):
    prim = col == VERM
    ax.plot([lo, hi], [y[i], y[i]], color=col, lw=1.5, solid_capstyle="butt", zorder=3)
    for b in (lo, hi):
        ax.plot([b, b], [y[i] - 0.12, y[i] + 0.12], color=col, lw=1.3, zorder=3)
    ax.plot([est], [y[i]], "o", color=col, ms=7.0 if prim else 5.2, zorder=4,
            mec="white", mew=0.8)
    ax.text(-0.035, y[i], labels[i], transform=tr, fontsize=FS_L, color=INK,
            va="center", ha="right", weight="bold" if prim else "normal")
    ax.text(1.035, y[i], f"{est:.2f} ({lo:.2f}–{hi:.2f})", transform=tr,
            fontsize=FS_N, color=col, va="center", ha="left",
            weight="bold" if prim else "normal")

ax.text(-0.035, y[0] + 0.80, "PUBLISHED ESTIMATES", transform=tr,
        fontsize=FS_N, color=GREY, weight="bold", va="bottom", ha="right")
ax.text(1.035, y[0] + 0.80, "different cohorts, covariates and outcomes",
        transform=tr, fontsize=FS_N - 0.5, color=GREY, va="bottom", ha="left")
ax.text(-0.035, y[n_pub] + 0.62,
        "SAME DATA — ONLY THE DESIGN CHANGES", transform=tr,
        fontsize=FS_N, color=INK, weight="bold", va="bottom", ha="right")
ax.text(1.035, y[n_pub] + 0.62, "n = 171,048 throughout", transform=tr,
        fontsize=FS_N - 0.5, color=INK, va="bottom", ha="left")
ax.plot([-0.40, 1.42], [y[n_pub] + 0.46] * 2, transform=tr, color=GREY,
        lw=0.6, clip_on=False)

ax.set_yticks([])
ax.set_ylim(-0.75, y[0] + 1.5)
ax.set_xlabel("Hazard ratio for mortality after loss to follow-up (log scale)",
              fontsize=FS_L, labelpad=3)
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8, 16])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8", "16"], fontsize=FS_L)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(GREY)
ax.tick_params(axis="y", length=0)
ax.text(0.45, -0.55, "← appears protective", fontsize=FS_N - 0.3, color=GREY,
        va="center", ha="center")
ax.text(4.5, -0.55, "appears harmful →", fontsize=FS_N - 0.3, color=GREY,
        va="center", ha="center")

fig.text(0.012, 0.055,
         "Hazard ratios only. Standardised mortality ratios against the general population "
         "(Kolappan 2006, Romanowski 2019) are excluded: they are not\n"
         "loss-to-follow-up-versus-comparator contrasts. The primary reported estimand is a "
         "60-month risk difference (+2.22 percentage points,\n"
         "1.88–2.52), which is not a hazard ratio and is not plotted here.",
         fontsize=FS_N - 0.6, color=GREY, va="bottom", ha="left", linespacing=1.5)

OUTDIR.mkdir(parents=True, exist_ok=True)
for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    p = OUTDIR / f"{STEM}.{ext}"
    fig.savefig(p, facecolor="white", **kw)
    print(f"wrote {p.relative_to(ROOT)}")
plt.close(fig)

out = pd.DataFrame({"row": labels, "hr": [x[0] for x in rows],
                    "ci_lo": [x[1] for x in rows], "ci_hi": [x[2] for x in rows],
                    "block": ["published"] * n_pub + ["reproduction"] * len(OURS),
                    "mechanism": [n.replace("\n", " ") for n in notes]})
out.to_csv(RESULTS / "fig5_design_forest_values.csv", index=False)
print(out[["row", "hr", "ci_lo", "ci_hi", "block"]].round(2).to_string(index=False))
print("\nMechanism notes are in the CSV for the caption; they do not fit on the plot.")
