#!/usr/bin/env python3
"""
73_make_fig4_revised.py — Figure 4 for the R1 revision.

LAYOUT (owner decision, 2026-08-20)
-----------------------------------
  A  Risk difference by month of disengagement   — ROLLING LANDMARK
  B  Risk ratio by subgroup                      — CCW
  C  Risk difference by subgroup                 — CCW

One design per panel. The previous Figure 4 put landmark hazard ratios and CCW
risk differences side by side inside one panel, which invites the reader to
compare two designs on two clocks — the risk the handoff flags as the main one
in this package. Panels B and C are both CCW, so the relative/absolute inversion
that Reviewer 1 asked about is visible without leaving a single design.

Panel A of the old figure (crude time-varying HR) is dropped.

WHAT IS NOT HERE
----------------
The by-month adjusted hazard ratios. Under this layout they leave the figure
entirely. "No safe month" rests on them — all six months exclude 1 under
Bonferroni — whereas on the RD scale month 6 crosses zero at both horizons
(2 y: -0.10, -0.70 to 0.61; 5 y: +0.94, -0.09 to 2.13). Owner's position
(2026-08-20) is that a CI grazing zero does not flip the message. Keep that in
mind when wording the claim: the figure will not show it.

SOURCES
-------
  A  ITT_Analysis/results/rolling_rd_by_month_boot.csv   (script 45b, B=300, M=5)
  B  CCW_analysis/results_v3/ccw_subgroups_h60_bootstrap.csv  (B=300, M=5, 60 mo)
  C  same file

Two subgroup cells are marked as not reportable in docs/number-registry.csv
because their intervals include zero (age >=65, homeless). They are DRAWN, in
grey, rather than silently dropped: a reader who sees five age bands and four
housing categories elsewhere would otherwise wonder what happened to them.

Usage:  python3 73_make_fig4_revised.py
        HZ=2 python3 73_make_fig4_revised.py     # panel A at 2 years instead of 5
"""
import os
import pathlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[2]
RES = ROOT / "ITT_Analysis" / "results"
CCW = ROOT / "CCW_analysis" / "results_v3"

HZ = float(os.environ.get("HZ", "5"))          # horizon for panel A

# Okabe-Ito, matching the design schematic. No red/green pairing (PLOS request).
BLUE, VERM, GREY = "#0072B2", "#D55E00", "#999999"

# ---------------------------------------------------------------- data ------
a = pd.read_csv(RES / "rolling_rd_by_month_boot.csv")
a = a[(a.dmon.notna()) & (a.time_y == HZ)].sort_values("dmon")

s = pd.read_csv(CCW / "ccw_subgroups_h60_bootstrap.csv")

# Order the subgroup rows for display, grouped by variable.
ORDER = [
    ("age_group", "15-24", "15–24"),
    ("age_group", "25-44", "25–44"),
    ("age_group", "45-64", "45–64"),
    ("age_group", "≥65", "≥65"),
    ("hiv_aids", "Positive", "HIV positive"),
    ("hiv_aids", "Negative", "HIV negative"),
    ("homelessness", "Yes", "Homeless"),
    ("homelessness", "No", "Housed"),
    ("sex", "Female", "Female"),
    ("sex", "Male", "Male"),
]
# Drug resistance and calendar period dropped (owner, 2026-08-21). The resistant
# cell rested on 167 lost-to-follow-up patients and "not tested" is not a clinical
# stratum; period is a health-system question, not effect modification.
rows = []
for sg, lv, lab in ORDER:
    m = s[(s.subgroup == sg) & (s.level == lv)]
    if not len(m):
        continue
    r = m.iloc[0]
    rows.append(dict(label=lab, group=sg, n=int(r["n"]),
                     rr=r["rr"], rr_lo=r["rr_lo"], rr_hi=r["rr_hi"],
                     rd=r["rd"], rd_lo=r["rd_lo"], rd_hi=r["rd_hi"]))
sub = pd.DataFrame(rows)
# Intervals spanning the null are NOT recoloured (owner, 2026-08-21): the
# interval already shows it, and a second visual channel implied a significance
# threshold the analysis does not use.

# ---------------------------------------------------------------- figure ----
fig = plt.figure(figsize=(15, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.05, 1.0], wspace=0.42)
# A and B are the two CCW forests and share the subgroup axis; C is the
# landmark timing panel and is independent.
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1], sharey=axA)
axC = fig.add_subplot(gs[0, 2])


def _panel_label(ax, letter, title):
    ax.set_title(title, fontsize=10.5, weight="bold", loc="left", pad=16)
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes,
            fontsize=14, weight="bold", va="bottom", ha="left")


# ---- A: RD by month, rolling landmark --------------------------------------
x = a.dmon.values
axC.axhline(0, color="black", lw=0.8, zorder=1)
axC.errorbar(x, a.rd, yerr=[a.rd - a.rd_lo, a.rd_hi - a.rd],
             fmt="o", ms=7, lw=0, elinewidth=1.6, capsize=4,
             color=VERM, ecolor=VERM, zorder=3)
for xi, rd, lo, hi, n in zip(x, a.rd, a.rd_lo, a.rd_hi, a.n_exposed):
    axC.annotate(f"n={n:,}", (xi, hi), textcoords="offset points",
                 xytext=(0, 7), ha="center", fontsize=7, color="#555555")
axC.set_xticks(range(1, 7))
axC.set_xlabel("Month of disengagement")
axC.set_ylabel(f"Risk difference at {int(HZ)} years (percentage points)")
_panel_label(axC, "C", "Rolling landmark: timing of disengagement")
axC.spines[["top", "right"]].set_visible(False)
axC.margins(x=0.10)

# ---- B and C: subgroups, CCW ------------------------------------------------
y = np.arange(len(sub))[::-1]


def forest(ax, est, lo, hi, null_at, xlabel, logx=False):
    ax.axvline(null_at, color="black", lw=0.8, zorder=1)
    for yi, e, l, h in zip(y, est, lo, hi):
        c = BLUE
        ax.plot([l, h], [yi, yi], color=c, lw=1.6, zorder=2)
        ax.plot([e], [yi], "o", ms=6, color=c, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(sub.label)
    ax.set_xlabel(xlabel)
    if logx:
        ax.set_xscale("log")
        ax.set_xticks([0.8, 1.0, 1.5, 2.0])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=0.03)
    # faint separators between subgroup variables
    for i in range(1, len(sub)):
        if sub.group.iloc[i] != sub.group.iloc[i - 1]:
            ax.axhline(y[i] + 0.5, color="#DDDDDD", lw=0.8, zorder=0)


forest(axA, sub.rr, sub.rr_lo, sub.rr_hi, 1.0,
       "Risk ratio at 60 months", logx=True)
_panel_label(axA, "A", "Clone-censor-weight: relative effect")

forest(axB, sub.rd, sub.rd_lo, sub.rd_hi, 0.0,
       "Risk difference at 60 months (percentage points)")
_panel_label(axB, "B", "Clone-censor-weight: absolute effect")
plt.setp(axB.get_yticklabels(), visible=False)

# The grey-interval sentence is gone with the grey. The two-clocks caveat is NOT
# optional -- CLAUDE.md forbids reading a landmark window against a CCW horizon --
# so it is kept, repointed at panel C.
fig.text(0.5, -0.02,
         "Panels A and B are measured from treatment start; panel C from each patient's "
         "loss-to-follow-up declaration date. The two clocks are not interchangeable and "
         "the magnitudes cannot be read against one another.",
         ha="center", fontsize=8, color="#555555", wrap=True)

fig.tight_layout()
for ext in ("png", "pdf"):
    out = RES / f"Figure_4_revised_{int(HZ)}y.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)
