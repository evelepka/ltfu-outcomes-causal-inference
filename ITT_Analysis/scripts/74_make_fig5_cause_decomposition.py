#!/usr/bin/env python3
"""Figure 5 — what the excess mortality after disengagement is made of.

Panel A  overall: the all-cause risk difference and its three components.
Panel B  by month of disengagement: the same decomposition, month by month.

WHY THE THREE PARTS ADD UP
--------------------------
The components come from 46d's Aalen-Johansen construction, not from 1 - KM per
cause. Censoring a competing-cause death at its own event time answers "risk of
dying of cause k in a world where nobody dies of anything else", which is a
hypothetical and whose parts do not sum to the all-cause risk. 46d assembles all
three cause-specific hazards into one survival function, so the parts sum to the
whole by construction -- the printed CHECK line in 46d is the gap, and it is
what makes this figure honest to add up by eye.

The unclassified cause is DRAWN, not folded away. Roughly an eighth of deaths
carry no usable cause; hiding them would make TB and non-TB look like a
partition of the deaths when they are not, and would leave a visible gap between
the stack and the all-cause bar with no explanation.

CAUTION ON READING PANEL B AGAINST FIGURE 4
-------------------------------------------
Both are on the same clock (each patient's own LTFU declaration date), so the
months correspond. But Figure 4 panel A comes from 45b, which fits ONE all-cause
Cox model, while the all-cause figure here is assembled from three cause-specific
models. Those are not algebraically identical, and the two differ slightly in the
middle months. If the owner takes option (b) in docs/status-2026-08-20-night.md,
Figure 4 panel A is redrawn from this same file and the difference disappears.

INPUT
-----
  rolling_cause_cif_boot.csv          full bootstrap (46d, B=180)
  rolling_cause_cif_boot_partial.csv  fallback, written by 46e from a partial run

Whichever is used, the replicate count is read from the file and printed into
the figure footnote rather than assumed.

Usage:  python3 74_make_fig5_cause_decomposition.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ITT_RESULTS_DIR                        # noqa: E402

RES = Path(ITT_RESULTS_DIR)
HZ = 5

# Okabe-Ito. No red/green pairing (PLOS request); grey reserved for "spans zero".
COL = {"tb": "#D55E00", "nontb": "#0072B2", "unclass": "#E69F00", "all": "#000000"}
# NOTE ON THE NULL CONVENTION -- deliberately different from Figure 4.
# Figure 4 greys out an interval spanning the null, which is fine there because
# every row is a separately labelled category. Here four series are interleaved
# at each month and greying one destroys the only cue to which cause it is. So
# the cause keeps its colour and the marker is drawn OPEN instead. If the two
# figures are ever harmonised, move Figure 4 to open markers, not this one to
# grey.
# Parallel "attributed to" wording, deliberately. An earlier draft used
# "Other causes" for the second category, which reads as "everything that is not
# tuberculosis" and therefore appears to swallow the third -- the categories
# looked overlapping when they are in fact a partition. The second category
# requires a POSITIVE classification as something other than tuberculosis; the
# third is the residue that cannot be classified at all.
LAB = {"tb": "Attributed to tuberculosis",
       "nontb": "Attributed to another cause",
       "unclass": "No cause attributed",
       "all": "All causes"}
# Panel B interleaves four series per month, so its legend has to stay narrow.
# Panel A spells the categories out in full a few centimetres to the left, and
# the colours are shared, so the legend can be terse without losing the reader.
SHORT = {"all": "All causes", "tb": "Tuberculosis",
         "nontb": "Another cause", "unclass": "Not attributed"}
GREY = "#999999"
ORDER = ["all", "tb", "nontb", "unclass"]


def load():
    for name in ("rolling_cause_cif_boot.csv", "rolling_cause_cif_boot_partial.csv"):
        p = RES / name
        if p.exists():
            d = pd.read_csv(p)
            if "n_reps" in d and d.n_reps.notna().any():
                # Report the RANGE, not the maximum. Replicates that could not
                # estimate a given month still contributed to others, so cells
                # differ in how many replicates stand behind them; quoting the
                # best-supported cell would overstate the thinnest one.
                n = d[d.time_y == HZ].n_reps.dropna()
                lo, hi = int(n.min()), int(n.max())
                return d, name, (lo if lo == hi else f"{lo}-{hi}")
    sys.exit("no bootstrapped CIF file yet -- run 46d, or 46e for a partial run")


def spans_zero(r):
    return pd.notna(r.rd_lo) and pd.notna(r.rd_hi) and r.rd_lo < 0 < r.rd_hi


def main() -> int:
    d, src, nrep = load()
    d = d[d.time_y == HZ]
    ov = d[d.dmon.isna()].set_index("cause")
    bym = d[d.dmon.notna()]

    fig = plt.figure(figsize=(13, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.35], wspace=0.28)
    axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # ---- A: overall -------------------------------------------------------
    y = np.arange(len(ORDER))[::-1]
    for yi, c in zip(y, ORDER):
        if c not in ov.index:
            continue
        r = ov.loc[c]
        null = spans_zero(r)
        axA.plot([r.rd_lo, r.rd_hi], [yi, yi], color=COL[c], lw=2.0,
                 alpha=0.45 if null else 1.0, zorder=2)
        axA.plot([r.rd], [yi], "o", ms=8 if c == "all" else 6.5,
                 mfc="white" if null else COL[c], mec=COL[c], mew=1.6, zorder=3)
    axA.axvline(0, color="black", lw=0.8, zorder=1)
    axA.axhline(y[0] - 0.5, color="#DDDDDD", lw=0.9, zorder=0)   # all vs parts
    axA.set_yticks(y)
    axA.set_yticklabels([LAB[c] for c in ORDER])
    axA.set_xlabel(f"Risk difference at {HZ} years (percentage points)")
    axA.set_title("Overall", fontsize=10.5, weight="bold", loc="left", pad=14)
    axA.text(-0.42, 1.06, "A", transform=axA.transAxes, fontsize=14,
             weight="bold", va="bottom", ha="left")
    axA.spines[["top", "right"]].set_visible(False)
    axA.margins(y=0.14)

    # ---- B: by month ------------------------------------------------------
    months = sorted(bym.dmon.unique())
    offs = {"all": -0.27, "tb": -0.09, "nontb": 0.09, "unclass": 0.27}
    for c in ORDER:
        for m in months:
            s = bym[(bym.dmon == m) & (bym.cause == c)]
            if not len(s):
                continue
            r = s.iloc[0]
            null = spans_zero(r)
            x = m + offs[c]
            axB.plot([x, x], [r.rd_lo, r.rd_hi], color=COL[c], lw=1.5,
                     alpha=0.45 if null else 1.0, zorder=2)
            axB.plot([x], [r.rd], "o", ms=5.5,
                     mfc="white" if null else COL[c], mec=COL[c], mew=1.4, zorder=3)
    for c in ORDER:                                   # legend proxies
        axB.plot([], [], "o-", color=COL[c], label=SHORT[c], ms=5.5, lw=1.5)
    axB.axhline(0, color="black", lw=0.8, zorder=1)
    axB.set_xticks(months)
    axB.set_xlabel("Month of disengagement")
    axB.set_ylabel(f"Risk difference at {HZ} years (pp)")
    axB.set_title("By month of disengagement", fontsize=10.5, weight="bold",
                  loc="left", pad=14)
    axB.text(-0.09, 1.06, "B", transform=axB.transAxes, fontsize=14,
             weight="bold", va="bottom", ha="left")
    # headroom so the legend clears the tallest interval
    lo, hi = axB.get_ylim()
    axB.set_ylim(lo, hi + 0.24 * (hi - lo))
    axB.legend(frameon=False, fontsize=8.5, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, 1.0), columnspacing=1.2, handletextpad=0.4)
    axB.spines[["top", "right"]].set_visible(False)
    axB.margins(x=0.06)

    fig.text(0.5, -0.04,
             "Cause-specific cumulative incidence (Aalen-Johansen) standardised to the disengaging "
             f"population, from each patient's own loss-to-follow-up declaration date. {nrep} cluster-bootstrap "
             "replicates resampling patients. The three cause categories are mutually exclusive and exhaustive: "
             "the second requires a positive attribution to something other than tuberculosis, the third is the "
             "residue that cannot be attributed at all. They sum to the all-cause estimate by construction. "
             "Open markers indicate an interval spanning zero.",
             ha="center", fontsize=8, color="#555555", wrap=True)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = RES / f"Figure_5_cause_decomposition.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print("wrote", out)
    print(f"  source: {src}  ({nrep} replicates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
