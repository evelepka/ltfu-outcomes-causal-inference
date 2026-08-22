#!/usr/bin/env python3
"""Figure 5 — excess mortality by cause and by month of disengagement.

LAYOUT FOLLOWS "Updated Figures Plan.docx"
------------------------------------------
  "Figure 5. Will present a 2 year (top panel) and 5 year (bottom panel) RD
   based on Landmark model showing TB and non-TB mortality, by month of
   disengagement."

Two stacked panels, one horizon each, month of disengagement on x. An earlier
draft had an overall panel and a single horizon; both were mine, neither is in
the plan, and both are gone.

THE ONE DEVIATION, AND WHY
--------------------------
The plan says "TB and non-TB". A third series is drawn anyway: deaths with no
attributable cause, roughly an eighth of the total. The plan was written before
46d had been run and before that share was known. TB and non-TB do NOT partition
the deaths -- in 46d, ev_nontb requires a positive attribution to something other
than tuberculosis, not merely the absence of a tuberculosis attribution. Drawing
only two series would invite the reader to add them and read the sum as the
all-cause excess, which it is not.

It is subordinate to the other two: smaller marker, thinner line. TB and non-TB
carry the message, as the plan intends.

INPUT
-----
  rolling_cause_cif_boot.csv          46d with REPORT=2,5 and B=180
  rolling_cause_cif_boot_partial.csv  fallback, written by 46e from a partial run

The replicate count is read from the file rather than assumed, and a horizon the
file does not carry is drawn as an explicit "missing" panel rather than silently
omitted.

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

# Okabe-Ito. No red/green pairing (PLOS request); grey reserved for "spans zero".
COL = {"tb": "#D55E00", "nottb": "#0072B2", "all": "#000000"}
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
       "nottb": "Not attributed to tuberculosis",
       "all": "All causes"}
# Panel B interleaves four series per month, so its legend has to stay narrow.
# Panel A spells the categories out in full a few centimetres to the left, and
# the colours are shared, so the legend can be terse without losing the reader.
SHORT = {"all": "All causes", "tb": "Attributed to tuberculosis",
         "nottb": "Not attributed to tuberculosis"}
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
                n = d.n_reps.dropna()
                lo, hi = int(n.min()), int(n.max())
                return d, name, (lo if lo == hi else f"{lo}-{hi}")
    sys.exit("no bootstrapped CIF file yet -- run 46d, or 46e for a partial run")


def spans_zero(r):
    return pd.notna(r.rd_lo) and pd.notna(r.rd_hi) and r.rd_lo < 0 < r.rd_hi


HORIZONS = [5]
# The plan asked for a 2-year panel above a 5-year one. Owner decision
# 2026-08-21: five years only, matching the single horizon the manuscript
# reports, and the 2-year bootstrap was stopped rather than finished.
# TWO classes, owner decision 2026-08-21. The methodology names two, and a third
# series implied a classification the paper never had. "Not attributed to
# tuberculosis" is everything else -- positively attributed to another cause,
# plus the deaths whose cause cannot be assigned to either. That is literally
# true of them, and the two parts still sum to the all-cause excess.
#
# The residual is NOT small: the unassignable deaths are ~11-12% of the total
# and are irreducible with these data. Using tuberculosis from any line of the
# certificate instead of the underlying cause alone moves only 89 of them
# (12.0% -> 11.4%), so folding them in is a presentation choice, not a
# measurement improvement. The footnote must say what the second class contains.
SERIES = ["tb", "nottb"]
OFFS = {"tb": -0.10, "nottb": 0.10}
STY = {"tb": dict(ms=6.5, lw=1.8), "nottb": dict(ms=6.5, lw=1.8)}


def add_not_tb(d):
    """not-TB = all-cause minus TB. Derived per bootstrap replicate where the
    draws allow it, never by subtracting one interval from another."""
    import numpy as np
    draws = RES / "rolling_cause_cif_draws.csv"
    bd = pd.read_csv(draws, on_bad_lines="skip") if draws.exists() else None
    rows = []
    for (dm, hz), g in d.groupby([d.dmon.fillna(-1), "time_y"]):
        g = g.set_index("cause")
        if not {"all", "tb"} <= set(g.index):
            continue
        r = g.loc["all"].copy()
        r["cause"] = "nottb"
        r["rd"] = g.loc["all"].rd - g.loc["tb"].rd
        lo = hi = np.nan
        if bd is not None:
            b = bd[(bd.time_y == hz) & (bd.dmon.isna() if dm == -1 else bd.dmon == dm)]
            w = b.pivot_table(index="rep", columns="cause", values="rd")
            if {"all", "tb"} <= set(w.columns):
                diff = (w["all"] - w["tb"]).dropna()
                if len(diff) >= 20:
                    lo, hi = diff.quantile(.025), diff.quantile(.975)
        r["rd_lo"], r["rd_hi"] = lo, hi
        r["dmon"] = np.nan if dm == -1 else dm
        rows.append(r)
    return pd.concat([d, pd.DataFrame(rows)], ignore_index=True)


def main() -> int:
    d, src, nrep = load()
    d = add_not_tb(d)
    have = sorted(d.time_y.unique())
    missing = [h for h in HORIZONS if h not in have]
    if missing:
        print(f"  NOTE: horizon(s) {missing} absent from {src}; those panels are blank.")
        print("        Re-run 46d with REPORT=2,5.")

    fig, axes = plt.subplots(len(HORIZONS), 1, figsize=(9.2, 5.0), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, hz, letter in zip(axes, HORIZONS, "AB"):
        sub = d[(d.time_y == hz) & d.dmon.notna()]
        ax.axhline(0, color="black", lw=0.8, zorder=1)
        if not len(sub):
            ax.text(0.5, 0.5, f"{hz}-year horizon not present in\n{src}\n"
                              "re-run 46d with REPORT=2,5",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="#B00020")
        for c in SERIES:
            for m in sorted(sub.dmon.unique()):
                r = sub[(sub.dmon == m) & (sub.cause == c)]
                if not len(r):
                    continue
                r = r.iloc[0]
                null = spans_zero(r)
                x = m + OFFS[c]
                ax.plot([x, x], [r.rd_lo, r.rd_hi], color=COL[c],
                        lw=STY[c]["lw"], alpha=0.45 if null else 1.0, zorder=2)
                ax.plot([x], [r.rd], "o", ms=STY[c]["ms"],
                        mfc="white" if null else COL[c], mec=COL[c], mew=1.4,
                        zorder=3)
        ax.set_ylabel(f"Risk difference at {hz} years\n(percentage points)")
        ax.set_title(f"Excess mortality by cause and month of disengagement, {hz} years",
                     fontsize=10.5, weight="bold", loc="left", pad=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(x=0.07)

    for c in SERIES:
        axes[0].plot([], [], "o-", color=COL[c], label=SHORT[c],
                     ms=STY[c]["ms"], lw=STY[c]["lw"])
    lo, hi = axes[0].get_ylim()
    axes[0].set_ylim(lo, hi + 0.22 * (hi - lo))
    axes[0].legend(frameon=False, fontsize=9, ncol=3, loc="upper center",
                   columnspacing=1.6, handletextpad=0.4)

    xs = sorted(d[d.dmon.notna()].dmon.unique())
    if xs:
            axes[-1].set_xticks(xs)
    axes[-1].set_xlabel("Month of disengagement")

    fig.text(0.5, -0.12,
             "Cause-specific cumulative incidence (Aalen-Johansen) standardised to the disengaging "
             f"population, from each patient's own loss-to-follow-up declaration date. {nrep} cluster-bootstrap "
             "replicates resampling patients. The second class is every death not attributed to tuberculosis: "
             "those positively attributed to another cause, and those whose cause cannot be assigned to either "
             "(about an eighth of deaths, chiefly respiratory and HIV codes carrying no mention of "
             "tuberculosis). The two classes sum to the all-cause excess. Open markers indicate an interval "
             "spanning zero.",
             ha="center", fontsize=8, color="#555555", wrap=True)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.20)
    for ext in ("png", "pdf"):
        out = RES / f"Figure_5_cause_decomposition.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print("wrote", out)
    print(f"  source: {src}  ({nrep} replicates)  horizons present: {have}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
