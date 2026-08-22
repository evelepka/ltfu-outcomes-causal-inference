#!/usr/bin/env python3
"""The primary result, displayed: standardised risk under each assigned strategy.

Nothing in the manuscript plotted the primary estimand. The clone-censor-weight
contrast existed only as numbers in the text while a whole figure went to crude
within-LTFU curves. This plots it, overall and for people living with HIV — the
subgroup with the largest absolute excess and the one Reviewer 1 named.

WHAT THESE CURVES ARE
---------------------
Weighted, model-standardised cumulative risk under two assigned strategies, NOT
Kaplan-Meier. Each arm is a reweighted set of clones with
inverse-probability-of-censoring weights correcting departure from the assigned
strategy. Say so in the caption — the figure itself carries no explanatory text,
by owner's instruction.

TWO THINGS THE READER WILL NOTICE
---------------------------------
1. The curves are indistinguishable for the first ~6 months. That is construction,
   not finding: both clones of every patient start at treatment initiation and a
   death before the strategies diverge counts in both arms. It is also why the
   risk DIFFERENCE rather than the ratio is the interpretable contrast.
2. Risk at month 1 is already ~4%. That is real -- 6,868 of 171,048 patients
   (4.0%) die within 30 days of starting treatment. The curves are anchored at
   (0, 0) so this reads as steep early mortality rather than a truncated axis.

The subgroup panel drops hiv_aids from the covariate set, matching how ccw_v3
computes subgroup contrasts.

Usage:  python3 ITT_Analysis/scripts/72_make_fig3_primary_curves.py
"""
import importlib.util
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
spec = importlib.util.spec_from_file_location("ccw", ROOT / "CCW_analysis/ccw_v3.py")
ccw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccw)

OUTDIR = ROOT / "Plos Medicine" / "R1" / "Figures" / "Draft"
RESULTS = ROOT / "ITT_Analysis" / "results"
STEM = "LTFU_primary_curves"

H, T = 60, 6
BLUE, VERM, GREY, INK = "#0072B2", "#D55E00", "#8A8A8A", "#222222"
# The difference is a CONTRAST between the two arms, not a property of either, so
# it gets a neutral colour rather than the exposed arm's vermillion. The shaded
# gap in the top row is the same quantity the bottom row plots, so it matches.
# Deep purple (owner's choice from the option sheet). The difference is a CONTRAST
# between the arms, not a property of either, so it gets a third hue rather than the
# exposed arm's vermillion. A warm hue was rejected on purpose: this colour also
# fills the gap between the arm curves, and a warm fill visually annexes the
# difference to the vermillion arm -- the confusion the recolour exists to remove.
DIFF = "#6A4C93"
ccw.HORIZON_M = H

lookup = ccw.build_cause_lookup(verbose=False)
TLS = [ccw.load_timeline(p, lookup, verbose=False)
       for p in sorted(ccw.MI_DIR.glob("imp_*.csv"))]
print(f"[72] {len(TLS)} imputations")


def curves(level=None):
    """Pooled arm curves. level=('hiv_aids','Positive') restricts and drops that
    covariate, exactly as ccw_v3 does for subgroup contrasts."""
    covs = ccw.COVS if level is None else [c for c in ccw.COVS if c != level[0]]
    got = {"disengage": [], "remain": []}
    for tl in TLS:
        s = tl if level is None else tl[tl[level[0]].to_numpy() == level[1]]
        s, X = ccw.attach_patterns(s.reset_index(drop=True), covs)
        for arm in ("disengage", "remain"):
            a = ccw.build_arm(s, arm, T=T, cause="all_cause")
            got[arm].append(ccw.weighted_risk(ccw.add_ipcw(ccw.expand(a), X,
                                                           verbose=False)))
    r1 = 100 * np.mean(np.vstack(got["disengage"]), axis=0)
    r0 = 100 * np.mean(np.vstack(got["remain"]), axis=0)
    # anchor at the origin: without a month-0 point the curve appears to start at 4%
    x = np.concatenate([[0], np.arange(1, H + 1)])
    return x, np.concatenate([[0.0], r1]), np.concatenate([[0.0], r0])


# ---------------------------------------------------------------- bootstrap --
# B=300 python3 72_... runs it; without B the figure uses any cached CSV.
# Resamples PATIENTS and rebuilds both arms per replicate, retaining the FULL
# 60-month risk-difference curve rather than only the endpoint -- which is what
# the previous bootstrap stored, and why a pointwise band was not available.
B = int(os.environ.get("B", "0"))
BAND_CSV = RESULTS / "fig3_primary_curve_bands.csv"


def boot_band(level, B, seed):
    rng = np.random.default_rng(seed)
    covs = ccw.COVS if level is None else [c for c in ccw.COVS if c != level[0]]
    base = [tl if level is None else tl[tl[level[0]].to_numpy() == level[1]]
            for tl in TLS]
    out = []
    t0 = time.time()
    for b in range(B):
        src = base[rng.integers(len(base))]
        idx = rng.integers(0, len(src), len(src))
        s = src.iloc[idx].reset_index(drop=True)
        s, X = ccw.attach_patterns(s, covs)
        try:
            cur = {}
            for arm in ("disengage", "remain"):
                a = ccw.build_arm(s, arm, T=T, cause="all_cause")
                cur[arm] = ccw.weighted_risk(ccw.add_ipcw(ccw.expand(a), X,
                                                          verbose=False))
            out.append(100 * (cur["disengage"] - cur["remain"]))
        except Exception:
            continue
        if (b + 1) % 25 == 0:
            el = time.time() - t0
            print(f"    rep {b+1}/{B}  {el/(b+1):.2f}s/rep  "
                  f"ETA {el/(b+1)*(B-b-1)/60:.1f} min", flush=True)
    return np.vstack(out)


bs = pd.read_csv(ROOT / "CCW_analysis/results_v3/ccw_bootstrap_h60.csv")
ov = bs[bs.estimate == "nested_T6"].iloc[0]
sg = pd.read_csv(ROOT / "CCW_analysis/results_v3/ccw_subgroups_h60_bootstrap.csv")
hp = sg[(sg.subgroup == "hiv_aids") & (sg.level == "Positive")].iloc[0]

PANELS = [
    ("All patients", None, (ov["rd_lo"], ov["rd_hi"]), 14.0),
    ("People living with HIV", ("hiv_aids", "Positive"),
     (hp["rd_lo"], hp["rd_hi"]), 38.0),
]

# ---- bootstrap bands, computed once and cached ------------------------------
bands = {}
if B > 0:
    for title, level, _, _ in PANELS:
        print(f"  bootstrapping {title} (B={B})")
        bands[title] = boot_band(level, B, seed=abs(hash(title)) % 10_000)
    rows_b = []
    for title, arr in bands.items():
        for m in range(arr.shape[1]):
            rows_b.append(dict(group=title, month=m + 1,
                               rd_lo=np.percentile(arr[:, m], 2.5),
                               rd_hi=np.percentile(arr[:, m], 97.5),
                               n_reps=arr.shape[0]))
    pd.DataFrame(rows_b).to_csv(BAND_CSV, index=False)
    print(f"wrote {BAND_CSV.relative_to(ROOT)}")
band_df = pd.read_csv(BAND_CSV) if BAND_CSV.exists() else None
if band_df is None:
    print("  NOTE: no cached bands; run with B=200 to draw the lower panels")

# ---------------------------------------------------------------- figure -----
fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0),
                         gridspec_kw={"height_ratios": [2.0, 1.15]})
rows, curve_rows = [], []
LETTERS = [("A", "C"), ("B", "D")]
for col, (title, level, (lo, hi), ymax) in enumerate(PANELS):
    axT, axD = axes[0, col], axes[1, col]
    x, r1, r0 = curves(level)
    rd_curve = r1 - r0
    rd = rd_curve[-1]
    print(f"  {title}: {r1[-1]:.2f}% vs {r0[-1]:.2f}%, RD {rd:+.2f} pp")
    rows.append(dict(group=title, risk_disengage=r1[-1], risk_remain=r0[-1],
                     rd=rd, rd_lo=lo, rd_hi=hi))
    # the per-month point curve, saved because the nadir and the crossing month
    # are quoted in the text and must come from a file, not from reading the plot
    for m in range(1, H + 1):
        curve_rows.append(dict(group=title, month=m, risk_disengage=r1[m],
                               risk_remain=r0[m], rd=rd_curve[m]))

    # ---- top: standardised risk under each strategy ----
    axT.fill_between(x, r0, r1, color=DIFF, alpha=0.13, lw=0, zorder=1)
    axT.plot(x, r1, color=VERM, lw=2.0, zorder=3,
             label=f"Disengaged within {T} months")
    axT.plot(x, r0, color=BLUE, lw=2.0, zorder=3, label="Remained in care")
    axT.annotate(f"{r1[-1]:.1f}%", (H, r1[-1]), xytext=(3, 1),
                 textcoords="offset points", fontsize=7.8, color=VERM,
                 weight="bold", va="center")
    axT.annotate(f"{r0[-1]:.1f}%", (H, r0[-1]), xytext=(3, -2),
                 textcoords="offset points", fontsize=7.8, color=BLUE,
                 weight="bold", va="center")
    axT.set_ylim(0, ymax)
    axT.set_ylabel("Cumulative mortality (%)", fontsize=8.2)
    axT.set_title(f"{LETTERS[col][0]}  {title}", fontsize=9, weight="bold",
                  loc="left", pad=6)

    # ---- bottom: the difference, with its band ----
    # A band on the risk curves themselves would overlap and invite the
    # "intervals overlap so it is null" misreading. The difference and its own
    # interval answer the question directly.
    axD.axhline(0, color=INK, lw=0.8, ls=(0, (4, 3)), zorder=2)
    if band_df is not None:
        bd = band_df[band_df.group == title].sort_values("month")
        axD.fill_between(bd.month, bd.rd_lo, bd.rd_hi, color=DIFF, alpha=0.22,
                         lw=0, zorder=1)
        ylo = min(bd.rd_lo.min(), -0.6) * 1.15
        yhi = max(bd.rd_hi.max(), rd) * 1.15
    else:
        ylo, yhi = -1.0, rd * 1.4
    axD.plot(x[1:], rd_curve[1:], color=DIFF, lw=1.8, zorder=3)
    axD.set_ylim(ylo, yhi)
    axD.set_ylabel("Risk difference (pp)", fontsize=8.2)
    axD.set_title(f"{LETTERS[col][1]}", fontsize=9, weight="bold", loc="left", pad=4)
    axD.annotate(f"{rd:+.2f}\n({lo:.2f} to {hi:.2f})", (H, rd),
                 xytext=(4, 0), textcoords="offset points", fontsize=7.6,
                 color=INK, weight="bold", va="center", linespacing=1.3)

    for ax_ in (axT, axD):
        ax_.set_xlim(0, H + 9)
        ax_.set_xticks([0, 12, 24, 36, 48, 60])
        ax_.tick_params(labelsize=7.8)
        for s_ in ("top", "right"):
            ax_.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"):
            ax_.spines[s_].set_color(GREY)
    axD.set_xlabel("Months since treatment initiation", fontsize=8.2)

axes[0, 0].legend(frameon=False, fontsize=7.8, loc="lower right",
                  bbox_to_anchor=(1.0, 0.02))

fig.tight_layout(pad=0.7)
OUTDIR.mkdir(parents=True, exist_ok=True)
# bbox_inches=tight: the bottom row's right-hand annotation is wider than the top
# row's, so plain tight_layout shifts those axes left and clips panel C's ylabel.
for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    q = OUTDIR / f"{STEM}.{ext}"
    fig.savefig(q, facecolor="white", bbox_inches="tight", pad_inches=0.12, **kw)
    print(f"wrote {q.relative_to(ROOT)}")
plt.close(fig)
pd.DataFrame(rows).to_csv(RESULTS / "fig3_primary_curves.csv", index=False)
pd.DataFrame(curve_rows).to_csv(RESULTS / "fig3_primary_rd_curve.csv", index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False))
