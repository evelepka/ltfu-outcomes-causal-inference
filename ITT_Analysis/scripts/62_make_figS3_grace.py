#!/usr/bin/env python3
"""Figure S3 - Sequential target-trial aHR by month of LTFU: no grace vs 30-day grace.

Two panels (early / late window). Y-axis is log-scaled (appropriate for hazard
ratios) but labelled with plain numbers (1, 2, 3 ...) rather than 10^0 scientific
notation.

Data: ITT_Analysis/results/target_trial_grace_vs_original.csv
  Panel A (early window, 6 mo):  model == 'early', cap == 0.5
  Panel B (late window, 6-24 mo): model == 'late',  cap == 2.0
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator
import pandas as pd

ROOT = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper"
CSV = os.path.join(ROOT, "ITT_Analysis/results/target_trial_grace_vs_original.csv")
OUT = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S3_grace.png"

BLUE, RED = "#1f77b4", "#d62728"
DODGE = 0.10

df = pd.read_csv(CSV)
df["m"] = df["Trial_Month"].str.replace("Month_", "").astype(int)

panels = [
    dict(ax=0, sel=(df.model == "early") & (df.cap == 0.5),
         title="A. Early window (first 6 months post-landmark)",
         yticks=[0.3, 0.5, 0.7, 1.0, 1.5, 2.0], ylim=(0.18, 2.6)),
    dict(ax=1, sel=(df.model == "late") & (df.cap == 2.0),
         title="B. Late window (6 to 24 months post-landmark)",
         yticks=[1.0, 1.5, 2.0, 2.5, 3.0, 3.5], ylim=(0.95, 3.9)),
]


def plain(v, _pos=None):
    s = f"{v:g}"
    return s


fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.95))

for p in panels:
    ax = axes[p["ax"]]
    d = df[p["sel"]].sort_values("m")
    x = d["m"].to_numpy()

    # no-grace series
    ax.errorbar(x - DODGE, d["HR_orig"],
                yerr=[d["HR_orig"] - d["CI_L_orig"], d["CI_H_orig"] - d["HR_orig"]],
                fmt="o", color=BLUE, markersize=7, capsize=3.5, lw=1.6,
                markeredgecolor="white", markeredgewidth=0.6,
                label="No grace period")
    # grace series
    ax.errorbar(x + DODGE, d["HR_grace"],
                yerr=[d["HR_grace"] - d["CI_L_grace"], d["CI_H_grace"] - d["HR_grace"]],
                fmt="s", color=RED, markersize=7, capsize=3.5, lw=1.6,
                markeredgecolor="white", markeredgewidth=0.6,
                label="Symmetric 30-day grace")

    ax.axhline(1.0, ls="--", color="0.55", lw=1.0, zorder=0)

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator(p["yticks"]))
    ax.yaxis.set_minor_locator(FixedLocator([]))      # no minor ticks/labels
    ax.yaxis.set_major_formatter(FuncFormatter(plain))
    ax.set_ylim(*p["ylim"])

    ax.set_xticks(range(1, 7))
    ax.set_xlim(0.5, 6.5)
    ax.set_title(p["title"], fontsize=12, fontweight="bold")
    ax.set_xlabel("Month of LTFU", fontsize=11)
    ax.set_ylabel("Adjusted hazard ratio (log scale)", fontsize=11)
    ax.grid(axis="y", color="0.88", lw=0.7, zorder=0)
    ax.tick_params(labelsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=10, loc="upper left")

fig.tight_layout(w_pad=2.5)
fig.savefig(OUT, dpi=360, facecolor="white", bbox_inches="tight", pad_inches=0.08)
print("wrote", OUT)
