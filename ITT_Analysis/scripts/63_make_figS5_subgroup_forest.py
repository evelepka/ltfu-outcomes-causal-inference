#!/usr/bin/env python3
"""Figure S5 - subgroup forest plot of late-window adjusted hazard ratios (cap = 2).

Primary defn-B + grace specification, MI-pooled. Sources (late window, cap = 2):
  age/sex/HIV/homelessness : target_trial_subgroup_interactions_grace_mi.csv
  drug resistance          : target_trial_resistance_grace_mi.csv
  calendar period          : target_trial_period_grace_mi.csv (stratified)

Fixes vs prior version: age categories in natural order; subgroups visually
grouped with bold headers; numeric aHR (95% CI) printed at right; log x-axis
with plain-number ticks.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import pandas as pd

R = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results"
OUT = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S5_forest.png"

sub  = pd.read_csv(os.path.join(R, "target_trial_subgroup_interactions_grace_mi.csv"))
res  = pd.read_csv(os.path.join(R, "target_trial_resistance_grace_mi.csv"))
per  = pd.read_csv(os.path.join(R, "target_trial_period_grace_mi.csv"))


def get(df, subgroup, level):
    r = df[(df.Subgroup == subgroup) & (df.Level == level) &
           (df.model == "late") & (df.cap == 2.0)].iloc[0]
    return r.HR, r.CI_L, r.CI_H


def get_per(period):
    r = per[(per.analysis == "stratified") & (per.period == period) &
            (per.cap == 2)].iloc[0]
    return r.HR, r.CI_L, r.CI_H


# (group header, [(label, HR, lo, hi), ...]) in display order (top -> bottom)
GROUPS = [
    ("Age, years", [
        ("15–24", *get(sub, "age_group", "15-24")),
        ("25–44", *get(sub, "age_group", "25-44")),
        ("45–64", *get(sub, "age_group", "45-64")),
        ("≥65",   *get(sub, "age_group", "≥65")),
    ]),
    ("Sex", [
        ("Male",   *get(sub, "sex", "Male")),
        ("Female", *get(sub, "sex", "Female")),
    ]),
    ("HIV status", [
        ("Negative", *get(sub, "hiv_aids", "Negative")),
        ("Positive", *get(sub, "hiv_aids", "Positive")),
    ]),
    ("Homelessness", [
        ("No",  *get(sub, "homelessness", "No")),
        ("Yes", *get(sub, "homelessness", "Yes")),
    ]),
    ("Drug resistance", [
        ("Drug-sensitive", *get(res, "resistance_clean", "Sensitive")),
        ("Drug-resistant", *get(res, "resistance_clean", "Resistant (Any)")),
        ("Not evaluated",  *get(res, "resistance_clean", "Not Evaluated")),
    ]),
    ("Calendar period", [
        ("Pre-COVID (2013–2019)",  *get_per("Pre-COVID (2013-2019)")),
        ("Post-COVID (2020–2023)", *get_per("Post-COVID (2020-2023)")),
    ]),
]

# Build row layout from the bottom up so the first group sits at the top.
rows = []           # (y, kind, payload)
y = 0.0
for gi, (header, items) in enumerate(reversed(GROUPS)):
    for label, hr, lo, hi in reversed(items):
        rows.append((y, "point", (label, hr, lo, hi)))
        y += 1.0
    rows.append((y, "header", header))
    y += 0.9
    if gi != len(GROUPS) - 1:
        y += 0.4   # extra gap between groups

ymax = y

fig, ax = plt.subplots(figsize=(8.6, 9.2))
ax.axvline(1.0, ls="--", color="0.55", lw=1.0, zorder=0)

XR = 6.6   # x position (in HR units, log) for the numeric column
for yy, kind, payload in rows:
    if kind == "header":
        ax.text(0.30, yy, payload, fontsize=10.5, fontweight="bold",
                va="center", ha="left")
    else:
        label, hr, lo, hi = payload
        ax.errorbar(hr, yy, xerr=[[hr - lo], [hi - hr]], fmt="o",
                    color="#1f1f1f", markersize=5.5, capsize=3, lw=1.4)
        ax.text(0.36, yy, label, fontsize=9.5, va="center", ha="left")
        ax.text(XR, yy, f"{hr:.2f} ({lo:.2f}–{hi:.2f})",
                fontsize=9, va="center", ha="left", family="DejaVu Sans")

ax.set_xscale("log")
ax.set_xlim(0.45, 6.3)
ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3, 5]))
ax.xaxis.set_minor_locator(FixedLocator([]))
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_ylim(-0.8, ymax)
ax.set_yticks([])
ax.set_xlabel("Late-window adjusted hazard ratio (log scale)", fontsize=11)
ax.set_xlim(0.32, 11.5)   # leave room for the numeric column on the right
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.text(XR, ymax - 0.3, "aHR (95% CI)", fontsize=9.5, fontweight="bold",
        va="center", ha="left")

fig.tight_layout()
fig.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.12)
print("wrote", OUT)
