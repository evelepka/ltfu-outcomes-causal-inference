#!/usr/bin/env python3
"""Figures S10 and S11 - severity-stratified bar charts.

S10: on-treatment cumulative mortality at 6 and 24 months, by HIV / hospitalised /
     smear status   (source: ontx_severity_stratified.csv)
S11: competing-risks cumulative incidence by month 3 (death-on-treatment vs LTFU),
     by HIV / hospitalised status   (source: competing_risks_by_severity.csv)

Fix vs prior version: axis labels no longer read "HIV-positive: Negative" /
"Smear-positive: Negative" (contradictory). Each stratum is labelled by its
variable name with the level beneath it.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results"

SHORT = {"HIV-positive": "HIV", "Hospitalised at diagnosis": "Hospitalised",
         "Smear-positive": "Smear"}


# ---------------- Figure S10 ----------------
def fig_s10():
    d = pd.read_csv(os.path.join(R, "ontx_severity_stratified.csv"))
    order = [("HIV-positive", "Negative"), ("HIV-positive", "Positive"),
             ("Hospitalised at diagnosis", "No"), ("Hospitalised at diagnosis", "Yes"),
             ("Smear-positive", "Negative"), ("Smear-positive", "Positive")]
    labels, m6, m24 = [], [], []
    for sv, lv in order:
        row = d[(d.stratum_var == sv) & (d.subgroup == lv)].iloc[0]
        labels.append(f"{SHORT[sv]}\n{lv}")
        m6.append(row.cum_mort_6mo_pct); m24.append(row.cum_mort_24mo_pct)

    x = np.arange(len(order)); w = 0.38
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    ax.bar(x - w/2, m6,  w, color="#8c8c8c", edgecolor="black", lw=0.6,
           label="6-month cumulative mortality")
    ax.bar(x + w/2, m24, w, color="white", edgecolor="black", lw=0.6,
           label="24-month cumulative mortality")
    for xi, v in zip(x - w/2, m6):
        ax.text(xi, v + 0.08, f"{v:.1f}", ha="center", va="bottom", fontsize=8.5)
    for xi, v in zip(x + w/2, m24):
        ax.text(xi, v + 0.08, f"{v:.1f}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Cumulative mortality (%)", fontsize=11)
    ax.set_ylim(0, max(m24) * 1.18)
    ax.legend(frameon=False, fontsize=10.5, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S10_bars.png", dpi=300, facecolor="white",
                bbox_inches="tight", pad_inches=0.1)
    print("wrote /Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S10_bars.png")


# ---------------- Figure S11 ----------------
def fig_s11():
    d = pd.read_csv(os.path.join(R, "competing_risks_by_severity.csv"))
    d = d[d.time_mo == 3]
    order = [("Hospitalised at diagnosis", "No"), ("Hospitalised at diagnosis", "Yes"),
             ("HIV-positive", "Negative"), ("HIV-positive", "Positive")]
    labels, death, ltfu = [], [], []
    for sv, lv in order:
        sub = d[(d.stratum_var == sv) & (d.subgroup == lv)]
        labels.append(f"{SHORT[sv]}\n{lv}")
        death.append(sub[sub.event == "death_on_tx"].cif_pct.iloc[0])
        ltfu.append(sub[sub.event == "ltfu"].cif_pct.iloc[0])

    x = np.arange(len(order)); w = 0.38
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    ax.bar(x - w/2, death, w, color="#c0392b", edgecolor="black", lw=0.6,
           label="Death on treatment")
    ax.bar(x + w/2, ltfu,  w, color="#2c7fb8", edgecolor="black", lw=0.6,
           label="Loss to follow-up")
    for xi, v in zip(x - w/2, death):
        ax.text(xi, v + 0.2, f"{v:.1f}", ha="center", va="bottom", fontsize=8.5)
    for xi, v in zip(x + w/2, ltfu):
        ax.text(xi, v + 0.2, f"{v:.1f}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Cumulative incidence by month 3 (%)", fontsize=11)
    ax.set_ylim(0, max(death) * 1.18)
    ax.legend(frameon=False, fontsize=10.5, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S11_bars.png", dpi=300, facecolor="white",
                bbox_inches="tight", pad_inches=0.1)
    print("wrote /Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S11_bars.png")


fig_s10()
fig_s11()
