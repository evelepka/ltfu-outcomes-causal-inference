#!/usr/bin/env python3
"""Figure S6 - complete-case (CC) vs multiply-imputed (MI) within-LTFU models.

Two panels: mortality (Cox) and retreatment (Fine-Gray). Each covariate shows
the MI-pooled and complete-case adjusted hazard ratio with 95% CI.

Fixes vs prior version: plain-number log x-axis (no 10^0 / 2x10^0 scientific
notation); legend moved below the panels so it no longer overlaps the row
labels; larger text.
"""
import os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import pandas as pd
import numpy as np

R = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results"
OUT = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S6_ccmi.png"

df = pd.read_csv(os.path.join(R, "multivariable_results_mi_cc.csv"))


def parse(est):
    m = re.match(r"\s*([\d.]+)\s*\(([\d.]+)[-–]([\d.]+)\)", str(est))
    if not m:
        return np.nan, np.nan, np.nan
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


LABELS = {
    "age_group25-44": "Age 25–44", "age_group45-64": "Age 45–64", "age_group≥65": "Age ≥65",
    "sexMale": "Male", "race_cleanBlack or Mixed": "Black or Mixed race",
    "race_cleanOther": "Other race", "edu_clean≤ 7 years": "Education ≤7 years",
    "edu_clean8 - 11 years": "Education 8–11 years", "hiv_aidsPositive": "HIV-positive",
    "diabetesYes": "Diabetes", "alcoholYes": "Alcohol use", "drug_useYes": "Drug use",
    "incarceratedYes": "Incarcerated", "homelessnessYes": "Homelessness",
    "hosp_admissionYes": "Hospitalised at diagnosis",
    "clinical_cleanExtrapulmonary": "Extrapulmonary TB",
    "clinical_cleanPulmonary and Extrapulmonary or disseminated": "Pulmonary + extrapulmonary",
    "dot_statusYes": "DOT", "tx_month_grp< 2 months": "LTFU <2 months",
    "tx_month_grp2 to <4 months": "LTFU 2 to <4 months",
}
# display order (top -> bottom) = the natural covariate order
ORDER = list(LABELS.keys())

MI_C, CC_C = "#1f1f1f", "#888888"
DODGE = 0.16

panels = [
    ("Cox_Death_Adjusted_MI", "Cox_Death_Adjusted_CC", "A. Mortality (Cox)"),
    ("FG_Retr_Adjusted_MI", "FG_Retr_Adjusted_CC", "B. Retreatment (Fine–Gray)"),
]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 8.2), sharey=True)

for ax, (mi_model, cc_model, title) in zip(axes, panels):
    mi = df[df.model == mi_model].set_index("term")
    cc = df[df.model == cc_model].set_index("term")
    terms = [t for t in ORDER if t in mi.index or t in cc.index]
    ypos = {t: i for i, t in enumerate(reversed(terms))}

    ax.axvline(1.0, ls="--", color="0.6", lw=1.0, zorder=0)
    for t in terms:
        y = ypos[t]
        if t in mi.index:
            hr, lo, hi = parse(mi.loc[t, "estimate"])
            ax.errorbar(hr, y + DODGE, xerr=[[hr - lo], [hi - hr]], fmt="o",
                        color=MI_C, markersize=5, capsize=2.5, lw=1.3,
                        label="Multiple imputation")
        if t in cc.index:
            hr, lo, hi = parse(cc.loc[t, "estimate"])
            ax.errorbar(hr, y - DODGE, xerr=[[hr - lo], [hi - hr]], fmt="s",
                        color=CC_C, markersize=5, capsize=2.5, lw=1.3,
                        label="Complete case")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([0.5, 0.7, 1, 1.5, 2, 3]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlim(0.45, 3.2)
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([LABELS[t] for t in reversed(terms)], fontsize=9.5)
    ax.set_title(title, fontsize=12, fontweight="bold", loc="left")
    ax.set_xlabel("Adjusted hazard ratio (log scale)", fontsize=11)
    ax.grid(axis="x", color="0.9", lw=0.6, zorder=0)
    ax.tick_params(labelsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# single shared legend below the panels (dedup handles)
handles, labels = axes[0].get_legend_handles_labels()
seen, H, L = set(), [], []
for h, l in zip(handles, labels):
    if l not in seen:
        seen.add(l); H.append(h); L.append(l)
fig.legend(H, L, loc="lower center", ncol=2, frameon=False, fontsize=11,
           bbox_to_anchor=(0.5, -0.01))

fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.1)
print("wrote", OUT)
