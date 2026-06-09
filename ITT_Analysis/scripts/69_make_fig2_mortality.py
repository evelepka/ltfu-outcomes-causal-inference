"""
69_make_fig2_mortality.py
Descriptive within-LTFU mortality figure (replacement for the stratified-
retreatment Figure 2). All panels characterize mortality AMONG patients who
experienced LTFU -- distinct from the causal contrast (Figure 4) and the
adjusted cause-specific HRs (Figure 5).

Panels (5):
  A  Cumulative mortality by month of disengagement (KM curves), within LTFU.
  B  24-mo mortality by HIV status.
  C  24-mo mortality by hospitalization at diagnosis.
  D  24-mo mortality by housing status.
  E  24-mo mortality by age group.
Crude (unadjusted) Kaplan-Meier; time measured from LTFU (time_d).
Output: Figure_2_mortality_descriptive.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

ROOT = os.path.expanduser("~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper")
d = pd.read_csv(os.path.join(ROOT, "ITT_Analysis/data/itt_cohort.csv"), low_memory=False,
               parse_dates=["best_start", "end_date"])

def _clean(t, e):
    t = pd.to_numeric(pd.Series(t).reset_index(drop=True), errors="coerce")
    e = pd.to_numeric(pd.Series(e).reset_index(drop=True), errors="coerce")
    m = t.notna() & e.notna() & (t > 0)
    return t[m].values, e[m].astype(int).values

def km_fit(t, e):
    t, e = _clean(t, e); return KaplanMeierFitter().fit(t, e)

def km_risk(t, e, h=2.0):
    return float(1 - km_fit(t, e).predict(h))

ltfu = d[(d.itt_group == "Loss to follow-up") & (d.time_d > 0)].copy()
ltfu["txm"] = ((ltfu.end_date - ltfu.best_start).dt.days - 30).clip(lower=1) / 30.4
ltfu["dmonth"] = np.floor(ltfu.txm).astype(int) + 1

# Panel A: KM by month of disengagement (M1-M6)
B_bins = [(m, f"M{m}") for m in range(1, 7)]
A_km = {lab: km_fit(ltfu[ltfu.dmonth == m].time_d, ltfu[ltfu.dmonth == m].event_d) for m, lab in B_bins}

# Panels B-E: within-LTFU 24-mo mortality by risk factor
def strat(col, mapping):
    return [(lab, km_risk(ltfu[ltfu[col] == k].time_d, ltfu[ltfu[col] == k].event_d)*100,
             int((ltfu[col] == k).sum())) for k, lab in mapping]
HIV  = strat("hiv_aids", [("Positive", "HIV+"), ("Negative", "HIV−")])
HOSP = strat("hosp_admission", [("Yes", "Hospitalized"), ("No", "Not hospitalized")])
HOUS = strat("homelessness", [("Yes", "Homeless"), ("No", "Housed")])
AGE  = strat("age_group", [("15-24", "15–24"), ("25-44", "25–44"), ("45-64", "45–64"), ("65+", "≥65")])
ALC  = strat("alcohol", [("Yes", "Alcohol use"), ("No", "No alcohol use")])
for nm, dd in [("HIV", HIV), ("HOSP", HOSP), ("HOUS", HOUS), ("AGE", AGE), ("ALC", ALC)]:
    print(nm, [(l, round(v, 1), n) for l, v, n in dd])

def barpanel(ax, dat, title):
    labs=[l for l,_,_ in dat]; vals=[v for _,v,_ in dat]; ns=[n for _,_,n in dat]
    cols = plt.cm.OrRd(np.linspace(0.45, 0.9, len(dat)))
    xx=np.arange(len(dat)); ax.bar(xx, vals, color=cols)
    for i in range(len(dat)):
        ax.text(xx[i], vals[i]+0.1, f"{vals[i]:.1f}%\n(n={ns[i]:,})", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xx); ax.set_xticklabels(labs); ax.set_ylabel("24-month mortality (%)")
    ax.set_title(title, fontsize=10, weight="bold"); ax.set_ylim(0, max(vals)*1.35)
    ax.spines[["top","right"]].set_visible(False)

fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
# A (x-axis in months)
ax = axes[0,0]; cmap = plt.cm.YlOrRd(np.linspace(0.4, 0.95, 6)); tl = np.linspace(0, 2, 100)
for i, (m, lab) in enumerate(B_bins):
    ax.plot(tl*12, (1-A_km[lab].predict(tl))*100, color=cmap[i], lw=2, label=lab)
ax.set_xlim(0, 24); ax.set_xticks([0, 6, 12, 18, 24])
ax.set_xlabel("Months since LTFU"); ax.set_ylabel("Cumulative mortality (%)")
ax.set_title("A. Mortality by month of disengagement", fontsize=10, weight="bold")
ax.legend(frameon=False, ncol=2, fontsize=8, title="Month of disengagement")
ax.spines[["top","right"]].set_visible(False)
barpanel(axes[0,1], HIV,  "B. HIV status")
barpanel(axes[0,2], HOSP, "C. Hospitalization at diagnosis")
barpanel(axes[1,0], HOUS, "D. Housing status")
barpanel(axes[1,1], AGE,  "E. Age at diagnosis (years)")
barpanel(axes[1,2], ALC,  "F. Alcohol use")
fig.suptitle("Mortality among patients lost to follow-up: by timing of disengagement (A) and risk factor (B–F)",
             fontsize=13, weight="bold")
fig.tight_layout(rect=[0,0,1,0.96])
out = os.path.join(ROOT, "ITT_Analysis/results/Figure_2_mortality_descriptive.png")
fig.savefig(out, dpi=300); print("wrote", out)
