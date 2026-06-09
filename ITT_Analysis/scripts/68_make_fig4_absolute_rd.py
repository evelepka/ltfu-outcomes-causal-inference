"""
68_make_fig4_absolute_rd.py
Absolute-scale companion panels for Figure 4:
  - Panel D: standardized 24-month all-cause mortality risk difference (RD) by
    month of disengagement (from 30o).
  - Panel C (right): subgroup late-window aHR (32d) paired with the subgroup
    standardized 24-month RD (30o), to show the relative-vs-absolute divergence.
Outputs: Figure_4D_absolute_rd.png and Figure_4C_subgroup_aHR_vs_RD.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.expanduser("~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper/ITT_Analysis/results")

bym = pd.read_csv(os.path.join(R, "target_trial_defnB_absolute_rd.csv"))
sg_rd = pd.read_csv(os.path.join(R, "target_trial_subgroup_absolute_rd.csv"))
sg_hr = pd.read_csv(os.path.join(R, "target_trial_subgroup_interactions_grace_mi.csv"))
sg_hr = sg_hr[(sg_hr.model == "late") & (sg_hr.cap == 2)]

# ---------- Panel D: by-month LATE-WINDOW (6-24mo) RD ----------
rd = bym[bym.metric == "rd_late"].copy()
rd["m"] = rd.Trial_Month.str.extract(r"(\d+)").astype(int)
rd = rd.sort_values("m")
fig, ax = plt.subplots(figsize=(5.2, 4.0))
ax.axhline(0, color="#999", lw=1, ls="--")
ax.errorbar(rd.m, rd.est*100, yerr=[(rd.est-rd.lo)*100, (rd.hi-rd.est)*100],
            fmt="o", color="#b2182b", ms=7, capsize=4, lw=1.6)
ax.set_xlabel("Month of disengagement (m)")
ax.set_ylabel("Late-window (6–24 mo) mortality\nrisk difference (percentage points)")
ax.set_title("D  Late-window absolute risk difference by month of disengagement", loc="left", fontsize=10, weight="bold")
ax.set_xticks(range(1, 7))
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(R, "Figure_4D_absolute_rd.png"), dpi=300)
plt.close(fig)

# ---------- Panel C right: subgroup aHR vs RD ----------
order = [("age_group", "15-24", "Age 15–24"), ("age_group", "25-44", "Age 25–44"),
         ("age_group", "45-64", "Age 45–64"), ("age_group", "≥65", "Age ≥65"),
         ("sex", "Female", "Female"), ("sex", "Male", "Male"),
         ("hiv_aids", "Negative", "HIV-negative"), ("hiv_aids", "Positive", "HIV-positive"),
         ("homelessness", "No", "Housed"), ("homelessness", "Yes", "Homeless")]
rows = []
for sub, lvl, lab in order:
    h = sg_hr[(sg_hr.Subgroup == sub) & (sg_hr.Level == lvl)]
    r = sg_rd[(sg_rd.Subgroup == sub) & (sg_rd.Level == lvl) & (sg_rd.metric == "rd_late")]
    if len(h) and len(r):
        rows.append(dict(lab=lab, hr=h.HR.values[0], hl=h.CI_L.values[0], hh=h.CI_H.values[0],
                         rd=r.est.values[0]*100, rl=r.lo.values[0]*100, rh=r.hi.values[0]*100))
df = pd.DataFrame(rows)
y = np.arange(len(df))[::-1]
fig, (axh, axr) = plt.subplots(1, 2, figsize=(8.4, 4.6), sharey=True,
                               gridspec_kw=dict(wspace=0.05))
# left: aHR (log scale)
axh.axvline(1, color="#999", lw=1, ls="--")
axh.errorbar(df.hr, y, xerr=[df.hr-df.hl, df.hh-df.hr], fmt="s", color="#2166ac", ms=6, capsize=3, lw=1.4)
axh.set_xscale("log"); axh.set_xticks([0.5, 1, 2, 4]); axh.set_xticklabels(["0.5", "1", "2", "4"])
axh.set_yticks(y); axh.set_yticklabels(df.lab)
axh.set_xlabel("Late-window aHR (relative)")
axh.set_title("C  Relative vs absolute, by subgroup", loc="left", fontsize=11, weight="bold")
axh.spines[["top", "right"]].set_visible(False)
# right: RD
axr.axvline(0, color="#999", lw=1, ls="--")
axr.errorbar(df.rd, y, xerr=[df.rd-df.rl, df.rh-df.rd], fmt="o", color="#b2182b", ms=6, capsize=3, lw=1.4)
axr.set_xlabel("Late-window (6–24 mo) risk difference (pp)")
axr.spines[["top", "left", "right"]].set_visible(False)
axr.tick_params(left=False)
fig.tight_layout()
fig.savefig(os.path.join(R, "Figure_4C_subgroup_aHR_vs_RD.png"), dpi=300)
plt.close(fig)
print("wrote Figure_4D_absolute_rd.png and Figure_4C_subgroup_aHR_vs_RD.png")
