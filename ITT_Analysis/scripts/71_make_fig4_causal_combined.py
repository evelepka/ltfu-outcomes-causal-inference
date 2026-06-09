"""
71_make_fig4_causal_combined.py  — single 4-panel Figure 4 (causal effect of LTFU).
  A  Adjusted time-varying HR over months (piecewise/step)
  B  By-month aHR (early & late windows)
  C  By-month absolute risk difference (early & late windows)
  D  Subgroup late-window aHR (left) + absolute RD (right)
Output: Figure_4_causal_combined.png  (no overall title; RD axes in percentage points)
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.expanduser("~/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper/ITT_Analysis/results")
red, blue, grey = "#b2182b", "#2166ac", "#999999"
red_l, blue_l = "#e8907f", "#7fb0d6"
months = list(range(1, 7))

HRT_INTERVAL = "2mo"   # "3mo", "2mo", or use Figure_3a_HR_over_time_defnB.csv for monthly
tv  = pd.read_csv(os.path.join(R, f"Figure_4a_HR_over_time_{HRT_INTERVAL}.csv"))
arr = pd.read_csv(os.path.join(R, "target_trial_defnB_mi_early_late_array.csv"))
rdm = pd.read_csv(os.path.join(R, "target_trial_defnB_bymonth_earlylate_rd.csv"))
hr  = pd.read_csv(os.path.join(R, "target_trial_subgroup_interactions_grace_mi.csv"))
hr  = hr[(hr.model == "late") & (hr.cap == 2)]
rds = pd.read_csv(os.path.join(R, "target_trial_subgroup_absolute_rd.csv"))

fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(2, 6, hspace=0.45, wspace=1.0)
axA = fig.add_subplot(gs[0, 0:2]); axB = fig.add_subplot(gs[0, 2:4]); axC = fig.add_subplot(gs[0, 4:6])
axDl = fig.add_subplot(gs[1, 0:3]); axDr = fig.add_subplot(gs[1, 3:6], sharey=axDl)

# A: time-varying HR as points with 95% CI (piecewise interval estimates)
mo = (tv.month_mid * 12).values
axA.axhline(1, color=grey, lw=1, ls="--")
axA.errorbar(mo, tv.HR, yerr=[tv.HR - tv.CI_L, tv.CI_H - tv.HR],
             fmt="o", color=blue, ms=5, capsize=2.5, lw=1.1, elinewidth=1.0)
axA.set_yscale("log"); axA.set_yticks([0.25,0.5,1,2,4]); axA.set_yticklabels(["0.25","0.5","1","2","4"])
axA.set_xlim(0,24); axA.set_xticks([0,6,12,18,24])
axA.set_xlabel("Months since disengagement"); axA.set_ylabel("Crude hazard ratio\n(LTFU vs on-treatment)")
axA.set_title("A  Crude time-varying HR", loc="left", fontsize=11, weight="bold")
axA.spines[["top","right"]].set_visible(False)

def by_series(ax, df_idx, color, lab, mk, dx, logy=True):
    hrs=[df_idx.loc[f"Month_{m}","est" if "est" in df_idx.columns else "HR"] for m in months]
    lo=[df_idx.loc[f"Month_{m}","lo" if "lo" in df_idx.columns else "CI_L"] for m in months]
    hi=[df_idx.loc[f"Month_{m}","hi" if "hi" in df_idx.columns else "CI_H"] for m in months]
    sc = 100 if not logy else 1
    ax.errorbar(np.array(months)+dx,[h*sc for h in hrs],
                yerr=[[(h-l)*sc for h,l in zip(hrs,lo)],[(h2-h)*sc for h,h2 in zip(hrs,hi)]],
                fmt=mk,color=color,ms=6,capsize=3,lw=1.3,label=lab)

# B: by-month aHR early & late
early=arr[arr.model=="early"].set_index("Trial_Month"); late=arr[(arr.model=="late")&(arr.cap==2)].set_index("Trial_Month")
axB.axhline(1,color=grey,lw=1,ls="--")
by_series(axB,early,blue_l,"Early (0–6 mo)","o",-0.12); by_series(axB,late,blue,"Late (6–24 mo)","s",0.12)
axB.set_yscale("log"); axB.set_yticks([0.5,1,2,4]); axB.set_yticklabels(["0.5","1","2","4"])
axB.set_xticks(months); axB.set_xticklabels([f"M{m}" for m in months])
axB.set_xlabel("Month of disengagement"); axB.set_ylabel("Adjusted hazard ratio")
axB.set_title("B  Relative effect by month", loc="left", fontsize=11, weight="bold")
axB.legend(frameon=False, fontsize=8); axB.spines[["top","right"]].set_visible(False)

# C: by-month RD early & late
ee=rdm[rdm.metric=="rd_early"].set_index("Trial_Month"); ll=rdm[rdm.metric=="rd_late"].set_index("Trial_Month")
axC.axhline(0,color=grey,lw=1,ls="--")
by_series(axC,ee,red_l,"Early (0–6 mo)","o",-0.12,logy=False); by_series(axC,ll,red,"Late (6–24 mo)","s",0.12,logy=False)
axC.set_xticks(months); axC.set_xticklabels([f"M{m}" for m in months])
axC.set_xlabel("Month of disengagement"); axC.set_ylabel("Risk difference\n(percentage points)")
axC.set_title("C  Absolute effect by month", loc="left", fontsize=11, weight="bold")
axC.legend(frameon=False, fontsize=8); axC.spines[["top","right"]].set_visible(False)

# D: subgroup aHR (left) + RD (right)
order=[("age_group","15-24","Age 15–24"),("age_group","25-44","Age 25–44"),("age_group","45-64","Age 45–64"),
       ("age_group","≥65","Age ≥65"),("sex","Female","Female"),("sex","Male","Male"),
       ("hiv_aids","Negative","HIV-negative"),("hiv_aids","Positive","HIV-positive"),
       ("homelessness","No","Housed"),("homelessness","Yes","Homeless")]
rdsl=rds[rds.metric=="rd_late"].set_index(["Subgroup","Level"])
rows=[]
for s,l,lab in order:
    h=hr[(hr.Subgroup==s)&(hr.Level==l)]; r=rdsl.loc[(s,l)]
    rows.append((lab,h.HR.values[0],h.CI_L.values[0],h.CI_H.values[0],r.est*100,r.lo*100,r.hi*100))
y=np.arange(len(rows))[::-1]
axDl.axvline(1,color=grey,lw=1,ls="--")
axDl.errorbar([r[1] for r in rows],y,xerr=[[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]],fmt="s",color=blue,ms=6,capsize=3,lw=1.3)
axDl.set_xscale("log"); axDl.set_xticks([0.5,1,2,4]); axDl.set_xticklabels(["0.5","1","2","4"])
axDl.set_yticks(y); axDl.set_yticklabels([r[0] for r in rows])
axDl.set_xlabel("Late-window adjusted hazard ratio (relative)")
axDl.set_title("D  Subgroup effect: relative (left) and absolute (right)", loc="left", fontsize=11, weight="bold")
axDl.spines[["top","right"]].set_visible(False)
axDr.axvline(0,color=grey,lw=1,ls="--")
axDr.errorbar([r[4] for r in rows],y,xerr=[[r[4]-r[5] for r in rows],[r[6]-r[4] for r in rows]],fmt="o",color=red,ms=6,capsize=3,lw=1.3)
axDr.set_xlabel("Late-window risk difference (percentage points)")
axDr.spines[["top","left","right"]].set_visible(False); axDr.tick_params(left=False)

out=os.path.join(R,"Figure_4_causal_combined.png")
fig.savefig(out,dpi=300,bbox_inches="tight"); print("wrote",out)
