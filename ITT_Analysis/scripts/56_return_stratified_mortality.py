"""
56. Late-window mortality among LTFU patients, stratified by 6-month return status
==================================================================================
Within the LTFU cohort (N = 21,619), apply a 6-month landmark from the LTFU date
to avoid immortal-time bias, then compare late-window (6–24 months from LTFU)
mortality between:
  • Returners by 6 months   — re-entered TB care within 6 months of LTFU
  • Non-returners by 6 months — no recorded re-entry by 6 months

Outputs:
  • ITT_Analysis/results/return_stratified_mortality.csv  (KM estimates + Cox HR)
  • ITT_Analysis/results/return_stratified_km.png/pdf      (KM curve plot)
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'return_stratified_mortality.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'return_stratified_km.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'return_stratified_km.pdf'

LANDMARK_YR = 0.5     # 6 months
MAX_FU_YR   = 2.0     # 24 months from LTFU

# ----- load cohort -----
d = pd.read_csv(COHORT, low_memory=False)
ltfu = d[d['itt_group']=='Loss to follow-up'].copy()
print(f'LTFU cohort N = {len(ltfu):,}')

# Returner by 6 months: re-notification event within 0.5 yr of LTFU
ltfu['returner_6mo'] = ((ltfu['event_rn']==1) & (ltfu['time_rn'] <= LANDMARK_YR)).astype(int)
print('Returner-by-6mo counts:'); print(ltfu['returner_6mo'].value_counts())

# Landmark eligibility: must be alive at LTFU + 6mo
#   - If died before LANDMARK_YR (event_d==1 and time_d <= LANDMARK_YR), excluded.
ltfu['died_before_landmark'] = ((ltfu['event_d']==1) & (ltfu['time_d'] <= LANDMARK_YR)).astype(int)
print(f'Died before 6-mo landmark, excluded from landmark cohort: '
      f'{ltfu["died_before_landmark"].sum():,}')

lm = ltfu[ltfu['died_before_landmark']==0].copy()
print(f'Landmark-eligible (alive at 6 mo) N = {len(lm):,}')
print('  Among landmark-eligible — returner counts:'); print(lm['returner_6mo'].value_counts())

# Late-window follow-up: from 6 mo to 24 mo, time_d already from LTFU date
lm['fu_yr']    = np.minimum(lm['time_d'], MAX_FU_YR) - LANDMARK_YR
lm['fu_event'] = ((lm['event_d']==1) & (lm['time_d'] <= MAX_FU_YR)).astype(int)
lm = lm[lm['fu_yr']>0].copy()   # drop zero-duration

# ----- KM estimates -----
kmf_R, kmf_N = KaplanMeierFitter(), KaplanMeierFitter()
R = lm[lm['returner_6mo']==1]
N = lm[lm['returner_6mo']==0]
kmf_R.fit(R['fu_yr'], R['fu_event'], label='Returner by 6mo')
kmf_N.fit(N['fu_yr'], N['fu_event'], label='Non-returner by 6mo')

# Late-window cumulative mortality at 6, 12, 18, 24 mo from LTFU
def cum_mortality(kmf, t_yr_from_ltfu):
    # KM is on fu_yr (since landmark); cumulative incidence = 1 - S(t - 0.5)
    t = t_yr_from_ltfu - LANDMARK_YR
    if t<=0: return 0.0, (0.0, 0.0)
    s_func = kmf.survival_function_
    ci_func = kmf.confidence_interval_
    idx = s_func.index.get_indexer([t], method='ffill')[0]
    if idx == -1: idx = 0
    s = float(s_func.iloc[idx, 0])
    lo = float(ci_func.iloc[idx, 0]); hi = float(ci_func.iloc[idx, 1])
    return 1-s, (1-hi, 1-lo)

rows=[]
for t_mo in (12, 18, 24):
    t_yr = t_mo/12.0
    for arm, kmf, group in [('Returner by 6mo', kmf_R, R), ('Non-returner by 6mo', kmf_N, N)]:
        m, (lo, hi) = cum_mortality(kmf, t_yr)
        rows.append({'metric':f'Cumulative mortality at {t_mo} mo from LTFU',
                     'arm': arm, 'n': len(group),
                     'estimate_pct': round(100*m, 2),
                     'ci_low_pct':   round(100*lo, 2),
                     'ci_high_pct':  round(100*hi, 2)})

# log-rank
lr = logrank_test(R['fu_yr'], N['fu_yr'], R['fu_event'], N['fu_event'])
print(f'\nlog-rank test p = {lr.p_value:.3e}')

# Cox HR (returner vs non-returner)
cox_df = lm[['fu_yr','fu_event','returner_6mo']].copy()
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='fu_yr', event_col='fu_event')
hr   = float(cph.hazard_ratios_['returner_6mo'])
ci   = cph.confidence_intervals_.loc['returner_6mo'].tolist()
hr_lo, hr_hi = np.exp(ci[0]), np.exp(ci[1])
p    = float(cph.summary.loc['returner_6mo','p'])
print(f'Crude Cox HR (returner vs non-returner) = {hr:.2f} (95% CI {hr_lo:.2f}–{hr_hi:.2f}), p={p:.3e}')
rows.append({'metric':'Crude Cox HR (returner vs non-returner), late window',
             'arm':'returner vs non-returner','n':len(lm),
             'estimate_pct': round(hr,3),
             'ci_low_pct':   round(hr_lo,3),
             'ci_high_pct':  round(hr_hi,3)})

# Adjusted Cox HR — control for baseline covariates available
adj_covars = ['age_tb','sex','hiv_aids','homelessness','alcohol','drug_use',
              'hosp_admission','diabetes','clinical_clean','resistance_clean',
              'lab_confirmed_stat','dot_status','incarcerated']
adj_df = lm[['fu_yr','fu_event','returner_6mo']+adj_covars].copy()
adj_df = pd.get_dummies(adj_df, columns=[c for c in adj_covars if c!='age_tb'],
                        drop_first=True, dummy_na=False)
adj_df = adj_df.dropna()
adj_df = adj_df.astype({c: float for c in adj_df.columns if adj_df[c].dtype==bool})
print(f'Adjusted Cox N = {len(adj_df):,}')
cph_adj = CoxPHFitter(penalizer=0.001)
cph_adj.fit(adj_df, duration_col='fu_yr', event_col='fu_event')
hr_a = float(cph_adj.hazard_ratios_['returner_6mo'])
ci_a = cph_adj.confidence_intervals_.loc['returner_6mo'].tolist()
hr_a_lo, hr_a_hi = np.exp(ci_a[0]), np.exp(ci_a[1])
p_a = float(cph_adj.summary.loc['returner_6mo','p'])
print(f'Adjusted Cox HR (returner vs non-returner) = {hr_a:.2f} (95% CI {hr_a_lo:.2f}–{hr_a_hi:.2f}), p={p_a:.3e}')
rows.append({'metric':'Adjusted Cox HR (returner vs non-returner), late window',
             'arm':'returner vs non-returner (adjusted)','n':len(adj_df),
             'estimate_pct': round(hr_a,3),
             'ci_low_pct':   round(hr_a_lo,3),
             'ci_high_pct':  round(hr_a_hi,3)})

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')

# ----- KM plot -----
fig, ax = plt.subplots(figsize=(8.5,5.5))
# Display the curves on a months-from-LTFU x-axis: add 6mo offset
def plot_arm(kmf, color, label):
    s = kmf.survival_function_
    ci = kmf.confidence_interval_
    x = (s.index.values + LANDMARK_YR)*12   # months from LTFU
    y_event = 1 - s.iloc[:,0].values
    ax.step(x, y_event*100, where='post', color=color, lw=2.0, label=label)
    ax.fill_between(x, (1-ci.iloc[:,1])*100, (1-ci.iloc[:,0])*100,
                    color=color, alpha=0.15, step='post')
plot_arm(kmf_R, '#2ca02c', f'Returner by 6 mo (n={len(R):,})')
plot_arm(kmf_N, '#d62728', f'Non-returner by 6 mo (n={len(N):,})')

ax.axvline(LANDMARK_YR*12, color='grey', ls=':', lw=0.8)
ax.text(LANDMARK_YR*12+0.2, 0.3, '6-mo landmark', fontsize=8.5, color='dimgrey')
ax.set_xlim(LANDMARK_YR*12, MAX_FU_YR*12)
ax.set_ylim(0, max(8, ax.get_ylim()[1]))
ax.set_xlabel('Months from LTFU date')
ax.set_ylabel('Cumulative mortality (%)')
ax.set_title('Late-window mortality among LTFU patients, by 6-month return status\n'
             '(6-month landmark; follow-up 6–24 months from LTFU)',
             loc='left', fontsize=11, fontweight='bold')
ax.legend(loc='upper left', frameon=False)
ax.grid(alpha=0.25)
for s in ('top','right'): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'wrote {OUT_PNG}')
