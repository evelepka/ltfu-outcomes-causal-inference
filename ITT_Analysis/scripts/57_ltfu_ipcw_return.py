"""
57. IPCW for informative censoring at return-to-care
=====================================================
Within the LTFU cohort, estimate the marginal late-window (6-24 mo from LTFU)
mortality under a hypothetical "no return to care" scenario, using
inverse-probability-of-censoring weighting where return is the informative
censoring event.

Reports three estimates side-by-side:
  (1) ITT-like crude mortality       — observed deaths, no return censoring
  (2) Naive return-censored mortality — return treated as censor, NO weighting
  (3) IPCW-adjusted mortality        — return treated as censor, IPCW

The contrast (1) vs (3) quantifies how much the ITT-like LTFU mortality is
driven by post-LTFU returners (per the within-LTFU stratified finding,
returners die more, especially of TB). The IPCW estimate is the cleanest
target for "persistent disengagement: what if no one returned?"

Weights are computed from a Cox model for the return-to-care hazard fit on
baseline covariates. We use stabilised weights and truncate at the 1st/99th
percentile to limit influence of extreme values.

Outputs:
  • ITT_Analysis/results/ltfu_ipcw_return.csv      — table of estimates
  • ITT_Analysis/results/ltfu_ipcw_return_km.png   — three KM curves overlaid
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'ltfu_ipcw_return.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'ltfu_ipcw_return_km.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'ltfu_ipcw_return_km.pdf'

LANDMARK_YR = 0.5
MAX_FU_YR   = 2.0

# ---------- LTFU cohort ----------
d = pd.read_csv(COHORT, low_memory=False)
ltfu = d[d['itt_group']=='Loss to follow-up'].copy()
print(f'LTFU N = {len(ltfu):,}')

ADJ = ['age_tb','sex','hiv_aids','homelessness','alcohol','drug_use',
       'hosp_admission','diabetes','clinical_clean','resistance_clean',
       'lab_confirmed_stat','dot_status','incarcerated','race_clean','edu_clean',
       'diagnosis_setting']

# Use the FULL LTFU cohort, time origin = LTFU date (event/return times already from LTFU date).
# We report cumulative mortality at 12, 18, 24 mo. No landmark restriction.
lm = ltfu.copy()
print(f'Full LTFU cohort N = {len(lm):,} (time origin = LTFU date)')

# Cap follow-up at 24 mo
lm['t_d']   = np.minimum(lm['time_d'], MAX_FU_YR)
lm['ev_d']  = ((lm['event_d']==1) & (lm['time_d']<=MAX_FU_YR)).astype(int)
lm['t_rn_capped'] = np.where((lm['event_rn']==1) & (lm['time_rn']<=MAX_FU_YR), lm['time_rn'], np.nan)
lm['returned_in_2y'] = lm['t_rn_capped'].notna().astype(int)

print('\nEvents within 24 mo from LTFU:')
print(f'  deaths (any return status): {lm.ev_d.sum():,}')
print(f'  returns (any death status): {lm.returned_in_2y.sum():,}')
print(f'  death-before-return:        {((lm.ev_d==1) & ((lm.returned_in_2y==0) | (lm.t_rn_capped>lm.t_d))).sum():,}')

# ---------- Fit return-to-care Cox model on baseline covariates ----------
# Event = return within 2 yr; duration = time_rn if returned, else min(time_d, 2yr).
ret_df = lm.copy()
ret_df['t_ret']  = np.where(ret_df['returned_in_2y']==1, ret_df['t_rn_capped'], ret_df['t_d'])
ret_df['ev_ret'] = ret_df['returned_in_2y']
ret_df = ret_df[ret_df['t_ret']>0].copy()
adj_df = ret_df[['t_ret','ev_ret']+ADJ].copy()
adj_df = pd.get_dummies(adj_df, columns=[c for c in ADJ if c!='age_tb'], drop_first=True, dummy_na=False)
adj_df = adj_df.dropna().astype({c: float for c in adj_df.columns if adj_df[c].dtype==bool})
print(f'\nReturn-model fitting N = {len(adj_df):,} (after dropping rows with missing covariates)')
cox_ret = CoxPHFitter(penalizer=0.001)
cox_ret.fit(adj_df, duration_col='t_ret', event_col='ev_ret')
print(f'  Return Cox concordance = {cox_ret.concordance_index_:.3f}')

# ---------- Compute IPCW weights ----------
# For each subject in adj_df: weight = 1 / S_ret(observed_duration | X).
# Stabilised version: numerator = marginal S_ret(observed_duration).
sf_func  = cox_ret.predict_survival_function
sf_per   = sf_func(adj_df.drop(columns=['t_ret','ev_ret']))   # rows = times, cols = subjects
# Marginal (unconditional) KM for return as numerator
km_ret_marg = KaplanMeierFitter().fit(adj_df['t_ret'], adj_df['ev_ret'])
def s_marg(t):
    s=km_ret_marg.survival_function_
    idx=s.index.get_indexer([t],method='ffill')[0]; idx=max(idx,0)
    return float(s.iloc[idx,0])
def s_cond(col, t):
    idx=sf_per.index.get_indexer([t],method='ffill')[0]; idx=max(idx,0)
    return float(sf_per.iloc[idx][col])

w_list=[]
# To make this efficient, evaluate at each unique observed duration once
for col,(t,ev) in zip(sf_per.columns, zip(adj_df['t_ret'].values, adj_df['ev_ret'].values)):
    s_c = s_cond(col, t)
    s_m = s_marg(t)
    w   = s_m / s_c if s_c>0 else np.nan
    w_list.append(w)
adj_df['w_ipcw'] = w_list
# Truncate at 1st/99th percentile
q1, q99 = adj_df['w_ipcw'].quantile([0.01, 0.99])
adj_df['w_trunc'] = adj_df['w_ipcw'].clip(lower=q1, upper=q99)
print(f'  IPCW weight summary (truncated): min={adj_df.w_trunc.min():.2f}, '
      f'median={adj_df.w_trunc.median():.2f}, max={adj_df.w_trunc.max():.2f}')

# Attach weights back to the mortality follow-up rows
lm_w = lm.loc[adj_df.index].copy()
lm_w['w_trunc'] = adj_df['w_trunc'].values
# Mortality follow-up under return-censoring:
#   t_mort = min(time_d capped at 2y, time_rn if returned-before-death within 2y)
#   ev_d_under_censor = death observed BEFORE return-to-care (or no return)
returned_before_death = (lm_w['returned_in_2y']==1) & (lm_w['t_rn_capped'] < lm_w['t_d'])
lm_w['t_mort'] = np.where(returned_before_death, lm_w['t_rn_capped'], lm_w['t_d'])
lm_w['ev_d_under_censor'] = ((lm_w['ev_d']==1) & ~returned_before_death).astype(int)

# ---------- Three cumulative-mortality estimates ----------
def km_estimate(df, t_col, ev_col, label, weights=None):
    kmf = KaplanMeierFitter()
    if weights is None:
        kmf.fit(df[t_col], df[ev_col], label=label)
    else:
        kmf.fit(df[t_col], df[ev_col], label=label, weights=df[weights])
    return kmf

# (1) ITT-like crude: no return censoring, observed mortality from LTFU date.
itt = lm_w.copy()
kmf_itt = km_estimate(itt[itt['t_d']>0], 't_d','ev_d','(1) ITT (no return censoring)')

# (2) Naive return-censored: censor at return, no IPCW
nv = lm_w[lm_w['t_mort']>0].copy()
kmf_nv = km_estimate(nv, 't_mort', 'ev_d_under_censor', '(2) Naive return-censored')

# (3) IPCW return-censored: same as (2) but weighted
kmf_ipcw = km_estimate(nv, 't_mort', 'ev_d_under_censor',
                       '(3) IPCW (return-censored + weighted)', weights='w_trunc')

def cum_at(kmf, t_yr_from_ltfu):
    t = t_yr_from_ltfu  # time scale is years from LTFU date directly
    s_f = kmf.survival_function_
    ci  = kmf.confidence_interval_
    idx = s_f.index.get_indexer([t], method='ffill')[0]; idx=max(idx,0)
    s   = float(s_f.iloc[idx,0])
    try:
        lo=float(ci.iloc[idx,0]); hi=float(ci.iloc[idx,1])
        return 1-s, (1-hi, 1-lo)
    except Exception:
        return 1-s, (np.nan, np.nan)

rows=[]
for t_mo in (12, 18, 24):
    for label, kmf, n in [
        ('(1) ITT (no return censoring)',         kmf_itt,  len(itt)),
        ('(2) Naive return-censored',             kmf_nv,   len(nv)),
        ('(3) IPCW return-censored + weighted',   kmf_ipcw, len(nv)),
    ]:
        m,(lo,hi) = cum_at(kmf, t_mo/12.0)
        rows.append({'metric':f'Cumulative mortality at {t_mo} mo from LTFU',
                     'estimate': label, 'n': n,
                     'pct': round(100*m,2),
                     'ci_lo_pct': round(100*lo,2) if not np.isnan(lo) else None,
                     'ci_hi_pct': round(100*hi,2) if not np.isnan(hi) else None})

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}\n')

# print summary
for r in rows:
    print(f'  {r["metric"]:<50} {r["estimate"]:<45} '
          f'{r["pct"]:>5.2f}%  ({r["ci_lo_pct"]:.2f}–{r["ci_hi_pct"]:.2f})')

# ---------- Plot: three overlaid KM curves ----------
fig, ax = plt.subplots(figsize=(9,5.6))
colors = {'itt':'#d62728', 'nv':'#888888', 'ipcw':'#1f77b4'}
def plot(kmf, color, label):
    s = kmf.survival_function_; ci = kmf.confidence_interval_
    x = s.index.values * 12   # time scale already years from LTFU date
    y = (1-s.iloc[:,0].values)*100
    ax.step(x, y, where='post', color=color, lw=2.2, label=label)
    try:
        ax.fill_between(x, (1-ci.iloc[:,1])*100, (1-ci.iloc[:,0])*100,
                        color=color, alpha=0.15, step='post')
    except Exception:
        pass
plot(kmf_itt,  colors['itt'],  '(1) ITT (no return censoring)')
plot(kmf_nv,   colors['nv'],   '(2) Naive return-censored (biased)')
plot(kmf_ipcw, colors['ipcw'], '(3) IPCW: persistent disengagement (target)')
ax.set_xlim(0, MAX_FU_YR*12)
ax.set_xlabel('Months from LTFU date')
ax.set_ylabel('Cumulative mortality (%)')
ax.set_title('Late-window mortality among LTFU patients:\n'
             'ITT vs. naive-censored vs. IPCW (persistent disengagement)',
             loc='left', fontsize=11, fontweight='bold')
ax.legend(loc='upper left', frameon=False)
ax.grid(alpha=0.25)
for s_ in ('top','right'): ax.spines[s_].set_visible(False)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'\nwrote {OUT_PNG}')
