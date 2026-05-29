"""
57b. LTFU vs on-treatment HR under ITT and IPCW (return-censoring)
====================================================================
Compare LTFU vs on-treatment 24-month mortality at the m=3 landmark
(middle of treatment; treatment_months in (2, 3] for LTFU, > 3 for on-tx).

Two estimates:
  • ITT-like      : observed mortality, no return censoring.
  • IPCW          : LTFU arm censored at return-to-care, weighted by IPCW;
                     on-tx arm unweighted.

The contrast between the two HRs quantifies how much of the LTFU vs on-tx
late-window effect is mediated through return-to-care.

Outputs:
  • ITT_Analysis/results/ipcw_ltfu_vs_ontx_m3.csv
  • ITT_Analysis/results/ipcw_ltfu_vs_ontx_m3_km.png/.pdf
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'ipcw_ltfu_vs_ontx_m3.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'ipcw_ltfu_vs_ontx_m3_km.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'ipcw_ltfu_vs_ontx_m3_km.pdf'
OUT_DATA_ITT  = BASE/'ITT_Analysis'/'data'/'analytic_itt_m3.csv'
OUT_DATA_IPCW = BASE/'ITT_Analysis'/'data'/'analytic_ipcw_m3.csv'

M_LANDMARK_YR = 3.0/12.0   # m=3 landmark
MAX_FU_YR     = 2.0

ADJ = ['age_tb','sex','hiv_aids','homelessness','alcohol','drug_use',
       'hosp_admission','diabetes','clinical_clean','resistance_clean',
       'lab_confirmed_stat','dot_status','incarcerated','race_clean','edu_clean',
       'diagnosis_setting']

# ---------- cohort ----------
d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'], errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],   errors='coerce')
d['tx_yrs']    =(d['end_date']-d['best_start']).dt.days/365.25

# LTFU arm at m=3: LTFU patients with tx_yrs in (2/12, 3/12]
# Their time origin = LTFU date (= end_date). time_d already from LTFU date.
# To align with m=3 landmark for the on-tx arm, treat their "follow-up time" simply
# as time_d (since LTFU happens close to m=3 anyway).
ltfu_m3 = d[(d['itt_group']=='Loss to follow-up')
            & (d['tx_yrs']>2.0/12.0) & (d['tx_yrs']<=3.0/12.0)].copy()
ltfu_m3['arm'] = 'LTFU'
print(f'LTFU arm at m=3: N={len(ltfu_m3):,}')

# Cap follow-up at 24 mo from LTFU date
ltfu_m3['t_d']   = np.minimum(ltfu_m3['time_d'], MAX_FU_YR)
ltfu_m3['ev_d']  = ((ltfu_m3['event_d']==1) & (ltfu_m3['time_d']<=MAX_FU_YR)).astype(int)
ltfu_m3['t_rn_capped'] = np.where((ltfu_m3['event_rn']==1) & (ltfu_m3['time_rn']<=MAX_FU_YR),
                                   ltfu_m3['time_rn'], np.nan)
ltfu_m3['returned_in_2y'] = ltfu_m3['t_rn_capped'].notna().astype(int)

# On-tx arm at m=3: tx_yrs > 3/12 AND was non-LTFU (continuing or completing treatment)
# Time origin = best_start + 3 months. Follow-up = (time_d_tx - 3/12), capped at 2 yr.
ontx = d[(d['itt_group']!='Loss to follow-up') & (d['tx_yrs']>3.0/12.0)].copy()
ontx['arm'] = 'OnTx'
# time_d_tx is years from best_start to event/censor
ontx['t_from_lm']    = ontx['time_d_tx'] - M_LANDMARK_YR
ontx = ontx[ontx['t_from_lm']>0]   # alive past m=3 landmark
ontx['t_d']   = np.minimum(ontx['t_from_lm'], MAX_FU_YR)
ontx['ev_d']  = ((ontx['event_d']==1) & (ontx['t_from_lm']<=MAX_FU_YR)).astype(int)
print(f'On-tx arm at m=3: N={len(ontx):,}')

# ---------- IPCW weights for LTFU arm: return-to-care censoring ----------
# Fit Cox model for return event among LTFU arm patients, then compute weights.
ret = ltfu_m3.copy()
ret['t_ret']  = np.where(ret['returned_in_2y']==1, ret['t_rn_capped'], ret['t_d'])
ret['ev_ret'] = ret['returned_in_2y']
ret = ret[ret['t_ret']>0]
ret_x = pd.get_dummies(ret[['t_ret','ev_ret']+ADJ],
                       columns=[c for c in ADJ if c!='age_tb'], drop_first=True, dummy_na=False)
ret_x = ret_x.dropna().astype({c:float for c in ret_x.columns if ret_x[c].dtype==bool})
print(f'Return-model fit N = {len(ret_x):,}')
cox_ret = CoxPHFitter(penalizer=0.001).fit(ret_x, duration_col='t_ret', event_col='ev_ret')
print(f'  return Cox concordance = {cox_ret.concordance_index_:.3f}')
# Stabilised weight at the LTFU subject's mortality-censoring time
sf_per = cox_ret.predict_survival_function(ret_x.drop(columns=['t_ret','ev_ret']))
km_marg = KaplanMeierFitter().fit(ret_x['t_ret'], ret_x['ev_ret'])
def s_marg(t):
    s=km_marg.survival_function_
    idx=s.index.get_indexer([t],method='ffill')[0]; idx=max(idx,0)
    return float(s.iloc[idx,0])
def s_cond(col, t):
    idx=sf_per.index.get_indexer([t],method='ffill')[0]; idx=max(idx,0)
    return float(sf_per.iloc[idx][col])
# mortality follow-up endpoint = min(t_d, t_rn if returned before death)
returned_before_death = (ret['returned_in_2y']==1) & (ret['t_rn_capped'] < ret['t_d'])
ret['t_mort'] = np.where(returned_before_death, ret['t_rn_capped'], ret['t_d'])
ret['ev_mort_under_censor'] = ((ret['ev_d']==1) & ~returned_before_death).astype(int)
ws=[]
for col, t in zip(sf_per.columns, ret['t_mort'].values):
    s_c=s_cond(col,t); s_m=s_marg(t)
    ws.append(s_m/s_c if s_c>0 else np.nan)
ret['w_ipcw'] = ws
q1,q99 = ret['w_ipcw'].quantile([0.01,0.99])
ret['w_trunc'] = ret['w_ipcw'].clip(lower=q1, upper=q99)
print(f'  IPCW weight (truncated): min={ret.w_trunc.min():.2f}, '
      f'median={ret.w_trunc.median():.2f}, max={ret.w_trunc.max():.2f}')

# ---------- Build pooled analytic dataset ----------
COMMON_COLS = ['arm','age_tb','sex','hiv_aids','homelessness','alcohol','drug_use',
               'hosp_admission','diabetes','clinical_clean','resistance_clean',
               'lab_confirmed_stat','dot_status','incarcerated','race_clean','edu_clean',
               'diagnosis_setting']

# (1) ITT pooled: LTFU uses observed t_d/ev_d (no return censoring)
itt_l = ltfu_m3[['t_d','ev_d']+COMMON_COLS].copy(); itt_l['w']=1.0
itt_o = ontx[['t_d','ev_d']+COMMON_COLS].copy();   itt_o['w']=1.0
itt = pd.concat([itt_l, itt_o], ignore_index=True)
itt['ltfu'] = (itt['arm']=='LTFU').astype(int)
itt = itt.drop(columns=['arm'])

# (3) IPCW pooled: LTFU uses t_mort/ev_mort_under_censor + IPCW weight; on-tx unweighted
ipw_l = ret[['t_mort','ev_mort_under_censor','w_trunc']+COMMON_COLS].rename(
    columns={'t_mort':'t_d','ev_mort_under_censor':'ev_d','w_trunc':'w'}).copy()
ipw_o = ontx[['t_d','ev_d']+COMMON_COLS].copy(); ipw_o['w']=1.0
ipw = pd.concat([ipw_l, ipw_o], ignore_index=True)
ipw['ltfu'] = (ipw['arm']=='LTFU').astype(int)
ipw = ipw.drop(columns=['arm'])

# ---------- Cox HRs ----------
def fit_cox(df, weight=False):
    df = df[df['t_d']>0].copy()
    X = pd.get_dummies(df, columns=[c for c in COMMON_COLS if c!='age_tb' and c!='arm'],
                        drop_first=True, dummy_na=False)
    X = X.dropna().astype({c:float for c in X.columns if X[c].dtype==bool})
    cph = CoxPHFitter(penalizer=0.001)
    kwargs = dict(duration_col='t_d', event_col='ev_d')
    if weight:
        kwargs['weights_col']='w'
        cph.fit(X, robust=True, **kwargs)
    else:
        X = X.drop(columns=['w'])
        cph.fit(X, **kwargs)
    hr = float(cph.hazard_ratios_['ltfu'])
    ci = cph.confidence_intervals_.loc['ltfu'].tolist()
    return hr, float(np.exp(ci[0])), float(np.exp(ci[1])), len(X)

# Unadjusted (crude) Cox for ITT and IPCW
def fit_crude(df, weight=False):
    df = df[df['t_d']>0].copy()
    cols = ['t_d','ev_d','ltfu'] + (['w'] if weight else [])
    sub = df[cols].copy()
    cph = CoxPHFitter()
    if weight:
        cph.fit(sub, duration_col='t_d', event_col='ev_d', weights_col='w', robust=True)
    else:
        cph.fit(sub, duration_col='t_d', event_col='ev_d')
    hr=float(cph.hazard_ratios_['ltfu'])
    ci=cph.confidence_intervals_.loc['ltfu'].tolist()
    return hr, float(np.exp(ci[0])), float(np.exp(ci[1])), len(sub)

rows=[]
for label, df, weight in [('ITT (no return censoring)',  itt, False),
                          ('IPCW (return-censored + weighted)', ipw, True)]:
    hr,lo,hi,n = fit_crude(df, weight=weight)
    print(f'\n[crude]    {label}:   HR={hr:.2f} ({lo:.2f}–{hi:.2f}),  N={n:,}')
    rows.append({'spec':'crude HR','estimate':label,
                 'HR':round(hr,3),'lo':round(lo,3),'hi':round(hi,3),'N':n})
    hr,lo,hi,n = fit_cox(df, weight=weight)
    print(f'[adjusted] {label}:   aHR={hr:.2f} ({lo:.2f}–{hi:.2f}),  N={n:,}')
    rows.append({'spec':'adjusted aHR','estimate':label,
                 'HR':round(hr,3),'lo':round(lo,3),'hi':round(hi,3),'N':n})

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')

# Save analytic datasets for downstream Bayesian fit in R
OUT_DATA_ITT.parent.mkdir(parents=True, exist_ok=True)
itt.to_csv(OUT_DATA_ITT, index=False)
ipw.to_csv(OUT_DATA_IPCW, index=False)
print(f'wrote {OUT_DATA_ITT} ({len(itt):,} rows)')
print(f'wrote {OUT_DATA_IPCW} ({len(ipw):,} rows)')

# ---------- KM curves (crude, unweighted) for visualisation ----------
def km(df, label, weights=None):
    k = KaplanMeierFitter()
    if weights is None:
        k.fit(df['t_d'], df['ev_d'], label=label)
    else:
        k.fit(df['t_d'], df['ev_d'], label=label, weights=df['w'])
    return k

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
def plot(ax, kmf_l, kmf_o, title):
    for kmf, color, lab in [(kmf_l,'#1f77b4', kmf_l._label),
                            (kmf_o,'#888888', kmf_o._label)]:
        s=kmf.survival_function_; ci=kmf.confidence_interval_
        x=s.index.values*12; y=(1-s.iloc[:,0].values)*100
        ax.step(x, y, where='post', color=color, lw=2.0, label=lab)
        try:
            ax.fill_between(x,(1-ci.iloc[:,1])*100,(1-ci.iloc[:,0])*100,
                            color=color, alpha=0.15, step='post')
        except Exception: pass
    ax.set_xlim(0, MAX_FU_YR*12); ax.set_xlabel('Months since m=3 landmark')
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.grid(alpha=0.25); ax.legend(loc='upper left', frameon=False, fontsize=9)
    for s_ in ('top','right'): ax.spines[s_].set_visible(False)

# Plot 1: ITT
l_itt = km(itt_l.rename(columns={'t_d':'t_d','ev_d':'ev_d'}).assign(label=None),
           f'LTFU at m=3 (n={len(ltfu_m3):,})')
o_itt = km(itt_o, f'On-treatment at m=3 (n={len(ontx):,})')
plot(axes[0], l_itt, o_itt, 'A. ITT (no return censoring)')

# Plot 2: IPCW
l_ipw = km(ret.assign(t_d=ret['t_mort'], ev_d=ret['ev_mort_under_censor'], w=ret['w_trunc']),
           f'LTFU at m=3, IPCW (n={len(ret):,})', weights='w')
o_ipw = km(ontx, f'On-treatment at m=3 (n={len(ontx):,})')
plot(axes[1], l_ipw, o_ipw, 'B. IPCW (LTFU return-censored + weighted)')

axes[0].set_ylabel('Cumulative mortality (%)')
fig.suptitle('LTFU vs on-treatment mortality at m=3 landmark (24-month follow-up)',
             fontsize=12.5, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'wrote {OUT_PNG}')
