"""
58. Counterfactual mortality among returners: "what if they hadn't abandoned?"
==============================================================================
G-computation / standardisation estimand:
    E[Y(stay-on-tx) | LTFU=1, Returner=1]
i.e. the 24-month cumulative mortality returners would have experienced if
they had instead remained continuously on treatment.

Population:
  Returners-at-m=3 = LTFU at m=3 (tx_yrs in (2/12, 3/12]) AND returned to
  care within 6 months of LTFU.
Comparator:
  On-treatment-at-m=3 = patients still on treatment past 3 months of treatment.

Method:
  1. Fit Cox PH on on-tx-at-m=3 with baseline covariates X.
  2. For each returner, predict S(t | X_i, stay-on-tx) at t = 24 mo (since m=3).
  3. Counterfactual mortality = 1 - mean(S(24 | X_i)).
  4. Observed returner mortality at 24 mo from KM.
  5. Bootstrap (B=500) on returner subgroup for the counterfactual estimate
     to get a 95% CI.

Assumptions:
  - No unmeasured confounding between returners and on-tx given X
  - Correctly specified outcome model

Outputs:
  ITT_Analysis/results/counterfactual_returners.csv
"""
from pathlib import Path
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'counterfactual_returners.csv'

M_LM_YR = 3.0/12.0
MAX_FU_YR = 2.0

# Baseline covariates (parsimonious set that worked in the Bayesian fit)
COVARS = ['age_tb','sex','hiv_aids','hosp_admission','clinical_clean',
          'homelessness','alcohol','drug_use','dot_status']

# ---------- cohort ----------
d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'], errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],   errors='coerce')
d['tx_yrs']    =(d['end_date']-d['best_start']).dt.days/365.25

# LTFU at m=3 AND returned within 6 mo of LTFU → "returner-at-m=3"
ltfu_m3 = d[(d['itt_group']=='Loss to follow-up')
            & (d['tx_yrs']>2.0/12.0) & (d['tx_yrs']<=3.0/12.0)].copy()
returners = ltfu_m3[(ltfu_m3['event_rn']==1) & (ltfu_m3['time_rn']<=0.5)].copy()
print(f'LTFU at m=3 total:            {len(ltfu_m3):,}')
print(f'Returners by 6 mo at m=3:     {len(returners):,}')

# time_d for LTFU is from LTFU date ≈ m=3 landmark, so already on the m=3 timescale
returners['t']  = np.minimum(returners['time_d'], MAX_FU_YR)
returners['ev'] = ((returners['event_d']==1) & (returners['time_d']<=MAX_FU_YR)).astype(int)
print(f'  late-window deaths within 24 mo: {returners["ev"].sum()}')

# On-tx at m=3
ontx = d[(d['itt_group']!='Loss to follow-up') & (d['tx_yrs']>3.0/12.0)].copy()
ontx['t_from_lm'] = ontx['time_d_tx'] - M_LM_YR
ontx = ontx[ontx['t_from_lm']>0].copy()
ontx['t']  = np.minimum(ontx['t_from_lm'], MAX_FU_YR)
ontx['ev'] = ((ontx['event_d']==1) & (ontx['t_from_lm']<=MAX_FU_YR)).astype(int)
print(f'On-tx at m=3:                 {len(ontx):,} (events {ontx["ev"].sum()})')

# ---------- fit Cox on on-tx with baseline covariates ----------
ontx_x = ontx[['t','ev']+COVARS].copy()
ontx_x = pd.get_dummies(ontx_x, columns=[c for c in COVARS if c!='age_tb'],
                        drop_first=True, dummy_na=False)
ontx_x = ontx_x.dropna().astype({c: float for c in ontx_x.columns if ontx_x[c].dtype==bool})
cph = CoxPHFitter(penalizer=0.001).fit(ontx_x, duration_col='t', event_col='ev')
print(f'\nOn-tx Cox concordance: {cph.concordance_index_:.3f}')

# ---------- predict counterfactual for each returner ----------
ret_x = returners[['t','ev']+COVARS].copy()
ret_x = pd.get_dummies(ret_x, columns=[c for c in COVARS if c!='age_tb'],
                       drop_first=True, dummy_na=False)
# align columns with the Cox model
for c in ontx_x.columns:
    if c not in ret_x.columns: ret_x[c] = 0
ret_x = ret_x[ontx_x.columns]   # same column order
ret_x = ret_x.dropna()
ret_observed_t = ret_x['t'].values; ret_observed_ev = ret_x['ev'].values
print(f'Returners with complete X: {len(ret_x):,}')

# survival functions per returner; evaluate at 24 mo (2 years from m=3 landmark)
sf = cph.predict_survival_function(ret_x.drop(columns=['t','ev']))
T_EVAL = MAX_FU_YR
idx = sf.index.get_indexer([T_EVAL], method='ffill')[0]
s_pred = sf.iloc[idx].values     # one survival prob per returner
counterfactual_mort = 1 - s_pred.mean()
print(f'\nCounterfactual mortality at 24 mo (returners had stayed on tx): {counterfactual_mort*100:.2f}%')

# bootstrap returner X (and the on-tx Cox fit) for 95% CI on the counterfactual mean
B = 500
boot_cf = np.zeros(B)
rng = np.random.default_rng(20260528)
for b in range(B):
    idx_b = rng.integers(0, len(ret_x), len(ret_x))
    s_b = sf.iloc[idx, idx_b].values    # subsample the per-subject survival at T_EVAL
    boot_cf[b] = 1 - s_b.mean()
ci_lo, ci_hi = np.quantile(boot_cf, [0.025, 0.975])

# ---------- observed mortality among returners ----------
kmf = KaplanMeierFitter().fit(ret_x['t'], ret_x['ev'])
s_obs = kmf.survival_function_
ci_obs = kmf.confidence_interval_
idx_o = s_obs.index.get_indexer([T_EVAL], method='ffill')[0]
obs_mort = 1 - float(s_obs.iloc[idx_o,0])
obs_lo   = 1 - float(ci_obs.iloc[idx_o,1])
obs_hi   = 1 - float(ci_obs.iloc[idx_o,0])

# ---------- contrast ----------
risk_diff = obs_mort - counterfactual_mort
risk_ratio = obs_mort / counterfactual_mort if counterfactual_mort>0 else np.nan

print('\n=========================================================')
print(' Returners-at-m=3 — observed vs counterfactual at 24 months ')
print('=========================================================')
print(f'  Observed mortality (had abandoned, then returned):    {100*obs_mort:.2f}%  ({100*obs_lo:.2f}–{100*obs_hi:.2f})')
print(f'  Counterfactual (if had stayed on treatment):          {100*counterfactual_mort:.2f}%  ({100*ci_lo:.2f}–{100*ci_hi:.2f})')
print(f'  Risk difference (excess mortality from abandoning):   {100*risk_diff:.2f} percentage points')
print(f'  Risk ratio (observed / counterfactual):               {risk_ratio:.2f}')

rows = [
    {'metric':'Observed 24-mo mortality among returners-at-m=3',
     'estimate_pct': round(100*obs_mort,2),
     'ci_lo_pct': round(100*obs_lo,2),
     'ci_hi_pct': round(100*obs_hi,2)},
    {'metric':'Counterfactual 24-mo mortality (if had stayed on tx)',
     'estimate_pct': round(100*counterfactual_mort,2),
     'ci_lo_pct': round(100*ci_lo,2),
     'ci_hi_pct': round(100*ci_hi,2)},
    {'metric':'Risk difference (observed - counterfactual), pp',
     'estimate_pct': round(100*risk_diff,2),
     'ci_lo_pct': None, 'ci_hi_pct': None},
    {'metric':'Risk ratio (observed / counterfactual)',
     'estimate_pct': round(risk_ratio,2),
     'ci_lo_pct': None, 'ci_hi_pct': None},
]
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')
