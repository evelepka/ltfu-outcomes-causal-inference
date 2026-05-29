"""
58b. Sensitivity to Analysis 7 — counterfactual using return-state covariates
==============================================================================
Analysis 7 estimated E[Y(stay-on-tx) | returner=1] using each returner's
BASELINE covariates. This implicitly assumed that returners and on-tx
patients with the same baseline X have the same potential outcome under
"continue treatment". Because returners were selected on unmeasured factors
(active TB recrudescence at the time of return) that we do not observe at
baseline, this assumption is plausibly violated in a direction that pushes
the counterfactual estimate downward.

Here we use the clinical characteristics recorded on the RETURN notification
itself as a partial proxy for those unmeasured factors. Demographic and
social variables (age, sex, race, education, homelessness, alcohol, drug
use, dot_status) are kept at their baseline values; clinical/severity
variables (hosp_admission, clinical form, bac1, sputum culture, lab
confirmation, drug resistance, HIV) are replaced with their values at the
return notification.

Two counterfactual estimates are reported:
  (A) Baseline-only:   uses first-notification covariates for returners — the
                       Analysis-7 estimand (= 3.63%).
  (B) Return-state:    uses return-notification clinical covariates +
                       baseline social/demographic covariates for returners.
The pair brackets the truth under standard interpretations: (A) is appropriate
under the strong exchangeability assumption; (B) controls for the part of
unmeasured baseline severity captured by the return notification.

Output:
  ITT_Analysis/results/counterfactual_returners_returnstate.csv
  ITT_Analysis/results/return_state_vs_baseline_table.csv
"""
from pathlib import Path
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
RAW    = BASE/'Data'/'Final_table_cleaned.csv'
OUT_CSV   = BASE/'ITT_Analysis'/'results'/'counterfactual_returners_returnstate.csv'
OUT_TABLE = BASE/'ITT_Analysis'/'results'/'return_state_vs_baseline_table.csv'

M_LM_YR = 3.0/12.0
MAX_FU_YR = 2.0

# --- cohort ---
print('Loading cohort...')
d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'],errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],  errors='coerce')
d['tx_yrs']    =(d['end_date']-d['best_start']).dt.days/365.25

# LTFU at m=3 AND returned by 6mo of LTFU → returners cohort (matches script 58)
ltfu_m3 = d[(d['itt_group']=='Loss to follow-up')
            & (d['tx_yrs']>2.0/12.0) & (d['tx_yrs']<=3.0/12.0)].copy()
returners = ltfu_m3[(ltfu_m3['event_rn']==1) & (ltfu_m3['time_rn']<=0.5)].copy()
returners['t']  = np.minimum(returners['time_d'], MAX_FU_YR)
returners['ev'] = ((returners['event_d']==1) & (returners['time_d']<=MAX_FU_YR)).astype(int)
print(f'Returners at m=3: N={len(returners):,}; late-window deaths={returners.ev.sum()}')

# --- raw data: pull each returner's return-notification clinical covariates ---
print('Loading raw clinical covariates...')
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=['sinan_clean','case_type','notification_date',
                           'hosp_admission','clinical_classif','lab_confirmed',
                           'sputum_culture','bac1','resistance','hiv'])
raw['notification_date']=pd.to_datetime(raw['notification_date'],errors='coerce')

# For each returner, find their RETURN notification = the first notification with
# notification_date > LTFU date (= end_date of first episode).
ret_ids = returners[['sinan_clean','end_date']].rename(columns={'end_date':'ltfu_date'})
m = raw.merge(ret_ids, on='sinan_clean', how='inner')
m = m[m['notification_date'] > m['ltfu_date']].copy()
m = m.sort_values(['sinan_clean','notification_date'])
return_notif = m.drop_duplicates('sinan_clean', keep='first').copy()
print(f'Return notifications found for {len(return_notif):,} of {len(returners):,} returners')

# --- map raw codes to cohort-clean labels (mirroring script 01 conventions) ---
def map_hosp(s): return 'Yes' if s=='S' else ('No' if s=='N' else np.nan)
def map_clinical(s):
    return {'Pul':'Pulmonary','Ext':'Extrapulmonary',
            'P+E':'Pulmonary and Extrapulmonary or disseminated',
            'Dissem':'Pulmonary and Extrapulmonary or disseminated'}.get(s, np.nan)
def map_lab(v):
    try: v=int(v); return 'Yes' if v==1 else 'No'
    except Exception: return np.nan
def map_smear(s):
    if s=='Pos': return 'Positive'
    if s=='Neg': return 'Negative'
    return 'Not Evaluated'
def map_sputum(s):
    if s=='Pos': return 'Positive'
    if s=='Neg': return 'Negative'
    return 'Not Evaluated'
def map_resist(s):
    if s=='SENS': return 'Sensitive'
    if isinstance(s,str) and s.startswith('TB'): return 'Resistant (Any)'
    return 'Not Evaluated'
def map_hiv(s):
    if s=='Pos': return 'Positive'
    if s=='Neg': return 'Negative'
    return 'Unknown'

return_notif['hosp_admission_ret']      = return_notif['hosp_admission'].map(map_hosp)
return_notif['clinical_clean_ret']      = return_notif['clinical_classif'].map(map_clinical)
return_notif['lab_confirmed_stat_ret']  = return_notif['lab_confirmed'].map(map_lab)
return_notif['bac1_clean_ret']          = return_notif['bac1'].map(map_smear)
return_notif['sputum_culture_clean_ret']= return_notif['sputum_culture'].map(map_sputum)
return_notif['resistance_clean_ret']    = return_notif['resistance'].map(map_resist)
return_notif['hiv_aids_ret']            = return_notif['hiv'].map(map_hiv)

ret_cols = ['hosp_admission_ret','clinical_clean_ret','lab_confirmed_stat_ret',
            'bac1_clean_ret','sputum_culture_clean_ret','resistance_clean_ret','hiv_aids_ret']

# Merge return-state covariates onto returner cohort
returners = returners.merge(return_notif[['sinan_clean']+ret_cols], on='sinan_clean', how='left')

# --- comparison: baseline (first notification) vs return-state distributions ---
def pct(s, lvl): return 100*(s==lvl).mean()
cmp_rows=[]
PAIRS = [
    ('Hospitalised at diagnosis (Yes)',      'hosp_admission',          'hosp_admission_ret',          'Yes'),
    ('Pulmonary TB',                          'clinical_clean',           'clinical_clean_ret',           'Pulmonary'),
    ('Extrapulmonary TB',                     'clinical_clean',           'clinical_clean_ret',           'Extrapulmonary'),
    ('Pulmonary + extrapulmonary',            'clinical_clean',           'clinical_clean_ret',           'Pulmonary and Extrapulmonary or disseminated'),
    ('Lab-confirmed (Yes)',                   'lab_confirmed_stat',       'lab_confirmed_stat_ret',       'Yes'),
    ('Smear-positive',                        'bac1_clean',               'bac1_clean_ret',               'Positive'),
    ('Sputum culture positive',               'sputum_culture_clean',     'sputum_culture_clean_ret',     'Positive'),
    ('Drug-resistant TB',                     'resistance_clean',         'resistance_clean_ret',         'Resistant (Any)'),
    ('HIV-positive',                          'hiv_aids',                 'hiv_aids_ret',                 'Positive'),
]
for label, base_col, ret_col, lvl in PAIRS:
    p_b = pct(returners[base_col], lvl)
    p_r = pct(returners[ret_col],  lvl)
    cmp_rows.append({'characteristic':label, 'baseline_pct':round(p_b,1),
                     'return_pct':round(p_r,1), 'shift_pp':round(p_r-p_b,1)})
cmp = pd.DataFrame(cmp_rows)
OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
cmp.to_csv(OUT_TABLE, index=False)
print(f'\nwrote {OUT_TABLE}')
print(cmp.to_string(index=False))

# --- fit on-tx Cox (same as script 58) ---
ontx = d[(d['itt_group']!='Loss to follow-up') & (d['tx_yrs']>3.0/12.0)].copy()
ontx['t_from_lm'] = ontx['time_d_tx'] - M_LM_YR
ontx = ontx[ontx['t_from_lm']>0]
ontx['t']  = np.minimum(ontx['t_from_lm'], MAX_FU_YR)
ontx['ev'] = ((ontx['event_d']==1) & (ontx['t_from_lm']<=MAX_FU_YR)).astype(int)

COVARS = ['age_tb','sex','hiv_aids','hosp_admission','clinical_clean',
          'homelessness','alcohol','drug_use','dot_status']
ontx_x = ontx[['t','ev']+COVARS].copy()
ontx_x = pd.get_dummies(ontx_x, columns=[c for c in COVARS if c!='age_tb'],
                        drop_first=True, dummy_na=False)
ontx_x = ontx_x.dropna().astype({c: float for c in ontx_x.columns if ontx_x[c].dtype==bool})
cph = CoxPHFitter(penalizer=0.001).fit(ontx_x, duration_col='t', event_col='ev')
print(f'\nOn-tx Cox concordance: {cph.concordance_index_:.3f}')

# --- prepare two prediction frames for returners ---
# (A) baseline X for clinical (HIV, hosp, clinical_clean) + baseline X for demographic/social
ret_A = returners[['t','ev']+COVARS].copy()

# (B) return-state X for clinical (HIV, hosp, clinical_clean) + baseline X for demographic/social
ret_B = returners[['t','ev']+COVARS].copy()
overrides = {
    'hiv_aids':          'hiv_aids_ret',
    'hosp_admission':    'hosp_admission_ret',
    'clinical_clean':    'clinical_clean_ret',
}
for tgt, src in overrides.items():
    ret_B[tgt] = returners[src].fillna(returners[tgt])

def predict_cf(df_):
    df_ = pd.get_dummies(df_, columns=[c for c in COVARS if c!='age_tb'],
                          drop_first=True, dummy_na=False)
    for c in ontx_x.columns:
        if c not in df_.columns: df_[c] = 0
    df_ = df_[ontx_x.columns]
    df_ = df_.dropna()
    sf = cph.predict_survival_function(df_.drop(columns=['t','ev']))
    idx = sf.index.get_indexer([MAX_FU_YR], method='ffill')[0]; idx=max(idx,0)
    return 1 - sf.iloc[idx].values, df_

cf_A_arr, df_A = predict_cf(ret_A)
cf_B_arr, df_B = predict_cf(ret_B)

# bootstrap CIs for the mean counterfactuals (subjects only; fixed Cox)
def boot_ci(arr, B=500):
    rng = np.random.default_rng(20260528)
    boot = np.array([arr[rng.integers(0,len(arr),len(arr))].mean() for _ in range(B)])
    return float(np.mean(arr)), tuple(np.quantile(boot,[0.025,0.975]))

mean_A, ci_A = boot_ci(cf_A_arr)
mean_B, ci_B = boot_ci(cf_B_arr)

# observed mortality KM
kmf = KaplanMeierFitter().fit(returners['t'], returners['ev'])
idx_o = kmf.survival_function_.index.get_indexer([MAX_FU_YR], method='ffill')[0]
obs = 1 - float(kmf.survival_function_.iloc[idx_o,0])
ci_obs = kmf.confidence_interval_
obs_lo = 1-float(ci_obs.iloc[idx_o,1]); obs_hi = 1-float(ci_obs.iloc[idx_o,0])

print('\n=================================================================')
print(' 24-mo mortality among returners — observed vs two counterfactuals ')
print('=================================================================')
print(f'  Observed:                                            {100*obs:.2f}%  ({100*obs_lo:.2f}–{100*obs_hi:.2f})')
print(f'  Counterfactual A  (baseline X for returners):        {100*mean_A:.2f}%  ({100*ci_A[0]:.2f}–{100*ci_A[1]:.2f})')
print(f'  Counterfactual B  (return-state clinical X):         {100*mean_B:.2f}%  ({100*ci_B[0]:.2f}–{100*ci_B[1]:.2f})')
print()
print(f'  Risk diff vs A:                                      {100*(obs-mean_A):.2f} pp     Risk ratio: {obs/mean_A:.2f}')
print(f'  Risk diff vs B:                                      {100*(obs-mean_B):.2f} pp     Risk ratio: {obs/mean_B:.2f}')

rows=[
    {'metric':'Observed 24-mo mortality (returners-at-m=3)',
     'estimate_pct': round(100*obs,2),
     'ci_lo': round(100*obs_lo,2), 'ci_hi': round(100*obs_hi,2)},
    {'metric':'Counterfactual A — baseline X for returners',
     'estimate_pct': round(100*mean_A,2),
     'ci_lo': round(100*ci_A[0],2), 'ci_hi': round(100*ci_A[1],2)},
    {'metric':'Counterfactual B — return-state clinical X for returners',
     'estimate_pct': round(100*mean_B,2),
     'ci_lo': round(100*ci_B[0],2), 'ci_hi': round(100*ci_B[1],2)},
    {'metric':'Risk difference vs A (pp)',
     'estimate_pct': round(100*(obs-mean_A),2), 'ci_lo': None, 'ci_hi': None},
    {'metric':'Risk difference vs B (pp)',
     'estimate_pct': round(100*(obs-mean_B),2), 'ci_lo': None, 'ci_hi': None},
    {'metric':'Risk ratio vs A',
     'estimate_pct': round(obs/mean_A,2), 'ci_lo': None, 'ci_hi': None},
    {'metric':'Risk ratio vs B',
     'estimate_pct': round(obs/mean_B,2), 'ci_lo': None, 'ci_hi': None},
]
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')
