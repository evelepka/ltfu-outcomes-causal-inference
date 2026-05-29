"""
55. Landmark SMD analysis
=========================
For landmark months m ∈ {1, 3, 6}, compare baseline characteristics between
the LTFU arm (disengaged during month m) and the on-treatment comparator
(still on treatment past month m). Standardised mean differences (SMDs) are
reported by panel:
  • TB-severity / early-mortality predictors
  • Social vulnerability

This is the apples-to-apples comparison the target trial emulates at each
landmark; it is NOT Table 1 (which compares overall LTFU vs non-LTFU at
treatment initiation).

Output: ITT_Analysis/results/landmark_baseline_smd.csv
"""
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT    = BASE/'ITT_Analysis'/'results'/'landmark_baseline_smd.csv'

SEVERITY = ['age_tb','hiv_aids','hosp_admission','bac1_clean','sputum_culture_clean',
            'lab_confirmed_stat','resistance_clean','clinical_clean','diabetes',
            'other_immuno_condition','diagnosis_setting']
SOCIAL   = ['sex','race_clean','edu_clean','homelessness','alcohol','drug_use',
            'mental_health','tobacco_use','incarcerated','dot_status']

# ----- load -----
d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'], errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],   errors='coerce')
d['tx_months'] = (d['end_date']-d['best_start']).dt.days / 30.4375
d = d.dropna(subset=['tx_months','itt_group'])
print(f'cohort N={len(d):,}; LTFU N={(d.itt_group=="Loss to follow-up").sum():,}')

# ----- SMD helpers -----
def smd_continuous(a, b):
    a,b = a.dropna(), b.dropna()
    if len(a)<2 or len(b)<2: return np.nan
    sa, sb = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt((sa+sb)/2)
    return (a.mean()-b.mean())/pooled if pooled>0 else np.nan

def smd_binary(p1, p0):
    pbar = (p1+p0)/2
    denom = np.sqrt(pbar*(1-pbar))
    return (p1-p0)/denom if denom>0 else np.nan

def summarise(arm_ltfu, arm_otx, var, kind):
    """Return list of rows: variable / level / pct_or_mean per arm / SMD."""
    rows=[]
    if kind=='continuous':
        m1,m0 = arm_ltfu[var].mean(), arm_otx[var].mean()
        s1,s0 = arm_ltfu[var].std(),  arm_otx[var].std()
        rows.append({'variable':var, 'level':'mean (SD)',
                     'ltfu': f'{m1:.1f} ({s1:.1f})',
                     'otx':  f'{m0:.1f} ({s0:.1f})',
                     'smd':  smd_continuous(arm_ltfu[var], arm_otx[var])})
    else:
        levels = sorted(pd.concat([arm_ltfu[var],arm_otx[var]]).dropna().unique().tolist())
        n1, n0 = len(arm_ltfu), len(arm_otx)
        for L in levels:
            p1 = (arm_ltfu[var]==L).mean()
            p0 = (arm_otx[var]==L).mean()
            rows.append({'variable':var, 'level':str(L),
                         'ltfu': f'{p1*100:.1f}%  ({(arm_ltfu[var]==L).sum()}/{n1})',
                         'otx':  f'{p0*100:.1f}%  ({(arm_otx[var]==L).sum()}/{n0})',
                         'smd':  smd_binary(p1,p0)})
    return rows

# ----- landmark arm definitions -----
def arms_at(m):
    """LTFU arm = disengaged in month m (m-1 < tx_months <= m);
       On-treatment arm = still on treatment past month m (tx_months > m)."""
    eligible = d[d['tx_months']>m-1]   # alive on treatment at least through month m-1
    ltfu = eligible[(eligible['itt_group']=='Loss to follow-up')
                    & (eligible['tx_months']> m-1)
                    & (eligible['tx_months']<= m)]
    otx  = eligible[eligible['tx_months']> m]
    return ltfu, otx

# ----- run -----
records=[]
for m in (1,3,6):
    ltfu, otx = arms_at(m)
    print(f'\nLandmark month {m}: LTFU N={len(ltfu):,}; on-tx N={len(otx):,}')
    for panel, vars_ in (('severity',SEVERITY),('social',SOCIAL)):
        for v in vars_:
            kind = 'continuous' if v=='age_tb' else 'categorical'
            for row in summarise(ltfu, otx, v, kind):
                row.update({'landmark_m':m, 'panel':panel,
                            'n_ltfu':len(ltfu), 'n_otx':len(otx)})
                records.append(row)

out = pd.DataFrame(records, columns=['landmark_m','panel','variable','level',
                                      'n_ltfu','ltfu','n_otx','otx','smd'])
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)
print(f'\nwrote {OUT}')

# ----- quick summary: |SMD| > 0.1 imbalances in severity panel, by landmark -----
print('\n\n==== |SMD| > 0.10 in SEVERITY panel ====')
sev = out[(out['panel']=='severity') & out['smd'].abs().gt(0.10)]
for m in (1,3,6):
    print(f'\n-- landmark m={m} --')
    sub = sev[sev['landmark_m']==m].copy()
    sub['dir'] = sub['smd'].apply(lambda x: '↑LTFU' if x>0 else '↓LTFU')
    sub = sub.sort_values('smd', key=lambda s: -s.abs())
    for _,r in sub.iterrows():
        print(f'  {r.dir}  SMD={r.smd:+.2f}  {r.variable}={r.level}   '
              f'LTFU {r.ltfu.split()[0]:>6} vs OnTx {r.otx.split()[0]:>6}')

print('\n==== |SMD| > 0.10 in SOCIAL panel ====')
soc = out[(out['panel']=='social') & out['smd'].abs().gt(0.10)]
for m in (1,3,6):
    print(f'\n-- landmark m={m} --')
    sub = soc[soc['landmark_m']==m].copy()
    sub['dir'] = sub['smd'].apply(lambda x: '↑LTFU' if x>0 else '↓LTFU')
    sub = sub.sort_values('smd', key=lambda s: -s.abs())
    for _,r in sub.iterrows():
        print(f'  {r.dir}  SMD={r.smd:+.2f}  {r.variable}={r.level}   '
              f'LTFU {r.ltfu.split()[0]:>6} vs OnTx {r.otx.split()[0]:>6}')
