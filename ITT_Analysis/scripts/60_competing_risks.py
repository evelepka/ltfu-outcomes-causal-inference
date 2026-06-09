"""
60. Competing risks framing — making the survival precondition visible
========================================================================
Direct visualisation of the structural argument underlying mechanism #1:
"death on treatment" and "LTFU" are mutually-exclusive competing events
from the moment treatment starts. Among patients with severe baseline TB,
the early cumulative incidence of death-on-treatment is substantial; some
fraction of those patients die on treatment before they ever have the
opportunity to disengage. The LTFU classification is therefore structurally
unavailable to them.

For each baseline-severity stratum, compute Aalen–Johansen cumulative
incidence functions for two competing events from best_start:
  • Death on treatment (case_outcome ∈ {Obito TB, Obito NTB})
  • LTFU              (case_outcome ∈ {Abandono, Abandono Primario, Faltoso})

Outputs:
  ITT_Analysis/results/competing_risks_by_severity.csv
  ITT_Analysis/results/competing_risks_by_severity.png/.pdf
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from lifelines import AalenJohansenFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
RAW    = BASE/'Data'/'Final_table_cleaned.csv'
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'competing_risks_by_severity.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'competing_risks_by_severity.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'competing_risks_by_severity.pdf'

MAX_MO = 24

# --------------- merge first-Novo case_outcome onto cohort ---------------
print('Loading raw...')
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=['sinan_clean','case_type','case_outcome','end_date'])
raw['end_date']=pd.to_datetime(raw['end_date'], format='%B %d, %Y', errors='coerce')
TRANSFER = {'Transf Outro Municipio','Transf Outro Estado/Pais'}
novo = raw[raw['case_type'].str.strip().str.lower().eq('novo')
           & raw['case_outcome'].notna()
           & raw['case_outcome'].str.strip().ne('')
           & raw['case_outcome'].ne('Mud Diag')
           & ~raw['case_outcome'].isin(TRANSFER)].copy()
novo = novo.sort_values('end_date')
first = novo.drop_duplicates('sinan_clean', keep='first')[['sinan_clean','case_outcome']]

d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'],errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],  errors='coerce')
d = d.merge(first, on='sinan_clean', how='left')
d['tx_mo'] = (d['end_date']-d['best_start']).dt.days/30.4375
d = d[(d['tx_mo']>=0) & d['case_outcome'].notna()].copy()
print(f'cohort with case_outcome: {len(d):,}')
print('first-Novo case_outcome distribution:')
print(d['case_outcome'].value_counts().head(10))

# --------------- define competing events ---------------
DEATH_OUTCOMES = {'Obito TB','Obito NTB'}
LTFU_OUTCOMES  = {'Abandono','Abandono Primario','Faltoso'}
# Map to event integer: 0 = censored (still on tx / cure / other), 1 = death on tx, 2 = LTFU
def event(o):
    if o in DEATH_OUTCOMES: return 1
    if o in LTFU_OUTCOMES:  return 2
    return 0
d['ev'] = d['case_outcome'].map(event)
d['t']  = np.minimum(d['tx_mo'], MAX_MO)
d.loc[d['tx_mo']>MAX_MO, 'ev'] = 0   # admin censor at MAX_MO

print('\nEvent counts (within 24 mo):')
print(pd.Series(d['ev']).value_counts())

# --------------- stratify by baseline severity ---------------
STRATA = [
    ('Hospitalised at diagnosis', 'hosp_admission', 'Yes', 'No'),
    ('HIV-positive',              'hiv_aids',      'Positive', 'Negative'),
]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)

records=[]
for ax, (label, col, sev_lvl, ref_lvl) in zip(axes, STRATA):
    sev = d[d[col]==sev_lvl]
    ref = d[d[col]==ref_lvl]

    for grp_label, grp, ls in [(f'{label}: yes', sev, '-'), (f'{label}: no', ref, '--')]:
        aj_death = AalenJohansenFitter(calculate_variance=False)
        aj_ltfu  = AalenJohansenFitter(calculate_variance=False)
        aj_death.fit(grp['t'], grp['ev'], event_of_interest=1)
        aj_ltfu .fit(grp['t'], grp['ev'], event_of_interest=2)
        # plot
        cif_d = aj_death.cumulative_density_
        cif_l = aj_ltfu .cumulative_density_
        ax.step(cif_d.index, cif_d.iloc[:,0]*100, where='post',
                color='#d62728', ls=ls, lw=2.0,
                label=f'Death on treatment — {grp_label} (n={len(grp):,})')
        ax.step(cif_l.index, cif_l.iloc[:,0]*100, where='post',
                color='#1f77b4', ls=ls, lw=2.0,
                label=f'LTFU — {grp_label} (n={len(grp):,})')
        # record CIF at 3, 6, 12, 24 months
        for t_mo in (3, 6, 12, 24):
            for kind, aj in [('death_on_tx',aj_death),('ltfu',aj_ltfu)]:
                cif = aj.cumulative_density_
                idx = cif.index.get_indexer([t_mo], method='ffill')[0]
                idx = max(idx, 0)
                records.append({'stratum_var':label,'subgroup':sev_lvl if grp is sev else ref_lvl,
                                'event':kind,'time_mo':t_mo,
                                'cif_pct':round(float(cif.iloc[idx,0])*100,2),
                                'n':len(grp)})
    ax.set_xlim(0, MAX_MO)
    ax.set_xlabel('Months since treatment start')
    ax.set_title(label, loc='left', fontsize=11, fontweight='bold')
    ax.grid(alpha=0.25)
    ax.legend(loc='upper left', frameon=False, fontsize=8)
    for s_ in ('top','right'): ax.spines[s_].set_visible(False)
axes[0].set_ylabel('Cumulative incidence (%)')

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'\nwrote {OUT_PNG}')

tbl = pd.DataFrame(records)
tbl.to_csv(OUT_CSV, index=False)
print(f'wrote {OUT_CSV}\n')
# Quick summary print
print('Cumulative incidence at 3 and 6 months (key for the structural argument):')
for stratum_var in tbl['stratum_var'].unique():
    sub = tbl[(tbl['stratum_var']==stratum_var) & (tbl['time_mo'].isin([3,6]))].copy()
    print(f'\n  {stratum_var}:')
    piv = sub.pivot_table(index=['subgroup','time_mo'], columns='event',
                          values='cif_pct', aggfunc='first')
    piv['death_vs_ltfu_ratio'] = (piv['death_on_tx']/piv['ltfu']).round(2)
    print(piv)
