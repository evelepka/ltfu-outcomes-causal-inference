"""
59. On-treatment mortality at the m=3 landmark, stratified by baseline severity
================================================================================
Tests the depletion-of-comparator mechanism: if early-window aHRs below 1.0
reflect acute TB deaths concentrated in severe-baseline patients on treatment
(who deplete the on-treatment comparator before their LTFU counterparts can
disengage), then within the on-treatment arm at the m=3 landmark we should see:
  (a) cumulative mortality that is concentrated early in the follow-up window;
  (b) much steeper early mortality among patients with severe baseline
      characteristics (hospitalised at diagnosis, HIV+, smear+, drug-resistant).

Outputs:
  ITT_Analysis/results/ontx_severity_stratified.csv
  ITT_Analysis/results/ontx_severity_stratified_km.png/.pdf
  ITT_Analysis/results/ontx_severity_hazard.png/.pdf
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'ontx_severity_stratified.csv'
OUT_KM  = BASE/'ITT_Analysis'/'results'/'ontx_severity_stratified_km.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'ontx_severity_stratified_km.pdf'

M_LM_YR = 3.0/12.0
MAX_FU_YR = 2.0

# --------------- cohort ---------------
d = pd.read_csv(COHORT, low_memory=False)
d['best_start']=pd.to_datetime(d['best_start'],errors='coerce')
d['end_date']  =pd.to_datetime(d['end_date'],  errors='coerce')
d['tx_yrs']    =(d['end_date']-d['best_start']).dt.days/365.25

# On-tx at m=3 (matches script 57b/58 definition)
ontx = d[(d['itt_group']!='Loss to follow-up') & (d['tx_yrs']>3.0/12.0)].copy()
ontx['t']    = np.minimum(ontx['time_d_tx'] - M_LM_YR, MAX_FU_YR)
ontx['ev']   = ((ontx['event_d']==1) & ((ontx['time_d_tx']-M_LM_YR)<=MAX_FU_YR)).astype(int)
ontx = ontx[ontx['t']>0].copy()
print(f'On-tx at m=3: N={len(ontx):,}; deaths={ontx.ev.sum():,}')

# --------------- strata definitions ---------------
STRATA = {
    'Hospitalised at diagnosis': ('hosp_admission', 'Yes', 'No'),
    'HIV-positive':              ('hiv_aids',      'Positive', 'Negative'),
    'Smear-positive':            ('bac1_clean',    'Positive', 'Negative'),
    'Drug-resistant TB':         ('resistance_clean', None, None),  # custom below
}

def label_resistance(s):
    s_str = str(s).strip()
    if s_str in ('Multidrug-Resistant TB','MDR-TB','RR-TB','XDR-TB',
                 'Rifampin-Resistant Tuberculosis','Multidrug Resistant',
                 'Polidrogo-resistente','Rifampin-Resistant','RR/MDR',
                 'Multidroga Resistente','Mono Rifampicin'):
        return 'Drug-resistant'
    if 'Resist' in s_str and 'Sensit' not in s_str: return 'Drug-resistant'
    if s_str == 'Sensitive': return 'Sensitive'
    return None  # unknown / not evaluated

# print resistance values to see what's actually in the column
print('\nresistance_clean values:')
print(ontx['resistance_clean'].value_counts(dropna=False).head(10))

# --------------- KM per stratum, both panels (early + full) ---------------
records=[]
km_results={}

# Show: for each stratum (severe vs non-severe), the KM curve over 0–24 months from m=3.
strata_plot = [
    ('Hospitalised at diagnosis', 'hosp_admission', 'Yes', 'No'),
    ('HIV-positive',              'hiv_aids',      'Positive', 'Negative'),
    ('Smear-positive',            'bac1_clean',    'Positive', 'Negative'),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), sharey=True)
colours = {'severe':'#d62728','non-severe':'#2ca02c'}

for ax, (lab, col, sev_lvl, ref_lvl) in zip(axes, strata_plot):
    sev = ontx[ontx[col]==sev_lvl]
    ref = ontx[ontx[col]==ref_lvl]
    for sub, color, name in [(sev, colours['severe'], f'{lab}: yes (n={len(sev):,})'),
                              (ref, colours['non-severe'], f'{lab}: no (n={len(ref):,})')]:
        if len(sub)<10: continue
        kmf = KaplanMeierFitter().fit(sub['t'], sub['ev'], label=name)
        s = kmf.survival_function_; ci = kmf.confidence_interval_
        x = s.index.values*12; y = (1-s.iloc[:,0].values)*100
        ax.step(x, y, where='post', color=color, lw=2.0, label=name)
        try:
            ax.fill_between(x,(1-ci.iloc[:,1])*100,(1-ci.iloc[:,0])*100,
                            color=color, alpha=0.15, step='post')
        except Exception: pass
        # quantitative summary at 6 and 24 mo
        for t_mo in (6, 24):
            idx = s.index.get_indexer([t_mo/12.0], method='ffill')[0]
            idx = max(idx, 0)
            cum_pct = float(1-s.iloc[idx,0])*100
            records.append({'stratum_var':lab,'subgroup':sev_lvl if sub is sev else ref_lvl,
                            'n':len(sub),'time_mo':t_mo,'cum_mortality_pct':round(cum_pct,2)})
    # mark early / late window boundary
    ax.axvline(6, color='grey', ls=':', lw=0.8)
    ax.text(6.2, 0.3, '6-mo boundary', fontsize=8, color='dimgrey')
    ax.set_xlim(0, MAX_FU_YR*12)
    ax.set_xlabel('Months from m=3 landmark')
    ax.set_title(lab, loc='left', fontsize=11, fontweight='bold')
    ax.legend(loc='upper left', frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    for s_ in ('top','right'): ax.spines[s_].set_visible(False)

axes[0].set_ylabel('Cumulative on-treatment mortality (%)')
fig.suptitle('On-treatment 24-month mortality at m=3 landmark, by baseline severity',
             fontsize=12.5, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(OUT_KM, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'\nwrote {OUT_KM}')

# --------------- save table ---------------
tbl = pd.DataFrame(records)
# Pivot to make readable comparison: cum_mortality at 6 mo vs 24 mo by stratum
piv = tbl.pivot_table(index=['stratum_var','subgroup'], columns='time_mo',
                      values='cum_mortality_pct', aggfunc='first').reset_index()
piv.columns.name=None
piv.rename(columns={6:'cum_mort_6mo_pct',24:'cum_mort_24mo_pct'}, inplace=True)
# share of 24-mo mortality occurring in the early window
piv['early_share_pct'] = (piv['cum_mort_6mo_pct']/piv['cum_mort_24mo_pct']*100).round(1)
piv.to_csv(OUT_CSV, index=False)
print(f'wrote {OUT_CSV}\n')
print(piv.to_string(index=False))
