"""
56b. Cause-specific late-window mortality among LTFU patients, by 6-month return status
=======================================================================================
Extends script 56 (all-cause) to TB-attributable vs non-TB death. Uses the
"hybrid" cause attribution matching the manuscript (script 30e):
  TB death     = SIM ICD-10 ^A1[5-9] | ^B90 | ^B200,  OR  TBweb case_outcome == 'Obito TB'
  Non-TB death = SIM ICD-10 known and not TB/respiratory/HIV-other,
                 OR TBweb case_outcome == 'Obito NTB'
  Unknown      = death without SIM code and without TBweb Obito attribution
(Respiratory and HIV-other are pulled out of the SIM bucket and treated as
unknown for the hybrid attribution, matching the primary specification.)

For cause-specific KM/Cox, deaths from the *other* cause are censored at the
death time, and unknown-cause deaths are censored at the death time.

Outputs:
  • ITT_Analysis/results/return_stratified_cause_specific.csv
  • ITT_Analysis/results/return_stratified_km_cause_specific.png/pdf
"""
from pathlib import Path
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter

# Project root: env override, else resolve relative to this file. The previous
# hardcoded absolute path pointed at 'Abandonment Paper', but the shared folder
# is now 'LTFU Paper', so the script could not open its own inputs. 22 of the
# 24 Python scripts in this directory still carry the stale path -- see
# docs/dead-ends.md.
import os
BASE = Path(os.environ.get("TB_ABANDONMENT_ROOT", "")) if os.environ.get(
    "TB_ABANDONMENT_ROOT") else Path(__file__).resolve().parents[2]
RAW    = BASE/'Data'/'Final_table_cleaned.csv'
COHORT = BASE/'ITT_Analysis'/'data'/'itt_cohort.csv'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'return_stratified_cause_specific_fixedattr.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'return_stratified_km_cause_specific_fixedattr.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'return_stratified_km_cause_specific_fixedattr.pdf'

LANDMARK_YR = 0.5
MAX_FU_YR   = 2.0

# ---------- cause attribution lookup (mirrors 30e classify_cod) ----------
print('Loading raw data for cause attribution...')
raw = pd.read_csv(RAW, low_memory=False,
                  usecols=['sinan_clean','case_type','case_outcome','end_date','dod','cause_of_death_code'])
raw['end_date'] = pd.to_datetime(raw['end_date'], format='%B %d, %Y', errors='coerce')
raw['dod']      = pd.to_datetime(raw['dod'],      format='%B %d, %Y', errors='coerce')

TRANSFER = {'Transf Outro Municipio','Transf Outro Estado/Pais'}
novo = raw[raw['case_type'].str.strip().str.lower().eq('novo')
           & raw['case_outcome'].notna()
           & raw['case_outcome'].str.strip().ne('')
           & raw['case_outcome'].ne('Mud Diag')
           & ~raw['case_outcome'].isin(TRANSFER)].copy()
novo = novo.sort_values('end_date')
first_outcome = novo.drop_duplicates('sinan_clean', keep='first')[['sinan_clean','case_outcome']]

# --- ADR-0003 FIX: take the Obito outcome from ANY episode ----------------
# Index-only lookup cannot see an LTFU patient's death, because their index
# episode closes as `Abandono`. 1,058 of 1,668 LTFU deaths (63.4%) are
# recorded on a retreatment episode instead. Verified same-death (median
# 0-day lag to death_date). This matters most HERE: without it, the
# returner stratum is precisely the stratum whose causes get discarded.
_ob = raw[raw['case_outcome'].astype(str).str.strip().isin(['Obito TB','Obito NTB'])]
_ob = _ob.sort_values('end_date').drop_duplicates('sinan_clean', keep='last')[
    ['sinan_clean','case_outcome']].rename(columns={'case_outcome':'_obito'})
first_outcome = first_outcome.merge(_ob, on='sinan_clean', how='outer')
first_outcome['case_outcome'] = first_outcome['_obito'].combine_first(
    first_outcome['case_outcome'])
first_outcome = first_outcome.drop(columns=['_obito'])
print(f"[ADR-0003] Obito outcome recovered from any episode for {len(_ob):,} individuals")

deathrec = raw[raw['dod'].notna() & raw['cause_of_death_code'].notna()
               & raw['cause_of_death_code'].str.strip().ne('')].copy()
deathrec['cause_of_death_code'] = deathrec['cause_of_death_code'].str.strip().str.upper()
deathrec = deathrec.sort_values('dod').drop_duplicates('sinan_clean', keep='last')[
    ['sinan_clean','cause_of_death_code']]

attr = pd.merge(first_outcome, deathrec, on='sinan_clean', how='outer')

def classify(cod, case_out):
    if pd.notna(cod) and cod:
        if re.match(r'^A1[5-9]', cod) or cod.startswith('B90') or cod.startswith('B200'):
            return 'tb_strict'
        if re.match(r'^J[0-9]', cod):
            return 'respiratory'
        if re.match(r'^B2[0-4]', cod) and not cod.startswith('B200'):
            return 'hiv_other'
        return 'non_tb'
    if case_out=='Obito TB':  return 'tb_via_tbweb'
    if case_out=='Obito NTB': return 'ntb_via_tbweb'
    return 'unknown'

attr['cod_class'] = [classify(c,co) for c,co in zip(attr['cause_of_death_code'], attr['case_outcome'])]
# hybrid (primary)
attr['tb_hybrid']    = attr['cod_class'].isin(['tb_strict','tb_via_tbweb'])
attr['nontb_hybrid'] = attr['cod_class'].isin(['non_tb','ntb_via_tbweb'])
print('cod_class distribution among individuals with attribution:')
print(attr['cod_class'].value_counts())

# ---------- LTFU cohort + landmark eligibility ----------
d = pd.read_csv(COHORT, low_memory=False)
ltfu = d[d['itt_group']=='Loss to follow-up'].copy()
ltfu = ltfu.merge(attr[['sinan_clean','tb_hybrid','nontb_hybrid','cod_class']],
                  on='sinan_clean', how='left')
ltfu['tb_hybrid']    = ltfu['tb_hybrid'].fillna(False)
ltfu['nontb_hybrid'] = ltfu['nontb_hybrid'].fillna(False)

ltfu['returner_6mo'] = ((ltfu['event_rn']==1) & (ltfu['time_rn']<=LANDMARK_YR)).astype(int)
ltfu['died_before_landmark'] = ((ltfu['event_d']==1) & (ltfu['time_d']<=LANDMARK_YR)).astype(int)
lm = ltfu[ltfu['died_before_landmark']==0].copy()
print(f'\nLandmark-eligible (alive at 6 mo) N={len(lm):,}')
print('Among deaths in landmark cohort during late window (time_d in (0.5, 2]):')
late_deaths = lm[(lm['event_d']==1) & (lm['time_d']>LANDMARK_YR) & (lm['time_d']<=MAX_FU_YR)]
print(f'  total late-window deaths: {len(late_deaths)}')
print(f'  TB-hybrid:    {late_deaths.tb_hybrid.sum()}')
print(f'  non-TB hybrid:{late_deaths.nontb_hybrid.sum()}')
print(f'  unknown (no attribution):{len(late_deaths)-late_deaths.tb_hybrid.sum()-late_deaths.nontb_hybrid.sum()}')

lm['fu_yr']    = np.minimum(lm['time_d'], MAX_FU_YR) - LANDMARK_YR
lm = lm[lm['fu_yr']>0].copy()

# Cause-specific event indicators: only the focal cause counts as event;
# deaths from any other cause (or unknown) are censored at death time.
lm['ev_overall'] = ((lm['event_d']==1) & (lm['time_d']<=MAX_FU_YR)).astype(int)
lm['ev_tb']      = (lm['ev_overall'].astype(bool) & lm['tb_hybrid']).astype(int)
lm['ev_ntb']     = (lm['ev_overall'].astype(bool) & lm['nontb_hybrid']).astype(int)

# ---------- KM + Cox per cause ----------
def cum_at(kmf, t_yr_from_ltfu):
    t = t_yr_from_ltfu - LANDMARK_YR
    s_func = kmf.survival_function_
    ci_func = kmf.confidence_interval_
    idx = s_func.index.get_indexer([t], method='ffill')[0]
    idx = max(idx, 0)
    s = float(s_func.iloc[idx,0])
    lo = float(ci_func.iloc[idx,0]); hi = float(ci_func.iloc[idx,1])
    return 1-s, (1-hi, 1-lo)

rows=[]; km_curves={}
for cause_label, ev_col, color_R, color_N in [
        ('TB-attributable',    'ev_tb',  '#1f77b4', '#aec7e8'),
        ('Non-TB',             'ev_ntb', '#8c564b', '#c49c94')]:
    R = lm[lm['returner_6mo']==1].copy()
    N = lm[lm['returner_6mo']==0].copy()
    kmf_R = KaplanMeierFitter().fit(R['fu_yr'], R[ev_col], label=f'{cause_label}: returner')
    kmf_N = KaplanMeierFitter().fit(N['fu_yr'], N[ev_col], label=f'{cause_label}: non-returner')
    km_curves[cause_label] = (kmf_R, kmf_N, color_R, color_N, len(R), len(N))
    for t_mo in (12, 18, 24):
        for arm, kmf, n in [('Returner by 6mo', kmf_R, len(R)),
                            ('Non-returner by 6mo', kmf_N, len(N))]:
            m,(lo,hi) = cum_at(kmf, t_mo/12.0)
            rows.append({'cause': cause_label, 'arm': arm, 'n': n,
                         'metric': f'Cumulative {cause_label.lower()} mortality at {t_mo} mo',
                         'estimate_pct': round(100*m,2),
                         'ci_low_pct':   round(100*lo,2),
                         'ci_high_pct':  round(100*hi,2)})
    # Cox: crude and adjusted
    crude = lm[['fu_yr',ev_col,'returner_6mo']].rename(columns={ev_col:'event'})
    cph = CoxPHFitter().fit(crude, duration_col='fu_yr', event_col='event')
    hr = float(cph.hazard_ratios_['returner_6mo'])
    lo,hi = [float(np.exp(x)) for x in cph.confidence_intervals_.loc['returner_6mo']]
    p = float(cph.summary.loc['returner_6mo','p'])
    rows.append({'cause':cause_label,'arm':'returner vs non-returner (crude)','n':len(lm),
                 'metric':'Crude Cox HR (cause-specific)',
                 'estimate_pct':round(hr,3),'ci_low_pct':round(lo,3),'ci_high_pct':round(hi,3)})
    print(f'\n{cause_label} crude HR (returner vs non-returner) = {hr:.2f} ({lo:.2f}–{hi:.2f}), p={p:.2e}')

    adj_covars = ['age_tb','sex','hiv_aids','homelessness','alcohol','drug_use',
                  'hosp_admission','diabetes','clinical_clean','resistance_clean',
                  'lab_confirmed_stat','dot_status','incarcerated']
    adj_df = lm[['fu_yr',ev_col,'returner_6mo']+adj_covars].rename(columns={ev_col:'event'})
    adj_df = pd.get_dummies(adj_df, columns=[c for c in adj_covars if c!='age_tb'],
                             drop_first=True, dummy_na=False)
    adj_df = adj_df.dropna().astype({c:float for c in adj_df.columns if adj_df[c].dtype==bool})
    cph_a = CoxPHFitter(penalizer=0.001).fit(adj_df, duration_col='fu_yr', event_col='event')
    hr_a = float(cph_a.hazard_ratios_['returner_6mo'])
    lo_a,hi_a = [float(np.exp(x)) for x in cph_a.confidence_intervals_.loc['returner_6mo']]
    p_a = float(cph_a.summary.loc['returner_6mo','p'])
    rows.append({'cause':cause_label,'arm':'returner vs non-returner (adjusted)','n':len(adj_df),
                 'metric':'Adjusted Cox HR (cause-specific)',
                 'estimate_pct':round(hr_a,3),'ci_low_pct':round(lo_a,3),'ci_high_pct':round(hi_a,3)})
    print(f'{cause_label} adjusted HR (returner vs non-returner) = {hr_a:.2f} ({lo_a:.2f}–{hi_a:.2f}), p={p_a:.2e}')

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')

# ---------- KM plot ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
def plot_curves(ax, cause):
    kmf_R, kmf_N, cR, cN, nR, nN = km_curves[cause]
    for kmf, color, lab in [(kmf_R, cR, f'Returner by 6 mo (n={nR:,})'),
                            (kmf_N, '#d62728', f'Non-returner by 6 mo (n={nN:,})')]:
        s = kmf.survival_function_
        ci = kmf.confidence_interval_
        x = (s.index.values + LANDMARK_YR)*12
        y = (1-s.iloc[:,0].values)*100
        ax.step(x, y, where='post', color=color, lw=2.0, label=lab)
        ax.fill_between(x, (1-ci.iloc[:,1])*100, (1-ci.iloc[:,0])*100,
                        color=color, alpha=0.15, step='post')
    ax.set_xlim(LANDMARK_YR*12, MAX_FU_YR*12)
    ax.set_xlabel('Months from LTFU date')
    ax.set_title(f'{cause} mortality', fontsize=11, loc='left', fontweight='bold')
    ax.grid(alpha=0.25); ax.legend(loc='upper left', frameon=False, fontsize=9)
    for s_ in ('top','right'): ax.spines[s_].set_visible(False)

plot_curves(axes[0], 'TB-attributable')
plot_curves(axes[1], 'Non-TB')
axes[0].set_ylabel('Cumulative mortality (%)')
fig.suptitle('Cause-specific late-window mortality in LTFU patients, by 6-month return status\n'
             '(hybrid TB attribution; 6-month landmark; follow-up 6–24 months from LTFU)',
             fontsize=12.5, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'wrote {OUT_PNG}')
