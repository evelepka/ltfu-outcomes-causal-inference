"""
55c. Propensity-score-weighted SMDs at landmark m=3 (MI-pooled)
================================================================
Companion to scripts 55/55b. The unadjusted SMDs compare the marginal
distribution of each covariate between the LTFU arm and the on-treatment
arm. Here we additionally report SMDs in the IPTW-weighted population, so
the "after" SMDs reflect how much imbalance remains once we condition on
the full baseline covariate set via the propensity score.

Pipeline (per imputed dataset, then averaged):
  1. Define LTFU arm at m=3 and on-tx arm at m=3 from `best_start` / `end_date`.
  2. Fit logistic regression for P(LTFU=1 at m=3 | X).
  3. Compute IPTW (ATE) weights: w = 1/PS for LTFU, w = 1/(1-PS) for on-tx.
  4. Compute the crude and weighted SMD for each (variable × level).
  5. Average SMDs across MI datasets (Rubin pooling of point estimates).

Outputs:
  ITT_Analysis/results/landmark_smd_psweighted.csv
  ITT_Analysis/results/landmark_smd_psweighted_loveplot.png/.pdf
"""
from pathlib import Path
import pandas as pd, numpy as np, glob
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
MI_DIR  = BASE/'ITT_Analysis'/'data'/'mi'
OUT_CSV = BASE/'ITT_Analysis'/'results'/'landmark_smd_psweighted.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'landmark_smd_psweighted_loveplot.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'landmark_smd_psweighted_loveplot.pdf'

M_LM_YR = 3.0/12.0

# Variables — full manuscript baseline set
SEVERITY = ['age_tb','hiv_aids','hosp_admission','bac1_clean','sputum_culture_clean',
            'lab_confirmed_stat','resistance_clean','clinical_clean','diabetes',
            'other_immuno_condition','diagnosis_setting']
SOCIAL   = ['sex','race_clean','edu_clean','homelessness','alcohol','drug_use',
            'mental_health','tobacco_use','incarcerated','dot_status']
ALL_VARS = SEVERITY + SOCIAL

def weighted_smd_binary(x, t, w):
    """SMD between treated (t==1) and untreated (t==0) with subject weights w."""
    t = np.asarray(t).astype(bool)
    w = np.asarray(w, dtype=float)
    w1, w0 = w[t], w[~t]
    x1, x0 = np.asarray(x)[t], np.asarray(x)[~t]
    if w1.sum()==0 or w0.sum()==0: return np.nan
    p1 = (x1*w1).sum()/w1.sum()
    p0 = (x0*w0).sum()/w0.sum()
    pbar = (p1+p0)/2
    return np.nan if pbar in (0,1) else (p1-p0)/np.sqrt(pbar*(1-pbar))

def weighted_smd_cont(x, t, w):
    t = np.asarray(t).astype(bool); w = np.asarray(w, dtype=float)
    x = np.asarray(x, dtype=float)
    w1, w0 = w[t], w[~t]; x1, x0 = x[t], x[~t]
    m1 = (x1*w1).sum()/w1.sum(); m0 = (x0*w0).sum()/w0.sum()
    v1 = (w1*(x1-m1)**2).sum()/w1.sum()
    v0 = (w0*(x0-m0)**2).sum()/w0.sum()
    pooled = np.sqrt((v1+v0)/2)
    return np.nan if pooled==0 else (m1-m0)/pooled

# --------------------- run per MI dataset ---------------------
imp_files = sorted(glob.glob(str(MI_DIR/'imp_*.csv')))
print(f'Found {len(imp_files)} MI datasets')

records = []   # one row per (mi, variable, level)

for mi_idx, path in enumerate(imp_files, start=1):
    print(f'\n--- MI dataset {mi_idx}: {Path(path).name} ---')
    d = pd.read_csv(path, low_memory=False)
    d['best_start']=pd.to_datetime(d['best_start'],errors='coerce')
    d['end_date']  =pd.to_datetime(d['end_date'],  errors='coerce')
    d['tx_yrs']    =(d['end_date']-d['best_start']).dt.days/365.25

    # LTFU at m=3: LTFU AND tx_yrs in (2/12, 3/12]
    # On-tx at m=3: non-LTFU AND tx_yrs > 3/12
    ltfu_mask = (d['itt_group']=='Loss to follow-up') & (d['tx_yrs']>2.0/12.0) & (d['tx_yrs']<=3.0/12.0)
    ontx_mask = (d['itt_group']!='Loss to follow-up') & (d['tx_yrs']>3.0/12.0)
    keep = d[ltfu_mask | ontx_mask].copy()
    keep['ltfu'] = ltfu_mask[keep.index].astype(int)
    print(f'  N total={len(keep):,}; LTFU={keep["ltfu"].sum():,}; on-tx={(keep["ltfu"]==0).sum():,}')

    # Build PS features (one-hot encode the categoricals)
    X_df = pd.get_dummies(keep[ALL_VARS],
                           columns=[v for v in ALL_VARS if v!='age_tb'],
                           drop_first=True, dummy_na=False)
    X_df = X_df.astype({c:float for c in X_df.columns if X_df[c].dtype==bool})
    # standardize age for stability
    sc = StandardScaler()
    X_arr = X_df.values
    X_arr[:,0] = sc.fit_transform(X_arr[:,[0]]).ravel()
    y_arr = keep['ltfu'].values

    # Logistic propensity model with light regularisation
    lr = LogisticRegression(max_iter=2000, C=1.0, solver='lbfgs')
    lr.fit(X_arr, y_arr)
    ps = lr.predict_proba(X_arr)[:,1]
    ps = np.clip(ps, 1e-3, 1-1e-3)
    # ATE IPTW weights
    w = np.where(y_arr==1, 1/ps, 1/(1-ps))
    # Truncate at 1/99 percentile
    q1, q99 = np.quantile(w, [0.01, 0.99])
    w = np.clip(w, q1, q99)
    print(f'  PS  range: {ps.min():.3f}–{ps.max():.3f}; weight median={np.median(w):.2f}')

    # Compute SMDs (crude and weighted) per variable level
    n1_raw, n0_raw = int((y_arr==1).sum()), int((y_arr==0).sum())
    for panel, vars_ in (('severity',SEVERITY),('social',SOCIAL)):
        for v in vars_:
            if v=='age_tb':
                smd_crude = weighted_smd_cont(keep[v].values, y_arr, np.ones_like(w))
                smd_adj   = weighted_smd_cont(keep[v].values, y_arr, w)
                records.append({'mi':mi_idx,'panel':panel,'variable':v,'level':'mean',
                                'smd_crude':smd_crude,'smd_adj':smd_adj})
            else:
                levels = sorted(keep[v].dropna().unique())
                for L in levels:
                    ind = (keep[v]==L).astype(int).values
                    smd_crude = weighted_smd_binary(ind, y_arr, np.ones_like(w))
                    smd_adj   = weighted_smd_binary(ind, y_arr, w)
                    records.append({'mi':mi_idx,'panel':panel,'variable':v,'level':str(L),
                                    'smd_crude':smd_crude,'smd_adj':smd_adj})

# --------------------- pool across MI ---------------------
df = pd.DataFrame(records)
pooled = (df.groupby(['panel','variable','level'], as_index=False)
            .agg(smd_crude=('smd_crude','mean'),
                 smd_adj  =('smd_adj','mean')))
pooled.to_csv(OUT_CSV, index=False)
print(f'\nwrote {OUT_CSV}')

# Headline: how many variable-levels were |SMD|>0.1 crude vs adjusted?
nc = (pooled['smd_crude'].abs()>0.10).sum()
na = (pooled['smd_adj'].abs()>0.10).sum()
print(f'\n|SMD|>0.10:  crude={nc}/{len(pooled)},  PS-adjusted={na}/{len(pooled)}')

# --------------------- Love plot crude vs adjusted ---------------------
LEVEL_PICK = {
    'age_tb':'mean', 'hiv_aids':'Positive', 'hosp_admission':'Yes',
    'bac1_clean':'Positive', 'sputum_culture_clean':'Positive',
    'lab_confirmed_stat':'No', 'resistance_clean':'Sensitive',
    'clinical_clean':'Extrapulmonary', 'diabetes':'Yes',
    'other_immuno_condition':'Yes',
    'diagnosis_setting':'Emergency / Inpatient',
    'sex':'Male', 'race_clean':'Black or Mixed', 'edu_clean':'≤ 7 years',
    'homelessness':'Yes', 'alcohol':'Yes', 'drug_use':'Yes',
    'mental_health':'Yes', 'tobacco_use':'Yes', 'incarcerated':'Yes',
    'dot_status':'No',
}
LABEL = {
    'age_tb':'Age (years, continuous)',
    'hiv_aids':'HIV-positive',
    'hosp_admission':'Hospitalised at diagnosis',
    'bac1_clean':'Smear-positive',
    'sputum_culture_clean':'Sputum culture positive',
    'lab_confirmed_stat':'Lab-unconfirmed',
    'resistance_clean':'Drug-susceptible TB',
    'clinical_clean':'Extrapulmonary TB',
    'diabetes':'Diabetes',
    'other_immuno_condition':'Other immunocompromise',
    'diagnosis_setting':'Emergency/inpatient dx setting',
    'sex':'Male',
    'race_clean':'Black or Mixed race',
    'edu_clean':'≤7 years education',
    'homelessness':'Homelessness',
    'alcohol':'Alcohol use',
    'drug_use':'Drug use',
    'mental_health':'Mental health condition',
    'tobacco_use':'Tobacco use',
    'incarcerated':'Incarcerated',
    'dot_status':'No DOT',
}
SEVERITY_ORDER = ['hosp_admission','diagnosis_setting','hiv_aids','clinical_clean',
                  'other_immuno_condition','diabetes','bac1_clean','sputum_culture_clean',
                  'lab_confirmed_stat','resistance_clean','age_tb']
SOCIAL_ORDER   = ['homelessness','drug_use','alcohol','tobacco_use','mental_health',
                  'incarcerated','dot_status','edu_clean','race_clean','sex']

def pick(df_, panel, order):
    rows=[]
    for v in order:
        lvl = LEVEL_PICK[v]
        sub = df_[(df_['panel']==panel) & (df_['variable']==v) & (df_['level']==lvl)]
        if not sub.empty: rows.append(sub.iloc[0])
    return pd.DataFrame(rows)

sev = pick(pooled, 'severity', SEVERITY_ORDER)
soc = pick(pooled, 'social',   SOCIAL_ORDER)

fig, axes = plt.subplots(1, 2, figsize=(12, 7.5), sharex=True,
                         gridspec_kw={'wspace':0.50})

def draw(ax, sub, order, title):
    ax.axvline(0, color='black', lw=0.6)
    for x in (-0.1, 0.1):
        ax.axvline(x, color='grey', lw=0.6, ls='--', alpha=0.6)
    y = {v:i for i,v in enumerate(reversed(order))}
    for v in order:
        r = sub[sub['variable']==v]
        if r.empty: continue
        y_pos = y[v]
        ax.scatter(r['smd_crude'].iloc[0], y_pos, color='#d62728', s=46,
                   edgecolor='black', lw=0.4, zorder=3, label='Crude')
        ax.scatter(r['smd_adj'].iloc[0],   y_pos, color='#1f77b4', s=46,
                   edgecolor='black', lw=0.4, zorder=3, label='PS-weighted')
        # connecting line
        ax.plot([r['smd_crude'].iloc[0], r['smd_adj'].iloc[0]],
                [y_pos, y_pos], color='lightgrey', lw=0.8, zorder=1)
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels([LABEL[v] for v in reversed(order)], fontsize=9)
    ax.set_xlabel('Standardised mean difference\n(LTFU − on-treatment, at m=3)', fontsize=10)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_xlim(-0.85, 0.85)
    ax.grid(axis='x', alpha=0.25)
    for s_ in ('top','right'): ax.spines[s_].set_visible(False)

draw(axes[0], sev, SEVERITY_ORDER, 'A. TB-severity / early-mortality predictors')
draw(axes[1], soc, SOCIAL_ORDER,   'B. Social vulnerability')

handles = [Line2D([0],[0], marker='o', color='w', label='Crude (unadjusted)',
                  markerfacecolor='#d62728', markeredgecolor='black', markersize=8),
           Line2D([0],[0], marker='o', color='w', label='IPTW propensity-weighted',
                  markerfacecolor='#1f77b4', markeredgecolor='black', markersize=8)]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=10)

fig.suptitle('',
             fontsize=12.5, fontweight='bold', y=1.0)
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'wrote {OUT_PNG}')
