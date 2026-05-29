"""
55b. Love plot of landmark SMDs
================================
Visualises ITT_Analysis/results/landmark_baseline_smd.csv as a Love plot:
one row per variable, x-axis is SMD (LTFU vs on-treatment), points
coloured by landmark month (1, 3, 6). Two panels (severity vs social).
Dashed vertical lines at ±0.1 mark the conventional imbalance threshold.

Positive SMD = LTFU arm has more of the indicated level (e.g. more "Yes"
to hospitalisation) than the on-treatment comparator at that landmark.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = Path('/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/'
            '.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper')
SRC = BASE/'ITT_Analysis'/'results'/'landmark_baseline_smd.csv'
OUT_PNG = BASE/'ITT_Analysis'/'results'/'landmark_smd_loveplot.png'
OUT_PDF = BASE/'ITT_Analysis'/'results'/'landmark_smd_loveplot.pdf'

df = pd.read_csv(SRC)

# Choose ONE informative level per variable for the Love plot.
# (For continuous age, use the synthetic "mean (SD)" row.)
LEVEL_PICK = {
    'age_tb'                 : 'mean (SD)',
    'hiv_aids'               : 'Positive',
    'hosp_admission'         : 'Yes',
    'bac1_clean'             : 'Positive',
    'sputum_culture_clean'   : 'Positive',
    'lab_confirmed_stat'     : 'No',
    'resistance_clean'       : 'Sensitive',
    'clinical_clean'         : 'Extrapulmonary',
    'diabetes'               : 'Yes',
    'other_immuno_condition' : 'Yes',
    'diagnosis_setting'      : 'Emergency / Inpatient',
    'sex'                    : 'Male',
    'race_clean'             : 'Black or Mixed',
    'edu_clean'              : '≤ 7 years',
    'homelessness'           : 'Yes',
    'alcohol'                : 'Yes',
    'drug_use'               : 'Yes',
    'mental_health'          : 'Yes',
    'tobacco_use'            : 'Yes',
    'incarcerated'           : 'Yes',
    'dot_status'             : 'No',
}
LABEL = {
    'age_tb'                 : 'Age (years, continuous)',
    'hiv_aids'               : 'HIV-positive',
    'hosp_admission'         : 'Hospitalised at diagnosis',
    'bac1_clean'             : 'Smear-positive',
    'sputum_culture_clean'   : 'Sputum culture positive',
    'lab_confirmed_stat'     : 'Lab-unconfirmed',
    'resistance_clean'       : 'Drug-susceptible TB',
    'clinical_clean'         : 'Extrapulmonary TB',
    'diabetes'               : 'Diabetes',
    'other_immuno_condition' : 'Other immunocompromise',
    'diagnosis_setting'      : 'Emergency/inpatient dx setting',
    'sex'                    : 'Male',
    'race_clean'             : 'Black or Mixed race',
    'edu_clean'              : '≤7 years education',
    'homelessness'           : 'Homelessness',
    'alcohol'                : 'Alcohol use',
    'drug_use'               : 'Drug use',
    'mental_health'          : 'Mental health condition',
    'tobacco_use'            : 'Tobacco use',
    'incarcerated'           : 'Incarcerated',
    'dot_status'             : 'No DOT',
}

# Order variables within each panel (top of plot = strongest claim of severity / vulnerability)
SEVERITY_ORDER = ['hosp_admission','diagnosis_setting','hiv_aids','clinical_clean',
                  'other_immuno_condition','diabetes','bac1_clean','sputum_culture_clean',
                  'lab_confirmed_stat','resistance_clean','age_tb']
SOCIAL_ORDER   = ['homelessness','drug_use','alcohol','tobacco_use','mental_health',
                  'incarcerated','dot_status','edu_clean','race_clean','sex']

def filter_picked(df_):
    rows=[]
    for v, lvl in LEVEL_PICK.items():
        sub = df_[(df_['variable']==v) & (df_['level']==lvl)]
        rows.append(sub)
    return pd.concat(rows, ignore_index=True)

picked = filter_picked(df)
sev = picked[picked['panel']=='severity'].copy()
soc = picked[picked['panel']=='social'].copy()

COLOURS = {1:'#1f77b4', 3:'#ff7f0e', 6:'#2ca02c'}

fig, axes = plt.subplots(1, 2, figsize=(11.5, 7.5), sharex=True,
                         gridspec_kw={'wspace':0.45})

def draw_panel(ax, panel_df, order, title):
    ax.axvline(0, color='black', lw=0.6)
    for x in (-0.1, 0.1):
        ax.axvline(x, color='grey', lw=0.6, ls='--', alpha=0.7)
    y_positions = {v:i for i,v in enumerate(reversed(order))}
    for v in order:
        for m, color in COLOURS.items():
            row = panel_df[(panel_df['variable']==v) & (panel_df['landmark_m']==m)]
            if row.empty: continue
            smd = row['smd'].iloc[0]
            ax.scatter(smd, y_positions[v], color=color, s=44, zorder=3, edgecolor='black', lw=0.4)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels([LABEL[v] for v in reversed(order)], fontsize=9)
    ax.set_xlabel('Standardised mean difference\n(LTFU − on-treatment, at landmark)', fontsize=10)
    ax.set_title(title, fontsize=11, loc='left', fontweight='bold')
    ax.set_xlim(-0.85, 0.85)
    ax.grid(axis='x', alpha=0.25)
    for spine in ('top','right'): ax.spines[spine].set_visible(False)

draw_panel(axes[0], sev, SEVERITY_ORDER,
           'A. TB-severity / early-mortality predictors')
draw_panel(axes[1], soc, SOCIAL_ORDER,
           'B. Social vulnerability')

# annotation on direction
axes[0].text(-0.83, len(SEVERITY_ORDER)-0.5, '← lower in LTFU', fontsize=8.5, color='dimgrey', ha='left')
axes[0].text( 0.83, len(SEVERITY_ORDER)-0.5, 'higher in LTFU →', fontsize=8.5, color='dimgrey', ha='right')

# legend
handles = [Line2D([0],[0], marker='o', color='w', label=f'Landmark month {m}',
                  markerfacecolor=c, markeredgecolor='black', markersize=8)
           for m,c in COLOURS.items()]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=10)

fig.suptitle('Baseline characteristics at landmark: LTFU arm vs on-treatment comparator',
             fontsize=12.5, fontweight='bold', y=0.99)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print('wrote', OUT_PNG)
print('wrote', OUT_PDF)
