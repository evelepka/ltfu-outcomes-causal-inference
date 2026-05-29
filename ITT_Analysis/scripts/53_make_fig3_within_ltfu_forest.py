"""Regenerate Figure 3 (within-LTFU multivariable forest plot) from new MI results.
Two-panel: Cox mortality (left) + Fine-Gray retreatment (right).
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re

CSV = '/Users/jasonandrews/Google Drive/My Drive/Abandonment Paper/ITT_Analysis/results/multivariable_results_mi_cc.csv'
OUT_PNG = '/Users/jasonandrews/Google Drive/My Drive/Abandonment Paper/ITT_Analysis/results/Figure_4_within_ltfu_forest.png'
OUT_PDF = '/Users/jasonandrews/Google Drive/My Drive/Abandonment Paper/ITT_Analysis/results/Figure_4_within_ltfu_forest.pdf'

df = pd.read_csv(CSV)

def parse_est(s):
    """Parse 'HR (lo-hi)' -> (HR, lo, hi)."""
    m = re.match(r'([\d.]+)\s*\(([\d.]+)[-–]([\d.]+)\)', str(s))
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3))
    return np.nan, np.nan, np.nan

# Map raw term names to display labels and groups
TERM_MAP = {
    "age_group25-44": ("Age 25-44 (ref. 15-24)", "Age"),
    "age_group45-64": ("Age 45-64", "Age"),
    "age_group≥65": ("Age ≥65", "Age"),
    "sexMale": ("Male (ref. female)", "Sex"),
    "race_cleanBlack or Mixed": ("Black or Mixed (ref. White)", "Race"),
    "race_cleanOther": ("Other/Indigenous", "Race"),
    "edu_clean≤ 7 years": ("≤7 years (ref. ≥12 years)", "Education"),
    "edu_clean8 - 11 years": ("8-11 years", "Education"),
    "hiv_aidsPositive": ("HIV positive (ref. negative)", "Clinical"),
    "diabetesYes": ("Diabetes", "Clinical"),
    "alcoholYes": ("Alcohol use", "Behavioural"),
    "drug_useYes": ("Drug use", "Behavioural"),
    "incarceratedYes": ("Incarcerated", "Social"),
    "homelessnessYes": ("Homeless", "Social"),
    "hosp_admissionYes": ("Hospitalised at diagnosis", "Disease"),
    "clinical_cleanExtrapulmonary": ("Extrapulmonary (ref. pulmonary)", "Disease"),
    "clinical_cleanPulmonary and Extrapulmonary or disseminated": ("Pulm+Extra/disseminated", "Disease"),
    "dot_statusYes": ("DOT (ref. no DOT)", "Treatment"),
    "tx_month_grp< 2 months": ("LTFU <2 mo (ref. ≥4 mo)", "Timing"),
    "tx_month_grp2 to <4 months": ("LTFU 2-4 mo", "Timing"),
}

# Get the two adjusted MI models
mort = df[df.model=='Cox_Death_Adjusted_MI'].copy()
retr = df[df.model=='FG_Retr_Adjusted_MI'].copy()

def expand(d):
    rows = []
    for _, r in d.iterrows():
        hr, lo, hi = parse_est(r['estimate'])
        rows.append({'term': r['term'], 'hr': hr, 'lo': lo, 'hi': hi, 'p': r['p_value']})
    return pd.DataFrame(rows)

mort = expand(mort)
retr = expand(retr)

# Order: by group then label
order = list(TERM_MAP.keys())
mort = mort.set_index('term').reindex(order).reset_index()
retr = retr.set_index('term').reindex(order).reset_index()

# Plot: tighter inter-column space, larger fonts
fig, axes = plt.subplots(1, 2, figsize=(13, 9), sharey=True,
                          gridspec_kw={'wspace': 0.05})
y = np.arange(len(order))[::-1]  # top to bottom
labels = [TERM_MAP[t][0] for t in order]

for ax, dat, title in [(axes[0], mort, 'A. Mortality (Cox MI-pooled)'),
                        (axes[1], retr, 'B. Retreatment (Fine-Gray MI-pooled)')]:
    ax.errorbar(dat['hr'], y,
                xerr=[dat['hr']-dat['lo'], dat['hi']-dat['hr']],
                fmt='o', color='black', ecolor='black', capsize=3, lw=1.5, markersize=6)
    ax.axvline(1, color='gray', ls='--', lw=1.0)
    ax.set_xscale('log')
    ax.set_xlim(0.2, 8)
    xticks = [0.25, 0.5, 1, 2, 4, 8]
    ax.set_xticks(xticks); ax.set_xticklabels([str(x) for x in xticks], fontsize=12)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('Adjusted hazard ratio (95% CI)', fontsize=13)
    ax.set_title(title, fontsize=14, loc='left', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=12)

plt.subplots_adjust(left=0.27, right=0.98, top=0.95, bottom=0.07, wspace=0.05)
plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_PDF, bbox_inches='tight', facecolor='white')
print(f'Wrote {OUT_PNG}')
print(f'Wrote {OUT_PDF}')
