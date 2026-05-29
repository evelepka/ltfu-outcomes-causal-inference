"""Final, correct multi-source DR-status classifier using only actual DST results.
Distinguishes RR/MDR vs INH-mono using rif and isoniazid signals."""
import pandas as pd
from statsmodels.stats.proportion import proportion_confint as ci

ROOT = '/Users/jasonandrews/Google Drive/My Drive/Abandonment Paper'

src = pd.read_csv(f'{ROOT}/Data/Final_table_cleaned.csv', low_memory=False)

def classify(row):
    res = str(row['resistance']).upper().strip() if pd.notna(row['resistance']) else None
    rif = str(row['rifasens']).strip() if pd.notna(row['rifasens']) else None
    inh = str(row['isonisens']).strip() if pd.notna(row['isonisens']) else None

    rif_resist = (rif == 'Resist')
    inh_resist = (inh == 'Resist')
    mdr_phen   = (res == 'TB MR')           # multidrug-resistant
    sdr_phen   = (res == 'TB R')            # single-drug-resistant (ambiguous: could be R or H)
    sens_phen  = (res == 'SENS')
    rif_sens   = (rif == 'Sens')
    inh_sens   = (inh == 'Sens')

    # Priority order:
    # 1) MDR if MDR phen, or both rif_resist AND inh_resist
    if mdr_phen or (rif_resist and inh_resist):
        return 'MDR-TB'
    # 2) RR (rifampin-resistant, non-MDR) if rif_resist and not INH-resistant
    if rif_resist:
        return 'RR-TB (rifampin-resistant)'
    # 3) INH-mono if isonisens Resist and not rif-resist
    if inh_resist:
        return 'INH-mono resistance'
    # 4) Single-drug resistant phenotypic with no INH/rif signals: ambiguous, group with RR/MDR as conservative
    if sdr_phen:
        return 'Other resistant'
    # 5) Sensitive only if explicit Sens
    if sens_phen or rif_sens or inh_sens:
        return 'Drug-sensitive'
    return 'Not Evaluated'

src['dr_status_v3'] = src.apply(classify, axis=1)
per_pt = src[['sinan_clean','dr_status_v3']].drop_duplicates('sinan_clean')

# Apply to new cohort
itt = pd.read_csv(f'{ROOT}/ITT_Analysis/data/itt_cohort.csv', low_memory=False, usecols=['sinan_clean','itt_group'])
new = itt.merge(per_pt, on='sinan_clean', how='left').fillna({'dr_status_v3':'Not Evaluated'})
new['is_ltfu'] = new.itt_group == 'Loss to follow-up'
print(f'NEW cohort N: {len(new):,}\n')

order = ['Drug-sensitive', 'RR-TB (rifampin-resistant)', 'MDR-TB', 'INH-mono resistance', 'Other resistant', 'Not Evaluated']
print('=== NEW cohort — corrected multi-source DR-status ===')
for grp in order:
    sub = new[new.dr_status_v3 == grp]
    if len(sub) == 0:
        print(f'  {grp:<32}: N=0')
        continue
    n_l = sub.is_ltfu.sum()
    p = 100 * n_l / len(sub)
    lo, hi = ci(n_l, len(sub), method='wilson')
    print(f'  {grp:<32}: N={len(sub):>7,} ({100*len(sub)/len(new):.1f}%)  LTFU={n_l:>5,} ({p:.1f}%, 95% CI {lo*100:.1f}-{hi*100:.1f})')

# Total DST performed
n_dst = (new.dr_status_v3 != 'Not Evaluated').sum()
n_resistant = new.dr_status_v3.isin(['RR-TB (rifampin-resistant)', 'MDR-TB', 'INH-mono resistance', 'Other resistant']).sum()
n_rr_mdr = new.dr_status_v3.isin(['RR-TB (rifampin-resistant)', 'MDR-TB']).sum()
n_inh_mono = (new.dr_status_v3 == 'INH-mono resistance').sum()
print()
print(f'DST yielded a result: {n_dst:,} ({100*n_dst/len(new):.1f}%)')
print(f'  Any resistance: {n_resistant:,} ({100*n_resistant/n_dst:.2f}% of tested)')
print(f'  RR or MDR-TB:   {n_rr_mdr:,} ({100*n_rr_mdr/n_dst:.2f}% of tested)')
print(f'  INH-mono:       {n_inh_mono:,} ({100*n_inh_mono/n_dst:.2f}% of tested)')

# Save lookup
out = per_pt.rename(columns={'dr_status_v3':'dr_status'})
out.to_csv(f'{ROOT}/ITT_Analysis/results/resistance/dr_status_corrected.csv', index=False)
print(f'\nSaved {ROOT}/ITT_Analysis/results/resistance/dr_status_corrected.csv')
