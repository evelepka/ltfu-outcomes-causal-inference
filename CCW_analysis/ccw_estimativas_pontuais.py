"""Recompute ALL point estimates with the CORRECT arm-specific IPCW (full covariates)."""
import pandas as pd, numpy as np, json
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
import warnings; warnings.filterwarnings('ignore')
SCRATCH = Path("/private/tmp/claude-501/-Users-evelynlepkadelima-repos-outcomes-after-tb-abandonment--claude-worktrees-serene-lewin-0c454f/ae5d64de-0097-404e-9411-b4f453f7faea/scratchpad")
tl0 = pd.read_csv(SCRATCH/"ccw_timeline.csv", low_memory=False)
MONTH,HORIZON=30.4,24
COVS=['age_group','sex','hiv_aids','homelessness','hosp_admission','drug_use','clinical_clean','dot_status']
def prep(tl):
    tl=tl.copy()
    tl['m_dis']=np.floor(np.maximum(tl['t_disengage'],0)/MONTH).clip(upper=HORIZON)
    tl['m_end']=np.floor(np.maximum(tl['t_txend'],0)/MONTH).clip(upper=HORIZON)
    tl['m_death']=np.where(tl['died']==1,np.floor(np.maximum(tl['t_death'],0)/MONTH),np.nan)
    return tl
def build(dat,kind):
    d=dat.copy(); md=d['m_death']; alive=lambda mev: md.isna()|(md>mev)
    if kind=='remain':
        d['m_dev']=np.where(d['is_ltfu'],np.where(alive(d['m_dis']),d['m_dis'],np.inf),np.inf)
    else:
        months=kind[1]; dis_in=d['is_ltfu']&d['m_dis'].isin(months)
        dev_c=np.where(alive(d['m_end']),d['m_end'],np.inf)
        dev_w=np.where(d['is_ltfu']&~dis_in&alive(d['m_dis']),d['m_dis'],np.inf)
        mdev=np.minimum(dev_c,dev_w); d['m_dev']=np.where(dis_in,np.inf,mdev)
    d['stop']=np.minimum(d['m_dev'],HORIZON); d['event_month']=np.where(md<d['stop'],md,np.inf)
    nm=np.where(d['event_month']<np.inf,d['event_month']+1,d['stop']); d['nmonths']=np.clip(nm,0,None).astype(int)
    return d[d['nmonths']>0]
def expand(c):
    n=c['nmonths'].values.astype(int); idx=np.repeat(np.arange(len(c)),n)
    L=c.iloc[idx].reset_index(drop=True); L['month']=np.concatenate([np.arange(k) for k in n])
    L['death']=((L['event_month']<np.inf)&(L['month']==L['event_month'])).astype(int)
    L['dev_next']=((L['m_dev']<HORIZON)&(L['month']==L['stop']-1)&(L['death']==0)).astype(int)
    L['rid']=idx; return L
def wr(L):  # arm-specific IPCW + weighted cumulative incidence
    for c in COVS: L[c]=L[c].astype(str)
    L['mS']=L['month'].astype(str)
    if L['dev_next'].sum()>5:
        enc=ColumnTransformer([('oh',OneHotEncoder(handle_unknown='ignore',sparse_output=True),COVS+['mS'])])
        encm=ColumnTransformer([('oh',OneHotEncoder(handle_unknown='ignore',sparse_output=True),['mS'])])
        Xd=enc.fit_transform(L); Xn=encm.fit_transform(L); y=L['dev_next'].values
        pdv=LogisticRegression(max_iter=150,C=1e6).fit(Xd,y).predict_proba(Xd)[:,1]
        pnv=LogisticRegression(max_iter=150,C=1e6).fit(Xn,y).predict_proba(Xn)[:,1]
        t=pd.DataFrame({'g':L['rid'].values,'d':1-pdv,'n':1-pnv})
        t['cd']=t.groupby('g')['d'].cumprod(); t['cn']=t.groupby('g')['n'].cumprod(); L['sw']=(t['cn']/t['cd']).values
    else: L['sw']=1.0
    lo,hi=L['sw'].quantile([0.005,0.995]); L['swt']=L['sw'].clip(lo,hi)
    g=L.groupby('month').apply(lambda x: np.average(x['death'],weights=x['swt'])); ci=1-(1-g).cumprod()
    r6=ci.loc[ci.index<=5].iloc[-1] if (ci.index<=5).any() else np.nan
    return ci.iloc[-1], r6

tl=prep(tl0)
# remain arm once (shared)
r24_0, r6_0 = wr(expand(build(tl,'remain')))
CONTRASTS={'primary_m0':[0],'intensive_m1_2':[1,2],'continuation_m3_6':[3,4,5,6],
           'grad_m1':[1],'grad_m2':[2],'grad_m3':[3],'grad_m4':[4],'grad_m5':[5],'grad_m6':[6]}
point={}
print("=== CORRECTED POINT ESTIMATES (arm-specific IPCW, full covs) ===\n")
print(f"  Reference (remain engaged): 6mo {r6_0*100:.2f}%  24mo {r24_0*100:.2f}%\n")
for name,mo in CONTRASTS.items():
    r24_1,r6_1=wr(expand(build(tl,('w',mo))))
    point[name]=dict(rr24=r24_1/r24_0, rd24=100*(r24_1-r24_0), rr6=r6_1/r6_0, rd6=100*(r6_1-r6_0),
                     risk24=r24_1*100, risk6=r6_1*100)
    print(f"  {name:18s}: 6mo RR={r6_1/r6_0:.2f} (RD {100*(r6_1-r6_0):+.2f}pp) | 24mo RR={r24_1/r24_0:.2f} (RD {100*(r24_1-r24_0):+.2f}pp)")
point['_remain']=dict(risk24=r24_0*100, risk6=r6_0*100)
json.dump(point, open(SCRATCH/"ccw_point_correct.json","w"), indent=2)
print("\nSaved corrected point estimates.")
