"""Optimized bootstrap: remain arm computed once per rep (shared), arm-specific IPCW.
3 primary contrasts (m0, m1-2, m3-6), 200 reps."""
import pandas as pd, numpy as np, json, time
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
import warnings; warnings.filterwarnings('ignore')

SCRATCH = Path("/private/tmp/claude-501/-Users-evelynlepkadelima-repos-outcomes-after-tb-abandonment--claude-worktrees-serene-lewin-0c454f/ae5d64de-0097-404e-9411-b4f453f7faea/scratchpad")
tl0 = pd.read_csv(SCRATCH / "ccw_timeline.csv", low_memory=False)
MONTH, HORIZON = 30.4, 24
COVS = ['age_group','hiv_aids','hosp_admission','homelessness']  # reduced set for weight model (speed)

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
    d['stop']=np.minimum(d['m_dev'],HORIZON)
    d['event_month']=np.where(md<d['stop'],md,np.inf)
    nm=np.where(d['event_month']<np.inf,d['event_month']+1,d['stop'])
    d['nmonths']=np.clip(nm,0,None).astype(int)
    return d[d['nmonths']>0]
def expand(c):
    n=c['nmonths'].values.astype(int); idx=np.repeat(np.arange(len(c)),n)
    L=c.iloc[idx].reset_index(drop=True); L['month']=np.concatenate([np.arange(k) for k in n])
    L['death']=((L['event_month']<np.inf)&(L['month']==L['event_month'])).astype(int)
    L['dev_next']=((L['m_dev']<HORIZON)&(L['month']==L['stop']-1)&(L['death']==0)).astype(int)
    L['rid']=idx; return L
_enc=None;_encm=None
def weighted_risk(L):
    """arm-specific IPCW + weighted cumulative incidence -> return 24mo & 6mo risk"""
    global _enc,_encm
    for c in COVS: L[c]=L[c].astype(str)
    L['mS']=L['month'].astype(str)
    if L['dev_next'].sum()>5:
        enc=ColumnTransformer([('oh',OneHotEncoder(handle_unknown='ignore',sparse_output=True),COVS+['mS'])])
        encm=ColumnTransformer([('oh',OneHotEncoder(handle_unknown='ignore',sparse_output=True),['mS'])])
        Xd=enc.fit_transform(L); Xn=encm.fit_transform(L); y=L['dev_next'].values
        pdv=LogisticRegression(max_iter=100,C=1e6).fit(Xd,y).predict_proba(Xd)[:,1]
        pnv=LogisticRegression(max_iter=100,C=1e6).fit(Xn,y).predict_proba(Xn)[:,1]
        t=pd.DataFrame({'g':L['rid'].values,'d':1-pdv,'n':1-pnv})
        t['cd']=t.groupby('g')['d'].cumprod(); t['cn']=t.groupby('g')['n'].cumprod()
        L['sw']=(t['cn']/t['cd']).values
    else: L['sw']=1.0
    lo,hi=L['sw'].quantile([0.005,0.995]); L['swt']=L['sw'].clip(lo,hi)
    g=L.groupby('month').apply(lambda x: np.average(x['death'],weights=x['swt']))
    ci=1-(1-g).cumprod()
    r6=ci.loc[ci.index<=5].iloc[-1] if (ci.index<=5).any() else np.nan
    return ci.iloc[-1], r6

PRIMARY={'primary_m0':[0],'intensive_m1_2':[1,2],'continuation_m3_6':[3,4,5,6]}
B=200
boot={k:{'rr24':[],'rr6':[]} for k in PRIMARY}
t_start=time.time()
for b in range(B):
    d=tl0.iloc[np.random.RandomState(2000+b).randint(0,len(tl0),len(tl0))].reset_index(drop=True)
    d=prep(d)
    # remain arm ONCE
    r24_0,r6_0 = weighted_risk(expand(build(d,'remain')))
    for k,mo in PRIMARY.items():
        r24_1,r6_1 = weighted_risk(expand(build(d,('w',mo))))
        boot[k]['rr24'].append(r24_1/r24_0); boot[k]['rr6'].append(r6_1/r6_0)
    if (b+1)%10==0:
        json.dump(boot,open(SCRATCH/"ccw_boot_primary.json","w"))
        el=time.time()-t_start
        print(f"rep {b+1}/{B} | {el/(b+1):.1f}s/rep | ETA {el/(b+1)*(B-b-1)/60:.0f} min", flush=True)

# summarize
print("\n=== BOOTSTRAP CI (200 reps) ===", flush=True)
point=json.load(open(SCRATCH/"ccw_point.json"))
summ={}
for k in PRIMARY:
    a=np.array(boot[k]['rr24']); a=a[np.isfinite(a)]
    a6=np.array(boot[k]['rr6']); a6=a6[np.isfinite(a6)]
    summ[k]=dict(rr24=point[k]['rr24'], rr24_lo=np.percentile(a,2.5), rr24_hi=np.percentile(a,97.5),
                 rr6=point[k]['rr6'], rr6_lo=np.percentile(a6,2.5), rr6_hi=np.percentile(a6,97.5),
                 rd24=point[k]['rd24'])
    print(f"  {k:18s}: 24mo RR {point[k]['rr24']:.2f} ({summ[k]['rr24_lo']:.2f}-{summ[k]['rr24_hi']:.2f}) | 6mo RR {point[k]['rr6']:.2f} ({summ[k]['rr6_lo']:.2f}-{summ[k]['rr6_hi']:.2f})", flush=True)
json.dump(summ,open(SCRATCH/"ccw_summary_primary.json","w"),indent=2)
print("BOOTSTRAP DONE", flush=True)
