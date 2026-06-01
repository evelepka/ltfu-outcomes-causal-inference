#!/usr/bin/env python3
"""Two-panel return-to-care figure.
A: risk-set matched cohort, mortality from TIME OF RETURN, returners vs matched disengaged controls (association).
B: g-computation counterfactual (§10), returners observed vs had-they-stayed-on-treatment (preventable cause).
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np
from lifelines import KaplanMeierFitter, CoxPHFitter
np.random.seed(42)

B="/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper"
d=pd.read_csv(B+"/ITT_Analysis/data/itt_cohort.csv",low_memory=False)
d['best_start']=pd.to_datetime(d.best_start,errors='coerce'); d['end_date']=pd.to_datetime(d.end_date,errors='coerce')
d['tx_yrs']=(d.end_date-d.best_start).dt.days/365.25
H=2.0
RET,CON,CF="#c0392b","#2c7fb8","#2ca02c"
fig,(axA,axB)=plt.subplots(1,2,figsize=(12,5.0))

# ================= Panel A: risk-set matched cohort (from time of return) =================
l=d[d.itt_group=='Loss to follow-up'].copy()
l['agegrp']=pd.cut(l.age_tb,[0,25,45,65,200],right=False).astype(str)
l['stratum']=l.agegrp+'|'+l.hiv_aids.astype(str)+'|'+l.hosp_admission.astype(str)
cols=['sinan_clean','time_d','event_d','time_rn','event_rn']
pools={k:{c:v[c].to_numpy() for c in cols} for k,v in l.groupby('stratum')}
rows=[]; K=5; setid=0
for t in l[l.event_rn==1].itertuples():
    tr=t.time_rn; p=pools[t.stratum]
    idx=np.where((p['time_d']>tr)&((p['event_rn']==0)|(p['time_rn']>tr))&(p['sinan_clean']!=t.sinan_clean))[0]
    if idx.size==0: continue
    sel=np.random.choice(idx,min(K,idx.size),replace=False)
    td=t.time_d if t.event_d==1 else np.inf
    rows.append((t.sinan_clean,1,max(min(td,tr+H)-tr,0.5/365.25),1 if td<=tr+H else 0))
    for j in sel:
        c_td=p['time_d'][j] if p['event_d'][j]==1 else np.inf
        c_tr=p['time_rn'][j] if p['event_rn'][j]==1 else np.inf
        if c_td<=c_tr and c_td<=tr+H: cev,cfu=1,c_td-tr
        else: cev,cfu=0,min(c_tr,tr+H)-tr
        rows.append((p['sinan_clean'][j],0,max(cfu,0.5/365.25),cev))
    setid+=1
m=pd.DataFrame(rows,columns=['subjid','grp','fu','event'])
kmf=KaplanMeierFitter()
for g,color,lab in [(1,RET,'Returned to care'),(0,CON,'Remained LTFU (matched)')]:
    sub=m[m.grp==g]; kmf.fit(sub.fu,sub.event); s=kmf.survival_function_
    y=(1-s.iloc[:,0].values)*100; y24=(1-kmf.predict(H))*100
    axA.step(s.index.values*12,y,where='post',color=color,lw=2.4,label=f'{lab} (n={len(sub):,})')
    axA.text(24.4,y24,f'{y24:.1f}%',color=color,va='center',fontsize=9.5,fontweight='bold')
cph=CoxPHFitter().fit(m[['fu','event','grp','subjid']],'fu','event',cluster_col='subjid')
hr=np.exp(cph.params_['grp']); lo,hiu=np.exp(cph.confidence_intervals_.loc['grp'].values)
axA.text(11,5.0,f'HR {hr:.2f}\n(95% CI {lo:.2f}–{hiu:.2f})',fontsize=11,fontweight='bold',color='#333333')
axA.set_xlim(0,27); axA.set_ylim(0,11); axA.set_xticks([0,6,12,18,24])
axA.set_xlabel('Months since return'); axA.set_ylabel('Cumulative all-cause mortality (%)')
axA.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(0,1.0))
axA.text(-0.02,1.04,'A',transform=axA.transAxes,fontsize=15,fontweight='bold')
for s_ in ('top','right'): axA.spines[s_].set_visible(False)

# ================= Panel B: g-computation counterfactual (from disengagement) =================
COV=['age_tb','sex','hiv_aids','hosp_admission','clinical_clean','homelessness','alcohol','drug_use','dot_status']
m3=d[(d.itt_group=='Loss to follow-up')&(d.tx_yrs>2/12)&(d.tx_yrs<=3/12)].copy()
ret=m3[(m3.event_rn==1)&(m3.time_rn<=0.5)].copy()
ret['t']=np.minimum(ret.time_d,H); ret['ev']=((ret.event_d==1)&(ret.time_d<=H)).astype(int)
ox=d[(d.itt_group!='Loss to follow-up')&(d.tx_yrs>3/12)].copy()
ox['tl']=ox.time_d_tx-3/12; ox=ox[ox.tl>0]; ox['t']=np.minimum(ox.tl,H); ox['ev']=((ox.event_d==1)&(ox.tl<=H)).astype(int)
oxd=pd.get_dummies(ox[['t','ev']+COV],columns=[c for c in COV if c!='age_tb'],drop_first=True,dummy_na=False)
oxd=oxd.dropna().astype({c:float for c in oxd.columns if oxd[c].dtype==bool})
cphB=CoxPHFitter(penalizer=0.001).fit(oxd,'t','ev')
rxd=pd.get_dummies(ret[['t','ev']+COV],columns=[c for c in COV if c!='age_tb'],drop_first=True,dummy_na=False)
for c in oxd.columns:
    if c not in rxd.columns: rxd[c]=0
rxd=rxd[oxd.columns].dropna()
sf=cphB.predict_survival_function(rxd.drop(columns=['t','ev'])); xc=sf.index.values; yc=(1-sf.mean(axis=1).values)*100
ko=KaplanMeierFitter().fit(rxd.t,rxd.ev); so=ko.survival_function_; obs_y=(1-so.iloc[:,0].values)*100
axB.step(so.index.values*12,obs_y,where='post',color=RET,lw=2.4,label=f'Observed (returned; n={len(rxd):,})')
axB.plot(xc*12,yc,color=CF,lw=2.4,ls='--',label='Counterfactual: had stayed on treatment')
axB.text(24.4,obs_y[-1],f'{obs_y[-1]:.1f}%',color=RET,va='center',fontsize=9.5,fontweight='bold')
axB.text(24.4,yc[-1],f'{yc[-1]:.1f}%',color=CF,va='center',fontsize=9.5,fontweight='bold')
axB.set_xlim(0,27); axB.set_ylim(0,11); axB.set_xticks([0,6,12,18,24])
axB.set_xlabel('Months since LTFU'); axB.set_ylabel('Cumulative all-cause mortality (%)')
axB.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(0,1.0))
axB.text(-0.02,1.04,'B',transform=axB.transAxes,fontsize=15,fontweight='bold')
for s_ in ('top','right'): axB.spines[s_].set_visible(False)

fig.tight_layout(w_pad=3)
fig.savefig(B+"/ITT_Analysis/results/Figure_3_return_pathway.png",dpi=200,bbox_inches='tight',facecolor='white')
print(f"A: returned 24mo {(1-kmf.fit(m[m.grp==1].fu,m[m.grp==1].event).predict(H))*100:.1f}% vs HR {hr:.2f}")
print(f"B: observed {obs_y[-1]:.1f}% vs counterfactual {yc[-1]:.1f}%")
print("wrote /tmp/fig_return_2panel.png")
