#!/usr/bin/env python3
"""Bootstrap confidence intervals for the CCW estimates at a 60-month horizon.

Supplies the intervals for the MAIN EFFECT (nested "disengage by month 6") and
for the per-month sequential trials -- handoff sections 2a and 7.

Resamples PATIENTS, not clone-months: a resample changes how often each patient
appears, then both arms are rebuilt and the IPCW refitted inside every replicate.
Bootstrapping the expanded person-month rows would treat correlated rows from one
patient as independent and understate the variance badly. An imputation is drawn
per replicate too, so imputation uncertainty is carried rather than ignored.

Also emits the PAIRED month-to-month differences. Overlapping marginal intervals
are not a test of whether the timing gradient exists; the months share their
comparator pool, so the difference has to be taken within replicate.

  B=3   HORIZON=60 python3 CCW_analysis/ccw_seq_bootstrap.py   # timing check
  B=300 HORIZON=60 python3 CCW_analysis/ccw_seq_bootstrap.py   # as run for the handoff
"""
import importlib.util, os, time
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(os.environ.get("TB_ABANDONMENT_ROOT",
                           Path(__file__).resolve().parent.parent))
spec = importlib.util.spec_from_file_location("ccw", BASE / "CCW_analysis/ccw_v3.py")
ccw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccw)

H = int(os.environ.get("HORIZON", 60))
B = int(os.environ.get("B", 3))
ccw.HORIZON_M = H
MD = ccw.MONTH_DAYS
MONTHS = [1, 2, 3, 4, 5, 6]
rng = np.random.default_rng(2026)

# A low-B timing check must NOT overwrite a real run's outputs: the draws file is
# what the CI tables and ccw_timing_heterogeneity.py read. Learned the hard way.
SUF = "" if B >= 50 else "_smoke"


def prep(tl):
    t_dis = tl["t_dis"].to_numpy(); t_death = tl["t_death"].to_numpy()
    t_admin = tl["t_admin"].to_numpy()
    is_t = tl[ccw.CAUSES["all_cause"][0]].to_numpy().astype(bool)
    anyd = tl["__any_death__"].to_numpy().astype(bool)
    t_comp = np.where(anyd & ~is_t, t_death, np.inf)
    t_ev = np.where(is_t, t_death, np.inf)
    with np.errstate(invalid="ignore"):
        return (np.where(np.isfinite(t_dis), np.floor(t_dis / MD), np.inf),
                np.where(np.isfinite(t_ev), np.floor(t_ev / MD), np.inf),
                np.where(np.isfinite(t_comp), np.floor(t_comp / MD), np.inf),
                np.floor(t_admin / MD))


def seq_rd(tl, Xpat, M, m_dis, m_ev, m_comp, m_adm):
    m0 = M - 1
    elig = (m_ev >= m0) & (m_comp >= m0) & (m_adm > m0) & (m_dis >= m0)
    dis_here = m_dis == m0
    out = {}
    for a in ("disengage", "remain"):
        if a == "disengage":
            r_dev = np.where(dis_here, np.inf, 1.0); mode, T_ = "window_end", 1
        else:
            r_dev = m_dis - m0; mode, T_ = "all", None
        rE, rC, rA = m_ev - m0, m_comp - m0, m_adm - m0
        r_cens = np.minimum(np.minimum(r_dev, rC), np.minimum(rA, float(H)))
        has = elig & (rE < r_cens) & (rE < float(H)) & np.isfinite(rE)
        nm = np.clip(np.where(has, rE + 1, r_cens), 0, H)
        nm = np.where(elig, nm, 0).astype(np.int32)
        dev = (elig & ~has & np.isfinite(r_dev) & (r_dev <= rA) & (r_dev <= rC)
               & (r_dev < float(H)))
        k = nm > 0
        if not k.any():
            return np.nan, np.nan
        d = {"nmonths": nm[k], "event_month": np.where(has, rE, np.inf)[k],
             "dev_binds": dev[k], "pat": tl["pat"].to_numpy()[k],
             "at_risk_mode": mode, "T": T_}
        rows = ccw.add_ipcw(ccw.expand(d), Xpat, verbose=False)
        out[a] = ccw.weighted_risk(rows)
    r1 = ccw.risk_at(out["disengage"], H); r0 = ccw.risk_at(out["remain"], H)
    return (r1 / r0 if r0 > 0 else np.nan), 100 * (r1 - r0)


def nested_rd(tl, Xpat, T=6):
    out = {}
    for a in ("disengage", "remain"):
        arm = ccw.build_arm(tl, a, T=T, cause="all_cause")
        rows = ccw.add_ipcw(ccw.expand(arm), Xpat, verbose=False)
        out[a] = ccw.weighted_risk(rows)
    r1 = ccw.risk_at(out["disengage"], H); r0 = ccw.risk_at(out["remain"], H)
    return (r1 / r0 if r0 > 0 else np.nan), 100 * (r1 - r0)


lookup = ccw.build_cause_lookup(verbose=False)
imps = sorted(ccw.MI_DIR.glob("imp_*.csv"))
base = ccw.load_timeline(imps[0], lookup, verbose=False)
n = len(base)

recs = []
t0 = time.time()
for b in range(B + 1):                       # b=0 is the point estimate
    imp = imps[0] if b == 0 else imps[rng.integers(len(imps))]
    tl = base if b == 0 else ccw.load_timeline(imp, lookup, verbose=False)
    if b > 0:
        idx = rng.integers(0, n, n)          # resample PATIENTS with replacement
        tl = tl.iloc[idx].reset_index(drop=True)
    tl, Xpat = ccw.attach_patterns(tl, ccw.COVS)
    md, me, mc, ma = prep(tl)
    rr, rd = nested_rd(tl, Xpat)
    recs.append({"b": b, "which": "nested_T6", "rr": rr, "rd": rd})
    for M in MONTHS:
        rr, rd = seq_rd(tl, Xpat, M, md, me, mc, ma)
        recs.append({"b": b, "which": f"seq_M{M}", "rr": rr, "rd": rd})
    if b == 0:
        print(f"  point estimates done in {time.time()-t0:.0f}s")
    elif b % 10 == 0 or b <= 2:
        el = time.time() - t0
        print(f"  replicate {b}/{B}  elapsed {el/60:.1f} min  "
              f"projected total {el/b*B/60:.0f} min", flush=True)

d = pd.DataFrame(recs)
d.to_csv(ccw.OUTDIR / f"ccw_bootstrap_draws_h{H}{SUF}.csv", index=False)
pt = d[d.b == 0].set_index("which")
bs = d[d.b > 0]
rows = []
for w, g in bs.groupby("which"):
    rows.append({
        "estimate": w,
        "rr": pt.loc[w, "rr"], "rr_lo": g.rr.quantile(.025), "rr_hi": g.rr.quantile(.975),
        "rd": pt.loc[w, "rd"], "rd_lo": g.rd.quantile(.025), "rd_hi": g.rd.quantile(.975),
        "B_ok": g.rd.notna().sum()})
out = pd.DataFrame(rows)
order = ["nested_T6"] + [f"seq_M{m}" for m in MONTHS]
out["o"] = out.estimate.map({k: i for i, k in enumerate(order)})
out = out.sort_values("o").drop(columns="o")
p = ccw.OUTDIR / f"ccw_bootstrap_h{H}{SUF}.csv"
out.to_csv(p, index=False)

print(f"\n=== CCW at {H} months, patient-level bootstrap, B={B} ===")
print(f"{'estimate':>11} {'RR':>6} {'95% CI':>16} {'RD pp':>8} {'95% CI':>18} {'B':>5}")
print("-" * 74)
for _, r in out.iterrows():
    print(f"{r.estimate:>11} {r.rr:>6.2f} {f'({r.rr_lo:.2f}-{r.rr_hi:.2f})':>16} "
          f"{r.rd:>+8.2f} {f'({r.rd_lo:+.2f} to {r.rd_hi:+.2f})':>18} {int(r.B_ok):>5}")
print(f"\nwrote {p.relative_to(ccw.ROOT)}   total {(time.time()-t0)/60:.1f} min")

# ---- PAIRED contrasts between months. Overlapping marginal CIs are NOT a test;
# the months share patients, so the difference must be taken within replicate.
w = bs.pivot(index="b", columns="which", values="rd")
PAIRS = [("seq_M5", "seq_M2"), ("seq_M5", "seq_M6"), ("seq_M5", "seq_M1"),
         ("seq_M4", "seq_M2"), ("seq_M3", "seq_M2"), ("seq_M6", "seq_M2"),
         ("seq_M2", "seq_M1")]
print("\n=== paired bootstrap differences in RD (pp), within replicate ===")
print(f"{'contrast':>21} {'diff':>7} {'95% CI':>18} {'2-sided p':>10}")
print("-" * 60)
prs = []
for a, b_ in PAIRS:
    dd = (w[a] - w[b_]).dropna()
    lo, hi = dd.quantile(.025), dd.quantile(.975)
    pv = 2 * min((dd <= 0).mean(), (dd >= 0).mean())
    pv = max(pv, 1 / len(dd))
    prs.append({"contrast": f"{a}-{b_}", "diff": pt.loc[a, "rd"] - pt.loc[b_, "rd"],
                "lo": lo, "hi": hi, "p": pv})
    print(f"{a+' vs '+b_:>21} {prs[-1]['diff']:>+7.2f} "
          f"{f'({lo:+.2f} to {hi:+.2f})':>18} {pv:>10.3f}")
pd.DataFrame(prs).to_csv(ccw.OUTDIR / f"ccw_bootstrap_paired_h{H}{SUF}.csv", index=False)

# does the gradient exist at all? spread across months 2-6 per replicate
sp = w[[f"seq_M{m}" for m in (2, 3, 4, 5, 6)]]
rng_pp = (sp.max(axis=1) - sp.min(axis=1)).dropna()
print(f"\nspread across months 2-6 (max-min RD): point "
      f"{sp.loc[:, :].max(axis=1).median() - sp.min(axis=1).median():+.2f} pp; "
      f"bootstrap 95% CI ({rng_pp.quantile(.025):+.2f} to {rng_pp.quantile(.975):+.2f})")
print(f"argmax month across replicates:")
am = sp.idxmax(axis=1).value_counts(normalize=True).sort_index()
for k, v in am.items():
    print(f"   {k}: {100*v:5.1f}% of replicates")
