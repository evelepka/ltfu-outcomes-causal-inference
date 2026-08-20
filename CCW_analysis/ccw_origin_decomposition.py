#!/usr/bin/env python3
"""Is the CCW's month-1 sign flip caused by the GRACE WINDOW (origin at day 0)
or by the WEIGHTING? Decompose by holding one fixed and varying the other.

This is the experiment behind handoff section 6, and it is what settled the
design split: CCW for the main effect and subgroups, rolling landmark for timing.

B and C are DIAGNOSTIC ONLY and must not be promoted into the manuscript. They
put the origin at a MONTH BOUNDARY rather than at each patient's own declaration
date, which is exactly the monthly-landmark defect ADR-0004 fixed -- for a
month-1 patient who disengaged on day 5, declaration is ~day 35 but the clock
starts at day 30.4. B/C is the superseded monthly landmark with weights bolted
on. Use 45b, which aligns on each patient's own declaration date.

  A  CCW as-is        origin = month boundary m0-1, 1-month grace, clones + IPCW
  B  landmark-aligned origin = month boundary m0,   no grace, per-protocol + IPCW+IPTW
  C  landmark-aligned origin = month boundary m0,   no grace, ITT + IPTW
  C0 same as C, unweighted (crude), to show how much the weights do

If B and C flip positive, origin alignment is the cause and a continuous CCW
would indeed remove the artifact -- by becoming the landmark.
Month 5 is run as a control: the mechanism predicts little change there.
"""
import importlib.util
import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

BASE = Path(os.environ.get("TB_ABANDONMENT_ROOT",
                           Path(__file__).resolve().parent.parent))
spec = importlib.util.spec_from_file_location("ccw", BASE / "CCW_analysis/ccw_v3.py")
ccw = importlib.util.module_from_spec(spec); spec.loader.exec_module(ccw)
H = 60; ccw.HORIZON_M = H; MD = ccw.MONTH_DAYS

lookup = ccw.build_cause_lookup(verbose=False)
tl = ccw.load_timeline(sorted(ccw.MI_DIR.glob("imp_*.csv"))[0], lookup, verbose=False)
tl, Xpat = ccw.attach_patterns(tl, ccw.COVS)
t_dis = tl["t_dis"].to_numpy(); t_d = tl["t_death"].to_numpy(); t_ad = tl["t_admin"].to_numpy()
is_t = tl[ccw.CAUSES["all_cause"][0]].to_numpy().astype(bool)
ad = tl["__any_death__"].to_numpy().astype(bool)
t_c = np.where(ad & ~is_t, t_d, np.inf); t_e = np.where(is_t, t_d, np.inf)
with np.errstate(invalid="ignore"):
    m_dis = np.where(np.isfinite(t_dis), np.floor(t_dis / MD), np.inf)
    m_ev = np.where(np.isfinite(t_e), np.floor(t_e / MD), np.inf)
    m_comp = np.where(np.isfinite(t_c), np.floor(t_c / MD), np.inf)
    m_adm = np.floor(t_ad / MD)
pats = tl["pat"].to_numpy()


def arm(mask, m0, per_protocol):
    rE, rC, rA = m_ev - m0, m_comp - m0, m_adm - m0
    r_dev = (m_dis - m0) if per_protocol else np.full(len(m_dis), np.inf)
    r_cens = np.minimum(np.minimum(r_dev, rC), np.minimum(rA, float(H)))
    has = mask & (rE < r_cens) & (rE < float(H)) & np.isfinite(rE)
    nm = np.clip(np.where(has, rE + 1, r_cens), 0, H)
    nm = np.where(mask, nm, 0).astype(np.int32)
    dev = (mask & ~has & np.isfinite(r_dev) & (r_dev <= rA) & (r_dev <= rC)
           & (r_dev < float(H)))
    k = nm > 0
    d = {"nmonths": nm[k], "event_month": np.where(has, rE, np.inf)[k],
         "dev_binds": dev[k], "pat": pats[k], "at_risk_mode": "all", "T": None}
    return d, k


def risk_landmark(m0, per_protocol, use_iptw=True):
    """Origin at month boundary m0; exposed = disengaged during month m0-1."""
    alive = (m_ev >= m0) & (m_comp >= m0) & (m_adm > m0)
    exp_m = alive & (m_dis == m0 - 1)          # already disengaged at the origin
    cmp_m = alive & (m_dis >= m0)              # still in care at the origin
    elig = exp_m | cmp_m
    A = exp_m[elig].astype(int)
    if use_iptw:
        X = Xpat[pats[elig]]
        ps = LogisticRegression(max_iter=200, C=1.0).fit(X, A).predict_proba(X)[:, 1]
        pm = A.mean()
        w_el = np.where(A == 1, pm / np.clip(ps, 1e-6, None),
                        (1 - pm) / np.clip(1 - ps, 1e-6, None))
    else:
        w_el = np.ones(elig.sum())
    w_all = np.zeros(len(elig)); w_all[np.where(elig)[0]] = w_el
    out = {}
    for lbl, mask, pp in (("exp", exp_m, False), ("cmp", cmp_m, per_protocol)):
        d, k = arm(mask, m0, pp)
        rows = ccw.add_ipcw(ccw.expand(d), Xpat, verbose=False)
        rows["swt"] = rows["swt"] * np.repeat(w_all[k], d["nmonths"])
        out[lbl] = ccw.weighted_risk(rows)
    r1 = ccw.risk_at(out["exp"], H); r0 = ccw.risk_at(out["cmp"], H)
    return r1, r0, int(exp_m.sum()), int(cmp_m.sum())


def risk_ccw(M):
    """The existing sequential CCW: origin at m0 = M-1, one-month grace."""
    m0 = M - 1
    elig = (m_ev >= m0) & (m_comp >= m0) & (m_adm > m0) & (m_dis >= m0)
    here = m_dis == m0
    o = {}
    for a in ("disengage", "remain"):
        if a == "disengage":
            r_dev = np.where(here, np.inf, 1.0); mode, T_ = "window_end", 1
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
        d = {"nmonths": nm[k], "event_month": np.where(has, rE, np.inf)[k],
             "dev_binds": dev[k], "pat": pats[k], "at_risk_mode": mode, "T": T_}
        o[a] = ccw.weighted_risk(ccw.add_ipcw(ccw.expand(d), Xpat, verbose=False))
    return (ccw.risk_at(o["disengage"], H), ccw.risk_at(o["remain"], H),
            int((elig & here).sum()), int(elig.sum()))


print(f"All contrasts at a {H}-month horizon, 1 imputation.\n")
hdr = (f"{'month':>5} {'variant':>34} {'risk exp%':>10} {'risk cmp%':>10} "
       f"{'RD pp':>8} {'RR':>6} {'n exp':>7}")
for M in (1, 5):
    print("=" * len(hdr)); print(hdr); print("-" * len(hdr))
    specs = [
        ("A  CCW as-is (grace, origin day 0)", lambda: risk_ccw(M)),
        ("B  landmark-aligned, per-protocol",  lambda: risk_landmark(M, True, True)),
        ("C  landmark-aligned, ITT",           lambda: risk_landmark(M, False, True)),
        ("C0 landmark-aligned, ITT, unweighted", lambda: risk_landmark(M, False, False)),
    ]
    for lbl, fn in specs:
        r1, r0, ne, _ = fn()
        print(f"{M:>5} {lbl:>34} {100*r1:>10.2f} {100*r0:>10.2f} "
              f"{100*(r1-r0):>+8.2f} {r1/r0 if r0>0 else float('nan'):>6.2f} {ne:>7,}")
    print()
