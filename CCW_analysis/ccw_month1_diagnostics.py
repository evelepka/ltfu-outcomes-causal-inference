#!/usr/bin/env python3
"""Why the CCW's month-1 contrast is not reported. Four diagnostics, one script.

These four ran separately while the 2026-08-19 design decision was being made and
are merged here because they are one argument. They supply handoff sections 4, 5,
the front-loading mechanism in section 3, and the "excluding the
definition-inconsistent records" entry under section 11.

The conclusion they support: month 1 is excluded from the CCW for reasons that do
not reference the sign of the estimate. That matters, because reporting the
landmark's month 1 while suppressing the CCW's would otherwise be choosing a
design by its answer.

  A. WHAT "MONTH M" INDEXES
     m_dis = floor(((closure - 30 d) - treatment start) / 30.4), floored at 0.
     So the index is the inferred DISENGAGEMENT month (last contact), NOT the
     declaration month, which follows ~30 days later. Also counts the records
     whose closure date is under 30 days -- for which the LTFU definition cannot
     have been met and t_dis was floored up from a negative value.

  B. THE IDENTIFIABILITY LIMIT
     Deaths before LTFU could be declared. A patient who stops attending on day 5
     and dies on day 10 is closed as Obito, never Abandono, so the effect of
     disengagement BEHAVIOUR including those deaths is not identifiable from
     TBweb. Neither design recovers it; the CCW merely imputes it through a
     baseline-covariate-only censoring model.

  C. WHERE THE MONTH-1 DEFICIT ARISES
     Arm risks at several horizons for the M=1 and M=5 sequential trials. The
     deficit is established in the first few months and then erodes, and the same
     shape appears at M=5 far smaller -- one mechanism (the comparator is in care,
     so it is still absorbing on-treatment TB deaths while the LTFU patient has
     left), not a month-1 quirk.

  D. DOES EXCLUDING THE DEFINITION-INCONSISTENT RECORDS HELP?
     No. It makes month 1 more negative, because those patients have HIGHER crude
     mortality than the rest of month 1. Months 2-6 do not move at all, since with
     m_dis = 0 they were already ineligible for every later trial. Recorded under
     "considered and rejected" so nobody retries it.

Estimator and weight machinery are imported from ccw_v3.py; the arm construction
matches ccw_sequential_prototype.py.

  python3 CCW_analysis/ccw_month1_diagnostics.py
  HORIZON=60 python3 CCW_analysis/ccw_month1_diagnostics.py
"""
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(os.environ.get("TB_ABANDONMENT_ROOT",
                           Path(__file__).resolve().parent.parent))
spec = importlib.util.spec_from_file_location("ccw", BASE / "CCW_analysis/ccw_v3.py")
ccw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccw)

H = int(os.environ.get("HORIZON", 60))
ccw.HORIZON_M = H
MD = ccw.MONTH_DAYS

lookup = ccw.build_cause_lookup(verbose=False)
imp = sorted(ccw.MI_DIR.glob("imp_*.csv"))[0]
tl_full = ccw.load_timeline(imp, lookup, verbose=False)
print(f"horizon {H} months, 1 imputation ({imp.name})\n")


def month_arrays(tl):
    """Integer-month clocks. Convert days -> months exactly once (see the month
    arithmetic note in the analysis-conventions rule)."""
    t_dis = tl["t_dis"].to_numpy()
    t_d = tl["t_death"].to_numpy()
    t_ad = tl["t_admin"].to_numpy()
    is_t = tl[ccw.CAUSES["all_cause"][0]].to_numpy().astype(bool)
    ad = tl["__any_death__"].to_numpy().astype(bool)
    t_c = np.where(ad & ~is_t, t_d, np.inf)
    t_e = np.where(is_t, t_d, np.inf)
    with np.errstate(invalid="ignore"):
        return (np.where(np.isfinite(t_dis), np.floor(t_dis / MD), np.inf),
                np.where(np.isfinite(t_e), np.floor(t_e / MD), np.inf),
                np.where(np.isfinite(t_c), np.floor(t_c / MD), np.inf),
                np.floor(t_ad / MD))


def seq_curves(tl, Xpat, M, m_dis, m_ev, m_comp, m_adm):
    """Both arms of the "disengage during month M" trial. Grace = month M itself,
    so a death in it is credited to BOTH arms."""
    m0 = M - 1
    elig = (m_ev >= m0) & (m_comp >= m0) & (m_adm > m0) & (m_dis >= m0)
    here = m_dis == m0
    out = {}
    for arm in ("disengage", "remain"):
        if arm == "disengage":
            r_dev = np.where(here, np.inf, 1.0)
            mode, T_ = "window_end", 1
        else:
            r_dev = m_dis - m0
            mode, T_ = "all", None
        rE, rC, rA = m_ev - m0, m_comp - m0, m_adm - m0
        r_cens = np.minimum(np.minimum(r_dev, rC), np.minimum(rA, float(H)))
        has = elig & (rE < r_cens) & (rE < float(H)) & np.isfinite(rE)
        nm = np.clip(np.where(has, rE + 1, r_cens), 0, H)
        nm = np.where(elig, nm, 0).astype(np.int32)
        dev = (elig & ~has & np.isfinite(r_dev) & (r_dev <= rA) & (r_dev <= rC)
               & (r_dev < float(H)))
        k = nm > 0
        if not k.any():
            return None, None, 0
        d = {"nmonths": nm[k], "event_month": np.where(has, rE, np.inf)[k],
             "dev_binds": dev[k], "pat": tl["pat"].to_numpy()[k],
             "at_risk_mode": mode, "T": T_}
        out[arm] = ccw.weighted_risk(ccw.add_ipcw(ccw.expand(d), Xpat, verbose=False))
    return out["disengage"], out["remain"], int((elig & here).sum())


# =============================================================================
# A. What "month M" indexes
# =============================================================================
print("=" * 78)
print("A. WHAT 'MONTH M' INDEXES")
print("=" * 78)
print("  t_dis = (closure date - 30 d) - treatment start,  floored at 0")
print("  m_dis = floor(t_dis / 30.4)   <- indexes the trial")
print("  => the inferred DISENGAGEMENT month, NOT the declaration month.\n")

raw = pd.read_csv(imp, low_memory=False,
                  usecols=["sinan_clean", "itt_group", "best_start", "end_date"])
for c in ("best_start", "end_date"):
    raw[c] = pd.to_datetime(raw[c], errors="coerce")
raw["closure_d"] = (raw.end_date - raw.best_start).dt.days
ltfu = raw[raw.itt_group.eq("Loss to follow-up")].copy()
ltfu["t_dis"] = np.maximum(ltfu.closure_d - 30, 0)
ltfu["m_dis"] = np.floor(ltfu.t_dis / MD).astype(int)
ltfu["declar_d"] = ltfu.closure_d          # closure IS the declaration
ltfu["m_declar"] = np.floor(ltfu.declar_d / MD).astype(int)

m1 = ltfu[ltfu.m_dis == 0]
print(f"  month 1 (m_dis = 0): n = {len(m1):,}")
print(f"    inferred disengagement day: median {m1.t_dis.median():.0f} "
      f"(IQR {m1.t_dis.quantile(.25):.0f}-{m1.t_dis.quantile(.75):.0f})")
print(f"    DECLARATION day:            median {m1.declar_d.median():.0f} "
      f"(IQR {m1.declar_d.quantile(.25):.0f}-{m1.declar_d.quantile(.75):.0f})")
for k, v in m1.m_declar.value_counts().sort_index().items():
    print(f"      declared in month {k + 1}: {v:>6,} ({100 * v / len(m1):5.1f}%)")
bad_n = int((m1.closure_d < 30).sum())
print(f"    closure < 30 d, so the definition cannot have been met: "
      f"{bad_n:,} ({100 * bad_n / len(m1):.1f}% of month 1) "
      f"-- t_dis floored up from a negative value")
print("\n  for contrast:")
for M in (2, 5):
    s = ltfu[ltfu.m_dis == M - 1]
    print(f"    month {M}: n={len(s):>6,}  disengagement median day "
          f"{s.t_dis.median():>3.0f}   declaration median day "
          f"{s.declar_d.median():>3.0f}   floored-negative "
          f"{int((s.closure_d - 30 < 0).sum())}")
print("\n  => label any figure axis 'month of disengagement (inferred as 30 days")
print("     before case closure)', NOT 'month of LTFU'.")

# =============================================================================
# B. The identifiability limit
# =============================================================================
print("\n" + "=" * 78)
print("B. DEATHS BEFORE LTFU COULD BE DECLARED")
print("=" * 78)
anyd = tl_full["__any_death__"].to_numpy().astype(bool)
td = tl_full["t_death"].to_numpy()
t_admin = tl_full["t_admin"].to_numpy()
is_l = tl_full["is_ltfu"].to_numpy().astype(bool)
n = len(tl_full)
print(f"  cohort {n:,}   recorded LTFU {is_l.sum():,} ({100 * is_l.mean():.1f}%)")
for cut in (30, 60, 90, 180):
    d = anyd & (td <= cut)
    print(f"    died by day {cut:>3}: {d.sum():>6,} ({100 * d.sum() / n:.2f}% of "
          f"cohort)   of whom recorded LTFU: {int((d & is_l).sum()):>5,}")
d30 = int((anyd & (td <= 30)).sum())
print(f"\n  Those {d30:,} deaths are excluded from BOTH arms in any")
print("  declaration-aligned design. An unknown fraction had already stopped")
print("  attending and would have been declared LTFU had they survived; TBweb")
print("  closes them as Obito either way. The effect of disengagement BEHAVIOUR")
print("  including them is therefore not identifiable from these data.")
print("  Direction is assertable: they are the frailest of the would-be exposed,")
print("  so excluding them biases toward the null.")

# =============================================================================
# C. Where the month-1 deficit arises
# =============================================================================
print("\n" + "=" * 78)
print("C. WHERE THE MONTH-1 DEFICIT ARISES")
print("=" * 78)
tl_p, Xpat = ccw.attach_patterns(tl_full, ccw.COVS)
md, me, mc, ma = month_arrays(tl_p)
GRID = [g for g in (3, 6, 12, 24, 36, 60) if g <= H]
for M in (1, 5):
    r1c, r0c, ne = seq_curves(tl_p, Xpat, M, md, me, mc, ma)
    if r1c is None:
        print(f"  month {M}: not estimable")
        continue
    print(f"\n  sequential trial, month {M} "
          f"(relative months from the trial origin, n exposed {ne:,})")
    print(f"  {'rel month':>10} {'disengage %':>12} {'remain %':>10} "
          f"{'RD pp':>8} {'RR':>6}")
    for g in GRID:
        r1 = ccw.risk_at(r1c, g)
        r0 = ccw.risk_at(r0c, g)
        print(f"  {g:>10} {100 * r1:>12.2f} {100 * r0:>10.2f} "
              f"{100 * (r1 - r0):>+8.2f} {r1 / r0 if r0 > 0 else np.nan:>6.2f}")
print("\n  => the deficit is established in the first few months and then erodes.")
print("     The same shape appears at month 5, far smaller, which is why it")
print("     crosses over there and never does at month 1. One mechanism: the")
print("     comparator is in care and still absorbing on-treatment TB deaths.")

# =============================================================================
# D. Does excluding the definition-inconsistent records help?
# =============================================================================
print("\n" + "=" * 78)
print("D. EXCLUDING THE RECORDS WITH CLOSURE < 30 DAYS")
print("=" * 78)
closure_d = (raw.end_date - raw.best_start).dt.days.to_numpy(dtype="float64")
is_ltfu_raw = raw.itt_group.eq("Loss to follow-up").to_numpy()
bad30 = is_ltfu_raw & np.isfinite(closure_d) & (closure_d < 30)
assert len(bad30) == len(tl_full), "raw extract and timeline are misaligned"
print(f"  definition-inconsistent LTFU: {bad30.sum():,} of "
      f"{is_ltfu_raw.sum():,} ({100 * bad30.sum() / is_ltfu_raw.sum():.1f}%)")

died60 = anyd & (td <= H * MD)
reach60 = t_admin >= H * MD
m_dis_f = np.where(np.isfinite(tl_full["t_dis"].to_numpy()),
                   np.floor(tl_full["t_dis"].to_numpy() / MD), np.inf)
in_m1 = m_dis_f == 0
print(f"\n  who are they? crude {H}-month mortality among month-1 LTFU:")
for lbl, msk in (("closure < 30 d", in_m1 & bad30), ("closure >= 30 d", in_m1 & ~bad30)):
    sub = msk & reach60
    print(f"    {lbl:>15}: n={msk.sum():>6,}  reaching {H} mo n={sub.sum():>6,}  "
          f"crude mortality {100 * died60[sub].mean():>5.2f}%")


def contrasts(tl):
    tl = tl.reset_index(drop=True)
    tl, Xp = ccw.attach_patterns(tl, ccw.COVS)
    md_, me_, mc_, ma_ = month_arrays(tl)
    res = {}
    o = {}
    for arm in ("disengage", "remain"):
        a = ccw.build_arm(tl, arm, T=6, cause="all_cause")
        o[arm] = ccw.weighted_risk(ccw.add_ipcw(ccw.expand(a), Xp, verbose=False))
    r1 = ccw.risk_at(o["disengage"], H)
    r0 = ccw.risk_at(o["remain"], H)
    res["nested_T6"] = (r1 / r0, 100 * (r1 - r0), np.nan)
    for M in range(1, 7):
        r1c, r0c, ne = seq_curves(tl, Xp, M, md_, me_, mc_, ma_)
        if r1c is None:
            res[f"seq_M{M}"] = (np.nan, np.nan, 0)
            continue
        r1 = ccw.risk_at(r1c, H)
        r0 = ccw.risk_at(r0c, H)
        res[f"seq_M{M}"] = (r1 / r0 if r0 > 0 else np.nan, 100 * (r1 - r0), ne)
    return res


full = contrasts(tl_full)
excl = contrasts(tl_full.loc[~bad30])
print("\n                   full cohort              excluding closure<30d")
hdr = (f"{'contrast':>11} | {'RR':>6} {'RD pp':>8} {'n_exp':>7} | "
       f"{'RR':>6} {'RD pp':>8} {'n_exp':>7} | {'dRD':>7}")
print(hdr)
print("-" * len(hdr))
rows = []
for k in full:
    a, b = full[k], excl[k]
    na = f"{a[2]:,}" if np.isfinite(a[2]) and a[2] else "-"
    nb = f"{b[2]:,}" if np.isfinite(b[2]) and b[2] else "-"
    print(f"{k:>11} | {a[0]:>6.2f} {a[1]:>+8.2f} {na:>7} | "
          f"{b[0]:>6.2f} {b[1]:>+8.2f} {nb:>7} | {b[1] - a[1]:>+7.2f}")
    rows.append({"contrast": k, "rr_full": a[0], "rd_full": a[1], "n_full": a[2],
                 "rr_excl": b[0], "rd_excl": b[1], "n_excl": b[2],
                 "rd_delta": b[1] - a[1]})
print("\n  => month 1 gets MORE negative, because those patients have higher")
print("     crude mortality than the rest of month 1. Months 2-6 do not move:")
print("     with m_dis = 0 they were already ineligible for every later trial.")
print("     Excluding them is not a fix. Rejected.")

out = ccw.OUTDIR / f"ccw_month1_diagnostics_h{H}.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\nwrote {out.relative_to(ccw.ROOT)}")
