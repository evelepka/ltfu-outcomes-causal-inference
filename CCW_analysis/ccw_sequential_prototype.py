#!/usr/bin/env python3
"""PROTOTYPE: sequential trial emulation with a grace period, to test whether a
cloning design can express a TIMING gradient of disengagement.

WHY
---
`ccw_v3.py --disjoint` uses the right strategy shape ("disengage DURING month T")
but both arms still start at day 0 with no eligibility restriction, so for T=6 the
arms accrue five months of IDENTICAL follow-up before the strategies diverge. That
dilution is what flattens the gradient. It is an artifact of carrying a nested-CCW
origin into a disjoint design, not evidence that disjoint designs cannot work.

This prototype instead emulates a SEPARATE trial per month:

  eligibility  alive AND still in care at the START of month M
  origin       that month boundary (NOT day 0)
  strategies   "disengage during month M" vs "remain engaged"
  grace        month M itself; a death during it is credited to BOTH clones,
               because at that point both strategies remain compatible with the
               observed data
  censoring    a clone is censored when it departs its assigned strategy,
               with stabilized IPCW for that artificial censoring
  follow-up    up to HORIZON_M months FROM THE ORIGIN

Because both arms are in care at the moment of assignment, this also removes the
asymmetry that drives the age >=65 inversion in the landmark (exposed defined by
having left, comparator by having stayed).

WHAT WE LEARNED (added 2026-08-19, after the design decision)
------------------------------------------------------------
The design works, and the timing question it was built to answer came back
qualified. Do not read the original framing above as still open.

  * IMPLEMENTATION VALIDATED. Nested T=1 and sequential M=1 are the same
    intervention on the same population and agree to within rounding.
  * THE GRADIENT IS NOT STATISTICALLY RESOLVABLE. An omnibus test that the risk
    difference is equal across disengagement months does not reject. See
    ccw_timing_heterogeneity.py. The apparent peak is within noise, and the
    landmark gives the same answer once the parameterisation is matched (42c).
    Do not report a peak month from this.
  * MONTH 1 IS NOT USABLE and is excluded. Its origin is day 0, up to 30 days
    before the exposure is realised; moving ONLY the origin flips its sign. See
    ccw_origin_decomposition.py and ccw_month1_diagnostics.py.
  * THE >=65 CLAIM BELOW HELD. CCW subgroups are monotone and positive there,
    which is why the 2026-08-19 decision puts subgroups on this design.
  * STATUS: superseded for TIMING by the rolling landmark (45b, 42c). Retained
    as the design demonstration and as the source of the per-month point
    estimates. Still 1 imputation by default, no CIs (see ccw_seq_bootstrap.py),
    and not wired into check.sh.

WHAT TO LOOK AT BEFORE TRUSTING IT
----------------------------------
Positivity is the honest risk: only ~2% of patients disengage in any given month,
so "everyone disengages during month M" is a weakly supported intervention. Judge
this prototype on (i) mean/max stabilized weights, (ii) how wide the per-trial
interval is, (iii) whether the point estimates order sensibly against the rolling
landmark's timing curve. Do not read the point estimates in isolation.

Estimator and weight machinery are imported from ccw_v3.py so this is not a
parallel implementation.

Usage:  python3 CCW_analysis/ccw_sequential_prototype.py            # months 2 and 5
        python3 CCW_analysis/ccw_sequential_prototype.py 2 3 4 5    # a gradient
        N_IMP=5 python3 CCW_analysis/ccw_sequential_prototype.py
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ccw", HERE / "ccw_v3.py")
ccw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccw)

# months of follow-up from each trial origin; HORIZON env var overrides.
# At long horizons ccw's internal HORIZON_M must match, since the IPCW cell
# encoding and the weighted-risk array length both key off it.
HORIZON_REL = int(os.environ.get("HORIZON", ccw.HORIZON_M))
ccw.HORIZON_M = HORIZON_REL
MD = ccw.MONTH_DAYS


def build_seq_arm(tl, arm, m0, cause="all_cause"):
    """One arm of the trial whose grace window is absolute month m0 (0-indexed).

    Returns the same dict shape ccw.expand() consumes, but on a RELATIVE month
    clock: r = absolute month - m0, so r=0 is the grace month.
    """
    t_dis = tl["t_dis"].to_numpy()
    t_death = tl["t_death"].to_numpy()
    t_admin = tl["t_admin"].to_numpy()
    ev_col = ccw.CAUSES[cause][0]
    is_target = tl[ev_col].to_numpy().astype(bool)
    any_death = tl["__any_death__"].to_numpy().astype(bool)

    t_compete = np.where(any_death & ~is_target, t_death, np.inf)
    t_event = np.where(is_target, t_death, np.inf)

    with np.errstate(invalid="ignore"):
        m_dis = np.where(np.isfinite(t_dis), np.floor(t_dis / MD), np.inf)
        m_event = np.where(np.isfinite(t_event), np.floor(t_event / MD), np.inf)
        m_compete = np.where(np.isfinite(t_compete), np.floor(t_compete / MD), np.inf)
        m_admin = np.floor(t_admin / MD)

    # ---- ELIGIBILITY: alive, uncensored and still in care at the START of m0 --
    elig = (m_event >= m0) & (m_compete >= m0) & (m_admin > m0) & (m_dis >= m0)

    # relative clock
    r_dis = m_dis - m0
    r_event = m_event - m0
    r_compete = m_compete - m0
    r_admin = m_admin - m0

    if arm == "disengage":
        # adherent = disengaged during the grace month (r_dis == 0). Everyone else
        # departs the strategy at the END of the grace month, i.e. r = 1.
        adherent = r_dis == 0
        r_dev = np.where(adherent, np.inf, 1.0)
        at_risk_mode, T = "window_end", 1      # expand() -> at_risk = (month == 0)
    elif arm == "remain":
        # departs whenever they disengage
        r_dev = r_dis
        at_risk_mode, T = "all", None
    else:
        raise ValueError(arm)

    r_cens = np.minimum(np.minimum(r_dev, r_compete),
                        np.minimum(r_admin, float(HORIZON_REL)))
    has_event = (elig & (r_event < r_cens) & (r_event < float(HORIZON_REL))
                 & np.isfinite(r_event))

    nmonths = np.where(has_event, r_event + 1, r_cens)
    nmonths = np.clip(nmonths, 0, HORIZON_REL)
    nmonths = np.where(elig, nmonths, 0).astype(np.int32)

    dev_binds = (elig & (~has_event) & np.isfinite(r_dev) & (r_dev <= r_admin)
                 & (r_dev <= r_compete) & (r_dev < float(HORIZON_REL)))

    keep = nmonths > 0
    return {"nmonths": nmonths[keep],
            "event_month": np.where(has_event, r_event, np.inf)[keep],
            "dev_binds": dev_binds[keep],
            "pat": tl["pat"].to_numpy()[keep],
            "at_risk_mode": at_risk_mode,
            "T": T,
            "n_elig": int(keep.sum()),
            "n_adherent": int((keep & (r_dis == 0)).sum()) if arm == "disengage" else None}


def one_trial(tl, Xpat, M, verbose=True):
    """M is 1-indexed as in the paper: M=2 is the second month of treatment."""
    m0 = M - 1
    out = {}
    for arm in ("disengage", "remain"):
        a = build_seq_arm(tl, arm, m0)
        if a["n_elig"] == 0:
            return None
        rows = ccw.expand(a)
        rows = ccw.add_ipcw(rows, Xpat, verbose=verbose, label=f"M{M} {arm}")
        ci = ccw.weighted_risk(rows, label=f"M{M} {arm}")
        out[arm] = {"ci": ci, "n": a["n_elig"], "adh": a["n_adherent"],
                    "w_mean": float(rows["sw"].mean()),
                    "w_max": float(rows["sw"].max()),
                    "w_trunc": float(np.mean(rows["swt"] != rows["sw"]))}
    return out


def main(months):
    imp_dir = ccw.MI_DIR
    imps = sorted(imp_dir.glob("imp_*.csv"))
    n_imp = int(os.environ.get("N_IMP", 1))
    imps = imps[:n_imp]
    print(f"[seq] sequential-trial prototype | {len(imps)} imputation(s) | "
          f"trials at month(s) {months} | horizon {HORIZON_REL} mo from each origin")

    lookup = ccw.build_cause_lookup(verbose=False)
    rows = []
    for p in imps:
        tl = ccw.load_timeline(p, lookup, verbose=False)
        tl, Xpat = ccw.attach_patterns(tl, ccw.COVS)
        for M in months:
            r = one_trial(tl, Xpat, M, verbose=(p is imps[0]))
            if r is None:
                print(f"  M{M}: not estimable")
                continue
            for h in sorted({6, 12, 24, HORIZON_REL}):
                r1 = ccw.risk_at(r["disengage"]["ci"], h)
                r0 = ccw.risk_at(r["remain"]["ci"], h)
                rows.append({"imp": p.name, "M": M, "horizon_rel_mo": h,
                             "risk_dis": 100 * r1, "risk_rem": 100 * r0,
                             "rr": r1 / r0 if r0 > 0 else np.nan,
                             "rd": 100 * (r1 - r0),
                             "n_elig": r["disengage"]["n"],
                             "n_adherent": r["disengage"]["adh"],
                             "w_mean_dis": r["disengage"]["w_mean"],
                             "w_max_dis": r["disengage"]["w_max"],
                             "w_trunc_dis": r["disengage"]["w_trunc"]})

    import pandas as pd
    d = pd.DataFrame(rows)
    if d.empty:
        print("nothing estimable")
        return 1
    g = (d.groupby(["M", "horizon_rel_mo"])
           .agg(risk_dis=("risk_dis", "mean"), risk_rem=("risk_rem", "mean"),
                rr=("rr", "mean"), rd=("rd", "mean"),
                n_elig=("n_elig", "max"), n_adherent=("n_adherent", "max"),
                w_mean=("w_mean_dis", "mean"), w_max=("w_max_dis", "max"),
                w_trunc=("w_trunc_dis", "mean"))
           .reset_index())
    print("\n=== SEQUENTIAL TRIALS: disengage during month M vs remain engaged ===")
    print("    (clock is months FROM THE TRIAL ORIGIN, not from treatment start)\n")
    hdr = (f"{'M':>2} {'horiz':>6} {'eligible':>9} {'adherent':>9} "
           f"{'risk_dis':>9} {'risk_rem':>9} {'RR':>6} {'RD pp':>7} "
           f"{'w_mean':>7} {'w_max':>8} {'trunc':>6}")
    print(hdr)
    print("-" * len(hdr))
    for _, r in g.iterrows():
        print(f"{int(r.M):>2} {int(r.horizon_rel_mo):>6} {int(r.n_elig):>9,} "
              f"{int(r.n_adherent):>9,} {r.risk_dis:>8.2f}% {r.risk_rem:>8.2f}% "
              f"{r.rr:>6.2f} {r.rd:>+7.2f} {r.w_mean:>7.3f} {r.w_max:>8.1f} "
              f"{r.w_trunc:>5.1%}")

    out = ccw.OUTDIR / "ccw_sequential_prototype.csv"
    g.to_csv(out, index=False)
    print(f"\n[seq] wrote {out.relative_to(ccw.ROOT)}")
    print("\nPROTOTYPE. Judge on weights and interval width, not point estimates.")
    return 0


if __name__ == "__main__":
    ms = [int(x) for x in sys.argv[1:]] or [2, 5]
    sys.exit(main(ms))
