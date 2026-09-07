#!/usr/bin/env python3
"""Cause-specific CCW at a 60-month horizon, with a patient-level bootstrap.

WHY THIS EXISTS
The CCW is the primary design, but until now the only cause-specific CCW output
was ccw_v3_shift30_main.csv, which stops at 24 months and carries no intervals.
The five-year cause split therefore had to be read off the rolling landmark --
a design whose justification is timing, not cause. Owner decision 2026-08-26:
the cause decomposition belongs on the primary design, so it is computed here.

Structure deliberately mirrors ccw_seq_bootstrap.py (same monkeypatched horizon,
same MI-pooled b=0, same patient-level resample) so the cause rows and the
all-cause row are produced by one pipeline and are directly comparable.

ADDITIVITY DOES NOT HOLD HERE, AND THAT IS NOT A BUG.
build_arm() censors a competing-cause death at its own event time, and unknown-
cause deaths are censored in BOTH the TB and the non-TB analysis (ccw_v3.py:315).
So tb_hybrid + nontb_hybrid need not equal all_cause, unlike the Aalen-Johansen
decomposition in the landmark (script 46d), where the parts sum to the whole by
construction. Any manuscript sentence claiming the components sum to the
all-cause excess, or reporting TB as a PERCENTAGE of it, is a statement about
the landmark decomposition and does NOT transfer to these numbers.

The paired TB-minus-non-TB difference is taken WITHIN replicate. Both causes are
computed on the same resample for exactly this reason: overlapping marginal
intervals are not a test of whether the excess is concentrated in TB deaths.

  B=3   HORIZON=60 python3 CCW_analysis/ccw_cause_bootstrap.py   # smoke test
  B=300 HORIZON=60 python3 CCW_analysis/ccw_cause_bootstrap.py   # real run
"""
import importlib.util
import os
import time
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
T_NESTED = 6
ccw.HORIZON_M = H

# all_cause is recomputed here rather than reused from ccw_bootstrap_h60.csv so
# that every row in the output comes from one pipeline and one set of replicates.
RUN = ["all_cause", "tb_hybrid", "nontb_hybrid"]

rng = np.random.default_rng(2026)

# A low-B smoke run must not overwrite a real run's outputs; same guard as
# ccw_seq_bootstrap.py, which learned it the hard way.
SUF = "" if B >= 50 else "_smoke"


def nested_rd(tl, Xpat, cause, T=T_NESTED):
    """Risk ratio and risk difference (pp) at H months for one cause."""
    out = {}
    for a in ("disengage", "remain"):
        arm = ccw.build_arm(tl, a, T=T, cause=cause)
        rows = ccw.add_ipcw(ccw.expand(arm), Xpat, verbose=False)
        out[a] = ccw.weighted_risk(rows)
    r1 = ccw.risk_at(out["disengage"], H)
    r0 = ccw.risk_at(out["remain"], H)
    return (r1 / r0 if r0 > 0 else np.nan), 100 * (r1 - r0)


lookup = ccw.build_cause_lookup(verbose=False)
imps = sorted(ccw.MI_DIR.glob("imp_*.csv"))
base = ccw.load_timeline(imps[0], lookup, verbose=False)
n = len(base)

print(f"cause-specific CCW  horizon={H}mo  T={T_NESTED}  B={B}  "
      f"imputations={len(imps)}  n={n:,}", flush=True)

recs = []
t0 = time.time()

# b=0 is the POINT estimate, MI-pooled across every imputation (Rubin). Using a
# single imputation here is what produced the old 2.22-versus-2.23 discrepancy.
pooled = {}
for imp in imps:
    tl_i = ccw.load_timeline(imp, lookup, verbose=False)
    tl_i, Xpat_i = ccw.attach_patterns(tl_i, ccw.COVS)
    for cause in RUN:
        rr, rd = nested_rd(tl_i, Xpat_i, cause)
        pooled.setdefault(cause, []).append((rr, rd))
for w, vals in pooled.items():
    recs.append({"b": 0, "which": w,
                 "rr": sum(v[0] for v in vals) / len(vals),
                 "rd": sum(v[1] for v in vals) / len(vals)})
print(f"  MI-pooled point estimates done in {time.time()-t0:.0f}s", flush=True)

for b in range(1, B + 1):
    imp = imps[rng.integers(len(imps))]
    tl = ccw.load_timeline(imp, lookup, verbose=False)
    idx = rng.integers(0, n, n)              # resample PATIENTS with replacement
    tl = tl.iloc[idx].reset_index(drop=True)
    tl, Xpat = ccw.attach_patterns(tl, ccw.COVS)
    for cause in RUN:                        # same resample for every cause
        rr, rd = nested_rd(tl, Xpat, cause)
        recs.append({"b": b, "which": cause, "rr": rr, "rd": rd})
    if b % 10 == 0 or b <= 2:
        el = time.time() - t0
        print(f"  replicate {b}/{B}  elapsed {el/60:.1f} min  "
              f"projected total {el/b*B/60:.0f} min", flush=True)

d = pd.DataFrame(recs)
d.to_csv(ccw.OUTDIR / f"ccw_cause_bootstrap_draws_h{H}{SUF}.csv", index=False)

pt = d[d.b == 0].set_index("which")
bs = d[d.b > 0]
rows = []
for w, g in bs.groupby("which"):
    rows.append({
        "cause": w, "label": ccw.CAUSES[w][1], "horizon_m": H, "T": T_NESTED,
        "rr": pt.loc[w, "rr"], "rr_lo": g.rr.quantile(.025), "rr_hi": g.rr.quantile(.975),
        "rd": pt.loc[w, "rd"], "rd_lo": g.rd.quantile(.025), "rd_hi": g.rd.quantile(.975),
        "B_ok": int(g.rd.notna().sum())})
out = pd.DataFrame(rows)
out["o"] = out.cause.map({k: i for i, k in enumerate(RUN)})
out = out.sort_values("o").drop(columns="o")
p = ccw.OUTDIR / f"ccw_cause_bootstrap_h{H}{SUF}.csv"
out.to_csv(p, index=False)

print(f"\n=== cause-specific CCW at {H} months, T={T_NESTED}, B={B} ===")
print(f"{'cause':>14} {'RR':>6} {'95% CI':>16} {'RD pp':>8} {'95% CI':>20} {'B':>5}")
print("-" * 78)
for _, r in out.iterrows():
    print(f"{r.cause:>14} {r.rr:>6.2f} {f'({r.rr_lo:.2f}-{r.rr_hi:.2f})':>16} "
          f"{r.rd:>+8.2f} {f'({r.rd_lo:+.2f} to {r.rd_hi:+.2f})':>20} {r.B_ok:>5}")

# Sum of the parts against the whole. Reported so the gap is visible rather than
# assumed away: unknown-cause deaths are censored in both cause analyses, so a
# shortfall here is expected and is the reason no percentage-of-total is emitted.
s = pt.loc["tb_hybrid", "rd"] + pt.loc["nontb_hybrid", "rd"]
print(f"\nadditivity: tb {pt.loc['tb_hybrid','rd']:+.3f} + nontb "
      f"{pt.loc['nontb_hybrid','rd']:+.3f} = {s:+.3f} pp   vs all-cause "
      f"{pt.loc['all_cause','rd']:+.3f} pp   (gap {s - pt.loc['all_cause','rd']:+.3f})")
print("  a nonzero gap is EXPECTED: unclassified deaths are censored in both "
      "cause analyses.\n  Do not report a TB share of the all-cause excess from "
      "these numbers.")

# Paired within-replicate TB minus non-TB. Marginal intervals overlapping is not
# a test of concentration; the two causes share patients within a replicate.
w = bs.pivot(index="b", columns="which", values="rd").dropna(
    subset=["tb_hybrid", "nontb_hybrid"])
diff = w["tb_hybrid"] - w["nontb_hybrid"]
pt_diff = pt.loc["tb_hybrid", "rd"] - pt.loc["nontb_hybrid", "rd"]
pd.DataFrame([{"contrast": "tb_minus_nontb_rd", "horizon_m": H,
               "estimate": pt_diff,
               "lo": diff.quantile(.025), "hi": diff.quantile(.975),
               "p_gt0": float((diff > 0).mean()), "B_ok": int(len(diff))}]
             ).to_csv(ccw.OUTDIR / f"ccw_cause_paired_h{H}{SUF}.csv", index=False)
print(f"\npaired TB - nonTB risk difference: {pt_diff:+.2f} pp "
      f"({diff.quantile(.025):+.2f} to {diff.quantile(.975):+.2f}); "
      f"share of replicates > 0: {(diff > 0).mean():.3f}")
print(f"\nwrote {p.relative_to(ccw.ROOT)}   total {(time.time()-t0)/60:.1f} min")
