#!/usr/bin/env python3
"""Appendix 1.1 on the PRIMARY estimand: the 60-month risk difference.

46g answered the exclusion question with a hazard ratio, which contradicts how
the paper reports itself -- response 1.17 states the primary estimand is the
absolute risk difference and that hazard ratios are secondary summaries. This
runs the clone-censor-weight contrast on both cohorts instead.

Both sides are unimputed, passed to ccw_v3 as a single pseudo-imputation, so the
estimator is identical and any difference is attributable to the cohort. ccw_v3
casts covariates with .astype(str), so a missing value becomes its own pattern
level rather than dropping the row -- the same handling on both sides.
"""
import sys
from pathlib import Path
import pandas as pd

CCW = Path("/Users/evelynlepkadelima/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com"
           "/My Drive/TB SP 2026/LTFU Paper/CCW_analysis")
SC = Path("/private/tmp/claude-501/-Users-evelynlepkadelima-repos-outcomes-after-tb-abandonment"
          "/2318c5ca-50f6-4c22-869e-8151617cd379/scratchpad")
sys.path.insert(0, str(CCW))
import ccw_subgroup_bootstrap as B          # noqa: E402  (repoints ccw_v3 paths)
C = B.C
HORIZON = 60
B.set_horizon(HORIZON)

lookup = C.build_cause_lookup(verbose=False)
rows = []
for tag, d in [("primary", SC / "ccw_primary"), ("with imputed start", SC / "ccw_variant")]:
    tl = C.load_timeline(d / "imp_01.csv", lookup, verbose=False)
    tl, Xp = C.attach_patterns(tl, C.COVS)
    ref = C.reference(tl, "all_cause", Xp)
    r = C.one_contrast(tl, C.PRIMARY_T, "all_cause", Xp, ref=ref)
    r = dict(r) if not isinstance(r, dict) else r
    r["cohort"] = tag
    r["n_patients"] = int(tl["pid"].nunique()) if "pid" in tl else len(tl)
    rows.append(r)
    keys = [k for k in r if any(s in str(k) for s in ("rd", "rr", "risk"))]
    print("  %-20s %s" % (tag, {k: (round(r[k], 4) if isinstance(r[k], float) else r[k])
                                for k in keys}))

df = pd.DataFrame(rows)
out = CCW / "results_v3" / "ccw_txstart_sensitivity_h60.csv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print("colunas:", list(df.columns))
print("wrote", out)
