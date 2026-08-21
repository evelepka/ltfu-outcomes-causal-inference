#!/usr/bin/env python3
"""Does the CCW timing gradient exist at all? Omnibus Wald test of H0: the risk
difference is equal across disengagement months 2-6, using the bootstrap
covariance of the five estimates. Reads the saved draws, so nothing is refitted.

This is the CCW half of the matched-parameterisation comparison in handoff
section 7. The landmark half is 42c. Both must use the SAME parameterisation --
free monthly bins -- or the comparison measures the smoothing basis rather than
the design.

A max-minus-min "spread" statistic is NOT a valid test here and must not be
substituted: max - min is non-negative in every replicate, so its bootstrap
interval can never contain zero.

  python3 CCW_analysis/ccw_timing_heterogeneity.py      # after ccw_seq_bootstrap.py
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(os.environ.get("TB_ABANDONMENT_ROOT",
                           Path(__file__).resolve().parent.parent))
P = BASE / "CCW_analysis/results_v3/ccw_bootstrap_draws_h60.csv"
d = pd.read_csv(P)
MS = [f"seq_M{m}" for m in (2, 3, 4, 5, 6)]
pt = d[d.b == 0].set_index("which").loc[MS, "rd"].to_numpy()
W = d[d.b > 0].pivot(index="b", columns="which", values="rd")[MS].dropna()
S = np.cov(W.to_numpy(), rowvar=False)

print(f"point RDs (months 2-6): {np.round(pt, 2)}")
print(f"point spread max-min:   {pt.max()-pt.min():+.2f} pp   "
      f"(descriptive only -- see the docstring: this is not a test)")
print(f"bootstrap SEs:          {np.round(np.sqrt(np.diag(S)), 2)}\n")

# successive-difference contrasts: 4 x 5
C = np.zeros((4, 5))
for i in range(4):
    C[i, i], C[i, i + 1] = -1.0, 1.0
c = C @ pt
V = C @ S @ C.T
stat = float(c @ np.linalg.solve(V, c))
p_chi = stats.chi2.sf(stat, df=4)
print("H0: RD equal across months 2-6")
print(f"  Wald chi2(4) = {stat:.2f}   p = {p_chi:.4f}")

# fully bootstrap version: recentre draws under H0, how often does the resampled
# statistic exceed the observed one?
Wc = W.to_numpy() - W.to_numpy().mean(axis=0)
null = []
for row in Wc:
    cc = C @ row
    try:
        null.append(float(cc @ np.linalg.solve(V, cc)))
    except np.linalg.LinAlgError:
        pass
null = np.array(null)
print(f"  bootstrap-calibrated p = {max((null >= stat).mean(), 1/len(null)):.4f} "
      f"(reference dist from {len(null)} recentred draws)")

print("\npairwise correlation of the monthly RDs across replicates:")
print(pd.DataFrame(S, index=MS, columns=MS)
      .pipe(lambda x: x.div(np.sqrt(np.diag(S)), axis=0).div(np.sqrt(np.diag(S)), axis=1))
      .round(2).to_string())
