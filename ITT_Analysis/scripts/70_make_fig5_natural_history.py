#!/usr/bin/env python3
"""Figure 5 — the natural history of loss to follow-up.

WHY THIS FIGURE
---------------
It replaces the cause-specific decomposition (demoted to a sentence) and the old
Figure 3 return-pathway figure (cut). Two reviewers could not follow that figure,
and on inspection its panels used three different time origins, three different
populations, and a control arm censored at the moment of return -- which removes
patients precisely when they become high-risk and is most of why its hazard ratio
was 6.93.

EVERY PANEL HERE USES ONE CLOCK: months since loss to follow-up. That single
constraint is the fix.

THE MESSAGE
-----------
Loss to follow-up is usually transient; the patients most likely to come back are
the ones at highest risk; and the deaths are concentrated among those who come
back. This is why the average effect is modest without the finding being
unimportant, and it is a characterisation of LTFU no prior paper offers.

WHAT IS AND IS NOT CLAIMED
--------------------------
Panel D is DESCRIPTIVE. It says where the deaths are, not what returning does.
The effect of return is not identifiable here: the appendix (section 9.3) shows a
1.57-fold residual hazard after baseline adjustment, and
`return_state_vs_baseline_table.csv` shows recorded clinical characteristics at
re-presentation are essentially unchanged from baseline (hospitalised 41.9% ->
41.9%, HIV 15.6% -> 15.2%). The deterioration that drives return is real but
invisible in TBweb, so no adjustment recovers it. Do not add a hazard ratio to
panel D.

Re-engagement is estimated with Aalen-Johansen, treating death before
re-notification as a competing event -- 1-KM would overstate it.

Usage:  python3 ITT_Analysis/scripts/70_make_fig5_natural_history.py
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
OUTDIR = ROOT / "Plos Medicine" / "R1" / "Figures" / "Draft"
RESULTS = ROOT / "ITT_Analysis" / "results"
STEM = "LTFU_natural_history"

BLUE, VERM, GREY, INK = "#0072B2", "#D55E00", "#8A8A8A", "#222222"
AMBER = "#E69F00"
FS_T, FS_L, FS_A = 8.5, 7.0, 6.4
HORIZON = 60.0          # months since LTFU
LANDMARK = 6.0          # months, for panel D


def aalen_johansen(t_ev, cause, grid):
    """CIF for cause 1 with cause 2 competing. cause: 1, 2, or 0 for censored.

    Discrete-time Aalen-Johansen: CIF_1(t) = sum_{s<=t} S(s-) * h_1(s), with
    S depleted by BOTH causes. Using 1-KM on cause 1 alone and censoring cause 2
    would assert the dead remain at risk of returning, and overstate the curve.
    """
    order = np.argsort(t_ev)
    t_ev, cause = t_ev[order], cause[order]
    n = len(t_ev)
    times = np.unique(t_ev[cause > 0])
    S, cif, out, k = 1.0, 0.0, [], 0
    at_risk = n
    ti = 0
    for g in grid:
        while ti < len(times) and times[ti] <= g:
            s = times[ti]
            at_risk = int(np.sum(t_ev >= s))
            d1 = int(np.sum((t_ev == s) & (cause == 1)))
            d2 = int(np.sum((t_ev == s) & (cause == 2)))
            if at_risk > 0:
                cif += S * d1 / at_risk
                S *= 1.0 - (d1 + d2) / at_risk
            ti += 1
        out.append(cif)
    return np.array(out)


def km_mortality(t_ev, ev, grid):
    """1-KM cumulative mortality (single event type)."""
    order = np.argsort(t_ev)
    t_ev, ev = t_ev[order], ev[order]
    times = np.unique(t_ev[ev == 1])
    S, out, ti = 1.0, [], 0
    for g in grid:
        while ti < len(times) and times[ti] <= g:
            s = times[ti]
            nr = int(np.sum(t_ev >= s))
            d = int(np.sum((t_ev == s) & (ev == 1)))
            if nr > 0:
                S *= 1.0 - d / nr
            ti += 1
        out.append(100 * (1 - S))
    return np.array(out)


# ---------------------------------------------------------------- data -------
d = pd.read_csv(ROOT / "ITT_Analysis/data/itt_cohort.csv", low_memory=False)
l = d[d.itt_group.eq("Loss to follow-up")].copy()
for c in ("time_rn", "event_rn", "time_d", "event_d"):
    l[c] = pd.to_numeric(l[c], errors="coerce")
l = l[l.time_rn.notna() & l.time_d.notna()]
print(f"[70] LTFU cohort: {len(l):,}")

# competing-risks structure, in MONTHS since LTFU
t_rn = np.where(l.event_rn.eq(1), l.time_rn * 12, np.inf)
t_dd = np.where(l.event_d.eq(1), l.time_d * 12, np.inf)
t_cn = np.maximum(l.time_rn, l.time_d).to_numpy() * 12      # administrative
T = np.minimum(np.minimum(t_rn, t_dd), t_cn)
CAUSE = np.where((t_rn <= t_dd) & np.isfinite(t_rn), 1,
        np.where(np.isfinite(t_dd), 2, 0))
l["T"], l["CAUSE"] = T, CAUSE
grid = np.linspace(0, HORIZON, 400)

# month of disengagement, from treatment start to LTFU declaration
l["best_start"] = pd.to_datetime(l.best_start, errors="coerce")
l["end_date"] = pd.to_datetime(l.end_date, errors="coerce")
dis_m = ((l.end_date - l.best_start).dt.days - 30) / 30.4
l["dis_grp"] = pd.cut(dis_m, [-99, 2, 4, 999],
                      labels=["<2 months", "2-4 months", "≥4 months"])

# ---------------------------------------------------------------- figure -----
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
(axA, axB), (axC, axD) = axes


def style(ax, ylab, ymax, xlab=False):
    ax.set_xlim(0, HORIZON)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 12, 24, 36, 48, 60])
    ax.tick_params(labelsize=FS_L, length=2.5, pad=1.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GREY)
    ax.set_ylabel(ylab, fontsize=FS_L)
    if xlab:
        ax.set_xlabel("Months since loss to follow-up", fontsize=FS_L)


def title(ax, letter, txt):
    ax.text(0.0, 1.05, letter, transform=ax.transAxes, fontsize=FS_T + 1,
            weight="bold", va="bottom")
    ax.text(0.055, 1.05, txt, transform=ax.transAxes, fontsize=FS_T,
            weight="bold", va="bottom")


rows = []

# ---- A: re-engagement overall ----------------------------------------------
cif = 100 * aalen_johansen(l["T"].to_numpy(), l["CAUSE"].to_numpy(), grid)
style(axA, "Re-engaged in care (%)", 60)
title(axA, "A", "Nearly half re-engage in care")
axA.plot(grid, cif, color=BLUE, lw=1.8)
for mo in (12, 60):
    v = np.interp(mo, grid, cif)
    axA.plot([mo, mo], [0, v], ls=":", lw=0.7, color=GREY)
    axA.text(mo - 1.2, v + 1.6, f"{v:.0f}%", fontsize=FS_A, color=BLUE,
             weight="bold", ha="right")
    rows.append(dict(panel="A", stratum="overall", month=mo, pct=v))
axA.text(HORIZON - 1, 4, f"n = {len(l):,} lost to follow-up", fontsize=FS_A,
         color=GREY, ha="right")

# ---- B: by HIV status -------------------------------------------------------
style(axB, "Re-engaged in care (%)", 60)
title(axB, "B", "Highest-risk patients re-engage most")
for lv, col in (("Positive", VERM), ("Negative", BLUE)):
    s = l[l.hiv_aids.astype(str).eq(lv)]
    if not len(s):
        continue
    c = 100 * aalen_johansen(s["T"].to_numpy(), s["CAUSE"].to_numpy(), grid)
    axB.plot(grid, c, color=col, lw=1.8)
    v24 = np.interp(24, grid, c)
    axB.text(HORIZON - 1, np.interp(HORIZON, grid, c),
             f"HIV {lv.lower()} ({v24:.1f}% at 24 mo)", fontsize=FS_A,
             color=col, ha="right", va="bottom")
    rows.append(dict(panel="B", stratum=f"HIV {lv}", month=24, pct=v24))

# ---- C: by month of disengagement ------------------------------------------
style(axC, "Re-engaged in care (%)", 60, xlab=True)
title(axC, "C", "Earlier disengagement, more re-engagement")
for lv, col in (("<2 months", VERM), ("2-4 months", AMBER), ("≥4 months", BLUE)):
    s = l[l.dis_grp.astype(str).eq(lv)]
    if not len(s):
        continue
    c = 100 * aalen_johansen(s["T"].to_numpy(), s["CAUSE"].to_numpy(), grid)
    axC.plot(grid, c, color=col, lw=1.8)
    v24 = np.interp(24, grid, c)
    axC.text(HORIZON - 1, np.interp(HORIZON, grid, c),
             f"LTFU {lv} ({v24:.1f}%)", fontsize=FS_A, color=col,
             ha="right", va="bottom")
    rows.append(dict(panel="C", stratum=f"LTFU {lv}", month=24, pct=v24))

# ---- D: mortality by re-engagement status, 6-month landmark ----------------
alive = (l.event_d.ne(1)) | (l.time_d * 12 > LANDMARK)
lm = l[alive].copy()
ret6 = lm.event_rn.eq(1) & (lm.time_rn * 12 <= LANDMARK)
style(axD, "Cumulative mortality (%)", 11, xlab=True)
title(axD, "D", "Deaths concentrate among those who return")
g2 = np.linspace(LANDMARK, HORIZON, 300)
for msk, col, lab in ((ret6, VERM, "re-engaged by 6 months"),
                      (~ret6, BLUE, "still disengaged at 6 months")):
    s = lm[msk]
    t = s.time_d.to_numpy() * 12
    ev = s.event_d.eq(1).to_numpy().astype(int)
    keep = t > LANDMARK
    y = km_mortality(t[keep], ev[keep], g2)
    axD.plot(g2, y, color=col, lw=1.8)
    axD.text(HORIZON - 1, min(y[-1], 10.2) - 0.55, f"{lab} (n={len(s):,})",
             fontsize=FS_A, color=col, ha="right", va="top")
    rows.append(dict(panel="D", stratum=lab, month=24,
                     pct=float(np.interp(24, g2, y))))
axD.axvline(LANDMARK, ls=":", lw=0.7, color=GREY)
axD.text(LANDMARK + 0.8, 10.3, "landmark: 6 months", fontsize=FS_A, color=INK)
axD.text(HORIZON - 1, 0.5, "descriptive: where the deaths are,\nnot the effect of returning",
         fontsize=FS_A, color=GREY, ha="right", style="italic")

fig.tight_layout(pad=1.1)
OUTDIR.mkdir(parents=True, exist_ok=True)
for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    p = OUTDIR / f"{STEM}.{ext}"
    fig.savefig(p, facecolor="white", **kw)
    print(f"wrote {p.relative_to(ROOT)}")
plt.close(fig)
out = RESULTS / "fig5_natural_history_values.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"wrote {out.relative_to(ROOT)}")
print(pd.DataFrame(rows).round(1).to_string(index=False))
