#!/usr/bin/env python3
"""Figure S1 - CONSORT-style cohort flow diagram for the TB abandonment paper.

Design goals (per reviewer feedback):
  * comfortable, un-squeezed text (no cramped letters);
  * clear vertical separation between the flow boxes;
  * boxes sized snugly to their content (no large empty margins inside boxes).

2026-10-16: the numbers are no longer hard-coded. They were written out once
and never rechecked, so the figure went on describing a cohort the pipeline had
left behind: age <15 read 6,238 against 6,178, data inconsistencies 111 against
1, pre-treatment deaths 2,063 against 549, end dates outside the window 17,624
against 17,834, and no recorded treatment start 1,394 against 2,778. The chain
still summed to 171,048 because every subtotal had been derived from the stale
counts, so nothing looked wrong. They are now derived from Final_table_cleaned.csv
with the same filter definitions as 01_itt_cohort_selection.py, in the order this
figure presents them, and the script asserts that the chain closes on the cohort.

2026-08-31: the final box read N = 171,069, twenty-one more than the cohort the
paper reports, because the flow omitted the "identified only post-mortem"
exclusion that appendix 1.2 lists. That step is now in the last exclusion box and
the total reads 171,048, matching the manuscript. The output path was also
absolute and pointed at the pre-rename "Abandonment Paper" tree, so the script
could not write its own figure; it now resolves the repository root from its own
location, as scripts 47 and 48 do.
"""
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm

# ---- styling ---------------------------------------------------------------
TITLE_FS = 11.5          # bold header line in each box
BODY_FS  = 10.5          # N line / exclusion items
LINE_H   = 0.345         # vertical space per text line (inches)
PAD_V    = 0.20          # internal top/bottom padding (inches)
PAD_H    = 0.22          # internal left/right padding (inches)

BLUE_FILL, BLUE_EDGE   = "#DCE6F2", "#4472C4"
FINAL_FILL             = "#B8CCE4"
PINK_FILL, PINK_EDGE   = "#FBE6E7", "#C0504D"
ARROW_C                = "#3F3F3F"

# ---- content: derived from the data, not typed ------------------------------
# Same source and same filter definitions as 01_itt_cohort_selection.py, applied
# in the order this figure presents them. Exclusion categories overlap, so the
# order decides which category a person is counted under; the convention here
# matches the Study population paragraph in the manuscript.
_ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
_df = pd.read_csv(_ROOT / "Data" / "Final_table_cleaned.csv", low_memory=False)
for _c in ("end_date", "dod", "tx_start"):
    _df[_c] = pd.to_datetime(_df[_c], errors="coerce")

# comprehensive death date: the `dod` field or an Obito treatment outcome
_a = _df.dropna(subset=["dod"]).groupby("sinan_clean")["dod"].max()
_b = (_df[_df.case_outcome.isin({"Obito TB", "Obito NTB"})].dropna(subset=["end_date"])
        .groupby("sinan_clean")["end_date"].max())
_death = pd.concat([_a, _b], axis=1).max(axis=1).rename("death_date")

n_total = _df.sinan_clean.nunique()
_novo = _df[_df.case_type.str.strip().str.title() == "Novo"]
n_nonfirst = n_total - _novo.sinan_clean.nunique()

_blank = _novo.case_outcome.isna() | (_novo.case_outcome.astype(str).str.strip() == "")
_mud = _novo.case_outcome.eq("Mud Diag")
_valid = _novo[~_blank & ~_mud]
_lost = set(_novo.sinan_clean) - set(_valid.sinan_clean)
_sub = _novo[_novo.sinan_clean.isin(_lost)]
_s_mud = set(_sub[_sub.case_outcome.eq("Mud Diag")].sinan_clean)
_s_blank = set(_sub[_blank.loc[_sub.index]].sinan_clean)
n_muddiag = len(_s_mud - _s_blank)
n_missing = len(_s_blank)                      # absorbs those meeting both criteria
n_firstep = _valid.sinan_clean.nunique()

_first = _valid.sort_values("end_date").groupby("sinan_clean", as_index=False).first()
_TR = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}
n_transfer = int(_first.case_outcome.isin(_TR).sum())
_c = _first[~_first.case_outcome.isin(_TR)].join(_death, on="sinan_clean")

def _drop(frame, keep):
    return int((~keep).sum()), frame[keep]

# The date checks are reported together, as the Study population section does:
# a record whose treatment start is missing and one whose end date precedes its
# start are both unusable for time-on-treatment alignment, for the same reason.
n_age,       _c = _drop(_c, _c.age_tb >= 15)
n_predeath,  _c = _drop(_c, ~(_c.death_date.notna() & (_c.death_date <= _c.tx_start)))
n_consistent = len(_c)
n_window,    _c = _drop(_c, _c.end_date.notna()
                        & (_c.end_date >= pd.Timestamp("2013-01-01"))
                        & (_c.end_date <= pd.Timestamp("2023-12-31")))
n_txdates,   _c = _drop(_c, _c.tx_start.notna() & ~(_c.end_date < _c.tx_start))
n_postmortem, _c = _drop(_c, _c.disease_discovery != "Descob. Apos Obito")
n_cohort = len(_c)

# the chain must close on the cohort the paper reports
assert n_total - n_nonfirst - n_muddiag - n_missing == n_firstep, "box 2 does not close"
assert n_firstep - n_transfer - n_age - n_predeath == n_consistent, "box 3"
assert n_consistent - n_window - n_txdates - n_postmortem == n_cohort, "box 4"
_expected = len(pd.read_csv(_ROOT / "ITT_Analysis" / "data" / "itt_cohort.csv",
                            usecols=["sinan_clean"], low_memory=False))
assert n_cohort == _expected, f"flow ends at {n_cohort:,}, cohort file has {_expected:,}"
print(f"[61] derived chain closes on {n_cohort:,}")

def _n(x): return f"{x:,}"

flow = [
    ("TBweb records, 2013 to 2023", f"N = {_n(n_total)} individuals"),
    ("First TB episode, eligible record", f"N = {_n(n_firstep)}"),
    ("Resident, adult, with consistent data", f"N = {_n(n_consistent)}"),
    ("Primary cohort: treatment within study\nperiod and with recorded start date",
     f"N = {_n(n_cohort)}"),
]

excl = [
    ["Excluded:",
     f"No first-episode (new case) record, n = {_n(n_nonfirst)}",
     f"Change of diagnosis, n = {_n(n_muddiag)}",
     f"Missing or incomplete outcome, n = {_n(n_missing)}"],
    ["Excluded:",
     f"Inter-state transfer during treatment, n = {_n(n_transfer)}",
     f"Age < 15 years, n = {_n(n_age)}",
     f"Pre-treatment deaths, n = {_n(n_predeath)}"],
    ["Excluded:",
     f"End dates outside 2013 to 2023, n = {_n(n_window)}",
     f"Missing or inconsistent treatment dates, n = {_n(n_txdates)}",
     f"Identified only post-mortem, n = {_n(n_postmortem)}"],
]

# ---- geometry (all in inches) ----------------------------------------------
FIG_W, FIG_H = 11.0, 9.0

# left (flow) column
L_X0, L_W = 0.45, 4.55          # left edge, width
L_CX = L_X0 + L_W / 2.0         # vertical-connector x

# right (exclusion) column
R_X0, R_W = 6.05, 4.55

TOP = FIG_H - 0.35              # top of first box
S   = 2.42                      # centre-to-centre spacing of flow boxes


def box_height(nlines, title_extra=0.0):
    return nlines * LINE_H + 2 * PAD_V + title_extra


# flow box line counts (title may wrap to 2 lines + N line)
flow_lines = [t[0].count("\n") + 2 for t in flow]
flow_h = [box_height(n) for n in flow_lines]

# centre y of each flow box
flow_cy = [TOP - flow_h[0] / 2.0]
for i in range(1, 4):
    flow_cy.append(flow_cy[0] - i * S)

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.invert_yaxis()              # so larger y plotted lower is unnecessary; keep normal
ax.set_axis_off()
# undo inversion: we want y increasing upward
ax.set_ylim(0, FIG_H)


def draw_box(cx, cy, w, h, lines, fill, edge, bold_first=True, align="center"):
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=0.10",
        linewidth=1.4, edgecolor=edge, facecolor=fill, mutation_aspect=1.0)
    ax.add_patch(box)
    # stack text lines from top of inner area
    n = len(lines)
    text_top = cy + h / 2 - PAD_V
    for i, (txt, fs, bold) in enumerate(lines):
        y = text_top - (i + 0.5) * LINE_H
        if align == "center":
            x = cx
            ha = "center"
        else:
            x = cx - w / 2 + PAD_H
            ha = "left"
        ax.text(x, y, txt, ha=ha, va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", color="#222222")


# ---- draw flow boxes -------------------------------------------------------
for i, (title, nline) in enumerate(flow):
    title_parts = title.split("\n")
    lines = [(p, TITLE_FS, True) for p in title_parts] + [(nline, BODY_FS, True)]
    # recompute height from actual line count
    h = box_height(len(lines))
    flow_h[i] = h
    fill = FINAL_FILL if i == len(flow) - 1 else BLUE_FILL
    draw_box(L_CX, flow_cy[i], L_W, h, lines, fill, BLUE_EDGE, align="center")

# recompute centres with final heights (kept uniform spacing S, fine)

# ---- vertical arrows between flow boxes ------------------------------------
for i in range(3):
    y_top = flow_cy[i] - flow_h[i] / 2
    y_bot = flow_cy[i + 1] + flow_h[i + 1] / 2
    ax.add_patch(FancyArrowPatch((L_CX, y_top), (L_CX, y_bot),
                 arrowstyle="-|>", mutation_scale=16, lw=1.6,
                 color=ARROW_C, shrinkA=0, shrinkB=0))

# ---- exclusion boxes + branch arrows ---------------------------------------
for i, items in enumerate(excl):
    cy = (flow_cy[i] + flow_cy[i + 1]) / 2.0          # midway in the gap
    h = box_height(len(items))
    lines = [(items[0], BODY_FS, True)] + [(t, BODY_FS, False) for t in items[1:]]
    draw_box(R_X0 + R_W / 2, cy, R_W, h, lines, PINK_FILL, PINK_EDGE,
             align="left")
    # branch: from vertical connector out to the exclusion box
    ax.add_patch(FancyArrowPatch((L_CX, cy), (R_X0, cy),
                 arrowstyle="-|>", mutation_scale=16, lw=1.6,
                 color=ARROW_C, shrinkA=0, shrinkB=0))

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ROOT = Path(__file__).resolve().parents[2]
out = str(ROOT / "ITT_Analysis" / "results" / "Figure_S1_consort.png")
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.10)
print("wrote", out)
