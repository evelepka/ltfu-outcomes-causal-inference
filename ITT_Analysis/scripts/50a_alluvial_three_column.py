"""
50a_alluvial_three_column.py
----------------------------
Three-column Sankey/alluvial of TB trajectories after loss to follow-up,
matching the style we used in earlier decks.

    Col 0: Index Loss to Follow Up (all LTFU individuals)
    Col 1: First outcome after LTFU
             - Died (without retreatment)
             - Retreatment
             - No further outcome available
    Col 2: Retreatment episode outcome (only for retreated branch)
             - Cured
             - Lost to follow up (abandoned again)
             - Died
             - Failure
             - Other/no further outcome available

Reads:
  - ITT_Analysis/data/itt_cohort.csv       (for LTFU membership + index end_date)
  - Data/Final_table_cleaned.csv           (to look up retreatment episode outcome)
Writes:
  - ITT_Analysis/results/Figure_1a_alluvial.png
  - ITT_Analysis/results/Figure_1a_alluvial.html   (interactive)
"""

import os
import sys
from pathlib import Path
from collections import defaultdict
import warnings

import numpy as np
import pandas as pd
import plotly.graph_objects as go


warnings.filterwarnings("ignore")


def _find_project_root():
    env = os.environ.get("TB_ABANDONMENT_ROOT")
    if env:
        p = Path(env).expanduser()
        if p.exists():
            return p
    for c in [
        Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
        Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
    ]:
        if c.exists():
            return c
    return Path(__file__).resolve().parents[2]


BASE = _find_project_root()
COHORT_CSV = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"
FULL_CSV = BASE / "Data" / "Final_table_cleaned.csv"
OUT_PNG = BASE / "ITT_Analysis" / "results" / "Figure_1a_alluvial.png"
OUT_HTML = BASE / "ITT_Analysis" / "results" / "Figure_1a_alluvial.html"

CENSOR_DATE = pd.Timestamp("2024-12-31")
MIN_FLOW = 5  # hide tiny flows

ABANDON_OUTCOMES = {"Abandono", "Abandono Primario", "Faltoso"}
OBITO_OUTCOMES = {"Obito TB", "Obito NTB"}

# Palette — muted pastels with strong accent for each outcome
NODE_COLORS = {
    "Index Loss to Follow Up":    "rgba(142,  68, 173, 0.90)",  # purple
    "Retreatment":                "rgba( 41, 128, 185, 0.90)",  # blue
    "Died (without retreatment)": "rgba( 44,  62,  80, 0.90)",  # dark
    "No further outcome available":"rgba(176, 190, 197, 0.90)", # grey
    "Lost to follow up":          "rgba(231,  76,  60, 0.90)",  # red
    "Cured":                      "rgba( 46, 204, 113, 0.90)",  # green
    "Died":                       "rgba( 44,  62,  80, 0.90)",  # dark
    "Failure":                    "rgba(230, 126,  34, 0.90)",  # orange
    "Other/no further outcome available": "rgba(176, 190, 197, 0.90)",
}
LINK_COLORS = {k: v.replace("0.90", "0.32") for k, v in NODE_COLORS.items()}


# ---------------------------------------------------------------------------
# Load cohort (LTFU subset) + raw dataset
# ---------------------------------------------------------------------------
print("Loading cohort ...")
cohort = pd.read_csv(COHORT_CSV, low_memory=False)
cohort["end_date"] = pd.to_datetime(cohort["end_date"], errors="coerce")
cohort["death_date"] = pd.to_datetime(cohort["death_date"], errors="coerce")
ltfu = cohort[cohort["itt_group"] == "Loss to follow-up"].copy()
print(f"  LTFU N = {len(ltfu):,}")

print("Loading raw dataset ...")
full = pd.read_csv(FULL_CSV, low_memory=False)
full["notification_date"] = pd.to_datetime(full["notification_date"], errors="coerce")
full["end_date"] = pd.to_datetime(full["end_date"], errors="coerce")
full = full.sort_values(["sinan_clean", "notification_date"])
# Filter to people in our cohort for speed
full = full[full["sinan_clean"].astype(str).isin(ltfu["sinan_clean"].astype(str))].copy()
print(f"  Raw rows for LTFU individuals: {len(full):,}")


# ---------------------------------------------------------------------------
# Identify retreatment episode + classify its outcome
# ---------------------------------------------------------------------------
# For each LTFU individual, find the FIRST notification AFTER their index end_date.
# If none: "No further outcome available" at Col 1 (unless they died first).
# If the death_date precedes any such notification: "Died (without retreatment)".
# Otherwise: "Retreatment", and classify its case_outcome.

def classify_retreat_outcome(outcome: str, ep_end, death_date) -> str:
    o = (str(outcome) if outcome is not None else "").strip()
    if o in ABANDON_OUTCOMES:
        return "Lost to follow up"
    if o == "Cura":
        return "Cured"
    if o in OBITO_OUTCOMES:
        return "Died"
    if o == "Falencia/Resistencia":
        return "Failure"
    # Transfer / other / NaN — fold death-during-episode into Died if applicable
    if pd.notna(death_date) and death_date <= CENSOR_DATE:
        ref = ep_end if pd.notna(ep_end) else pd.NaT
        if pd.notna(ref) and death_date > ref:
            return "Died"
    return "Other/no further outcome available"


full_by_id = {str(sid): grp for sid, grp in full.groupby("sinan_clean")}


def col1_and_retx_outcome(row):
    sid = str(row["sinan_clean"])
    idx_end = row["end_date"]
    death = row["death_date"]
    grp = full_by_id.get(sid, pd.DataFrame())
    later = grp[grp["notification_date"] > idx_end].sort_values("notification_date") if len(grp) else pd.DataFrame()
    if len(later) == 0:
        # No retreatment episode on file
        if pd.notna(death) and death <= CENSOR_DATE and death > idx_end:
            return pd.Series(["Died (without retreatment)", None])
        return pd.Series(["No further outcome available", None])
    # Retreatment exists — but did they die BEFORE retreatment notification?
    first_retx_date = later.iloc[0]["notification_date"]
    if pd.notna(death) and death <= CENSOR_DATE and death <= first_retx_date:
        return pd.Series(["Died (without retreatment)", None])
    ep = later.iloc[0]
    retx_out = classify_retreat_outcome(ep.get("case_outcome"), ep.get("end_date"), death)
    return pd.Series(["Retreatment", retx_out])


print("Classifying flows ...")
ltfu[["col1", "col2"]] = ltfu.apply(col1_and_retx_outcome, axis=1)

# Sanity print
print("\nCol 1 distribution:")
for v, c in ltfu["col1"].value_counts().items():
    print(f"  {v:<35}: {c:>6,}  ({c/len(ltfu)*100:.1f}%)")
retreated = ltfu[ltfu["col1"] == "Retreatment"]
print("\nCol 2 distribution (within retreated):")
for v, c in retreated["col2"].value_counts().items():
    print(f"  {v:<38}: {c:>6,}  ({c/len(retreated)*100:.1f}%)")


# ---------------------------------------------------------------------------
# Build Sankey
# ---------------------------------------------------------------------------
TOTAL = len(ltfu)

ALL_NODES_ORDERED = [
    (0, "Index Loss to Follow Up"),
    (1, "Retreatment"),
    (1, "Died (without retreatment)"),
    (1, "No further outcome available"),
    (2, "Lost to follow up"),
    (2, "Cured"),
    (2, "Died"),
    (2, "Failure"),
    (2, "Other/no further outcome available"),
]
node_idx = {n: i for i, n in enumerate(ALL_NODES_ORDERED)}
node_size = defaultdict(int)
flows = defaultdict(int)

for _, row in ltfu.iterrows():
    node_size[(0, "Index Loss to Follow Up")] += 1
    node_size[(1, row["col1"])] += 1
    flows[((0, "Index Loss to Follow Up"), (1, row["col1"]))] += 1

for _, row in retreated.iterrows():
    node_size[(2, row["col2"])] += 1
    flows[((1, "Retreatment"), (2, row["col2"]))] += 1

X_POSITIONS = {0: 0.01, 1: 0.5, 2: 0.99}

node_labels, node_colors_arr, node_x = [], [], []
for col, lab in ALL_NODES_ORDERED:
    count = node_size.get((col, lab), 0)
    if count == 0:
        node_labels.append("")
    else:
        pct = count / TOTAL * 100
        if col == 2:
            denom = node_size.get((1, "Retreatment"), TOTAL)
            pct2 = count / denom * 100 if denom else 0
            node_labels.append(f"<b>{lab}</b><br>{count:,} ({pct2:.1f}%)")
        elif col == 1:
            node_labels.append(f"<b>{lab}</b><br>{count:,} ({pct:.1f}%)")
        else:
            node_labels.append(f"<b>{lab}</b><br>N = {count:,}")
    node_colors_arr.append(NODE_COLORS.get(lab, "rgba(150,150,150,0.85)"))
    node_x.append(X_POSITIONS[col])

link_src, link_dst, link_val, link_col, link_lab = [], [], [], [], []
for (src_node, dst_node), cnt in flows.items():
    if cnt < MIN_FLOW:
        continue
    link_src.append(node_idx[src_node])
    link_dst.append(node_idx[dst_node])
    link_val.append(cnt)
    link_col.append(LINK_COLORS.get(dst_node[1], "rgba(150,150,150,0.25)"))
    link_lab.append(f"{src_node[1]} → {dst_node[1]}<br>{cnt:,} individuals ({cnt/TOTAL*100:.1f}%)")

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        pad=22, thickness=24,
        line=dict(color="white", width=0.8),
        label=node_labels,
        color=node_colors_arr,
        x=node_x,
        hovertemplate="%{label}<extra></extra>",
    ),
    link=dict(
        source=link_src, target=link_dst, value=link_val,
        color=link_col, label=link_lab,
        hovertemplate="%{label}<extra></extra>",
    ),
))

# Column headers
headers = {
    0: "Index Loss to Follow Up",
    1: "First Outcome After Loss to Follow Up",
    2: "Retreatment Episode Outcome",
}
annotations = [
    dict(x=X_POSITIONS[c], y=1.06, xref="paper", yref="paper",
         text=f"<b>{t}</b>", showarrow=False,
         font=dict(size=13, color="#1a1a2e", family="Arial"),
         xanchor="center", yanchor="bottom")
    for c, t in headers.items()
]
annotations.append(dict(
    x=0.0, y=-0.08, xref="paper", yref="paper",
    text=(f"<i>N = {TOTAL:,} individuals with loss to follow-up on their first new (Novo) TB episode. "
          f"Retreatment episode outcomes draw from the subsequent notification in SINAN "
          f"(case_outcome field). Follow-up censored {CENSOR_DATE.date()}. "
          f"Flows &lt; {MIN_FLOW} hidden.</i>"),
    showarrow=False, font=dict(size=10, family="Arial", color="#555"),
    xanchor="left", yanchor="top"
))

fig.update_layout(
    title=dict(
        text=f"<b>TB trajectories after loss to follow-up</b>"
             f"<br><sup>São Paulo | LTFU cohort N = {TOTAL:,} | follow-up through {CENSOR_DATE.date()}</sup>",
        x=0.5, xanchor="center",
        font=dict(size=18, family="Arial", color="#1a1a2e")
    ),
    font=dict(size=11, family="Arial", color="#2c2c2c"),
    paper_bgcolor="#F8F9FA",
    plot_bgcolor="#F8F9FA",
    height=650, width=1500,
    margin=dict(l=30, r=30, t=120, b=140),
    annotations=annotations,
)

fig.write_html(str(OUT_HTML), include_plotlyjs="cdn", full_html=True)
print(f"\nWrote {OUT_HTML}")
try:
    fig.write_image(str(OUT_PNG), width=1600, height=800, scale=2)
    print(f"Wrote {OUT_PNG}")
except Exception as e:
    print(f"PNG export skipped: {e}")
