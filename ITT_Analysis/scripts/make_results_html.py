"""
make_results_html.py
--------------------
Build a single self-contained HTML report summarizing all figures and tables
produced for the ITT analysis. Embeds PNGs as base64 so the file can be
shared without any asset paths.

Output: ITT_Analysis/results/results_summary.html
"""

import base64
import os
from pathlib import Path
from datetime import datetime

import pandas as pd


# ---------------------------------------------------------------------------
# Project root resolution
# ---------------------------------------------------------------------------
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
RESULTS = BASE / "ITT_Analysis" / "results"
DATA = BASE / "Data"
COHORT_CSV = BASE / "ITT_Analysis" / "data" / "itt_cohort.csv"

OUT_HTML = RESULTS / "results_summary.html"


def embed_png(path: Path) -> str:
    """Return data: URI for a PNG."""
    if not path.exists():
        return ""
    b64 = base64.b64encode(path.read_bytes()).decode()
    return f"data:image/png;base64,{b64}"


def df_to_html(df: pd.DataFrame, classes: str = "tbl") -> str:
    if df is None or df.empty:
        return '<p class="muted">(no data)</p>'
    return df.to_html(classes=classes, index=False, border=0, escape=False,
                      float_format=lambda x: f"{x:.3g}")


def read_csv_safe(path: Path):
    if not path.exists():
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Pull together the pieces
# ---------------------------------------------------------------------------
cohort = pd.read_csv(COHORT_CSV, low_memory=False, usecols=["itt_group"])
N_total = len(cohort)
itt_counts = cohort["itt_group"].value_counts().to_dict()
N_ltfu = itt_counts.get("Loss to follow-up", 0)
N_nonltfu = itt_counts.get("Non-LTFU", 0)

flowchart = read_csv_safe(DATA / "exclusion_flowchart.csv")
fig1_stats = read_csv_safe(RESULTS / "Figure_1_caption_stats.csv")
fig2_vals = read_csv_safe(RESULTS / "Figure_2_cif_values_24mo.csv")
tt_array = read_csv_safe(RESULTS / "target_trial_mi_6mo_array_hr.csv")
tt_subgroups = read_csv_safe(RESULTS / "target_trial_subgroup_interactions_mi.csv")
mi_cc = read_csv_safe(RESULTS / "multivariable_results_mi_cc.csv")

fig1_src = embed_png(RESULTS / "Figure_1_descriptive.png")
fig2_src = embed_png(RESULTS / "Figure_2_stratified_retreatment_24mo.png")
fig3_src = embed_png(RESULTS / "Figure_3_causal_mortality.png")

# Quick summary stats for the top banner
if fig2_vals is not None:
    fig2_vals_fmt = fig2_vals.copy()
    for c in ["month_6", "month_12", "month_18", "month_24"]:
        if c in fig2_vals_fmt.columns:
            fig2_vals_fmt[c] = fig2_vals_fmt[c].apply(lambda x: f"{x:.1f}%")
else:
    fig2_vals_fmt = None

if tt_array is not None:
    tt_array_fmt = tt_array.copy()
    tt_array_fmt["HR (95% CI)"] = tt_array_fmt.apply(
        lambda r: f"{r['HR']:.2f} ({r['CI_Lower']:.2f}–{r['CI_Upper']:.2f})", axis=1
    )
    tt_array_fmt["p-value"] = tt_array_fmt["P_Value"].apply(
        lambda p: "<0.001" if p < 0.001 else f"{p:.3f}"
    )
    tt_array_fmt = tt_array_fmt[["Trial_Month", "HR (95% CI)", "p-value"]]
else:
    tt_array_fmt = None

if tt_subgroups is not None:
    tt_sub_fmt = tt_subgroups.copy()
    tt_sub_fmt["HR (95% CI)"] = tt_sub_fmt.apply(
        lambda r: f"{r['HR']:.2f} ({r['CI_L']:.2f}–{r['CI_H']:.2f})", axis=1
    )
    tt_sub_fmt["p-value"] = tt_sub_fmt["P_Value"].apply(
        lambda p: "<0.001" if p < 0.001 else f"{p:.3f}"
    )
    tt_sub_fmt = tt_sub_fmt[["Subgroup", "Level", "HR (95% CI)", "p-value", "N_Imp"]]
else:
    tt_sub_fmt = None

# Split the mi_cc table by model for display
mi_cc_panels = {}
if mi_cc is not None:
    for model in mi_cc["model"].unique():
        mi_cc_panels[model] = mi_cc[mi_cc["model"] == model][["term", "estimate", "p_value"]]

generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
CSS = """
:root {
  --ink: #1f2937;
  --muted: #6b7280;
  --accent: #2c3e50;
  --accent-2: #e74c3c;
  --bg-panel: #f8fafc;
  --border: #e5e7eb;
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  max-width: 1200px; margin: 0 auto; padding: 2rem;
  color: var(--ink); line-height: 1.55; background: #fff;
}
h1 { font-size: 1.75rem; margin: 0 0 .25rem; }
h2 { font-size: 1.35rem; margin-top: 2.5rem; padding-bottom: .35rem;
     border-bottom: 2px solid var(--accent); color: var(--accent); }
h3 { font-size: 1.1rem; margin-top: 1.5rem; color: var(--accent); }
.meta { color: var(--muted); margin-bottom: 2rem; }
.banner {
  background: var(--bg-panel); border: 1px solid var(--border);
  border-radius: 8px; padding: 1rem 1.25rem; margin: 1rem 0 2rem;
}
.banner .stat { display: inline-block; margin-right: 2.5rem; }
.banner .stat .num { font-size: 1.5rem; font-weight: 700; color: var(--accent-2); }
.banner .stat .lbl { color: var(--muted); font-size: .85rem; display: block; }
figure { margin: 1.5rem 0; }
figure img {
  width: 100%; height: auto; display: block;
  border: 1px solid var(--border); border-radius: 6px;
  background: #fff;
}
figcaption {
  text-align: center; color: var(--muted); font-size: .9rem;
  margin-top: .4rem;
}
.tbl {
  border-collapse: collapse; width: 100%; margin: .75rem 0 1.5rem;
  font-size: .9rem;
}
.tbl th, .tbl td {
  text-align: left; padding: 6px 10px;
  border-bottom: 1px solid var(--border);
}
.tbl th {
  background: var(--bg-panel); font-weight: 600;
  border-bottom: 2px solid var(--accent);
}
.tbl tr:hover td { background: #fafbfc; }
.muted { color: var(--muted); }
.note {
  font-size: .85rem; color: var(--muted); padding: .5rem .75rem;
  background: #fffbeb; border-left: 3px solid #f59e0b; border-radius: 3px;
  margin: .75rem 0;
}
details { margin: .5rem 0; }
details summary {
  cursor: pointer; font-weight: 600; color: var(--accent);
  padding: .35rem 0;
}
code { background: var(--bg-panel); padding: 1px 5px; border-radius: 3px; font-size: .9em; }
nav { background: var(--bg-panel); padding: .75rem 1.25rem;
      border-radius: 6px; margin: 1rem 0 2rem;
      border: 1px solid var(--border); }
nav a { margin-right: 1.25rem; color: var(--accent); text-decoration: none;
        font-weight: 500; }
nav a:hover { text-decoration: underline; }
.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
@media (max-width: 800px) { .two-col { grid-template-columns: 1fr; } }
"""

def section(title, anchor, html_inner):
    return f'<section id="{anchor}"><h2>{title}</h2>\n{html_inner}\n</section>'


html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>TB Abandonment — Results Summary</title>
<style>{CSS}</style>
</head>
<body>

<h1>Outcomes After TB Treatment Abandonment</h1>
<p class="meta">Results summary &middot; generated {generated_at} &middot;
cohort N = {N_total:,} ({N_ltfu:,} LTFU, {N_nonltfu:,} Non-LTFU)</p>

<div class="banner">
  <div class="stat"><span class="num">{N_total:,}</span><span class="lbl">Total ITT cohort</span></div>
  <div class="stat"><span class="num">{N_ltfu:,}</span><span class="lbl">Loss to follow-up</span></div>
  <div class="stat"><span class="num">{N_nonltfu:,}</span><span class="lbl">Non-LTFU (maintained care)</span></div>
  <div class="stat"><span class="num">{100*N_ltfu/N_total:.1f}%</span><span class="lbl">LTFU proportion</span></div>
</div>

<nav>
  <a href="#cohort">Cohort</a>
  <a href="#fig1">Figure 1 — Descriptive</a>
  <a href="#fig2">Figure 2 — Retreatment</a>
  <a href="#fig3">Figure 3 — Mortality</a>
  <a href="#mi">MI multivariable</a>
</nav>

{section("Cohort & exclusion flowchart", "cohort",
    '<p>Single cohort file at <code>ITT_Analysis/data/itt_cohort.csv</code>. '
    'See <code>COHORT_SELECTION.md</code> for the full definition. '
    'Inclusion/exclusion chain below is auto-generated by <code>01_itt_cohort_selection.py</code>.</p>'
    + df_to_html(flowchart))}

{section("Figure 1 — Descriptive", "fig1", f'''
<figure>
  <img src="{fig1_src}" alt="Figure 1">
  <figcaption>Figure 1. A — Clinical flow after loss to follow-up. B — Event-timing rainclouds.</figcaption>
</figure>
<h3>Summary statistics</h3>
{df_to_html(fig1_stats)}
''')}

{section("Figure 2 — Stratified cumulative incidence of retreatment (24 months)", "fig2", f'''
<figure>
  <img src="{fig2_src}" alt="Figure 2">
  <figcaption>Figure 2. Cumulative incidence of retreatment after abandonment, up to 24 months. Aalen–Johansen estimates treating death as a competing risk. Gray's test p-values annotated on each panel.</figcaption>
</figure>
<h3>Cumulative incidence values by horizon (%)</h3>
{df_to_html(fig2_vals_fmt)}
''')}

{section("Figure 3 — Causal mortality panels", "fig3", f'''
<figure>
  <img src="{fig3_src}" alt="Figure 3">
  <figcaption>Figure 3. A — Crude time-varying cumulative hazard illustrating immortal-time bias. B — MI-pooled sequential target-trial hazard ratio by month of abandonment. C — MI-pooled subgroup-stratified target-trial hazard ratios (forest plot).</figcaption>
</figure>
<div class="two-col">
  <div>
    <h3>B. Target-trial HR by month of abandonment (MI-pooled)</h3>
    {df_to_html(tt_array_fmt)}
  </div>
  <div>
    <h3>C. Subgroup-stratified target-trial HR (MI-pooled)</h3>
    {df_to_html(tt_sub_fmt)}
    <p class="note">These MI-pooled estimates differ from earlier complete-case numbers reported in <code>README.md</code>. MI pulls estimates toward the null with tighter CIs; complete-case dropped ~45% of rows with any missingness. Both are legitimate depending on assumed missingness mechanism.</p>
  </div>
</div>
''')}
""" + section(
    "Within-LTFU multivariable models (MI-pooled + complete-case sensitivity)",
    "mi",
    '<p>Pooled Cox (mortality) and Fine–Gray (retreatment with death as competing risk) '
    'within the LTFU subgroup (N = 21,619). Imputation on the full cohort via '
    '<code>miceforest</code>; pooling via Rubin\'s rules in R.</p>\n'
    + "".join(
        f'<details><summary>{model} ({len(tbl)} terms)</summary>\n{df_to_html(tbl)}\n</details>'
        for model, tbl in sorted(mi_cc_panels.items())
    )
    if mi_cc_panels
    else '<p class="muted">(no MI results on disk yet)</p>'
) + """

<p class="meta" style="margin-top:3rem">
Generated from the ITT pipeline at <code>github.com/jasonandr/outcomes-after-tb-abandonment</code>.
</p>

</body></html>
"""

OUT_HTML.write_text(html, encoding="utf-8")
print(f"Wrote {OUT_HTML}")
print(f"Size: {OUT_HTML.stat().st_size / 1024:.1f} KB")
