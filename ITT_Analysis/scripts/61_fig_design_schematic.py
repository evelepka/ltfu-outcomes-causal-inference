#!/usr/bin/env python3
"""
Design schematic for the LTFU paper (R1 comment 15, R2 comment 1).

One figure, five blocks:

  A  what the surveillance data actually record, on an expanded 0-8 month axis:
     four patient tracks, including the patient who dies inside the 30-day
     absence window and is therefore recorded as an on-treatment death (R3's
     competing-exposure objection) and the comparator patient who is in care at
     someone else's landmark and is lost later (R2's point about the reference
     group).
  B  comparator = treatment completers             -> healthy-survivor comparator
  C  exposure classified from treatment initiation -> immortal time
  D  rolling landmark at each patient's own LTFU declaration date  [PRIMARY]
  E  clone-censor-weight from day 0                -> the bias analysis

Panels B-E share one x-axis (months since treatment initiation) so that the
whole point of the figure -- these designs start the clock in different places
-- is visible without reading any text. Each panel names its own bias in the
plot and its time origin and estimand in the key beside it, which is what R1's
pivotal comment 2 asked for.

The curves are SCHEMATIC. They illustrate the direction of each bias and are
not fitted estimates; nothing here is read from a results CSV, so this script
cannot drift from the analysis. Point estimates are OFF by default and, when
switched on with --numbers, are quoted only from docs/number-registry.csv ids
(see ESTIMATES below). Per ADR-0005 the early/late split is NOT presented as a
reported estimand: the windows appear only to explain bias.

Canvas is 7.5 x 8.75 in, the PLOS Medicine maximum, rendered at 300 dpi.
Palette is Okabe-Ito; PLOS asked us to avoid red and green.

Usage
-----
    python3 ITT_Analysis/scripts/61_fig_design_schematic.py
    python3 ITT_Analysis/scripts/61_fig_design_schematic.py --numbers

Writes SVG (editable), PDF, PNG and 300-dpi TIFF to
Plos Medicine/R1/Figures/Draft/.
"""

from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

# ---------------------------------------------------------------- paths -----
ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
OUTDIR = ROOT / "Plos Medicine" / "R1" / "Figures" / "Draft"
STEM = "LTFU_design_schematic"

# --------------------------------------------------------------- palette ----
BLUE = "#0072B2"  # remained in care / comparator
VERM = "#D55E00"  # lost to follow-up / exposed
PURP = "#A0498C"  # the misclassification message
AMBER = "#E69F00"  # bias shading
AMBER_INK = "#7A5000"
BLUE_INK = "#004E7A"
INK = "#1A1A1A"
GREY = "#767676"

# ----------------------------------------------------------------- type -----
for family in ("Helvetica", "Arial", "DejaVu Sans"):
    if any(f.name == family for f in matplotlib.font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = family
        break
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "svg.fonttype": "none",  # keep text editable in the SVG
        "pdf.fonttype": 42,
    }
)

FS_TITLE = 8.0
FS_BODY = 6.5  # right-hand keys
FS_ANNOT = 6.2  # in-panel annotation
FS_TICK = 6.4
LINE = 1.40  # line-height multiplier for hand-stacked text

FIG_W, FIG_H = 7.5, 10.0

# --------------------------------------------------------------- geometry ---
A_BOX = [0.078, 0.777, 0.907, 0.173]  # panel A axes
ROW_TOPS = [0.742, 0.604, 0.466, 0.328, 0.190]  # top of each of panels B-F
PLOT_W, PLOT_H = 0.505, 0.0875
KEY_X, KEY_W, KEY_H = 0.605, 0.380, 0.118
KEY_WRAP = 62
KEY_LABEL_W = 0.30  # field-name column, as a fraction of the key axes width
KEY_VAL_WRAP = 42
# Opaque backing for any label that sits over a curve, a rule or a shaded band.
TBOX = dict(facecolor="white", edgecolor="none", pad=0.9)

# --------------------------------------------------- timeline event times ---
GAP = 30 / 30.4375  # the 30 days the LTFU definition requires, in months
DECL = 3.4  # LTFU declaration date (observed in TBweb)
LAST_DOSE = DECL - GAP  # inferred by subtracting 30 days
DEATH_IN_GAP = 3.0
P4_STOP = 5.0
PLANNED_END = 6.0
T_DETAIL = 8.0  # right edge of the timeline tracks in panel A
T_MAX = 24.0  # right edge of panels B-F
T_CLONE = 4.0  # the "become LTFU by month T" horizon in panel E
YMAX = 1.75  # panel B-F y-limit; curves stay under 0.95, text sits above 1.0

# ------------------------------------------------------------- estimates ----
# Quoted ONLY with --numbers. Bracketed keys are docs/number-registry.csv ids.
# UPDATED 2026-08-20 for the 2026-08-19 design decision: CCW carries the main
# effect and the subgroups, the rolling landmark carries timing. The PRIMARY /
# BIAS labels on panels D and E were the other way round and are now swapped.
# Panel E moves to the 60-month horizon the decision reports.
# 2026-08-21: panels B-D are the three prior-literature designs, each REPRODUCED
# in this cohort by 49_literature_design_reproduction.R, so the spread is shown on
# our own data rather than by cross-study arithmetic. The previous version had no
# panel for the initiation-origin/completer-comparator design -- the one that
# produces the published 4-5 fold estimates -- and mis-cited roll_early 0.89 for
# the all-non-LTFU design, whose reproduced value is 0.35.
ESTIMATES = {
    "B": "aHR 4.97 (4.48-5.51), reproduced here",
    "C": "aHR 3.55 (3.22-3.92), reproduced here",
    "D": "aHR 0.35 (0.32-0.37), reproduced here",
    "E": "aHR 2.42 (2.18-2.69); RD +2.61 pp at 5 y",
    "F": "RD +2.22 pp (1.88-2.52); RR 1.24 (1.20-1.27)",
}


# ---------------------------------------------------------------- helpers ---
def cuminc(x, x0, plateau, rate):
    """Schematic cumulative-incidence curve rising from origin x0."""
    return plateau * (1.0 - np.exp(-rate * np.maximum(np.asarray(x, float) - x0, 0.0)))


def panel_title(ax, letter, title, tag=None):
    ax.text(
        0.0, 1.04, letter, transform=ax.transAxes, fontsize=FS_TITLE + 1,
        weight="bold", va="bottom", ha="left",
    )
    ax.text(
        0.032, 1.04, title, transform=ax.transAxes, fontsize=FS_TITLE,
        weight="bold", va="bottom", ha="left",
    )
    if tag:
        ax.text(
            1.0, 1.04, tag, transform=ax.transAxes, fontsize=FS_ANNOT,
            weight="bold", color=VERM, va="bottom", ha="right",
        )


def curve_axis(ax, xlabel=False):
    ax.set_xlim(-0.7, T_MAX + 0.5)
    ax.set_ylim(0, YMAX)
    ax.add_patch(
        Rectangle((0, 0), PLANNED_END, YMAX, facecolor="#F1F1F1", edgecolor="none", zorder=0)
    )
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=FS_TICK, length=2.5, pad=1.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(GREY)
    ax.spines["bottom"].set_color(GREY)
    ax.set_ylabel("Cumulative\nmortality risk", fontsize=FS_TICK, labelpad=3)
    if xlabel:
        ax.set_xlabel("Months since treatment initiation", fontsize=FS_TICK, labelpad=2)
    ax.annotate(  # direction-only arrow: the vertical scale is schematic
        "",
        xy=(-0.45, 0.86),
        xytext=(-0.45, 0.06),
        arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.7, mutation_scale=6),
        annotation_clip=False,
    )


def bias_note(ax, text, color, width=46):
    """The panel's bias statement, in the band reserved above the curves."""
    ax.text(
        T_MAX + 0.4,
        YMAX * 0.995,
        "\n".join(textwrap.wrap(text, width)),
        fontsize=FS_ANNOT,
        color=color,
        ha="right",
        va="top",
        linespacing=1.35,
        weight="bold",
        bbox=TBOX,
        zorder=6,
    )


def landmark_line(ax, x, label, top=0.99, label_y=None):
    ax.plot([x, x], [0, top], ls=(0, (1.6, 1.6)), lw=0.85, color=INK, zorder=5)
    ax.text(
        x + 0.25, (top if label_y is None else label_y) + 0.03, label,
        fontsize=FS_ANNOT, color=INK,
        ha="left", va="bottom", weight="bold", linespacing=1.3, bbox=TBOX, zorder=6,
    )


def series_label(ax, y, text, color, va="bottom"):
    """Label a curve just clear of it, not on top of it.

    Sitting the text ON the line and giving it a white background punches a
    visible gap in the curve, which reads as a rendering fault. Instead the label
    is nudged off the line and pulled left of the right edge, into the space
    between the two diverging curves, so its backing only covers empty ground.
    """
    dy = 0.055 if va == "bottom" else -0.055
    ax.text(
        19.8, y + dy, text, fontsize=FS_ANNOT, color=color, ha="right", va=va,
        bbox=TBOX, zorder=6,
    )


def key_block(ax, rows):
    """Right-hand key as a two-column table: field name | value.

    Columns sit at fixed x in every panel and each panel's block opens with a
    rule, so the five keys read as one table banded by panel rather than as five
    paragraphs of prose. Values are phrases, not sentences -- anything needing a
    sentence belongs in the caption.
    """
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    lh = (FS_BODY * LINE / 72.0) / (KEY_H * FIG_H)  # one line, in axes fraction
    ax.plot([0, 1], [1.0, 1.0], color=GREY, lw=0.5, clip_on=False, zorder=1)
    y = 1.0 - lh * 0.55
    for field, value in rows:
        body = textwrap.wrap(value, KEY_VAL_WRAP)
        ax.text(0, y, field, fontsize=FS_BODY, weight="bold", color=INK, va="top", ha="left")
        ax.text(
            KEY_LABEL_W, y, "\n".join(body), fontsize=FS_BODY, color=INK,
            va="top", ha="left", linespacing=LINE,
        )
        y -= lh * (len(body) + 0.30)


# ================================================================ panel A ===
ROW = {"p1": 8.8, "p2": 6.6, "p3": 4.4, "p4": 2.2}
BAR_H = 0.42
LH_A = (FS_ANNOT * 1.35 / 72.0) / (A_BOX[3] * FIG_H) * 10.0  # one line, in data units


def draw_timeline(ax):
    ax.set_xlim(-3.5, 11.6)
    ax.set_ylim(0, 11.3)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, T_DETAIL + 1, 1))
    ax.tick_params(axis="x", labelsize=FS_TICK, length=2.5, pad=1.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_bounds(0, T_DETAIL)
    ax.spines["bottom"].set_color(GREY)
    ax.text(
        (T_DETAIL / 2 + 3.5) / 15.1, -0.105, "Months since treatment initiation",
        transform=ax.transAxes, fontsize=FS_TICK, color=INK, ha="center", va="top",
    )
    ax.add_patch(
        Rectangle((0, 0), PLANNED_END, 9.35, facecolor="#F1F1F1", edgecolor="none", zorder=0)
    )
    # the declaration date, carried down across every track
    ax.plot(
        [DECL, DECL], [0.55, 9.30], ls=(0, (1.0, 2.0)), lw=0.7, color=GREY, alpha=0.9, zorder=1
    )

    def oncare(y, x0, x1):
        ax.add_patch(
            Rectangle((x0, y - BAR_H / 2), x1 - x0, BAR_H, facecolor=BLUE,
                      edgecolor="none", zorder=2)
        )

    def lost(y, x0, x1):
        ax.add_patch(
            Rectangle((x0, y - BAR_H / 2), x1 - x0, BAR_H, facecolor=VERM,
                      edgecolor="none", zorder=2)
        )
        ax.add_patch(
            FancyArrowPatch((x1, y), (x1 + 0.55, y), arrowstyle="-|>",
                            mutation_scale=7, lw=1.5, color=VERM)
        )

    def absence(y, x0, x1):
        ax.add_patch(
            Rectangle(
                (x0, y - BAR_H / 2), x1 - x0, BAR_H, facecolor="white",
                edgecolor=GREY, lw=0.6, hatch="/////", zorder=3,
            )
        )

    def rowlabel(y, text):
        ax.text(
            -0.35, y, text, fontsize=FS_ANNOT, color=INK, weight="bold",
            ha="right", va="center", linespacing=1.3,
        )

    # -- patient 1 -----------------------------------------------------------
    y = ROW["p1"]
    rowlabel(y, "Patient 1\ncompletes treatment")
    oncare(y, 0, PLANNED_END)
    ax.plot([PLANNED_END], [y], marker="D", ms=3.8, color=INK, zorder=6)
    ax.text(
        PLANNED_END + 0.35, y, "treatment completed", fontsize=FS_ANNOT, color=INK, va="center",
    )

    # -- patient 2 -----------------------------------------------------------
    y = ROW["p2"]
    rowlabel(y, "Patient 2\nlost to follow-up")
    oncare(y, 0, LAST_DOSE)
    absence(y, LAST_DOSE, DECL)
    lost(y, DECL, T_DETAIL)
    ax.plot([LAST_DOSE], [y], marker="o", ms=3.8, mfc="white", mec=INK, mew=0.9, zorder=6)
    ax.plot([DECL], [y], marker="D", ms=3.8, color=INK, zorder=6)
    ax.annotate(
        "last dose (inferred: declaration - 30 d)",
        xy=(LAST_DOSE, y + BAR_H / 2),
        xytext=(LAST_DOSE - 0.30, y + 0.46),
        fontsize=FS_ANNOT, color=INK, ha="center", va="bottom", linespacing=1.3,
        arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY, shrinkB=1.5),
    )
    ax.annotate(
        "LTFU declared (the date TBweb records)",
        xy=(DECL, y + BAR_H / 2),
        xytext=(DECL + 2.6, y + 0.46),
        fontsize=FS_ANNOT, color=INK, ha="center", va="bottom", weight="bold", linespacing=1.3,
        arrowprops=dict(arrowstyle="-", lw=0.6, color=GREY, shrinkB=1.5),
    )

    # -- patient 3 -----------------------------------------------------------
    y = ROW["p3"]
    rowlabel(y, "Patient 3\nstops attending,\ndies before declaration")
    oncare(y, 0, LAST_DOSE)
    absence(y, LAST_DOSE, DEATH_IN_GAP)
    ax.plot([LAST_DOSE], [y], marker="o", ms=3.8, mfc="white", mec=INK, mew=0.9, zorder=6)
    ax.text(
        LAST_DOSE - 0.15, y + 0.40, "stops attending", fontsize=FS_ANNOT, color=INK,
        ha="right", va="bottom",
    )
    ax.plot([DEATH_IN_GAP], [y], marker="X", ms=5.4, color=INK, zorder=6)
    ax.plot([DECL], [y], marker="D", ms=3.8, mfc="white", mec=GREY, mew=0.9, zorder=6)
    ax.text(
        DEATH_IN_GAP - 0.1, y + 0.40, "death", fontsize=FS_ANNOT, color=INK,
        ha="center", va="bottom", weight="bold",
    )
    ax.text(
        DECL + 0.30,
        y,
        "\n".join(
            textwrap.wrap(
                "Never meets the definition, so recorded as an on-treatment death. "
                "Death competes with the exposure.",
                44,
            )
        ),
        fontsize=FS_ANNOT, color=PURP, ha="left", va="center", linespacing=1.35,
    )

    # -- patient 4 -----------------------------------------------------------
    y = ROW["p4"]
    rowlabel(y, "Patient 4\ncomparator, then\nlost to follow-up")
    oncare(y, 0, P4_STOP)
    absence(y, P4_STOP, PLANNED_END)
    lost(y, PLANNED_END, T_DETAIL)
    ax.plot([PLANNED_END], [y], marker="D", ms=3.8, color=INK, zorder=6)
    ax.plot([DECL], [y], marker="|", ms=10, mew=1.4, color=INK, zorder=6)
    ax.text(
        0.0,
        y - 0.72,
        "\n".join(
            textwrap.wrap(
                "In care on Patient 2's declaration date, so eligible as a comparator "
                "there, and stays in that arm after being lost.",
                104,
            )
        ),
        fontsize=FS_ANNOT, color=BLUE_INK, ha="left", va="top", linespacing=1.35,
    )

    # -- key -----------------------------------------------------------------
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=BLUE, edgecolor="none", label="in care"),
        Rectangle((0, 0), 1, 1, facecolor=VERM, edgecolor="none", label="lost to follow-up"),
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=GREY, lw=0.6, hatch="/////",
                  label="treatment missed, definition not yet met"),
        Line2D([], [], marker="D", ls="none", color=INK, ms=3.8, label="observed date"),
        Line2D([], [], marker="o", ls="none", mfc="white", mec=INK, mew=0.9, ms=3.8,
               label="inferred date"),
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.005),
        ncol=5,
        frameon=False,
        fontsize=FS_ANNOT,
        handlelength=1.2,
        handleheight=0.9,
        columnspacing=1.4,
        labelspacing=0.45,
        borderpad=0.0,
    )


# ================================================================ panel B ===
def draw_najera(ax):
    """Origin = treatment initiation, comparator = CURED only.

    The design missing from the previous version of this figure, and the one that
    produces the published 4-5 fold estimates. To be classified cured a patient
    must survive AND complete therapy, so the comparator contributes ~6 months in
    which it cannot die, while the exposed arm need only survive to its own
    default date. Comparator immortal time therefore EXCEEDS exposed immortal
    time and the net bias inflates -- which is the opposite of the intuition that
    immortal time is always null-ward.
    """
    curve_axis(ax)
    x = np.linspace(0, T_MAX, 400)
    ax.plot(x, cuminc(x, 0, 0.95, 0.13), color=VERM, lw=1.6, solid_capstyle="round")
    xc = np.linspace(PLANNED_END, T_MAX, 300)
    ax.plot(xc, cuminc(xc, PLANNED_END, 0.26, 0.12), color=BLUE, lw=1.6, solid_capstyle="round")
    ax.add_patch(
        Rectangle((0, 0), PLANNED_END, YMAX, facecolor=AMBER, alpha=0.20,
                  edgecolor="none", zorder=1)
    )
    ax.plot([0, PLANNED_END], [0, 0], color=BLUE, lw=3.0, solid_capstyle="butt", zorder=4)
    ax.annotate(
        "\n".join(
            textwrap.wrap(
                "cured patients must survive to month 6, so the COMPARATOR "
                "cannot die here",
                38,
            )
        ),
        xy=(PLANNED_END * 0.55, 0.03),
        xytext=(0.25, 1.02),
        fontsize=FS_ANNOT, color=AMBER_INK, ha="left", va="bottom", linespacing=1.35,
        arrowprops=dict(arrowstyle="-|>", lw=0.7, color=AMBER_INK, mutation_scale=5),
    )
    landmark_line(ax, 0, "origin: day 0", label_y=0.74)
    series_label(ax, 0.93, "ever lost to follow-up", VERM)
    series_label(ax, 0.21, "cured / completed therapy", BLUE)
    bias_note(
        ax,
        "Bias: comparator immortal time plus healthy-survivor selection. Both inflate.",
        INK,
    )


# ================================================================ panel C ===
def draw_completers(ax):
    curve_axis(ax)
    x = np.linspace(PLANNED_END, T_MAX, 300)
    ax.plot(x, cuminc(x, PLANNED_END, 0.92, 0.16), color=VERM, lw=1.6, solid_capstyle="round")
    ax.plot(x, cuminc(x, PLANNED_END, 0.30, 0.13), color=BLUE, lw=1.6, solid_capstyle="round")
    ax.add_patch(
        Rectangle((0, 0), PLANNED_END, YMAX, facecolor=AMBER, alpha=0.20,
                  edgecolor="none", zorder=1)
    )
    ax.text(
        PLANNED_END / 2,
        0.45,
        "\n".join(
            textwrap.wrap("both arms must reach month 6", 16)
        ),
        fontsize=FS_ANNOT, color=AMBER_INK, ha="center", va="center", linespacing=1.35, zorder=3,
    )
    landmark_line(ax, PLANNED_END, "landmark: month 6")
    series_label(ax, 0.90, "lost to follow-up", VERM)
    series_label(ax, 0.25, "completed therapy", BLUE)
    bias_note(
        ax,
        "Bias: healthy-survivor comparator. Inflates the contrast.",
        INK,
    )


# ================================================================ panel D ===
def draw_initiation(ax):
    curve_axis(ax)
    x = np.linspace(0, T_MAX, 400)
    ax.plot(x, cuminc(x, 0, 0.64, 0.14), color=BLUE, lw=1.6, solid_capstyle="round")
    xl = np.linspace(DECL, T_MAX, 300)
    ax.plot(xl, cuminc(xl, DECL, 0.50, 0.13), color=VERM, lw=1.6, solid_capstyle="round")
    ax.add_patch(
        Rectangle((0, 0), DECL, YMAX, facecolor=AMBER, alpha=0.20, edgecolor="none", zorder=1)
    )
    ax.plot([0, DECL], [0, 0], color=VERM, lw=3.0, solid_capstyle="butt", zorder=4)
    ax.annotate(
        "\n".join(
            textwrap.wrap(
                "the EXPOSED arm cannot die here: immortal time",
                38,
            )
        ),
        xy=(DECL * 0.5, 0.03),
        xytext=(0.25, 1.02),
        fontsize=FS_ANNOT, color=AMBER_INK, ha="left", va="bottom", linespacing=1.35,
        arrowprops=dict(arrowstyle="-|>", lw=0.7, color=AMBER_INK, mutation_scale=5),
    )
    landmark_line(ax, 0, "origin: day 0", label_y=0.74)
    series_label(ax, 0.67, "never lost to follow-up", BLUE)
    series_label(ax, 0.43, "ever lost to follow-up", VERM, va="top")
    bias_note(
        ax,
        "Bias: exposed-arm immortal time; comparator absorbs on-treatment deaths. "
        "Null-ward, so not what inflates the literature.",
        INK,
    )


# ================================================================ panel E ===
def draw_rolling(ax):
    curve_axis(ax)
    x = np.linspace(DECL, T_MAX, 300)
    ax.plot(x, cuminc(x, DECL, 0.92, 0.15), color=VERM, lw=1.6, solid_capstyle="round")
    ax.plot(x, cuminc(x, DECL, 0.42, 0.11), color=BLUE, lw=1.6, solid_capstyle="round")
    landmark_line(ax, DECL, "origin: each patient's own LTFU declaration date")
    series_label(ax, 0.90, "lost to follow-up", VERM)
    series_label(ax, 0.35, "in care at that same date", BLUE)
    bias_note(
        ax,
        "Residual bias: deaths before declaration cannot be classified as exposed.",
        PURP,
    )


# ================================================================ panel F ===
def draw_ccw(ax):
    curve_axis(ax, xlabel=True)
    xs = np.linspace(0, T_CLONE, 200)
    shared = cuminc(xs, 0, 0.55, 0.14)
    y_t = float(shared[-1])
    ax.plot(xs, shared, color=BLUE, lw=2.6, solid_capstyle="round", zorder=3)
    ax.plot(xs, shared, color=VERM, lw=1.3, ls=(0, (2.2, 1.9)), zorder=4)
    xa = np.linspace(T_CLONE, T_MAX, 300)
    ax.plot(xa, y_t + cuminc(xa, T_CLONE, 0.62, 0.17), color=VERM, lw=1.6, solid_capstyle="round")
    ax.plot(xa, y_t + cuminc(xa, T_CLONE, 0.34, 0.12), color=BLUE, lw=1.6, solid_capstyle="round")
    ax.add_patch(
        Rectangle((0, 0), T_CLONE, 0.55, facecolor=BLUE, alpha=0.13, edgecolor="none", zorder=1)
    )
    ax.plot([0], [0], marker="o", ms=4.2, color=INK, zorder=6)
    landmark_line(ax, 0, "origin: day 0, both clones")
    landmark_line(
        ax, T_CLONE, "strategies diverge:\nbecome LTFU by month T", top=0.99, label_y=0.56
    )
    ax.annotate(
        "\n".join(
            textwrap.wrap(
                "deaths here count in BOTH arms",
                40,
            )
        ),
        xy=(T_CLONE * 0.45, y_t * 0.45),
        xytext=(T_CLONE + 1.8, 0.015),
        fontsize=FS_ANNOT, color=BLUE_INK, ha="left", va="bottom", linespacing=1.35,
        arrowprops=dict(arrowstyle="-|>", lw=0.7, color=BLUE_INK, mutation_scale=5),
    )
    series_label(ax, 0.90, "became LTFU by T", VERM)
    series_label(ax, 0.42, "remained in care", BLUE, va="top")
    bias_note(
        ax,
        "No immortal time, no survivor conditioning. Ratios diluted: read the RD.",
        BLUE_INK,
    )


# ============================================================== assembly ====
# Panels B-D are the three prior-literature designs, each reproduced in this cohort
# (49_literature_design_reproduction.R). Keys are deliberately parallel -- same three
# fields, short values -- so the panels can be read across rather than one by one.
# Values are PHRASES. Anything that needs a sentence belongs in the caption, not
# here: five panels of prose in this column is what made the previous version
# unreadable. Fields are identical across panels so the column reads downward.
SPECS = [
    (
        "B",
        "Comparator = cured patients, followed from treatment initiation",
        None,
        draw_najera,
        [
            ("Origin", "Day 0"),
            ("Exposure", "Ever LTFU, from the treatment outcome"),
            ("Comparator", "Recorded as cured"),
        ],
    ),
    (
        "C",
        "Symmetric end-of-treatment landmark",
        None,
        draw_completers,
        [
            ("Origin", "Month 6; both arms must reach it"),
            ("Exposure", "Ever LTFU"),
            ("Comparator", "Recorded as cured"),
        ],
    ),
    (
        "D",
        "Comparator = everyone not lost to follow-up",
        None,
        draw_initiation,
        [
            ("Origin", "Day 0"),
            ("Exposure", "Ever LTFU"),
            ("Comparator", "All others, including on-treatment deaths"),
        ],
    ),
    (
        "E",
        "Rolling landmark at the LTFU declaration date",
        "TIMING ANALYSIS",
        draw_rolling,
        [
            ("Origin", "Each patient's own declaration date"),
            ("Exposure", "LTFU declared on that date"),
            ("Comparator", "In care that date; may be lost later"),
        ],
    ),
    (
        "F",
        "Clone, censor and weight from day 0",
        "PRIMARY ANALYSIS",
        draw_ccw,
        [
            ("Origin", "Day 0, both clones"),
            ("Exposure", "Strategy: become LTFU by month T"),
            ("Comparator", "Same patients, strategy: stay in care"),
        ],
    ),
]


def build(show_numbers: bool):
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    ax_a = fig.add_axes(A_BOX)
    panel_title(ax_a, "A", "Timeline of events and exposure assignment")
    draw_timeline(ax_a)

    for top, (letter, title, tag, drawer, rows) in zip(ROW_TOPS, SPECS):
        ax_p = fig.add_axes([A_BOX[0], top - 0.013 - PLOT_H, PLOT_W, PLOT_H])
        panel_title(ax_p, letter, title, tag)
        drawer(ax_p)
        ax_k = fig.add_axes([KEY_X, top - KEY_H, KEY_W, KEY_H])
        table = list(rows)
        if show_numbers:
            table.append(("Estimate", ESTIMATES[letter]))
        key_block(ax_k, table)

    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--numbers",
        action="store_true",
        help="add registry-sourced point estimates to the right-hand keys",
    )
    ap.add_argument("--stem", default=STEM)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig = build(args.numbers)
    stem = args.stem + ("_with_numbers" if args.numbers else "")
    written = []
    for ext, kw in (
        ("svg", {}),
        ("pdf", {}),
        ("png", {"dpi": 300}),
        ("tif", {"dpi": 300, "pil_kwargs": {"compression": "tiff_lzw"}}),
    ):
        path = OUTDIR / f"{stem}.{ext}"
        fig.savefig(path, facecolor="white", **kw)
        if ext == "tif":
            # PLOS / NAAS want RGB without an alpha channel
            from PIL import Image

            with Image.open(path) as im:
                im.convert("RGB").save(path, compression="tiff_lzw", dpi=(300, 300))
        written.append(path)
    plt.close(fig)
    for p in written:
        print(f"wrote {p.relative_to(ROOT)}  ({p.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
