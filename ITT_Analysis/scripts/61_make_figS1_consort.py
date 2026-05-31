#!/usr/bin/env python3
"""Figure S1 - CONSORT-style cohort flow diagram for the TB abandonment paper.

Design goals (per reviewer feedback):
  * comfortable, un-squeezed text (no cramped letters);
  * clear vertical separation between the flow boxes;
  * boxes sized snugly to their content (no large empty margins inside boxes).

Numbers are hard-coded and verified to add up at each step.
"""
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

# ---- content ---------------------------------------------------------------
# Main flow boxes: (bold title, N line)
flow = [
    ("TBweb notifications, 2013 to 2023", "N = 235,629"),
    ("First TB episode, eligible record", "N = 200,300"),
    ("Resident, adult, with consistent data", "N = 190,087"),
    ("Primary cohort: treatment within study\nperiod and with recorded start date",
     "N = 171,069"),
]

# Exclusion boxes sit in the gap between successive flow boxes.
excl = [
    ["Excluded:",
     "Subsequent (non-first) notifications, n = 16,347",
     "Change of diagnosis, n = 6,552",
     "Missing or incomplete outcome, n = 12,430"],
    ["Excluded:",
     "Inter-state transfer during treatment, n = 1,801",
     "Age < 15 years, n = 6,238",
     "Data inconsistencies, n = 111",
     "Pre-treatment deaths, n = 2,063"],
    ["Excluded:",
     "End dates outside 2013 to 2023, n = 17,624",
     "No recorded treatment start date, n = 1,394"],
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
out = "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper/ITT_Analysis/results/Figure_S1_consort.png"
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.10)
print("wrote", out)
