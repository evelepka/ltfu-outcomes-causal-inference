#!/usr/bin/env python3
"""Export the submission TIFFs (Plos Medicine/<round>/Figures/TIFF_<round>/FigN.tif).

WHY THIS EXISTS
---------------
The R2 TIFFs were exported by hand, so after a figure was regenerated nothing
said how to remake them. This reproduces that export: PLOS's maximum width of
7.5 in at 300 dpi (2,250 px), aspect ratio kept, RGB with no alpha, LZW.

Only figures whose source PNG is named in MAP are exported; any other TIFF in
the folder is left alone. An existing TIFF is moved to old/ before it is
replaced, never overwritten in place.

Usage:  python3 ITT_Analysis/scripts/78_export_plos_tiffs.py            # R2, all in MAP
        ROUND=R2 FIGS=4 python3 ITT_Analysis/scripts/78_export_plos_tiffs.py
"""
import os
import shutil
from datetime import date
from pathlib import Path

from PIL import Image

ROOT = Path(os.environ.get("TB_ABANDONMENT_ROOT", Path(__file__).resolve().parents[2]))
ROUND = os.environ.get("ROUND", "R2")
FIGDIR = ROOT / "Plos Medicine" / ROUND / "Figures"
OUT = FIGDIR / f"TIFF_{ROUND}"
WIDTH_PX, DPI = 2250, 300

# figure number -> source PNG, relative to FIGDIR unless absolute (script that draws it)
MAP = {
    2: ROOT / "ITT_Analysis" / "results" / "Figure_1_descriptive.png",  # 50a, then 50_make_fig1_descriptive.R
    4: "LTFU_primary_curves.png",     # 72_make_fig3_primary_curves.py
    5: "Figure_4_revised_5y.png",     # 73_make_fig4_revised.py
}

wanted = [int(x) for x in os.environ.get("FIGS", ",".join(map(str, MAP))).split(",")]
OUT.mkdir(parents=True, exist_ok=True)
for n in wanted:
    src = FIGDIR / MAP[n]
    dst = OUT / f"Fig{n}.tif"
    im = Image.open(src)
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, "white")
        im = im.convert("RGBA")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    h = round(im.height * WIDTH_PX / im.width)
    im = im.resize((WIDTH_PX, h), Image.LANCZOS)
    if dst.exists():
        old = OUT / "old"
        old.mkdir(exist_ok=True)
        shutil.move(str(dst), str(old / f"Fig{n}_{date.today():%Y_%m_%d}_replaced.tif"))
    im.save(dst, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    print(f"Fig{n}: {src.name} -> {dst.relative_to(ROOT)}  {im.size[0]}x{im.size[1]} px")
