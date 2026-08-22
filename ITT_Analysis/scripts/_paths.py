"""Shared project path resolution for Python scripts in this repo.

The R side has had `_paths.R` for a while; the Python side grew a copy-pasted
`_find_project_root()` in ~16 scripts instead. New Python code should import
from here:

    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from _paths import PROJECT_ROOT, COHORT_CSV, ITT_RESULTS_DIR

Resolution order:
  1. TB_ABANDONMENT_ROOT environment variable
  2. Known Google Drive mounts for current collaborators

A candidate only counts if it contains ITT_Analysis/data. There is deliberately
NO repo-relative fallback: the old fallback resolved to the git checkout (which
holds no data) and turned a missing-mount problem into a confusing
file-not-found several hundred lines later.
"""
from __future__ import annotations

import os
from pathlib import Path

_CANDIDATES = [
    Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/TB SP 2026/LTFU Paper",
    Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
    Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper",
]


def _is_project_root(p: Path) -> bool:
    return (p / "ITT_Analysis" / "data").is_dir()


def find_project_root() -> Path:
    env = os.environ.get("TB_ABANDONMENT_ROOT")
    if env:
        p = Path(env).expanduser()
        if not _is_project_root(p):
            raise FileNotFoundError(
                f"TB_ABANDONMENT_ROOT is set to '{p}' but that is not a project root "
                f"(no ITT_Analysis/data inside it)."
            )
        return p.resolve()

    for c in _CANDIDATES:
        if _is_project_root(c):
            return c

    tried = "\n  - ".join(str(c) for c in _CANDIDATES)
    raise FileNotFoundError(
        "Could not locate the project root (the Google Drive folder holding Data/ "
        f"and ITT_Analysis/).\nTried:\n  - {tried}\n"
        "Mount Google Drive, or set TB_ABANDONMENT_ROOT to the folder containing "
        "ITT_Analysis/data."
    )


PROJECT_ROOT = find_project_root()
DATA_DIR = PROJECT_ROOT / "Data"
ITT_DATA_DIR = PROJECT_ROOT / "ITT_Analysis" / "data"
ITT_RESULTS_DIR = PROJECT_ROOT / "ITT_Analysis" / "results"
ITT_MI_DIR = ITT_DATA_DIR / "mi"

COHORT_CSV = ITT_DATA_DIR / "itt_cohort.csv"
FLOWCHART_CSV = DATA_DIR / "exclusion_flowchart.csv"
FINAL_TABLE_CSV = DATA_DIR / "Final_table_cleaned.csv"
