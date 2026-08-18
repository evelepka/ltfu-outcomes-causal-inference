#!/usr/bin/env python3
"""Refuse to publish patient data. Run before every commit/push to this repo.

Individual-level TBweb/SIM records are shared via Google Drive ONLY. They carry
sinan_clean (SINAN notification ID) with dob, dod, tx_city and address_type
alongside HIV status, homelessness and incarceration history -- directly
identifying, and covered by an ethics approval for secondary use that does not
extend to publication. GitHub retains deleted content in caches and forks, so a
mistake here is not reversible.

Exit 0 = safe to push. Exit 1 = STOP.
"""
import csv, os, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

IDENT = {"sinan_clean", "sinan_padded", "sinan_original", "dob", "dod",
         "death_date", "tx_city", "address_type", "notification_date",
         "best_start", "end_date"}
MAX_CSV_ROWS = 1000          # aggregate result tables are tens of rows
MAX_NONFIG_BYTES = 2_000_000


# ---------------------------------------------------------------------------
# The CLAUDE.md in THIS repo is the public one: repo layout and the
# data-never-here rule. The project's working CLAUDE.md lives in Google Drive and
# carries the journal submission ID, reviewer state and unpublished estimates.
# Copying the private one over this one would publish all of that, so fail loudly.
# ---------------------------------------------------------------------------
PRIVATE_MARKERS = [
    (r"PMEDICINE-D-\d", "journal submission identifier"),
    (r"\bReviewer\s*[#0-9]", "reviewer references (peer review is confidential)"),
    (r"\baHR\s*\d\.\d", "adjusted hazard ratio estimate"),
    (r"\blate aHR\b", "unpublished effect estimate"),
    (r"Plos Medicine/R1", "revision working directory"),
]


def check_public_claude_md():
    p = ROOT / "CLAUDE.md"
    if not p.exists():
        return []
    t = p.read_text(errors="replace")
    out = []
    for pat, why in PRIVATE_MARKERS:
        m = re.search(pat, t)
        if m:
            out.append(f"CLAUDE.md contains {why} ({m.group(0)!r}). This file is "
                       f"PUBLIC. The working CLAUDE.md belongs in Google Drive only.")
    return out


def main():
    files = subprocess.run(["git", "ls-files"], capture_output=True,
                           text=True).stdout.split()
    fails = []
    for f in files:
        if not os.path.isfile(f):
            continue
        if f.lower().endswith((".docx", ".doc", ".gdoc")):
            fails.append(f"manuscript/report binary staged: {f}")
        if (os.path.getsize(f) > MAX_NONFIG_BYTES
                and not f.startswith("figures/")):
            fails.append(f"oversized non-figure ({os.path.getsize(f)/1e6:.1f} MB): {f}")
        if not f.lower().endswith(".csv"):
            continue
        try:
            with open(f, newline="", errors="replace") as fh:
                hdr = next(csv.reader(fh))
                rows = sum(1 for _ in fh)
        except Exception as e:                                  # noqa: BLE001
            fails.append(f"unreadable CSV {f}: {e}")
            continue
        bad = {c.strip().strip('"').lower() for c in hdr} & IDENT
        if bad:
            fails.append(f"IDENTIFIER COLUMN in {f}: {sorted(bad)} ({rows:,} rows)")
        if rows > MAX_CSV_ROWS:
            fails.append(f"{f} has {rows:,} rows (>{MAX_CSV_ROWS}); "
                         f"aggregate tables only")
    fails.extend(check_public_claude_md())
    print(f"check_no_patient_data: {len(files):,} tracked files")
    for x in fails:
        print(f"  FAIL  {x}")
    print("SAFE TO PUSH" if not fails else f"DO NOT PUSH ({len(fails)} problem(s))")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main())
