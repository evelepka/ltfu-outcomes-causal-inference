#!/usr/bin/env python3
"""Flag deaths with tuberculosis anywhere on the death certificate.

The analysis dataset carries one cause field, `cause_of_death_code`, which is
verified here to be identical to SIM's CAUSABAS -- the UNDERLYING cause. The
Methods say so and call the tuberculosis estimate conservative for that reason.

The full certificate is nevertheless available: `Banco de dados/LINKAGE SIM (1).xlsx`,
sheet "Limpo", carries LINHAA-LINHAD and LINHAII alongside CAUSABAS. Tuberculosis
recorded as a contributing condition is therefore recoverable, and it is not
rare: it adds roughly a quarter again on top of the underlying-cause count, and
it sits mostly under HIV, COPD and lung cancer -- precisely the deaths the
cause-specific analysis has to hold out because they might be tuberculosis.

This writes the flag; nothing consumes it unless TB_ANY_LINE=1 is set.

Output: Data/tb_any_line_flag.csv  (sinan_clean, tb_anyline, tb_underlying)
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT, DATA_DIR                    # noqa: E402

SRC = Path(PROJECT_ROOT) / "Banco de dados" / "LINKAGE SIM (1).xlsx"
LINES = ["LINHAA", "LINHAB", "LINHAC", "LINHAD", "LINHAII"]
# same codes the main analysis calls tuberculosis
TB = re.compile(r"\b(A1[5-9]\d?|B90\d?|B200)")


def has_tb(s):
    return (s.fillna("").astype(str).str.upper()
             .str.replace(r"[^A-Z0-9]", " ", regex=True).str.contains(TB))


def main() -> int:
    if not SRC.exists():
        sys.exit(f"missing {SRC}")
    d = pd.read_excel(SRC, sheet_name="Limpo",
                      usecols=["SINAN", "CAUSABAS"] + LINES, dtype=str)
    key = (d.SINAN.fillna("").astype(str).str.strip()
             .str.replace(r"\.0$", "", regex=True))

    underlying = has_tb(d.CAUSABAS)
    anyline = underlying.copy()
    for c in LINES:
        anyline |= has_tb(d[c])

    out = (pd.DataFrame({"sinan_clean": key,
                         "tb_underlying": underlying.astype(int),
                         "tb_anyline": anyline.astype(int)})
             .query("sinan_clean != ''")
             .groupby("sinan_clean", as_index=False).max())

    # the join key and the cause field are both verifiable against the analysis
    # extract; if either drifts this is the place it will show
    raw = pd.read_csv(Path(DATA_DIR) / "Final_table_cleaned.csv", dtype=str,
                      usecols=["sinan_clean", "cause_of_death_code"], low_memory=False)
    raw = raw[raw.cause_of_death_code.notna() & (raw.cause_of_death_code.str.strip() != "")]
    chk = raw.merge(d.assign(k=key)[["k", "CAUSABAS"]].drop_duplicates("k"),
                    left_on=raw.sinan_clean.str.strip(), right_on="k", how="inner")
    agree = (chk.cause_of_death_code.str.upper().str.strip()
             == chk.CAUSABAS.str.upper().str.strip()).mean()
    if agree < 0.99:
        sys.exit(f"CAUSABAS matches cause_of_death_code in only {agree:.1%} "
                 f"of {len(chk):,} joined deaths -- the key or the field has changed")
    print(f"  key check: CAUSABAS == cause_of_death_code in {agree:.1%} of {len(chk):,} deaths")

    dst = Path(DATA_DIR) / "tb_any_line_flag.csv"
    out.to_csv(dst, index=False)
    n_u, n_a = int(out.tb_underlying.sum()), int(out.tb_anyline.sum())
    print(f"wrote {dst}  ({len(out):,} people)")
    print(f"  tuberculosis as underlying cause : {n_u:,}")
    print(f"  tuberculosis anywhere            : {n_a:,}   (+{n_a - n_u:,}, "
          f"{100 * (n_a - n_u) / max(n_u, 1):.0f}% more)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
