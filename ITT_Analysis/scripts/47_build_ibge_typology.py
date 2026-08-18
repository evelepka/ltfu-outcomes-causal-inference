#!/usr/bin/env python3
"""47. Municipality -> IBGE rural/urban typology crosswalk (Reviewer 1 comment 11).

SOURCE
------
IBGE, "Classificacao e caracterizacao dos espacos rurais e urbanos do Brasil:
uma primeira aproximacao" (2017). Table `Tipologia_municipal_rural_urbano.xlsx`:
  geoftp.ibge.gov.br/organizacao_do_territorio/tipologias_do_territorio/
  classificacao_rural_e_urbana/
  classificacao_e_caracterizacao_dos_espacos_rurais_e_urbanos_do_brasil_2017/tabelas/
Five national levels: Urbano, IntermediarioAdjacente, IntermediarioRemoto,
RuralAdjacente, RuralRemoto. Sao Paulo State has only THREE (no Remoto).

WHAT THIS WRITES
----------------
`ITT_Analysis/external/municipality_typology_sp.csv` -- a MUNICIPALITY-level
crosswalk only. Deliberately no `sinan_clean`: this file carries no patient data,
so it is safe for the public code mirror.

Keyed on the RAW TBweb `tx_city` string, verbatim, for every distinct value present
in `Data/Final_table_cleaned.csv`. That is deliberate: accent/case normalisation is
locale-dependent and easy to get subtly wrong in R, so it is done ONCE here in
Python and `_rolling.R` only ever does an exact string match. The script asserts
that every distinct `tx_city` in the data is mapped, so a silent coverage loss is
impossible.

Columns:
  municipality  RAW TBweb `tx_city` value, exactly as stored
  ibge_tipo     the IBGE level, verbatim ("NotAMunicipality" for custody)
  geo_class     descriptive analysis variable, Sao Paulo city split out
  geo4          the MODEL variable: geo_class with custody folded into the
                reference level, because `geo_class == "Prison"` is 100% collinear
                with the `incarcerated` covariate (verified 2026-08-18: all 20,669
                are incarcerated; 99.97% of incarcerated are custody, 7 exceptions
                of 20,676). Including it as its own level triggers "Loglik
                converged before variable 14; coefficient may be infinite" and
                contaminates the exposure estimate.
  cd_gcmun      IBGE municipality geocode

WHY `geo_class` IS NOT JUST `ibge_tipo`
---------------------------------------
Measured 2026-08-18: 97.0% of matched patients live in `Urbano` municipalities, and
crude mortality is nearly flat across IBGE levels (Urbano 11.3%, RuralAdjacente
11.4%, IntermediarioAdjacente 10.0%). So `ibge_tipo` alone has almost no
discriminating power here and cannot address the confounding that actually exists,
which is WITHIN `Urbano`: Sao Paulo city holds 33.5% of the cohort but 41.9% of the
exposed, and the n-weighted ecological correlation between a municipality's LTFU
rate and its mortality rate is +0.61.

`geo_class` therefore splits Sao Paulo city out as its own level and uses the IBGE
typology for everything else. Both parts are externally defined -- no cutpoint is
chosen from our own outcome data.

NAME RECONCILIATION
-------------------
TBweb and IBGE disagree on three spellings; all three were checked by hand against
the IBGE SP list (645 municipalities) and are unambiguous:
  MOGI-MIRIM        -> Moji Mirim      (IBGE keeps the archaic 'j')
  FLORINEA          -> Florinia        (IBGE spelling differs)
  SAO JOSE DO RIO*  -> Urbano          (TBweb truncation; BOTH candidates,
                                        Rio Preto and Rio Pardo, are Urbano,
                                        so the typology is unambiguous)
`DETENTO` is not a municipality -- TBweb records it in `tx_city` for people
treated in custody (~12% of the cohort). It gets its own `geo_class` level; note
`incarcerated` is already a model covariate, so this is belt-and-braces.

Usage:  python3 ITT_Analysis/scripts/47_build_ibge_typology.py
"""
import sys
import unicodedata
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "ITT_Analysis" / "external"
XLSX = EXT / "Tipologia_municipal_rural_urbano.xlsx"
OUT = EXT / "municipality_typology_sp.csv"
URL = ("https://geoftp.ibge.gov.br/organizacao_do_territorio/tipologias_do_territorio/"
       "classificacao_rural_e_urbana/"
       "classificacao_e_caracterizacao_dos_espacos_rurais_e_urbanos_do_brasil_2017/"
       "tabelas/Tipologia_municipal_rural_urbano.xlsx")

ALIASES = {"MOGI-MIRIM": "MOJI MIRIM", "MOGI MIRIM": "MOJI MIRIM",
           "FLORINEA": "FLORINIA"}
TRUNCATED = {"SAO JOSE DO RIO*": "Urbano"}   # both candidates are Urbano
SAO_PAULO = "SAO PAULO"
PRISON = "DETENTO"
GEO_REF = "Urbano"        # reference level for geo4; custody folds in here
RAW_CSV = ROOT / "Data" / "Final_table_cleaned.csv"


def norm(s):
    if pd.isna(s):
        return None
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
    return " ".join(s.upper().replace("'", " ").replace("`", " ").split())


def main():
    EXT.mkdir(parents=True, exist_ok=True)
    if not XLSX.exists():
        print(f"[47] downloading {URL}")
        urllib.request.urlretrieve(URL, XLSX)
    ib = pd.read_excel(XLSX, sheet_name="Tipologia_munic_rural_urbano")
    sp = ib[ib.SIG_UF == "SP"].copy()
    sp["tx_city_key"] = sp.NM_MUN.map(norm)
    print(f"[47] IBGE Sao Paulo municipalities: {len(sp)}")
    print(sp.TIPO.value_counts().to_string())

    rows = []
    for _, r in sp.iterrows():
        key = r.tx_city_key
        geo = "Sao Paulo city" if key == SAO_PAULO else r.TIPO
        rows.append({"tx_city_key": key, "ibge_tipo": r.TIPO,
                     "geo_class": geo, "cd_gcmun": r.CD_GCMUN})
    for alias, target in ALIASES.items():
        t = sp[sp.tx_city_key == target]
        if t.empty:
            print(f"[47] WARNING alias target missing: {target}")
            continue
        r = t.iloc[0]
        rows.append({"tx_city_key": alias, "ibge_tipo": r.TIPO,
                     "geo_class": r.TIPO, "cd_gcmun": r.CD_GCMUN})
    for key, tipo in TRUNCATED.items():
        rows.append({"tx_city_key": key, "ibge_tipo": tipo,
                     "geo_class": tipo, "cd_gcmun": pd.NA})
    rows.append({"tx_city_key": PRISON, "ibge_tipo": "NotAMunicipality",
                 "geo_class": "Prison", "cd_gcmun": pd.NA})

    xw = pd.DataFrame(rows).drop_duplicates("tx_city_key")

    # --- key on the RAW TBweb string, for every value present in the data ------
    raw = pd.read_csv(RAW_CSV, low_memory=False, usecols=["tx_city"])
    distinct = raw.tx_city.dropna().drop_duplicates().to_frame()
    distinct["tx_city_key"] = distinct.tx_city.map(norm)
    out = distinct.merge(xw, on="tx_city_key", how="left")

    unmapped = out[out.geo_class.isna()].tx_city.tolist()
    if unmapped:
        raise SystemExit(f"[47] {len(unmapped)} tx_city values unmapped: "
                         f"{unmapped[:10]} -- add an alias in ALIASES/TRUNCATED")
    out["geo4"] = out.geo_class.where(out.geo_class != "Prison", GEO_REF)
    out = out.rename(columns={"tx_city": "municipality"})
    out = out[["municipality", "ibge_tipo", "geo_class", "geo4", "cd_gcmun"]]
    assert "sinan_clean" not in out.columns, "crosswalk must carry no patient data"
    assert out.municipality.is_unique, "one row per raw municipality string"
    out.to_csv(OUT, index=False)
    print(f"\n[47] {len(distinct)} distinct tx_city values in the data, all mapped")
    print(f"[47] wrote {len(out)} rows -> {OUT.relative_to(ROOT)}")
    print("\ngeo_class (descriptive):"); print(out.geo_class.value_counts().to_string())
    print("\ngeo4 (model variable, custody folded):")
    print(out.geo4.value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
