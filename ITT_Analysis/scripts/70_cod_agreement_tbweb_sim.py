"""
70_cod_agreement_tbweb_sim.py
-----------------------------
Agreement between TBweb treatment-outcome death classification (Obito TB vs
Obito NTB) and SIM underlying-cause-of-death ICD-10 coding, for the cohort.
Addresses reviewer comment #1 (check agreement in causes of death between
TBweb and SIM; compare with our hybrid definition).

Outputs:
  ITT_Analysis/results/cod_agreement_tbweb_sim.csv   (cross-tab + metrics)
Prints a summary suitable for the appendix.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd


def root():
    for c in [Path.home() / "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper",
              Path.home() / "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper"]:
        if c.exists():
            return c
    raise SystemExit("GDrive root not found")


BASE = root()
raw = pd.read_csv(BASE / "Data" / "Final_table_cleaned.csv", low_memory=False)
cohort = pd.read_csv(BASE / "ITT_Analysis" / "data" / "itt_cohort.csv", low_memory=False)
cohort_ids = set(cohort["sinan_clean"])

TRANSFER = {"Transf Outro Municipio", "Transf Outro Estado/Pais"}

# ---- SIM cause classification (mirror 30i classify_cod, ICD-only) ----------
def sim_class(cod):
    if pd.isna(cod) or not str(cod).strip():
        return np.nan  # no SIM cause
    c = str(cod).upper().strip()
    if c.startswith(("A15", "A16", "A17", "A18", "A19", "B90")) or c.startswith("B200"):
        return "TB"
    return "non-TB"


# Death record per individual: latest dod with a cause code
raw["dod_parsed"] = pd.to_datetime(raw["dod"], errors="coerce")
deathrec = raw[raw["dod_parsed"].notna()].sort_values("dod_parsed")
sim = (deathrec.dropna(subset=["cause_of_death_code"])
       .groupby("sinan_clean")["cause_of_death_code"].last().reset_index())
sim["sim_cause"] = sim["cause_of_death_code"].apply(sim_class)

# TBweb death classification: Obito TB / Obito NTB among the first Novo episode
novo = raw[raw["case_type"].astype(str).str.strip().str.title() == "Novo"]
novo = novo[novo["case_outcome"].notna() & (novo["case_outcome"].astype(str).str.strip() != "")
            & (novo["case_outcome"] != "Mud Diag") & (~novo["case_outcome"].isin(TRANSFER))]
novo = novo.sort_values("end_date").groupby("sinan_clean", as_index=False).first()
novo["tbweb_cause"] = novo["case_outcome"].map({"Obito TB": "TB", "Obito NTB": "non-TB"})
tbweb = novo[["sinan_clean", "case_outcome", "tbweb_cause"]]

m = tbweb.merge(sim[["sinan_clean", "sim_cause"]], on="sinan_clean", how="outer")
m = m[m["sinan_clean"].isin(cohort_ids)]

# Coverage
n_tbweb_death = m["tbweb_cause"].notna().sum()
n_sim_cause = m["sim_cause"].notna().sum()
both = m[m["tbweb_cause"].notna() & m["sim_cause"].notna()]
print(f"Cohort N = {len(cohort_ids):,}")
print(f"Individuals with a TBweb death classification (Obito TB/NTB): {n_tbweb_death:,}")
print(f"Individuals with a SIM ICD cause of death:                    {n_sim_cause:,}")
print(f"Both sources available (agreement denominator):              {len(both):,}")

# Cross-tab among both available
ct = pd.crosstab(both["tbweb_cause"], both["sim_cause"], rownames=["TBweb"], colnames=["SIM"])
print("\nCross-tab (TBweb rows x SIM cols), individuals with both:")
print(ct)

# Agreement metrics (treat as binary TB vs non-TB)
a = both.copy()
a["agree"] = a["tbweb_cause"] == a["sim_cause"]
po = a["agree"].mean()
# Cohen's kappa
p_tb_tbweb = (a["tbweb_cause"] == "TB").mean()
p_tb_sim = (a["sim_cause"] == "TB").mean()
pe = p_tb_tbweb * p_tb_sim + (1 - p_tb_tbweb) * (1 - p_tb_sim)
kappa = (po - pe) / (1 - pe)
# Sensitivity/specificity of TBweb vs SIM (SIM as reference)
tp = ((a["tbweb_cause"] == "TB") & (a["sim_cause"] == "TB")).sum()
fn = ((a["tbweb_cause"] == "non-TB") & (a["sim_cause"] == "TB")).sum()
fp = ((a["tbweb_cause"] == "TB") & (a["sim_cause"] == "non-TB")).sum()
tn = ((a["tbweb_cause"] == "non-TB") & (a["sim_cause"] == "non-TB")).sum()
sens = tp / (tp + fn) if (tp + fn) else np.nan
spec = tn / (tn + fp) if (tn + fp) else np.nan
print(f"\nPercent agreement: {po*100:.1f}%")
print(f"Cohen's kappa:     {kappa:.3f}")
print(f"TBweb vs SIM (SIM ref): sensitivity for TB = {sens*100:.1f}%, specificity = {spec*100:.1f}%")

# Save
out = BASE / "ITT_Analysis" / "results" / "cod_agreement_tbweb_sim.csv"
with open(out, "w") as f:
    f.write(f"cohort_N,{len(cohort_ids)}\n")
    f.write(f"n_tbweb_death_class,{n_tbweb_death}\n")
    f.write(f"n_sim_cause,{n_sim_cause}\n")
    f.write(f"n_both,{len(both)}\n")
    f.write(f"percent_agreement,{po}\n")
    f.write(f"cohen_kappa,{kappa}\n")
    f.write(f"tbweb_sens_for_TB,{sens}\n")
    f.write(f"tbweb_spec,{spec}\n\n")
    f.write("crosstab_TBweb_rows_SIM_cols\n")
    ct.to_csv(f)
print(f"\nWrote {out}")
