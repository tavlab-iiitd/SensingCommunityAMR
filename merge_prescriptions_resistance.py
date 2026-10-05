"""Link monthly urine resistance to monthly prescriptions of the same antibiotic.

The linkage is at platform level (antibiotic and month), not patient level.
Combination products only match an identically composed laboratory agent.

Input : data/ast_clean.csv, data/prescriptions_clean.csv
Output: data/prescription_resistance_monthly.csv
"""

import numpy as np
import pandas as pd

SPELLING = {
    "amoxycillin": "amoxicillin",
    "clavulanate": "clavulanic acid",
    "tetractcline": "tetracycline",
    "ceoxitin": "cefoxitin",
    "cephotaxime": "cefotaxime",
    "cotrimoxazole": "sulfamethoxazole/trimethoprim",
}
COUNTS = ["SENSITIVE", "RESISTANT", "INTERMEDIATE"]


def join_key(name):
    """Lower-case, component-sorted drug name, e.g. 'sulfamethoxazole/trimethoprim'."""
    if pd.isna(name):
        return np.nan
    text = str(name).lower().strip().replace("`", "")
    text = text.replace("+", "/").replace("&", "/").replace("-", " ")
    for wrong, right in SPELLING.items():
        text = text.replace(wrong, right)
    parts = [" ".join(p.split()) for p in text.split("/") if p.split()]
    return "/".join(sorted(parts)) if parts else np.nan


# Monthly resistance per organism and antibiotic, urine cultures only
ast = pd.read_csv("data/ast_clean.csv", low_memory=False)
ast["organism_clean"] = ast["organism_clean"].astype(str).str.strip()
ast = ast[(ast["test_name"] == "URINE") & ~ast["organism_clean"].str.lower().isin(["no growth", "nan"])].copy()
ast["salt_name"] = ast["antibiotic_clean"].astype(str).str.strip().str.lower()
ast["join_key"] = ast["antibiotic_clean"].map(join_key)
ast = ast[ast["join_key"].notna()]

resistance = ast.groupby(["month_year", "join_key", "salt_name", "organism_clean", "aware_2025",
                          "who_class_2025"], dropna=False, as_index=False)[COUNTS].sum()
resistance["total_tests"] = resistance[COUNTS].sum(axis=1)
resistance["pct_resistant"] = (resistance["RESISTANT"] / resistance["total_tests"] * 100).round(2)

# Monthly prescriptions per product in the 25 states/UTs used for state-level analyses
rx = pd.read_csv("data/prescriptions_clean.csv", low_memory=False)
rx = rx[rx["state_included"]].copy()
rx["join_key"] = rx["product_clean"].map(join_key)

prescriptions = rx.groupby(["month_year", "join_key", "product_clean"], as_index=False).agg(
    total_prescriptions_count=("Rx_Count", "sum"))
antibiotic_totals = rx.groupby("month_year", as_index=False)["Rx_Count"].sum().rename(
    columns={"Rx_Count": "total_ab_prescriptions"})

# Keep antibiotics present in both datasets
merged = (resistance.merge(prescriptions, on=["month_year", "join_key"], how="left")
          .merge(antibiotic_totals, on="month_year", how="left"))
merged = merged[merged["join_key"].isin(set(prescriptions["join_key"]))].copy()
merged["pres_normalized"] = (merged["total_prescriptions_count"]
                             / merged["total_ab_prescriptions"] * 100).round(4)

merged = merged[["month_year", "salt_name", "product_clean", "join_key", "organism_clean",
                 "aware_2025", "who_class_2025", "pct_resistant", "RESISTANT", "SENSITIVE",
                 "INTERMEDIATE", "total_tests", "total_prescriptions_count",
                 "total_ab_prescriptions", "pres_normalized"]]
merged = merged.sort_values(["month_year", "salt_name", "organism_clean"])
merged.to_csv("data/prescription_resistance_monthly.csv", index=False)

print(f"{len(merged):,} rows; {merged['salt_name'].nunique()} antibiotics in both datasets")
