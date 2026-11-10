"""National month–drug linkage used for the exploratory lead–lag analysis."""

import pandas as pd

from linkage_names import join_key

COUNTS = ["SENSITIVE", "RESISTANT", "INTERMEDIATE"]

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


rx = pd.read_csv("data/prescriptions_clean.csv", low_memory=False)
rx = rx[rx["state_included"]].copy()
rx["join_key"] = rx["product_clean"].map(join_key)

prescriptions = rx.groupby(["month_year", "join_key", "product_clean"], as_index=False).agg(
    total_prescriptions_count=("Rx_Count", "sum"))
antibiotic_totals = rx.groupby("month_year", as_index=False)["Rx_Count"].sum().rename(
    columns={"Rx_Count": "total_ab_prescriptions"})


merged = (resistance.merge(prescriptions, on=["month_year", "join_key"], how="left", validate="many_to_one")
          .merge(antibiotic_totals, on="month_year", how="left", validate="many_to_one"))
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
