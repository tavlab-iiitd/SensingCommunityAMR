"""Clean the prescription data and classify every product with WHO AWaRe 2025.

Input : data/Prescription_Data_Pan_India_20260728.csv, data/aware-2025.xlsx
Output: data/prescriptions_clean.csv
"""

import pandas as pd

from name_cleaning import aware_group, load_who, product_label

NOT_ANTIBACTERIAL = {"quiniodochlor", "piperazine+ephedrine+guaifenesin+alcohol",
                     "sodium aminosalicylate"}
UNKNOWN_STATE = {"", "nan", "none", "null", "unknown", "not available", "others", "other"}
MIN_STATE_PRESCRIPTIONS = 1000

rx = pd.read_csv("data/Prescription_Data_Pan_India_20260728.csv", encoding="latin1", low_memory=False)
rx.columns = rx.columns.str.strip()

rx["date_of_rx"] = pd.to_datetime(rx["date_of_rx"], format="%Y-%m-%d", errors="coerce")
rx = rx[rx["date_of_rx"].between("2023-01-01", "2025-12-31")].copy()
rx["month_year"] = rx["date_of_rx"].dt.strftime("%Y-%m")
rx["Rx_Count"] = pd.to_numeric(rx["Rx_Count"]).astype(int)
rx["patient_rx_pincode"] = pd.to_numeric(rx["patient_rx_pincode"], errors="coerce")
for column in ["salt_name", "district", "patient_rx_state"]:
    rx[column] = rx[column].astype(str).str.strip()

salt = rx["salt_name"].str.lower().str.replace(r"\s*\+\s*", "+", regex=True)
rx = rx[~salt.isin(NOT_ANTIBACTERIAL)].copy()

# AWaRe 2025: the whole product must match a WHO entry, otherwise it is Unclassified
who = load_who()
products = rx["salt_name"].unique()
groups = {p: aware_group(p, who) for p in products}
rx["aware_2025"] = rx["salt_name"].map(lambda p: groups[p][0])
rx["who_class_2025"] = rx["salt_name"].map(lambda p: groups[p][1])
rx["product_clean"] = rx["salt_name"].map({p: product_label(p, who) for p in products})

# State-level analyses use states/UTs with at least 1,000 antibiotic prescriptions
state = rx["patient_rx_state"]
rx["patient_rx_state"] = state.where(~state.str.lower().isin(UNKNOWN_STATE), "Not available")
state_totals = rx[rx["patient_rx_state"] != "Not available"].groupby("patient_rx_state")["Rx_Count"].sum()
rx["state_included"] = rx["patient_rx_state"].isin(state_totals[state_totals >= MIN_STATE_PRESCRIPTIONS].index)

columns = ["month_year", "salt_name", "product_clean", "aware_2025", "who_class_2025",
           "patient_rx_state", "district", "patient_rx_pincode", "Rx_Count", "state_included"]
rx[columns].to_csv("data/prescriptions_clean.csv", index=False)

print(f"{len(rx):,} rows, {rx['Rx_Count'].sum():,} antibiotic prescriptions")
print(rx.groupby("aware_2025")["Rx_Count"].sum())
print(f"{rx.loc[rx['state_included'], 'patient_rx_state'].nunique()} states/UTs with >= 1,000 prescriptions")
