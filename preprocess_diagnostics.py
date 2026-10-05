"""Clean the antimicrobial susceptibility testing (AST) data.

Standardises specimen, organism and antibiotic names and links each antibiotic to
WHO AWaRe 2025.

Input : data/AMR_Pan_India_Diagnostic_Antibiotic.csv, data/aware-2025.xlsx
Output: data/ast_clean.csv
"""

import pandas as pd

from name_cleaning import clean_antibiotic, clean_organism, lab_aware_group, load_who

SPECIMEN = {
    "URINE CULTURE AND SENSITIVITY": "URINE",
    "SPUTUM CULTURE AND SENSITIVITY (AEROBIC)": "SPUTUM",
    "CULTURE AEROBIC BLOOD (AUTOMATED), ADULT": "BLOOD",
    "CULTURE AEROBIC BLOOD (AUTOMATED), PEDIATRIC": "BLOOD",
}
NOT_ANTIBACTERIAL = {"Amphotericin B", "Caspofungin", "Fluconazole", "Flucytosine",
                     "Micafungin", "Voriconazole", "Optochin"}

ast = pd.read_csv("data/AMR_Pan_India_Diagnostic_Antibiotic.csv", low_memory=False)
ast["end_state_date"] = pd.to_datetime(ast["end_state_date"], format="%Y-%m-%d")
ast = ast[ast["end_state_date"].between("2023-01-01", "2025-12-31")].copy()
ast["month_year"] = ast["end_state_date"].dt.strftime("%Y-%m")
ast["end_state_date"] = ast["end_state_date"].dt.strftime("%Y-%m-%d")
ast["test_name"] = ast["test_name"].replace(SPECIMEN)

organisms = {name: clean_organism(name) for name in ast["organism_isolated"].unique()}
antibiotics = {name: clean_antibiotic(name) for name in ast["antibiotic"].unique()}
ast["organism_clean"] = ast["organism_isolated"].map(organisms)
ast["antibiotic_clean"] = ast["antibiotic"].map(antibiotics)
ast = ast[~ast["antibiotic_clean"].isin(NOT_ANTIBACTERIAL)].copy()

who = load_who()
groups = {a: lab_aware_group(a, who) for a in ast["antibiotic_clean"].dropna().unique()}
ast["aware_2025"] = ast["antibiotic_clean"].map(lambda a: groups.get(a, (None, None))[0])
ast["who_class_2025"] = ast["antibiotic_clean"].map(lambda a: groups.get(a, (None, None))[1])

columns = ["end_state_date", "month_year", "order_state", "Order_Pincode", "test_name",
           "organism_clean", "antibiotic_clean", "aware_2025", "who_class_2025",
           "RESISTANT", "INTERMEDIATE", "SENSITIVE"]
ast[columns].to_csv("data/ast_clean.csv", index=False)

print(f"{len(ast):,} AST rows, {ast[['RESISTANT', 'INTERMEDIATE', 'SENSITIVE']].sum().sum():,} interpretable results")
