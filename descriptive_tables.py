"""Table 1 and Supplementary Tables 1 and 3 (descriptive summaries).

Input : data/prescriptions_clean.csv, data/State_Month_Prescription_updated.csv,
        data/AMR_Pan_India_Aggregated_Diagnostic.csv, data/ast_clean.csv
Output: results/tables/table1.csv, supp_table1_diagnostic_dataset.csv,
        supp_table3_prescription_dataset.csv, supp_table3_unclassified_products.csv
"""

import os

import pandas as pd

from name_cleaning import clean_organism

os.makedirs("results/tables", exist_ok=True)
GROUPS = ["Access", "Watch", "Reserve", "Not recommended", "Unclassified"]


def n(x):
    return f"{int(x):,}"


def pct(x, total):
    return f"{100 * x / total:.2f}"


def block(title, counts, total):
    """Table rows (title, name, count, %) for a series of counts."""
    return [[title if i == 0 else "", name, count, pct(count, total)]
            for i, (name, count) in enumerate(counts.items())]


# ---------------------------------------------------------------- prescriptions
rx = pd.read_csv("data/prescriptions_clean.csv", low_memory=False)
medicines = pd.read_csv("data/State_Month_Prescription_updated.csv")
medicines = medicines[medicines["month_year"].between("2023-01", "2025-12")]

total_rx = rx["Rx_Count"].sum()
total_medicines = medicines["presc"].sum()
monthly_share = 100 * rx.groupby("month_year")["Rx_Count"].sum() / medicines.groupby("month_year")["presc"].sum()

located = rx[rx["patient_rx_state"] != "Not available"].copy()
located["state"] = located["patient_rx_state"].str.title()
state_totals = located.groupby("state")["Rx_Count"].sum().sort_values(ascending=False)
included = rx.loc[rx["state_included"], "Rx_Count"].sum()
small_states = located[~located["state_included"]].groupby("state")["Rx_Count"].sum().sort_values(ascending=False)
no_state = rx.loc[rx["patient_rx_state"] == "Not available", "Rx_Count"].sum()

state_months = located.groupby(["state", "month_year"])["Rx_Count"].sum()   # active state-months only
possible_state_months = located["state"].nunique() * 36
q1, median, q3 = state_months.quantile([0.25, 0.5, 0.75])

aware = rx.groupby("aware_2025")["Rx_Count"].sum()
classes = rx.groupby("who_class_2025", dropna=False)["Rx_Count"].sum().sort_values(ascending=False)
classes.index = classes.index.fillna("Unclassified")

# ---------------------------------------------------------------- diagnostics
dx = pd.read_csv("data/AMR_Pan_India_Aggregated_Diagnostic.csv", encoding="ISO-8859-1", low_memory=False)
dx["end_state_date"] = pd.to_datetime(dx["end_state_date"], format="%Y-%m-%d")
dx = dx[dx["end_state_date"].between("2023-01-01", "2025-12-31")].copy()
dx["month_year"] = dx["end_state_date"].dt.strftime("%Y-%m")
dx["specimen"] = dx["test_name"].map({
    "URINE CULTURE AND SENSITIVITY": "Urine",
    "SPUTUM CULTURE AND SENSITIVITY (AEROBIC)": "Sputum",
    "CULTURE AEROBIC BLOOD (AUTOMATED), ADULT": "Blood",
    "CULTURE AEROBIC BLOOD (AUTOMATED), PEDIATRIC": "Blood",
})
dx["organism"] = dx["organism_isolated"].map({o: clean_organism(o) for o in dx["organism_isolated"].unique()})
# Culture positive = any culture not reported as "No growth" (includes 67 reports with no readable organism name)
dx["positive"] = ~(dx["organism"] == "No growth")

samples = dx["Total_Test"].sum()
positives = dx.loc[dx["positive"], "Total_Test"].sum()
urine = dx.loc[dx["specimen"] == "Urine", "Total_Test"].sum()
organisms = dx[dx["positive"]].groupby("organism")["Total_Test"].sum().sort_values(ascending=False)
dx_possible = dx["order_state"].nunique() * dx["month_year"].nunique()
dx_active = len(dx[["order_state", "month_year"]].drop_duplicates())

ast = pd.read_csv("data/ast_clean.csv", low_memory=False)
ast = ast[ast["test_name"] == "URINE"].copy()
ast["tested"] = ast[["RESISTANT", "INTERMEDIATE", "SENSITIVE"]].sum(axis=1)
urine_resistance = ast.groupby(["organism_clean", "antibiotic_clean"])[["RESISTANT", "tested"]].sum()

# ---------------------------------------------------------------- Table 1
table1 = [
    ["Dataset size", "Total medicines recorded in all prescriptions", n(total_medicines), "—"],
    ["", "Antibiotic prescriptions; in states/UTs with >=1,000", f"{n(total_rx)}; {n(included)}", "—"],
    ["", "Antibiotics as a share of all prescribed medicines", pct(total_rx, total_medicines) + "%", "—"],
    ["", "Monthly antibiotic share (range)", f"{monthly_share.min():.2f}%-{monthly_share.max():.2f}%", "—"],
    ["Geographic coverage", "States/UTs observed; state-specific eligible; below threshold",
     f"{located['state'].nunique()}; {located.loc[located['state_included'], 'state'].nunique()}; "
     f"{len(small_states)} ({n(small_states.sum())} prescriptions)", str(dx["order_state"].nunique())],
    ["", "Districts represented", str(rx["district"].nunique()), str(dx["district"].nunique())],
    ["", "Pincodes represented", n(rx["patient_rx_pincode"].nunique()), n(dx["Order_Pincode"].nunique())],
    ["Temporal coverage", "Study period", "2023-2025", "2023-2025"],
    ["", "Zero-activity state/UT-months",
     f"{possible_state_months - len(state_months)}/{n(possible_state_months)}; "
     f"{pct(possible_state_months - len(state_months), possible_state_months)}%",
     f"{dx_possible - dx_active}/{dx_possible}"],
    ["", "Active state/UT-months",
     f"{n(len(state_months))}/{n(possible_state_months)}; {pct(len(state_months), possible_state_months)}%",
     f"{dx_active}/{dx_possible}; {pct(dx_active, dx_possible)}%"],
    ["State-month intensity", "Median prescriptions per active state/UT-month (IQR)",
     f"{n(median)} ({n(q1)}-{n(q3)})", "—"],
]
for group in GROUPS:
    table1.append(["AWaRe composition" if group == "Access" else "", group,
                   f"{n(aware[group])}; {pct(aware[group], total_rx)}%", "—"])
access = 100 * aware["Access"] / total_rx
table1 += [
    ["Gap from WHO Access targets", "Former >=60%; 2030 >=70% (percentage points)",
     f"-{60 - access:.2f}; -{70 - access:.2f}", "—"],
    ["", "Access-to-Watch ratio", f"{aware['Access'] / aware['Watch']:.3f}", "—"],
]
for rank, (name, count) in enumerate(classes.head(3).items(), start=1):
    table1.append(["Dominant antibiotic classes" if rank == 1 else "", f"Rank {rank}",
                   f"{name}: {n(count)}; {pct(count, total_rx)}%", "—"])
table1 += [
    ["Diagnostic context", "Samples analysed", "—", n(samples)],
    ["", "Culture-positive samples", "—", f"{n(positives)}; {pct(positives, samples)}%"],
    ["", "Urine among diagnostic samples", "—", f"{n(urine)}; {pct(urine, samples)}%"],
    ["", "Dominant culture-positive organism", "—",
     f"{organisms.index[0]}: {n(organisms.iloc[0])}; {pct(organisms.iloc[0], positives)}%"],
]
for organism, short in [("Escherichia coli", "E. coli"), ("Klebsiella pneumoniae", "K. pneumoniae")]:
    for antibiotic, group in [("Amoxicillin/Clavulanic acid", "Access"), ("Nitrofurantoin", "Access"),
                              ("Amikacin", "Access"), ("Ceftriaxone", "Watch"),
                              ("Ciprofloxacin", "Watch"), ("Meropenem", "Watch")]:
        r, tested = urine_resistance.loc[(organism, antibiotic)]
        first = organism == "Escherichia coli" and antibiotic == "Amoxicillin/Clavulanic acid"
        table1.append(["Urine isolates" if first else "",
                       f"{short}-{antibiotic.lower()} ({group})", "—",
                       f"{n(r)}/{n(tested)}; {pct(r, tested)}%"])

pd.DataFrame(table1, columns=["Domain", "Measure", "Prescription dataset", "Diagnostic dataset"]).to_csv(
    "results/tables/table1.csv", index=False)

# ---------------------------------------------------------------- Supplementary Table 1
pos = dx[dx["positive"]]
s1 = [
    ["Geographic coverage", "States represented", dx["order_state"].nunique(), ""],
    ["", "Districts represented", dx["district"].nunique(), ""],
    ["", "Pincodes represented", dx["Order_Pincode"].nunique(), ""],
    ["Temporal coverage", "Study period", "2023-2025", ""],
    ["", "Total months", dx["month_year"].nunique(), ""],
    ["", "Total possible state-month cells", dx_possible, ""],
    ["", "Observed state-month cells", dx_active, ""],
    ["", "Active state-month proportion (%)", pct(dx_active, dx_possible), ""],
    ["Culture outcome", "Culture positive", positives, pct(positives, samples)],
    ["", "No growth", samples - positives, pct(samples - positives, samples)],
]
s1 += block("Sample type (overall)", dx.groupby("specimen")["Total_Test"].sum().sort_values(ascending=False), samples)
s1 += block("Sample type (culture positive)", pos.groupby("specimen")["Total_Test"].sum().sort_values(ascending=False), positives)
s1 += block("Top contributing states overall", dx.groupby("order_state")["Total_Test"].sum().nlargest(5), samples)
s1 += block("Top contributing states for positive cultures", pos.groupby("order_state")["Total_Test"].sum().nlargest(5), positives)
s1 += block("Most frequently isolated organisms overall", organisms.head(3), positives)
for specimen, top in [("Urine", 3), ("Blood", 2), ("Sputum", 2)]:
    subset = pos[pos["specimen"] == specimen]
    counts = subset.groupby("organism")["Total_Test"].sum().nlargest(top)
    s1 += block(f"Most frequently isolated organisms in {specimen}", counts, subset["Total_Test"].sum())

pd.DataFrame(s1, columns=["Category", "Subcategory", "N", "%"]).to_csv(
    "results/tables/supp_table1_diagnostic_dataset.csv", index=False)

# ---------------------------------------------------------------- Supplementary Table 3
s3 = [
    ["Total antibiotics prescribed", "", total_rx, pct(total_rx, total_medicines)],
    ["Geographic coverage", "States and UTs represented", located["state"].nunique(), ""],
    ["", "Districts represented", rx["district"].nunique(), ""],
    ["", "Pincodes represented", rx["patient_rx_pincode"].nunique(), ""],
    ["", "States retained for state-specific analyses (>=1,000 prescriptions)",
     located.loc[located["state_included"], "state"].nunique(), pct(included, total_rx)],
    ["Temporal coverage", "Study period", "2023-2025", ""],
    ["Median", "Prescription units per state-month", int(median), ""],
    ["IQR", "Prescription units per state-month (Q1-Q3)", f"{int(q1)}-{int(q3)}", ""],
]
s3 += block("Top contributing states", state_totals.head(5), total_rx)
s3 += block("AWaRe classification", aware[GROUPS], total_rx)
s3 += block("Antibiotic classes", classes, total_rx)
s3 += block("States below 1,000 prescriptions (excluded from state analyses)", small_states, total_rx)
s3 += [["No recorded state (national summaries only)", "", no_state, pct(no_state, total_rx)]]

pd.DataFrame(s3, columns=["Category", "Subcategory", "N", "%"]).to_csv(
    "results/tables/supp_table3_prescription_dataset.csv", index=False)

unclassified = (rx[rx["state_included"] & (rx["aware_2025"] == "Unclassified")]
                .groupby("product_clean")["Rx_Count"].sum().sort_values(ascending=False))
unclassified.rename("prescriptions").to_csv("results/tables/supp_table3_unclassified_products.csv")

print(pd.DataFrame(table1).to_string(index=False, header=False))
print(f"\n{len(unclassified)} unclassified products")
