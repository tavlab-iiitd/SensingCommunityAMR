"""Descriptive monthly trends for the ten states with the largest diagnostic volumes."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results/statewise/trend_outputs"
OUT.mkdir(parents=True, exist_ok=True)

START, END = "2023-01", "2025-12"
MONTHS = pd.date_range(START, END, freq="MS")
N_STATES = 10
RANK_BY = "diagnostic_samples"
MIN_SHARED_MONTHS = 1
N_BOOT, SEED = 1000, 42
TOP_PRODUCTS = 10
TOP_DRUGS = 10
MIN_TESTED = 1

GROUPS = ["Access", "Watch", "Reserve", "Not recommended"]
GROUP_COLORS = dict(zip(GROUPS, ["steelblue", "darkorange", "firebrick", "seagreen"]))
ORGANISMS = ["Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa"]
COUNTS = ["RESISTANT", "INTERMEDIATE", "SENSITIVE"]
EXCLUSIONS = {"Klebsiella pneumoniae": ["Ampicillin"]}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


from name_cleaning import clean_organism

def state_name(x):
    if pd.isna(x):
        return pd.NA
    s = " ".join(str(x).strip().lower().replace("&", "and").split())
    aliases = {"pondicherry": "puducherry", "nct delhi": "delhi", "new delhi": "delhi",
               "orissa": "odisha", "telengana": "telangana", "chattisgarh": "chhattisgarh",
               "dadra and nagar haveli": "dadra and nagar haveli and daman and diu",
               "daman and diu": "dadra and nagar haveli and daman and diu"}
    return pd.NA if s in {"", "nan", "none", "unknown", "not available"} else aliases.get(s, s).title()

def pins(s):
    s = s.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)
    return s.where(s.str.fullmatch(r"[1-9]\d{5}", na=False), "Missing")

rx = pd.read_csv(DATA / "prescriptions_clean.csv", usecols=["month_year", "patient_rx_state",
                 "patient_rx_pincode", "Rx_Count", "aware_2025", "product_clean"], low_memory=False)
ast = pd.read_csv(DATA / "ast_clean.csv", usecols=["month_year", "order_state", "Order_Pincode",
                  "test_name", "organism_clean", "antibiotic_clean", *COUNTS], low_memory=False)
dx = pd.read_csv(ROOT / "data/AMR_Pan_India_Aggregated_Diagnostic.csv", encoding="latin1",
                 usecols=["end_state_date", "order_state", "test_name", "organism_isolated", "Total_Test"])
dx["month_year"] = pd.to_datetime(dx["end_state_date"]).dt.strftime("%Y-%m")
for frame, state_col in [(rx, "patient_rx_state"), (ast, "order_state"), (dx, "order_state")]:
    frame["state"] = frame[state_col].map(state_name)
    frame["month"] = pd.to_datetime(frame["month_year"], format="%Y-%m")
rx = rx.loc[rx["month"].isin(MONTHS) & rx["state"].notna()].copy()
ast = ast.loc[ast["month"].isin(MONTHS) & ast["state"].notna()].copy()
dx = dx.loc[dx["month"].isin(MONTHS) & dx["state"].notna()].copy()
rx["pin"] = pins(rx["patient_rx_pincode"])
ast["pin"] = pins(ast["Order_Pincode"])
rx["aware_2025"] = rx["aware_2025"].fillna("Unclassified")
assert rx["Rx_Count"].notna().all() and rx["Rx_Count"].ge(0).all()
assert ast[COUNTS].notna().all().all() and ast[COUNTS].ge(0).all().all()
assert dx["Total_Test"].notna().all() and dx["Total_Test"].ge(0).all()
ast["tested"] = ast[COUNTS].sum(axis=1)
ast = ast.loc[ast["test_name"].eq("URINE") & ast["organism_clean"].isin(ORGANISMS)
              & ast["antibiotic_clean"].notna() & ast["tested"].gt(0)].copy()
labels = {s: clean_organism(s) for s in dx["organism_isolated"].dropna().unique()}
dx["organism"] = dx["organism_isolated"].map(labels)

rx_month = rx.groupby(["state", "month"])["Rx_Count"].sum().rename("antibiotic_records")
dx_month = dx.groupby(["state", "month"])["Total_Test"].sum().rename("diagnostic_samples")
ast_month = ast.groupby(["state", "month"])["tested"].sum().gt(0).astype(int).rename("has_urine_ast")
coverage_month = pd.concat([rx_month, dx_month, ast_month], axis=1).fillna(0)
coverage_month["shared"] = coverage_month["antibiotic_records"].gt(0) & coverage_month["has_urine_ast"].eq(True)
coverage = coverage_month.groupby("state").agg(
    antibiotic_records=("antibiotic_records", "sum"), diagnostic_samples=("diagnostic_samples", "sum"),
    rx_months=("antibiotic_records", lambda x: x.gt(0).sum()),
    diagnostic_months=("diagnostic_samples", lambda x: x.gt(0).sum()),
    urinary_ast_months=("has_urine_ast", "sum"), shared_months=("shared", "sum"))
urine = dx["test_name"].eq("URINE CULTURE AND SENSITIVITY")
coverage["urine_ecoli"] = dx.loc[urine & dx["organism"].eq("Escherichia coli")].groupby("state")["Total_Test"].sum()
coverage = coverage.fillna(0).astype(int).reset_index()
eligible = coverage.loc[coverage["antibiotic_records"].ge(1000)
                         & coverage["shared_months"].ge(MIN_SHARED_MONTHS)].copy()
assert RANK_BY in ["diagnostic_samples", "antibiotic_records", "urine_ecoli"]
selected = eligible.sort_values([RANK_BY, "state"], ascending=[False, True]).head(N_STATES).copy()
if len(selected) < N_STATES:
    raise ValueError(f"Only {len(selected)} states meet the coverage rule; reduce N_STATES or MIN_SHARED_MONTHS.")
selected.insert(0, "rank", range(1, len(selected) + 1))
STATES = selected["state"].tolist()
coverage["selected"] = coverage["state"].isin(STATES)
coverage.to_csv(OUT / "all_state_coverage.csv", index=False)
selected.to_csv(OUT / "top10_states.csv", index=False)
coverage_month.reset_index().to_csv(OUT / "state_month_coverage.csv", index=False)
print(selected)
print("Ranking:", RANK_BY)
print("States with fewer than 30 overlapping months:", selected.loc[selected["shared_months"].lt(30), "state"].tolist())
rx = rx.loc[rx["state"].isin(STATES)].copy()
ast = ast.loc[ast["state"].isin(STATES)].copy()

def rolling(values):
    """Centred, unweighted three-month mean."""
    return pd.DataFrame(values).rolling(3, center=True, min_periods=1).mean().to_numpy()

def monthly_trend(panel, rng, smooth=False):
    """panel: month, pin, count, denominator; one row per pincode-month."""
    result = panel.groupby("month")[["count", "denominator"]].sum().reindex(MONTHS)
    observed = result["denominator"].gt(0).to_numpy()
    result["share_pct"] = 100 * result["count"] / result["denominator"].where(result["denominator"].gt(0))
    result["trend_pct"] = rolling(result[["share_pct"]].to_numpy()).ravel() if smooth else result["share_pct"]
    result.loc[~observed, "trend_pct"] = np.nan
    draws = np.full((len(MONTHS), N_BOOT), np.nan)
    result["n_pincodes"] = 0
    result["unlocated_denominator"] = 0.0
    for i, month in enumerate(MONTHS):
        rows = panel.loc[panel["month"].eq(month) & panel["denominator"].gt(0)]
        located = rows.loc[rows["pin"].ne("Missing"), ["count", "denominator"]].to_numpy()
        fixed = rows.loc[rows["pin"].eq("Missing"), ["count", "denominator"]].sum().to_numpy()
        result.loc[month, "n_pincodes"] = len(located)
        result.loc[month, "unlocated_denominator"] = fixed[1]
        if len(located) < 2:
            continue
        sample = rng.integers(0, len(located), size=(N_BOOT, len(located)))
        totals = located[sample].sum(axis=1) + fixed
        draws[i] = 100 * totals[:, 0] / totals[:, 1]
    usable = np.isfinite(draws).all(axis=1)
    if smooth:
        draws = rolling(draws)

        incomplete = observed & ~usable
        usable &= pd.Series(incomplete).rolling(3, center=True, min_periods=1).sum().eq(0).to_numpy()
    result["ci_low"], result["ci_high"] = np.nan, np.nan
    if usable.any():
        result.loc[usable, ["ci_low", "ci_high"]] = np.percentile(draws[usable], [2.5, 97.5], axis=1).T
    result["n_pincodes"] = result["n_pincodes"].astype(int)
    return result.rename_axis("month").reset_index()

rx_den = rx.groupby(["state", "month", "pin"], as_index=False)["Rx_Count"].sum().rename(columns={"Rx_Count": "denominator"})
product_totals = rx.groupby(["aware_2025", "product_clean"], as_index=False)["Rx_Count"].sum()
top_products = (product_totals.loc[product_totals["aware_2025"].isin(GROUPS[:3])]
                .sort_values(["aware_2025", "Rx_Count", "product_clean"], ascending=[True, False, True])
                .groupby("aware_2025", sort=False).head(TOP_PRODUCTS))
top_products.to_csv(OUT / "selected_prescription_products.csv", index=False)

def prescription_trends(column, series):
    counts = rx.groupby(["state", "month", "pin", column], as_index=False)["Rx_Count"].sum()
    rng, output = np.random.default_rng(SEED), []
    for state in STATES:
        denominator = rx_den.loc[rx_den["state"].eq(state), ["month", "pin", "denominator"]]
        for name in series:
            numerator = counts.loc[counts["state"].eq(state) & counts[column].eq(name), ["month", "pin", "Rx_Count"]]
            panel = denominator.merge(numerator, on=["month", "pin"], how="left", validate="one_to_one")
            panel["count"] = panel["Rx_Count"].fillna(0)
            output.append(monthly_trend(panel, rng, smooth=True).assign(state=state, series=name))
    return pd.concat(output, ignore_index=True)

rx_aware = prescription_trends("aware_2025", GROUPS)
rx_products = prescription_trends("product_clean", top_products["product_clean"].tolist())
rx_products = rx_products.merge(top_products[["product_clean", "aware_2025"]],
                                left_on="series", right_on="product_clean", validate="many_to_one").drop(columns="product_clean")
rx_aware.to_csv(OUT / "prescription_aware_monthly.csv", index=False)
rx_products.to_csv(OUT / "prescription_product_monthly.csv", index=False)
print(rx_aware.head())

ast_plot = ast.copy()
for organism, excluded in EXCLUSIONS.items():
    ast_plot = ast_plot.loc[~(ast_plot["organism_clean"].eq(organism) & ast_plot["antibiotic_clean"].isin(excluded))]
drug_totals = ast_plot.groupby(["organism_clean", "antibiotic_clean"], as_index=False)["tested"].sum()
top_drugs = (drug_totals.sort_values(["organism_clean", "tested", "antibiotic_clean"], ascending=[True, False, True])
             .groupby("organism_clean", sort=False).head(TOP_DRUGS))
top_drugs.to_csv(OUT / "selected_resistance_drugs.csv", index=False)
ast_panel = ast_plot.groupby(["state", "organism_clean", "antibiotic_clean", "month", "pin"], as_index=False)[["RESISTANT", "tested"]].sum()
ast_panel = ast_panel.rename(columns={"RESISTANT": "count", "tested": "denominator"})
rng, results = np.random.default_rng(SEED), []
for organism in ORGANISMS:
    drugs = top_drugs.loc[top_drugs["organism_clean"].eq(organism), "antibiotic_clean"]
    for state in STATES:
        subset = ast_panel.loc[ast_panel["state"].eq(state) & ast_panel["organism_clean"].eq(organism)]
        for drug in drugs:
            panel = subset.loc[subset["antibiotic_clean"].eq(drug)]
            results.append(monthly_trend(panel, rng).assign(state=state, organism=organism, series=drug))
resistance = pd.concat(results, ignore_index=True).rename(
    columns={"count": "R", "denominator": "N", "share_pct": "R_percent"})
resistance.to_csv(OUT / "resistance_monthly.csv", index=False)
resistance_coverage = resistance.groupby(["state", "organism", "series"], as_index=False).agg(
    tested=("N", "sum"), resistant=("R", "sum"),
    months_tested=("N", lambda x: x.gt(0).sum()),
    months_ge30=("N", lambda x: x.ge(30).sum()),
    min_tested_observed=("N", "min"), max_tested_observed=("N", "max"))
resistance_coverage.to_csv(OUT / "resistance_coverage.csv", index=False)
print(resistance_coverage.head(10))

def plot_states(data, series_order, title, ylabel, filename, colors=None, resistance_plot=False):
    nrows = (len(STATES) + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(15, 3.6 * nrows), sharex=True, sharey=True, squeeze=False)
    palette = colors or {name: plt.cm.tab10.colors[i % 10] for i, name in enumerate(series_order)}
    for ax, state in zip(axes.flat, STATES):
        state_data = data.loc[data["state"].eq(state)].copy()
        if resistance_plot:
            state_data.loc[state_data["N"].fillna(0).lt(MIN_TESTED), ["trend_pct", "ci_low", "ci_high"]] = np.nan
        for name in series_order:
            d = state_data.loc[state_data["series"].eq(name)].sort_values("month").copy()
            color = palette[name]
            ax.plot(d["month"].to_numpy(), d["trend_pct"].to_numpy(), color=color, lw=1.3,
                    marker="o" if resistance_plot else None, markersize=2.4, label=name)
            ax.fill_between(d["month"].to_numpy(), d["ci_low"].to_numpy(), d["ci_high"].to_numpy(),
                            color=color, alpha=0.10, lw=0)
        n_months = state_data.groupby("month")["trend_pct"].count().gt(0).sum()
        ax.set_title(f"{state} | {n_months}/{len(MONTHS)} months", loc="left", fontsize=11)
        if not state_data["trend_pct"].notna().any():
            ax.text(0.5, 0.5, "No results meeting the plotting threshold", transform=ax.transAxes, ha="center", color="0.4")
        ax.set_ylim(-2, 102 if resistance_plot else max(5, data["ci_high"].max() * 1.08, data["trend_pct"].max() * 1.08))
        ax.set_xlim(MONTHS[0], MONTHS[-1])
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
        ax.grid(alpha=0.18)
        ax.tick_params(axis="x", rotation=35)
    for ax in axes.flat[len(STATES):]:
        ax.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2 if len(series_order) > 4 else 4,
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.015))
    if resistance_plot and MIN_TESTED > 1:
        title += f" (N ≥ {MIN_TESTED})"
    fig.suptitle(title, fontsize=17, y=0.993)
    fig.supylabel(ylabel, x=0.008)
    fig.subplots_adjust(top=0.965, bottom=0.12 if len(series_order) > 4 else 0.07,
                        left=0.06, right=0.99, hspace=0.27, wspace=0.11)
    fig.savefig(OUT / f"{filename}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

plot_states(rx_aware, GROUPS, "AWaRe prescription shares | centred 3-month mean",
            "% of all antibiotic prescription records", "prescription_aware_top10", GROUP_COLORS)

for group in GROUPS[:3]:
    products = top_products.loc[top_products["aware_2025"].eq(group), "product_clean"].tolist()
    plot_states(rx_products.loc[rx_products["aware_2025"].eq(group)], products,
                f"{group} products | centred 3-month mean", "% of all antibiotic prescription records",
                f"prescription_products_{group.lower()}_top10")

for organism in ORGANISMS:
    drugs = top_drugs.loc[top_drugs["organism_clean"].eq(organism), "antibiotic_clean"].tolist()
    plot_states(resistance.loc[resistance["organism"].eq(organism)], drugs,
                f"{organism} | urinary resistance, observed monthly percentages",
                "Resistance: 100 × R / (R + I + S)", f"resistance_{organism.replace(' ', '_').lower()}_top10",
                resistance_plot=True)

med = pd.read_csv(ROOT / "data/State_Month_Prescription_updated.csv")
med["state"] = med["Rx_State"].map(state_name)
med["month"] = pd.to_datetime(med["month_year"], format="%Y-%m")
med = med.loc[med["state"].isin(STATES) & med["month"].isin(MONTHS)]
med = med.groupby(["state", "month"])["presc"].sum().rename("all_medicines")
grid = pd.MultiIndex.from_product([STATES, MONTHS], names=["state", "month"])
all_med = pd.concat([rx_month.rename("antibiotic_records"), med], axis=1).reindex(grid).reset_index()
all_med["share_pct"] = 100 * all_med["antibiotic_records"] / all_med["all_medicines"].where(all_med["all_medicines"].gt(0))
all_med["trend_pct"] = all_med.groupby("state", sort=False)["share_pct"].transform(
    lambda s: s.rolling(3, center=True, min_periods=1).mean())
all_med.loc[all_med["share_pct"].isna(), "trend_pct"] = np.nan
all_med["ci_low"], all_med["ci_high"] = np.nan, np.nan
all_med["series"] = "All antibiotics"
all_med.to_csv(OUT / "prescription_share_all_medicines.csv", index=False)
plot_states(all_med, ["All antibiotics"], "Antibiotics as a share of all medicines | centred 3-month mean",
            "% of all medicine prescription records", "prescription_all_medicines_top10")

assert len(STATES) == N_STATES and len(set(STATES)) == N_STATES
assert selected["shared_months"].ge(MIN_SHARED_MONTHS).all()
assert (rx_aware["count"] <= rx_aware["denominator"]).all()
assert rx_aware["trend_pct"].dropna().between(0, 100).all()
assert resistance["R_percent"].dropna().between(0, 100).all()
assert resistance.loc[resistance["N"].isna(), "R_percent"].isna().all()
assert all_med["share_pct"].dropna().between(0, 100).all()


direct_rx = rx.groupby(["state", "month", "aware_2025"])["Rx_Count"].sum()
check_rx = rx_aware.set_index(["state", "month", "series"])
np.testing.assert_allclose(check_rx["count"], direct_rx.reindex(check_rx.index).fillna(0))
direct_den = rx.groupby(["state", "month"])["Rx_Count"].sum()
pd.testing.assert_series_equal(rx_aware.groupby(["state", "month"])["denominator"].first().sort_index(),
                               direct_den.sort_index(), check_names=False, check_dtype=False)
direct_ast = ast.groupby(["state", "organism_clean", "antibiotic_clean", "month"])[COUNTS].sum()
check_ast = resistance.set_index(["state", "organism", "series", "month"])
matched = direct_ast.reindex(check_ast.index)
np.testing.assert_allclose(check_ast["R"], matched["RESISTANT"], equal_nan=True)
np.testing.assert_allclose(check_ast["N"], matched.sum(axis=1, min_count=1), equal_nan=True)
summary = resistance_coverage.groupby(["state", "organism"])["months_tested"].agg(["min", "max"]).reindex(
    pd.MultiIndex.from_product([STATES, ORGANISMS], names=["state", "organism"]))
print(summary.rename(columns={"min": "fewest_months_for_a_drug", "max": "most_months_for_a_drug"}))
print("Checks passed. Outputs:", OUT)
print("CSV files:", ", ".join(p.name for p in sorted(OUT.glob("*.csv"))))
