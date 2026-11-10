"""Monthly urinary resistance with state bootstrap intervals."""

import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)
COUNTS = ["RESISTANT", "INTERMEDIATE", "SENSITIVE"]
ORGANISMS = {"escherichia coli": 5, "klebsiella pneumoniae": 6, "pseudomonas aeruginosa": 7}
PLOT_EXCLUSIONS = {"klebsiella pneumoniae": {"ampicillin"}}

ast = pd.read_csv("data/ast_clean.csv", low_memory=False)
ast = ast[ast["test_name"] == "URINE"].copy()
ast["organism"] = ast["organism_clean"].str.lower()
ast["antibiotic"] = ast["antibiotic_clean"].str.lower()
ast["tested"] = ast[COUNTS].sum(axis=1)


def state_bootstrap(resistant, tested, rng, n_boot=1000):
    """Resample whole states; return the 2.5th and 97.5th percentiles of % resistant."""
    draws = rng.integers(0, len(tested), size=(n_boot, len(tested)))
    r, n = resistant[draws].sum(axis=1), tested[draws].sum(axis=1)
    boot = r[n > 0] / n[n > 0] * 100
    return np.percentile(boot, 2.5), np.percentile(boot, 97.5)


rng = np.random.default_rng(42)
replacement_rng = np.random.default_rng(42)
results = []
for organism, figure_number in ORGANISMS.items():
    d = ast[ast["organism"] == organism]
    ranked = d.groupby("antibiotic")["tested"].sum().sort_values(ascending=False)
    original_top10 = ranked.head(10).index.tolist()
    excluded = PLOT_EXCLUSIONS.get(organism, set())
    top10 = ranked.loc[~ranked.index.isin(excluded)].head(10).index.tolist()


    d = d[d["antibiotic"].isin(set(original_top10) | set(top10))]
    state_counts = d.groupby(["month_year", "antibiotic", "order_state"], dropna=False)[COUNTS].sum().reset_index()
    state_counts["tested"] = state_counts[COUNTS].sum(axis=1)

    rows = []
    for (month, antibiotic), g in state_counts.groupby(["month_year", "antibiotic"]):
        states = g.groupby("order_state")[["RESISTANT", "tested"]].sum()
        resistant, tested = states["RESISTANT"].to_numpy(), states["tested"].to_numpy()
        if tested.sum() == 0:
            continue
        bootstrap_rng = rng if antibiotic in original_top10 else replacement_rng
        ci_low, ci_high = state_bootstrap(resistant, tested, bootstrap_rng)
        if antibiotic not in top10:
            continue
        rows.append({"organism": organism, "antibiotic": antibiotic, "month_year": month,
                     "R_percent": resistant.sum() / tested.sum() * 100, "ci_low": ci_low, "ci_high": ci_high,
                     "N": int(tested.sum()), "R": int(resistant.sum()), "n_states": len(states)})
    monthly = pd.DataFrame(rows).sort_values(["antibiotic", "month_year"])
    results.append(monthly)

    months = sorted(monthly["month_year"].unique())
    fig, ax = plt.subplots(figsize=(14, 6))
    for i, antibiotic in enumerate(top10):
        a = monthly[monthly["antibiotic"] == antibiotic]
        x = [months.index(m) for m in a["month_year"]]
        color = plt.cm.tab10.colors[i % 10]
        ax.plot(x, a["R_percent"], marker="o", label=antibiotic.title(), color=color, linewidth=1.8)
        ax.fill_between(x, a["ci_low"], a["ci_high"], color=color, alpha=0.15)
    ax.set_xticks(range(len(months)))
    ax.set_xticklabels(months, rotation=60, ha="right")
    ax.set_ylim(0, 105)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_ylabel("Resistance %")
    ax.set_xlabel("Month-Year")
    ax.set_title(f"Monthly Resistance Trends: {organism.title()} — URINE", fontsize=11)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize="small", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"results/figures/supp_fig{figure_number}_{organism.replace(' ', '_')}_urine_resistance.png",
                bbox_inches="tight", dpi=300)
    if organism == "klebsiella pneumoniae":
        plt.savefig(f"results/figures/supp_fig{figure_number}_{organism.replace(' ', '_')}_urine_resistance.pdf",
                    bbox_inches="tight")
    plt.close()

pd.concat(results).to_csv("results/figures/supp_fig5-7_resistance_monthly_data.csv", index=False)
