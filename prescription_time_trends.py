"""Monthly prescription trends: Figure 2 and Supplementary Figures 2 and 4.

Shares are calculated in the 25 states/UTs with >=1,000 antibiotic prescriptions and smoothed
with a centred 3-month rolling mean. 95% intervals: 1,000 bootstrap resamples of states within
each month, smoothed in the same way.

Figure 2   : AWaRe groups and top products as % of all antibiotic prescriptions
Supp Fig 2 : antibiotics as % of all medicines (a) and per 100 non-antibiotic medicines (b)
Supp Fig 4 : as Figure 2, but as % of all medicine prescriptions

Input : data/prescriptions_clean.csv, data/State_Month_Prescription_updated.csv
Output: results/figures/fig2_*, supp_fig2_*, supp_fig4_* (png, pdf and data csv)
"""

import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)
GROUPS = ["Access", "Watch", "Reserve", "Not recommended"]
COLORS = {"Access": "steelblue", "Watch": "darkorange", "Reserve": "firebrick", "Not recommended": "seagreen"}
TOP_N = 10

rx = pd.read_csv("data/prescriptions_clean.csv", low_memory=False)
rx = rx[rx["state_included"]].copy()
rx["state"] = rx["patient_rx_state"].str.title()
rx["month"] = pd.to_datetime(rx["month_year"], format="%Y-%m")

medicines = pd.read_csv("data/State_Month_Prescription_updated.csv")
medicines = medicines[medicines["month_year"].between("2023-01", "2025-12")].copy()
medicines["state"] = medicines["Rx_State"].str.strip().str.title()
medicines["month"] = pd.to_datetime(medicines["month_year"], format="%Y-%m")
medicines = medicines.groupby(["month", "state"], as_index=False)["presc"].sum()
medicines = medicines[medicines["state"].isin(rx["state"])].rename(columns={"presc": "denominator"})

antibiotics = rx.groupby(["month", "state"], as_index=False)["Rx_Count"].sum().rename(
    columns={"Rx_Count": "denominator"})


def rolling(x):
    return pd.Series(x).rolling(3, center=True, min_periods=1).mean().to_numpy()


def monthly_trend(panel, rng, n_boot=1000):
    """Smoothed monthly share (%) of `count` in `denominator`, with a bootstrap over states.

    `panel` has one row per state-month. Draw order (replicate, then month) is kept fixed.
    """
    trend = panel.groupby("month", as_index=False)[["count", "denominator"]].sum()
    trend["share_pct"] = 100 * trend["count"] / trend["denominator"]
    trend["rolling_mean"] = rolling(trend["share_pct"])

    by_month = {m: d[["count", "denominator"]].to_numpy() for m, d in panel.groupby("month")}
    boot = np.full((n_boot, len(trend)), np.nan)
    for b in range(n_boot):
        share = np.full(len(trend), np.nan)
        for i, month in enumerate(trend["month"]):
            rows = by_month[month]
            sample = rows[rng.integers(0, len(rows), size=len(rows))]
            if sample[:, 1].sum() > 0:
                share[i] = 100 * sample[:, 0].sum() / sample[:, 1].sum()
        boot[b] = rolling(share)
    trend["ci_low"] = np.nanpercentile(boot, 2.5, axis=0)
    trend["ci_high"] = np.nanpercentile(boot, 97.5, axis=0)
    return trend


def group_trends(denominator):
    """Panel a: AWaRe groups. Panels b-d: top products of Access, Watch and Reserve."""
    counts = rx.groupby(["month", "state", "aware_2025"], as_index=False)["Rx_Count"].sum()
    panel = (pd.DataFrame({"aware_2025": GROUPS}).merge(denominator, how="cross")
             .merge(counts, on=["month", "state", "aware_2025"], how="left"))
    panel["count"] = panel["Rx_Count"].fillna(0)
    rng = np.random.default_rng(42)
    groups = pd.concat([monthly_trend(panel[panel["aware_2025"] == g], rng).assign(series=g, panel="a")
                        for g in GROUPS])

    totals = rx.groupby(["aware_2025", "product_clean"], as_index=False)["Rx_Count"].sum()
    totals = totals.sort_values(["aware_2025", "Rx_Count"], ascending=[True, False])
    top = totals.groupby("aware_2025").head(TOP_N)
    counts = rx.groupby(["month", "state", "product_clean"], as_index=False)["Rx_Count"].sum()
    panel = (top[["aware_2025", "product_clean"]].merge(denominator, how="cross")
             .merge(counts, on=["month", "state", "product_clean"], how="left"))
    panel["count"] = panel["Rx_Count"].fillna(0)
    rng = np.random.default_rng(42)
    products = []
    for letter, group in zip("bcd", GROUPS[:3]):
        for product in top.loc[top["aware_2025"] == group, "product_clean"]:
            trend = monthly_trend(panel[panel["product_clean"] == product], rng)
            products.append(trend.assign(series=product, panel=letter))
    return pd.concat([groups] + products, ignore_index=True)


def draw(ax, data, title, ylabel):
    for series, d in data.groupby("series", sort=False):
        line = ax.plot(d["month"], d["rolling_mean"], marker="o", markersize=4, linewidth=1.6,
                       label=series, color=COLORS.get(series))[0]
        ax.fill_between(d["month"], d["ci_low"], d["ci_high"], alpha=0.2, linewidth=0, color=line.get_color())
    ax.set_title(title, fontsize=15)
    ax.set_xlabel("Month", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.tick_params(axis="x", rotation=45, labelsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=12, frameon=False)


def save_trend_figure(data, name, denominator_label):
    ylabel = f"% of {denominator_label.lower()}"
    titles = {"a": f"AWaRe Category Share of {denominator_label} (3-month rolling)",
              "b": f"Share of {denominator_label} (3-month rolling): Access",
              "c": f"Share of {denominator_label} (3-month rolling): Watch",
              "d": f"Share of {denominator_label} (3-month rolling): Reserve"}
    fig, axes = plt.subplots(2, 2, figsize=(26, 16))
    for ax, letter in zip(axes.flat, "abcd"):
        draw(ax, data[data["panel"] == letter], titles[letter], ylabel)
        ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=22, fontweight="bold")
    fig.subplots_adjust(hspace=0.75, wspace=0.18)
    for ext in ["png", "pdf"]:
        fig.savefig(f"results/figures/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    data.to_csv(f"results/figures/{name}_data.csv", index=False)


# Figure 2: denominator = all antibiotic prescriptions
save_trend_figure(group_trends(antibiotics), "fig2_prescription_trends", "Total Antibiotic Prescriptions")

# Supplementary Figure 4: denominator = all medicine prescriptions
save_trend_figure(group_trends(medicines), "supp_fig4_prescription_trends_all_medicines",
                  "All Medicine Prescriptions")

# Supplementary Figure 2: antibiotics among all medicines (a) and per 100 non-antibiotic medicines (b)
panel = antibiotics.rename(columns={"denominator": "count"}).merge(medicines, on=["month", "state"])
share = monthly_trend(panel, np.random.default_rng(42))
panel["denominator"] = panel["denominator"] - panel["count"]
per_100 = monthly_trend(panel, np.random.default_rng(42))

fig, axes = plt.subplots(2, 1, figsize=(12, 10))
for ax, data, title, ylabel in [
        (axes[0], share, "a  Antibiotic share of all medicine prescriptions (3-month rolling average)",
         "Antibiotic prescriptions as a\npercentage of all medicines"),
        (axes[1], per_100, "b  Antibiotic prescriptions per 100 non-antibiotic medicine prescriptions",
         "Antibiotic prescriptions per 100\nnon-antibiotic medicine prescriptions")]:
    ax.plot(data["month"], data["rolling_mean"], marker="o", linewidth=1.5, label="Antibiotic prescriptions")
    ax.fill_between(data["month"], data["ci_low"], data["ci_high"], alpha=0.2, linewidth=0, label="95% CI")
    ax.set_title(title, loc="left")
    ax.set_xlabel("Month")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=60)
    ax.grid(alpha=0.3)
    ax.legend()
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"results/figures/supp_fig2_antibiotic_share.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
pd.concat([share.assign(panel="a"), per_100.assign(panel="b")]).to_csv(
    "results/figures/supp_fig2_antibiotic_share_data.csv", index=False)

print(f"Monthly antibiotic share in the 25 states/UTs: {share['share_pct'].min():.2f}%-{share['share_pct'].max():.2f}%")
