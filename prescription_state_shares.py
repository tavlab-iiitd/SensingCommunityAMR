"""Figure 1 (Access share by state) and Supplementary Figure 3 (Watch share by state).

Share of all antibiotic prescriptions in each of the 25 states/UTs with >=1,000
antibiotic prescriptions. 95% intervals: 2,000 bootstrap resamples of pincodes within each state.

Input : data/prescriptions_clean.csv
Output: results/figures/fig1_state_access.(png/pdf), supp_fig3_state_watch.(png/pdf) and their data
"""

import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)

rx = pd.read_csv("data/prescriptions_clean.csv", low_memory=False)
rx = rx[rx["state_included"]].copy()
rx["state"] = rx["patient_rx_state"].str.title()
rx["pincode"] = rx["patient_rx_pincode"].astype("Int64").astype(str).replace("<NA>", "Missing pincode")


def state_share(group, n_boot=2000, seed=42):
    clusters = rx.groupby(["state", "pincode"], as_index=False)["Rx_Count"].sum()
    in_group = rx[rx["aware_2025"] == group].groupby(["state", "pincode"], as_index=False)["Rx_Count"].sum()
    clusters = clusters.merge(in_group, on=["state", "pincode"], how="left", suffixes=("_total", "_group"))
    clusters["Rx_Count_group"] = clusters["Rx_Count_group"].fillna(0)

    rng = np.random.default_rng(seed)
    rows = []
    for state, d in clusters.groupby("state"):
        group_rx, total_rx = d["Rx_Count_group"].to_numpy(), d["Rx_Count_total"].to_numpy()
        draws = rng.integers(0, len(d), size=(n_boot, len(d)))
        boot = group_rx[draws].sum(axis=1) / total_rx[draws].sum(axis=1) * 100
        rows.append({"state": state, "prescriptions": int(total_rx.sum()), "pincodes": len(d),
                     "share_pct": 100 * group_rx.sum() / total_rx.sum(),
                     "ci_low": np.percentile(boot, 2.5), "ci_high": np.percentile(boot, 97.5)})
    return pd.DataFrame(rows).sort_values("share_pct").reset_index(drop=True)


# Figure 1: Access, with the former WHO 60% benchmark and the 2030 70% target
access = state_share("Access")
access.to_csv("results/figures/fig1_state_access_data.csv", index=False)

fig, ax = plt.subplots(figsize=(11, 11))
y = access.index
ax.hlines(y, 0, access["share_pct"], color="#d9d9d9", linewidth=1.4)
ax.hlines(y, access["ci_low"], access["ci_high"], color="#4a7fb5", linewidth=3.2, alpha=0.4)
ax.scatter(access["share_pct"], y, color="#1f4e79", s=85, edgecolors="white", linewidths=0.6, zorder=3)
ax.axvline(60, color="#777777", linestyle=":", linewidth=1.3)
ax.axvline(70, color="#c0392b", linestyle="--", linewidth=1.6)
ax.text(59.2, len(access) - 0.1, "Former WHO ≥60%", fontsize=8, ha="right", fontweight="bold", color="#666666")
ax.text(70.8, len(access) - 0.1, "UNGA 2030 ≥70%", fontsize=8, ha="left", fontweight="bold", color="#c0392b")
for i, row in access.iterrows():
    ax.text(76, i, f"{row.share_pct:.1f}  [{row.ci_low:.1f}–{row.ci_high:.1f}]", va="center", fontsize=8)
ax.set_yticks(y)
ax.set_yticklabels(access["state"], fontsize=9.5)
ax.set_ylim(-1.4, len(access) + 0.6)
ax.set_xlim(0, 95)
ax.set_xlabel("Access antibiotics as a share of total antibiotic prescriptions, %")
ax.grid(axis="x", alpha=0.28, linestyle=":")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("State-level variation in Access antibiotic use across the digital\n"
             "pharmacy platform, India, 2023–2025", fontsize=13, fontweight="bold", pad=16, loc="left")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"results/figures/fig1_state_access.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)

# Supplementary Figure 3: Watch
watch = state_share("Watch")
watch.to_csv("results/figures/supp_fig3_state_watch_data.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 12))
y = watch.index
ax.hlines(y, 0, watch["share_pct"], color="#cccccc", linewidth=1.3)
ax.hlines(y, watch["ci_low"], watch["ci_high"], color="#4a7fb5", linewidth=3, alpha=0.4)
ax.scatter(watch["share_pct"], y, color="#1f4e79", s=85, edgecolors="white", linewidths=0.5, zorder=3)
for i, row in watch.iterrows():
    ax.text(row.ci_high + 1.2, i, f"{row.share_pct:.1f}%  [{row.ci_low:.1f}–{row.ci_high:.1f}]",
            va="center", fontsize=7, color="#444444")
ax.set_yticks(y)
ax.set_yticklabels(watch["state"], fontsize=9)
ax.set_xlim(0, watch["ci_high"].max() + 22)
ax.set_xlabel("Watch antibiotics as a share of total antibiotic prescriptions, %")
ax.grid(axis="x", alpha=0.3, linestyle=":")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("State-level Watch antibiotic share, India, 2023–2025", fontsize=12, fontweight="bold", pad=12)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"results/figures/supp_fig3_state_watch.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Access {access['share_pct'].min():.1f}%-{access['share_pct'].max():.1f}%; "
      f"Watch {watch['share_pct'].min():.1f}%-{watch['share_pct'].max():.1f}%")
