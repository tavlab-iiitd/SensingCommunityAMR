"""Supplementary Figure 1: geographic and temporal footprint of the two data streams.

(a) antibiotic prescriptions and (b) diagnostic culture samples by state/UT (log colour scale),
(c) months with any recorded activity in each state/UT. These are platform-recorded counts,
not population rates.

Input : data/prescriptions_clean.csv, data/AMR_Pan_India_Aggregated_Diagnostic.csv,
        geodata/geoBoundaries-IND-ADM1_simplified.geojson (geoBoundaries, CC BY 2.5 IN)
Output: results/figures/supp_fig1_geographic_footprint.(png/pdf) and _data.csv
"""

import json
import math
import os
import re
import unicodedata

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FuncFormatter, LogLocator

os.makedirs("results/figures", exist_ok=True)

# Spelling variants -> the state name used in the boundary file
ALIASES = {
    "andaman nicobar islands": "andaman and nicobar islands",
    "andaman and nicobar": "andaman and nicobar islands",
    "chattisgarh": "chhattisgarh",
    "dadra nagar haveli daman diu": "dadra and nagar haveli and daman and diu",
    "dadra nagar haveli": "dadra and nagar haveli and daman and diu",
    "dadra and nagar haveli": "dadra and nagar haveli and daman and diu",
    "dadra and nagar haveli daman and diu": "dadra and nagar haveli and daman and diu",
    "dadra and nagar haveli and daman diu": "dadra and nagar haveli and daman and diu",
    "daman diu": "dadra and nagar haveli and daman and diu",
    "daman and diu": "dadra and nagar haveli and daman and diu",
    "jammu kashmir": "jammu and kashmir",
    "nct delhi": "delhi",
    "nct of delhi": "delhi",
    "new delhi": "delhi",
    "orissa": "odisha",
    "pondicherry": "puducherry",
    "telengana": "telangana",
    "uttaranchal": "uttarakhand",
}


def state_key(name):
    if pd.isna(name):
        return None
    text = unicodedata.normalize("NFKD", str(name)).encode("ascii", "ignore").decode()
    text = re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()
    if text in {"", "nan", "none", "unknown", "not available"}:
        return None
    return ALIASES.get(text, text)


rx = pd.read_csv("data/prescriptions_clean.csv", usecols=["month_year", "patient_rx_state", "Rx_Count",
                                                           "state_included"], low_memory=False)
rx["state"] = rx["patient_rx_state"].map(state_key)
rx = rx[rx["state"].notna()]

dx = pd.read_csv("data/AMR_Pan_India_Aggregated_Diagnostic.csv", encoding="ISO-8859-1",
                 usecols=["end_state_date", "order_state", "Total_Test"], low_memory=False)
dx["month_year"] = pd.to_datetime(dx["end_state_date"]).dt.strftime("%Y-%m")
dx = dx[dx["month_year"].between("2023-01", "2025-12")]
dx["state"] = dx["order_state"].map(state_key)
dx = dx[dx["state"].notna()]

geo = json.load(open("geodata/geoBoundaries-IND-ADM1_simplified.geojson", encoding="utf-8"))["features"]
names = {state_key(f["properties"]["shapeName"]): f["properties"]["shapeName"] for f in geo}
assert set(rx["state"]) | set(dx["state"]) <= set(names), "state names missing from the boundary file"

summary = pd.DataFrame({"state": sorted(names)})
summary["state_name"] = summary["state"].map(names)
summary = summary.merge(rx.groupby("state").agg(prescriptions=("Rx_Count", "sum"),
                                                included=("state_included", "any")).reset_index(),
                        on="state", how="left")
summary = summary.merge(dx.groupby("state").agg(diagnostic_samples=("Total_Test", "sum")).reset_index(),
                        on="state", how="left")
summary[["prescriptions", "diagnostic_samples"]] = summary[["prescriptions", "diagnostic_samples"]].fillna(0)
summary["included"] = summary["included"].eq(True)
summary.to_csv("results/figures/supp_fig1_geographic_footprint_data.csv", index=False)

# Monthly presence: 1 = prescriptions only, 2 = diagnostics only, 3 = both
months = sorted(set(rx["month_year"]) | set(dx["month_year"]))
active = summary[(summary["prescriptions"] > 0) | (summary["diagnostic_samples"] > 0)]
rx_active = set(zip(rx.loc[rx["Rx_Count"] > 0, "state"], rx.loc[rx["Rx_Count"] > 0, "month_year"]))
dx_active = set(zip(dx.loc[dx["Total_Test"] > 0, "state"], dx.loc[dx["Total_Test"] > 0, "month_year"]))
presence = np.array([[((s, m) in rx_active) + 2 * ((s, m) in dx_active) for m in months] for s in active["state"]])


def polygons(geometry):
    shapes = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]
    return [np.asarray(shape[0]) for shape in shapes]


def choropleth(ax, column, title, cmap_name, outline_included=False):
    values = summary.set_index("state")[column]
    included = summary.set_index("state")["included"]
    norm = LogNorm(vmin=max(1.0, values[values > 0].min()), vmax=values.max())
    cmap = plt.get_cmap(cmap_name)
    for feature in geo:
        key = state_key(feature["properties"]["shapeName"])
        bold = outline_included and included[key]
        for ring in polygons(feature["geometry"]):
            ax.add_patch(Polygon(ring, closed=True, rasterized=True,
                                 facecolor=cmap(norm(values[key])) if values[key] > 0 else "#eeeeee",
                                 edgecolor="#202020" if bold else "#777777", linewidth=0.75 if bold else 0.28))
    xy = np.vstack([ring for f in geo for ring in polygons(f["geometry"])])
    (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
    ax.set_xlim(x0 - 0.02 * (x1 - x0), x1 + 0.02 * (x1 - x0))
    ax.set_ylim(y0 - 0.02 * (y1 - y0), y1 + 0.02 * (y1 - y0))
    ax.set_aspect(1 / math.cos(math.radians((y0 + y1) / 2)))
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=7)
    bar = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal", fraction=0.035, pad=0.015)
    bar.set_label("Platform-recorded events (log scale)", fontsize=8)
    bar.locator = LogLocator(base=10, subs=(1.0,))
    bar.formatter = FuncFormatter(lambda v, _: f"{v / 1e6:g}M" if v >= 1e6 else f"{v / 1e3:g}k" if v >= 1e3 else f"{v:g}")
    bar.update_ticks()
    bar.ax.tick_params(labelsize=7, length=2)
    if outline_included:
        ax.legend(handles=[Patch(facecolor="none", edgecolor="#202020", linewidth=1.2, label="Included in state analysis")],
                  loc="lower left", frameon=False, fontsize=7.5, bbox_to_anchor=(0.0, 0.005))


fig = plt.figure(figsize=(13.2, 13.6))
grid = fig.add_gridspec(2, 2, height_ratios=[1.02, 1.15], hspace=0.18, wspace=0.08)
choropleth(fig.add_subplot(grid[0, 0]), "prescriptions", "a  Antibiotic prescription footprint", "Blues", True)
choropleth(fig.add_subplot(grid[0, 1]), "diagnostic_samples", "b  Diagnostic culture footprint", "Oranges")

ax = fig.add_subplot(grid[1, :])
colors = ListedColormap(["#f0f0f0", "#4c78a8", "#f58518", "#54a24b"])
ax.imshow(presence, aspect="auto", cmap=colors, norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4), interpolation="nearest")
ax.set_yticks(range(len(active)))
ax.set_yticklabels([f"● {n}" if inc else f"  {n}" for n, inc in zip(active["state_name"], active["included"])], fontsize=6.9)
ax.set_xticks(range(0, len(months), 3))
ax.set_xticklabels(months[::3], rotation=45, ha="right", fontsize=7)
for i, month in enumerate(months):
    if month.endswith("-01") and i > 0:
        ax.axvline(i - 0.5, color="white", linewidth=1.2)
ax.tick_params(length=0)
ax.set_title("c  Monthly presence of platform-recorded activity", loc="left", fontsize=12, fontweight="bold", pad=8)
ax.set_xlabel("Calendar month", fontsize=9)
ax.set_ylabel("State/UT (● included in prescription state analysis)", fontsize=9)
ax.legend(handles=[Patch(facecolor="#f0f0f0", edgecolor="#cccccc", label="No recorded activity"),
                   Patch(facecolor="#4c78a8", label="Prescription only"),
                   Patch(facecolor="#f58518", label="Diagnostic only"),
                   Patch(facecolor="#54a24b", label="Both streams")],
          ncol=4, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.11))

fig.suptitle("Geographic and temporal footprint of the online pharmacy platform",
             x=0.055, y=0.985, ha="left", fontsize=15, fontweight="bold")
fig.text(0.055, 0.962, "Platform-recorded events, not population rates", ha="left", fontsize=10, color="#555555")
fig.subplots_adjust(left=0.12, right=0.98, top=0.945, bottom=0.085)
for ext in ["png", "pdf"]:
    fig.savefig(f"results/figures/supp_fig1_geographic_footprint.{ext}", dpi=450, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"Prescription states/UTs: {(summary['prescriptions'] > 0).sum()} map features "
      f"({summary['included'].sum()} included); diagnostic states/UTs: {(summary['diagnostic_samples'] > 0).sum()}")
