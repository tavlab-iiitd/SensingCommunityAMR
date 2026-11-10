"""Export the geographic footprint as three separate panels."""

from pathlib import Path
import hashlib
import json
import math
import os
import re
import unicodedata

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
FIGURES = HERE / "results/figures"
OUT = FIGURES / "figure_01_separate_panels"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FuncFormatter, LogLocator

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})


ALIASES = {
    "andaman nicobar islands": "andaman and nicobar islands",
    "andaman and nicobar": "andaman and nicobar islands",
    "chattisgarh": "chhattisgarh",
    **{name: "dadra and nagar haveli and daman and diu" for name in [
        "dadra nagar haveli daman diu", "dadra nagar haveli", "dadra and nagar haveli",
        "dadra and nagar haveli daman and diu", "dadra and nagar haveli and daman diu",
        "daman diu", "daman and diu"]},
    "jammu kashmir": "jammu and kashmir", "nct delhi": "delhi",
    "nct of delhi": "delhi", "new delhi": "delhi", "orissa": "odisha",
    "pondicherry": "puducherry", "telengana": "telangana", "uttaranchal": "uttarakhand",
}


def state_key(name):
    if pd.isna(name):
        return None
    text = unicodedata.normalize("NFKD", str(name)).encode("ascii", "ignore").decode()
    text = re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()
    if text in {"", "nan", "none", "unknown", "not available"}:
        return None
    return ALIASES.get(text, text)


summary = pd.read_csv(FIGURES / "supp_fig1_geographic_footprint_data.csv")
summary["included"] = summary["included"].eq(True)
geo = json.loads((HERE / "geodata/geoBoundaries-IND-ADM1_simplified.geojson").read_text())["features"]
assert set(summary["state"]) == {state_key(f["properties"]["shapeName"]) for f in geo}


rx = pd.read_csv(DATA / "prescriptions_clean.csv",
                 usecols=["month_year", "patient_rx_state", "Rx_Count"], low_memory=False)
rx["state"] = rx["patient_rx_state"].map({s: state_key(s) for s in rx["patient_rx_state"].unique()})
rx = rx[rx["state"].notna()]
dx = pd.read_csv(DATA / "AMR_Pan_India_Aggregated_Diagnostic.csv", encoding="ISO-8859-1",
                 usecols=["end_state_date", "order_state", "Total_Test"], low_memory=False)
dx["month_year"] = pd.to_datetime(dx["end_state_date"]).dt.strftime("%Y-%m")
dx = dx[dx["month_year"].between("2023-01", "2025-12")].copy()
dx["state"] = dx["order_state"].map({s: state_key(s) for s in dx["order_state"].unique()})
dx = dx[dx["state"].notna()]
months = sorted(set(rx["month_year"]) | set(dx["month_year"]))
active = summary[(summary["prescriptions"] > 0) | (summary["diagnostic_samples"] > 0)]
rx_active = set(zip(rx.loc[rx["Rx_Count"] > 0, "state"], rx.loc[rx["Rx_Count"] > 0, "month_year"]))
dx_active = set(zip(dx.loc[dx["Total_Test"] > 0, "state"], dx.loc[dx["Total_Test"] > 0, "month_year"]))
presence = np.array([[int((s, m) in rx_active) + 2 * int((s, m) in dx_active)
                      for m in months] for s in active["state"]])


def polygons(geometry):
    shapes = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]
    return [np.asarray(shape[0]) for shape in shapes]


def save_panel(fig, name):

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"{name}.{extension}", dpi=600, bbox_inches="tight",
                    pad_inches=.12, facecolor="white")
    plt.close(fig)


def map_panel(column, cmap_name, colorbar_label, filename, outline_included=False):
    values = summary.set_index("state")[column]
    included = summary.set_index("state")["included"]
    norm = LogNorm(vmin=max(1.0, values[values > 0].min()), vmax=values.max())
    cmap = plt.get_cmap(cmap_name)
    fig = plt.figure(figsize=(6.8, 8.2))
    ax = fig.add_axes([.025, .19, .95, .80])
    for feature in geo:
        key = state_key(feature["properties"]["shapeName"])
        bold = outline_included and included[key]
        for ring in polygons(feature["geometry"]):
            ax.add_patch(Polygon(ring, closed=True,
                facecolor=cmap(norm(values[key])) if values[key] > 0 else "#eeeeee",
                edgecolor="#202020" if bold else "#777777", linewidth=.75 if bold else .28))
    xy = np.vstack([ring for f in geo for ring in polygons(f["geometry"])])
    (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
    ax.set_xlim(x0 - .02 * (x1 - x0), x1 + .02 * (x1 - x0))
    ax.set_ylim(y0 - .02 * (y1 - y0), y1 + .02 * (y1 - y0))
    ax.set_aspect(1 / math.cos(math.radians((y0 + y1) / 2)))
    ax.axis("off")


    colorbar_axis = fig.add_axes([.17, .115, .66, .023])
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_axis, orientation="horizontal")
    bar.set_label(colorbar_label, fontsize=11, labelpad=7)
    bar.locator = LogLocator(base=10, subs=(1.0,))
    bar.formatter = FuncFormatter(lambda v, _: f"{v / 1e6:g}M" if v >= 1e6 else
                                 f"{v / 1e3:g}k" if v >= 1e3 else f"{v:g}")
    bar.update_ticks()
    bar.ax.tick_params(labelsize=10, length=3)
    handles = [Patch(facecolor="#eeeeee", edgecolor="#777777", linewidth=.5,
                     label="No recorded activity")]
    if outline_included:
        handles.append(Patch(facecolor="none", edgecolor="#202020", linewidth=1.2,
                             label="Included in state analysis"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .012),
               ncol=len(handles), frameon=False, fontsize=10, handlelength=1.5, columnspacing=1.5)
    save_panel(fig, filename)


map_panel("prescriptions", "Blues", "Antibiotic prescriptions (log scale)",
          "Figure_1a_prescription_footprint", outline_included=True)
map_panel("diagnostic_samples", "Oranges", "Diagnostic culture samples (log scale)",
          "Figure_1b_diagnostic_footprint")


fig, ax = plt.subplots(figsize=(14, 10))
colors = ListedColormap(["#f0f0f0", "#4c78a8", "#f58518", "#54a24b"])
norm = BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], 4)
ax.pcolormesh(np.arange(len(months) + 1) - .5, np.arange(len(active) + 1) - .5,
              presence, cmap=colors, norm=norm, shading="flat", edgecolors="none")
ax.set_xlim(-.5, len(months) - .5)
ax.set_ylim(len(active) - .5, -.5)
ax.set_yticks(range(len(active)))
ax.set_yticklabels([f"● {name}" if inc else f"  {name}"
                    for name, inc in zip(active["state_name"], active["included"])], fontsize=10)
ax.set_xticks(range(0, len(months), 3))
ax.set_xticklabels(months[::3], rotation=45, ha="right", fontsize=11)
for i, month in enumerate(months):
    if month.endswith("-01") and i > 0:
        ax.axvline(i - .5, color="white", linewidth=1.2)
ax.tick_params(length=0)
ax.set_xlabel("Calendar month", fontsize=13, labelpad=9)
ax.set_ylabel("State/UT", fontsize=13, labelpad=12)
handles = [Patch(facecolor="#f0f0f0", edgecolor="#cccccc", label="No recorded activity"),
           Patch(facecolor="#4c78a8", label="Prescription only"),
           Patch(facecolor="#f58518", label="Diagnostic only"),
           Patch(facecolor="#54a24b", label="Both streams"),
           Line2D([], [], color="black", marker="o", linestyle="none", markersize=6,
                  label="Included in state analysis")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.63, .018),
           ncol=3, frameon=False, fontsize=11, columnspacing=1.7, handlelength=1.6)
fig.subplots_adjust(left=.29, right=.985, top=.99, bottom=.23)
save_panel(fig, "Figure_1c_monthly_presence")


summary.to_csv(OUT / "Figure_1ab_source_data.csv", index=False)
heatmap_data = pd.DataFrame([(s, m, int(presence[i, j]))
    for i, s in enumerate(active["state"]) for j, m in enumerate(months)],
    columns=["state_key", "month", "coverage_code"])
heatmap_data.to_csv(OUT / "Figure_1c_source_data.csv", index=False)
metadata = {"map_data": "Exact verified Figure 1 state totals",
            "heatmap_cells": int(presence.size),
            "state_count": len(active), "month_count": len(months),
            "map_color_scales": {column: {"minimum": float(max(1, summary.loc[summary[column] > 0, column].min())),
                                          "maximum": float(summary[column].max())}
                                 for column in ["prescriptions", "diagnostic_samples"]},
            "presentation": "Separate panels; no titles, subtitles or explanatory annotations; axis labels and legends retained.",
            "source_summary_sha256": hashlib.sha256((FIGURES / "supp_fig1_geographic_footprint_data.csv").read_bytes()).hexdigest()}
(OUT / "validation.json").write_text(json.dumps(metadata, indent=2) + "\n")
print("Saved three panels as PNG, PDF and SVG in:", OUT)
