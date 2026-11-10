"""NCDC, ICMR and Tata 1mg geographic footprints in 2024."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import unicodedata
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Polygon


HERE = Path(__file__).resolve().parents[1]
NO_ACTIVITY = "#E8EBED"
PANEL_INK = "#202A32"
PANELS = (
    ("a", "ncdc_reporting_labs", "#3E648B"),
    ("b", "icmr_listed_institutions", "#7A6A9E"),
    ("c", "prescription_records", "#B66D24"),
    ("d", "diagnostic_samples", "#167B72"),
)

COLOUR_SCALES = {
    "prescription_records": (10, 1_000_000, ("#F4E7D5", "#B66D24", "#6E390C")),
    "diagnostic_samples": (10, 10_000, ("#DDEEEA", "#167B72", "#084A46")),
}

SMALL_FEATURES = {
    "delhi", "chandigarh", "puducherry", "lakshadweep",
    "dadra and nagar haveli and daman and diu",
}


def state_key(name):
    """Match accented GeoJSON names to the harmonised source CSV."""
    name = unicodedata.normalize("NFKD", str(name)).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", " ", name.lower()).strip()


def exterior_rings(geometry):
    """Yield exterior polygon rings, including offshore islands."""
    if geometry["type"] == "Polygon":
        polygons = [geometry["coordinates"]]
    elif geometry["type"] == "MultiPolygon":
        polygons = geometry["coordinates"]
    else:
        raise ValueError(f"Unsupported geometry: {geometry['type']}")
    for polygon in polygons:
        points = np.asarray(polygon[0], dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
            raise ValueError("Invalid boundary coordinates")
        yield points


def load_inputs(csv_path, boundary_path):
    columns = [panel[1] for panel in PANELS]
    values = {}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = {"state_key", *columns} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing CSV columns: {sorted(missing)}")
        for row in reader:
            key = row["state_key"]
            if key in values:
                raise ValueError(f"Duplicate state in CSV: {key}")
            values[key] = {column: int(row[column]) for column in columns}
            if any(count < 0 for count in values[key].values()):
                raise ValueError(f"Negative counts for {key}")

    geojson = json.loads(boundary_path.read_text(encoding="utf-8"))
    features = []
    boundary_keys = set()
    for feature in geojson["features"]:
        key = state_key(feature["properties"]["shapeName"])
        if key in boundary_keys:
            raise ValueError(f"Duplicate boundary: {key}")
        boundary_keys.add(key)
        features.append((key, list(exterior_rings(feature["geometry"]))))
    missing = set(values) - boundary_keys
    if missing:
        raise ValueError(f"States without matching boundaries: {sorted(missing)}")
    return values, features


def ring_area(points):
    x, y = points[:, 0], points[:, 1]
    return abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def draw_panel(ax, letter, column, colour, values, features):
    scale = COLOUR_SCALES.get(column)
    if scale:
        lower, upper, colours = scale
        norm = LogNorm(vmin=lower, vmax=upper, clip=True)
        cmap = LinearSegmentedColormap.from_list(column, colours)

    for key, boundary_rings in features:
        count = values.get(key, {}).get(column, 0)
        fill = NO_ACTIVITY if count == 0 else cmap(norm(count)) if scale else colour
        for points in boundary_rings:
            ax.add_patch(Polygon(points, closed=True, facecolor=fill,
                                 edgecolor="white", linewidth=0.38, zorder=1))
        if key in SMALL_FEATURES:
            centre = max(boundary_rings, key=ring_area).mean(axis=0)
            ax.scatter(*centre, s=12, facecolor=fill, edgecolor="white",
                       linewidth=0.35, zorder=3)


    ax.set_xlim(67.0, 98.5)
    ax.set_ylim(6.0, 37.7)
    ax.set_aspect(1 / math.cos(math.radians(22)))
    ax.set_axis_off()
    ax.text(0.0, 0.995, letter, transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold", color=PANEL_INK)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "results/surveillance/source_data/geographic_footprints_2024.csv")
    parser.add_argument("--boundaries", type=Path, default=HERE / "geodata/geoBoundaries-IND-ADM1_simplified.geojson")
    parser.add_argument("--output-stem", type=Path,
                        default=HERE / "results/surveillance/Figure_surveillance_four_maps_2024")
    parser.add_argument("--width", type=float, default=7.2, help="Figure width in inches")
    parser.add_argument("--height", type=float, default=7.8, help="Figure height in inches")
    parser.add_argument("--dpi", type=int, default=600, help="PNG and TIFF resolution")
    parser.add_argument("--font", default="Arial", help="Panel-letter font")
    args = parser.parse_args()
    if args.width <= 0 or args.height <= 0 or args.dpi <= 0:
        parser.error("Width, height and dpi must be positive")

    values, features = load_inputs(args.data, args.boundaries)
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [args.font, "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 2, figsize=(args.width, args.height), facecolor="white")
    fig.subplots_adjust(left=0.025, right=0.985, bottom=0.015, top=0.985,
                        wspace=0.04, hspace=0.04)
    for ax, panel in zip(axes.flat, PANELS):
        draw_panel(ax, *panel, values, features)


    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png", "tiff"):
        output = args.output_stem.with_suffix("." + extension)
        options = {"dpi": args.dpi} if extension in {"png", "tiff"} else {}
        if extension == "tiff":
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(output, **options)
        print(output)
    plt.close(fig)


if __name__ == "__main__":
    main()
