"""Pooled urinary resistance by state; percentages hidden below 30 tests."""

from pathlib import Path
import csv
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/statewise"
OUT = BASE / "heatmap_outputs"
AST = ROOT / "data/ast_clean.csv"
SELECTION = BASE / "trend_outputs/top10_states.csv"
START, END = "2023-01", "2025-12"
MIN_TESTED = 30
COUNTS = ["RESISTANT", "INTERMEDIATE", "SENSITIVE"]


ENTEROBACTERALES = [
    "Amoxicillin/Clavulanic acid", "Ceftriaxone", "Ciprofloxacin",
    "Trimethoprim/Sulfamethoxazole", "Nitrofurantoin", "Gentamicin",
    "Amikacin", "Piperacillin/Tazobactam", "Ertapenem", "Meropenem",
]
PANELS = {
    "Escherichia coli": ENTEROBACTERALES,
    "Klebsiella pneumoniae": ENTEROBACTERALES,
    "Pseudomonas aeruginosa": [
        "Ceftazidime", "Cefepime", "Aztreonam", "Ciprofloxacin",
        "Levofloxacin", "Amikacin", "Piperacillin/Tazobactam",
        "Imipenem", "Meropenem",
    ],
}
LABELS = {
    "Amoxicillin/Clavulanic acid": "Amoxicillin/\nclavulanic acid",
    "Trimethoprim/Sulfamethoxazole": "Trimethoprim/\nsulfamethoxazole",
    "Piperacillin/Tazobactam": "Piperacillin/\ntazobactam",
}
SLUGS = {
    "Escherichia coli": "ecoli",
    "Klebsiella pneumoniae": "klebsiella",
    "Pseudomonas aeruginosa": "pseudomonas",
}
LOW_N = "#e2e5e9"
NO_TESTS = "#fafafa"


def state_name(value):
    if pd.isna(value):
        return None
    name = " ".join(str(value).strip().lower().replace("&", "and").split())
    aliases = {
        "pondicherry": "puducherry", "nct delhi": "delhi", "new delhi": "delhi",
        "orissa": "odisha", "telengana": "telangana", "chattisgarh": "chhattisgarh",
        "dadra and nagar haveli": "dadra and nagar haveli and daman and diu",
        "daman and diu": "dadra and nagar haveli and daman and diu",
    }
    if name in {"", "nan", "none", "null", "unknown", "not available", "other", "others"}:
        return None
    return aliases.get(name, name).title()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def prepare():
    selection = pd.read_csv(SELECTION).sort_values("rank")
    states = selection["state"].tolist()
    assert len(states) == 10 and len(set(states)) == 10
    assert selection["antibiotic_records"].ge(1000).all()
    assert selection["shared_months"].ge(1).all()

    ast = pd.read_csv(AST, usecols=["month_year", "order_state", "Order_Pincode",
                                  "test_name", "organism_clean", "antibiotic_clean", *COUNTS])
    ast = ast.loc[ast["test_name"].eq("URINE")
                  & ast["month_year"].between(START, END)
                  & ast["organism_clean"].isin(PANELS)].copy()
    assert ast[COUNTS].notna().all().all()
    assert ast[COUNTS].ge(0).all().all()
    assert np.equal(ast[COUNTS], np.floor(ast[COUNTS])).all().all()
    ast["state"] = ast["order_state"].map(state_name)
    ast["tested"] = ast[COUNTS].sum(axis=1)
    ast = ast.loc[ast["tested"].gt(0)].copy()
    unknown_rows = int(ast["state"].isna().sum())
    ast = ast.loc[ast["state"].notna()].copy()
    pins = ast["Order_Pincode"].astype("string").str.strip().str.replace(r"\.0$", "", regex=True)
    ast["pin"] = pins.where(pins.str.fullmatch(r"[1-9]\d{5}", na=False))
    ast = ast.loc[[drug in PANELS[org] for org, drug in
                   zip(ast["organism_clean"], ast["antibiotic_clean"])]].copy()
    keys = ["state", "organism_clean", "antibiotic_clean"]
    pooled = ast.groupby(keys, as_index=False).agg(
        resistant=("RESISTANT", "sum"), intermediate=("INTERMEDIATE", "sum"),
        sensitive=("SENSITIVE", "sum"), tested=("tested", "sum"),
        months_tested=("month_year", "nunique"), first_month=("month_year", "min"),
        last_month=("month_year", "max"), valid_pincodes=("pin", "nunique"),
    )
    pooled["resistance_pct"] = 100 * pooled["resistant"] / pooled["tested"]
    assert pooled["resistance_pct"].between(0, 100).all()
    assert pooled["tested"].eq(pooled[["resistant", "intermediate", "sensitive"]].sum(axis=1)).all()
    assert not pooled.duplicated(keys).any()

    grid = pd.DataFrame([(state, org, drug) for org, drugs in PANELS.items()
                         for state in states for drug in drugs], columns=keys)
    shown = grid.merge(pooled, on=keys, how="left", validate="one_to_one")
    count_fields = ["resistant", "intermediate", "sensitive", "tested", "months_tested", "valid_pincodes"]
    shown[count_fields] = shown[count_fields].fillna(0).astype(int)
    shown["display_status"] = np.select(
        [shown["tested"].eq(0), shown["tested"].lt(MIN_TESTED)],
        ["No tests", "Fewer than 30 tests"], default="Displayed",
    )
    shown["display_resistance_pct"] = shown["resistance_pct"].where(shown["tested"].ge(MIN_TESTED))
    shown = shown.merge(selection[["state", "rank", "diagnostic_samples"]], on="state", validate="many_to_one")
    return states, selection, pooled, shown, unknown_rows


def verify(shown):
    """Reconcile all plotted R/I/S counts by streaming the source independently."""
    keys = ["state", "organism_clean", "antibiotic_clean"]
    expected = {tuple(row): [0, 0, 0] for row in shown[keys].itertuples(index=False, name=None)}
    with AST.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["test_name"] != "URINE" or not START <= row["month_year"] <= END:
                continue
            key = (state_name(row["order_state"]), row["organism_clean"], row["antibiotic_clean"])
            if key in expected:
                for i, col in enumerate(COUNTS):
                    expected[key][i] += int(row[col])
    for row in shown.itertuples(index=False):
        key = (row.state, row.organism_clean, row.antibiotic_clean)
        assert expected[key] == [row.resistant, row.intermediate, row.sensitive], key
        assert row.tested == sum(expected[key]), key
    assert shown.loc[shown["tested"].lt(MIN_TESTED), "display_resistance_pct"].isna().all()
    assert shown.loc[shown["tested"].eq(0), "resistance_pct"].isna().all()


    historical = pd.read_csv(BASE / "results/urinary_ast_coverage_by_state_organism_agent.csv")
    matched = shown.merge(historical, on=keys, suffixes=("", "_prior"), validate="one_to_one")
    for col in ["tested", "resistant", "intermediate", "months_tested"]:
        np.testing.assert_array_equal(matched[col], matched[col + "_prior"])
    prior_monthly = pd.read_csv(BASE / "trend_outputs/resistance_monthly.csv")
    prior_monthly = prior_monthly.rename(columns={"organism": "organism_clean", "series": "antibiotic_clean"})
    prior_totals = prior_monthly.groupby(keys, as_index=False)[["R", "N"]].sum()
    joined = shown.merge(prior_totals, on=keys, validate="one_to_one")
    np.testing.assert_array_equal(joined["resistant"], joined["R"])
    np.testing.assert_array_equal(joined["tested"], joined["N"])
    return {
        "status": "PASS", "cells_checked_directly_against_source": len(shown),
        "cells_matching_existing_coverage_audit": len(matched),
        "cells_matching_existing_monthly_outputs": len(joined),
        "display_status_counts": shown["display_status"].value_counts().to_dict(),
        "suppressed_percentages_are_missing": True,
        "no_tests_are_distinct_from_zero_resistance": True,
    }


def draw_panel(organism, states, shown):
    drugs = PANELS[organism]
    panel = shown.loc[shown["organism_clean"].eq(organism)].set_index(["state", "antibiotic_clean"])
    grid = pd.MultiIndex.from_product([states, drugs], names=["state", "antibiotic_clean"])
    panel = panel.reindex(grid)
    values = panel["display_resistance_pct"].to_numpy().reshape(len(states), len(drugs))
    numbers = panel["tested"].to_numpy().reshape(values.shape)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad(LOW_N)
    fig = plt.figure(figsize=(13.6, 8.8), facecolor="white")
    ax = fig.add_axes([0.135, 0.245, 0.765, 0.565])
    im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(np.arange(len(drugs)))
    ax.set_xticklabels([LABELS.get(d, d) for d in drugs], rotation=40, ha="right", fontsize=10)
    ax.set_yticks(np.arange(len(states)))
    ax.set_yticklabels(states, fontsize=11)
    ax.tick_params(axis="both", which="both", length=0, pad=9)
    ax.set_xticks(np.arange(-0.5, len(drugs), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(states), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i in range(len(states)):
        for j in range(len(drugs)):
            n = int(numbers[i, j])
            if n == 0:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=NO_TESTS,
                                       edgecolor="#d4d8dd", linewidth=0.4, hatch="///"))
                ax.text(j, i, "No tests", ha="center", va="center", fontsize=9, color="#69717a")
            elif n < MIN_TESTED:
                ax.text(j, i, f"n={n:,}", ha="center", va="center", fontsize=10, color="#515b65")
            else:
                value = values[i, j]
                color = "white" if value >= 60 else "#20262c"
                rate_label = "<0.1%" if 0 < value < 0.1 else f"{value:.1f}%"
                ax.text(j, i - 0.13, rate_label, ha="center", va="center",
                        fontsize=11, fontweight="semibold", color=color)
                ax.text(j, i + 0.20, f"n={n:,}", ha="center", va="center", fontsize=8.6, color=color)

    cbax = fig.add_axes([0.925, 0.30, 0.016, 0.46])
    colorbar = fig.colorbar(im, cax=cbax, ticks=[0, 20, 40, 60, 80, 100])
    colorbar.set_label("Resistant (%)", fontsize=11, labelpad=10)
    colorbar.ax.tick_params(labelsize=10, length=3)
    colorbar.outline.set_visible(False)
    fig.text(0.135, 0.945, "State-wise urinary resistance, 2023–2025", fontsize=20, fontweight="semibold")
    fig.text(0.135, 0.899, organism, fontsize=17, fontstyle="italic")
    fig.text(0.135, 0.856, "Pooled resistance percentage and number tested (n)", fontsize=11, color="#505b65")
    fig.legend(handles=[Patch(facecolor=LOW_N, label=f"Fewer than {MIN_TESTED} tests: percentage hidden"),
                        Patch(facecolor=NO_TESTS, edgecolor="#d4d8dd", hatch="///", label="No tests")],
               loc="lower left", bbox_to_anchor=(0.13, 0.078), frameon=False, ncol=2, fontsize=10)
    fig.text(0.135, 0.055, "Resistance = 100 × ΣR / Σ(R + I + S). Counts are summed across observed months.",
             fontsize=10, color="#444e58")
    fig.text(0.135, 0.027, "Same ten states, ordered by diagnostic volume. Descriptive, unadjusted platform estimates.",
             fontsize=10, color="#444e58")
    return fig


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "axes.unicode_minus": False})
    states, selection, pooled, shown, unknown_rows = prepare()
    validation = verify(shown)
    pooled.to_csv(OUT / "pooled_resistance_selected_agents_all_states.csv", index=False)
    shown.to_csv(OUT / "heatmap_source_data.csv", index=False)
    selection.to_csv(OUT / "state_selection.csv", index=False)
    pdf_path = OUT / "statewise_pooled_resistance_heatmaps.pdf"
    with PdfPages(pdf_path, metadata={"Title": "State-wise pooled urinary resistance, 2023–2025",
                                     "Subject": "Three organism-specific descriptive heatmaps"}) as pdf:
        for organism in PANELS:
            fig = draw_panel(organism, states, shown)
            stem = "statewise_pooled_resistance_" + SLUGS[organism]
            fig.savefig(OUT / f"{stem}.png", dpi=220, facecolor="white")
            fig.savefig(OUT / f"{stem}.svg", facecolor="white")
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)
    manifest = {
        "period": [START, END], "specimen": "URINE", "min_tested_for_display": MIN_TESTED,
        "formula": "100 * sum(RESISTANT) / sum(RESISTANT + INTERMEDIATE + SENSITIVE)",
        "state_order": states, "drug_panels": PANELS,
        "state_selection": "Existing top ten by all diagnostic-sample volume, >=1000 antibiotic records and >=1 shared urinary-AST month",
        "unknown_state_source_rows_excluded": unknown_rows,
        "source_files": {str(p.relative_to(ROOT)): digest(p) for p in [AST, SELECTION, Path(__file__)]},
        "runtime": {"python": platform.python_version(), "executable": sys.executable,
                    "pandas": pd.__version__, "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "validation": validation,
    }
    (OUT / "validation_and_sources.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    print(f"Saved heatmaps, source counts and manifest to {OUT}")


if __name__ == "__main__":
    main()
