"""Link urine AST and prescriptions by pincode, month and whole-product drug key."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from linkage_names import join_key as norm_key

ROOT = Path(__file__).resolve().parents[1]
RX_PATH = ROOT / "data/prescriptions_clean.csv"
AST_PATH = ROOT / "data/ast_clean.csv"
OUT = ROOT / "results/pincode_merge"
COUNTS = ["SENSITIVE", "INTERMEDIATE", "RESISTANT"]
KEYS = ["pincode", "month_year", "antibiotic_key"]
SERIES = ["pincode", "organism_clean", "antibiotic_key"]
CELL = KEYS + ["organism_clean"]


def pins(s):
    clean = s.astype("string").str.strip().str.replace(r"\.0+$", "", regex=True)
    return clean.where(clean.str.fullmatch(r"[1-9][0-9]{5}", na=False))


def labels(s):
    return "; ".join(sorted(set(s.dropna().astype(str))))


def state_names(s):
    return s.astype("string").str.strip().str.lower().replace({
        "pondicherry": "puducherry", "orissa": "odisha",
        "daman & diu": "dadra and nagar haveli and daman and diu",
        "dadra & nagar haveli": "dadra and nagar haveli and daman and diu",
        "jammu & kashmir": "jammu and kashmir",
        "andaman & nicobar islands": "andaman and nicobar islands",
        "not available": pd.NA, "unknown": pd.NA, "nan": pd.NA,
    })


def validate_counts(d, cols):
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="raise")
        if d[c].isna().any() or (d[c] < 0).any() or (d[c] % 1 != 0).any():
            raise ValueError(f"Missing, negative, or non-integer counts: {c}")
        d[c] = d[c].astype("int64")


def longest_run(months):
    values = sorted(set(pd.Period(m, freq="M").ordinal for m in months))
    best = current = 0
    previous = None
    for value in values:
        current = current + 1 if previous is not None and value == previous + 1 else 1
        best = max(best, current)
        previous = value
    return best


def digest(path):
    with path.open("rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def main():
    audit = {"generated_at_utc": datetime.now(timezone.utc).isoformat(),
             "study_period": ["2023-01", "2025-12"],
             "join_keys": KEYS,
             "purpose": "Pincode-month-drug linkage for local descriptive comparisons.",
             "inputs": {str(p.relative_to(ROOT)): {"sha256": digest(p)}
                        for p in [RX_PATH, AST_PATH]}}
    rx = pd.read_csv(RX_PATH, usecols=["month_year", "patient_rx_pincode",
        "patient_rx_state", "product_clean", "Rx_Count", "state_included"],
        dtype={"patient_rx_pincode": "string"}, low_memory=False)
    rx = rx[rx.month_year.between("2023-01", "2025-12") &
            rx.state_included.astype(str).str.lower().eq("true")].copy()
    validate_counts(rx, ["Rx_Count"])
    audit["eligible_prescription_records_before_pin_filter"] = int(rx.Rx_Count.sum())
    rx["pincode"] = pins(rx.patient_rx_pincode)
    audit["prescription_records_excluded_invalid_pincode"] = int(rx.loc[rx.pincode.isna(), "Rx_Count"].sum())
    rx = rx[rx.pincode.notna() & rx.Rx_Count.gt(0)].copy()
    rx["antibiotic_key"] = rx.product_clean.map(norm_key)
    if rx.antibiotic_key.isna().any() or rx.antibiotic_key.eq("").any():
        raise ValueError("Prescription drug key is missing")
    rx["state"] = state_names(rx.patient_rx_state)
    rx_states = rx.groupby("pincode").state.agg(labels).rename("prescription_states")
    denominator = rx.groupby(["pincode", "month_year"], as_index=False).Rx_Count.sum().rename(
        columns={"Rx_Count": "total_ab_prescriptions"})
    products = rx.groupby(KEYS, as_index=False).Rx_Count.sum().rename(
        columns={"Rx_Count": "recorded_product_count"})
    product_names = rx.groupby("antibiotic_key").product_clean.agg(labels).rename("prescription_products")
    assert products.recorded_product_count.sum() == rx.Rx_Count.sum() == denominator.total_ab_prescriptions.sum()

    ast = pd.read_csv(AST_PATH, usecols=["month_year", "Order_Pincode", "order_state",
        "test_name", "organism_clean", "antibiotic_clean", "aware_2025", *COUNTS],
        dtype={"Order_Pincode": "string"}, low_memory=False)
    ast = ast[ast.month_year.between("2023-01", "2025-12") &
              ast.test_name.str.strip().str.upper().eq("URINE")].copy()
    ast["organism_clean"] = ast.organism_clean.astype("string").str.strip()
    ast = ast[ast.organism_clean.notna() & ~ast.organism_clean.str.lower().isin(
        ["", "no growth", "nan", "none", "unknown"])].copy()
    validate_counts(ast, COUNTS)
    ast["total_tests"] = ast[COUNTS].sum(axis=1)
    ast = ast[ast.total_tests.gt(0)].copy()
    ast["pincode"] = pins(ast.Order_Pincode)
    audit["urine_ast_results_before_pin_filter"] = int(ast.total_tests.sum())
    audit["ast_results_excluded_invalid_pincode"] = int(ast.loc[ast.pincode.isna(), "total_tests"].sum())
    ast = ast[ast.pincode.notna()].copy()
    ast["antibiotic_key"] = ast.antibiotic_clean.map(norm_key)
    common = ast.antibiotic_key.isin(products.antibiotic_key.unique())
    audit["ast_results_excluded_drug_not_in_prescription_stream"] = int(ast.loc[~common, "total_tests"].sum())
    ast = ast[common].copy()
    ast["state"] = state_names(ast.order_state)
    ast_states = ast.groupby("pincode").state.agg(labels).rename("diagnostic_states")
    names = ast.groupby("antibiotic_key").agg(
        laboratory_antibiotics=("antibiotic_clean", labels), aware_2025=("aware_2025", labels))
    resistance = ast.groupby(CELL, as_index=False)[COUNTS].sum()
    resistance["total_tests"] = resistance[COUNTS].sum(axis=1)
    resistance["pct_resistant"] = 100 * resistance.RESISTANT / resistance.total_tests


    state_map = pd.concat([rx[["pincode", "state"]], ast[["pincode", "state"]]]).drop_duplicates()
    conflicts = state_map.groupby("pincode").state.nunique().gt(1).rename("pincode_state_conflict")
    merged = (resistance.merge(products, on=KEYS, how="left", validate="many_to_one")
        .merge(denominator, on=["pincode", "month_year"], how="left", validate="many_to_one")
        .merge(names, on="antibiotic_key", how="left", validate="many_to_one")
        .merge(product_names, on="antibiotic_key", how="left", validate="many_to_one")
        .merge(rx_states, on="pincode", how="left", validate="many_to_one")
        .merge(ast_states, on="pincode", how="left", validate="many_to_one")
        .merge(conflicts, on="pincode", how="left", validate="many_to_one"))
    recorded = merged.recorded_product_count.notna()
    observed = merged.total_ab_prescriptions.notna()
    merged["match_status"] = np.select([recorded, observed],
        ["same_drug_recorded", "rx_month_observed_no_product"], default="no_rx_activity_in_pincode_month")
    merged["total_prescriptions_count"] = merged.recorded_product_count
    merged.loc[observed & ~recorded, "total_prescriptions_count"] = 0
    merged["pres_normalized"] = 100 * merged.total_prescriptions_count / merged.total_ab_prescriptions
    merged = merged.sort_values(SERIES + ["month_year"]).reset_index(drop=True)
    assert not merged.duplicated(CELL).any()
    assert len(merged) == len(resistance)
    assert merged[COUNTS].sum().equals(ast[COUNTS].sum())
    assert merged.pct_resistant.between(0, 100).all()
    assert merged.pres_normalized.dropna().between(0, 100).all()
    assert merged.loc[merged.match_status.eq("no_rx_activity_in_pincode_month"), "pres_normalized"].isna().all()

    merged["paired"] = merged.pres_normalized.notna()
    merged["positive_match"] = merged.match_status.eq("same_drug_recorded")
    merged["paired_tests_ge30"] = merged.paired & merged.total_tests.ge(30)
    coverage = merged.groupby(SERIES, as_index=False).agg(
        ast_months=("month_year", "nunique"), paired_months=("paired", "sum"),
        months_with_same_drug_recorded=("positive_match", "sum"),
        paired_months_tested_ge30=("paired_tests_ge30", "sum"),
        total_ast_tests=("total_tests", "sum"), median_tests_per_observed_month=("total_tests", "median"),
        pincode_state_conflict=("pincode_state_conflict", "max"))
    runs = merged[merged.paired].groupby(SERIES).month_year.agg(longest_run).rename("longest_consecutive_paired_months")
    coverage = coverage.merge(runs, on=SERIES, how="left", validate="one_to_one")
    coverage["longest_consecutive_paired_months"] = coverage.longest_consecutive_paired_months.fillna(0).astype(int)
    coverage["meets_24_paired_months"] = coverage.paired_months.ge(24)
    coverage["meets_24_months_with_same_drug_recorded"] = coverage.months_with_same_drug_recorded.ge(24)
    coverage = coverage.sort_values(["paired_months", "total_ast_tests"], ascending=False)
    merged = merged.drop(columns=["paired", "positive_match", "paired_tests_ge30"])

    audit.update({"merged_ast_cells": len(merged), "pincodes_with_eligible_ast": int(merged.pincode.nunique()),
        "organisms": int(merged.organism_clean.nunique()), "antibiotics": int(merged.antibiotic_key.nunique()),
        "match_status_counts": {k: int(v) for k, v in merged.match_status.value_counts().items()},
        "series_total": len(coverage), "series_ge24_paired_months": int(coverage.meets_24_paired_months.sum()),
        "series_ge24_positive_rx_months": int(coverage.meets_24_months_with_same_drug_recorded.sum()),
        "median_paired_months": float(coverage.paired_months.median()),
        "pincodes_with_conflicting_state_labels": int(merged.loc[merged.pincode_state_conflict, "pincode"].nunique()),
        "validation": "Counts conserved; unique AST cell keys; many-to-one joins; percentages bounded; missing exposure preserved."})
    audit["main_organisms"] = {}
    for org in ["Escherichia coli", "Klebsiella pneumoniae"]:
        d = coverage[coverage.organism_clean.eq(org)]
        audit["main_organisms"][org] = {"series": len(d), "series_ge24_paired_months": int(d.meets_24_paired_months.sum()),
            "median_paired_months": float(d.paired_months.median()) if len(d) else None,
            "series_ge24_months_tested_ge30": int(d.paired_months_tested_ge30.ge(24).sum())}

    OUT.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT / "prescription_resistance_pincode_monthly.csv", index=False, float_format="%.8g")
    coverage.to_csv(OUT / "pincode_longitudinal_coverage.csv", index=False)
    (OUT / "merge_summary.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: v for k, v in audit.items() if k != "inputs"}, indent=2))
    print(f"Outputs: {OUT}")


if __name__ == "__main__":
    main()
