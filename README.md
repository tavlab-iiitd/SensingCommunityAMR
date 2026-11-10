# SensingCommunityAMR

Code for *Antibiotic prescribing and resistance beyond hospital surveillance networks in India*,
including the Nature Health analyses through 10 November 2026.

## Run

Use Python 3.10+ and R. Place the five study inputs listed in [data/README.md](data/README.md)
in `data/`, then run from this directory:

```sh
python -m pip install -r requirements.txt
Rscript R/install_packages.R
python run_all.py
```

The full run prepares the data, fits the models and runs the additional analyses.
Individual stages can also be run:

```sh
python run_all.py --section core
python run_all.py --section models
python run_all.py --section additional
```

Run `core` first. The other two stages use its cleaned data. Scripts can also be run
individually, for example `python python/merge_pincode.py` or `Rscript R/resistance_trends.R`.
All generated tables, figures and fitted models go to `results/`; cleaned data go to `data/`.
A rerun replaces the generated outputs. Study data and generated results are excluded by `.gitignore`.

## Analyses

| Code | Output |
|---|---|
| `python/preprocess_prescriptions.py`, `python/preprocess_diagnostics.py` | Cleaned records and WHO AWaRe 2025 assignments |
| `python/descriptive_tables.py` | Cohort, diagnostic and prescribing summaries |
| `python/prescription_state_shares.py` | State Access and Watch shares, 1,000 pincode bootstrap samples |
| `python/prescription_time_trends.py` | Monthly prescription shares, both denominators |
| `python/resistance_monthly.py` | Monthly resistance, three organisms |
| `python/geographic_footprint.py`, `python/geographic_panels.py` | Geographic coverage, monthly activity and separate figure panels |
| `R/01_resistance_profiles.R` | Pooled resistance profiles |
| `R/prescription_trends.R`, `R/resistance_trends.R` | Corrected monthly effects, confidence intervals and BH-adjusted p-values |
| `python/merge_prescriptions_resistance.py`, `R/leadlag.R` | National exploratory lead–lag analyses |
| `python/merge_pincode.py` | Pincode–month–organism–antibiotic linkage and coverage |
| `python/pincode_changes.py` | Annual pincode comparisons, direction counts and coverage sensitivities |
| `python/surveillance_comparison.py` | 2024 NCDC and ICMR comparisons and network rosters |
| `python/surveillance_maps.py`, `python/four_network_maps.py` | Surveillance coverage maps |
| `python/state_coverage.py` | State, district and pincode coverage summaries |
| `python/state_trends.py` | Descriptive trends for the ten largest diagnostic states |


## Pincode analysis

`merge_pincode.py` links urine AST to prescriptions using pincode, month and the complete
normalised drug product. Prescription counts are summed to unique keys before joining;
resistance remains organism-specific. State-label conflicts are retained and flagged.
All antibiotic products contribute to the prescription denominator in the 25 eligible states/UTs.

A missing drug count is assigned zero only when other antibiotic prescriptions were recorded
in that pincode-month. No recorded prescription activity leaves exposure missing.
Resistance is `100 × R / (R + I + S)`. Prescription share is `100 × drug count / all antibiotic counts`.
Prescription counts repeat across organisms and must not be summed across organism rows.

`pincode_changes.py` pools annual counts from all observed months, then compares 2023 with
2025 and successive years. The primary view applies no additional coverage cutoffs or
state-conflict exclusions. It retains K. pneumoniae–ampicillin with an explicit flag.
The earlier coverage-filtered and organism/geography-restricted views are sensitivity outputs.
Missing denominators remain missing; no significance tests are fitted.

The main outputs in `results/pincode_changes/` are:

- `all_candidate_comparisons_and_flags.csv`: annual rates, changes, missingness and sensitivity flags.
- `counts_by_view.csv` and `counts_by_organism_and_drug.csv`: evaluable comparisons and co-increases.
- `change_directions.csv`: both increased, both decreased, opposite directions or at least one unchanged.
- `annual_pincode_counts.csv`: annual numerators and denominators.

The submitted primary analysis has 11,493 evaluable comparisons across 648 pincodes;
921 comparisons show both measures increasing. The 2023–2025 successive-year result is
78 comparisons across 63 pincodes. These are descriptive ecological comparisons.

## Trend models and figures

`R/monthly_marginal_effect.R` calculates the equally weighted mean of the conditional
probability derivative over observed state-months. It includes fitted state intercepts
and slopes and propagates their joint prediction covariance with the delta method.
The p-value tests that monthly effect. BH adjustment is within AWaRe group and denominator
for prescriptions, and within organism for resistance. The original model eligibility
rules and fallback order are retained. The primary tables end in `_monthly.csv`.

The final combined figure reproduces the saved 8 November artwork: ten E. coli drugs and
nine K. pneumoniae drugs after removing ampicillin. The separate monthly resistance figure
selects ten eligible K. pneumoniae drugs, including ceftazidime. These selections are intentional
and produce different panels; the combined-figure caption should describe nine K. pneumoniae drugs.

National lead–lag results use the national month–drug input, independently of the pincode
analysis. Missing product-month exposures retain the original missing-value convention.
Peaks are exploratory, with nominal p-values after searching lags and pairs. The figure exports
include the original constant bound and the lag-specific nominal bound.

## Surveillance and state analyses

The main 2024 comparison uses NCDC hospital-outpatient urinary E. coli resistance percentages.
NCDC drug-specific outpatient denominators were not reported and are left missing.
The additional ICMR comparisons use non-susceptibility derived from susceptible/tested counts;
this is distinct from resistance. Network estimates are not pooled or tested against each other.
The scripts retain source URLs, page references and the published precision.

State trends use the ten states with the largest diagnostic volumes among those with at
least 1,000 antibiotic prescriptions and one shared prescription/urinary-AST month.
Bootstrap intervals resample pincodes within each state-month. Pooled heatmaps suppress
percentages below 30 tests. These analyses do not fit state-specific associations.


## Software and boundaries

Reference runtime: Python 3.10, pandas 2.3.3, numpy 2.2.6, matplotlib 3.10.8, openpyxl 3.1.5;
R 4.1.3, data.table 1.15.4, zoo 1.8.12, glmmTMB 1.1.7, dplyr 1.1.4, tidyr 1.3.1,
forecast 8.21. Model fits can vary with package versions.

Map boundaries are from geoBoundaries gbOpen, CC BY 2.5 IN
