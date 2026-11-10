"""Check Not recommended assignments against the WHO workbook."""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/aware_check'
OUT.mkdir(parents=True, exist_ok=True)
ALIASES = {
    'amoxycillin': 'amoxicillin', 'tazobactum': 'tazobactam',
    'potassium clavulanate': 'clavulanic acid', 'serratiopeptidase': 'serrapeptase',
    'tobramicin': 'tobramycin', 'benzyl penicillin': 'benzylpenicillin',
    'penicillin g': 'benzylpenicillin', 'lactobacillus sporogenes': 'bacillus coagulans',
    'bromelain': 'bromelains',
}

def key(value):
    parts = re.split(r'[+/]', str(value).casefold().strip())
    return tuple(sorted(ALIASES.get(' '.join(p.split()), ' '.join(p.split()))
                        for p in parts if p.strip()))

who = pd.read_excel(ROOT / 'data/aware-2025.xlsx', sheet_name='Not recommended', header=None)
assert who.iat[3, 0] == 'Antibiotic'
entries = {}
for i in range(4, len(who)):
    if pd.notna(who.iat[i, 0]):
        name = str(who.iat[i, 0]).strip()
        entries.setdefault(key(name), []).append((i + 1, name))
rx = pd.read_csv(ROOT / 'data/prescriptions_clean.csv', low_memory=False)
lookup = {name: key(name) in entries for name in rx.salt_name.unique()}
rx['independent_nr'] = rx.salt_name.map(lookup)
assert rx.independent_nr.eq(rx.aware_2025.eq('Not recommended')).all()
products = rx.groupby(['salt_name', 'aware_2025'], as_index=False).Rx_Count.sum()
products['pct_all_antibiotics'] = 100 * products.Rx_Count / rx.Rx_Count.sum()
products.to_csv(OUT / 'all_product_assignments.csv', index=False)
nr = products[products.aware_2025.eq('Not recommended')].copy()
nr['WHO_name'] = nr.salt_name.map(lambda s: '; '.join(n for _, n in entries[key(s)]))
nr['WHO_excel_rows'] = nr.salt_name.map(lambda s: '; '.join(str(i) for i, _ in entries[key(s)]))
nr.to_csv(OUT / 'not_recommended_product_matches.csv', index=False)
rows = []
for scope, frame in [('National cohort', rx), ('Eligible states', rx[rx.state_included])]:
    numerator = int(frame.loc[frame.independent_nr, 'Rx_Count'].sum())
    denominator = int(frame.Rx_Count.sum())
    rows.append({'scope': scope, 'not_recommended_prescriptions': numerator,
                 'antibiotic_prescriptions': denominator, 'percentage': 100 * numerator / denominator})
pd.DataFrame(rows).to_csv(OUT / 'cohort_denominators.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
