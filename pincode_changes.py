"""Annual pincode prescribing and resistance changes, with coverage sensitivity flags."""

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/pincode_changes'
from linkage_names import join_key as norm_key
from merge_pincode import pins

KEYS = ['pincode', 'organism_clean', 'antibiotic_key']
YEARS = [2023, 2024, 2025]
ORGS = ['Escherichia coli', 'Klebsiella pneumoniae']
MERGED = ROOT / 'results/pincode_merge/prescription_resistance_pincode_monthly.csv'
RX = ROOT / 'data/prescriptions_clean.csv'

def sha(path):
    with path.open('rb') as f:
        return hashlib.sha256(f.read()).hexdigest()

def total_summary(d, label, previous_year_filters=False):
    e = d[d.endpoint_rates_available]
    r = d[d.resistance_endpoints_available]
    p = d[d.prescription_endpoints_available]
    q = d[d.all_three_years_available & (d.previous_coverage_eligible_2024 if previous_year_filters else True)]
    def nn(frame, flag):
        v = frame[frame[flag]]
        return len(v), v.pincode.nunique()
    row = {'view': label, 'candidate_comparisons': len(d), 'candidate_pincodes': d.pincode.nunique(),
           'resistance_evaluable_comparisons': len(r), 'resistance_evaluable_pincodes': r.pincode.nunique(),
           'prescription_evaluable_comparisons': len(p), 'prescription_evaluable_pincodes': p.pincode.nunique(),
           'joint_evaluable_comparisons': len(e), 'joint_evaluable_pincodes': e.pincode.nunique(),
           'joint_evaluable_positive_same_drug_both_years': int(e.same_drug_positive_both_years.sum()),
           'three_year_evaluable_comparisons': len(q), 'three_year_evaluable_pincodes': q.pincode.nunique()}
    for flag in ['resistance_increased', 'share_increased', 'count_increased',
                 'both_share_resistance_increased', 'both_count_resistance_increased',
                 'all_three_increased', 'both_share_resistance_increased_each_year']:
        source = q if flag.endswith('each_year') else e
        row[flag + '_comparisons'], row[flag + '_pincodes'] = nn(source, flag)
    row['resistance_increased_among_resistance_evaluable'] = int(r.resistance_increased.sum())
    row['both_share_resistance_increased_positive_same_drug_both_years'] = int((e.both_share_resistance_increased & e.same_drug_positive_both_years).sum())
    return row

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in [MERGED, RX]}
    raw = pd.read_csv(MERGED, dtype={'pincode': str}, low_memory=False)
    assert not raw.duplicated(KEYS + ['month_year']).any()
    d = raw[raw.organism_clean.isin(ORGS)].copy()
    d['year'] = d.month_year.str[:4].astype(int)
    d['month'] = d.month_year.str[5:7].astype(int)
    d['quarter'] = (d.month - 1) // 3 + 1
    flags = d.groupby(KEYS).agg(state_label_conflict=('pincode_state_conflict','max'),
                                recorded_states=('diagnostic_states','first'))
    index = pd.MultiIndex.from_frame(d[KEYS].drop_duplicates().merge(pd.DataFrame({'year': YEARS}), how='cross'))
    annual = d.groupby(KEYS + ['year']).agg(resistant=('RESISTANT','sum'), tested=('total_tests','sum'),
        ast_months=('month','nunique'), ast_quarters=('quarter','nunique')).reindex(index).reset_index()
    for col in ['resistant', 'tested', 'ast_months', 'ast_quarters']:
        annual[col] = annual[col].fillna(0).astype('int64')
    assert annual.tested.sum() == d.total_tests.sum()
    assert annual.resistant.sum() == d.RESISTANT.sum()
    annual['resistance_pct'] = 100 * annual.resistant / annual.tested.replace(0, np.nan)

    rx = pd.read_csv(RX, usecols=['patient_rx_pincode', 'month_year', 'Rx_Count', 'product_clean', 'state_included'],
                     dtype={'patient_rx_pincode':'string'}, low_memory=False)
    in_period = rx.month_year.between('2023-01','2025-12')
    state_included = rx.state_included.astype(str).str.lower().eq('true')
    upstream_excluded_state_records = int(rx.loc[in_period & ~state_included,'Rx_Count'].sum())
    rx = rx[in_period & state_included].copy()
    rx['pincode'] = pins(rx.patient_rx_pincode)
    rx = rx[rx.pincode.isin(d.pincode.unique()) & rx.Rx_Count.gt(0)].copy()
    rx['antibiotic_key'] = rx.product_clean.map(norm_key)
    rx['year'] = rx.month_year.str[:4].astype(int)
    den = rx.groupby(['pincode','year'], as_index=False).agg(total_ab_rx=('Rx_Count','sum'), rx_months=('month_year','nunique'))
    num = rx.groupby(['pincode','year','antibiotic_key'], as_index=False).Rx_Count.sum().rename(columns={'Rx_Count':'drug_rx'})
    annual = annual.merge(den,on=['pincode','year'],how='left',validate='many_to_one').merge(
        num,on=['pincode','year','antibiotic_key'],how='left',validate='many_to_one')
    observed = annual.total_ab_rx.gt(0)
    annual.loc[observed & annual.drug_rx.isna(),'drug_rx'] = 0
    annual['prescription_share_pct'] = 100 * annual.drug_rx / annual.total_ab_rx
    annual['rates_available'] = annual.tested.gt(0) & observed
    columns = ['resistant','tested','ast_months','ast_quarters','resistance_pct','drug_rx',
               'total_ab_rx','rx_months','prescription_share_pct','rates_available']
    w = annual.set_index(KEYS+['year'])[columns].unstack('year')
    w.columns = [f'{col}_{year}' for col,year in w.columns]
    w = w.join(flags).reset_index()
    w['previous_intrinsic_combination_exclusion'] = w.organism_clean.eq('Klebsiella pneumoniae') & w.antibiotic_key.eq('ampicillin')
    w['previous_geographic_and_combination_scope'] = ~w.state_label_conflict & ~w.previous_intrinsic_combination_exclusion
    w['resistance_endpoints_available'] = w.tested_2023.gt(0) & w.tested_2025.gt(0)
    w['prescription_endpoints_available'] = w.total_ab_rx_2023.gt(0) & w.total_ab_rx_2025.gt(0)
    w['endpoint_rates_available'] = w.resistance_endpoints_available & w.prescription_endpoints_available
    w['same_drug_positive_both_years'] = w.drug_rx_2023.gt(0) & w.drug_rx_2025.gt(0)
    w['all_three_years_available'] = w.endpoint_rates_available & w.rates_available_2024

    for before, after in [(2023,2025),(2023,2024),(2024,2025)]:
        avail = w[f'rates_available_{before}'] & w[f'rates_available_{after}']
        w[f'res_up_{before}_{after}'] = w[f'resistant_{after}']*w[f'tested_{before}'] > w[f'resistant_{before}']*w[f'tested_{after}']
        w[f'share_up_{before}_{after}'] = w[f'drug_rx_{after}']*w[f'total_ab_rx_{before}'] > w[f'drug_rx_{before}']*w[f'total_ab_rx_{after}']
        w[f'joint_up_{before}_{after}'] = avail & w[f'res_up_{before}_{after}'] & w[f'share_up_{before}_{after}']
    w['resistance_increased'] = w.resistance_endpoints_available & w.res_up_2023_2025
    w['share_increased'] = w.prescription_endpoints_available & w.share_up_2023_2025
    w['count_increased'] = w.prescription_endpoints_available & w.drug_rx_2025.gt(w.drug_rx_2023)
    w['both_share_resistance_increased'] = w.joint_up_2023_2025
    w['both_count_resistance_increased'] = w.endpoint_rates_available & w.resistance_increased & w.count_increased
    w['all_three_increased'] = w.both_share_resistance_increased & w.count_increased
    w['both_share_resistance_increased_each_year'] = w.all_three_years_available & w.joint_up_2023_2024 & w.joint_up_2024_2025
    w['delta_resistance_pp'] = w.resistance_pct_2025-w.resistance_pct_2023
    w['delta_prescription_share_pp'] = w.prescription_share_pct_2025-w.prescription_share_pct_2023
    rules = [('tested',30),('ast_months',6),('ast_quarters',4),('rx_months',12),('total_ab_rx',100),('drug_rx',10)]
    failcols = []
    for col, cutoff in rules:
        fail = f'failed_previous_{col}_minimum_{cutoff}'
        failcols.append(fail)
        w[fail] = ~(w[f'{col}_2023'].ge(cutoff) & w[f'{col}_2025'].ge(cutoff))
    w['previous_coverage_eligible'] = ~w[failcols].any(axis=1)
    w['previous_coverage_eligible_2024'] = np.logical_and.reduce([w[f'{col}_2024'].ge(cutoff) for col,cutoff in rules])
    w['previous_all_criteria_eligible'] = w.previous_geographic_and_combination_scope & w.previous_coverage_eligible
    w['unavailable_reason'] = [
        '; '.join([label for label,missing in [
            ('no AST denominator in 2023', n23<=0), ('no AST denominator in 2025', n25<=0),
            ('no prescription denominator in 2023', pd.isna(p23) or p23<=0),
            ('no prescription denominator in 2025', pd.isna(p25) or p25<=0)] if missing])
        for n23,n25,p23,p25 in zip(w.tested_2023,w.tested_2025,w.total_ab_rx_2023,w.total_ab_rx_2025)]
    w['previous_exclusion_reasons'] = [
        '; '.join(([c for c in failcols if row[c]]) +
                  (['state_label_conflict'] if row.state_label_conflict else []) +
                  (['K_pneumoniae_ampicillin'] if row.previous_intrinsic_combination_exclusion else []))
        for _,row in w.iterrows()]

    views = {'no_coverage_or_additional_exclusions': w,
             'no_coverage_thresholds_previous_geographic_combination_scope': w[w.previous_geographic_and_combination_scope],
             'previous_all_filters': w[w.previous_all_criteria_eligible]}
    summaries, by_drug = [], []
    for label, frame in views.items():
        previous_year_filters = label=='previous_all_filters'
        summaries.append(total_summary(frame,label+' | both organisms',previous_year_filters))
        for org in ORGS:
            subset=frame[frame.organism_clean.eq(org)]
            summaries.append(total_summary(subset,label+' | '+org,previous_year_filters))
            for drug,g in subset.groupby('antibiotic_key'):
                by_drug.append({'organism':org,'antibiotic':drug,**total_summary(g,label,previous_year_filters)})
    summary = pd.DataFrame(summaries)
    drugs = pd.DataFrame(by_drug)


    filter_audit = []
    for label,base in [('both organisms',w[w.previous_geographic_and_combination_scope & w.endpoint_rates_available]),
                       ('E_coli_amoxicillin_clavulanic_acid',w[w.previous_geographic_and_combination_scope & w.endpoint_rates_available & w.organism_clean.eq(ORGS[0]) & w.antibiotic_key.eq('amoxicillin/clavulanic acid')])]:
        remaining = pd.Series(True,index=base.index)
        for i,fail in enumerate(failcols):
            only=base[fail] & ~base[[c for c in failcols if c!=fail]].any(axis=1)
            newly=remaining & base[fail]
            remaining &= ~base[fail]
            filter_audit.append({'scope':label,'rule':fail,'base_evaluable_comparisons':len(base),
                'base_pincodes':base.pincode.nunique(),'fails_rule_alone_comparisons':int(base[fail].sum()),
                'affected_pincodes_rule_alone':base.loc[base[fail],'pincode'].nunique(),
                'co_increases_failing_rule_alone':int((base[fail]&base.both_share_resistance_increased).sum()),
                'fails_only_this_rule_comparisons':int(only.sum()),'cumulative_step':i+1,
                'newly_excluded_at_step':int(newly.sum()),'remaining_after_step':int(remaining.sum()),
                'remaining_pincodes_after_step':base.loc[remaining,'pincode'].nunique()})
    audit_table=pd.DataFrame(filter_audit)
    oldscope=views['no_coverage_thresholds_previous_geographic_combination_scope']
    excluded=oldscope[oldscope.endpoint_rates_available & ~oldscope.previous_coverage_eligible]
    conflict_summary=w[w.state_label_conflict].groupby('pincode',as_index=False).agg(
        candidate_comparisons=('antibiotic_key','size'),evaluable_comparisons=('endpoint_rates_available','sum'),
        co_increases=('both_share_resistance_increased','sum'))
    directions = []
    for label, frame in [('Both organisms', w), *[(org, w[w.organism_clean.eq(org)]) for org in ORGS]]:
        e = frame[frame.endpoint_rates_available].copy()
        res_sign = np.sign(e.resistant_2025 * e.tested_2023 - e.resistant_2023 * e.tested_2025)
        rx_sign = np.sign(e.drug_rx_2025 * e.total_ab_rx_2023 - e.drug_rx_2023 * e.total_ab_rx_2025)
        groups = {
            'both_increased': (res_sign > 0) & (rx_sign > 0),
            'both_decreased': (res_sign < 0) & (rx_sign < 0),
            'opposite_directions': res_sign * rx_sign < 0,
            'at_least_one_unchanged': (res_sign == 0) | (rx_sign == 0),
        }
        assert sum(int(mask.sum()) for mask in groups.values()) == len(e)
        for direction, mask in groups.items():
            directions.append({'organism': label, 'direction': direction,
                'comparisons': int(mask.sum()), 'pincodes': int(e.loc[mask, 'pincode'].nunique()),
                'antibiotics': int(e.loc[mask, 'antibiotic_key'].nunique())})
    pd.DataFrame(directions).to_csv(OUT / 'change_directions.csv', index=False)
    annual.to_csv(OUT / 'annual_pincode_counts.csv', index=False)

    files={'all_candidate_comparisons_and_flags':w,'counts_by_view':summary,'counts_by_organism_and_drug':drugs,
           'previous_filter_impact':audit_table,'comparisons_previously_excluded_by_coverage':excluded,
           'co_increases_no_coverage_filters':w[w.both_share_resistance_increased],
           'amoxicillin_clavulanic_acid_ecoli_all_pincodes':w[w.organism_clean.eq(ORGS[0]) & w.antibiotic_key.eq('amoxicillin/clavulanic acid')],
           'previous_state_conflict_exclusions':conflict_summary}
    for name,frame in files.items():
        frame.to_csv(OUT/(name+'.csv'),index=False)
    metadata={'inputs':inputs,'study_period':YEARS,'organisms':ORGS,
        'full_merge_pincodes_all_organisms':int(raw.pincode.nunique()),
        'two_organism_pincodes_before_additional_exclusions':int(d.pincode.nunique()),
        'upstream_scope':'Urine AST with valid six-digit pincodes and whole-product matches; prescriptions from the 25 eligible states/UTs.',
        'upstream_prescription_records_excluded_by_state_flag':upstream_excluded_state_records,
        'restored_state_conflict_pincodes':sorted(w.loc[w.state_label_conflict,'pincode'].unique()),
        'restored_K_pneumoniae_ampicillin_candidates':int(w.previous_intrinsic_combination_exclusion.sum()),
        'restored_K_pneumoniae_ampicillin_evaluable':int((w.previous_intrinsic_combination_exclusion&w.endpoint_rates_available).sum()),
        'restored_K_pneumoniae_ampicillin_co_increases':int((w.previous_intrinsic_combination_exclusion&w.both_share_resistance_increased).sum()),
        'missingness':'No AST denominator or no total-antibiotic prescription denominator makes the relevant rate unavailable. Missing data are not assigned zero.',
        'zero_prescriptions':'Preserves existing convention: no same-drug records within an observed prescription pincode-year is zero recorded activity, not evidence of zero use. Positive-same-drug-in-both-years is separately flagged.',
        'calculations':'Ratios of summed annual counts. Exact directional comparisons by integer cross multiplication. No models or significance tests.',
        'validation':'Annual resistant and tested counts conserved; unique pincode-organism-drug-month keys.'}
    for path,digest in inputs.items():
        assert sha(ROOT/path)==digest
    (OUT/'provenance_and_checks.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(summary.to_string(index=False))

if __name__=='__main__':
    main()
