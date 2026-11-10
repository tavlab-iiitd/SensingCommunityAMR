"""State, district and pincode coverage of prescriptions and diagnostic samples."""

from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/statewise'
for folder in ['results', 'figures']:
    (OUT / folder).mkdir(exist_ok=True, parents=True)
SOURCE = ROOT / 'python/name_cleaning.py'
import name_cleaning as cleaning
INPUTS = {
    'prescriptions': ROOT / 'data/prescriptions_clean.csv',
    'diagnostic_samples': ROOT / 'data/AMR_Pan_India_Aggregated_Diagnostic.csv',
    'ast': ROOT / 'data/ast_clean.csv',
    'name_cleaning': SOURCE,
}
MONTHS = pd.period_range('2023-01', '2025-12', freq='M').astype(str)
QUARTERS = pd.period_range('2023Q1', '2025Q4', freq='Q').astype(str)
ORGS = {'Escherichia coli': 'urine_ecoli', 'Klebsiella pneumoniae': 'urine_kp',
        'Pseudomonas aeruginosa': 'urine_pa'}
UNKNOWN = {'', 'nan', 'none', 'null', 'unknown', 'not available', 'others', 'other'}

def state_name(x):
    s = str(x).strip().lower()
    if s in UNKNOWN:
        return 'Unknown state'
    return {'pondicherry': 'Puducherry', 'daman & diu': 'Dadra and Nagar Haveli and Daman and Diu',
            'dadra & nagar haveli': 'Dadra and Nagar Haveli and Daman and Diu',
            'jammu & kashmir': 'Jammu and Kashmir',
            'andaman & nicobar islands': 'Andaman and Nicobar Islands'}.get(s, s.title())

def district_name(x):
    s = ' '.join(str(x).strip().upper().split())
    return pd.NA if s.lower() in UNKNOWN else s

def pincode(x):
    s = str(x).strip().removesuffix('.0')
    return s if len(s) == 6 and s.isdigit() and s[0] != '0' else pd.NA

def save(df, name):
    df.to_csv(OUT / 'results' / (name + '.csv'), index=False)

def add_time(d):
    d['quarter'] = pd.PeriodIndex(d['month_year'], freq='M').asfreq('Q').astype(str)
    return d[d.month_year.between('2023-01', '2025-12')].copy()

rx = pd.read_csv(INPUTS['prescriptions'], low_memory=False)
rx['state'] = rx.patient_rx_state.map(state_name)
rx['district_clean'] = rx.district.map(district_name)
rx['pin'] = rx.patient_rx_pincode.map(pincode)
rx = add_time(rx)
dx = pd.read_csv(INPUTS['diagnostic_samples'], encoding='latin1', low_memory=False)
dx['month_year'] = pd.to_datetime(dx.end_state_date).dt.strftime('%Y-%m')
dx['state'] = dx.order_state.map(state_name)
dx['district_clean'] = dx.district.map(district_name)
dx['pin'] = dx.Order_Pincode.map(pincode)
dx['organism'] = dx.organism_isolated.map({x: cleaning.clean_organism(x) for x in dx.organism_isolated.unique()})
dx['urine'] = dx.test_name.eq('URINE CULTURE AND SENSITIVITY')
dx['no_growth'] = dx.organism.eq('No growth')

dx['culture_positive_convention'] = ~dx.no_growth
dx = add_time(dx)
ast = pd.read_csv(INPUTS['ast'], low_memory=False)
ast['state'] = ast.order_state.map(state_name)
ast['pin'] = ast.Order_Pincode.map(pincode)
ast['tested'] = ast[['RESISTANT', 'INTERMEDIATE', 'SENSITIVE']].sum(axis=1)
ast = add_time(ast)
ast = ast[(ast.test_name == 'URINE') & (ast.tested > 0)].copy()

assert rx.Rx_Count.sum() == 2784503
assert dx.Total_Test.sum() == 155764
assert dx.loc[dx.urine, 'Total_Test'].sum() == 146273
assert dx.loc[dx.culture_positive_convention, 'Total_Test'].sum() == 35657
assert (ast[['RESISTANT', 'INTERMEDIATE', 'SENSITIVE']] >= 0).all().all()


def coverage(keys):
    r = rx.groupby(keys, dropna=False).agg(antibiotic_records=('Rx_Count','sum'),
        rx_months=('month_year','nunique'), rx_pincodes=('pin','nunique'))
    d = dx.groupby(keys, dropna=False).agg(diagnostic_samples=('Total_Test','sum'),
        diagnostic_months=('month_year','nunique'), diagnostic_pincodes=('pin','nunique'))
    result = r.join(d, how='outer', sort=False)
    filters = {'culture_positive': dx.culture_positive_convention, 'urine_samples': dx.urine,
        'urine_positive': dx.urine & dx.culture_positive_convention}
    filters.update({short: dx.urine & dx.organism.eq(org) for org, short in ORGS.items()})
    for label, mask in filters.items():
        result[label] = dx.loc[mask].groupby(keys, dropna=False).Total_Test.sum()
    rm = rx.groupby(keys + ['month_year'], dropna=False).Rx_Count.sum().rename('rx')
    dm = dx.groupby(keys + ['month_year'], dropna=False).Total_Test.sum().rename('dx')
    both = pd.concat([rm, dm], axis=1).dropna()
    result['both_stream_months'] = both.reset_index().groupby(keys, dropna=False).size()
    eq = dx[dx.urine & dx.organism.eq('Escherichia coli')].groupby(keys + ['quarter'], dropna=False).Total_Test.sum()
    for threshold in [10,30,50]:
        result[f'ecoli_quarters_ge{threshold}'] = eq[eq >= threshold].reset_index().groupby(keys, dropna=False).size()
    return result.fillna(0).astype(int).reset_index()

state = coverage(['state']).sort_values('diagnostic_samples', ascending=False)
district = coverage(['state', 'district_clean'])
pin = coverage(['state', 'pin'])
save(district, 'district_counts_raw_labels')
save(pin, 'pincode_counts')


shared = pin[(pin.antibiotic_records>0) & (pin.diagnostic_samples>0) & pin.pin.notna()]
for col, value in [('shared_pincodes', shared.groupby('state').size()),
                   ('rx_in_shared_pincodes', shared.groupby('state').antibiotic_records.sum()),
                   ('dx_in_shared_pincodes', shared.groupby('state').diagnostic_samples.sum())]:
    state[col] = state.state.map(value).fillna(0).astype(int)
for stream, num, den in [('rx','rx_in_shared_pincodes','antibiotic_records'),('dx','dx_in_shared_pincodes','diagnostic_samples')]:
    state[f'{stream}_percent_in_shared_pincodes'] = 100*state[num]/state[den].replace(0,np.nan)
state['rx_percent_national'] = 100*state.antibiotic_records/rx.Rx_Count.sum()
state['diagnostic_percent_national'] = 100*state.diagnostic_samples/dx.Total_Test.sum()


within = []
for s, ds in district[district.diagnostic_samples>0].groupby('state'):
    ps = pin[(pin.state==s) & pin.pin.notna()]
    good = ds[(ds.district_clean.notna()) & (ds.urine_ecoli>=100) &
              (ds.ecoli_quarters_ge10>=8) & (ds.antibiotic_records>=1000)]
    goodpin = ps[(ps.urine_ecoli>=100) & (ps.ecoli_quarters_ge10>=8) & (ps.antibiotic_records>=1000)]
    top = ds.sort_values('urine_ecoli',ascending=False).iloc[0]
    within.append({'state':s, 'diagnostic_district_labels':int(ds.district_clean.nunique()),
        'top_ecoli_district':top.district_clean, 'top_district_ecoli':int(top.urine_ecoli),
        'top_district_ecoli_percent':100*top.urine_ecoli/max(ds.urine_ecoli.sum(),1),
        'districts_passing_exploratory_screen':len(good),
        'pincodes_passing_exploratory_screen':len(goodpin),
        'screened_district_names':'; '.join(good.sort_values('urine_ecoli',ascending=False).district_clean)})
within = pd.DataFrame(within)
save(within, 'within_state_feasibility')


maps=[]
for d, label in [(rx,'prescription'),(dx,'diagnostic')]:
    m=d[['state','pin','district_clean']].dropna().drop_duplicates()
    m['source']=label
    maps.append(m)
maps=pd.concat(maps)
conflict=maps.groupby(['state','pin']).district_clean.nunique()
bad=conflict[conflict>1].reset_index()[['state','pin']]
save(maps.merge(bad,on=['state','pin']), 'pincode_district_label_conflicts')
multi_state=maps.groupby('pin').state.nunique()
save(maps[maps.pin.isin(multi_state[multi_state>1].index)].sort_values(['pin','state']),
     'pincodes_recorded_under_multiple_states')

suspect=dx[(dx.state=='Delhi') & dx.district_clean.isin(['BUDAUN','GAUTAM BUDDHA NAGAR','HISAR','SONIPAT'])]
save(suspect.groupby(['state','district_clean','pin'],dropna=False).Total_Test.sum().rename('diagnostic_samples').reset_index(),
     'delhi_apparent_cross_state_district_labels')
save(dx.groupby(['state','organism'],dropna=False).Total_Test.sum().rename('diagnostic_records').reset_index(),
     'diagnostic_organism_counts_all_specimens')


agents=['Ciprofloxacin','Nitrofurantoin','Amikacin','Ceftriaxone','Ertapenem',
        'Amoxicillin/Clavulanic acid','Trimethoprim/Sulfamethoxazole','Meropenem','Fosfomycin']
aa=ast[ast.organism_clean.isin(ORGS) & ast.antibiotic_clean.isin(agents)]
keys=['state','organism_clean','antibiotic_clean']
ag=aa.groupby(keys).agg(tested=('tested','sum'),resistant=('RESISTANT','sum'),
    intermediate=('INTERMEDIATE','sum'),months_tested=('month_year','nunique'),pincodes=('pin','nunique'))
q=aa.groupby(keys+['quarter']).tested.sum()
for threshold in [10,30,50]:
    ag[f'quarters_ge{threshold}']=q[q>=threshold].reset_index().groupby(keys).size()
ag=ag.fillna(0).reset_index()
ag['resistance_pct']=100*ag.resistant/ag.tested
save(ag, 'urinary_ast_coverage_by_state_organism_agent')
save(q.rename('tested').reset_index(), 'urinary_ast_quarterly_tested')
state['recommended_scope']='Prescription description only'
state.loc[state.diagnostic_samples>0,'recommended_scope']='Limited diagnostic description'
state.loc[state.urine_ecoli>=100,'recommended_scope']='Pooled E. coli profile; agent coverage review'
state.loc[(state.urine_ecoli>=400)&(state.ecoli_quarters_ge30>=10)&(state.both_stream_months>=30),'recommended_scope']='Secondary longitudinal candidate'
state.loc[(state.urine_ecoli>=1000)&(state.ecoli_quarters_ge30==12)&(state.both_stream_months==36),'recommended_scope']='Primary longitudinal candidate'
save(state, 'state_counts_and_coverage')

smrx=rx.groupby(['state','month_year']).Rx_Count.sum().rename('antibiotic_records')
smdx=dx.groupby(['state','month_year']).Total_Test.sum().rename('diagnostic_samples')
smec=dx[dx.urine & dx.organism.eq('Escherichia coli')].groupby(['state','month_year']).Total_Test.sum().rename('urine_ecoli')
grid=pd.MultiIndex.from_product([state.state, MONTHS],names=['state','month_year'])
sm=pd.concat([smrx,smdx,smec],axis=1).reindex(grid).fillna(0).astype(int).reset_index()
save(sm,'state_monthly_counts')
aw=rx.groupby(['state','aware_2025']).Rx_Count.sum().rename('antibiotic_records').reset_index()
aw['percent_of_state_antibiotic_records']=100*aw.antibiotic_records/aw.groupby('state').antibiotic_records.transform('sum')
save(aw,'state_aware_composition')

assert state.antibiotic_records.sum()==2784503
assert state.diagnostic_samples.sum()==155764
assert state.urine_ecoli.sum()==21549
assert state.urine_kp.sum()==5759
assert state.urine_pa.sum()==1629
assert district.diagnostic_samples.sum()==155764
assert pin.diagnostic_samples.sum()==155764
assert (state.shared_pincodes<=state.diagnostic_pincodes).all()


plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
show=state[state.diagnostic_samples>0].sort_values('diagnostic_samples')
fig,axes=plt.subplots(1,3,figsize=(14,10),sharey=True,layout='constrained')
for ax,col,title,color in zip(axes,['antibiotic_records','diagnostic_samples','urine_ecoli'],
        ['Outpatient antibiotic records','Home-collected diagnostic samples','Urinary E. coli records'],['#2166ac','#b35806','#237a57']):
    ax.barh(show.state,show[col],color=color)
    ax.set_xscale('log');ax.set_title(title,fontsize=12);ax.set_xlabel('Count, logarithmic scale')
    ax.set_xlim(1,show[col].max()*5)
    for i,v in enumerate(show[col]):
        ax.text(max(v,1)*1.08,i,f'{v:,}',va='center',fontsize=8)
fig.suptitle('Tata 1mg state-level observation volume, 2023–2025',fontsize=16)
for ext in ['png','pdf']: fig.savefig(OUT/'figures'/f'state_counts.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)

primary=state[state.recommended_scope.isin(['Primary longitudinal candidate','Secondary longitudinal candidate'])].sort_values('urine_ecoli',ascending=False).state.tolist()
fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained',gridspec_kw={'width_ratios':[1.3,1]})
mat=ag[(ag.organism_clean=='Escherichia coli') & ag.state.isin(primary)].pivot(index='state',columns='antibiotic_clean',values='quarters_ge30').reindex(primary).reindex(columns=agents).fillna(0)
im=axes[0].imshow(mat,vmin=0,vmax=12,cmap='Blues',aspect='auto')
axes[0].set_xticks(range(len(agents)),agents,rotation=55,ha='right',fontsize=8)
axes[0].set_yticks(range(len(primary)),primary)
for i in range(len(mat)):
    for j in range(len(agents)):
        v=mat.iloc[i,j];axes[0].text(j,i,str(int(v)),ha='center',va='center',color='white' if v>7 else 'black')
axes[0].set_title('E. coli: quarters with at least 30 tests (of 12)',fontsize=11)
w=within.set_index('state').reindex(primary)
axes[1].barh(primary,w.top_district_ecoli_percent,color='#798c9c')
axes[1].invert_yaxis();axes[1].set_xlim(0,105);axes[1].set_xlabel('Urinary E. coli records in the largest district (%)')
axes[1].set_title('Within-state concentration',fontsize=11)
for i,v in enumerate(w.top_district_ecoli_percent):axes[1].text(v+1,i,f'{v:.1f}%',va='center',fontsize=9)
fig.suptitle('Feasibility of state and within-state analyses',fontsize=15)
for ext in ['png','pdf']:fig.savefig(OUT/'figures'/f'analysis_feasibility.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)

manifest={'study_period':'2023-01 to 2025-12','python':platform.python_version(),
    'pandas':pd.__version__,'numpy':np.__version__,'matplotlib':matplotlib.__version__,
    'inputs':{k:{'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in INPUTS.items()},
    'checks':{'antibiotic_records':int(rx.Rx_Count.sum()),'diagnostic_samples':int(dx.Total_Test.sum()),
    'urine_samples':int(state.urine_samples.sum()),'urine_ecoli':int(state.urine_ecoli.sum()),
    'urine_kp':int(state.urine_kp.sum()),'urine_pa':int(state.urine_pa.sum()),
    'invalid_or_missing_rx_pincode_records':int(rx.loc[rx.pin.isna(),'Rx_Count'].sum()),
    'invalid_or_missing_dx_pincode_samples':int(dx.loc[dx.pin.isna(),'Total_Test'].sum()),
    'state_pincodes_with_conflicting_district_labels':len(bad),
    'pincodes_recorded_under_multiple_states':int((multi_state>1).sum()),
    'delhi_samples_with_apparent_cross_state_district_labels':int(suspect.Total_Test.sum())},
    'scope':'Descriptive counts and feasibility only. No trends, causal estimates or patient linkage.'}
(OUT/'results/validation_and_sources.json').write_text(json.dumps(manifest,indent=2))
cols=['state','antibiotic_records','diagnostic_samples','urine_ecoli','urine_kp','urine_pa','both_stream_months','recommended_scope']
print(state[cols].to_string(index=False))
print('\nWITHIN STATE\n'+within[within.state.isin(primary)].to_string(index=False))
print('\nVALIDATION\n'+json.dumps(manifest['checks'],indent=2))
