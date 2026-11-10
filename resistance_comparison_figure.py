"""NCDC outpatient comparison and monthly E. coli and K. pneumoniae resistance."""


from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/figures/resistance_comparison'
OUT.mkdir(parents=True, exist_ok=True)
comparison = pd.read_csv(ROOT / 'results/surveillance/source_data/resistance_comparison_2024.csv')
data = comparison[comparison.organism.eq('Escherichia coli')].dropna(subset=['tata_R_percent', 'ncdc_OPD_R_percent'])
data = data.sort_values(['ncdc_OPD_R_percent', 'antibiotic'], ascending=[False, True])
data = data.rename(columns={'tata_R': 'tata_resistant', 'tata_N': 'tata_tested',
    'tata_R_percent': 'tata_resistance_percent', 'ncdc_OPD_R_percent': 'ncdc_outpatient_resistance_percent'})
data = data[['antibiotic', 'tata_resistant', 'tata_tested',
    'tata_resistance_percent', 'ncdc_outpatient_resistance_percent']]
assert len(data) == 16
assert np.allclose(data.tata_resistance_percent, 100 * data.tata_resistant / data.tata_tested)
data.to_csv(OUT / 'panel_a_source_data.csv', index=False)
monthly = pd.read_csv(ROOT / 'results/figures/supp_fig5-7_resistance_monthly_data.csv')
monthly = monthly.rename(columns={'antibiotic': 'antibiotic_lower'})
ast = pd.read_csv(ROOT / 'data/ast_clean.csv', usecols=['test_name', 'organism_clean',
    'antibiotic_clean', 'RESISTANT', 'INTERMEDIATE', 'SENSITIVE'])
ast = ast[ast.test_name.eq('URINE')].copy()
ast['tested'] = ast[['RESISTANT', 'INTERMEDIATE', 'SENSITIVE']].sum(axis=1)
ast['antibiotic'] = ast.antibiotic_clean.str.lower()
orders = {}
for organism in ['escherichia coli', 'klebsiella pneumoniae']:
    counts = ast[ast.organism_clean.str.lower().eq(organism)].groupby('antibiotic').tested.sum()
    orders[organism] = counts.sort_values(ascending=False).head(10).index.tolist()
colours = dict(zip(orders['escherichia coli'], matplotlib.colors.TABLEAU_COLORS.values()))
for panel, organism in [('b', 'escherichia coli'), ('c', 'klebsiella pneumoniae')]:
    drugs = [d for d in orders[organism] if not (panel == 'c' and d == 'ampicillin')]
    frame = monthly[monthly.organism.eq(organism) & monthly.antibiotic_lower.isin(drugs)]
    assert frame.antibiotic_lower.nunique() == len(drugs)
    frame.to_csv(OUT / f'panel_{panel}_source_data.csv', index=False)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
    'axes.labelsize': 8, 'axes.titlesize': 9, 'axes.linewidth': .55,
    'xtick.labelsize': 6.8, 'ytick.labelsize': 7,
    'xtick.major.width': .5, 'ytick.major.width': .5,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
fig = plt.figure(figsize=(7.5, 10.5), facecolor='white')


data = pd.read_csv(OUT / 'panel_a_source_data.csv')
ax = fig.add_axes([.285, .720, .690, .250])
y = np.arange(len(data))
for i in range(0, len(data), 2):
    ax.axhspan(i-.48, i+.48, color='#F5F6F7', zorder=0)
for offset, col, colour, label, fmt in [
    (-.18, 'tata_resistance_percent', '#167D9A', 'Tata 1mg', '.1f'),
    (.18, 'ncdc_outpatient_resistance_percent', '#D48739', 'NCDC outpatients', '.0f')]:
    bars = ax.barh(y+offset, data[col], height=.32, color=colour, label=label, zorder=3)
    ax.bar_label(bars, labels=[format(v, fmt) for v in data[col]], padding=2, fontsize=5.8)
ax.set_yticks(y, [s[0].upper()+s[1:].lower() for s in data.antibiotic])
ax.set_ylim(len(data)-.5, -.5)
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.set_xlabel('Resistance (%)', labelpad=4)
ax.set_ylabel('Antibiotic', labelpad=6)
ax.set_title(r'$\it{Escherichia\ coli}$ — 2024', loc='left', pad=8)
ax.tick_params(axis='y', length=0, pad=5, labelsize=6.8)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.grid(axis='x', color='#E4E7E9', linewidth=.4)
ax.set_axisbelow(True)
ax.legend(loc='lower right', bbox_to_anchor=(1.01, 1.015), ncol=2,
          frameon=False, fontsize=6.8, handlelength=1.3,
          columnspacing=1.2, handletextpad=.5, borderaxespad=0)
fig.text(.023, .985, 'a', fontsize=12, weight='bold', va='top')

legends = []
for panel, bottom, legend_y, upper, species in [
    ('b', .440, .386, 100, r'$\it{Escherichia\ coli}$'),
    ('c', .100, .058, 90, r'$\it{Klebsiella\ pneumoniae}$')]:
    df = pd.read_csv(OUT / f'panel_{panel}_source_data.csv')
    months = sorted(df.month_year.unique())
    x = np.arange(len(months))
    ax = fig.add_axes([.095, bottom, .880, .200])
    for drug in colours:
        if drug not in set(df.antibiotic_lower):
            continue
        d = df[df.antibiotic_lower == drug].set_index('month_year').loc[months]
        label = f'{drug[0].upper()+drug[1:]} (n={int(d.N.sum()):,})'
        ax.fill_between(x, d.ci_low, d.ci_high, color=colours[drug], alpha=.15, linewidth=0)
        ax.plot(x, d.R_percent, color=colours[drug], linewidth=.9,
                marker='o', markersize=2.4, markeredgewidth=0, label=label)
    ax.set_xlim(-.6, 35.6)
    ax.set_ylim(0, upper)
    ax.set_yticks(np.arange(0, upper+1, 10))
    ticks = list(range(0, 36, 3)) + [35]
    tick_labels = [pd.Timestamp(months[i]).strftime('%b\n%Y') for i in ticks]
    ax.set_xticks(ticks, tick_labels)
    ax.set_xticks(x, minor=True)
    ax.tick_params(axis='x', which='minor', length=1.5, width=.4)
    ax.tick_params(axis='x', which='major', pad=3)
    ax.set_xlabel('Month–year', labelpad=3)
    ax.set_ylabel('Resistance (%)', labelpad=5)
    ax.set_title(species + ' — 2023–2025', loc='left', pad=8)
    ax.spines[['top', 'right']].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i].casefold())
    legend = fig.legend([handles[i] for i in order], [labels[i] for i in order],
        loc='upper center', bbox_to_anchor=(.53, legend_y), ncol=3,
        frameon=False, fontsize=6.3, handlelength=1.6, handletextpad=.45,
        columnspacing=1.25, labelspacing=.55, borderaxespad=0)
    legends.append(legend)
    fig.text(.023, bottom+.218, panel, fontsize=12, weight='bold', va='top')
    assert df.ci_high.max() < upper
    assert len(ax.lines) == (10 if panel == 'b' else 9)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for legend in [fig.axes[0].get_legend()] + legends:
    bbox = legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    assert 0 < bbox.x0 < bbox.x1 < 1 and 0 < bbox.y0 < bbox.y1 < 1, bbox
    print('legend bounds', tuple(round(v, 4) for v in bbox.bounds))
    assert all('\n' not in t.get_text() for t in legend.get_texts())
    boxes = [t.get_window_extent(renderer) for t in legend.get_texts()]
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i+1:])
stem = OUT / 'resistance_comparison_and_trends'
fig.savefig(stem.with_suffix('.png'), dpi=600)
fig.savefig(stem.with_suffix('.svg'))
fig.savefig(stem.with_suffix('.pdf'))
fig.savefig(OUT / 'preview.png', dpi=160)
plt.close(fig)
metadata = {
    'study_period': ['2023-01', '2025-12'], 'comparison_year': 2024,
    'drug_colours': colours,
    'selection': 'Ten most-tested E. coli drugs; original ten K. pneumoniae drugs with ampicillin removed, leaving nine.',
    'counts': 'Antibiotic-specific R+I+S results, not unique patients.',
    'comparison': 'NCDC outpatient percentages at published precision; drug-specific N unavailable.',
}
(OUT / 'figure_details.json').write_text(json.dumps(metadata, indent=2) + '\n')
print(stem)
