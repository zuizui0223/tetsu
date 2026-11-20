"""Render only aggregate predictive gains; no source site/taxon data."""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT/'data/derived/v304_incremental_foreign_rail.json').read_text(encoding='utf-8'))
    dest = ROOT/'docs/manuscript_assets'
    dest.mkdir(exist_ok=True)
    national = [r for r in data['national'] if r['penalty'] == .01 and r['cz_reference'] == 'urban'
                and not r['small_sample_external_diagnostic']]
    local = [r for r in data['local_kashmir']['panels'] if r['penalty'] == .01
             and r['panel'] == 'foreign_two_study_support']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw={'width_ratios': [1, 1.35]})
    colors = ['#267f83' if r['foreign_increment']['log_loss'] > 0 else '#b7574f' for r in national]
    values = [r['foreign_increment']['log_loss'] for r in national]
    axes[0].barh(np.arange(len(national)), values, color=colors)
    axes[0].set_yticks(np.arange(len(national)),
                      [f"{r['held_rail_study']}  ({r['target_taxa']} taxa)" for r in national])
    axes[0].invert_yaxis()
    axes[0].set_title('A. Entire railway study withheld\nDomestic reference + foreign rail vs reference', fontsize=11)
    axes[0].set_xlabel('Log-loss reduction (nat per species)')
    for j, v in enumerate(values):
        axes[0].annotate(f'{v:+.4f}', (v, j), xytext=(4 if v >= 0 else -4, 0),
                         textcoords='offset points', ha='left' if v >= 0 else 'right', va='center', fontsize=9)
    axes[0].set_xlim(-.035, .052)
    contrasts = ['foreign_beyond_local', 'foreign_beyond_local_and_track', 'matched_foreign_rail_vs_road']
    labels = ['Add foreign rail to local vegetation', 'Add foreign rail after other domestic rail sites',
              'Foreign rail vs foreign road (same 2 countries)']
    rows = []
    for r, color, offset in zip(local, ['#267f83', '#b16a2b'], [-.13, .13]):
        for j, key in enumerate(contrasts):
            v = r['comparisons'][key]['log_loss']
            lo, hi = v['block_bootstrap_95_percentile']
            axes[1].plot([lo, hi], [j+offset]*2, color=color, lw=2)
            axes[1].plot(v['gain'], j+offset, 'o', color=color,
                         label=f"{r['year']}: {r['sites']} sites" if j == 0 else None)
            rows.append({'analysis': 'Kashmir_58_taxa', 'group': r['year'], 'contrast': key,
                         'gain': v['gain'], 'lower': lo, 'upper': hi})
    axes[1].set_yticks(range(3), labels, fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set_title('B. Kashmir: 5 spatial blocks withheld\nSame-site, same-year, equal-area local references', fontsize=11)
    axes[1].set_xlabel('Log-loss reduction (nat per site-species)')
    axes[1].xaxis.set_major_locator(MaxNLocator(5))
    axes[1].legend(loc='lower right', frameon=False)
    for ax in axes:
        ax.axvline(0, color='#555555', lw=.8)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', alpha=.15)
    fig.suptitle('Foreign railway frequency adds information, but specificity is conditional', fontsize=13)
    fig.text(.04, .025,
             'Positive = better prediction. A: national references are not local matched controls.\n'
             'B: 58 externally selected taxa; bars are 999 resamples of 5 spatial blocks, conditional on fitted predictions.\n'
             'Same Kashmir study in two years; no global or causal confidence interval.', fontsize=8, color='#444444')
    fig.tight_layout(rect=[0, .15, 1, .93], w_pad=3)
    for suffix in ('png', 'svg'):
        fig.savefig(dest/f'v304_foreign_rail_increment.{suffix}', dpi=180, bbox_inches='tight')
    svg = dest/'v304_foreign_rail_increment.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',
                   encoding='utf-8')
    for r in national:
        rows.append({'analysis': 'national_study_holdout', 'group': r['held_rail_study'],
                     'contrast': 'foreign_beyond_domestic', 'gain': r['foreign_increment']['log_loss'],
                     'lower': '', 'upper': ''})
    with (dest/'v304_foreign_rail_increment.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
