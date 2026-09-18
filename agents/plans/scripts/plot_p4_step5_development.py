#!/usr/bin/env python3
"""Plot response support, full-image source wings, and frozen held-out null scores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main() -> None:
    """Write a standalone scientific figure from archived, already-calculated results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--development', type=Path, required=True)
    parser.add_argument('--nulls', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    development = json.loads(args.development.read_text())
    nulls = json.loads(args.nulls.read_text())
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), layout='constrained')
    support = development['response_support']['10-20']['fractions']
    radii = [f['radius_pixels'] for f in support]
    axes[0, 0].plot(radii, [f['energy_fraction_p05_median_p95'][1] for f in support], 'o-', label='all signed response')
    axes[0, 0].plot(radii, [f['negative_energy_fraction_p05_median_p95'][1] for f in support], 's--', label='negative lobes only')
    axes[0, 0].axvline(3.6, color='0.6', lw=1, label='one λ/D')
    axes[0, 0].set(title='Stored response support\nSeparation 10–20 pixels',
                   xlabel='Circle radius [pixels]', ylabel='Median fraction of 11×11 stamp energy', ylim=(0, 1.05))
    axes[0, 0].legend(fontsize=8)
    for record in development['full_image_leakage']:
        radius = ((record['trial']['row']-127.5)**2 + (record['trial']['column']-127.5)**2)**0.5
        axes[0, 1].plot([m['radius_pixels'] for m in record['masks']],
                        [m['difference_energy_fraction_inside_circle'] for m in record['masks']],
                        'o-', label=f'separation {radius:.1f} px')
    axes[0, 1].axvline(7.3, color='0.6', lw=1, ls='--', label='trial mask: 7.3 px')
    axes[0, 1].set(title='Full-image finite-source support', xlabel='Circle radius [pixels]',
                   ylabel='Fraction of full difference-image energy', ylim=(0, 1.05))
    axes[0, 1].legend(fontsize=8)
    titles = {'identity': 'Identity', 'diagonal': 'Diagonal', 'pca0': 'Local scale + mean (zero modes)', 'pca': 'PCA: three modes'}
    for axis, model in zip((axes[0, 2], *axes[1]), titles):
        threshold = nulls['models'][model]['threshold']
        for role, color, marker in (('calibration', '0.55', 'o'), ('evaluation', '#5b46a8', '^')):
            rows = [r for r in nulls['trials'] if r['role'] == role and r['common_eligible']]
            axis.scatter([r['nominal_radius'] for r in rows],
                         [r['models'][model]['search_score']/threshold for r in rows],
                         color=color, marker=marker, s=24, alpha=0.8, label=role)
        axis.axhline(1, color='#b33a30', ls='--', lw=1)
        count = nulls['models'][model]['evaluation_exceedances']
        axis.set(title=f'{titles[model]}\n{count}/28 evaluation exceedances', xlabel='Separation [pixels]',
                 ylabel='Search score / frozen threshold')
        axis.legend(fontsize=8)
    for axis in axes.ravel():
        axis.grid(alpha=0.2)
    fig.suptitle('Step 5 development and null calibration\nOverlapping spatial trials; no detection-gain claim', fontsize=12)
    fig.savefig(args.output, dpi=180)


if __name__ == '__main__':
    main()
