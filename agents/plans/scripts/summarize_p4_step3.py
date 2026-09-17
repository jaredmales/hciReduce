#!/usr/bin/env python3
"""Export compact tables and scientific figures from completed P4 Step-3 trials.

Keep all position/brightness trials in the exported CSV, including failed fits.
The plots describe conditional response calibration, not detection completeness.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from run_p4_response_convergence import fingerprint


def distribution(rows: list[dict], key: str) -> dict:
    """Summarize every finite observation while reporting missing results explicitly."""
    values = np.asarray([row[key] for row in rows if row[key] is not None], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f'nonfinite values in {key}')
    return {'count': len(values), 'missing': len(rows)-len(values),
            'minimum': float(values.min()) if len(values) else None,
            'median': float(np.median(values)) if len(values) else None,
            'maximum': float(values.max()) if len(values) else None}


def injection_figure(rows: list[dict], destination: Path) -> None:
    """Show all six analytic-template recoveries as brightness changes."""
    combinations = [name for name in ('mean', 'sigmaMean') if any(row['combination'] == name for row in rows)]
    figure, axes = plt.subplots(3, len(combinations), figsize=(5*len(combinations), 9),
                               sharex=True, sharey='row', squeeze=False, layout='constrained')
    metrics = [('fixed_position_bias_fraction', 'Fixed-position contrast bias (%)', 100),
               ('position_error_pixels', 'Fitted position error (pixels)', 1),
               ('cosine', 'Template / injected-response cosine', 1)]
    colors = plt.get_cmap('tab10').colors
    reference = sorted({row['amplitude'] for row in rows})[1]
    for column, combination in enumerate(combinations):
        for position in sorted({row['position'] for row in rows}):
            trials = sorted([row for row in rows if row['combination'] == combination and
                             row['position'] == position and row['template'] == 'analytic'],
                            key=lambda row: row['amplitude'])
            amplitudes = [row['amplitude']/reference for row in trials]
            label = f"P{position}: ({trials[0]['row']}, {trials[0]['column']})"
            for row_index, (key, title, scale) in enumerate(metrics):
                values = [float('nan') if row[key] is None else scale*row[key] for row in trials]
                axes[row_index, column].plot(amplitudes, values, 'o-', color=colors[position],
                                             markersize=4, linewidth=1.2, label=label)
                axes[row_index, column].set_ylabel(title)
                axes[row_index, column].grid(alpha=.2)
        axes[0, column].set_title('Arithmetic mean' if combination == 'mean' else 'Sigma-clipped mean')
        axes[0, column].axhline(0, color='black', linewidth=.6)
        axes[2, column].axhline(1, color='black', linewidth=.6)
        axes[2, column].set_xscale('log', base=4)
        axes[2, column].set_xticks([.25, 1, 4], ['0.25', '1', '4'])
        axes[2, column].set_xlabel('Injected contrast / AF Lep reference contrast')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=3 if len(combinations) == 2 else 2, frameon=False)
    title = ('P4 analytic response: conditional injection calibration\n'
             '621 frames; six positions; baseline-subtracted residuals') if len(combinations) == 2 else (
             'P4 analytic response: injection calibration\n621 frames; six positions\n'
             'Baseline-subtracted residuals')
    figure.suptitle(title, fontsize=13)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def timing_figure(rows: list[dict], controls: list[dict], destination: Path) -> None:
    """Plot median runtime and memory for controlled two-worker subset trials."""
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    variants = ('uncached_dense', 'cached_dense', 'cached_sparse')
    labels = ('No source cache', 'Source cache', 'Cache + sparse product')
    for column, frames in enumerate((24, 96)):
        for batch, offset, color in ((1, -.18, '#527eaa'), (32, .18, '#db9654')):
            trials = [next(row for row in rows if row['frames'] == frames and
                           row['variant'] == variant and row['requested_batch'] == batch)
                      for variant in variants]
            label = 'Batch 1' if batch == 1 else f"Batch {trials[0]['realized_batch']}"
            for row_index, (key, scale, title) in enumerate((
                    ('median_wall_seconds', 1, 'Wall time (seconds)'),
                    ('median_peak_rss_kib', 1/1024, 'Peak resident memory (MiB)'))):
                axes[row_index, column].bar(np.arange(3)+offset,
                    [trial[key]*scale for trial in trials], width=.36, label=label, color=color)
                axes[row_index, column].set_xticks(np.arange(3), labels, rotation=15, ha='right')
                axes[row_index, column].set_ylabel(title)
                axes[row_index, column].grid(axis='y', alpha=.2)
        axes[0, column].set_title(f'{frames} frames')
        for precision, style in (('P4-M32D64', '--'), ('P4-D64', ':')):
            control = next(row for row in controls if (row['frames'], row['precision']) == (frames, precision))
            label = 'Paired refit: ' + precision.removeprefix('P4-')
            axes[0, column].axhline(control['median_wall_seconds'], color='black', linestyle=style,
                                   linewidth=1, label=label)
            axes[1, column].axhline(control['median_peak_rss_kib']/1024, color='black', linestyle=style, linewidth=1)
        axes[0, column].legend(frameon=False)
    figure.suptitle('P4 response products: controlled serial trials\n'
                    'Median of three trials; 2 OpenMP workers; 1 BLAS thread', fontsize=13)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    """Export reproducible summaries after all full-data and controlled trials finish."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.experiment.resolve(), args.output.resolve()
    science_path = root / 'scientific_comparison/summary.json'
    injections_path = root / 'injections/summary.json'
    timing_path = root / 'reuse_benchmark/summary.json'
    control_path = root / 'reuse_benchmark/refit_summary.json'
    science = json.loads(science_path.read_text())
    injections = json.loads(injections_path.read_text())
    timings = json.loads(timing_path.read_text())
    controls = json.loads(control_path.read_text())
    if len(injections) != 72 or len(timings) != 12 or len(controls) != 4:
        raise ValueError('expected 72 injection comparisons, 12 analytic timing groups and 4 refit controls')
    metrics = ('fixed_position_bias_fraction', 'contrast_bias_fraction', 'position_error_pixels',
               'cosine', 'relative_template_error', 'best_scaled_template_error',
               'paired_relative_template_error', 'one_sided_vs_paired_relative_error',
               'raw_contrast_bias_fraction', 'raw_position_error_pixels')
    groups = []
    for combination in ('mean', 'sigmaMean'):
        for template in ('analytic', 'refit'):
            for amplitude in sorted({row['amplitude'] for row in injections}):
                rows = [row for row in injections if (row['combination'], row['template'], row['amplitude']) ==
                        (combination, template, amplitude)]
                if len(rows) != 6:
                    raise ValueError('missing injection positions')
                groups.append({'combination': combination, 'template': template, 'amplitude': amplitude,
                    'fit_status_counts': {status: sum(row['fit_status'] == status for row in rows)
                                         for status in sorted({row['fit_status'] for row in rows})},
                    'raw_fit_status_counts': {status: sum(row['raw_fit_status'] == status for row in rows)
                                             for status in sorted({row['raw_fit_status'] for row in rows})},
                    **{key: distribution(rows, key) for key in metrics}})
    output.mkdir(parents=True, exist_ok=False)
    for source, name in ((science_path, 'science.json'), (injections_path, 'injections.json'),
                         (root / 'injections/summary.csv', 'injections.csv'),
                         (timing_path, 'timing.json'), (control_path, 'refit_timing.json')):
        shutil.copy2(source, output / name)
    summary = {'schema': 1, 'injection_groups': groups,
               'provenance': [fingerprint(path) for path in
                              (science_path, injections_path, timing_path, control_path, Path(__file__).resolve())]}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    injection_figure(injections, output / 'injections.png')
    timing_figure(timings, controls, output / 'timing.png')
    print(output)


if __name__ == '__main__':
    main()
