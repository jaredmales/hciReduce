#!/usr/bin/env python3
"""Audit completed Step-5 products and archive recovery/photometry review figures.

This reads the frozen experiment without rerunning reductions or changing its
thresholds. Paired amplitude increments are secondary diagnostics only; detection,
raw contrast error, and interval coverage retain the original positive images.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from astropy.io import fits

MODELS = ('identity', 'diagonal', 'pca0', 'pca')
LABELS = ('Identity', 'Diagonal', 'Zero-mode PCA', 'PCA: three modes')
COLORS = ('#606060', '#c47b16', '#22866f', '#7052a3')
LEVELS = (0.5, 1.0, 2.0)
OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))


def fingerprint(path: Path) -> dict:
    """Record exact file content independently of the experiment runner."""
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def read_json(path: Path):
    """Read one frozen record."""
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    """Stop the review on an inconsistent product, including under Python optimization."""
    if not condition:
        raise RuntimeError(message)


def audit(root: Path) -> tuple[dict, dict, dict]:
    """Recompute every recorded measurement from FITS and verify frozen file hashes."""
    checked = {}

    def verify(expected: dict) -> None:
        """Verify a file once while rejecting conflicting recorded fingerprints."""
        path = expected['path']
        if path not in checked:
            checked[path] = fingerprint(Path(path))
        require(checked[path] == expected, 'fingerprint mismatch: ' + path)

    manifest = read_json(root / 'manifest.json')
    reduction_manifest = read_json(root / 'reductions/manifest.json')
    complete = read_json(root / 'complete.json')
    for record in (manifest['frozen_records'] + reduction_manifest['frozen_records'] +
                   reduction_manifest['input_records'] + [manifest['holdout_null_results'], complete['results']]):
        verify(record)
    require(read_json(root / 'state.json')['status'] == 'complete', 'analysis unfinished')
    require(read_json(root / 'reductions/state.json')['status'] == 'complete', 'reductions unfinished')
    require(complete['frozen_inputs_unchanged'], 'runner reported changed inputs')
    jobs = read_json(root / 'jobs.json')
    design = read_json(root / 'injection_design.json')
    nulls = read_json(root / 'null_results.json')
    protocol = read_json(root / 'protocol.json')
    data = read_json(root / 'results.json')
    require(len(jobs) == complete['reductions'] == 18, 'unexpected trial count')
    require(len(data['measurements']) == complete['model_measurements'] == 72, 'unexpected measurement count')
    require(jobs == reduction_manifest['trials'], 'reduction jobs differ from evaluation jobs')
    require(data['null_models'] == nulls['models'], 'threshold records changed')
    require(len(reduction_manifest['input_records']) == manifest['source_frame_count'] == 621, 'frame count changed')
    require({(j['nominal_radius'], j['angular_block']) for j in jobs} ==
            {tuple(site) for site in design['sites']}, 'injection sites changed')
    times, science_hashes = [], {}
    for job in jobs:
        directory = root / 'reductions' / job['name']
        finished = read_json(directory / 'complete.json')
        require(finished['trial'] == job, 'completed reduction has different trial parameters')
        for record in finished['products']:
            verify(record)
        science = directory / 'finim.fits'
        science_hashes[job['name']] = checked[str(science)]['sha256']
        image, header = fits.getdata(science, header=True)
        require(image.squeeze().shape == (256, 256), 'not a full-image reduction')
        require(header['COMBINATION METHOD'].strip() == 'mean' and
                int(header['P4 LOCAL STAMP SIZE']) == 0, 'wrong reduction contract')
        times.append(finished['elapsed_seconds'])
        baseline = next(t for t in nulls['trials'] if t['name'] == job['null_trial'])
        require(baseline['common_eligible'] and baseline['role'] == 'evaluation', 'ineligible injection site')
        require((job['row'], job['column']) == (baseline['row'], baseline['column']), 'injection position changed')
        scale = nulls['models']['pca0']['threshold'] * baseline['models']['pca0']['center_conditional_sigma']
        require(job['reference_contrast_scale'] == scale and
                job['contrast'] == scale * job['brightness_multiplier'], 'brightness definition changed')
    require({(r['trial']['name'], r['model']) for r in data['measurements']} ==
            {(j['name'], m) for j in jobs for m in MODELS}, 'duplicate or missing measurements')
    require({j['brightness_multiplier'] for j in jobs} == set(design['brightness_multipliers']) == set(LEVELS),
            'brightness levels changed')
    for record in data['measurements']:
        job, model = record['trial'], record['model']
        require(job in jobs, 'measurement trial differs from frozen jobs')
        directory = root / 'analysis' / job['name'] / model
        require(read_json(directory / 'measurement.json') == record, 'individual measurement differs from summary')
        for product in record['products']:
            verify(product)
        require(checked[str(directory / 'science.fits')]['sha256'] == science_hashes[job['name']],
                'analysis science differs from full positive image')
        maps = {}
        for role in ('psf_amplitude', 'psf_sigma', 'psf_score', 'noise_status', 'noise_samples'):
            image, header = fits.getdata(directory / f'science_{role}.fits', header=True)
            maps[role] = image.squeeze()
            require(header['HCIA NOISE PRODUCT'] == role and header['HCIA NOISE CALIBRATED'] == 0,
                    'wrong conditional-map metadata')
            require(header['HCIA NOISE MODEL'] == ('pca' if model == 'pca0' else model) and
                    header['HCIA NOISE MAX MODES'] == (0 if model == 'pca0' else 3) and
                    header['HCIA NOISE MIN SAMPLES'] == 8 and header['HCIA NOISE FLOOR FRACTION'] == 0.1 and
                    header['HCIA NOISE ARC STEP'] == 5 and header['HCIA NOISE GUARD'] == 0 and
                    header['HCIA NOISE EXCLUSION'] == 'NONZERO_STENCIL_PIXEL_CENTERS', 'noise policy changed')
        command = read_json(directory / 'command.json')
        for index, key in enumerate(('Rows', 'Columns', 'Radii')):
            option = '--noise.exclude' + key + '=' + ','.join(str(c[index]) for c in protocol['extra_exclusions'])
            require(option in command, 'held-out exclusions changed')
        x, y = job['row'], job['column']
        statuses = [float(maps['noise_status'][y + dy, x + dx]) for dx, dy in OFFSETS]
        require(all(s == 0 for s in statuses), 'this completed pilot unexpectedly contains an invalid search')
        amplitude = float(maps['psf_amplitude'][y, x])
        sigma = float(maps['psf_sigma'][y, x])
        score = max(float(maps['psf_score'][y + dy, x + dx]) for dx, dy in OFFSETS)
        threshold = nulls['models'][model]['threshold']
        expected = {'statuses': statuses, 'search_valid': True, 'center_amplitude': amplitude,
            'center_conditional_sigma': sigma, 'search_score': score, 'threshold': threshold,
            'detected': score > threshold, 'center_training_samples': int(maps['noise_samples'][y, x]),
            'contrast_bias_fraction': amplitude / job['contrast'] - 1,
            'inside_one_conditional_sigma': abs(amplitude - job['contrast']) <= sigma}
        require(all(record[k] == v for k, v in expected.items()), 'FITS measurement mismatch')
    for group in data['groups']:
        rows = [r for r in data['measurements'] if r['model'] == group['model'] and
                r['trial']['brightness_multiplier'] == group['brightness_multiplier']]
        require(group['trials'] == len(rows) == 6 and group['invalid_searches'] == 0 and
                group['detections'] == sum(r['detected'] for r in rows) and
                group['inside_one_conditional_sigma'] == sum(r['inside_one_conditional_sigma'] for r in rows) and
                group['median_contrast_bias_fraction'] == np.median([r['contrast_bias_fraction'] for r in rows]),
                'aggregate result mismatch')
    verification = {'unique_files_hashed': len(checked), 'source_frames': 621,
        'full_image_reductions_checked': len(jobs), 'measurements_recomputed_from_fits': len(data['measurements']),
        'all_checks_passed': True,
        'provenance': [fingerprint(root / name) for name in ('manifest.json', 'reductions/manifest.json',
            'complete.json', 'jobs.json', 'protocol.json', 'injection_design.json', 'null_results.json', 'results.json')],
        'reduction_seconds': {'sum': sum(times), 'median': float(np.median(times)), 'minimum': min(times), 'maximum': max(times)},
        'timing_scope': 'Recorded per-reduction wall times; excludes supervisor and subsequent filter analysis.'}
    return data, nulls, verification


def summarize(data: dict, nulls: dict, verification: dict) -> dict:
    """Describe raw errors and secondary paired increments without changing primary outcomes."""
    rows = []
    by_name = {t['name']: t for t in nulls['trials']}
    for record in data['measurements']:
        job, model = record['trial'], record['model']
        baseline = by_name[job['null_trial']]['models'][model]
        rows.append({'trial': job['name'], 'model': model,
            'standardized_error': (record['center_amplitude'] - job['contrast']) / record['center_conditional_sigma'],
            'baseline_center_amplitude': baseline['center_amplitude'],
            'baseline_search_score': baseline['search_score'],
            'baseline_detected': baseline['search_score'] > record['threshold'],
            'paired_increment_error_fraction': (record['center_amplitude'] - baseline['center_amplitude']) / job['contrast'] - 1})
    groups = []
    for group in data['groups']:
        chosen = [r for r in data['measurements'] if r['model'] == group['model'] and
                  r['trial']['brightness_multiplier'] == group['brightness_multiplier']]
        paired = [d['paired_increment_error_fraction'] for d in rows if d['model'] == group['model'] and
                  d['trial'] in {r['trial']['name'] for r in chosen}]
        groups.append({**group,
            'raw_error_range': [min(r['contrast_bias_fraction'] for r in chosen), max(r['contrast_bias_fraction'] for r in chosen)],
            'paired_increment_error_min_median_max': [min(paired), float(np.median(paired)), max(paired)]})
    models = {}
    for model in MODELS:
        chosen = [r for r in rows if r['model'] == model]
        centers = [t['models'][model] for t in nulls['trials'] if t['role'] == 'evaluation' and t['common_eligible']]
        models[model] = {'positive_inside_one_conditional_sigma': sum(abs(r['standardized_error']) <= 1 for r in chosen),
            'positive_trials': len(chosen), 'median_absolute_standardized_error': float(np.median([abs(r['standardized_error']) for r in chosen])),
            'null_centers_inside_one_conditional_sigma': sum(abs(r['center_amplitude']) <= r['center_conditional_sigma'] for r in centers),
            'null_trials': len(centers),
            'preexisting_detected_sites': sorted({r['trial'].rsplit('_l', 1)[0] for r in chosen if r['baseline_detected']})}
    return {'verification': verification, 'groups': groups, 'models': models, 'diagnostics': rows,
        'diagnostic_caveat': 'Paired increments subtract separately measured baseline amplitudes. Learned weights adapt to each full positive image. This is not detection completeness or raw photometric interval coverage.'}


def plot(data: dict, review: dict, output: Path) -> None:
    """Show individual site recovery and photometric errors with their conditioning explicit."""
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10})
    rows = data['measurements']
    diagnostics = {(r['trial'], r['model']): r for r in review['diagnostics']}
    sites = sorted({(r['trial']['nominal_radius'], r['trial']['angular_block']) for r in rows})
    site_colors = plt.get_cmap('tab10').colors[:len(sites)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True, layout='constrained')
    for axis, model, label in zip(axes.flat, MODELS, LABELS):
        for site, color in zip(sites, site_colors):
            selected = sorted([r for r in rows if r['model'] == model and
                (r['trial']['nominal_radius'], r['trial']['angular_block']) == site],
                key=lambda r: r['trial']['brightness_multiplier'])
            baseline = diagnostics[(selected[0]['trial']['name'], model)]['baseline_search_score']
            values = [baseline / selected[0]['threshold']] + [r['search_score'] / r['threshold'] for r in selected]
            axis.plot([0, *LEVELS], values, 'o-', color=color, markersize=4, lw=1.4, label=f'r={site[0]} px, block {site[1]}')
        counts = [next(g['detections'] for g in data['groups'] if g['model'] == model and g['brightness_multiplier'] == level) for level in LEVELS]
        null = data['null_models'][model]
        axis.axhline(1, color='black', ls='--', lw=1)
        axis.set(title=f'{label}: recovered {counts[0]}, {counts[1]}, {counts[2]} of 6\nNull exceedances: {null["evaluation_exceedances"]}/28',
                 xticks=[0, *LEVELS], xlabel='Injected contrast / site-specific reference', ylabel='Search score / frozen threshold')
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8, ncol=2, loc='upper left')
    fig.suptitle('Step 5: six fixed sites, three positive brightness levels\nZero denotes the existing baseline; detection requires strict exceedance of one', fontsize=12)
    fig.savefig(output / 'recovery.png', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2), layout='constrained')
    for index, (model, color) in enumerate(zip(MODELS, COLORS)):
        values = [[], [], []]
        for site_index, site in enumerate(sites):
            for level, marker in zip(LEVELS, ('v', 'o', '^')):
                record = next(r for r in rows if r['model'] == model and r['trial']['brightness_multiplier'] == level and
                              (r['trial']['nominal_radius'], r['trial']['angular_block']) == site)
                diagnostic = diagnostics[(record['trial']['name'], model)]
                measures = (100 * record['contrast_bias_fraction'], diagnostic['standardized_error'],
                            100 * diagnostic['paired_increment_error_fraction'])
                x = index + (site_index - 2.5) * 0.075
                for axis, collection, value in zip(axes, values, measures):
                    axis.scatter(x, value, marker=marker, color=color, s=28, alpha=0.8)
                    collection.append(value)
        for axis, collection in zip(axes, values):
            axis.plot([index - 0.26, index + 0.26], [np.median(collection)] * 2, color='black', lw=1.5)
    axes[0].set(title='Raw contrast error', ylabel='100 × (measured / injected − 1) [%]')
    axes[1].axhspan(-1, 1, color='#d3e5d3', zorder=0)
    axes[1].set(title='Conditional intervals under-cover\nIdentity uses unit covariance', ylabel='(Measured − injected) / conditional sigma')
    axes[2].set(title='Secondary paired-increment diagnostic', ylabel='100 × [(positive − baseline) / injected − 1] [%]')
    for axis in axes:
        axis.axhline(0, color='0.5', lw=0.7)
        axis.set(xticks=range(4), xticklabels=('Identity', 'Diagonal', 'PCA 0', 'PCA 3'))
        axis.grid(axis='y', alpha=0.2)
    axes[0].legend(handles=[Line2D([], [], marker=m, linestyle='none', color='0.3', label=f'{level:g} × reference')
                           for level, m in zip(LEVELS, ('v', 'o', '^'))], fontsize=8)
    fig.suptitle('Step 5 photometry: 18 measurements per method, reusing six spatial sites\nPoints show individual measurements; black bars show pooled medians', fontsize=12)
    fig.savefig(output / 'photometry.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Verify a completed frozen experiment before producing its review archive."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.input.resolve(), args.output.resolve()
    require(root != output and root not in output.parents, 'review must be outside the frozen experiment')
    data, nulls, verification = audit(root)
    review = summarize(data, nulls, verification)
    review['verification']['review_script'] = fingerprint(Path(__file__))
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(root / 'results.json', output / 'results.json')
    (output / 'review.json').write_text(json.dumps(review, indent=2, allow_nan=False) + '\n')
    plot(data, review, output)
    print(json.dumps(verification, indent=2))


if __name__ == '__main__':
    main()
