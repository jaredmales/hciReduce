#!/usr/bin/env python3
"""Add fixed 3.6-pixel Gaussian references to the completed Step-5 image experiment.

Reuse the original null split, common eligibility, five-pixel searches, and positive
images. Freeze each reference threshold on calibration nulls only. This is an
extension requested after inspecting the original evaluation, not a new blind test.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.signal import convolve2d

from run_p4_step5_full_injections import fingerprint, write_json

REFERENCES = ('gaussian_raw', 'gaussian_snr', 'identity_snr')
OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))


def read_json(path: Path):
    """Read a saved experiment record."""
    return json.loads(path.read_text())


def verify(record: dict) -> None:
    """Refuse changed source inputs or frozen software."""
    if fingerprint(Path(record['path'])) != record:
        raise RuntimeError('input changed: ' + record['path'])


def score(image: np.ndarray, trial: dict) -> dict:
    """Measure exactly the existing five-pixel signed search maximum."""
    x, y = trial['row'], trial['column']
    values = [float(image[y + dy, x + dx]) for dx, dy in OFFSETS]
    if not all(np.isfinite(values)):
        raise RuntimeError('reference has invalid support at ' + trial['name'])
    peak = int(np.argmax(values))
    return {'search_score': values[peak], 'center_value': values[0],
            'peak_row_column': [x + OFFSETS[peak][0], y + OFFSETS[peak][1]]}


def gaussian_oracle(image: np.ndarray) -> np.ndarray:
    """Independently calculate the documented 15-square, mask-normalized Gaussian in FP64."""
    sigma = 3.6 / math.sqrt(8 * math.log(2))
    yy, xx = np.indices((15, 15))
    kernel = np.exp(-((xx - 7)**2 + (yy - 7)**2) / (2 * sigma**2))
    valid = np.isfinite(image)
    numerator = convolve2d(np.where(valid, image, 0), kernel, mode='same')
    denominator = convolve2d(valid.astype(float), kernel, mode='same')
    result = np.full(image.shape, np.nan)
    np.divide(numerator, denominator, out=result, where=valid & (denominator > 0))
    return result


def analyze_image(directory: Path, source: Path, software: Path, response: Path, environment: dict) -> dict:
    """Run the production Gaussian and ordinary SNR paths, checking the exported smoothing exactly."""
    directory.mkdir()
    science = directory / 'science.fits'
    shutil.copy2(source, science)
    commands = []

    def run(command: list[str], name: str) -> None:
        """Preserve every production command and its log."""
        commands.append(command)
        write_json(directory / 'commands.json', commands)
        with (directory / (name + '.log')).open('w') as log:
            subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)

    gaussian = directory / 'gaussian.fits'
    run(['taskset', '-c', '12,13', str(software / 'hciGaussianReference'), str(science), str(gaussian), '3.6'], 'smooth')
    base = ['taskset', '-c', '12,13', str(software / 'hciAnalyze'), '--lambdaD=3.6',
        '--planet.sep=11.782', '--planet.PA=262.051', '--planet.R=7.3', '--snr.apertureR=1',
        '--snr.minRad=0', '--snr.maxRad=60', '--filter.hpfGaussFW=0', '--noise.model=identity',
        '--noise.outputDiagnostics=false', '--noise.only=false']
    run([*base, f'--file={science}', '--filter.psfResponse=', '--filter.lpfGaussFW=3.6'], 'gaussian_snr')
    shutil.move(directory / 'science_snr.fits', directory / 'gaussian_snr.fits')
    replay = directory / 'replay.fits'
    shutil.copy2(gaussian, replay)
    run([*base, f'--file={replay}', '--filter.psfResponse=', '--filter.lpfGaussFW=0'], 'replay')
    actual = fits.getdata(directory / 'gaussian_snr.fits').squeeze()
    repeated = fits.getdata(directory / 'replay_snr.fits').squeeze()
    if not np.array_equal(actual, repeated, equal_nan=True):
        raise RuntimeError('exported Gaussian differs from production CLI smoothing')
    run([*base, f'--file={science}', f'--filter.psfResponse={response}', '--filter.lpfGaussFW=0'], 'identity_snr')
    shutil.move(directory / 'science_snr.fits', directory / 'identity_snr.fits')
    image = fits.getdata(science).squeeze()
    smooth = fits.getdata(gaussian).squeeze()
    expected = gaussian_oracle(image)
    finite = np.isfinite(expected)
    error = float(np.max(np.abs(smooth[finite] - expected[finite])))
    scale = float(np.max(np.abs(expected[finite])))
    if not np.array_equal(np.isfinite(smooth), finite) or error > 3e-6 * scale:
        raise RuntimeError('Gaussian export disagrees with independent mask-aware convolution')
    gaussian_header = fits.getheader(directory / 'gaussian_snr.fits')
    identity_header = fits.getheader(directory / 'identity_snr.fits')
    if (not np.isclose(gaussian_header['LPFGFW'], 3.6) or gaussian_header['HCIAPSF'].strip() or
        gaussian_header['SNRSMALL'] != 1 or identity_header['LPFGFW'] != 0 or
        identity_header['HCIAPSF'] != str(response)):
        raise RuntimeError('unexpected filtering metadata')
    result = {'source': fingerprint(source), 'gaussian_cli_replay_bitwise_equal': True,
        'independent_gaussian_max_absolute_error': error, 'independent_gaussian_relative_to_peak': error / scale,
        'products': [fingerprint(p) for p in sorted(directory.glob('*.fits'))]}
    write_json(directory / 'complete.json', result)
    return result


def maps(directory: Path) -> dict:
    """Load three different detection statistics without assigning them common units."""
    return {name: fits.getdata(directory / filename).squeeze() for name, filename in
            zip(REFERENCES, ('gaussian.fits', 'gaussian_snr.fits', 'identity_snr.fits'))}


def plot_comparison(result: dict, output: Path) -> None:
    """Compare Gaussian and identity recovery at their calibration-only thresholds."""
    models = ('gaussian_raw', 'identity', 'gaussian_snr', 'identity_snr')
    labels = ('Gaussian 3.6 px: smoothed intensity', 'Identity: conditional matched-filter score',
              'Gaussian 3.6 px: application annular SNR', 'Identity: application annular SNR')
    sites = sorted({(r['trial']['nominal_radius'], r['trial']['angular_block']) for r in result['injections']})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True, layout='constrained')
    for axis, model, label in zip(axes.flat, models, labels):
        for site, color in zip(sites, plt.get_cmap('tab10').colors):
            selected = sorted([r for r in result['injections'] if r['model'] == model and
                (r['trial']['nominal_radius'], r['trial']['angular_block']) == site],
                key=lambda r: r['trial']['brightness_multiplier'])
            baseline = next(t for t in result['null_trials'] if t['name'] == selected[0]['trial']['null_trial'])
            threshold = result['models'][model]['threshold']
            axis.plot([0, 0.5, 1, 2], [baseline['models'][model]['search_score'] / threshold] +
                [r['search_score'] / threshold for r in selected], 'o-', color=color, markersize=4,
                label=f'r={site[0]} px, block {site[1]}')
        counts = result['models'][model]['detections_by_brightness']
        axis.axhline(1, color='black', ls='--', lw=1)
        axis.set(title=f'{label}\nRecovery: {counts}; nulls: {result["models"][model]["evaluation_exceedances"]}/28',
                 xticks=[0, 0.5, 1, 2], xlabel='Injected contrast / original site-specific reference',
                 ylabel='Search statistic / calibration threshold')
        axis.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8, ncol=2, loc='upper left')
    fig.suptitle('Step 5 extension: Gaussian reference at the user-specified FWHM\nSame images and trials; added after inspecting the original evaluation', fontsize=12)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def main() -> None:
    """Freeze the added reference, run saved-image analyses, and archive the comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation', type=Path, required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    args = parser.parse_args()
    original, root = args.evaluation.resolve(), args.output.resolve()
    protocol = read_json(original / 'protocol.json')
    nulls = read_json(original / 'null_results.json')
    previous = read_json(original / 'results.json')
    jobs = read_json(original / 'jobs.json')
    original_manifest = read_json(original / 'manifest.json')
    inputs = [protocol['science'], *protocol['source_response_inputs'],
              *original_manifest['frozen_records'], read_json(original / 'complete.json')['results']]
    inputs += [read_json(original / 'reductions' / j['name'] / 'complete.json')['products'][0] for j in jobs]
    for record in inputs:
        verify(record)
    root.mkdir(parents=True, exist_ok=False)
    software = root / 'software'
    software.mkdir()
    for name in ('hciAnalyze', 'libhcireduce.so', 'libmxlib.so', 'liblapack.so.3', 'libopenblas.so.0'):
        shutil.copy2(original / 'software' / name, software / name)
    shutil.copy2(args.helper, software / 'hciGaussianReference')
    for name in (Path(__file__).name, 'run_p4_step5_full_injections.py'):
        shutil.copy2(Path(__file__).with_name(name), software / name)
    repository = Path(__file__).resolve().parents[3]
    for name in ('benchmarks/hciGaussianReference.cpp', 'src/apps/hciAnalyze.cpp'):
        shutil.copy2(repository / name, software / Path(name).name)
    response = Path(next(r['path'] for r in protocol['source_response_inputs'] if r['path'].endswith('_manifest.fits')))
    frozen = [fingerprint(p) for p in sorted(software.iterdir())]
    manifest = {'purpose': 'post-review fixed Gaussian reference; no FWHM tuning and no new reductions',
        'fwhm_pixels': 3.6, 'kernel_size_pixels': 15, 'references': list(REFERENCES),
        'search_offsets_row_column': OFFSETS, 'common_eligibility': 'original 28 calibration and 28 evaluation searches',
        'threshold_rule': 'ceil((n+1)*0.95) order statistic; strict exceedance',
        'snr_reference': 'ordinary application annular mean/stddev and small-sample correction; known source masked at radius 7.3; trial neighborhoods are NOT excluded from this legacy normalization',
        'raw_reference': 'production Gaussian-smoothed intensity, no data-fitted noise normalization; raw scale is arbitrary',
        'inputs': inputs, 'software': frozen, 'original_protocol': fingerprint(original / 'protocol.json'),
        'original_nulls': fingerprint(original / 'null_results.json'), 'original_results': fingerprint(original / 'results.json')}
    write_json(root / 'manifest.json', manifest)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LD_LIBRARY_PATH=str(software))
    environment.pop('HCIANALYZE_GLOBAL_CONFIG', None)
    audits = [analyze_image(root / 'baseline', Path(protocol['science']['path']), software, response, environment)]
    baseline_maps = maps(root / 'baseline')
    # Freeze thresholds before processing any positive images.
    trials = []
    for trial in nulls['trials']:
        record = {**trial, 'models': dict(trial['models'])}
        for model in REFERENCES:
            record['models'][model] = score(baseline_maps[model], trial) if trial['common_eligible'] else None
        trials.append(record)
    summaries = dict(nulls['models'])
    for model in REFERENCES:
        calibration = [t['models'][model]['search_score'] for t in trials if t['role'] == 'calibration' and t['common_eligible']]
        rank = math.ceil((len(calibration) + 1) * 0.95)
        threshold = sorted(calibration)[rank - 1]
        evaluation = [t for t in trials if t['role'] == 'evaluation' and t['common_eligible']]
        summaries[model] = {'threshold': threshold, 'threshold_rank': rank, 'calibration_trials': len(calibration),
            'calibration_exceedances': sum(s > threshold for s in calibration), 'evaluation_trials': len(evaluation),
            'evaluation_exceedances': sum(t['models'][model]['search_score'] > threshold for t in evaluation)}
    write_json(root / 'thresholds.json', summaries)
    threshold_hash = fingerprint(root / 'thresholds.json')
    injections = [{k: r[k] for k in ('trial', 'model', 'search_score', 'detected')} for r in previous['measurements']]
    for job in jobs:
        directory = root / job['name']
        audits.append(analyze_image(directory, original / 'reductions' / job['name'] / 'finim.fits', software, response, environment))
        for model, image in maps(directory).items():
            measurement = score(image, job)
            injections.append({'trial': job, 'model': model, **measurement,
                'detected': measurement['search_score'] > summaries[model]['threshold']})
        print(job['name'], flush=True)
    for model in summaries:
        summaries[model] = {**summaries[model], 'detections_by_brightness': [sum(r['detected'] for r in injections
            if r['model'] == model and r['trial']['brightness_multiplier'] == level) for level in (0.5, 1, 2)]}
    for record in [*inputs, *frozen, threshold_hash]:
        verify(record)
    result = {'manifest': fingerprint(root / 'manifest.json'), 'models': summaries, 'null_trials': trials,
              'injections': injections, 'production_audits': audits,
              'caveat': 'All comparisons are descriptive on the same correlated field. Legacy annular SNR normalizes using trial neighborhoods too; it is not the held-out noise estimator used for covariance.'}
    write_json(root / 'results.json', result)
    write_json(root / 'complete.json', {'images': len(audits), 'new_injection_measurements': len(jobs) * len(REFERENCES),
        'all_gaussian_cli_replays_bitwise_equal': True, 'frozen_inputs_unchanged': True, 'results': fingerprint(root / 'results.json')})
    args.archive.mkdir(parents=True, exist_ok=False)
    for name in ('manifest.json', 'results.json', 'complete.json', 'thresholds.json'):
        shutil.copy2(root / name, args.archive / name)
    plot_comparison(result, args.archive / 'comparison.png')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
